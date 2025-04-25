#include <iostream>
#include <iomanip>

#include <chrono>

#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <vector>
//---------------------------------------------------------------
constexpr auto N = 500;         // 粒子数
//---------------------------------------------------------------
constexpr auto deg = 3;         // 空間次元
constexpr auto dt  = 5e-3;      // 時間刻み
constexpr auto N_A = N*4/5;     // A粒子の数（80:20のKob-Andersen mixture）
//---------------------------------------------------------------
constexpr auto rho = 1.2;       // 数密度
const double Lbox = std::pow(N/rho, 1.0/deg);
                                // シミュレーションボックスの一辺の大きさ
const double Linv = 1.0/Lbox;   // Lboxの逆数
//---------------------------------------------------------------
double conf[N][deg], velo[N][deg], force[N][deg];
                                // 粒子の位置、粒子の速度、粒子にかかる力を保持する配列
enum {X, Y, Z};
//---------------------------------------------------------------
// カットオフ距離からのマージン
constexpr auto margin = 0.7;
// マージンの距離の2乗
constexpr auto SKIN2  = (margin*0.5) * (margin*0.5);
// 相互作用ペアリスト（余裕を見て、1粒子当たり100粒子まで相互作用しても足りるように）
int list[N*100];
// 各粒子は list のどこからどこまでを見れば良いか、先頭インデックスを記録しておく配列
int point[N];
// 隣接リスト構築時点での配置を保存しておく配列
double NL_config[N][deg];
//---------------------------------------------------------------
void init_lattice() {
    const auto ln   = std::ceil(std::pow(N, 1.0/deg));
    const auto haba = Lbox/ln;

    for (int i=0; i<N; i++) {
        const int iz = std::floor(i/(ln*ln));
        const int iy = std::floor((i - iz*ln*ln)/ln);
        const int ix = i - iz*ln*ln - iy*ln;

        // 正方格子に粒子を配置
        conf[i][X] = haba*0.5 + haba * ix;
        conf[i][Y] = haba*0.5 + haba * iy;
        conf[i][Z] = haba*0.5 + haba * iz;

        // minimum-image convention
        for (int d=0; d<deg; d++) {
            conf[i][d] -= Lbox * std::floor(conf[i][d] * Linv + 0.5);
        }
    }
}
void init_species(std::mt19937 &mt) {
    // 0,1,...,N-1が格納されたstd::vector配列を初期化
    std::vector<int> v(N);
    std::iota(v.begin(), v.end(), 0);

    // その配列をシャッフル
    std::shuffle(v.begin(), v.end(), mt);

    for (int i=0; i<N; i+=2) {
        // シャッフルされた配列内で隣接する粒子インデックスペアを使って
        // 粒子位置を交換する
        const int id0 = v[i];
        const int id1 = v[i+1];

        // 交換前の位置を一時的に保存しておく        
        const double position0_X = conf[id0][X];
        const double position0_Y = conf[id0][Y];
        const double position0_Z = conf[id0][Z];

        const double position1_X = conf[id1][X];
        const double position1_Y = conf[id1][Y];
        const double position1_Z = conf[id1][Z];

        // 交換する
        conf[id0][X] = position1_X;
        conf[id0][Y] = position1_Y;
        conf[id0][Z] = position1_Z;
        conf[id1][X] = position0_X;
        conf[id1][Y] = position0_Y;
        conf[id1][Z] = position0_Z;
    }
}
//---------------------------------------------------------------
inline void remove_drift() {
    double vel1 = 0.0, vel2 = 0.0, vel3 = 0.0;

    // 系全体の速度を計算する
    for (int i=0; i<N; i++) {
        vel1 += velo[i][X];
        vel2 += velo[i][Y];
        vel3 += velo[i][Z];
    }
    vel1 /= N;
    vel2 /= N;
    vel3 /= N;
    // 各粒子の速度から、系全体の速度/Nを引いておく
    for (int i=0; i<N; i++) {
        velo[i][X] -= vel1;
        velo[i][Y] -= vel2;
        velo[i][Z] -= vel3;
    }
}
void init_vel_MB(const double T_targ, std::mt19937 &mt) {
    std::normal_distribution<double> dist_trans(0.0, std::sqrt(T_targ));
    for (int i=0; i<N; i++) {
        velo[i][X] = dist_trans(mt);
        velo[i][Y] = dist_trans(mt);
        velo[i][Z] = dist_trans(mt);
    }
    remove_drift();
}
//---------------------------------------------------------------
// LJパラメータを保存しておくルックアップ・テーブル
constexpr double MBLJ_sij1[2][2] = {
    {1.0, 0.8},
    {0.8, 0.88}
};
constexpr double MBLJ_energy[2][2] = {
    {1.0, 1.5},
    {1.5, 0.5}
};
//---------------------------------------------------------------
// LJポテンシャルとその一階微分
double LJpotential(const double rij1, const double sij1) {
    const double rij2 = rij1 * rij1;
    const double rij6 = rij2 * rij2 * rij2;
    const double sij2 = sij1 * sij1;
    const double sij6 = sij2 * sij2 * sij2;
    return 4.0 * sij6 * (sij6 - rij6)/(rij6 * rij6);
}
double deriv_1st_LJpotential(const double rij1, const double sij1) {
    const double rij2 = rij1 * rij1;
    const double rij6 = rij2 * rij2 * rij2;
    const double sij2 = sij1 * sij1;
    const double sij6 = sij2 * sij2 * sij2;
    return -24.0/rij1 * sij6 * (2.0 * sij6 - rij6)/(rij6 * rij6);
}
//---------------------------------------------------------------
void generate_NL() {
    // 系全体でのijペアの個数を数える変数
    auto nlist = -1;

    for (int i=0; i<N; i++) {
        // 粒子iのijペア番号の先頭を記録しておく
        point[i] = nlist+1;

        const int si = i>=N_A;
        for (int j=i+1; j<N; j++) {
            const int sj = j>=N_A;

            // 距離 rij を計算
            double dx = conf[i][X] - conf[j][X];
            double dy = conf[i][Y] - conf[j][Y];
            double dz = conf[i][Z] - conf[j][Z];
            dx -= Lbox * std::floor(dx * Linv + 0.5);
            dy -= Lbox * std::floor(dy * Linv + 0.5);
            dz -= Lbox * std::floor(dz * Linv + 0.5);
            const double rij2 = dx*dx + dy*dy + dz*dz;

            // カットオフ距離を計算
            const double rc1 = (si==sj ? 1.5 : 2.0);

            // カットオフ距離 + marginの距離を計算
            const double rlist1 = rc1 + margin;
            const double rlist2 = rlist1 * rlist1;

            // 距離 rij が閾値 rlist 以下かを判定
            if (rij2 < rlist2) {
                // 現在の粒子ペアijを相互作用ペアとしてカウント
                nlist++;
                // jのインデックスを相方粒子として保存する
                list[nlist] = j;
            }
        }
    }

    // 隣接リスト構築時点の配置 conf を NL_config で記憶しておく
    std::copy(*conf, *conf + N*deg, *NL_config);
}
//---------------------------------------------------------------
void calc_force() {
    // 力を計算する前に force をゼロ埋め
    std::fill(*force, *force + N*deg, 0.0);

    // 力の計算のループ
    for (int i=0; i<N-1; i++) {
        // 粒子インデックスi+1に対して、相方の粒子jの先頭インデックスを取得
        int pend = point[i+1];
        // 粒子i+1の相方と粒子iの相方が同じ（つまり、粒子iと相互作用する粒子jがいない）とき、
        // 粒子iをスキップ
        if (pend == point[i]) continue;

        // 粒子 i の粒子種を判定
        const int si = i>=N_A;

        for (int p=point[i]; p<pend; p++) {
            // 粒子jのインデックスを取得
            const int j = list[p];

            // 粒子 j の粒子種を判定
            const int sj = j>=N_A;

            // ij 間の距離の2乗を計算
            double dx = conf[i][X] - conf[j][X];
            double dy = conf[i][Y] - conf[j][Y];
            double dz = conf[i][Z] - conf[j][Z];
            dx -= Lbox * std::floor(dx * Linv + 0.5);
            dy -= Lbox * std::floor(dy * Linv + 0.5);
            dz -= Lbox * std::floor(dz * Linv + 0.5);
            const double rij2 = dx*dx + dy*dy + dz*dz;

            // カットオフ距離を計算
            // 粒子iと粒子jが同種粒子なら1.5
            // 粒子iと粒子jが異種粒子なら2.0
            const double rc1 = (si==sj ? 1.5 : 2.0);

            // 粒子iと粒子jが相互作用するか判定
            if (rij2 < rc1 * rc1) {
                const double rij1 = std::sqrt(rij2);
                const double sij1 = MBLJ_sij1[si][sj];

                double deriv_1st = deriv_1st_LJpotential(rij1, sij1) - deriv_1st_LJpotential(rc1, sij1);
                deriv_1st *= MBLJ_energy[si][sj];

                force[i][X] -= deriv_1st * dx / rij1;
                force[i][Y] -= deriv_1st * dy / rij1;
                force[i][Z] -= deriv_1st * dz / rij1;
                force[j][X] += deriv_1st * dx / rij1;
                force[j][Y] += deriv_1st * dy / rij1;
                force[j][Z] += deriv_1st * dz / rij1;
            }
        }
    }
}
//---------------------------------------------------------------
double calc_potential() {
    double ans = 0.0;
    for (int i=0; i<N-1; i++) {
        const int si = i>=N_A;
        const int pend = point[i+1];

        for (int p=point[i]; p<pend; p++) {
            const int j = list[p];
            const int sj = j>=N_A;

            double dx = conf[i][X] - conf[j][X];
            double dy = conf[i][Y] - conf[j][Y];
            double dz = conf[i][Z] - conf[j][Z];
            dx -= Lbox * std::floor(dx * Linv + 0.5);
            dy -= Lbox * std::floor(dy * Linv + 0.5);
            dz -= Lbox * std::floor(dz * Linv + 0.5);

            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double rc1 = (si==sj ? 1.5 : 2.0);
            if (rij2 < rc1 * rc1) {
                const double rij1 = std::sqrt(rij2);
                const double sij1 = MBLJ_sij1[si][sj];
                double potential = LJpotential(rij1, sij1) - LJpotential(rc1, sij1) - deriv_1st_LJpotential(rc1, sij1) * (rij1 - rc1);
                ans += potential * MBLJ_energy[si][sj];
            }
        }
    }
    return ans;
}
//---------------------------------------------------------------
inline void velocity_update() {
    for (int i=0; i<N; i++) {
        velo[i][X] += 0.5 * dt * force[i][X];
        velo[i][Y] += 0.5 * dt * force[i][Y];
        velo[i][Z] += 0.5 * dt * force[i][Z];
    }
}
inline void position_update() {
    for (int i=0; i<N; i++) {
        conf[i][X] += dt * velo[i][X];
        conf[i][Y] += dt * velo[i][Y];
        conf[i][Z] += dt * velo[i][Z];
    }
}
inline void NL_check() {
    // 前回、隣接リストを構築した配置と比べて
    // 変位が最大の2粒子を探す
    double max1st = 0.0, max2nd = 0.0;
    for (int i=0; i<N; i++) {
        double dx = conf[i][X] - NL_config[i][X];
        double dy = conf[i][Y] - NL_config[i][Y];
        double dz = conf[i][Z] - NL_config[i][Z];
        dx -= Lbox * std::floor(dx * Linv + 0.5);
        dy -= Lbox * std::floor(dy * Linv + 0.5);
        dz -= Lbox * std::floor(dz * Linv + 0.5);

        const double rij2 = dx*dx + dy*dy + dz*dz;
        if (rij2 > max1st) {
            max2nd = max1st;
            max1st = rij2;
        }
    }

    // 変位の和がmarginを越えていたら隣接リストを作り直す
    if (max1st + max2nd + 2*std::sqrt(max1st*max2nd) > margin*margin) {
        generate_NL();
    }
}
//---------------------------------------------------------------
void print_energies(const long t) {
    // ポテンシャルエネルギーを計算
    const double U = calc_potential();

    // 運動エネルギーを計算
    double K = 0.0;
    for (int i=0; i<N; i++) {
        K += 0.5 * (velo[i][X]*velo[i][X]
                    + velo[i][Y]*velo[i][Y]
                    + velo[i][Z]*velo[i][Z]);
    }

    // 時刻、1粒子当たりの運動エネルギー、1粒子当たりのポテンシャルエネルギー、1粒子当たりの全エネルギーを出力
    std::cout << std::setprecision(15) << std::scientific
              << dt * t << ","
              << K/N << ","
              << U/N << ","
              << (K + U)/N << std::endl;
}
//---------------------------------------------------------------
void NVE(const double tsim) {
    calc_force();

    // 周期境界条件の下で、何個目の箱のミラーに位置しているのか保存する配列
    int box[N][deg];
    std::fill(*box, *box + N*deg, 0);

    // log間隔で出力するために使う変数たち //////////
    const auto logbin = std::pow(10.0, 1.0/9);
    int counter = 5;
    auto checker = 1e-3 * std::pow(logbin, counter);
    ////////////////////////////////////////////

    long t = 0;
    print_energies(t);

    const long steps = tsim/dt;
    while (t < steps) {
        velocity_update();
        position_update();

        // PBC
        for (int i=0; i<N; i++) {
            double xbox = std::floor(conf[i][X] * Linv + 0.5);
            double ybox = std::floor(conf[i][Y] * Linv + 0.5);
            double zbox = std::floor(conf[i][Z] * Linv + 0.5);
            box[i][X] += xbox;
            box[i][Y] += ybox;
            box[i][Z] += zbox;
            conf[i][X] -= Lbox * xbox;
            conf[i][Y] -= Lbox * ybox;
            conf[i][Z] -= Lbox * zbox;
        }

        // 更新した位置においてマージンを飛び出した粒子がいないか判定
        NL_check();

        calc_force();
        velocity_update();

        t++;
        if (dt*t > checker) {
            checker *= logbin;
            print_energies(t);
        }
    }
}
//---------------------------------------------------------------
int main() {
    // 疑似乱数生成器を適当に初期化
    std::mt19937 mt(123456789);

    // 正方格子に粒子を配置する
    init_lattice();

    // A粒子とB粒子の初期位置はランダムに混ぜておく
    init_species(mt);

    // 初期化した粒子位置で隣接リストを構築する
    generate_NL();

    // 初期速度は、温度T = 1.0のマクスウェル・ボルツマン分布から引っ張ってくる
    init_vel_MB(1.0, mt);

    // t = 1e5の長さのNVEシミュレーションを実行（時間を測りながら）
    auto start = std::chrono::system_clock::now();
    NVE(1e2);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    // 実行時間を出力
    std::cout << elapsed << std::endl;

    // 終配置を出力
    for (int i=0; i<N; i++) {
        for (int d=0; d<deg; d++) {
            std::cout << std::setprecision(15) << std::scientific
                      << conf[i][d] << std::endl;
        }
    }
    
}
