分子動力学シミュレーションの参考になりそうなコード。
目的は、[LAMMPS](https://github.com/lammps/lammps)とか[HOOMD-blue](https://github.com/glotzerlab/hoomd-blue)とか、オープンソースだけど巨大なコードを解読するのはツライ、という初学者の修養の用に供することである。

# 系

なるべくプログラムを単純にするため、このリポジトリでは一つの系に限定したコードを書く。
系の相互作用はLennard-Jonesポテンシャルとする。
クーロン力や重力などの長距離相互作用は難しいので、カットオフ付きの短距離のLJ相互作用を扱う。

粒子 $i$ と粒子 $j$ の距離を $r_{ij}$ とすると、この2粒子間の相互作用は以下の式で表される：

$$
V(r_{ij}) = \phi(r_{ij}) - \phi(r_{ij}^\mathrm{cut}) - \phi^\prime(r_{ij}^\mathrm{cut}) (r_{ij} - r^\mathrm{cut}_{ij}).
$$

但し、

$$
\phi(r_{ij}) = 4\epsilon_{ij} \left[ \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{12} - \left( \frac{\sigma_{ij}}{r_{ij}} \right)^6 \right]
$$

である。

3次元系を考える。系に含まれる粒子数は $N$ とし、シミュレーションの箱の大きさ $L$ は数密度 $\rho = N/L^3 = 1.2$ より定める。
各方向で周期境界条件を採用する（下の説明も参照）。

この系には、LJパラメータ $\sigma_{ij}$ / $\epsilon_{ij}$ の異なる2成分の粒子が含まれている。
それぞれA粒子とB粒子と呼び、粒子数は80:20の割合である。
どちらも粒子の質量は $m$ とする。 各LJパラメータは以下の通り：

  |            | AA  | AB  | BB   |
  |------------|-----|-----|------|
  | $\sigma$   | 1.0 | 0.8 | 0.88 |
  | $\epsilon$ | 1.0 | 1.5 | 0.5  |

上の表式を見れば分かるように、このポテンシャルはカットオフ距離 $r_{ij}^\mathrm{cut}$ でポテンシャルと力（ポテンシャルの一階微分）が連続になるようにシフトされている。
カットオフ距離は、

-   AAとBBのとき $1.5\sigma_{AA} = 1.5$
-   ABのとき $2.5\sigma_{AB} = 2.0$

とする。 各粒子種ペアのポテンシャルをグラフにすると、こんな感じ
![](./potential.png)

質量、長さ、エネルギー、時間の単位は、それぞれ $m$ 、 $\sigma_{AA}$ 、 $\epsilon_{AA}$ 、 $(m\sigma_{AA}^2/\epsilon_{AA})^{1/2}$ とする。

ちなみに、ここで導入した[2成分Lennard-Jones粒子系](https://doi.org/10.1063/5.0004093)は、Kob-Andersenモデルと呼ばれる過冷却液体の最も基本的なモデルの修正版である。
[オリジナルのKob-Andersenモデル](https://doi.org/10.1103/physreve.51.4626)とは若干異なるが、オリジナル版は低温で結晶化しやすいという難点があるため、個人的好みから修正版を用いることにする。

# 参考文献
- Allen, Tildesley, 'Computer Simulation of Liquids', Oxford University Press
- Frenkel, Smit, 'Understanding Molecular Simulation', Academic Press
