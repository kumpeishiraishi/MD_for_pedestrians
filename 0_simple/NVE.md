分子動力学法（molecular dynamics; MD）は、計算機で物理現象を解析する強力な手法である。ここでは、MDでよく用いられるアルゴリズムを導出してみよう。

粒子の1次元の運動を考える（質量は $m$ ）。時刻 $t$ での粒子の位置を $x(t)$ 、速度を $v(t) = \dot{x}(t)$ 、加速度を $a(t) = \ddot{x}(t)$ と表記する。 $x(t + \Delta t)$ を $t$ の近くでテイラー展開すると

$$
x(t + \Delta t) = x(t) + \Delta t v(t) + \frac{\Delta t^2}{2} a(t) + \cdots
$$

となる（以下、 $\Delta t^3$ 以上の項は考えないことにする）。同様に、 $x(t)$ を $t + \Delta t$ の近くで展開することも出来る。すると、

$$
x(t) = x(t + \Delta t) - \Delta t v(t + \Delta t) + \frac{\Delta t^2}{2} a(t + \Delta t)
$$

と書くことができる。これらの両辺を足し合わせて整理すると、

$$
v(t + \Delta t) = v(t) + \frac{\Delta t}{2} \left[a(t) + a(t + \Delta t)\right]
$$

となる。

ところで、運動方程式 $ma(t) = F(t)$ から、加速度 $a(t)$ は粒子にかかる力 $F(t)$ を用いて

$$
a(t) = \frac{F(t)}{m}
$$

と表すことができる。これを使って、上のテイラー展開を書き直すと

$$
x(t + \Delta t) = x(t) + \Delta t v(t) + \frac{\Delta t^2}{2m} F(t)
$$

$$
v(t + \Delta t) = v(t) + \frac{\Delta t}{2m} \left[F(t) + F(t + \Delta t)\right]
$$

と書ける。少しだけ式変形すると

$$
x(t + \Delta t) = x(t) + \Delta t \left[ v(t) + \frac{\Delta t}{2m} F(t) \right]
$$

$$
v(t + \Delta t) = \left[ v(t) + \frac{\Delta t}{2m} F(t) \right] + \frac{\Delta t}{2m} F(t + \Delta t)
$$

と、鉤括弧内に両式で共通した部分が現れる。

これが、速度べルレ法（velocity Verlet）というアルゴリズムでよく使われる形である。改めてアルゴリズムを整理すると、速度ベルレ法は

1.  $v\left(t + \frac{\Delta t}{2}\right) = v(t) + \frac{\Delta t}{2m}F(t)$ で速度を更新する
2.  $x(t + \Delta t) = x(t) + \Delta t v\left(t + \frac{\Delta t}{2}\right)$ で位置を更新する
3.  新しい時刻での位置 $x(t + \Delta t)$ で力 $F(t + \Delta t)$ を計算する
4.  $v(t + \Delta t) = v\left(t + \frac{\Delta t}{2}\right) + \frac{\Delta t}{2m}F(t + \Delta t)$ で速度を更新する
5.  ステップ1に戻る

と、まとめられる。これまで単一の粒子による1次元の運動を考えてきたが、3次元空間中の $N$ 粒子の運動（物質の研究をするときに一番自然な設定）に変更するのは簡単なはずだ。

さて、これは運動方程式を解いていることに他ならない。そのときには当然エネルギー保存則が成り立つので、速度ベルレ法によって数値積分された系はミクロカノニカル分布に従った状態を取っていく。つまり、粒子数 $N$ 、体積 $V$ 、エネルギー $E$ が一定の状態を時々刻々発生させていく[^1]ため、速度ベルレ法は「 *NVE* ダイナミクス」と呼ばれることも多い。

[^1]: このような一連のミクロ状態の時系列を「トラジェクトリー」や「（相空間中の）軌道」などと呼ぶ
