# DCBR v2
タンパク質の残基間の速度相関関数を解析するためのプログラムである。

## 速度相関関数の定義
<img src="https://latex.codecogs.com/gif.latex?C_{i,j}(n):=\left\langle\cfrac{\vec{v}_i(t)\cdot\vec{v}_j(t+n\Delta t)}{|\vec{v}_i(t)| |\vec{v}_j(t+n\Delta t)|}\right\rangle_t"/>

- <img src="https://latex.codecogs.com/gif.latex?\vec{v}_i:=\vec{u}_i/\Delta t"/>
- <img src="https://latex.codecogs.com/gif.latex?\vec{u}_i:=\vec{p}_i(t-\Delta t)-\vec{p}_i(t)"/>
- <img src="https://latex.codecogs.com/gif.latex?\vec{p}_i"/> : <img src="https://latex.codecogs.com/gif.latex?i"/> 番目の残基のうち<img src="https://latex.codecogs.com/gif.latex?C_\alpha"/>を除いた全ての原子の位置
- <img src="https://latex.codecogs.com/gif.latex?t"/> : サンプリング対象の時刻
- <img src="https://latex.codecogs.com/gif.latex?\Delta t"/> : 入力ファイル内の最小の時間間隔
- <img src="https://latex.codecogs.com/gif.latex?\langle\cdots\rangle_t"/> : <img src="https://latex.codecogs.com/gif.latex?t"/> で平均
- <img src="https://latex.codecogs.com/gif.latex?|\vec{a}|"/> : <img src="https://latex.codecogs.com/gif.latex?\vec{a}"/> のノルム(L2)

## 使用
dcbr

in=入力.pdb

out=出力ファイル名

step=サンプリング開始時刻-終了時刻

target=対象残基は何番めの残基か？(i)

delta=n_max (n=[0:n_max]で解析する)

オプション --silence

を１行で実行すること。

## ビルド
例: clang++-5.0 -std=c++1z -stdlib=libc++ main.cpp -O2 -o dcbr

また、コードは boost を利用しているため、設定していないのなら、それについてパスを通す必要がある。

boost を macports でインストールしているのなら、-I/opt/local/include/ を上記のものに追加する。
