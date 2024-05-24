from CST import CST  # CSTモジュールからCSTクラスをインポート

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def _update_left(i, ax_left, res, airfoil):
    """
    左側のプロットを更新する関数
    Args:
    - i: 現在のフレーム番号
    - ax_left: 左側のプロットの軸オブジェクト
    - res: 最適化結果のデータフレーム
    - airfoil: 翼型形状データ
    """
    ax_left.clear()  # 現在描写されているグラフを消去
    ax_left.set_title("cd = {:.4f} cm = {:.4f}".format(res['cd'].iloc[i], -res['cm'].iloc[i]))  # タイトルに現在のcdとcmの値を表示
    ax_left.set_xlim(0, 1)  # x軸の範囲を設定
    ax_left.set_ylim(-0.5, 0.5)  # y軸の範囲を設定
    ax_left.plot(airfoil.x[i], airfoil.y[i], color='blue')  # 翼型形状をプロット

def _update_right(i, ax_right, res):
    """
    右側のプロットを更新する関数
    Args:
    - i: 現在のフレーム番号
    - ax_right: 右側のプロットの軸オブジェクト
    - res: 最適化結果のデータフレーム
    """
    ax_right.clear()  # 現在描写されているグラフを消去
    ax_right.set_title('cd vs cm')  # タイトルを設定
    ax_right.set_ylim(0, 0.012)  # y軸の範囲を設定
    ax_right.set_xlim(-0.1, 0.02)  # x軸の範囲を設定
    ax_right.set_xticks([-0.20, -0.16, -0.12, -0.08, -0.04, 0.00, 0.04])  # x軸のメモリを設定
    ax_right.set_ylabel('cd')  # y軸のラベルを設定
    ax_right.set_xlabel('cm')  # x軸のラベルを設定
    ax_right.grid()  # グリッドを表示
    ax_right.plot(-res['cm'], res['cd'], '.', color='blue')  # cd vs cmのプロット
    ax_right.plot(-res['cm'].values[i], res['cd'].values[i], '*', color='red', markersize=10)  # 現在のフレームのデータを強調表示

def gif_animation(res, airfoil):
    """
    GIFアニメーションを作成する関数
    Args:
    - res: 最適化結果のデータフレーム
    - airfoil: 翼型形状データ
    """
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12, 6))  # 左右に2つのプロットを持つ図を作成

    def _update(i):
        _update_left(i, ax_left, res, airfoil)  # 左側のプロットを更新
        _update_right(i, ax_right, res)  # 右側のプロットを更新

    ani = animation.FuncAnimation(fig, _update, interval=100, frames=len(res), repeat_delay=2000)  # アニメーションを作成

    # 保存
    ani.save("./MDO/sample.gif", writer="pillow")  # アニメーションをGIFとして保存

def objective_space(res):
    """
    目的空間のプロットを表示する関数
    Args:
    - res: 最適化結果のデータフレーム
    """
    fig = plt.figure(figsize=(6, 6))  # 図のサイズを指定して作成

    plt.title('cd vs cm')  # タイトルを設定
    plt.ylabel('cd')  # y軸のラベルを設定
    plt.xlabel('cm')  # x軸のラベルを設定
    plt.ylim(0, 0.012)  # y軸の範囲を設定
    plt.xlim(-0.20, 0.04)  # x軸の範囲を設定
    plt.xticks([-0.20, -0.16, -0.12, -0.08, -0.04, 0.00, 0.04])  # x軸の目盛りを設定
    plt.grid()  # グリッドを表示
    plt.plot(-res['cm'], res['cd'], '.', color='blue')  # cd vs cmのプロット

    # 表示
    plt.show()  # プロットを表示

# 最適化結果の読み込み
res_X = pd.read_csv('./MDO/res_X.csv')
res_F = pd.read_csv('./MDO/res_F.csv')

# 結果を結合し、cdでソート
res = pd.concat([res_X, res_F], axis=1)
res = res.sort_values('cd')
res_np = res.to_numpy()

# 翼型形状を生成
airfoil = CST()
n_wu = int(res_X.shape[1]/2)
wu = res_np[:, :n_wu].tolist()
wl = res_np[:, n_wu:].tolist()
dz = np.zeros(len(wu)).tolist()
airfoil.create_airfoil(wu=wu, wl=wl, dz=dz, Node=100)

# 目的空間のプロットを表示
objective_space(res)
# GIFアニメーションを作成
gif_animation(res, airfoil)