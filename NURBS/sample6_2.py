import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

from mybspline import MyBSpline


if __name__=='__main__':
    np.set_printoptions(linewidth=999)

    # --- Bスプラインの次数と制御点数の設定 ---
    k = 3           # Bスプラインの次数（4次＝クインティック）
    n_cp = 13       # 制御点数（曲線の自由度）
    n_sample = 100

    # --- キンク（曲率変化点）の指定 ---
    # 曲線の形状を変化させたい位置を指定（例：0.3）
    kink_list = []
    kink_list = [(0.3, k)]

    # --- 新しいBスプラインオブジェクトの初期化 ---
    bspline = MyBSpline()

    # --- サンプル点の定義（y軸に沿って √y の曲線） ---
    x = np.linspace(0, 1, 31)
    y = x ** 2
    z = 0.1 * (y ** 2)
    sample_points = np.column_stack([x, y, z])  # 3次元座標

    # --- 補間点の定義（√yの代表点） ---
    u_interp = np.array([0, 0.3, 1.0])
    interpolate_points = np.column_stack([np.sqrt(u_interp), u_interp, 0.1 * u_interp])

    # --- ノットベクトルの生成（キンク位置を反映） ---
    knot_vector = bspline._calc_knot_vector(n_cp, k, kink_list=kink_list)

    # --- ノットベクトルからキンク（重複）を除去 ---
    knot_vector, knot_to_add = bspline._remove_kink_from_knot(knot_vector, k)
    print('knot_vector = ', '['+', '.join([f'{val:.3f}' for val in knot_vector])+']')

    # --- 制御点数をノットベクトルから再計算 ---
    n_cp = len(knot_vector) - (k + 1)

    # --- Bスプラインのフィッティング（補間＋近似） ---
    control_points = bspline.fit_Bspline(
        sample_points=sample_points,
        interpolate_points=interpolate_points,
        # u_interp=u_interp,
        n_cp=n_cp,
        k=k,
        knot_vector=knot_vector,
    )

    # --- ノットベクトルの確認出力 ---
    print('knot_vector = ', '['+', '.join([f'{val:.3f}' for val in bspline.t])+']')

    # --- フィッティング結果を使ってBスプラインオブジェクトを再構築 ---
    bspline = MyBSpline(knot_vector, control_points, k)

    # --- 除去されたキンクノットを再挿入（必要に応じて） ---
    for knot in knot_to_add:
        counter = Counter(knot_to_add)
        bspline = bspline.insert_knot(knot, 1)

    # --- 最終的なノットベクトルの確認出力 ---
    print('knot_vector = ', '['+', '.join([f'{val:.3f}' for val in bspline.t])+']')

    # --- Bspline曲線の可視化 ---
    bspline.plot_bspline(ref_points=None)
    bspline.plot_basis_function()
    bspline.plot_ddt()
    bspline.plot_ddx(axis=1)

    plt.show()