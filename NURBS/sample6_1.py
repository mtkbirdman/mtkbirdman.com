import numpy as np
import matplotlib.pyplot as plt

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

    # --- Bスプラインオブジェクトの初期化（空の状態） ---
    bspline = MyBSpline()

    # --- 3次元空間上の形状定義（root → mid → tip） ---
    root = np.array([0, 0, 0])         # 始点
    mid = np.array([0, 0.3, 0.03])     # 中間点（キンク位置）
    tip = np.array([0.7, 1.0, 0.10])   # 終点

    # --- サンプル点の生成（root→mid→tipを直線で接続） ---
    sample_points = []
    # root→mid間を線形補間
    for val in np.linspace(0, 1, int(0.3 * n_sample) + 1)[:-1]:
        sample_points.append(root * (1 - val) + mid * val)
    # mid→tip間を線形補間
    for val in np.linspace(0, 1, int(0.7 * n_sample) + 1)[:-1]:
        sample_points.append(mid * (1 - val) + tip * val)
    # tipを追加
    sample_points.append(tip)

    # --- 補間点の定義（曲線が厳密に通過すべき点） ---
    interpolate_points = np.column_stack([root, mid, tip]).T
    u_interp = np.array([0, 0.3, 1.0])  # 各補間点に対応するパラメータ値

    # --- ノットベクトルの生成（キンク位置を反映） ---
    knot_vector = bspline._calc_knot_vector(n_cp, k, kink_list=kink_list)
    print('knot_vector = ', '['+', '.join([f'{val:.3f}' for val in knot_vector])+']')

    # --- Bスプラインのフィッティング（補間＋近似） ---
    control_points = bspline.fit_Bspline(
        sample_points=sample_points,
        interpolate_points=interpolate_points,
        # u_interp=u_interp,
        n_cp=n_cp,
        k=k,
        knot_vector=knot_vector,
    )

    # --- フィッティング結果を使ってBスプラインオブジェクトを再構築 ---
    bspline = MyBSpline(knot_vector, control_points, k)

    # --- ノットベクトルの確認出力 ---
    print('knot_vector = ', '['+', '.join([f'{val:.3f}' for val in bspline.t])+']')

    # --- Bspline曲線の可視化 ---
    bspline.plot_bspline(ref_points=None)
    bspline.plot_basis_function()
    bspline.plot_ddt()
    bspline.plot_ddx(axis=1)

    plt.show()