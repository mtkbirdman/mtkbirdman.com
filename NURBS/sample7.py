import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid

from mybspline import MyBSpline


def fx(x):
    x = np.asarray(x)
    y = np.where(
        x <= 0.5,
        np.sqrt(x),
        np.sqrt(0.5) + (1/2)*(0.5**(-1/2))*(x - 0.5)
    )
    return y


def fy(y):
    y = np.asarray(y)
    y_split = np.sqrt(0.5)
    a = 1 / (2 * np.sqrt(0.5))  # 傾き

    x = np.where(
        y <= y_split,
        y**2,
        0.5 + (y - y_split) / a
    )
    return x


def equal_arc_length_points(fx, a, b, n=100):
    # 解像度（累積長さの評価用のサンプル数）
    num_samples = 1000

    # xのサンプル点
    x_dense = np.linspace(a, b, num_samples)
    y_dense = fx(x_dense)

    # 曲線の勾配（dy/dx）→ dx, dy による距離要素を計算
    dx = np.gradient(x_dense)
    dy = np.gradient(y_dense)

    if y_dense.ndim == 1:
        ds = np.sqrt(dx**2 + dy**2)
    else:
        dy_norm = np.linalg.norm(dy, axis=0)
        ds = np.sqrt(dx**2 + dy_norm**2)

    # 累積距離（アーク長）
    arc_length = cumulative_trapezoid(ds, x_dense, initial=0)

    # アーク長の最大値をn等分する
    target_lengths = np.linspace(0, arc_length[-1], n + 1)

    # アーク長に対するxの補間関数を作成
    interp_func = interp1d(arc_length, x_dense)

    # 等アーク長に対応するx座標を取得
    x_equal = interp_func(target_lengths)

    return x_equal


if __name__=='__main__':
    np.set_printoptions(linewidth=999)

    # ===== サンプル関数の作成 =====

    n_sample = 999
    x = equal_arc_length_points(fx, a=0.0, b=1.0, n=n_sample)
    y = fx(x)
    sample_points = np.column_stack((x, y))
    print('sample_points.shape = ', sample_points.shape)

    n_cp = 10

    # --- 補間点の定義（曲線が厳密に通過すべき点） ---
    n = 2
    x = np.concatenate([equal_arc_length_points(fx, a=0.0, b=0.5, n=n), equal_arc_length_points(fx, a=0.5, b=1.0, n=n)[1:]])
    y = fx(x)
    interpolate_points = np.column_stack((x, y))
    print('interpolate_points.shape = ', interpolate_points.shape)

    fig, ax = plt.subplots()
    ax.plot(sample_points[:, 0], sample_points[:, 1], 'o-', markerfacecolor='none')
    ax.plot(interpolate_points[:, 0], interpolate_points[:, 1], '*', color='k', markerfacecolor='yellow', markersize=10)


    # ===== k=3のBSpline曲線 =====

    print('\nk=3')
    bspline = MyBSpline()
    k = 3      # Bスプラインの次数（3次＝キュービック）
    constraints = [
        (2, [0.5, fx(0.5)]), (3, [0.5, fx(0.5)]), 
        (2, [1.0, fx(1.0)]), 
        # (3, [1.0, fx(1.0)]),
    ]
    kink_points = [[0.5, fx(0.5)]]
    u_kink = bspline._get_u_interp(sample_points, np.linspace(0, 1, len(sample_points)+1), kink_points)[0]
    kink_list = [(u_kink, 1)]  # 位置と重複度
    knot_vector = bspline._calc_knot_vector(n_cp=n_cp, k=k, kink_list=kink_list)
    print('u_kink', u_kink)
    print('len(knot_vector)', len(knot_vector))
    print('knot_vector = ['+', '.join([f'{val:.3f}' for val in knot_vector])+']')
    control_points = bspline.fit_Bspline(
        sample_points=sample_points,
        n_cp=n_cp,
        k=k,
        interpolate_points=interpolate_points,
        knot_vector=knot_vector,
        constraints=constraints,
    )
    bspline = MyBSpline(knot_vector, control_points, k)
    bspline.plot_bspline(ref_points=interpolate_points)
    bspline.plot_ddt()
    bspline.plot_ddx()
    bspline.plot_basis_function()
    sample3 = bspline(np.linspace(0, 1, n_sample+1))


    # ===== k=4のBSpline曲線 =====

    print('\nk=4')
    bspline = MyBSpline()
    k = 4      # Bスプラインの次数（3次＝キュービック）
    constraints = [
        (2, [0.5, fx(0.5)]), (3, [0.5, fx(0.5)]), 
        (2, [1.0, fx(1.0)]), (3, [1.0, fx(1.0)]),
    ]
    kink_points = [[0.5, fx(0.5)]]
    u_kink = bspline._get_u_interp(sample_points, np.linspace(0, 1, len(sample_points)+1), kink_points)[0]
    kink_list = [(u_kink, 2)]  # 位置と重複度
    knot_vector = bspline._calc_knot_vector(n_cp=n_cp+2, k=k, kink_list=kink_list)
    print('u_kink', u_kink)
    print('len(knot_vector)', len(knot_vector))
    print('knot_vector = ['+', '.join([f'{val:.3f}' for val in knot_vector])+']')
    control_points = bspline.fit_Bspline(
        sample_points=sample_points,
        n_cp=n_cp+2,
        k=k,
        interpolate_points=interpolate_points,
        knot_vector=knot_vector,
        constraints=constraints,
    )
    bspline = MyBSpline(knot_vector, control_points, k)
    bspline.plot_bspline(ref_points=interpolate_points)
    bspline.plot_ddt()
    bspline.plot_ddx()
    bspline.plot_basis_function()
    sample4 = bspline(np.linspace(0, 1, n_sample+1))


    # ===== k=5のBSpline曲線 =====

    print('\nk=5')
    bspline = MyBSpline()
    k = 5      # Bスプラインの次数（3次＝キュービック）
    constraints = [
        (2, [0.5, fx(0.5)]), (3, [0.5, fx(0.5)]), (4, [0.5, fx(0.5)]),
        (2, [1.0, fx(1.0)]), (3, [1.0, fx(1.0)]), (4, [1.0, fx(1.0)]),
    ]
    kink_points = [[0.5, fx(0.5)]]
    u_kink = bspline._get_u_interp(sample_points, np.linspace(0, 1, len(sample_points)+1), kink_points)[0]
    kink_list = [(u_kink, 2)]  # 位置と重複度
    knot_vector = bspline._calc_knot_vector(n_cp=n_cp+4, k=k, kink_list=kink_list)
    print('u_kink', u_kink)
    print('len(knot_vector)', len(knot_vector))
    print('knot_vector = ['+', '.join([f'{val:.3f}' for val in knot_vector])+']')
    control_points = bspline.fit_Bspline(
        sample_points=sample_points,
        n_cp=n_cp+4,
        k=k,
        interpolate_points=interpolate_points,
        knot_vector=knot_vector,
        constraints=constraints,
    )
    bspline = MyBSpline(knot_vector, control_points, k)
    bspline.plot_bspline(ref_points=interpolate_points)
    bspline.plot_ddt()
    bspline.plot_ddx()
    bspline.plot_basis_function()
    sample5 = bspline(np.linspace(0, 1, n_sample+1))

    # ===== 比較プロット1 =====
    fig, ax = plt.subplots()
    ax.plot(sample_points[:, 0], sample_points[:, 1], '-', markerfacecolor='none', label='sample')
    ax.plot(sample3[:, 0], sample3[:, 1], '-', markerfacecolor='none', label='k=3')
    ax.plot(sample4[:, 0], sample4[:, 1], '-', markerfacecolor='none', label='k=4')
    ax.plot(sample5[:, 0], sample5[:, 1], '-', markerfacecolor='none', label='k=5')

    # ===== 比較プロット2 =====
    fig, ax = plt.subplots()

    # 基準ライン（ゼロ誤差ライン）
    ax.plot(sample_points[:, 0], np.zeros_like(sample_points[:, 0]), '-', label='sample')

    # 各補間結果との距離を計算して描画
    for sample, label in zip([sample3, sample4, sample5], ['k=3', 'k=4', 'k=5']):
        interp_y = np.interp(sample_points[:, 0], sample[:, 0], sample[:, 1])
        # 距離 = 垂直方向の差（1次元の場合）
        distances = interp_y - sample_points[:, 1]
        ax.plot(sample_points[:, 0], distances, '-', markerfacecolor='none', label=label)

    ax.set_ylabel('Vertical Distance')  # 距離ラベルに変更
    ax.legend()
    plt.show()
