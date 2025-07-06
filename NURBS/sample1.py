import numpy as np
import matplotlib.pyplot as plt

from mybspline import MyBSpline


if __name__=='__main__':

    # --- 制御点数と次数の設定 ---
    n_cp = 8         # 制御点の数（Bスプラインの自由度を決定）
    k = 3            # Bスプラインの次数（3次＝キュービック）

    # --- 制御点の座標を生成 ---
    # x座標：0 〜 √2π の範囲を等間隔に分割
    x = np.linspace(0, np.sqrt(2) * np.pi, n_cp)

    # y座標：xに対応するsin関数の値（滑らかな波形）
    y = np.sin(x)

    # 制御点を2次元座標としてまとめる（shape: [n_cp, 2]）
    control_points = np.column_stack((x, y))

    # --- ノットベクトルの設定 ---
    # 初期化
    knot_vector = None

    # ノットベクトルの例①：0.500が2回重複（t=0.500で2階微分が不連続）
    knot_vector = [0.000, 0.000, 0.000, 0.000, 0.250, 0.500, 0.500, 0.750, 1.000, 1.000, 1.000, 1.000]

    # ノットベクトルの例②：0.333が3回重複（t=0.333で1階部分が不連続：キンク）
    # knot_vector = [0.000, 0.000, 0.000, 0.000, 0.333, 0.333, 0.333, 0.667, 1.000, 1.000, 1.000, 1.000]

    # ノットベクトルの例③：0.500が4回重複（やりすぎ）
    # knot_vector = [0.000, 0.000, 0.000, 0.000, 0.500, 0.500, 0.500, 0.500, 1.000, 1.000, 1.000, 1.000]

    # --- Bスプラインオブジェクトの生成 ---
    # 指定したノットベクトル・制御点・次数を用いて MyBSpline を初期化
    bspline = MyBSpline(knot_vector, control_points, k)

    # --- ノットベクトルの確認出力 ---
    # bspline.t に格納されたノットベクトルを取得
    knot_vector = bspline.t

    # 制御点数とノットベクトル長を表示
    print(f'n={n_cp}')  # 制御点数
    print(f'm={len(knot_vector)}', '[' + ', '.join([f'{v:.3f}' for v in knot_vector]) + ']')  # ノットベクトルの中身

    bspline.plot_bspline(ref_points=None)
    bspline.plot_basis_function()
    bspline.plot_ddt()
    bspline.plot_ddx()

    plt.show()