import numpy as np
import matplotlib.pyplot as plt

from mybspline import MyBSpline


if __name__=='__main__':
    # --- 空のBスプラインオブジェクトを初期化 ---
    # この段階ではノットベクトルや制御点は未設定
    bspline = MyBSpline()

    # --- 補間点の定義 ---
    n_cp = 6  # 制御点数（補間点数と一致させることで完全補間を実現）
    k = 3            # Bスプラインの次数（3次＝キュービック）

    # x座標：0 〜 √2π の範囲を等間隔に分割
    x = np.linspace(0, np.sqrt(2) * np.pi, n_cp)

    # y座標：xに対応するsin関数の値（滑らかな波形）
    y = np.sin(x)

    # 補間点を2次元座標としてまとめる（shape: [n_cp, 2]）
    interpolate_points = np.column_stack((x, y))

    # --- Bスプラインのフィッティング（補間条件付き） ---
    # sample_points=None → 補間点のみを使ってフィッティング
    # u_interp=None → パラメータ値は自動的に均等割り当て
    # n_cp → 制御点数（補間点数と一致）
    # degree=k → 指定された次数（例：3次＝キュービック）

    control_points = bspline.fit_Bspline(
        sample_points=None,
        interpolate_points=interpolate_points,
        # u_interp=None,
        n_cp=n_cp,
        k=k,
    )

    # --- ノットベクトルの再計算 ---
    # 補間点に対応するノットベクトルを生成（open uniform knot）
    knot_vector = bspline._calc_knot_vector(n_cp, k)

    # --- フィッティング結果を使って新たにBスプラインオブジェクトを生成 ---
    # ノットベクトルと制御点を明示的に指定して初期化
    bspline = MyBSpline(knot_vector, control_points, k)

    # --- Bスプライン曲線の可視化 ---
    bspline.plot_bspline(ref_points=None)
    bspline.plot_basis_function()
    bspline.plot_ddt()
    bspline.plot_ddx()

    plt.show()