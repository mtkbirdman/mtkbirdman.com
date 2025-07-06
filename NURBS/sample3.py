import numpy as np
import matplotlib.pyplot as plt

from mybspline import MyBSpline


if __name__=='__main__':


    # --- サンプル点の定義 ---
    # x座標：0 〜 √2π の範囲を100分割した等間隔の点
    x = np.linspace(0, np.sqrt(2) * np.pi, 100)

    # y座標：xに対応するsin関数の値（滑らかな波形）
    y = np.sin(x)

    # サンプル点を2次元座標としてまとめる（shape: [100, 2]）
    sample_points = np.column_stack((x, y))

    # --- 制御点数の設定 ---
    n_cp = 5         # 制御点数（少ないほど近似の自由度が低くなる）
    k = 4            # Bスプラインの次数（3次＝キュービック）

    # --- Bスプラインオブジェクトの初期化（空の状態） ---
    bspline = MyBSpline()

    # --- Bスプラインのフィッティング（近似） ---
    # sample_points → 近似対象の点群
    # interpolate_points=None → 補間条件なし（近似のみ）
    # u_interp=None → パラメータ値は自動的に均等割り当て
    # n_cp → 制御点数（自由度の調整）
    # degree=k → 指定された次数（例：3次＝キュービック）

    control_points = bspline.fit_Bspline(
        sample_points=sample_points,
        n_cp=n_cp,
        k=k,
    )

    # --- ノットベクトルの生成（open uniform knot） ---
    # ノットベクトルの長さは m = n_cp + k + 1
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