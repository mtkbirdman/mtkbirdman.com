import numpy as np
import matplotlib.pyplot as plt

from mybspline import MyBSpline


if __name__=='__main__':


    # --- 空のBスプラインオブジェクトを初期化 ---
    # この段階ではノットベクトルや制御点は未設定
    bspline = MyBSpline()

    # --- フィッティングに使用するサンプル点の定義 ---
    n_cp = 10  # 制御点数（近似の自由度を決定）
    k = 3      # Bスプラインの次数（3次＝キュービック）

    # x座標：0 〜 √2π の範囲を100分割した等間隔の点
    x = np.linspace(0, np.sqrt(2) * np.pi, 100)

    # y座標：xに対応するsin関数の値（滑らかな波形）
    y = np.sin(x)

    # サンプル点を2次元座標としてまとめる（shape: [100, 2]）
    sample_points = np.column_stack((x, y))

    # --- 補間点の定義（曲線が厳密に通過すべき点） ---
    # 少数の代表点を選び、補間条件として与える
    x = np.linspace(0, np.sqrt(2) * np.pi, 4)
    y = np.sin(x)
    interpolate_points = np.column_stack((x, y))

    # --- Bスプラインのフィッティング（補間＋近似＋周期条件） ---
    # sample_points → 近似対象の点群
    # interpolate_points → 曲線が厳密に通過すべき点
    # periodic=True → 曲線の始点と終点で1階・2階微分が一致する周期条件を課す

    control_points = bspline.fit_Bspline(
        sample_points=sample_points,
        n_cp=n_cp,
        k=k,
        interpolate_points=interpolate_points,
        periodic=True,
    )

    # --- ノットベクトルの生成（open uniform knot） ---
    # ノットベクトルの長さは m = n_cp + k + 1
    # 端点に k 個の重複ノットを配置し、内部は均等分割
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