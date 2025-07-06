import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline
import matplotlib.cm as cm
import warnings
import math

from collections import Counter


class MyBSpline(BSpline):


    def __init__(self, knot_vector=None, coeff=None, degree=3):
        """
        MyBSpline クラスの初期化メソッド

        Parameters:
            knot_vector (array-like, optional): ノットベクトル。指定しない場合は自動生成される。
            coeff (array-like, optional): Bスプラインの係数（制御点）。指定しない場合は等間隔で初期化される。
            degree (int): Bスプラインの次数（デフォルトは3次＝キュービック）
        """

        # 係数が指定されていない場合、次数+1個の等間隔の値で初期化
        # これは、最低限の制御点数を確保するための処理（ベジェ曲線）
        if coeff is None:
            coeff = np.linspace(0, 1, degree + 1)

        # ノットベクトルが指定されていない場合、自動生成する
        # ノットベクトルの長さは「制御点数 + 次数 + 1」で決まる
        if knot_vector is None:
            m = len(coeff) + degree + 1  # ノットベクトルの長さ
            # 制御点数と次数をもとに、uniform knot を生成
            knot_vector = self._calc_knot_vector(n_cp=len(coeff), k=degree)

        # 親クラス BSpline の初期化を呼び出す
        # ここで、ノットベクトル・係数・次数を設定する
        super().__init__(knot_vector, coeff, degree)


    def basis_matrix(self, u_values, nu=0):
        """
        指定されたパラメータ値におけるBスプライン基底関数（またはその微分）の値を行列として返す関数

        Parameters:
            u_values (array-like): 評価対象となるパラメータ値の配列（例：np.linspace(0, 1, 100)）
            nu (int): 微分階数（0なら通常の基底関数、1以上ならその階数の微分）

        Returns:
            basis_matrix (ndarray): shape=(len(u_values), n_basis) の基底関数行列
                                    各行は u_values の各点における基底関数の値（または微分値）を表す
        """

        # パラメータ値を NumPy 配列に変換（リストなども受け入れ可能にするため）
        u_values = np.asarray(u_values)

        # 基底関数の数（制御点の数）を計算
        # BSpline の定義より、n_basis = len(t) - k - 1
        # ここで t はノットベクトル、k は次数
        n_basis = len(self.t) - (self.k + 1)

        # 出力行列を初期化（ゼロで埋める）
        # 行数はパラメータ値の数、列数は基底関数の数
        basis_matrix = np.zeros((len(u_values), n_basis))

        # 各基底関数についてループ
        # Bスプラインの基底関数は、1つの係数だけが1で他が0のときに得られる
        for i in range(n_basis):
            # i番目の基底関数を抽出するための係数ベクトルを作成
            # 例：i=2 の場合 → [0, 0, 1, 0, 0, ...]
            coeff = np.zeros(n_basis)
            coeff[i] = 1

            # BSplineオブジェクトを作成（i番目の基底関数）
            # self.t: ノットベクトル, coeff: 係数, self.k: 次数
            basis = BSpline(self.t, coeff, self.k)

            # 必要に応じて微分を計算（nu=0ならそのまま、nu>0なら微分）
            basis = basis.derivative(nu=nu)

            # u_values における基底関数（またはその微分）の値を評価し、行列に格納
            # basis(u_values) は shape=(len(u_values),) のベクトル
            basis_matrix[:, i] = basis(u_values)

        # 最終的な基底関数行列を返す
        return basis_matrix


    def dim(self):
        """
        Bスプラインの制御点が属する空間の次元数を返す関数

        Returns:
            int: 制御点の次元（例：2次元なら2、3次元なら3）
        """

        # self.c は Bスプラインの制御点（係数）を格納する配列
        # 形状は (制御点数, 空間次元) となっている
        # 例：2次元空間なら [[x0, y0], [x1, y1], ...]
        #     3次元空間なら [[x0, y0, z0], [x1, y1, z1], ...]

        # shape[1] を取得することで、各制御点が持つ成分数（＝空間次元）を取得できる
        return self.c.shape[1]


    def ddx(self, r, u_values, axis=0):
        """
        指定された軸（通常はx軸）に関する r 階微分を、他の成分に対して計算する関数

        Parameters:
            r (int): 微分階数（1階、2階など）
            u_values (array-like): 微分を評価するパラメータ値（例：np.linspace(0, 1, 100)）
            axis (int): 微分対象の独立変数の軸（通常は x=0）

        Returns:
            dict: 各成分の r 階微分結果を格納した辞書（key: 軸番号, value: 微分値の配列）
                  ただし、axis 自身は除外される（他の成分のみ）
        """

        # 0階から r階までの微分をすべて計算してリストに格納
        # derivatives[i] は i階微分の結果（shape=(len(u_values), dim)）
        derivatives = [self.derivative(nu=i)(u_values) for i in range(r + 1)]

        # 微分対象軸（例：x軸）の1階微分値を取得（連鎖律の分母に使う）
        dx = derivatives[1][:, axis]

        # 結果を格納する辞書（key: 軸番号, value: 微分値の配列）
        results = {}

        # 全成分に対してループ（ただし微分対象軸は除外）
        for target_axis in range(self.dim()):
            if target_axis == axis:
                continue  # 微分対象軸はスキップ

            # --- 1階微分の場合（連鎖律の基本形） ---
            if r == 1:
                # dy/dx = (dy/dt) / (dx/dt)
                results[target_axis] = derivatives[1][:, target_axis] / dx
                continue

            # --- 高階微分の場合（連鎖律の再帰的適用） ---
            # まず、r階微分の対象成分を取得
            dk_target = derivatives[r][:, target_axis]

            # 最初の項：dy^{(r)}/dt^{(r)} × (dx/dt)^{r-1}
            terms = dk_target * dx**(r - 1)

            # 連鎖律の補正項を計算（j=1〜r-1）
            # これは多変数関数の高階微分における一般化された連鎖律に基づく
            for j in range(1, r):
                coeff = math.comb(r, j)  # 二項係数（rCj）を計算
                # 補正項：rCj × dy^{(j)}/dt^{(j)} × dx^{(r-j)}/dt^{(r-j)} × (dx/dt)^{r-j-1}
                terms -= coeff * derivatives[j][:, target_axis] * derivatives[r - j][:, axis] * dx**(r - j - 1)

            # 最終的な微分値を dx^{2r - 1} で割って正規化
            results[target_axis] = terms / dx**(2 * r - 1)

        # 各成分の微分結果を辞書で返す
        return results


    def fit_Bspline(
        self,
        sample_points=None,       # フィッティング対象のサンプル点（近似用）
        interpolate_points=None,  # 補間すべき点（厳密に通過させたい点）
        n_cp=None,                # 制御点の数（指定しない場合は自動決定）
        k=3,                 # Bスプラインの次数（デフォルトは3次＝キュービック）
        periodic=False,           # 周期条件を課すかどうか（Trueなら周期的な曲線）
        knot_vector=None,         # ノットベクトル（指定しない場合は自動生成）
        constraints=[],
        ):

        if interpolate_points is None:
            n_cp_min = k + 1
            if n_cp is None:
                n_cp = n_cp_min
                print(f'n_cp is set to {n_cp}.')

            if sample_points is None:
                raise ValueError("At least one of sample_points or interpolate_points must be provided.")
            else:
                if n_cp < n_cp_min:
                    print(f'''
                        n_cp must be >= k + 1 = {k + 1}. 
                        Modified n_cp = {n_cp} -> {k+1}.
                    ''')
                    n_cp = n_cp_min
        else:
            n_cp_min = len(interpolate_points) + 2 * periodic + len(constraints) 
            if n_cp is None:
                n_cp = n_cp_min
                print(f'n_cp is set to {n_cp}.')
            
            if sample_points is None:
                if n_cp != n_cp_min:
                    print(f'''
Warning: n_cp must be = len(interpolate_points) + 2 * periodic + len(constraints) = {n_cp_min}.
Modified to n_cp = {n_cp} -> {n_cp_min}.
                    ''')
                elif len(interpolate_points) < k + 1:
                    raise ValueError(f'Warning: len(interpolate_points) must be >= k + 1 = {k + 1}. Current = {len(interpolate_points)}')
            else:
                if n_cp < n_cp_min:
                    print(f'''
Warning: n_cp must be >= len(interpolate_points) + 2 * periodic + len(constraints) = {n_cp_min}.
Modified to n_cp = {n_cp} -> {n_cp_min}.
                    ''')
                    n_cp = n_cp_min
        
        # print('n_cp', n_cp)
        
        if sample_points is not None:
            sample_points = np.asarray(sample_points)
            n_samples, n_dim = sample_points.shape
            u_sample = np.linspace(0, 1, n_samples)
        
        if interpolate_points is not None:
            _, n_dim = interpolate_points.shape
            if sample_points is not None:
                u_interp = self._get_u_interp(sample_points, u_sample, interpolate_points)
                # print('u_interp', u_interp)
            else:
                interpolate_points = np.asarray(interpolate_points)
                u_interp = np.linspace(0, 1, len(interpolate_points))

        if knot_vector is None:
            # print(n_cp)
            knot_vector = self._calc_knot_vector(n_cp=n_cp, k=k)

        self.t = knot_vector
        self.k = k

        if sample_points is not None:

            # --- 基底関数行列 A を構築（サンプル点に対する評価） ---
            # 各サンプル点におけるBスプライン基底関数の値を行列形式で取得
            # A.shape = (n_samples, n_control_points)
            A = self.basis_matrix(u_sample)

            # --- 目的関数の右辺（サンプル点） ---
            # 最小二乗法の目的：Ax ≈ b を満たすような x（制御点）を求める
            b = sample_points

            # --- 補間条件または周期条件がある場合の制約処理 ---
            constraints_rows = []  # 制約行列 C の各行を格納するリスト
            rhs_rows = []          # 制約の右辺ベクトル d の各行を格納するリスト

            # --- 補間点の制約（指定された点を厳密に通過） ---
            if interpolate_points is not None:

                # 各補間点に対して、基底関数の値を1行ずつ構築
                for idx, u in enumerate(u_interp):
                    row = np.zeros((n_cp,))  # 制約行の初期化（長さ = 制御点数）
                    for i in range(n_cp):
                        coeff = np.zeros(n_cp)
                        coeff[i] = 1  # i番目の基底関数だけ1にする
                        basis = BSpline(knot_vector, coeff, k)  # 単一基底関数を構築
                        row[i] = basis(u)  # uにおける基底関数の値を評価
                    constraints_rows.append(row)  # 制約行を追加
                    rhs_rows.append(interpolate_points[idx])  # 対応する補間点を右辺に追加

            # --- 周期条件の制約（1階・2階微分が一致） ---
            if periodic:
                for r in [1, 2]:  # 1階・2階微分に対して制約を課す
                    # 各制御点に対して単位ベクトルを使って微分行列を構築
                    deriv = BSpline(knot_vector, np.eye(n_cp), k).derivative(nu=r)
                    row_start = deriv(0.0)  # 始点での微分値
                    row_end = deriv(1.0)    # 終点での微分値
                    constraints_rows.append(row_start - row_end)  # 差分が0になるように制約
                    rhs_rows.append(np.zeros(n_dim))  # 右辺はゼロベクトル（等しいことを意味）

            if constraints:
                # 各制約対象点に対応するパラメータを求める
                constraint_points = np.array([pt for _, pt in constraints])
                u_constraints = self._get_u_interp(sample_points, u_sample, constraint_points)

                for (r, pt), u in zip(constraints, u_constraints):
                    # r階微分の基底関数行列の1行を構築
                    # print(f'\nr = {r}, u = {u}')
                    deriv_basis = BSpline(knot_vector, np.eye(n_cp), k).derivative(nu=r)
                    row = deriv_basis(u)
                    constraints_rows.append(row)
                    rhs_rows.append(np.zeros(n_dim))  # 微分=0の制約
                    # print(f'row = ', row)
                    # print(f'np.zeros(n_dim) = ', np.zeros(n_dim))
                    
            if len(constraints_rows) > 0:
                # --- 制約行列 C と右辺 d の構築 ---
                C = np.vstack(constraints_rows)  # 制約行列（shape: [n_constraints, n_control_points]）
                d = np.vstack(rhs_rows)          # 制約右辺ベクトル（shape: [n_constraints, n_dim]）

                # --- KKTシステムの構築（制約付き最小二乗法） ---
                # 通常の最小二乗法の正規方程式：AᵀA x = Aᵀb
                ATA = A.T @ A
                ATb = A.T @ b
                CT = C.T

                # KKT行列（拡張連立方程式）を構築
                # [ AᵀA   Cᵀ ] [x]   = [Aᵀb]
                # [  C     0 ] [λ]     [ d ]
                KKT_matrix = np.block([
                    [ATA, CT],
                    [C, np.zeros((C.shape[0], C.shape[0]))]
                ])
                rhs_full = np.vstack([ATb, d])  # 右辺ベクトルを連結

                # --- KKTシステムを解いて制御点を取得 ---
                # 解ベクトル [x, λ] を求める（x: 制御点, λ: ラグランジュ乗数）
                solution = np.linalg.solve(KKT_matrix, rhs_full)
                self.c = solution[:n_cp]  # 制御点のみを抽出して保存
            else:
                # --- 通常の最小二乗法（制約なし） ---
                # Ax ≈ b を最小二乗で解く（np.linalg.lstsq は安定性が高い）
                self.c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

        else:
            
            A = self.basis_matrix(u_interp)
            b = interpolate_points
            
            constraints_rows = []  # 制約行列 C の各行を格納するリスト
            rhs_rows = []          # 制約の右辺ベクトル d の各行を格納するリスト

            # --- 周期条件の制約（1階・2階微分が一致） ---
            if periodic:
                for r in [1, 2]:  # 1階・2階微分に対して制約を課す
                    # 各制御点に対して単位ベクトルを使って微分行列を構築
                    deriv = BSpline(knot_vector, np.eye(n_cp), k).derivative(nu=r)
                    row_start = deriv(0.0)  # 始点での微分値
                    row_end = deriv(1.0)    # 終点での微分値
                    constraints_rows.append(row_start - row_end)  # 差分が0になるように制約
                    rhs_rows.append(np.zeros(n_dim))  # 右辺はゼロベクトル（等しいことを意味）

            if constraints:
                # 各制約対象点に対応するパラメータを求める
                constraint_points = np.array([pt for _, pt in constraints])
                u_constraints = self._get_u_interp(interpolate_points, u_interp, constraint_points)
                # print('u_constraints', u_constraints)

                for (r, pt), u in zip(constraints, u_constraints):
                    # r階微分の基底関数行列の1行を構築
                    deriv_basis = BSpline(knot_vector, np.eye(n_cp), k).derivative(nu=r)
                    row = deriv_basis(u)
                    constraints_rows.append(row)
                    rhs_rows.append(np.zeros(n_dim))  # 微分=0の制約

            if len(constraints_rows) > 0:

                # print('n_cp', n_cp)
                # print('A', A.shape)
                # print('b', b.shape)

                # --- 制約行列 C と右辺 d の構築 ---
                C = np.vstack(constraints_rows)  # 制約行列（shape: [n_constraints, n_control_points]）
                d = np.vstack(rhs_rows)          # 制約右辺ベクトル（shape: [n_constraints, n_dim]）
                
                # print('C', C.shape)
                # print('d', d.shape)

                A = np.vstack([A, C])
                b = np.vstack([b, d])  # 右辺ベクトルを連結

                # print('A', A.shape)
                # print('b', b.shape)

            solution = np.linalg.solve(A, b)
            self.c = solution[:n_cp]  # 制御点のみを抽出して保存

        # --- フィッティングされた制御点を返す ---
        return self.c


    def _get_u_interp(self, sample_points, u_sample, interpolate_points):
        """
        各補間点に対して、最近傍2点のサンプル点から距離に基づき線形内挿で u_interp を計算する関数。

        Parameters:
            sample_points (ndarray): shape=(N, D) のサンプル点（N個のD次元ベクトル）
            u_sample (ndarray): shape=(N,) の各サンプル点に対応するパラメータ値
            interpolate_points (ndarray): shape=(M, D) の補間点

        Returns:
            u_interp (ndarray): shape=(M,) の各補間点に対応するパラメータ値（線形内挿された u 値）
        """
        sample_points = np.asarray(sample_points)
        u_sample = np.asarray(u_sample)
        interpolate_points = np.asarray(interpolate_points)

        # print('u_sample', u_sample)

        u_interp = np.zeros(len(interpolate_points))

        for i, p in enumerate(interpolate_points):
            # 各補間点に対する距離を計算
            distances = np.linalg.norm(sample_points - p, axis=1)
            # 最も近い2点のインデックスを取得
            nearest_indices = np.argsort(distances)[:2]
            i1, i2 = nearest_indices

            # サンプル点・パラメータ値・距離を取得
            p1, p2 = sample_points[i1], sample_points[i2]
            u1, u2 = u_sample[i1], u_sample[i2]
            d1, d2 = distances[i1], distances[i2]
            # print('p1, p2 = ', p1, p2)
            # print('u1, u2 = ', u1, u2)
            # print('d1, d2 = ', d1, d2)
        
            if d1 + d2 == 0:
                # 補間点がサンプル点と一致している場合（距離ゼロ）
                u_interp[i] = u1  # または u2（同じ）
            else:
                # 線形補間: 重み付き平均
                w1 = d2 / (d1 + d2)
                w2 = d1 / (d1 + d2)
                # print('w1, w2 = ', w1, w2)
                u_interp[i] = w1 * u1 + w2 * u2
        

        return u_interp


    def _remove_kink_from_knot(self, knot_vector, k):
        """
        ノットベクトルから不要な重複（キンク）を除去する関数。
        ただし、先頭と末尾の k 個のノットは保持する（open uniform knot の境界条件を維持するため）。

        Parameters:
            knot_vector (list or array): 元のノットベクトル
            k (int): Bスプラインの次数（境界部のノット数を決定する）

        Returns:
            cleaned_knot_vector (list): 重複を除去したノットベクトル
            duplicates (list): 除去された重複ノットの一覧
        """

        # --- 境界部のノットを保持（open uniform knot の条件） ---
        # 先頭 k 個と末尾 k 個のノットはそのまま残す
        head = knot_vector[:k]      # 例： [0, 0, 0]（k=3 の場合）
        tail = knot_vector[-k:]     # 例： [1, 1, 1]

        # --- 中央部分のノットを抽出（重複除去対象） ---
        middle = knot_vector[k:-k]  # 境界部を除いたノット列

        # --- 重複除去処理（順序を保持しながら重複を除去） ---
        seen = set()                # すでに出現したノット値を記録する集合
        cleaned_middle = []         # 重複を除去した中央ノット列
        duplicates = []             # 除去された重複ノットの記録

        for val in middle:
            if val not in seen:
                seen.add(val)           # 初出の値は記録して保持
                cleaned_middle.append(val)
            else:
                duplicates.append(val)  # 重複していた値は記録のみ（除去）

        # --- 最終的なノットベクトルを構築 ---
        # 境界部 + 重複除去済み中央部 + 境界部（末尾）
        cleaned_knot_vector = head + cleaned_middle + tail

        # --- 結果を返す ---
        return cleaned_knot_vector, duplicates


    def _calc_knot_vector(self, n_cp, k, kink_list=[]):
        """
        Bスプラインのノットベクトルを生成する関数。
        任意の位置に指定した重複度でノット（キンク）を挿入可能。

        Parameters:
            n_cp (int): 制御点の数
            k (int): Bスプラインの次数
            kink_list (list of tuple): キンク位置とその重複度のリスト（例：[(0.3, 2), (0.7, 3)]）

        Returns:
            knot_vector (list of float): 生成されたノットベクトル
        """

        # キンク位置のリストと、その重複度の合計を取り出す
        kink_positions = [t for t, _ in kink_list]
        kink_multiplicities = [m for _, m in kink_list]
        kink_total_mult = sum(kink_multiplicities)

        # ノットベクトルの全長（= m）を計算：n_cp + k + 1
        m = n_cp + k + 1

        # ノットの始点と終点（0と1）とキンクを含めた全ノット区間
        knot_value = [0] + kink_positions + [1]

        # キンク区間以外に割り当てる自由なノット数を計算
        # 両端の固定ノット 2*(k+1)、キンク重複分のノットは除外
        m_add = m - 2 * (k + 1) - kink_total_mult

        # 区間に比例したノット数の割り当て
        interval_lengths = np.diff(np.asarray(knot_value))
        m_add_list = [int(val) for val in interval_lengths * m_add]

        # ノット数の補正（誤差分を調整）
        total = sum(m_add_list)
        if total < m_add:
            sorted_idx = np.argsort(-interval_lengths)
            for i in range(m_add - total):
                m_add_list[sorted_idx[i % len(sorted_idx)]] += 1

        # ノットベクトルの構築
        knot_vector = [0.0] * (k + 1)  # 始端の重複ノット（k+1個）

        for i in range(len(knot_value) - 1):
            t0, t1 = knot_value[i], knot_value[i + 1]
            seg_count = m_add_list[i]

            # 開区間にノットを均等配置（両端除く）
            if seg_count > 0:
                interior_knots = np.linspace(t0, t1, seg_count + 2)[1:-1]
                knot_vector.extend(np.round(interior_knots, 6).tolist())

            # キンク位置（t1）に重複ノットを追加（指定されている場合）
            for t, mult in kink_list:
                if np.isclose(t, t1):
                    knot_vector.extend([t] * mult)
                    break

        knot_vector.extend([1.0] * (k + 1))  # 終端の重複ノット

        return knot_vector


    def plot_basis_function(self, show=False):
        """
        Bスプラインの基底関数とその微分をプロットする関数。
        各微分階数（0階〜最大階数）に対して、基底関数の形状を可視化する。
        """

        # ノットベクトルと次数を取得
        knot_vector = self.t
        degree = self.k

        # 微分可能な最大階数を取得（通常は次数と同じだが、ノットの重複などで制限される場合あり）
        r_max = self._max_derivative_order()

        # プロット領域を準備（微分階数ごとにサブプロットを並べる）
        fig, axes = plt.subplots(figsize=(3 * (r_max + 1), 3), ncols=r_max + 1)

        # 各微分階数（nu）に対してループ
        for nu in range(r_max + 1):
            ax = axes.flatten()[nu]  # 対応するサブプロットを取得

            # 各ノット区間 [t0, t1] に対してループ（次数に基づいて有効な区間を選ぶ）
            for j, (t0, t1) in enumerate(zip(knot_vector[degree:-(degree + 1)], knot_vector[(degree + 1):-degree])):
                # パラメータ値を区間内で細かく分割（評価点を生成）
                u_values = np.linspace(t0, t1, 101)[:-1]  # 最後の点は重複を避けるため除外

                # 指定された微分階数で基底関数行列を計算（shape: [len(u_values), n_basis]）
                basis_matrix = self.basis_matrix(u_values, nu=nu)

                # 各基底関数（列）に対してループしてプロット
                for i, row in enumerate(basis_matrix.T):
                    color = cm.tab10(i % 10)  # 色をインデックスに基づいて選択
                    # 最後の区間の終点でラベルを表示（凡例用）
                    label = f'P{i}' if t1 == knot_vector[-(degree + 1)] else None
                    ax.plot(u_values, row, label=label, color=color)

            # ノット位置に縦線を描画（基底関数の切り替わり位置を示す）
            for knot in knot_vector:
                ax.axvline(knot, linestyle=':', color='k')

            # 軸ラベルとグリッドを設定
            ax.set_xlabel('parameter t')
            ax.set_ylabel(f'basis function')  # 微分階数を表示
            ax.set_title(f'B(t) r={nu}')  # 微分階数を表示
            ax.grid(True)

        # 凡例を右上に表示（最後のサブプロットにまとめて表示）
        ax.legend(bbox_to_anchor=(1, 1), loc='upper left')

        # レイアウト調整して表示
        plt.tight_layout()
        if show:
            plt.show()


    def plot_bspline(self, ref_points=None, show=False):
        """
        Bスプライン曲線と制御点、ノット位置を可視化する関数。
        必要に応じて、参照点（補間点やサンプル点）も重ねて表示可能。

        Parameters:
            ref_points (array-like, optional): 参照点（補間点やサンプル点など）をプロットに追加する
        """

        # --- Bスプラインの基本情報を取得 ---
        knot_vector = self.t         # ノットベクトル
        control_points = self.c      # 制御点（係数）
        k = self.k                   # Bスプラインの次数

        # --- プロット領域の準備 ---
        fig, ax = plt.subplots(figsize=(4, 3))

        # --- 制御点のプロット（点と線で接続） ---
        ax.plot(control_points[:, 0], control_points[:, 1],
                'o:', markerfacecolor='none', color='k', label='control point')

        # --- Bスプライン曲線の描画 ---
        # 有効なノット区間ごとに曲線を描画（次数に基づいて区間を選定）
        for i, (t0, t1) in enumerate(zip(knot_vector[k:-(k+1)], knot_vector[(k+1):-k])):
            color = cm.tab10(i)  # 色をインデックスに基づいて選択

            # パラメータ値を区間内で細かく分割
            u = np.linspace(t0, t1, 999)

            # Bスプライン曲線を評価（各 u に対する座標を取得）
            curve = self(u)  # __call__ メソッドで曲線評価

            # 曲線を描画
            ax.plot(curve[:, 0], curve[:, 1], '-', color=color)

            # --- ノット位置のラベル表示 ---
            for t in [t0, t1]:
                point = self(t)  # ノット位置における曲線座標
                ax.plot(point[0], point[1], 'x', color=color)  # ノット位置をマーク
                ax.text(point[0] + 0.1, point[1], f't={t:.3f}', fontsize=8,
                        color='k', ha='left', va='center')  # ノット値をラベル表示

        # --- 参照点の描画（補間点やサンプル点など） ---
        if ref_points is not None:
            ref_points = np.asarray(ref_points)
            ax.plot(ref_points[:, 0], ref_points[:, 1], '.', color='r', alpha=0.5)

        # --- 凡例と軸ラベルの設定 ---
        ax.legend()
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid(True)

        # --- プロットを表示 ---
        if show:
            plt.show()


    def plot_ddt(self, show=False):
        """
        各成分の微分（1階〜k階）をプロットする関数。
        Bスプライン曲線の各空間成分（x, y, zなど）に対して、時間（パラメータt）に関する微分を可視化する。
        """

        # --- Bスプラインの基本情報を取得 ---
        knot_vector = self.t         # ノットベクトル
        k = self.k                   # Bスプラインの次数
        control_points = self.c      # 制御点（係数）
        n_dim = self.dim()           # 空間次元（2D, 3Dなど）

        # --- 微分可能な最大階数を取得（通常は次数と同じ） ---
        r_max = self._max_derivative_order()

        # --- プロット領域の準備 ---
        # 各成分 × 各微分階数 のサブプロットを作成
        fig, axes = plt.subplots(figsize=(3 * r_max, 3 * n_dim), ncols=r_max, nrows=n_dim)
        axis_label = ['x', 'y', 'z']  # 軸ラベル（最大3次元まで対応）

        # --- 各空間成分（x, y, zなど）に対してループ ---
        for i in range(n_dim):
            # --- 各微分階数（1階〜最大階数）に対してループ ---
            for r in range(1, r_max + 1):
                # 対応するサブプロットを取得
                ax = axes.flatten()[i * r_max + r - 1]

                # r階微分のBスプライン関数を取得
                bspline_deriv = self.derivative(nu=r)

                # y軸の表示範囲を初期化（後で自動調整）
                vmin, vmax = np.inf, -np.inf

                # --- 有効なノット区間ごとにループ ---
                for j, (t0, t1) in enumerate(zip(knot_vector[k:-(k+1)], knot_vector[(k+1):-k])):
                    # パラメータ値を区間内で細かく分割
                    u_values = np.linspace(t0, t1, 100)[:-1]  # 最後の点は重複を避ける
                    color = cm.tab10(j)  # 色をインデックスに基づいて選択

                    # 各 u に対する微分値を評価（shape: [len(u_values), n_dim]）
                    deriv_vals = np.round(np.array([bspline_deriv(u) for u in u_values]), 6)

                    # i番目の成分（x, y, zなど）をプロット
                    ax.plot(u_values, deriv_vals[:, i], color=color)

                    # y軸の表示範囲を更新
                    vmin = min(vmin, np.min(deriv_vals[:, i]))
                    vmax = max(vmax, np.max(deriv_vals[:, i]))

                # --- ノット位置に縦線を描画（関数の切り替わり点） ---
                for knot in knot_vector:
                    ax.axvline(knot, linestyle=':', color='k')

                # --- 軸ラベルとタイトルの設定 ---
                ax.set_xlabel('parameter t')
                ax.set_ylabel('derivative')
                ax.set_title(rf'$\frac{{d^{r}{axis_label[i]}}}{{dt^{r}}}$')  # 数式形式で表示
                ax.grid(True)

        # --- レイアウト調整して表示 ---
        plt.tight_layout()
        if show:
            plt.show()


    def plot_ddx(self, axis=0, show=False):
        """
        指定された軸（通常はx軸）に関する他成分の微分（連鎖律に基づく）をプロットする関数。
        各成分の1階〜最大階数の微分を、指定軸に対して可視化する。

        Parameters:
            axis (int): 微分対象の独立変数の軸（通常は x=0）
        """

        # --- Bスプラインの基本情報を取得 ---
        knot_vector = self.t         # ノットベクトル
        k = self.k                   # Bスプラインの次数
        n_dim = self.dim()           # 空間次元（2D, 3Dなど）

        # --- 微分可能な最大階数を取得（通常は次数と同じ） ---
        r_max = self._max_derivative_order()

        # --- 微分対象軸以外の成分インデックスを取得 ---
        dims = [x for x in range(n_dim)]
        axis_label = ['x', 'y', 'z']  # 軸ラベル（最大3次元まで対応）
        dims.remove(axis)             # 微分対象軸は除外

        # --- プロット領域の準備 ---
        # 各成分 × 各微分階数 のサブプロットを作成
        fig, axes = plt.subplots(figsize=(3 * r_max, 3 * (n_dim - 1)), ncols=r_max, nrows=n_dim - 1)

        # --- 各成分（微分対象軸以外）に対してループ ---
        for j, dim in enumerate(dims):
            # --- 各微分階数（1階〜最大階数）に対してループ ---
            for r in range(1, r_max + 1):
                # 対応するサブプロットを取得（次元数に応じて1D/2D対応）
                if r_max == 1 and (n_dim - 1) == 1:
                  ax = axes
                elif (n_dim - 1) > 1:
                  ax = axes.flatten()[j * r_max + r - 1]
                else:
                  ax = axes[r - 1]

                # y軸の表示範囲を初期化（後で自動調整）
                vmin, vmax = np.inf, -np.inf

                # --- 有効なノット区間ごとにループ ---
                for i, (t0, t1) in enumerate(zip(knot_vector[k:-(k+1)], knot_vector[(k+1):-k])):
                    color = cm.tab10(i)  # 色をインデックスに基づいて選択

                    # パラメータ値を区間内で細かく分割
                    u_values = np.linspace(t0, t1, 101)[1:-1]

                    # 指定軸に関する r階微分を連鎖律で計算（target_axis=dim）
                    result = np.round(self.ddx(r, u_values, axis=axis)[dim], 6)

                    # 微分結果をプロット
                    ax.plot(u_values, result, label=f'$(d/d{axis_label[axis]})^{r}$', color=color)

                    # y軸の表示範囲を更新
                    vmin = min(vmin, np.min(result))
                    vmax = max(vmax, np.max(result))

                # --- ノット位置に縦線を描画（関数の切り替わり点） ---
                for knot in knot_vector:
                    ax.axvline(knot, linestyle=':', color='k')

                # --- 軸ラベルとタイトルの設定 ---
                ax.set_xlabel('parameter t')
                ax.set_ylabel('derivatives')
                ax.set_title(rf'$\left(\frac{{d{axis_label[dim]}}}{{d{axis_label[axis]}}}\right)^{r}$')
                ax.grid(True)

        # --- レイアウト調整して表示 ---
        plt.tight_layout()
        if show:
            plt.show()


    def _max_derivative_order(self):
        """
        Bスプラインが微分可能な最大階数を返す関数。
        通常は次数（self.k）と一致するが、ノットの重複などにより制限される場合がある。

        Returns:
            result (int): 微分可能な最大階数（0以上の整数）
        """

        # 初期値として、次数（k）を最大微分階数の候補とする
        result = self.k

        # --- 微分可能かどうかを確認するループ ---
        # 高階微分が ValueError を出す場合があるため、例外処理で安全に判定する
        while result > 0:
            try:
                # 指定階数で微分を試みる（成功すればその階数は有効）
                self.derivative(nu=result)
                break  # 微分成功 → ループ終了
            except ValueError:
                # 微分失敗 → 階数を1つ下げて再試行
                result -= 1
        # 最終的に判定された最大微分階数を返す
        return result


if __name__=='__main__':
    print('see ')