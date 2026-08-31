import math
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import BSpline
# from scipy.optimize import minimize_scalar
# from collections import Counter


class MyBSpline(BSpline):

    def __init__(self, t=None, c=None, k=4):
        
        # BSpline の初期化
        if t is not None and c is not None:
            super().__init__(t, c, k)

    def basis_matrix(self, u_values, nu=0, w=None):
        u_values = np.asarray(u_values)

        # t, c, k = self.tck
        # t = self.t
        # c = self.c
        # k = self.k
        # n = len(c)
        # basis_matrix = np.zeros((len(u_values), n))

        # for i in range(n):
        #     c_i = np.zeros(n)
        #     c_i[i] = 1
        #     spline = BSpline(t, c_i, k)
        #     basis_matrix[:, i] = spline(u_values, nu=nu)
        
        t, c, k = self.tck
        c_eye = np.eye(len(c), len(c))
        bspline = BSpline(t, c_eye, k)
        basis_matrix = bspline(u_values, nu=nu)

        # 重み付き処理（Rational B-Spline的な操作）
        if w is not None:
            w = np.asarray(w)
            if len(w) != basis_matrix.shape[1]:
                raise ValueError("The weight vector w must have the same length as the number of basis functions")
            weighted_basis = basis_matrix * w
            denominator = np.sum(weighted_basis, axis=1, keepdims=True)
            basis_matrix = weighted_basis / denominator

        return basis_matrix
        
    def dim(self):
        return self.c.shape[1]  # c: (n_ctrl_points, n_dim)

    def ddx(self, r, u_values, axis=0):
        u_values = np.asarray(u_values)

        # 0階〜r階の導関数の行列を計算（各階ごとに basis_matrix × control_points）
        derivatives = [
            self.basis_matrix(u_values, nu=i) @ self.c
            for i in range(r + 1)
        ]

        # dx/dt（r=1のときのaxis成分）
        dx = derivatives[1][:, axis]

        # ゼロ除算を防ぐために0→NaNに
        safe_dx = np.where(dx == 0, np.nan, dx)

        results = {}

        for target_axis in range(self.dim()):
            if target_axis == axis:
                continue

            # 1階導関数ならシンプルな商の微分
            if r == 1:
                results[target_axis] = derivatives[1][:, target_axis] / safe_dx
                continue

            # r階微分（連鎖律の一般化）
            dk_target = derivatives[r][:, target_axis]
            terms = dk_target * safe_dx ** (r - 1)

            for j in range(1, r):
                c = math.comb(r, j)
                dyj_dtj = derivatives[j][:, target_axis]
                dxrj_dtrj = derivatives[r - j][:, axis]
                terms -= c * dyj_dtj * dxrj_dtrj * safe_dx ** (r - j - 1)

            results[target_axis] = terms / safe_dx ** (2 * r - 1)

        return results

    def fit_Bspline(
        self,
        sample_points=None,
        interpolate_points=None,
        u_sample=None,
        u_interp=None,
        n_cp=None,
        k=4,
        periodic=False,
        t=None,
        constraints=[]
    ):
        if interpolate_points is None:
            n_cp_min = k + 1
            if n_cp is None:
                n_cp = n_cp_min
                print(f'n_cp is set to {n_cp}.')
            if sample_points is None:
                raise ValueError("At least one of sample_points or interpolate_points must be provided.")
            elif n_cp < n_cp_min:
                print(f'n_cp must be >= k + 1 = {k + 1}. Modified n_cp = {n_cp} -> {k + 1}.')
                n_cp = n_cp_min
        else:
            interpolate_points = np.asarray(interpolate_points)
            n_cp_min = len(interpolate_points) + 2 * periodic + len(constraints)
            if n_cp is None:
                n_cp = n_cp_min
                print(f'n_cp is set to {n_cp}.')
            if sample_points is None:
                if n_cp != n_cp_min:
                    print(f'Warning: n_cp must be = len(interpolate_points) + 2 * periodic + len(constraints) = {n_cp_min}. '
                        f'Modified n_cp = {n_cp} -> {n_cp_min}.')
                    n_cp = n_cp_min
                elif len(interpolate_points) < k + 1:
                    raise ValueError(f'Warning: len(interpolate_points) must be >= k + 1 = {k + 1}. '
                                    f'Current = {len(interpolate_points)}')
            else:
                if n_cp < n_cp_min:
                    print(f'Warning: n_cp must be >= len(interpolate_points) + 2 * periodic. '
                        f'Modified n_cp = {n_cp} -> {n_cp_min}.')
                    n_cp = n_cp_min

        # サンプル点とパラメータ
        if sample_points is not None:
            sample_points = np.asarray(sample_points)
            n_samples, n_dim = sample_points.shape
            if u_sample is None:
                u_sample = np.linspace(0, 1, n_samples)

        if interpolate_points is not None:
            interpolate_points = np.asarray(interpolate_points)
            _, n_dim = interpolate_points.shape
            if u_interp is None:
                if sample_points is not None:
                    u_interp = self._get_u_mapping(sample_points, u_sample, interpolate_points)
                else:
                    u_interp = np.linspace(0, 1, len(interpolate_points))

        # ノットベクトル生成
        if t is None:
            t = self._calc_knot_vector(n_cp=n_cp, k=k)
        self.t = np.ascontiguousarray(t)
        self.k = k
        self.c = np.eye(len(t)-(k+1), len(t)-(k+1))

        if sample_points is not None:
            A = self.basis_matrix(u_sample)
            b = sample_points
            constraints_rows = []
            rhs_rows = []

            if interpolate_points is not None:
                for idx, u in enumerate(u_interp):
                    row = self.basis_matrix([u], nu=0)[0]
                    constraints_rows.append(row)
                    rhs_rows.append(interpolate_points[idx])

            if periodic:
                for r in [1, 2]:
                    basis_start = self.basis_matrix([0.0], nu=r)[0]
                    basis_end = self.basis_matrix([1.0], nu=r)[0]
                    row = basis_start - basis_end
                    constraints_rows.append(row)
                    rhs_rows.append(np.zeros(n_dim))

            if constraints:
                constraint_points = np.array([pt for _, pt in constraints])
                u_constraints = self._get_u_mapping(sample_points, u_sample, constraint_points)
                for (r, pt), u in zip(constraints, u_constraints):
                    row = self.basis_matrix([u], nu=r)[0]
                    constraints_rows.append(row)
                    rhs_rows.append(np.zeros(n_dim))

            if len(constraints_rows) > 0:
                C = np.vstack(constraints_rows)
                d = np.vstack(rhs_rows)
                ATA = A.T @ A
                ATb = A.T @ b
                CT = C.T
                KKT_matrix = np.block([
                    [ATA, CT],
                    [C, np.zeros((C.shape[0], C.shape[0]))]
                ])
                rhs_full = np.vstack([ATb, d])
                solution = np.linalg.solve(KKT_matrix, rhs_full)
                self.c = solution[:n_cp]
                return self.c, KKT_matrix, rhs_full
            else:
                self.c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
                return self.c, A, b

        else:
            # 補間点のみ与えられているケース
            A = self.basis_matrix(u_interp)
            b = interpolate_points
            constraints_rows = []
            rhs_rows = []

            if periodic:
                for r in [1, 2]:
                    basis_start = self.basis_matrix([0.0], nu=r)[0]
                    basis_end = self.basis_matrix([1.0], nu=r)[0]
                    row = basis_start - basis_end
                    constraints_rows.append(row)
                    rhs_rows.append(np.zeros(n_dim))

            if constraints:
                constraint_points = np.array([pt for _, pt in constraints])
                u_constraints = self._get_u_mapping(interpolate_points, u_interp, constraint_points)
                for (r, pt), u in zip(constraints, u_constraints):
                    row = self.basis_matrix([u], nu=r)[0]
                    constraints_rows.append(row)
                    rhs_rows.append(np.zeros(n_dim))

            if len(constraints_rows) > 0:
                C = np.vstack(constraints_rows)
                d = np.vstack(rhs_rows)
                A = np.vstack([A, C])
                b = np.vstack([b, d])
            solution = np.linalg.solve(A, b)
            self.c = solution[:n_cp]
            return self.c, A, b
        
    def dt(self, u_param, r):
        basis_row = self.basis_matrix([u_param], nu=r)[0]
        return basis_row @ self.c

    def curvature(self, u_param):
        r1 = self.dt(u_param, 1)
        r2 = self.dt(u_param, 2)
        cross_norm = np.linalg.norm(np.cross(r1, r2))
        r1_norm = np.linalg.norm(r1)
        if r1_norm == 0:
            return 0.0
        kappa = cross_norm / (r1_norm**3)
        return kappa

    def get_parameters_from_points(self, target_points, u_min=0.0, u_max=1.0):
        """
        任意の点列 `target_points` に対して、その点を最も近く再現するような
        曲線上のパラメータ `u` を最適化により求める。
        """
        target_points = np.asarray(target_points)
        parameters = []

        for pt in target_points:
            def objective(u):
                curve_pt = self(u)  # 曲線上の点
                return np.linalg.norm(curve_pt - pt)  # 目的関数：ユークリッド距離

            result = minimize_scalar(objective, bounds=(u_min, u_max), method='bounded')
            parameters.append(result.x)

        return np.array(parameters)

    def _get_u_mapping(self, sample_points, u_sample, interpolate_points):
        """
        既知のサンプル点 `sample_points` とその対応パラメータ `u_sample` に基づき、
        任意の `interpolate_points` に対応するパラメータ `u` を線形補間で推定。
        """
        sample_points = np.asarray(sample_points)
        u_sample = np.asarray(u_sample)
        interpolate_points = np.asarray(interpolate_points)
        u_interp = np.zeros(len(interpolate_points))

        for i, p in enumerate(interpolate_points):
            distances = np.linalg.norm(sample_points - p, axis=1)
            nearest_indices = np.argsort(distances)[:2]

            i1, i2 = nearest_indices
            p1, p2 = sample_points[i1], sample_points[i2]
            u1, u2 = u_sample[i1], u_sample[i2]
            d1, d2 = distances[i1], distances[i2]

            if d1 + d2 == 0:
                # p がちょうど sample_points の1点に一致する場合
                u_interp[i] = u1
            else:
                # 距離の逆数を重みとした線形補間
                w1 = d2 / (d1 + d2)
                w2 = d1 / (d1 + d2)
                u_interp[i] = w1 * u1 + w2 * u2

        return u_interp

    def _calc_knot_vector(self, n_cp, k, kink_list=[]):
        """
        ノットベクトルを計算する。
        - n_cp: 制御点の数
        - k: 曲線次数
        - kink_list: [(knot_position, multiplicity), ...] 特異点ノットとその重複度のリスト

        ノットベクトルは区間[0,1]内に生成される。
        kink_listが指定されている場合、kink位置での重複ノットを設定。
        """
        # kink_listから位置と重複度を分離
        kink_positions = [t for t, _ in kink_list]
        kink_multiplicities = [m for _, m in kink_list]
        kink_total_mult = sum(kink_multiplicities)

        m = n_cp + k + 1  # ノットベクトルの総数

        # 基本ノット値は [0] + kink_positions + [1]
        knot_values = [0.0] + kink_positions + [1.0]

        # 追加ノット数を計算（内側のノット数）
        m_add = m - 2 * (k + 1) - kink_total_mult

        # 区間長を計算
        interval_lengths = np.diff(np.asarray(knot_values))

        # 区間ごとに分割ノット数を割合で配分（整数化）
        rn_add_list = [int(val * m_add) for val in interval_lengths]

        total = sum(rn_add_list)

        # 合計がm_addに満たない場合は長い区間からノットを追加
        if total < m_add:
            sorted_idx = np.argsort(-interval_lengths)
            for i in range(m_add - total):
                rn_add_list[sorted_idx[i % len(sorted_idx)]] += 1

        # ノットベクトル構築開始（端点ノットはk+1個ずつ）
        t = [0.0] * (k + 1)

        for i in range(len(knot_values) - 1):
            t0, t1 = knot_values[i], knot_values[i + 1]
            seg_count = rn_add_list[i]

            if seg_count > 0:
                # 区間内に等間隔で内部ノットを追加（端点は含まない）
                interior_knots = np.linspace(t0, t1, seg_count + 2)[1:-1]
                t.extend(np.round(interior_knots, 6).tolist())

            # kink_listに該当するノットは重複度分追加
            for knot, mult in kink_list:
                if np.isclose(knot, t1):
                    t.extend([knot] * mult)
                    break

        # 終端ノットもk+1個追加
        t.extend([1.0] * (k + 1))

        return t

    def _remove_kink_from_knot(self, t, k):
        head = t[:k]
        tail = t[-k:]
        middle = t[k:-k]

        seen = set()
        cleaned_middle = []
        duplicates = []

        for val in middle:
            if val not in seen:
                seen.add(val)
                cleaned_middle.append(val)
            else:
                duplicates.append(val)

        cleaned_knot_vector = head + cleaned_middle + tail
        return cleaned_knot_vector, duplicates

    def plot_bspline(self, ref_points=None, show=False, legend=True, figsize=(4, 3), seg=True):
        t = self.t
        c = self.c
        k = self.k

        fig, ax = plt.subplots(figsize=figsize)
        
        # 制御点を点線でプロット
        ax.plot(c[:, 0], c[:, 1], 'o:', markerfacecolor='none', color='k', label='$\mathbf{p}_{i}$')
        
        if seg:
            # 重複を除いたノット値のリスト
            t_unique = np.array(sorted(set(t.tolist())))
            
            # 各ノット区間ごとに曲線を描画
            for i, (t0, t1) in enumerate(zip(t_unique[:-1], t_unique[1:])):
                color = cm.tab10(i % 10)  # 10色のカラーマップをループ使用
                u = np.linspace(t0, t1, 21)
                curve = self(u)  # 曲線上の点群
                
                ax.plot(curve[:, 0], curve[:, 1], color=color)
                
                # 区間端点の描画とラベル
                point_start = self(t0)
                ax.plot(point_start[0], point_start[1], 'x', color=color)
                ax.text(point_start[0] + 0.1, point_start[1], f't: {t0:.3f}', fontsize=8, color='k', ha='left', va='center')
            point_end = self(t1)
            ax.plot(point_end[0], point_end[1], 'x', color=color)

            ax.text(point_end[0] + 0.1, point_end[1], f't: {t1:.3f}', fontsize=8, color='k', ha='left', va='center')
        else:
            t0, t1 = np.min(t), np.max(t)
            u = np.linspace(t0, t1, 101)
            curve = self(u)  # 曲線上の点群
            
            ax.plot(curve[:, 0], curve[:, 1], label='$\mathbf{c}(u)$')
        
        # 参照点がある場合は描画
        if ref_points is not None:
            ref_points = np.asarray(ref_points)
            ax.plot(ref_points[:, 0], ref_points[:, 1], 'o', label='reference points')
        
        if legend:
            ax.legend()
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        ax.grid(True)
        
        if show:
            plt.show()

    def plot_basis_function(self, show=False, legend=True, r_max=None):
        t = self.t
        k = self.k
        if r_max is None:
            r_max = k  # 最大微分次数（ここでは次数kと同じと仮定）

        # fig, axes = plt.subplots(nrows=1, ncols=r_max + 1, figsize=(4 * (r_max + 1) + 1, 3))

        t_unique = np.array(sorted(set(t.tolist())))

        for nu in range(r_max + 1):
            # if r_max == 0:
            #     ax = axes
            # else:
            #     ax = axes.flatten()[nu]
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4 + 1, 3))

            for j, (t0, t1) in enumerate(zip(t_unique[:-1], t_unique[1:])):
                # u_values = np.linspace(t0, t1, 41)[:-1]  # 端点を除く
                u_values = t0 + 0.5 * (1 - np.cos(np.pi * np.linspace(0, 1, 102)[:-1])) * (t1 - t0)
                    
                basis_matrix = self.basis_matrix(u_values, nu=nu)
                for i, row in enumerate(basis_matrix.T):
                    color = cm.tab10(i % 10)
                    if t1 == t[-(k + 1)]:
                        label = f'$N_{{{i}, {k}}}^{nu}(u)$' if nu != 0 else f'$N_{{{i}, {k}}}(u)$' 
                    else:
                        label = None
                    ax.plot(u_values, row, label=label, color=color)

            for knot in t:
                ax.axvline(knot, linestyle=':', color='k', alpha=0.5)

            ax.set_xlabel('parameter u')
            y_label = f'$N_{{i, k}}^{nu}(u)$' if nu != 0 else f'$N_{{i, k}}(u)$' 
            ax.set_ylabel(y_label)
            # ax.set_title(f'$N_{{i, {k}}}^{nu}$')
            ax.grid(True)
            if legend:
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()
        if show:
            plt.show()

    def plot_ddt(self, show=False):
        t = self.t
        k = self.k
        c = self.c
        n_dim = self.dim()
        r_max = k

        fig, axes = plt.subplots(figsize=(3 * r_max, 3 * n_dim), ncols=r_max, nrows=n_dim, squeeze=False)
        
        axis_label = ['x', 'y', 'z']

        for i in range(n_dim):
            for r in range(1, r_max + 1):
                ax = axes[i, r - 1]
                t_unique = np.array(sorted(set(t.tolist())))
                
                for j, (t0, t1) in enumerate(zip(t_unique[:-1], t_unique[1:])):
                    # チェビシェフ点風の分布をパラメータuに設定
                    u_values = t0 + 0.5 * (1 - np.cos(np.pi * np.linspace(0, 1, 22)[:-1])) * (t1 - t0)
                    color = cm.tab10(j % 10)
                    basis_matrix = self.basis_matrix(u_values, nu=r)
                    deriv_vals = basis_matrix @ c[:, i]
                    ax.plot(u_values, deriv_vals, color=color)

                for knot in t:
                    ax.axvline(knot, linestyle=':', color='k', alpha=0.1)
                
                ax.set_xlabel('parameter t')
                ax.set_ylabel('derivative')
                ax.set_title(rf'$\frac{{d^{r}{axis_label[i]}}}{{du^{{{r}}}}}$')
                ax.grid(True)

        plt.tight_layout()
        if show:
            plt.show()

    def plot_ddx(self, axis=0, show=False):
        t = self.t
        k = self.k
        n_dim = self.dim()
        r_max = k  # 最大微分次数（仮定）

        axis_label = ['x', 'y', 'z'][:n_dim]
        # dims = axis_label.copy()
        dims = [x for x in range(n_dim)]
        dims.remove(axis)  # 描画対象の軸以外をリストアップ

        nrows = n_dim - 1
        ncols = r_max
        figsize = (3 * ncols, 3 * nrows)
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False)

        for j, dim in enumerate(dims):
            for r in range(1, r_max + 1):
                ax = axes[j, r - 1]
                vmin, vmax = np.inf, -np.inf

                t_unique = np.array(sorted(set(t.tolist())))
                for i, (t0, t1) in enumerate(zip(t_unique[:-1], t_unique[1:])):
                    # チェビシェフ点風の分布をパラメータuに設定
                    u_values = t0 + 0.5 * (1 - np.cos(np.pi * np.linspace(0, 1, 22)[:-1])) * (t1 - t0)
                    result = np.round(self.ddx(r, u_values, axis=axis)[dim], 6)

                    color = cm.tab10(i % 10)
                    ax.plot(u_values, result, label=f'$(d/d{axis_label[axis]})^{r}$ to {dim}', color=color)

                    vmin = min(vmin, np.min(result))
                    vmax = max(vmax, np.max(result))

                for knot in t:
                    ax.axvline(knot, linestyle=':', color='k', alpha=0.1)

                ax.set_xlabel('parameter t')
                ax.set_ylabel('derivatives')
                ax.set_title(rf'$\left(\frac{{d{{{axis_label[dim]}}}}}{{d{{{axis_label[axis]}}}}}\right)^{{{r}}}$')
                ax.grid(True)

        plt.tight_layout()
        if show:
            plt.show()

    def plot_curvature(self, show=False):
        t = self.t
        c = self.c
        k = self.k

        fig, ax = plt.subplots(figsize=(4, 3))
        t_unique = np.array(sorted(set(t.tolist())))

        for i, (t0, t1) in enumerate(zip(t_unique[:-1], t_unique[1:])):
            color = cm.tab10(i % 10)
            u_values = np.linspace(t0, t1, 101)[:-1]
            curvature_values = np.array([self.curvature(u) for u in u_values])
            ax.plot(u_values, curvature_values, label='Curvature K(t)', color=color)

        for knot in t:
            ax.axvline(knot, linestyle=':', color='k', alpha=0.1)

        ax.set_xlabel('parameter t')
        ax.set_ylabel('Curvature')
        ax.set_title('Curvature along B-spline curve')
        ax.grid(True)
        plt.tight_layout()
        if show:
            plt.show()

    # def max_derivative_order(self):
    #     result = self.k
    #     while result > 0:
    #         try:
    #             self.derivative(nu=result)
    #             break
    #         except ValueError:
    #             result -= 1
    #     return result
