# 必要なモジュールのインポート
from NLFS import Plane  # PlaneクラスをNLFSモジュールからインポート
import time  # 時間計測のためのtimeモジュール
import numpy as np  # 数値計算用のnumpyモジュール
import pandas as pd  # データ操作用のpandasモジュール
import math  # 数学関数用のmathモジュール
from pymoo.core.problem import Problem  # PymooのProblemクラス
from pymoo.algorithms.moo.nsga2 import NSGA2  # NSGA-IIアルゴリズム
from pymoo.operators.sampling.lhs import LHS  # ラテン超方格サンプリング
from pymoo.operators.crossover.sbx import SBX  # シミュレーテッド・バイナリー・クロスオーバー
from pymoo.operators.mutation.pm import PolynomialMutation  # 多項式突然変異
from pymoo.optimize import minimize  # 最適化関数
from pymoo.visualization.scatter import Scatter  # 散布図の可視化
from pymoo.util.display.column import Column  # 表示用カラム
from pymoo.util.display.output import Output  # 表示用出力
from pymoo.config import Config  # Pymooの設定

# Pymooの警告を無効化
Config.warnings['not_compiled'] = False

# 問題設定用のクラス
class MyProblem(Problem):
    def __init__(self, u0=5.5, gamma0=3.5, uvw_gE=np.array([0, 0, 0])):
        self.u0 = u0  # 初期速度
        self.w0 = u0 * math.tan(math.radians(gamma0))  # 初期速度
        self.uvw_gE = uvw_gE  # 対地風速
        # 最適化変数の下限と上限
        xl = np.array([-30, -2, -10,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -5])
        xu = np.array([ -2, -1,   0, 10, 60, 60, 60, +1, +1, +1, +1, +1, +1, +5])
        super().__init__(n_var=14, n_obj=3, n_constr=4, xl=xl, xu=xu)  # 問題の初期化

    def _evaluate(self, X, out, *args, **kwargs):
        # Xを丸める
        X = np.round(X,5)
        # 目的関数のリスト
        f1_list = []
        f2_list = []
        f3_list = []
        g3_list = []
        g4_list = []
        for i in range(X.shape[0]):
            theta0 = X[i, 13]  # 初期ピッチ角
            initial_state = [0, 0, 0, self.u0, 0, self.w0, 0, theta0, 0, 0, 0, 0, 0, 0]  # 初期状態
            airplane = Plane(uvw_gE=self.uvw_gE, X=X[i, :13])  # Planeクラスのインスタンス作成
            solution = airplane.simulate(initial_state=initial_state)  # シミュレーション実行
            XE, YE, ZE = solution.y[:3]  # 結果の取得
            Distance = np.sqrt(XE[-1]**2 + YE[-1]**2)  # 飛行距離

            # 時間差分と速度成分の計算
            dt = np.diff(solution.t)
            uE = np.diff(XE) / dt
            wE = np.diff(ZE) / dt
            gamma = np.degrees(np.arctan(wE / uE))  # 経路角の計算

            # 目的関数1の計算
            target = np.where(solution.t <= X[i, 3], X[i, 0], np.where(solution.t >= X[i, 4], 0, X[i, 1]))
            diff1 = np.sum(np.abs(gamma - target[1:])) / solution.t[-1]

            # 目的関数2の計算
            phi = solution.y[6]
            target = np.where(solution.t <= X[i, 5], phi, np.where(solution.t >= X[i, 6], 0, X[i, 2]))
            diff2 = np.sum(np.abs(phi - target))/ np.abs(solution.t[-1]-X[i,5])

            # 目的関数リストに追加
            f1_list.append(diff1)
            f2_list.append(diff2)
            f3_list.append(-Distance)
            g3_list.append(X[i,4]-solution.t[-1]) # time_flare < flight time
            g4_list.append(X[i,6]-solution.t[-1]) # time_break < flight time

        # 目的関数の結果をnumpy配列に変換
        f1 = np.array(f1_list)
        f2 = np.array(f2_list)
        f3 = np.array(f3_list)

        # 制約条件の計算
        g1 = X[:, 3] - X[:, 4]
        g2 = X[:, 5] - X[:, 6]
        g3 = np.array(g3_list)
        g4 = np.array(g4_list)

        # 出力として目的関数と制約条件を設定
        out["F"] = np.column_stack([f1, f2, f3])
        out["G"] = np.column_stack([g1, g2, g3, g4])

# 出力表示用のカスタムクラス
class MyOutput(Output):
    def __init__(self):
        super().__init__()
        self.Dist_max = Column("Dist_max", width=13)
        self.Dist_ave = Column("Dist_ave", width=13)
        self.diff1 = Column("diff1", width=13)
        self.diff2 = Column("diff2", width=13)
        self.RunTime = Column("RunTime", width=13)
        self.columns += [self.Dist_max, self.Dist_ave, self.diff1, self.diff2, self.RunTime]
        self.start_time = time.time()  # 開始時刻を記録

    def update(self, algorithm):
        super().update(algorithm)
        self.diff1.set(np.min(algorithm.pop.get("F"), axis=0)[0])
        self.diff2.set(np.min(algorithm.pop.get("F"), axis=0)[1])
        self.Dist_max.set(-np.min(algorithm.pop.get("F"), axis=0)[2])
        self.Dist_ave.set(-np.average(algorithm.pop.get("F"), axis=0)[2])
        self.RunTime.set(time.time()-self.start_time)

if __name__ == '__main__':
    start_time = time.time()  # 開始時刻を記録

    # 最適化問題の設定
    problem = MyProblem(
        u0=5.5,
        gamma0=3.5,
        uvw_gE=np.array([0, 2*math.sqrt(2), 0]),
    )

    # NSGA-IIアルゴリズムの設定
    algorithm = NSGA2(
        pop_size=100,
        sampling=LHS(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PolynomialMutation(eta=20),
        eliminate_duplicates=True,
    )

    # 最適化の実行
    res = minimize(
        problem,
        algorithm,
        ("n_gen", 100),
        seed=1,
        save_history=False,
        verbose=True,
        output=MyOutput(),
    )
    
    res_X = pd.DataFrame(np.round(res.X,5), columns=[f'X{i}' for i in range(14)])
    res_F = pd.DataFrame(res.F, columns=['diff1', 'diff2', 'Distance'])
    res_G = pd.DataFrame(res.G, columns=['Time_Cruise', 'Time_Turn'])

    # 結果をCSVファイルに保存
    res = pd.concat([res_X, res_F, res_G], axis=1)  # 変数データと目的関数データを結合
    res = res.sort_values('Distance')  # 距離でソート
    res[[f'X{i}' for i in range(14)]].to_csv('./res_X.csv', index=False)
    res[['Distance','diff1', 'diff2']].to_csv('./res_F.csv', index=False)
    res[['Time_Cruise', 'Time_Turn', 'Time_Flare', 'Time_NoBank']].to_csv('./res_G.csv', index=False)

    # 終了時刻を記録し、実行時間を計算
    end_time = time.time()
    print(f"Run time: {end_time - start_time} s")
