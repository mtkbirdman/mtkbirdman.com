# 必要なモジュールのインポート
from NLFS import Plane  # NLFSモジュールからPlaneクラスをインポート
import math  # 数学関数用のmathモジュール
import numpy as np  # 数値計算用のnumpyモジュール
import pandas as pd  # データ操作用のpandasモジュール
from scipy.spatial.transform import Rotation  # 座標変換のためのライブラリ
import matplotlib.pyplot as plt  # プロット用のmatplotlibモジュール
import matplotlib.gridspec as gridspec  # サブプロットのレイアウト用のgridspecモジュール

# 3Dプロットを行う関数
def plot_3D(solution):
    fig = plt.figure(figsize=(16/2, 9/2))  # 図のサイズを指定して作成
    ax = fig.add_subplot(111, projection='3d')  # 3Dプロット用のサブプロットを追加
    ax.plot(solution.y[0], solution.y[1], 10 - solution.y[2], color='blue')  # 3Dプロット
    ax.set_xlabel('X')  # X軸のラベルを設定
    ax.set_ylabel('Y')  # Y軸のラベルを設定
    ax.set_zlabel('hE')  # Z軸のラベルを設定
    xmax = max(np.max(solution.y[0]), np.max(np.abs(solution.y[1])))  # X軸の最大値を計算
    ymin = min(-np.max(solution.y[0]) / 2, np.min(solution.y[1]))  # Y軸の最小値を計算
    ymax = ymin + xmax  # Y軸の最大値を計算
    ax.set_xlim(0, max(np.max(solution.y[0]), np.max(np.abs(solution.y[1]))))  # X軸の範囲を設定
    ax.set_ylim(ymin, ymax)  # Y軸の範囲を設定
    ax.set_zlim(0, 20)  # Z軸の範囲を設定
    ax.invert_yaxis()  # Y軸を反転
    ax.view_init(elev=60, azim=-180)  # プロットの視点を設定
    plt.show()  # 図を表示

# 2Dプロットを行う関数
def plot_2D(solutions):
    fig = plt.figure()  # 図を作成
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])  # サブプロットのレイアウトを設定

    # 左側のサブプロットを設定
    ax2 = fig.add_subplot(gs[0])
    ax2.set_ylabel('X')  # Y軸のラベルを設定
    ax2.set_xlabel('Y')  # X軸のラベルを設定
    ax2.set_ylim(0, 500)  # Y軸の範囲を設定
    ax2.set_xlim(-400, 100)  # X軸の範囲を設定
    ax2.set_yticks(np.linspace(0, 500, 6))  # Y軸の目盛りを設定
    ax2.set_xticks(np.linspace(-400, 100, 6))  # X軸の目盛りを設定
    ax2.set_aspect(aspect='equal')  # アスペクト比を1:1に設定

    # 右側のサブプロットを設定
    ax1 = fig.add_subplot(gs[1])
    ax1.set_ylabel('X')  # Y軸のラベルを設定
    ax1.set_xlabel('hE')  # X軸のラベルを設定
    ax1.set_ylim(0, 500)  # Y軸の範囲を設定
    ax1.set_xlim(0, 10)  # X軸の範囲を設定
    ax1.set_yticks(np.linspace(0, 500, 6))  # Y軸の目盛りを設定
    ax1.set_xticks(np.linspace(0, 10, 3))  # X軸の目盛りを設定
    ax1.set_aspect(aspect=0.1)  # アスペクト比を1:10に設定
    ax1.invert_xaxis()  # X軸を反転

    for solution in solutions:
        ax1.plot(10 - solution.y[2], solution.y[0])  # hE対Xのプロット
        ax2.plot(solution.y[1], solution.y[0])  # Y対Xのプロット
    
    # レイアウトの調整
    ax1.grid()  # グリッドを表示
    ax2.grid()  # グリッドを表示
    plt.tight_layout()  # レイアウトを自動調整

    plt.savefig('graph_position.png')
    plt.show()  # 図を表示

def plot_2Ds(solutions):
    solutions = []
    for N in range(11):
        # 初期条件の設定
        X = res_X.values[N]
        theta0 = X[13]  # 最適化結果から初期ピッチ角を取得
        initial_state = [0, 0, 0, u0, 0, w0, 0, theta0, 0, 0, 0, 0, 0, 0]  # 初期状態

        # 最も飛距離に優れた個体のシミュレーションの実行
        plane = Plane(uvw_gE=uvw_gE, X=X)  # Planeクラスのインスタンスを作成
        solution = plane.simulate(initial_state=initial_state)  # シミュレーションを実行
        XE, YE, ZE = solution.y[:3,-1]  # 結果の取得
        print('Distance {:6.3f} Time {:6.3f}'.format(math.sqrt(XE**2 + YE**2), solution.t[-1]))
        solutions.append(solution)

    plot_2D(solutions=solutions)

def plot_results(df):
    xlabel = 't' 
    ylabels =['hE', 'V','alpha', 'beta', 'p', 'q','r', 'gamma', 'VE', 'el', 'rd', 'phi', 'theta', 'psi'] 
    cols=2
    rows=len(ylabels)//cols 

    fig, axes = plt.subplots(rows, cols,figsize=(16,9), sharex=True) # 770 EX 
    for n_graph,ylabel in enumerate (ylabels): # ZYKL 
        i=n_graph%rows 
        j=n_graph//rows 
        ax = axes[i,j] # y軸が複数の場合、各y 軸に対応するサブプロットを使用する 
        ax.plot(df[xlabel].astype(float), df[ylabel].astype (float), label=ylabel) 
        if i == rows-1: 
            ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.minorticks_on() # 補助目盛を表示 
        ax.grid(axis='both', which='both')
        dX = ax.get_xticks() [1] - ax.get_xticks() [0] # X 
        dY = ax.get_yticks() [1] - ax.get_yticks() [0] # Y 
        xminorticks_interval = dX / 10
        yminorticks_interval = dY / 10
        ax.xaxis.set_minor_locator(plt.MultipleLocator(xminorticks_interval)) 
        ax.yaxis.set_minor_locator(plt.MultipleLocator(yminorticks_interval))
        ax.tick_params (axis='both', which='major')
        ax.tick_params (axis='both', which='minor', width=0.2,grid_linewidth=0.2) 
    plt.tight_layout()
    plt.savefig('graph_result.png')
    plt.show() #グラフを表示 

# CSVファイルからデータを読み込み
res_X = pd.read_csv('./res_X.csv')  # 最適化結果の変数データ
res_F = pd.read_csv('./res_F.csv')  # 最適化結果の目的関数データ

# 初期条件の設定
N = 0
u0 = 5.5  # 初期速度
w0 = u0 * math.tan(math.radians(3.5))  # 初期降下率
uvw_gE=np.array([2*math.sqrt(2), 0, 0])  # 風速
X = res_X.values[N]
theta0 = X[13]  # 最適化結果から初期ピッチ角を取得
initial_state = [0, 0, 0, u0, 0, w0, 0, theta0, 0, 0, 0, 0, 0, 0]  # 初期状態

# 最も飛距離に優れた個体のシミュレーションの実行
plane = Plane(uvw_gE=uvw_gE, X=X)  # Planeクラスのインスタンスを作成
solution = plane.simulate(initial_state=initial_state)  # シミュレーションを実行
XE, YE, ZE = solution.y[:3,-1]  # 結果の取得

# シミュレーション結果の表示
print('XE {:6.3f} YE {:6.3f} ZE {:6.3f}'.format(XE, YE, 10 - ZE, solution.t[-1]))
print('Distance {:6.3f} Time {:6.3f}'.format(math.sqrt(XE**2 + YE**2), solution.t[-1])+'\n')
print('dive   {:6.3f} deg from {:6.3f} to {:6.3f}'.format(X[0], 0, X[3]))
print('cruise {:6.3f} deg from {:6.3f} to {:6.3f}'.format(X[1], X[3], X[4]))
print('trun   {:6.3f} deg from {:6.3f} to {:6.3f}'.format(X[2], X[5], X[6]))

# シミュレーション結果をCSVファイルに保存
columns = ['t', 'XE', 'YE', 'ZE', 'u', 'v', 'w', 'phi', 'theta', 'psi', 'p', 'q', 'r', 'el', 'rd']
df = pd.DataFrame(np.concatenate([solution.t.reshape(1, -1).T, solution.y.T], axis=1), columns=columns)
df['hE'] = 10-df['ZE'].values 
uvw_g = []
for i in range(len(df['hE'])):
    euler = solution.y[6:9,i]
    rot = Rotation.from_euler('XYZ', -euler, degrees=True)
    uvw_g.append(rot.apply(uvw_gE*((df['hE'].iloc[i]/10)**(1/7)))) # 風速勾配 1/7乗則
ug, vg, wg = np.array(uvw_g).T
df['V'] = np.sqrt((df['u'].values+ug)**2+(df['v'].values+vg)**2+(df['w'].values+wg)**2) 
df['alpha'] = np.degrees(np.arctan((df['w'].values+wg)/(df['u'].values+ug))) 
df['beta'] = np. degrees (np.arctan((df['v'].values+vg)/(df['V'].values))) 
df['gamma'] = df['theta'] - df['alpha']
dt=np.diff(df['t'].values) 
uE=np.diff (df['XE'].values)/dt 
vE=np.diff (df['YE'].values)/dt 
wE=np.diff (df['ZE'].values)/dt 
df['VE'] = np.concatenate ((np.sqrt(uE**2+vE**2+wE**2), np.sqrt(uE**2+vE**2+wE**2) [-1].reshape(-1)), axis=0) 
df['ug'] = ug
df['vg'] = vg
df['wg'] = wg
df.to_csv('./result.csv')

# 2Dプロットの表示
plot_2D(solutions=[solution])
plot_results(df=df)
