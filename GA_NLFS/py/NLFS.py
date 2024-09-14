import time  # 時間管理のためのライブラリ
import math  # 数学関数のためのライブラリ
import numpy as np  # 数値計算のためのライブラリ
from scipy.integrate import solve_ivp  # 常微分方程式の数値解法
from scipy.spatial.transform import Rotation  # 座標変換のためのライブラリ
import matplotlib.pyplot as plt  # グラフ描画のためのライブラリ
from mpl_toolkits.mplot3d import Axes3D  # 3Dプロットのためのライブラリ

# 飛行機クラスの定義
class Airplane:
    def __init__(self, uvw_gE=np.array([0, 0, 0]), X=np.zeros(13)):
        # 機体の基本的なパラメータ
        self.mass = 100  # 期待質量 [kg]
        self.airspeed0 = 10  # 釣り合い対気速度 [m/s]
        self.alpha0 = 1.45  # 釣り合い迎え角 [度]
        self.CDp0 = 0.015  # 抗力係数
        self.Cmw0 = -0.115  # モーメント係数
        self.CLMAX = 1.700  # 最大揚力係数
        self.Sw = 18  # 翼面積 [m^2]
        self.bw = 25  # 翼幅 [m]
        self.cMAC = 0.75  # 平均翼弦長 [m]
        self.aw = 0.105  # 揚力傾斜
        self.hw = (0.333-0.250)  # 重心位置
        self.ew = 0.985  # 飛行機効率
        self.AR = (self.bw*self.bw)/self.Sw  # アスペクト比
        self.Downwash = False  # ダウンウォッシュ効果
        self.St = 1.500  # 垂直尾翼面積 [m^2]
        self.at = 0.080  # 垂直尾翼揚力傾斜
        self.lt = 3.200  # 垂直尾翼のアーム長 [m]
        self.dh_max = 0.5  # 最大重心移動量 [-]
        self.el_max = 10  # 最大エレベータ舵角 [deg]
        self.tau = 1.000  # エレベータ舵角効率
        self.VH = (self.St*self.lt)/(self.Sw*self.cMAC)  # 垂直尾翼容積比
        self.rd_max = 15  # 最大ラダー舵角 [deg]
        self.CGEMIN = 0.25  # 最小地面効果係数
        self.Cyb = -0.0036
        self.Cyp = -0.4500
        self.Cyr = 0.1420 
        self.Cydr = 0.0018
        self.Clb = -0.0040
        self.Clp = -0.8200
        self.Clr = 0.2200 
        self.Cldr = 0.0000
        self.Cnb = -0.0005
        self.Cnp = -0.1300
        self.Cnr = -0.0025
        self.Cndr = -0.0003
        self.Ixx = 1000  # 慣性モーメント [kg*m^2]
        self.Iyy = 70  # 慣性モーメント [kg*m^2]
        self.Izz = 1000  # 慣性モーメント [kg*m^2]
        self.Izx = -8  # 慣性モーメント [kg*m^2]
        self.Ixz = 8  # 慣性モーメント [kg*m^2]
        self.epsilon0 = 0  # 釣り合い吹き下ろし角
        self.gravity = 9.81  # 重力加速度 [m/s^2]
        self.rho = 1.225  # 空気密度 [kg/m^3]
        self.uvw_gE = uvw_gE  # 地表風ベクトル
        self.PID_Lon = np.zeros(3)  # 縦方向PID制御の係数
        self.PID_lat = np.zeros(3)  # 横方向PID制御の係数
        self.t_old = 0  # 前回の時間ステップ
        self.dh = 0  # 重心移動量

        self.gamma = None
        self.phi = None
        self.gamma_dive = X[0]  # ダイブ角
        self.gamma_cruise = X[1]  # 巡航角
        self.phi_turn = X[2]  # 旋回角
        self.time_pullup = X[3]  # 引き起こし時刻 [s]
        self.time_flare = X[4]  # フレア時刻 [s]
        self.time_entry = X[5]  # 旋回開始時刻 [s]
        self.time_break = X[6]  # 旋回終了時刻 [s]
        self.PID_param_lon = X[7:10]  # 縦方向PIDパラメータ
        self.PID_param_lat = X[10:13]  # 横方向PIDパラメータ

    # 飛行機の運動方程式
    def motion_equations(self, t, state):
        # 状態変数の展開
        XE, YE, ZE, u, v, w, phi, theta, psi, p, q, r , el, rd= state
        uvw = np.array([u, v, w], dtype=float)
        hE = max([10+ZE, 0.001])  # 地表高度の計算
        euler = np.array([phi, theta, psi], dtype=float)
        self.phi = phi
        
        # 機体座標系から地表座標系への変換
        rot = Rotation.from_euler('ZYX', np.flip(euler), degrees=True)
        uE, vE, wE = rot.apply(uvw)
        
        # 地表座標系から機体座標系への変換
        rot = Rotation.from_euler('XYZ', -euler, degrees=True)
        ug, vg, wg = rot.apply(self.uvw_gE*((hE/10)**(1/7))) # 風速勾配 1/7乗則
        
        # 対気速度の計算
        Airspeed = np.sqrt((u+ug)**2+(v+vg)**2+(w+wg)**2)
        
        # 迎え角と経路角の計算
        if abs((w+wg)/(u+ug)) < 1:
            alpha = math.degrees(math.asin((w+wg)/(u+ug)))
        else:
            alpha = np.sign((w+wg)/(u+ug))*90
        self.gamma = theta-alpha 

        # 横滑り角の計算
        beta = math.degrees(math.asin((v+vg)/Airspeed))

        # 揚力係数の計算
        CGE = (self.CGEMIN+33*(hE/self.bw)**1.5)/(1+33*(hE/self.bw)**1.5)
        CL0 = (self.mass*self.gravity)/(0.5*self.rho*self.airspeed0**2*self.Sw)
        CLt0 = (self.Cmw0+CL0*self.hw)/(self.VH+(self.St/self.Sw)*self.hw)
        CLw0 = CL0-(self.St/self.Sw)*CLt0
        if self.Downwash:
            self.epsilon0 = (CL0/(math.pi*self.ew*self.AR))*math.degrees(1)
        CLw = CLw0+self.aw*(alpha-self.alpha0)
        CLt = (CLt0+self.at*((alpha-self.alpha0)+(1-CGE*(CLw/CLw0))*self.epsilon0+el*self.tau+((self.lt-self.dh*self.cMAC)/Airspeed)*q)) 
        
        # 最大揚力係数の制限
        if abs(CLw) > self.CLMAX: 
            CLw = (CLw/abs(CLw))*self.CLMAX # Stall 
        if abs(CLt) > self.CLMAX: 
            CLt = (CLt/abs(CLt))*self.CLMAX # Stall 
        CL = CLw+(self.St/self.Sw)*CLt # CL 

        # 抗力係数の計算
        CD = (self.CDp0*(math.cos(math.radians(alpha-self.alpha0))**2+(4*2.5)/(math.pi*0.6)*(1/0.04)*math.sin(math.radians(alpha-self.alpha0))**2)
              +((CL*CL)/(math.pi*self.ew*self.AR))*CGE) # CD 
        
        # 動微係数を安定軸から機体軸へ変換
        cosA = math.cos(math.radians(alpha))
        sinA = math.sin(math.radians(alpha))
        Cyp = self.Cyp*cosA-self.Cyr*sinA
        Cyr = -self.Cyp*sinA+self.Cyr*cosA
        Clp = self.Clp*cosA*cosA-self.Cnp*sinA*cosA-self.Clr*sinA*cosA+self.Cnr*sinA*sinA
        Cnp = self.Clp*sinA*cosA+self.Cnp*cosA*cosA-self.Clr*sinA*sinA-self.Cnr*sinA*cosA
        Clr = self.Clp*sinA*cosA-self.Cnp*sinA*sinA+self.Clr*cosA*cosA-self.Cnr*sinA*cosA
        Cnr = self.Clp*sinA*sinA+self.Cnp*sinA*cosA+self.Clr*sinA*cosA+self.Cnr*cosA*cosA

        # 機体軸の空力係数の計算
        Cx = CL*sinA-CD*cosA 
        Cy = (self.Cyb*beta+Cyp*(1/math.degrees(1))*((p*self.bw)/(2*Airspeed))+Cyr*(1/math.degrees(1))*((r*self.bw)/(2*Airspeed))+ self.Cydr*rd) # Cy 
        Cz = -CL*cosA-CD*sinA # Cz 
        Cl = (self.Clb*beta+Clp*(1/math.degrees(1))*((p*self.bw)/(2*Airspeed))+Clr*(1/math.degrees(1))*((r*self.bw)/(2*Airspeed))+ self.Cldr*rd) # CL 
        Cm = self.Cmw0+CLw*self.hw-self.VH*CLt+CL*self.dh # Cm 
        Cn = (self.Cnb*beta+ Cnp*(1/math.degrees(1))*((p*self.bw)/(2*Airspeed))+Cnr*(1/math.degrees(1))*((r*self.bw)/(2*Airspeed))+ self.Cndr*rd) # Cn 
        
        # 機体にはたらく力とモーメントの計算
        X = 0.5*self.rho*Airspeed**2*self.Sw*Cx 
        Y = 0.5*self.rho*Airspeed**2*self.Sw*Cy 
        Z = 0.5*self.rho*Airspeed**2*self.Sw*Cz 
        L = ((self.Iyy-self.Izz)*math.radians(q*r)+self.Ixz*math.radians(p*q)+math.degrees(0.5*self.rho*Airspeed**2*self.Sw*self.bw*Cl)) 
        M = ((self.Izz-self.Ixx)*math.radians(r*p)+self.Ixz*math.radians(r**2-p**2)+math.degrees(0.5*self.rho*Airspeed**2*self.Sw*self.cMAC*Cm))
        N = ((self.Ixx-self.Iyy)*math.radians(p*q)-self.Ixz*math.radians(q*r)+math.degrees(0.5*self.rho*Airspeed**2*self.Sw*self.bw*Cn)) 
        
        # 絶対座標系の速度、機体座標系の加速度、姿勢角速度、角加速度の計算
        rad_phi = math.radians(phi) 
        rad_theta = math.radians(theta) 
        dphi = (p+(r*math.cos(rad_phi)+q*math.sin(rad_phi))*math.tan(rad_theta)) 
        dtheta = (q*math.cos(rad_phi)-r*math.sin(rad_phi)) 
        dpsi = ((r*math.cos(rad_phi)+q*math.sin(rad_phi))/math.cos(rad_theta)) 
        du = (-math.radians(q)*w+math.radians(r)*v-self.gravity*math.sin(rad_theta)+(X/self.mass)) 
        dv = (-math.radians(r)*u+math.radians(p)*w+self.gravity*math.cos(rad_theta)*math.sin(rad_phi)+(Y/self.mass)) 
        dw = (-math.radians(p)*v+math.radians(q)*u+self.gravity*math.cos(rad_theta)*math.cos(rad_phi)+(Z/self.mass)) 
        dp = ((L/self.Ixx)+(self.Ixz/self.Ixx)*(N/self.Izz))/(1-(self.Ixz**2)/(self.Izz*self.Ixx)) 
        dq = M/self.Iyy 
        dr = ((N/self.Izz)+(self.Ixz/self.Izz)*(L/self.Ixx))/(1-(self.Ixz**2)/(self.Ixx*self.Izz)) 
        
        # 入力の計算
        d_el, d_rd = self._calculate_inputs(t, state)
        # 最大舵角の制限 
        if np.abs(el) > self.el_max: 
            d_el = 0
        if np.abs(rd) > self.rd_max: 
            d_rd = 0
        
        # 時間の更新
        self.t_old = t         
        
        return [uE, vE, wE, du, dv, dw, dphi, dtheta, dpsi, dp, dq, dr, d_el, d_rd] 
    
    def _calculate_inputs(self, t, state): 
        # 操舵量の計算する
        if t-self.t_old <=  0: 
            return 0, 0
        else: 
            dt = t-self.t_old 
            
            # 縦方向の制御
            pullup_flag = True 
            flare_flag = True 
            if t < self.time_pullup: # 引き起こし
                gamma_target = self.gamma_dive
            elif t < self.time_flare: # 定常滑空
                if pullup_flag: 
                    self.PID_Lon[1] = 0 # Iをリセット
                    pullup_flag = False 
                gamma_target = self.gamma_cruise
            else: # フレア 
                if flare_flag: 
                    self.PID_Lon[1] = 0 # Iをリセット 
                    flare_flag = False 
                gamma_target = 0

            self.PID_Lon[2] = ((self.gamma-gamma_target)-self.PID_Lon[0])/dt 
            self.PID_Lon[0] = self.gamma-gamma_target
            self.PID_Lon[1] +=  self.PID_Lon[0]*dt 

            # 横方向の制御 
            turn_flag = True
            if t < self.time_entry: # 操縦しない            
                phi_target = self.phi
            elif t < self.time_break: # 旋回中
                phi_target = self.phi_turn
            else: # 水平飛行
                if turn_flag: 
                    self.PID_lat[1] = 0 # Iをリセット
                    turn_flag = False
                phi_target = 0
                
            self.PID_lat[2] = ((self.phi-phi_target)-self.PID_lat[0])/dt 
            self.PID_lat[0] = (self.phi-phi_target)
            self.PID_lat[1] +=  self.PID_lat[0]*dt 

            # 操舵変化量の計算 
            if t < self.time_pullup:
                d_el = np.dot(self.PID_param_lon, self.PID_Lon)*self.el_max 
            elif t < self.time_flare:
                d_el = np.dot(self.PID_param_lon, self.PID_Lon)*self.el_max 
            else: 
                d_el = np.dot(self.PID_param_lon, self.PID_Lon)*self.el_max 

            if t < self.time_entry or t > self.time_break: 
                d_rd = np.dot(self.PID_param_lat, self.PID_lat)*self.rd_max 
            else: 
                d_rd = np.dot(self.PID_param_lat, self.PID_lat)*self.rd_max 
            
            # 操舵変化量の制限
            d_el = min(max(d_el,-self.el_max/1.0), self.el_max/1.0) 
            d_rd = min(max(d_rd,-self.rd_max/1.0), self.rd_max/1.0) 
            
            return d_el, d_rd
        
    def simulate(self, initial_state, t_start = 0, t_end = 100, first_step = 1/30): 
        # 解析時間
        t_span = [t_start, t_end] 
        # 終了条件
        landing = Landing() 
        stall = Stall() 
        reverse = Reverse() 
        overbank = OverBank() 
        landing.terminal = True 
        stall.terminal = True 
        reverse.terminal = True 
        overbank.terminal = True 
        # 運動方程式を解く 
        solution = solve_ivp(self.motion_equations , t_span, initial_state , events = [landing, stall, reverse, overbank] , first_step = first_step) 
        return solution

class Landing: # 着水
    def __init__(self): 
        pass 
    def __call__(self, t, state): 
        XE, YE, ZE, u, v, w, phi, theta, psi, p, q, r, el, rd = state 
        return 10-ZE 

class Stall: # 失速
    def _init__(self): 
        pass 
    def __call__(self, t, state): 
        XE, YE, ZE, u, v, w, phi, theta, psi, p, q, r, el, rd = state 
        return w/u-math.sin(math.radians(30)) 

class Reverse: # 逆走
    def __init__(self): 
        pass 
    def __call__(self, t, state): 
        XE, YE, ZE, u, v, w, phi, theta, psi, p, q, r, el, rd = state 
        uvw = np.array([u, v, w], dtype=float)
        euler = np.array([phi, theta, psi], dtype=float)
        
        # 機体座標系から地表座標系への変換
        rot = Rotation.from_euler('ZYX', np.flip(euler), degrees=True)
        uE, vE, wE = rot.apply(uvw)
        return uE

class OverBank: #過度のバンク 
    def __init__(self): 
        pass 
    def __call__(self, t, state): 
        XE, YE, ZE, u, v, w, phi, theta, psi, p, q, r, el, rd = state 
        return abs(phi)-45 
    
if __name__ == "__main__":
    # シミュレーションの開始時刻を記録
    start_time = time.time()
    
    # Airplaneクラスのインスタンスを作成し、外部風速をゼロに設定
    plane = Airplane(uvw_gE=np.array([0, 0, 0]))
    
    # シミュレーションを実行し、初期状態を設定
    # 初期状態: [XE, YE, ZE, u, v, w, phi, theta, psi, p, q, r, el, rd]
    # XE, YE, ZE: 位置座標, u, v, w: 機体軸速度成分, phi, theta, psi: 姿勢オイラー角, p, q, r: 角速度成分, el, rd: エレベータ舵角、ラダー舵角
    solution = plane.simulate(initial_state=[0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        
    # シミュレーション結果の最終時刻における位置座標を出力
    print(solution.y[0, -1], solution.y[1, -1], 10-solution.y[2, -1])
    
    # シミュレーションの終了時刻を記録し、実行時間を計算
    end_time = time.time()
    print(f"Run time: {end_time-start_time} s")
    
    # 結果をプロットするための図を設定
    fig = plt.figure(figsize=(16/2, 9/2))
    ax = fig.add_subplot(111, projection='3d')
    
    # シミュレーション結果の3D散布図を作成
    ax.scatter(solution.y[0], solution.y[1], 10-solution.y[2], color='gray')
    
    # 軸ラベルを設定
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # x軸とy軸の表示範囲を設定
    xmax = max(np.max(solution.y[0]), np.max(np.abs(solution.y[1])))
    ymin = min(-np.max(solution.y[0])/2, np.min(solution.y[1]))
    ymax = ymin+xmax
    ax.set_xlim(0, max(np.max(solution.y[0]), np.max(np.abs(solution.y[1]))))
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(0, 20)
    
    # y軸を反転して表示
    ax.invert_yaxis()
    
    # 3Dプロットの視点を設定
    ax.view_init(elev=60, azim=-180)
    
    # プロットを表示
    plt.show()