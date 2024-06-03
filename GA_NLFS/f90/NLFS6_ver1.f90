program Analysis
    implicit none
    integer::type !どの解析を行うか

    write(*,*) "What Analysis?"
    write(*,*) "0: NLFS6"
    write(*,*) "1: GA_NLFS6"
    read(*,*) type

    if(type==0) call Test_NLFS()
    if(type==1) call GA_NLFS()

end program 

subroutine Test_NLFS()
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan !プログラム中で使用する関数の名前を宣言
    real(8)::Distance !飛行距離 [m]
    integer,parameter::genom_length = 21 !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8),dimension(0:14,0:2)::NLFS_txt !NLFS.txtの中身
    real(8),dimension(0:11,0:100)::Plane_txt !NLFS.txtの中身
    integer,parameter::FLightLog_data = 47 !飛行解析のデータ数
    real(8),dimension(:,:),allocatable::FlightLog !飛行解析のデータ
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    integer::iteration_max !最大反復回数 [s]
    integer::i_max_copy !最大反復回数のコピー [s]
    integer::pilot_method !操縦方式．0：エレベータ，1：重心移動
    !カウンター
    integer::i,num
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min,CPU_time_h !initial time, finish time, time rate
    
    call system_clock(initial_time) !開始時間の読み込み

    !遺伝子情報の読み込み
    open(10,file="LIST.txt") !ファイルを開く
    read(10,*) genom_list(0),genom_list(1),genom_list(2)
    do i=3,8
        read(10,*) genom_list(i) !1行ずつ値を読み込む
    end do
    do i=0,3
        read(10,*) genom_list(9+i*3),genom_list(10+i*3),genom_list(11+i*3)
    end do
    close(10) !ファイルを閉じる

    write(*,*) "-----------------------------------------------------------------------------------"
    write(*,'(A15)') "GENOM"
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "phi","=",genom_list(0),"theta","=",genom_list(1),"psi","=",genom_list(2)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_pullup","=",genom_list(3),"dive_angle","=",genom_list(6),"cluse_angle","=",genom_list(7)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_bank","=",genom_list(4),"time_bank_end","=",genom_list(5),"bank_angle","=",genom_list(8)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[0]","=",genom_list(9),"PID_lon2[0]","=",genom_list(12)&
    ,"PID_lat1[0]","=",genom_list(15),"PID_lat2[0]","=",genom_list(18)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(10),"PID_lon2[1]","=",genom_list(13)&
    ,"PID_lat1[1]","=",genom_list(16),"PID_lat2[1]","=",genom_list(19)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(11),"PID_lon2[2]","=",genom_list(14)&
    ,"PID_lat1[2]","=",genom_list(17),"PID_lat2[2]","=",genom_list(20)
    write(*,*)

    open(10,file="NLFS.txt") !ファイルを開く
        do i=0, 14
            read(10,*) NLFS_txt(i,0:2)
        end do
        pilot_method=int(NLFS_txt(0,0)) !操縦方式
        time_max=NLFS_txt(1,0) !最大飛行時間 [s]
        time_step=NLFS_txt(2,0) !タイムステップ [s]
    close(10) !ファイルを閉じる
    !機体データ読み込み
    Plane_txt=0.0D0
    open(10,file="Plane.txt") !ファイルを開く
        do i=0, 11
            if(i<5 .or. 7<i) then
                read(10,*) Plane_txt(i,0:16)
            else
                num=int((Plane_txt(4,0)+1)*(Plane_txt(4,1)+1)*(Plane_txt(4,2)+1))
                read(10,*) Plane_txt(i,0:num-1)
            end if
        end do
    close(10) !ファイルを閉じる

    iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
    allocate(FlightLog(0:iteration_max,0:FLightLog_data))

    i_max_copy=iteration_max
    call NLFS(Distance,genom_length,genom_list,i_max_copy,pilot_method,FlightLog_data,FlightLog,NLFS_txt,Plane_txt)
    
    call Output_FlightLog(i_max_copy,iteration_max,FlightLog_data,FlightLog,genom_length,genom_list) !飛行解析データをcsvファイルで書き出し

    write(*,*) 
    write(*,'("Flight Distance =",F8.3,"[m]")') Distance
    write(*,*) 
    
    write(*,*) "------------------------------------------------------------"
    call system_clock(finish_time,time_rate) !終了時間の読み込み
    CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
    CPU_time_min=0; CPU_time_h=0
    if(CPU_time_s>60) then
        CPU_time_min=int(CPU_time_s/60.0D0)
        CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        if(CPU_time_min>60) then
            CPU_time_h=int(CPU_time_min/60)
            CPU_time_min=CPU_time_min-60*CPU_time_h
        end if
    end if
    write(*,'(A24,I12,A12,I12,A12,F12.3,A12)') "Computation Time",CPU_time_h,"[h]",CPU_time_min,"[min]",CPU_time_s,"[s]"

end subroutine  

subroutine NLFS(Distance,genom_length,genom_list,iteration_max,pilot_method,FlightLog_data,FlightLog,NLFS_txt,Plane_txt)
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan,mod !プログラム中で使用する関数の名前を宣言
    !パラメーター
    real(8),parameter::relaxation_coef=1.00D0 !操舵の緩和係数．舵角の急激な変化を抑える
    real(8),parameter::weighting_factor1=0.000D0 !飛距離と偏差の総和の重みづけ．大きければ大きいほど偏差の総和の少なさが優先される．2.5
    real(8),parameter::weighting_factor2=1.00D0 !縦と横・方向の偏差の総和の重みづけ．大きければ大きいほど横・方向の少なさが優先される．10
    real(8),parameter::weighting_factor3=0.000D0 !飛距離に対する入力の時間変化率の総和の重みづけ．大きければ大きいほど入力が滑らかになることが優先される．
    !引数
    real(8),intent(OUT)::Distance !飛行距離
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1),intent(IN)::genom_list !遺伝子情報
    integer,intent(INOUT)::iteration_max !最大反復回数
    integer,intent(INOUT)::pilot_method !操縦方式．0：エレベータ，1：重心移動
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    real(8),dimension(0:iteration_max,0:FlightLog_data),intent(OUT)::FlightLog !飛行解析のデータ
    real(8),dimension(0:14,0:2),intent(IN)::NLFS_txt !NLFS.txtの中身
    real(8),dimension(0:11,0:100),intent(IN)::Plane_txt !Plane.txtの中身
    !非線形解析用
    real(8)::Lift !揚力 [N]
    real(8)::Drag !抗力 [N]
    real(8),dimension(0:2)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    real(8),dimension(0:2)::Position !位置ベクトル [xE,yE,hE] [m]
    real(8),dimension(0:2)::Velocity !速度ベクトル [dxE/dt,dyE/dt,dhE/dt] [m/s]
    real(8),dimension(0:2)::Velocity_old !速度ベクトル [dxE/dt,dyE/dt,dhE/dt] [m/s]
    real(8),dimension(0:2)::Altitude_angle !姿勢角 [Φ,θ,Ψ] [deg]
    real(8),dimension(0:2)::Altitude_angular_velocity !姿勢角速度 [dΦ/dt,dθ/dt,dΨ/dt] [deg]
    real(8),dimension(0:2)::Alt_ang_vel_old !姿勢角速度 [dΦ/dt,dθ/dt,dΨ/dt] [deg]
    real(8),dimension(0:2)::Ground_speed !機体軸における対地速度ベクトル [u,v,w] [m/s]
    real(8),dimension(0:2)::Acceleration !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2)::Acceleration_old !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2)::Angular_velocity !機体軸における角速度ベクトル [p,q,r] [deg/s]
    real(8),dimension(0:2)::Angular_acceleration !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8),dimension(0:2)::Angular_acceleration_old !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8),dimension(0:2)::Air_angle !迎角，横滑り角，経路角 [α,β,γ] [deg]
    real(8),dimension(0:2)::Input !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2)::Input_old !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2)::Wind !対地風ベクトル [ugE,vgE,wgE] [m/s] @プラホ高度における
    real(8),dimension(0:2)::gust !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    real(8),dimension(0:1)::Speed !対地速度，対気速度 [VE,V] [m/s]
    real(8),dimension(0:2)::Trim_speed !トリム速度ベクトル [u0,v0,w0] [m/s]
    real(8),dimension(0:2)::Trim_angle !トリム角度 [Φ0,θ0,Ψ0] [m/s]
    real(8),dimension(0:2)::Target !目標値．ダイブ角 [deg]，定常滑空角 [deg]，バンク角 [deg]
    real(8),dimension(0:2)::PID_para_lon1,PID_para_lon2 !縦の制御パラメーター(ダイブ),(定常滑空) [-]
    real(8),dimension(0:2)::PID_para_lat1,PID_para_lat2 !横・方向の制御パラメーター(定常滑空),(バンク) [m/s]
    real(8),dimension(0:2)::PID_lon1,PID_lon2 !縦の制御パラメーター(ダイブ)，(定常滑空)
    real(8),dimension(0:2)::PID_lat1,PID_lat2 !横・方向の制御パラメーター(定常滑空),(バンク)
    real(8)::Trim_velocity !トリム速度 V [m/s]
    real(8)::Mass,Weight !機体質量 [kg]，機体重量 [N]
    real(8)::nx,ny,nz !荷重倍数
    !着水判定用
    real(8)::alpha_stall !失速角 [deg]
    real(8)::hE_water !着水高度
    real(8),parameter::limit_load_factor =3.0D0 !制限荷重倍数 [G]
    !操舵の制限
    real(8)::de_max,dh_max,dr_max !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg]
    !評価関数
    real(8)::sum_error,sum_D_input !偏差の総和，入力の時間変化率の総和
    !時間
    real(8)::time !時間 [s]
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    real(8)::time_pullup !引き起こし終了時刻 [s]
    real(8)::time_bank !バンク開始時刻 [s]
    real(8)::time_bank_end !バンク終了時刻 [s]
    !その他
    real(8)::rho,gravity !空気密度 [kg/m^3]，重力加速度 [m/s^2]
    real(8)::phi,theta,psi
    real(8)::p,q,r
    real(8)::u,v,w,ug,vg,wg
    real(8)::xE,yE,hE,hE_0
    real(8)::alpha,beta,gamma
    real(8)::tmp,de,dh,dr
    real(8),dimension(:,:),allocatable::matrix1,matrix2
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    !カウンター
    integer::i
    integer::iteration
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min !initial time, finish time, time rate

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換
    
    call system_clock(initial_time) !開始時間の読み込み

    !配列の初期化
    Lift=0.0D0;         Drag=0.0D0
    Force=0.0D0;        Moment=0.0D0
    Position=0.0D0;     Velocity=0.0D0; Velocity_old=0.0D0 !速度ベクトル [m]
    Altitude_angle=0.0D0;               Altitude_angular_velocity=0.0D0 !姿勢角 [deg]
    Alt_ang_vel_old=0.0D0 !姿勢角 [deg]
    Ground_speed=0.0D0 !機体軸における対地速度ベクトル [m/s]
    Acceleration=0.0D0;                 Acceleration_old=0.0D0
    Angular_velocity=0.0D0;             Angular_acceleration=0.0D0
    Angular_acceleration_old=0.0D0
    Air_angle=0.0D0;    Input=0.0D0;    Input_old=0.0D0 !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg]
    Wind=0.0D0;         gust=0.0D0;     speed=0.0D0
    PID_lon1=0.0D0;     PID_lon2=0.0D0; PID_lat1=0.0D0;     PID_lat2=0.0D0
    FlightLog=0.0D0;    sum_error=0.0D0;sum_D_input=0.0D0;  Distance=0.0D0
    CPU_time_min=0; time=0.0D0

    !---------------------------------------------------!
    !--------------------値の読み込み--------------------!
    !---------------------------------------------------!

    do i=0,2
        Altitude_angle(i)=genom_list(i)
    end do
    time_pullup=genom_list(3) !引き起こし時間 [s]
    time_bank=genom_list(4) !バンク開始時間 [s]
    time_bank_end=genom_list(5) !バンク開始時間 [s]
    do i=0,2
        Target(i)=genom_list(6+i)
        PID_para_lon1(i)=genom_list(9+i)
        PID_para_lon2(i)=genom_list(12+i)
        PID_para_lat1(i)=genom_list(15+i)
        PID_para_lat2(i)=genom_list(18+i)
    end do

    pilot_method=int(NLFS_txt(0,0)) !操縦方式
    time_max=NLFS_txt(1,0) !最大飛行時間 [s]
    time_step=NLFS_txt(2,0) !タイムステップ [s]
    alpha_stall=NLFS_txt(3,0) !失速角 [s]
    hE_water=NLFS_txt(4,0) !着水高度 [s]
    Wind(0:2)=NLFS_txt(5,0:2) !対地風ベクトル [m/s]
    Trim_velocity=NLFS_txt(6,0) !トリム速度V [m/s]
    Trim_speed(0:2)=NLFS_txt(7,0:2)!トリム速度u,v,w
    Trim_angle(0:2)=NLFS_txt(8,0:2)!トリム姿勢角
    !初期条件
    Position(0:2)=NLFS_txt(9,0:2) !位置
    speed(0)=NLFS_txt(10,0) !プラホからの飛び出し速度 [m/s]
    Air_angle(2)=NLFS_txt(11,0) !経路角 [deg] プラホなら -4 [deg]
    !舵角の制限
    de_max=NLFS_txt(12,0)
    dh_max=NLFS_txt(13,0)
    dr_max=NLFS_txt(14,0)
    
    hE_0=Position(2) !プラホ高度 [m]

    !機体データ読み込み
    rho=Plane_txt(0,0); gravity=Plane_txt(0,1) !空気密度，重力加速度
    Mass=Plane_txt(1,0) !機体重量 [kg]
    Weight=Mass*gravity !機体重量を計算 [N]

    !---------------------------------------------------!
    !--------------------初期値の計算--------------------!
    !---------------------------------------------------!

    !Φ，θ，ψ [rad]を格納
    phi=rad*Altitude_angle(0);    theta=rad*altitude_angle(1);    psi=rad*altitude_angle(2)

    !ug,vg,wg [m/s]を計算
    allocate(matrix1(0:2,0:0),matrix2(0:2,0:0))
    do i=0,2
        matrix1(i,0)=wind(i) !対地風ベクトルをmatrix1に格納
    end do
    call Transform_coordinate(phi,theta,psi,matrix1,matrix2) !絶対座標系から機体軸座標に変換する
    do i=0,2
        gust(i)=matrix2(i,0)
    end do
    deallocate(matrix1,matrix2)
    ug=gust(0);     vg=gust(1);     wg=gust(2)

    !u,v,w [m/s]を計算
    Ground_speed(0)=speed(0)*cos(theta-Air_angle(2)*(pi/180))
    Ground_speed(1)=0.0D0
    Ground_speed(2)=speed(0)*sin(theta-Air_angle(2)*(pi/180))
    u=Ground_speed(0);  v=Ground_speed(1);  w=Ground_speed(2);  

    !対気速度VEの計算
    speed(1)=sqrt((u+ug)*(u+ug)+(v+vg)*(v+vg)+(w+wg)*(w+wg))

    !迎角，横滑り角，経路角の初期値の計算
    Air_angle(0)=deg*atan(w/u)
    Air_angle(1)=deg*atan(v/speed(0))
    Air_angle(2)=Altitude_angle(1)-Air_angle(0)

    !---------------------------------------------------!
    !--------------------反復計算開始--------------------!
    !---------------------------------------------------!

    do iteration=0,iteration_max

        u=Ground_speed(0);          w=Ground_speed(2);              v=Ground_speed(1)           !u,v,w [m/s]
        phi=rad*Altitude_angle(0);  theta=rad*altitude_angle(1);    psi=rad*altitude_angle(2)   !Φ，θ，ψ [rad]を格納
        p=Angular_velocity(0);      q=Angular_velocity(1);          r=Angular_velocity(2)       !p,q,r [deg]
        xE=Position(0);             yE=Position(1);                 hE=Position(2)              !位置xE,yE,hE [m]を格納       
        alpha=Air_angle(0);         beta=Air_angle(1);              gamma=air_angle(2)          !迎角α,横滑り角β [deg] を格納       
        de=Input(0);                dh=Input(1);                    dr=Input(2)                 !操舵を格納

        call system_clock(finish_time,time_rate) !終了時間の読み込み
        CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
        if(CPU_time_s>60) then
            CPU_time_min=int(CPU_time_s/60.0D0)
            CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        end if
        

        !-------------------------------------------------!
        !--------------------時間の更新--------------------!
        !-------------------------------------------------!

        time=time_step*iteration

        !-----------------------------------------------------------------!
        !--------------------対地風ベクトルを機体軸に変換--------------------!
        !-----------------------------------------------------------------!
        
        !対地風ベクトルをmatrix1に格納
        allocate(matrix1(0:2,0:0),matrix2(0:2,0:0))
        do i=0,2
            matrix1(i,0)=wind(i)
        end do
        call Transform_coordinate(phi,theta,psi,matrix1,matrix2) !絶対座標系から機体軸座標に変換する
        do i=0,2
            gust(i)=matrix2(i,0)
        end do
        deallocate(matrix1,matrix2)
        !1/6乗則によって地面付近の風を弱くする．
        gust=gust*((hE/hE_0)**(1.0D0/6.0D0))

        !-----------------------------------------------!
        !--------------------PID制御--------------------!
        !-----------------------------------------------!

        !P,I,Dを計算
        !縦
        if(time<time_pullup) then !ダイブ中
            PID_lon1(2)=((Air_angle(2)-Target(0))-PID_lon1(0))/time_step ![deg/s]
            PID_lon1(0)=(Air_angle(2)-Target(0)) ![deg]
            PID_lon1(1)=PID_lon1(1)+PID_lon1(0)*time_step ![deg*s]
        else !定常滑空
            PID_lon2(2)=((Air_angle(2)-Target(1))-PID_lon2(0))/time_step ![deg/s]
            PID_lon2(0)=(Air_angle(2)-Target(1)) ![deg]
            PID_lon2(1)=PID_lon2(1)+PID_lon2(0)*time_step ![deg*s]
        end if
        !横・方向
        if(time<time_bank .or. time>time_bank_end) then !定常滑空(目標バンク角0 [deg])
            PID_lat1(2)=((Altitude_angle(0)-0.0D0)-PID_lat1(0))/time_step ![deg/s]
            PID_lat1(0)=(Altitude_angle(0)-0.0D0) ![deg]
            PID_lat1(1)=PID_lat1(1)+PID_lat1(0)*time_step ![deg*s]
        elseif(time>time_bank .and. time<time_bank_end) then !バンク中
            PID_lat2(2)=((Altitude_angle(0)-Target(2))-PID_lat2(0))/time_step ![deg/s]
            PID_lat2(0)=(Altitude_angle(0)-Target(2)) ![deg]
            PID_lat2(1)=PID_lat2(1)+PID_lat2(0)*time_step ![deg*s]
        end if

        !操舵量を計算．エレベータ舵角も重心移動距離も負のピッチングモーメントを発生させる向きが正．
        !縦
        Input=0.0D0
        if(time<time_pullup) then !ダイブ中
            do i=0,2
                Input(pilot_method)=Input(pilot_method)+PID_para_lon1(i)*PID_lon1(i)
            end do
        else !定常滑空
            do i=0,2
                Input(pilot_method)=Input(pilot_method)+PID_para_lon2(i)*PID_lon2(i)
            end do
        end if
        IF(pilot_method==0) then !エレベータで操縦するなら
            Input(0)=Input(0)*(de_max/dh_max) !PIDの値を調整
        end if
        !横・方向
        if(time<time_bank .or. time>time_bank_end) then !定常滑空(目標バンク角0 [deg])
            do i=0,2
                Input(2)=Input(2)+PID_para_lat1(i)*PID_lat1(i)*(dr_max/dh_max)
            end do
        elseif(time>time_bank .and. time<time_bank_end) then !バンク中
            do i=0,2
                Input(2)=Input(2)+PID_para_lat2(i)*PID_lat2(i)*(dr_max/dh_max)
            end do
        end if
        !緩和係数をかける
        do i=0,2
            Input(i)=Input_old(i)+relaxation_coef*(Input(i)-input_old(i))
        end do 
        !舵角を制限
        if(abs(Input(0))>de_max) then
            Input(0)=(Input(0)/abs(Input(0)))*de_max
        end if
        if(abs(Input(1))>dh_max) then
            Input(1)=(Input(1)/abs(Input(1)))*dh_max
        end if
        if(abs(input(2))>dh_max) then
            Input(2)=(Input(2)/abs(Input(2)))*dr_max
        end if
        !Input_oldの更新はこのループの最後で

        !------------------------------------------------!
        !--------------------全機計算--------------------!
        !------------------------------------------------!

        call Plane_calculation &
            (Altitude_angle,Ground_speed,Angular_velocity,gust,Input &
            ,hE,Lift,Drag,Force,Moment,Acceleration,Angular_acceleration,Plane_txt)

        !---------------------------------------------------!
        !--------------------加速度を積分--------------------!
        !---------------------------------------------------!

        !Adams-Bashforth法で積分
        Ground_speed=Ground_speed+((3.0D0*Acceleration-Acceleration_old)/2.0D0)*time_step
        Angular_velocity=Angular_velocity+((3.0D0*Angular_acceleration-Angular_acceleration_old)/2.0D0)*time_step

        !oldを更新
        Acceleration_old=Acceleration
        Angular_acceleration_old=Angular_acceleration

        !---------------------------------------------------!
        !--------------------角速度を計算--------------------!
        !---------------------------------------------------!
        
        p=Angular_velocity(0)
        q=Angular_velocity(1)
        r=Angular_velocity(2)
        
        !機体軸p,q,rから慣性系dΦ/dt,dθ/dt,dΨ/dtへ
        Altitude_angular_velocity(0)=p+(r*cos(phi)+q*sin(phi))*tan(theta)
        Altitude_angular_velocity(1)=q*cos(phi)-r*sin(phi)
        Altitude_angular_velocity(2)=(r*cos(phi)+q*sin(phi))/cos(theta)

        !-------------------------------------------------!
        !--------------------速度を計算--------------------!
        !-------------------------------------------------!
        
        !速度を計算
        ug=gust(0)
        vg=gust(1)
        wg=gust(2)
        u=Ground_speed(0)
        w=Ground_speed(2)
        v=Ground_speed(1)
        speed(0)=sqrt(u*u+v*v+w*w)
        speed(1)=sqrt((u+ug)**2+(v+vg)**2+(w+wg)**2)

        !機体軸u,v,wから慣性系dxE/dt,dyE/dt,dhE/dtへ
        allocate(matrix1(0:2,0:0),matrix2(0:2,0:0))
        do i=0,2
            matrix1(i,0)=Ground_speed(i)
        end do
        call Inverse_Transform_coordinate(phi,theta,psi,matrix1,matrix2) !機体座標系から絶対軸座標に変換する
        do i=0,2
            Velocity(i)=matrix2(i,0)
        end do
        Velocity(2)=-Velocity(2) !dhE/dtは上向きが正
        deallocate(matrix1,matrix2)

        !-----------------------------------------------------------------!
        !--------------------位置，姿勢角,荷重倍数を計算--------------------!
        !-----------------------------------------------------------------!

        !Adams-Bashforth法で積分
        Altitude_angle=Altitude_angle+((3.0D0*Altitude_angular_velocity-Alt_ang_vel_old)/2.0D0)*time_step
        Position=Position+((3.0D0*Velocity-Velocity_old)/2.0D0)*time_step
        
        !dx/dt_oldを更新
        Alt_ang_vel_old=Altitude_angular_velocity
        Velocity_old=Velocity
        
        !迎角，横滑り角，経路角の初期値の計算
        Air_angle(0)=deg*atan((w+wg)/(u+ug))
        Air_angle(1)=deg*atan((v+vg)/speed(0))
        Air_angle(2)=Altitude_angle(1)-deg*atan(w/u) !対地．gustを含めない

        !荷重倍数の計算
        nx=Force(0)/Weight
        ny=Force(1)/Weight
        nz=-Force(2)/Weight

        !-----------------------------------------------!
        !--------------------着水判定--------------------!
        !-----------------------------------------------!

        tmp=0.0D0
        if(iteration==iteration_max) then
            tmp=tmp+1.0D0
        elseif(alpha>alpha_stall) then !Stalled
            tmp=tmp+1.0D0
        elseif(position(2)<hE_water) then !Landing
            tmp=tmp+1.0D0
        elseif(Velocity(0)<0.0D0) then !Return Flight"
            tmp=tmp+1.0D0
        end if
        if(tmp>0.0D0) then
            iteration_max=iteration
            exit
        end if

        !偏差の総和を計算する
        sum_error=sum_error+(abs(PID_lon1(0))+abs(PID_lon2(0))+(abs(PID_lat1(0))+abs(PID_lat2(0)))*weighting_factor2)*time_step

        !入力の時間変化の総和を計算する
        sum_D_input=sum_D_input+(abs(Input(0)-Input_old(0))/de_max)
        sum_D_input=sum_D_input+(abs(Input(1)-Input_old(1))/dh_max)
        sum_D_input=sum_D_input+(abs(Input(2)-Input_old(2))/dr_max)

        !Input_oldを更新
        Input_old=Input

        !--------------------------------------------------------------!
        !--------------------FlightLogにデータを入力--------------------!
        !--------------------------------------------------------------!

        FlightLog(iteration,0)=time
        do i=0,2
            FlightLog(iteration,1+i)=Position(i) !位置ベクトル [m]
            FlightLog(iteration,4+i)=Velocity(i) !速度ベクトル [m/s]
            FlightLog(iteration,7+i)=Altitude_angle(i) !姿勢角 [deg]
            FlightLog(iteration,10+i)=Altitude_angular_velocity(i) !姿勢角速度 [deg]
            FlightLog(iteration,13+i)=Ground_speed(i) !機体軸における対地速度ベクトル [m/s]
            FlightLog(iteration,16+i)=Angular_velocity(i) !機体軸における角速度ベクトル [deg/s]
            FlightLog(iteration,19+i)=Air_angle(i) !迎角，横滑り角，経路角 [deg]
            FlightLog(iteration,22+i)=Input(i) !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg]
            FlightLog(iteration,25+i)=Wind(i) !対地風ベクトル [m/s]
            FlightLog(iteration,28+i)=gust(i) !突風ベクトル [m/s]
        end do
        do i=0,1
            FlightLog(iteration,31+i)=speed(i) !対地速度，対気速度 [m/s]
        end do
        FlightLog(iteration,33)=Lift !揚力 [N]
        FlightLog(iteration,34)=Drag !抗力 [N]
        do i=0,2
            FlightLog(iteration,35+i)=Force(i) !機体にはたらく外力 [X,Y,Z] [N]
            FlightLog(iteration,38+i)=Moment(i) !機体にはたらくモーメント [L,M,N] [N*m]
        end do
        FlightLog(iteration,41)=nx !荷重倍数 [-]
        FlightLog(iteration,42)=ny !荷重倍数 [-]
        FlightLog(iteration,43)=nz !荷重倍数 [-]
        FlightLog(iteration,44)=Lift/Drag !揚抗比
        FlightLog(iteration,45)=sqrt(ug*ug+vg*vg+wg*wg) !現在の高度での風速
        if(time<time_pullup) then !ダイブ中
            FlightLog(iteration,46)=Target(0) !目標経路角
        else !定常滑空
            FlightLog(iteration,46)=Target(1) !目標経路角
        end if
        if(time<time_bank .or. time>time_bank_end) then !定常滑空(目標バンク角0 [deg])
            FlightLog(iteration,47)=0 !目標バンク角
        elseif(time>time_bank .and. time<time_bank_end) then !バンク中
            FlightLog(iteration,47)=Target(2) !目標バンク角
        end if

    end do

    Distance=sqrt(Position(0)**2+Position(1)**2)
    !飛距離から偏差の総和を重みづけして差をとる
    Distance=Distance-sum_error*weighting_factor1-sum_D_input*weighting_factor3
    if(isnan(Distance)) Distance=0.0D0

end subroutine

subroutine Plane_calculation &
    (Altitude_angle,Ground_speed,Angular_velocity,gust,Input &
    ,hE,Lift,Drag,Force,Moment,Acceleration,Angular_acceleration,Plane_txt)
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan,sqrt,sin,cos,tan !プログラム中で使用する関数の名前を宣言
    !引数
    real(8),dimension(0:2),intent(IN)::Altitude_angle !姿勢角 [Φ,θ,Ψ] [deg]
    real(8),dimension(0:2),intent(INOUT)::Ground_speed !機体軸における対地速度ベクトル [u,v,w] [m/s]
    real(8),dimension(0:2),intent(INOUT)::Angular_velocity !機体軸における角速度ベクトル [p,q,r] [deg/s]
    real(8),dimension(0:2),intent(IN)::Input !エレベータ舵角 [deg]，重心移動距離 [-]，ラダー舵角 [deg] [δe,δh,δr]
    real(8),dimension(0:2),intent(IN)::gust !機体軸における突風ベクトル [ug,vg,wg] [m/s]
    real(8),dimension(0:2),intent(OUT)::Acceleration !機体軸における対地加速度ベクトル [du/dt,dv/dt,dw/dt] [m/s]
    real(8),dimension(0:2),intent(OUT)::Angular_acceleration !機体軸における角加速度ベクトル [dp/dt,dq/dt,dr/dt] [deg/s]
    real(8)::hE !高度 [m]
    real(8),intent(OUT)::Lift !揚力 [N]
    real(8),intent(OUT)::Drag !抗力 [N]
    real(8),dimension(0:2),intent(OUT)::Force !機体にはたらく外力 [X,Y,Z] [N]
    real(8),dimension(0:2),intent(OUT)::Moment !機体にはたらくモーメント [L,M,N] [N*m]
    !機体データ
    real(8),dimension(0:11,0:100),intent(IN)::Plane_txt !NLFS.txtの中身
    real(8)::rho, gravity!空気密度，重力加速度
    real(8)::Mass !機体重量 [kg]
    real(8)::S, b, c_mac !主翼面積，スパン，平均空力翼弦長
    real(8)::Ixx,Iyy,Izz,Ixz !慣性モーメント [kg*m^4]
    real(8),dimension(0:100)::Cx_coef, Cz_coef, Cm_coef !Cx，Cz，Cmを多項式近似したときの係数
    real(8)::Cx, Cy, Cz, Cl, Cm, Cn
    real(8)::Cyb, Cydr, Clb, Clp, Clr, Cldr, Cmq, Cma_, Cmde, Cmdh, Cnb, Cnp, Cnr, Cndr
    !その他
    real(8)::Velocity
    real(8)::X,Y,Z,L,M,N
    real(8)::phi,theta,psi
    real(8)::de,dh,dr
    real(8)::p,q,r
    real(8)::u,v,w,ug,vg,wg
    real(8)::alpha,beta,gamma
    real(8)::pi,deg,rad !円周率，degへの変換，radへの変換
    !カウンター
    integer::i,j,k,num
    integer::i_max,j_max,k_max

    pi=4.0D0*atan(1.0D0) !円周率
    deg=(180.0D0/pi) !degへの変換
    rad=(pi/180.0D0) !radへの変換

    Lift=0.0D0 !揚力 [N]
    Drag=0.0D0 !抗力 [N]
    Force=0.0D0 !機体にはたらく外力 [X,Y,Z] [N]
    Moment=0.0D0 !機体にはたらくモーメント [L,M,N] [N*m]
    Acceleration=0.0D0
    Angular_acceleration=0.0D0

    !機体データ読み込み
    open(10,file="Plane.txt") !ファイルを開く
        rho=Plane_txt(0,0); gravity=Plane_txt(0,1) !空気密度，重力加速度
        Mass=Plane_txt(1,0) !機体重量 [kg]
        S=Plane_txt(2,0);b=Plane_txt(2,1);c_mac=Plane_txt(2,2) !主翼面積，スパン，平均空力翼弦長
        Ixx=Plane_txt(3,0);Iyy=Plane_txt(3,1);Izz=Plane_txt(3,2);Ixz=Plane_txt(3,3) !慣性モーメント [kg*m^4]
        i_max=int(Plane_txt(4,0));j_max=int(Plane_txt(4,1));k_max=int(Plane_txt(4,2));
        Cx_coef(0:100)=Plane_txt(5,0:100)
        Cz_coef(0:100)=Plane_txt(6,0:100)
        Cm_coef(0:100)=Plane_txt(7,0:100)
        Cyb=Plane_txt(8,0); Cydr=Plane_txt(8,1)
        Clb=Plane_txt(9,0); Clp=Plane_txt(9,1); Clr=Plane_txt(9,2); Cldr=Plane_txt(9,3)
        Cmq=Plane_txt(10,0); Cma_=Plane_txt(10,1); Cmde=Plane_txt(10,2); Cmdh=Plane_txt(10,3)
        Cnb=Plane_txt(11,0); Cnp=Plane_txt(11,1); Cnr=Plane_txt(11,2); Cndr=Plane_txt(11,3)
    close(10) !ファイルを閉じる

    !わかりやすいように変数変換
    u=Ground_speed(0);          w=Ground_speed(2);              v=Ground_speed(1)           !u,v,w [m/s]
    ug=gust(0);                 wg=gust(2);                     vg=gust(1)                  !u,v,w [m/s]
    phi=Altitude_angle(0);      theta=Altitude_angle(1);        psi=Altitude_angle(2)   !Φ，θ，ψ [rad]を格納
    p=Angular_velocity(0);      q=Angular_velocity(1);          r=Angular_velocity(2)       !p,q,r [deg]       
    de=Input(0);                dh=Input(1);                    dr=Input(2)                 !操舵を格納

    !速度V [m/s]，全機迎角α [deg]，横滑り角β [deg]を計算する
    Velocity=sqrt((u+ug)**2+(v+vg)**2+(w+wg)**2)
    alpha=deg*asin((w+wg)/(u+ug))
    beta=deg*asin((v+vg)/Velocity)
    gamma=theta-alpha

    !空力係数を計算
    Cx=0;Cz=0;Cm=0
    num=0
    do i=0,i_max
        do j=0,j_max
            do k=0,k_max
                Cx=Cx+Cx_coef(num)*(Velocity**i)*(alpha**j)*(hE**k)
                Cz=Cz+Cz_coef(num)*(Velocity**i)*(alpha**j)*(hE**k)
                Cm=Cm+Cm_coef(num)*(Velocity**i)*(alpha**j)*(hE**k)
                num=num+1
            end do
        end do
    end do
    Cy=Cyb*beta+Cydr*dr
    Cl=Clb*beta+(b/(2*Velocity))*(Clp*p+Clr*r)*rad+Cldr*dr
    Cm=Cm+(c_mac/(2*Velocity))*(Cmq*q)*rad+Cmde*de+Cmdh*dh
    Cn=Cnb*beta+(b/(2*Velocity))*(Cnp*p+Cnr*r)*rad+Cndr*dr

    !機体にはたらく力，モーメントを計算
    Force(0)=(1.0D0/2.0D0)*rho*Velocity*Velocity*S*Cx
    Force(1)=(1.0D0/2.0D0)*rho*Velocity*Velocity*S*Cy
    Force(2)=(1.0D0/2.0D0)*rho*Velocity*Velocity*S*Cz
    Moment(0)=(1.0D0/2.0D0)*rho*Velocity*Velocity*S*b*Cl
    Moment(1)=(1.0D0/2.0D0)*rho*Velocity*Velocity*S*c_mac*Cm
    Moment(2)=(1.0D0/2.0D0)*rho*Velocity*Velocity*S*b*Cn

    X=Force(0)
    Y=Force(1)
    Z=Force(2)
    L=(Iyy-Izz)*rad*q*r+Ixz*rad*p*q+deg*Moment(0)
    M=(Izz-Ixx)*rad*r*p+Ixz*rad*(r*r-p*p)+deg*Moment(1)
    N=(Ixx-Iyy)*rad*p*q-Ixz*rad*q*r+deg*Moment(2)
    Lift=X*sin(rad*alpha)-Z*cos(rad*alpha)
    Drag=-X*cos(rad*alpha)-Z*sin(rad*alpha)
    
    !非線形6自由度運動方程式を解く．青本p.16
    !並進運動方程式
    Acceleration(0)=-rad*q*w+rad*r*v-gravity*sin(rad*theta)+(X/Mass)
    Acceleration(1)=-rad*r*u+rad*p*w+gravity*cos(rad*theta)*sin(rad*phi)+(Y/Mass)
    Acceleration(2)=-rad*p*v+rad*q*u+gravity*cos(rad*theta)*cos(rad*phi)+(Z/Mass)
    !回転運動方程式
    Angular_acceleration(0)=((L/Ixx)+(Ixz/Ixx)*(N/Izz))/(1.0D0-(Ixz**2)/(Izz*Ixx))
    Angular_acceleration(1)=M/Iyy
    Angular_acceleration(2)=((N/Izz)+(Ixz/Izz)*(L/Ixx))/(1.0D0-(Ixz**2)/(Ixx*Izz))

end subroutine

subroutine Output_FlightLog(i_max_copy,iteration_max,FlightLog_data,FlightLog,genom_length,genom_list)
    !フライトログをcsvファイルで出力するための脳筋subroutine
    implicit none !暗黙の変数宣言を無効にする
    integer,intent(IN)::i_max_copy !実際に反復した回数
    integer,intent(IN)::iteration_max !最大反復回数
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    real(8),dimension(0:iteration_max,0:FlightLog_data),intent(IN)::FlightLog !飛行解析のデータ
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1),intent(IN)::genom_list !遺伝子情報
    integer::i,j

    open (10, file='FLIGHTLOG.txt',status='replace')
        do i = 0, i_max_copy
            do j=0,FlightLog_data
                write(10,'(F0.8,X)',advance='NO') FlightLog(i,j)
            end do
            write(10,*)
        end do
    close (10)

    open (10, file='LIST.txt', status='replace')
        write(10,'(F0.8,X,F0.8,X,F0.8)') genom_list(0),genom_list(1),genom_list(2)
        do i=3,8
            write(10,'(F0.8)') genom_list(i) !1行ずつ値を読み込む
        end do
        do i=0,3
            write(10,'(F0.8,X,F0.8,X,F0.8)') genom_list(9+i*3),genom_list(10+i*3),genom_list(11+i*3)
        end do
    close(10) !ファイルを閉じる

end subroutine

subroutine GA_NLFS
    !遺伝的アルゴリズムによって飛行経路を最適化する
    !符号:実数値エンコーディング
    !選択:エリート主義
    !交叉:二点交叉(or 一様交叉)
    !変位:摂動
    !世代:定常状態モデル
    implicit none !暗黙の変数宣言を無効にする
    integer,parameter::genom_length=21 !遺伝子情報の長さ
    integer,parameter::max_genom_list=60 !遺伝子集団の大きさ
    integer,parameter::select_genom=20 !遺伝子選択数
    integer,parameter::max_generation=90 !繰り返す世代
    real(8),parameter::individual_mutation=0.020D0 !個体突然変異率
    real(8),parameter::genom_mutation=0.050D0 !遺伝子突然変異率
    integer,parameter::FlightLog_data=47 !飛行解析のデータ数
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8),dimension(0:genom_length-1,0:max_genom_list-1)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1,0:select_genom-1)::elite_class !エリート
    real(8),dimension(0:genom_length-1,0:select_genom-1)::progeny_class !子孫
    real(8),dimension(0:max_genom_list-1)::Distance !飛行距離 [m]
    real(8),dimension(0:max_generation)::Distance_Max !各世代の最大飛行距離 [m]
    !カウンター
    integer::i,j,generation
    !計算時間計測
    real(8)::CPU_time_s
    integer::initial_time,finish_time,time_rate,CPU_time_min,CPU_time_h !initial time, finish time, time rate
    
    call system_clock(initial_time) !開始時間の読み込み

    !------------------------------------------------------!
    !--------------------遺伝子集団を作る--------------------!
    !------------------------------------------------------!

    do i=0,max_genom_list-1
        call create_genom(genom_length,genom_list)
        do j=0,genom_length-1
            genom_class(j,i)=genom_list(j)
        end do
    end do
    
    !------------------------------------------------------!
    !--------------------遺伝子を評価する--------------------!
    !------------------------------------------------------!

    call evalution(genom_length,max_genom_list,genom_class,Distance,FlightLog_data)

    !-------------------------------------------------------!
    !--------------------遺伝子の並び替え--------------------!
    !-------------------------------------------------------!

    call sort_genom(max_genom_list,genom_length,Distance,genom_class)

    !---------------------------------------------------!
    !--------------------世代のループ--------------------!
    !---------------------------------------------------!

    do generation=0,max_generation

        write(*,*) "------------------------------------------------------------"
        write(*,'(A24,I12,A12,I12)') "Generaion",generation,"/",max_generation 
        call system_clock(finish_time,time_rate) !終了時間の読み込み
        CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
        CPU_time_min=0.0D0
        if(CPU_time_s>60) then
            CPU_time_min=int(CPU_time_s/60.0D0)
            CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        end if
        write(*,'(A24,I12,A12,F12.3,A12)') "Computation Time",CPU_time_min,"[min]",CPU_time_s,"[s]"
        write(*,*) "------------------------------------------------------------"
        !if(CPU_time_min>=95) exit !97分で計算打ち切り．CPUの都合上

        !-----------------------------------------------------!
        !--------------------エリートの選択--------------------!
        !-----------------------------------------------------!

        call select_elite(select_genom,max_genom_list,genom_length,genom_class,elite_class)

        !-------------------------------------------!
        !--------------------交叉--------------------!
        !-------------------------------------------!

        call cross_over(select_genom,genom_length,max_genom_list,genom_class,elite_class,progeny_class)
        
        !-----------------------------------------------!
        !--------------------突然変異--------------------!
        !-----------------------------------------------!

        call mutation(individual_mutation,genom_mutation,genom_length,max_genom_list,genom_class)

        !------------------------------------------------------!
        !--------------------遺伝子を評価する--------------------!
        !------------------------------------------------------!

        call evalution(genom_length,max_genom_list,genom_class,Distance,FlightLog_data)

        !-------------------------------------------------------!
        !--------------------遺伝子の並び替え--------------------!
        !-------------------------------------------------------!

        call sort_genom(max_genom_list,genom_length,Distance,genom_class)

        !-------------------------------------------------!
        !--------------------結果を表示--------------------!
        !-------------------------------------------------!

        if(mod(generation,1)==0 .OR. generation==max_generation)  then
            call output_genom(genom_length,max_genom_list,genom_class,Distance,generation,FlightLog_data)
        end if
        Distance_Max(generation)=Distance(0)

    end do

    write(*,*) "------------------------------------------------------------"
    call system_clock(finish_time,time_rate) !終了時間の読み込み
    CPU_time_s=(finish_time-initial_time)/dble(time_rate) ![s]
    if(CPU_time_s>60) then
        CPU_time_min=int(CPU_time_s/60.0D0)
        CPU_time_s=CPU_time_s-60.D0*CPU_time_min
        if(CPU_time_min>60) then
            CPU_time_h=int(CPU_time_min/60)
            CPU_time_min=CPU_time_min-60*CPU_time_h
        end if
    end if
    write(*,'(A24,I12,A12,I12,A12,F12.3,A12)') "Computation Time",CPU_time_h,"[h]",CPU_time_min,"[min]",CPU_time_s,"[s]"
    write(*,*) "------------------------------------------------------------"
    !call write1D(max_generation,Distance_Max)
    
end subroutine

subroutine create_genom(genom_length,genom_list)
    !ランダムに遺伝子を生成する
    implicit none !暗黙の変数宣言を無効にする
    !引数
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    real(8),dimension(0:genom_length-1),intent(OUT)::genom_list !遺伝子情報
    !変数
    real(8),parameter::phi_max =0.0D0 !最大初期バンク角 [deg]
    real(8),parameter::theta_max =4.0D0 !最大初期ピッチ角 [deg]
    real(8),parameter::psi_max =0.0D0 !最大初期ヨー角 [deg]
    real(8)::time_max !最大飛行時間 [s]
    real(8),parameter::max_dive_angle =-20.0D0 !最大ダイブ角 [deg]
    real(8),parameter::glide_angle =-2.000D0 !定常滑空角 [deg]
    real(8),parameter::max_bank_angle =10.000D0 !最大バンク角 [deg]
    real(8)::rnd !0~1の乱数
    real(8)::tmp
    !カウンター
    integer::i

    !操縦方式，飛行時間，タイムステップの読み込み
    open(10,file="NLFS.txt") !ファイルを開く
    read(10,*)  !操縦方式
    read(10,*) time_max !最大飛行時間
    read(10,*)  !タイムステップ
    close(10) !ファイルを閉じる
    
    call random(rnd)
    genom_list(0)=phi_max*rnd !初期バンク角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のバンク角にする
        genom_list(0)=-genom_list(0)
    end if

    call random(rnd)
    genom_list(1)=theta_max*rnd !初期ピッチ角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のピッチ角にする
        genom_list(1)=-genom_list(1)
    end if

    call random(rnd)
    genom_list(2)=-psi_max*rnd !初期ヨー角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のヨー角にする
        genom_list(2)=-genom_list(2)
    end if
    genom_list(2)=psi_max !初期ヨー角 [deg]

    call random(rnd)
    genom_list(3)=time_max*rnd*0.25D0 !引き起こし時刻 [s]

    call random(rnd)
    genom_list(4)=time_max*rnd*0.5D0 !バンク開始時刻 [s]

    call random(rnd)
    genom_list(5)=time_max*rnd !バンク終了時刻 [s]
    if(genom_list(5)<genom_list(4)) then !開始時刻＜終了時刻となるように並び替え
        tmp=genom_list(4)
        genom_list(4)=genom_list(5)
        genom_list(5)=tmp
    end if

    call random(rnd)
    genom_list(6)=glide_angle+max_dive_angle*rnd !ダイブ角 [deg]
    !genom_list(6)=max_dive_angle*rnd !ダイブ角 [deg]

    call random(rnd)
    !genom_list(7)=glide_angle !定常滑空角 [deg]
    !genom_list(7)=max_dive_angle*rnd !定常滑空角 [deg]
    genom_list(7)=glide_angle*rnd !定常滑空角 [deg]

    call random(rnd)
    genom_list(8)=max_bank_angle*rnd !バンク角 [deg]
    call random(rnd)
    if(rnd<0.5D0) then !半分の確率で負のバンク角にする
        genom_list(8)=-genom_list(8)
    end if

    do i=0,11
        call random(rnd)
        genom_list(9+i)=rnd !PIDパラメータ
    end do

end subroutine

subroutine evalution(genom_length,max_genom_list,genom_class,Distance,FlightLog_data)
    !各遺伝子を評価する．
    implicit none !暗黙の変数宣言を無効にする
    intrinsic atan !プログラム中で使用する関数の名前を宣言
    !引数
    integer,intent(IN)::genom_length !遺伝子情報の長さ
    integer,intent(IN)::max_genom_list !遺伝子集団の大きさ
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(IN)::genom_class !遺伝子集団
    real(8),dimension(0:max_genom_list-1),intent(OUT)::Distance !飛行距離 [m]
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    !変数
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8),dimension(0:14,0:2)::NLFS_txt !NLFS.txtの中身
    real(8),dimension(0:11,0:100)::Plane_txt !NLFS.txtの中身
    real(8),dimension(:,:),allocatable::FlightLog !飛行解析のデータ
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    integer::iteration_max !最大反復回数 [s]
    integer::i_max_copy !最大反復回数のコピー [s]
    integer::pilot_method !操縦方式．0：エレベータ，1：重心移動
    !カウンター
    integer::i,j,num


    !解析条件読み込み
    open(10,file="NLFS.txt") !ファイルを開く
        do i=0, 14
            !write(*,*) i
            read(10,*) NLFS_txt(i,0:2)
        end do
        pilot_method=int(NLFS_txt(0,0)) !操縦方式
        time_max=NLFS_txt(1,0) !最大飛行時間 [s]
        time_step=NLFS_txt(2,0) !タイムステップ [s]
    close(10) !ファイルを閉じる
    !機体データ読み込み
    Plane_txt=0.0D0
    open(10,file="Plane.txt") !ファイルを開く
        do i=0, 11
            if(i<5 .or. 7<i) then
                read(10,*) Plane_txt(i,0:16)
            else
                num=int((Plane_txt(4,0)+1)*(Plane_txt(4,1)+1)*(Plane_txt(4,2)+1))
                read(10,*) Plane_txt(i,0:num-1)
            end if
        end do
    close(10) !ファイルを閉じる

    iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
    allocate(FlightLog(0:iteration_max,0:FlightLog_data))
    
    num=0
    iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
    do i=0,max_genom_list-1

        !遺伝子情報の読み込み
        do j=0,genom_length-1
            genom_list(j)=genom_class(j,i) !1行ずつ値を読み込む
        end do

        i_max_copy=iteration_max
        call NLFS(Distance(i),genom_length,genom_list,i_max_copy,pilot_method,FlightLog_data,FlightLog,NLFS_txt,Plane_txt)
        !write(*,'(A24,I12,A12,I12,A12,F12.5)') "Genom",num,"/",max_genom_list,"Distance",Distance(i)
        num=num+1
        if(isnan(Distance(i))) Distance(i)=0.0D0
    
    end do

end subroutine

subroutine sort_genom(max_genom_list,genom_length,Distance,genom_class)
    !ゲノムを優秀な順番に並び替える
    implicit none
    integer,intent(IN):: max_genom_list
    integer,intent(IN):: genom_length
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:max_genom_list-1),intent(INOUT):: Distance
    integer::first,last

    first=0;last=max_genom_list-1
    call quicksort1(Distance,genom_class,first,last,max_genom_list-1,genom_length-1)

end subroutine

recursive subroutine quicksort1(a,b,first,last,n,m)
    !
    ! Original
    ! https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !
    implicit none
    integer,intent(IN)::first,last,n,m
    real(8),dimension(0:n),intent(INOUT)::a
    real(8),dimension(0:m,0:n),intent(INOUT)::b
    real(8),dimension(0:m)::c
    real(8)::x,t
    integer::i,j

    x = a((first+last)/2) !基準を設定
    i = first !境界の始まり
    j = last !境界の終わり
    do
        do while (a(i) > x) !基準より小さな値を最初から探すループ
        i=i+1 !基準より小さな値の要素番号
        end do
        do while (x > a(j)) !基準より大きな値を最後から探すループ
        j=j-1 !基準より大きな値の要素番号
        end do
        if (i >= j) exit !基準より小さな値が大きな値より最後側にあれば終了
        !基準より小さな値と大きな値を入れ替える
        t=a(i);    c(0:m)=b(0:m,i)
        a(i)=a(j); b(0:m,i)=b(0:m,j)
        a(j)=t;    b(0:m,j)=c(0:m)
        i=i+1
        j=j-1
    end do
    if (first < i-1) call quicksort1(a,b,first,i-1,n,m)
    if (j+1 < last)  call quicksort1(a,b,j+1,last,n,m)

    return
end subroutine quicksort1

subroutine select_elite(sselect_genom,max_genom_list,genom_length,genom_class,elite_class)
    !エリートを選択する
    implicit none
    integer,intent(IN)::sselect_genom !遺伝子選択数
    integer,intent(IN):: max_genom_list
    integer,intent(IN):: genom_length
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1,0:sselect_genom-1),intent(OUT)::elite_class !エリート
    integer::i,j

    do i=0,sselect_genom-1
        do j=0,genom_length-1
            elite_class(j,i)=genom_class(j,i)
        end do
    end do

end subroutine

subroutine cross_over(sselect_genom,genom_length,max_genom_list,genom_class,elite_class,progeny_class)
    !ランダムに2つのエリートを選択し，一様交叉を行う
    implicit none
    integer,intent(IN)::sselect_genom !遺伝子選択数
    integer,intent(IN)::genom_length
    integer,intent(IN)::max_genom_list
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1,0:sselect_genom-1),intent(IN)::elite_class !エリート
    real(8),dimension(0:genom_length-1,0:sselect_genom-1),intent(OUT)::progeny_class !子孫
    real(8),dimension(0:genom_length-1)::genom_list1,genom_list2 !遺伝子情報
    integer::tmp,genom1,genom2
    real(8)::rnd
    integer::i,j

    !select_genomの数だけ子孫を作る
    do i=0,sselect_genom-1
        !ランダムにエリートを2つ選択する
        call random(rnd)
        do j=0,genom_length-1
            genom_list1(j)=elite_class(j,int(sselect_genom*rnd))
        end do
        call random(rnd)
        do j=0,genom_length-1
            genom_list2(j)=elite_class(j,int(sselect_genom*rnd))
        end do
        !一様交叉を行う
        !do j=0,genom_length-1
        !    call random(rnd)
        !    If(rnd<0.5D0) then
        !        progeny_class(j,i)=genom_list1(j)
        !    else
        !        progeny_class(j,i)=genom_list2(j)
        !    end if
        !二点交叉を行う
        !ランダムに2つの遺伝子を選択する
        call random(rnd)
        genom1=int(genom_length*rnd)
        call random(rnd)
        genom2=int(genom_length*rnd)
        if(genom1>genom2) then !genom1<genom2になるよう並び替え
            tmp=genom1
            genom1=genom2
            genom2=tmp
        end if
        do j=0,genom_length-1 !まずgenom_list1を格納し
            progeny_class(j,i)=genom_list1(j)
        end do
        do j=genom1,genom2 !指定した範囲だけgenom_list2と入れ替える
            progeny_class(j,i)=genom_list2(j)
        end do 
    end do

    !現在の世代の評価の低い個体と子孫を入れ替える
    do i=0,sselect_genom-1
        do j=0,genom_length-1
            genom_class(j,max_genom_list-1-i)=progeny_class(j,i)
        end do
    end do

end subroutine

subroutine mutation(individual_mutation,genom_mutation,genom_length,max_genom_list,genom_class)
    !突然変異を行う
    implicit none
    real(8),intent(IN)::individual_mutation !固体突然変異率
    real(8),intent(IN)::genom_mutation !遺伝子突然変異率
    integer,intent(IN)::genom_length
    integer,intent(IN)::max_genom_list
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(INOUT)::genom_class !遺伝子集団
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8)::time_max !最大飛行時間 [s]
    real(8),parameter::max_dive_angle =-60.0D0 !最大ダイブ角 [deg]
    real(8),parameter::max_glide_angle =-2.0D0 !最大定常滑空角 [deg]
    real(8),parameter::max_bank_angle =10.0D0 !最大バンク角 [deg]
    real(8)::rnd
    real(8)::tmp
    integer::i,j

    !操縦方式，飛行時間，タイムステップの読み込み
    open(10,file="NLFS.txt") !ファイルを開く
    read(10,*)  !操縦方式
    read(10,*) time_max !最大飛行時間
    read(10,*)  !タイムステップ
    close(10) !ファイルを閉じる

    !個体のループ
    do i=0,max_genom_list-1
        call random(rnd)
        if(rnd<individual_mutation) then !個体突然変異
            call create_genom(genom_length,genom_list)
            do j=0,genom_length-1
                genom_class(j,i)=genom_list(j)
            end do
        end if
        
        !遺伝子のループ
        call create_genom(genom_length,genom_list)
        do j=0,genom_length-1
            call random(rnd)
            if(rnd<genom_mutation) then !遺伝子突然変異率
                call random(tmp)
                call random(rnd)
                rnd=(rnd-0.5D0)/abs(rnd-0.5D0) !1.0か-1.0を作る
                genom_class(j,i)=genom_class(j,i)*(1.0D0+rnd*tmp) !遺伝子を0.0~2.0倍にする
            end if
        end do
        
        !バンク開始時刻＜終了時刻となるように並び替え
        if(genom_class(5,i)<genom_class(4,i)) then 
            tmp=genom_Class(4,i)
            genom_class(4,i)=genom_class(4,i)
            genom_class(5,i)=tmp
        end if
        
        genom_class(2,i)=genom_list(2) !Ψを初期値に固定
    end do

end subroutine

subroutine output_genom(genom_length,max_genom_list,genom_class,Distance,generation,FlightLog_data)
    !結果の出力
    implicit none
    integer,intent(IN)::genom_length
    integer,intent(IN)::max_genom_list
    real(8),dimension(0:genom_length-1,0:max_genom_list-1),intent(IN)::genom_class !遺伝子集団
    real(8),dimension(0:max_genom_list-1),intent(IN)::Distance !飛行距離 [m]
    integer,intent(IN)::generation
    integer,intent(IN)::FlightLog_data !FLightLogのデータ数
    real(8),dimension(0:genom_length-1)::genom_list !遺伝子情報
    real(8),dimension(0:14,0:2)::NLFS_txt !NLFS.txtの中身
    real(8),dimension(0:11,0:100)::Plane_txt !NLFS.txtの中身
    real(8)::Distance_Max,Distance_Min,Distance_Ave
    real(8),dimension(:,:),allocatable::FlightLog !飛行解析のデータ
    real(8)::time_max !最大飛行時間 [s]
    real(8)::time_step !タイムステップ [s]
    integer::iteration_max !最大反復回数 [s]
    integer::i_max_copy !最大反復回数のコピー [s]
    integer::pilot_method !操縦方式．0：エレベータ，1：重心移動
    real(8)::tmp
    integer::i,num

    Distance_Max=Distance(0)
    Distance_Min=Distance(max_genom_list-1)
    do i=0,max_genom_list-1
        Distance_Ave=Distance_Ave+Distance(i)
    end do
    Distance_Ave=Distance_Ave/max_genom_list
    do i=0,genom_length-1
        genom_list(i)=genom_class(i,0)
    end do

    open(10,file="NLFS.txt") !ファイルを開く
        do i=0, 14
            read(10,*) NLFS_txt(i,0:2)
        end do
        pilot_method=int(NLFS_txt(0,0)) !操縦方式
        time_max=NLFS_txt(1,0) !最大飛行時間 [s]
        time_step=NLFS_txt(2,0) !タイムステップ [s]
    close(10) !ファイルを閉じる
    !機体データ読み込み
    Plane_txt=0.0D0
    open(10,file="Plane.txt") !ファイルを開く
        do i=0, 11
            if(i<5 .or. 7<i) then
                read(10,*) Plane_txt(i,0:16)
            else
                num=int((Plane_txt(4,0)+1)*(Plane_txt(4,1)+1)*(Plane_txt(4,2)+1))
                read(10,*) Plane_txt(i,0:num-1)
            end if
        end do
    close(10) !ファイルを閉じる

    iteration_max=int(time_max/time_step)+1 !最大反復回数の読み込み
    allocate(FlightLog(0:iteration_max,0:FlightLog_data))

    i_max_copy=iteration_max
    call NLFS(tmp,genom_length,genom_list,i_max_copy,pilot_method,FlightLog_data,FlightLog,NLFS_txt,Plane_txt)

    call Output_FlightLog(i_max_copy,iteration_max,FlightLog_data,FlightLog,genom_length,genom_list) !飛行解析データをcsvファイルで書き出し

    write(*,*) "------------------------------------------------------------"
    write(*,'(A15,A1,I3)') "generation","=",generation
    write(*,'(A15,A1,F8.3)') "MAX","=",Distance_Max
    write(*,'(A15,A1,F8.3)') "MIN","=",Distance_Min
    write(*,'(A15,A1,F8.3)') "AVERAGE","=",Distance_Ave
    write(*,'(A15)') "BEST GENOM"
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "phi","=",genom_list(0),"theta","=",genom_list(1),"psi","=",genom_list(2)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_pullup","=",genom_list(3),"dive_angle","=",genom_list(6),"cluse_angle","=",genom_list(7)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "time_bank","=",genom_list(4),"time_bank_end","=",genom_list(5),"bank_angle","=",genom_list(8)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[0]","=",genom_list(9),"PID_lon2[0]","=",genom_list(12)&
    ,"PID_lat1[0]","=",genom_list(15),"PID_lat2[0]","=",genom_list(18)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(10),"PID_lon2[1]","=",genom_list(13)&
    ,"PID_lat1[1]","=",genom_list(16),"PID_lat2[1]","=",genom_list(19)
    write(*,'(A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3,3X,A15,A1,F8.3)') &
    "PID_lon1[1]","=",genom_list(11),"PID_lon2[2]","=",genom_list(14)&
    ,"PID_lat1[2]","=",genom_list(17),"PID_lat2[2]","=",genom_list(20)

end subroutine

subroutine random(rnd)
    !疑似的な乱数(0~1)を生成する
    !https://qiita.com/ocian/items/e5eabe60c6b31cb48fac
    implicit none
    integer::seedsize
    integer,allocatable::seed(:)
    real(8),intent(OUT)::rnd
    integer::c !時間を入れる

    call system_clock(count=c) !時間を取得
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    call random_seed(get=seed)
    seed=c !時間を全部に代入
    call random_number(rnd) !rndに乱数をセット
end subroutine

subroutine Rotate_bector(axis,theta,matrix1,matrix2)
    !3行1列のベクトルを指定された軸まわりに指定された角度だけ回転させる．
    implicit none
    integer,intent(IN)::axis !0:x軸，1：y軸，2：z軸
    real(8),intent(IN)::theta
    real(8),dimension(0:2,0:2)::rtt
    real(8),dimension(0:2,0:0),intent(IN)::matrix1
    real(8),dimension(0:2,0:0),intent(OUT)::matrix2
    real(8),dimension(0:2,0:0)::matrix
    real(8)::l1,m1
    !回転行列の定義
    l1 = cos(theta)
    m1 = sin(theta)

    if(axis==0) then !x軸まわりの回転行列
        rtt = 0.0D0
        rtt(0,0) = 1.0D0
        rtt(1,1) = l1
        rtt(2,1) = -m1
        rtt(1,2) = m1
        rtt(2,2) = l1
    elseif(axis==1) then !y軸まわりの回転行列
        rtt = 0.0D0
        rtt(0,0) = l1
        rtt(2,0) = m1
        rtt(1,1) = 1.0D0
        rtt(0,2) = -m1
        rtt(2,2) = l1
    elseif(axis==2) then !z軸まわりの回転行列
        rtt = 0.0D0
        rtt(0,0) = l1
        rtt(1,0) = -m1
        rtt(0,1) = m1
        rtt(1,1) = l1
        rtt(2,2) = 1.0D0
    end if

    matrix = matrix1
    call matrix_multiplication(3,3,rtt,3,1,matrix,matrix2) !軸まわりに回転

end subroutine

subroutine Transform_coordinate(phi,theta,psi,matrix1,matrix2)
    !3行1列のベクトルを絶対座標系から機体軸座標に変換する
    implicit none
    real(8),intent(IN)::phi,theta,psi
    real(8),dimension(0:2,0:0),intent(IN)::matrix1
    real(8),dimension(0:2,0:0),intent(OUT)::matrix2
    real(8),dimension(0:2,0:0)::matrix

    matrix=matrix1
    call rotate_bector(2,psi,matrix,matrix2) !ヨー方向に回転
    matrix=matrix2
    call rotate_bector(1,theta,matrix,matrix2) !ピッチ方向に回転
    matrix=matrix2
    call rotate_bector(0,phi,matrix,matrix2) !ロール方向に回転

end subroutine

subroutine Inverse_Transform_coordinate(phi,theta,psi,matrix1,matrix2)
    !3行1列のベクトルを機体座標系から絶対座標系に変換する
    implicit none
    real(8),intent(IN)::phi,theta,psi
    real(8),dimension(0:2,0:0),intent(IN)::matrix1
    real(8),dimension(0:2,0:0),intent(OUT)::matrix2
    real(8),dimension(0:2,0:0)::matrix

    matrix=matrix1
    call rotate_bector(0,-phi,matrix,matrix2) !ロール方向に逆回転
    matrix=matrix2
    call rotate_bector(1,-theta,matrix,matrix2) !ピッチ方向に逆回転
    matrix=matrix2
    call rotate_bector(2,-psi,matrix,matrix2) !ヨー方向に逆回転

end subroutine

subroutine matrix_multiplication(m,n,matrix1,p,q,matrix2,matrix3)
    !行列の積の計算
    implicit none
    integer,intent(IN)::m,n,p,q
    real(8),dimension(0:m-1,0:n-1),intent(IN)::matrix1
    real(8),dimension(0:p-1,0:q-1),intent(IN)::matrix2
    real(8),dimension(0:m-1,0:q-1),intent(OUT)::matrix3
    real(8)::a
    integer::i,j,k
    matrix3 = 0.0D0

    if(n /= p) then !行と列の数の不一致を知らせる
        write(*,*) "IMPOSSIBLE!"
    else !行列の積
        do i = 0, m-1
            do j = 0, q-1
                a = 0.0D0
                do k = 0, n-1
                a = a+matrix1(i,k)*matrix2(k,j)
                end do
                matrix3(i,j) = a
            end do
        end do
    end if

end subroutine matrix_multiplication

subroutine Ans(n,A,b,x)
    !ガウスの消去法によって，連立方程式Ax=bを解く．
    !Aはn×n行列，x，bはn×1行列
    !参考文献：「解析塾秘伝！有限要素法のつくり方」　p.82-
    !引数
    integer,intent(IN)::n
    real(8),dimension(0:n-1,0:n-1),intent(INOUT)::A
    real(8),dimension(0:n-1),intent(INOUT)::b
    real(8),dimension(0:n-1),intent(OUT)::x(0:n-1) !解を入れる配列
    !変数の宣言
    real(8)::pivot !マトリックスの体格成分
    real(8)::p !計算に使用するマトリックスの成分
    !カウンター
    integer::r,c,rr,cc

    !前進消去
    do r = 0,n - 1
        !対角成分をpivotに代入
        pivot = A(r, r)
        do c = r,N - 1    
            A(r, c) = A(r, c) / pivot
        end do
        b(r) = b(r) / pivot
        do rr = r + 1,n - 1
            p = A(rr, r)
            do cc = r,N - 1
                A(rr, cc) = A(rr, cc) - p * A(r, cc)
            end do
            b(rr) = b(rr) - p * b(r)
        end do
    end do

    !後退代入
    do r = N - 1,0,-1
        x(r) = b(r)
        do c = r + 1,N - 1
            x(r) = x(r) - A(r, c) * x(c)
        end do
    end do

End subroutine

subroutine write1D(n,x) !一次元配列の表示
    implicit none
    integer,intent(IN)::n
    real(8),dimension(0:n),intent(IN)::x
    integer::i

    do i = 0, n
    write(*,*) i,x(i)
    end do

end subroutine
