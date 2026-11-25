# 計算理論と数値モデル — Glider Landing Planner v0.1.16

**対象：** 航空工学・飛行力学・数値計算の検証者。**目的：** 現行コードが実際に評価する量・数式・仮定・失敗条件を、対応する関数とテストに結び付けて説明すること。一般のグライダー飛行性能や実際の安全性を保証する文書ではありません。

**適用対象は同梱v0.1.16のTypeScript実装です。** 旧設計書に登場する固定L/Dの「モデルA」は現行コードで**削除済み**です。以下のB/Cと5次Bezierの数式はソースに照合しています。操作は[README](../README.md)、開発手順は[DEVELOPMENT](DEVELOPMENT.md)、旧設計の経緯は[CHANGELOG](../CHANGELOG.md)を参照してください。

## 1. 計算対象と成立性の意味

計算するのは、一定の入力風と指定対気速度を仮定した**水平航跡**（局所ENU座標）、飛行時間、移動距離、沈下率に基づく高度、5次Bezierの接続曲線に必要なバンク角、チェックポイント、診断です。高度は離脱側から前進積分した**予想高度**と、ファイナル開始の目標高度から逆方向に積分した**必要高度**を別の値として保持します。

ここでの`feasible`は接続曲線の**数値的・幾何学的な成立性**、`valid`は生成された診断に`error`がないことを示します。実機の飛行可能性、失速余裕、法規・空域、地形・障害物クリアランスを判定したものではありません。最大バンク制限の計算上の成立と、飛行の安全性は同義ではありません。対応：`planning/plan.ts::buildPlan()`、`planning/connector.ts::evaluateConnector()`。

## 2. 座標系、高度、単位

### 2.1 緯度・経度と局所ENU

入出力の水平地点はWGS84楕円体の緯度・経度（十進度）です。`geometry/geo.ts::LocalFrame`は**着陸地点**を原点に取り、地心地球固定座標（ECEF）から局所東E・北Nへの回転変換を行います。WGS84のパラメータは

\[
a=6\,378\,137\ \mathrm{m},\qquad f=\frac{1}{298.257223563},\qquad e^2=f(2-f)
\]

で、緯度\(\varphi\)・経度\(\lambda\)、計算上の楕円体高\(h_e=0\)から

\[
N(\varphi)=\frac{a}{\sqrt{1-e^2\sin^2\varphi}},\quad
\begin{bmatrix}X\\Y\\Z\end{bmatrix}
=\begin{bmatrix}
(N+h_e)\cos\varphi\cos\lambda\\
(N+h_e)\cos\varphi\sin\lambda\\
(N(1-e^2)+h_e)\sin\varphi
\end{bmatrix}
\]

を求めます。着陸地点との差\((\Delta X,\Delta Y,\Delta Z)\)を使い、原点緯度\(\varphi_0\)・経度\(\lambda_0\)で

\[
\begin{aligned}
E&=-\sin\lambda_0\,\Delta X+\cos\lambda_0\,\Delta Y,\\
N&=-\sin\varphi_0\cos\lambda_0\,\Delta X
   -\sin\varphi_0\sin\lambda_0\,\Delta Y
   +\cos\varphi_0\,\Delta Z.
\end{aligned}
\]

高度座標（Up）はこの局所変換では直接利用せず、**MSL高度を独立したスカラー**として計算します。逆変換`toGeo()`はENUの水平平面からECEFへ戻し、楕円体緯度を**8回反復**して復元します。このため長距離・大高度差での測地精度を保証するものではありません。検証：`source/tests/engine.test.mjs`の`T02`。

### 2.2 高度の基準

- **MSL高度** \(H_{\mathrm{MSL}}\)：離脱高度・着陸標高・内部サンプルの高度。内部単位はm。
- **滑空場基準高度**：\(H_{\mathrm{rel}}=H_{\mathrm{MSL}}-H_{\mathrm{landing,MSL}}\)。着陸地点の標高をゼロとする相対値。
- **AGL高度**：各位置の地表面からの高さ。**地形標高を参照していないため計算しません。** 滑空場基準高度をAGLと読んではいけません。

表示時のみ\(1\ \mathrm{ft}=0.3048\ \mathrm{m}\)、\(1\ \mathrm{km/h}=1/3.6\ \mathrm{m/s}\)で変換します。高度の数値入力も画面単位に合わせて内部mへ戻します。対応：`ui/height-units.ts`、`ui/checkpoints.ts::landingRelativeHeight()`、`ui/chart.ts::fmtAlt()`。検証：`engine.test.mjs::T01`、`heights.test.mjs`。

新規計画の仮値は、妻沼滑空場の**代表座標**（36.211388889°N、139.418888889°E）、着陸標高89 ft MSL、離脱地点は代表点から真方位320°・850 m、開始機首方位と進入真方位140°、離脱高度は着陸標高＋1200 ft、基準対気速度90 km/hです。代表座標は正式な滑走路端・接地点を意味せず、性能設定も実機を同定するものではありません。対応：`model/defaults.ts::createDemo()`。初期代表位置・標高の出典として、旧READMEは[日本学生航空連盟](https://www.jsal.or.jp/)の妻沼滑空場資料を挙げています。初期最大滑空比34の参考として[Alexander Schleicher](https://www.alexander-schleicher.de/)のASK 21資料を挙げていますが、**初期速度90 km/hでの実機性能を検証する実測データは同梱していません。**

## 3. 速度、風、機首方位、対地航跡

全方位は北を0°・東を90°とする**真方位**で、風向\(\theta_w\)は風が**吹いてくる方向**です。入力する\(V\)は対気速度（m/s）、\(W\)は風速（m/s）。局所ENUでの北基準の単位ベクトルを

\[
\boldsymbol{u}(\theta)=\begin{bmatrix}\sin\theta\\\cos\theta\end{bmatrix}
\]

（式内角はラジアン）とすると、風の対地速度ベクトルと機体の対地速度は

\[
\boldsymbol{v}_w=-W\boldsymbol{u}(\theta_w),\qquad
\boldsymbol{v}_g=V\boldsymbol{u}(\psi)+\boldsymbol{v}_w
\]

です。\(\psi\)は機首真方位、\(\|\boldsymbol{v}_g\|\)が対地速度、\(\operatorname{atan2}(v_{g,E},v_{g,N})\)が対地航跡角です。対応：`flight/core.ts::windVector()`、`airGround()`。

**機首方位指定（`heading`）：** 指定真方位を機首方位として使用するため、横風があれば地上航跡は偏流します。

**対地航跡指定（`ground-track`）：** 指定地上航跡角\(\chi\)を維持するため、\(\boldsymbol{t}=\boldsymbol{u}(\chi)\)、\(w_\parallel=\boldsymbol{v}_w\cdot\boldsymbol{t}\)、\(w_\perp^2=\|\boldsymbol{v}_w\|^2-w_\parallel^2\)として

\[
G=w_\parallel+\sqrt{V^2-\max(0,w_\perp^2)},\qquad
\boldsymbol{v}_a=G\boldsymbol{t}-\boldsymbol{v}_w
\]

から機首方位を逆算します。平方根の中身が実装許容値以下、または\(G\)がほぼゼロ以下の場合は`headingForTrack()`が`null`を返し、`WIND_TRACK_INFEASIBLE`等の診断となります。これは風補正の運動学的解であり、実際の操縦可能性を判定しません。検証：`engine.test.mjs::T06`〜`T08`。

地図補助方位線にのみ使う磁気偏角\(D\)（東偏を正）について\(\theta_{\mathrm{true}}=\theta_{\mathrm{mag}}+D\)とします。**風・機首・進入経路の計算は磁気偏角で補正しません。** 対応：`geometry/geo.ts::magneticToTrue()`、`ui/map.ts::redraw()`。

## 4. 直線・旋回・順逆積分

### 4.1 状態と1レグ積分

各`State`は東北座標\((E,N)\)、MSL高度\(H\)、機首方位\(\psi\)、航跡角、符号付きバンク角、対気／対地速度、積算時間\(t\)、経路距離\(s\)、高度計算の有効性を持ちます。`FlightLeg`は直線（`lengthM`）または旋回（方向・角度・バンク）です。各レグの速度は`speedOverrideMps ?? defaultSpeedMps`を使用します。

直線では対地速度\(G\)を求め、指定距離\(L\)から\(\Delta t=L/G\)を計算します。方向付き積分では\(\varepsilon=+1\)（離脱側順方向）／\(-1\)（着陸側逆方向）を用い、

\[
\Delta\boldsymbol{x}=\varepsilon\,\boldsymbol{v}_g\Delta t,\quad
\Delta H=-\varepsilon\,w(V,0)\Delta t,\quad
\Delta s=\varepsilon L
\]

として数値サンプルを作ります。直線長からの分割数は\(\max(1,\min(2000,\lceil L/20\rceil))\)です。速度変更時は**瞬時の切替**とし、加減速・運動エネルギー変化を計算せず`SPEED_STEP_IDEALIZED`を通知します。接続する直線の目標機首方位が不連続なら`BOUNDARY_DISCONTINUITY`を通知します。対応：`flight/core.ts::simulateLeg()`。検証：`engine.test.mjs::T03`、`T13`。

### 4.2 バンク一定旋回

重力加速度\(g=9.80665\ \mathrm{m/s^2}\)、符号付きバンク角\(\phi\)（右旋回＋、左旋回−）を用い、機首真方位の旋回率を

\[
\dot\psi=\omega=\frac{g\tan\phi}{V}\quad[\mathrm{rad/s}]
\]

とします。無風なら旋回半径は\(R=V^2/(g|\tan\phi|)\)、角度\(\alpha\)の所要時間は\(T=\alpha/|\omega|\)です。これらはバンク一定・速度一定・協調旋回を仮定した関係であり、ロールイン／ロールアウトは含みません。速度ベクトルから位置を4次Runge–Kutta相当の重み付き時間積分で進め、方位は上記\(\omega\)から更新します。高度は区間中の沈下率\(w(V,|\phi|)\)を用いて積分します。

旋回サンプリングの角度幅は\(\max(0.25^\circ,\alpha/50000)\)で、\(2^\circ\)を超える必要がある場合は`TURN_COMPLEXITY_LIMIT`となります。経路全体でも12万サンプル以上は`PATH_COMPLEXITY_LIMIT`となります。旋回角を単に720°で切り捨てる処理はありません。対応：`flight/core.ts::advanceTurn()`、`simulateLeg()`、`simulatePath()`。検証：`engine.test.mjs::T04`、`T05`、`T24`。

**機首方位指定の旋回：** 指定角\(\alpha\)の**機首方位変化**を積分します。**対地航跡指定の旋回：** 一定の地上円弧半径を維持する制御ではなく、サンプルごとの**対地航跡角の累積変化**が指定角へ達するまで積分します。終端は最大25回の二分探索で調整します。風で航跡角変化が非単調、対地速度がほぼゼロ、または積分点数を超えた場合は`TURN_TARGET_NOT_REACHED`です。

### 4.3 旋回角指定／目標方位指定

「旋回角指定」は保存値`angleDeg`を使用します。「目標方位指定」は指定した方向（機首または航跡の真方位）に到達する**正の旋回角**を逆算します。離脱側は旋回後、着陸側は旋回前の方位を境界条件とします。右旋回を\(+1\)、左旋回を\(-1\)とする符号\(\sigma\)、順逆の向きを考慮すると、角度の余りは

\[
\alpha_0=\operatorname{mod}_{360}\!\left(\sigma\bigl[\psi_{\mathrm{target}}-\psi_{\mathrm{boundary}}\bigr]\right)
\]

（着陸側逆方向では差の符号を反転）です。余りが0°なら最小の正の解を**360°**とし、\(\alpha_0+360^\circ k\ (k\ge0)\)のうち直前の`angleDeg`に最も近い正の値を選びます。同距離なら小さい方を採ります。着陸側で解いた角度は、そのまま順方向の高度計算へ渡します。対応：`flight/core.ts::closestTurnAngle()`、`simulateLeg()`、`planning/plan.ts::rememberResolvedTurnAngles()`。検証：`source/tests/v013-turn-target.test.mjs`。

### 4.4 離脱側と着陸側

`planning/plan.ts::buildPlan()`は離脱地点から`departureLegs`を順方向（\(\varepsilon=+1\)）、ファイナル開始点から`arrivalLegs`を**配列逆順で時間逆方向**（\(\varepsilon=-1\)）に積分します。着陸側の配列は**実際の飛行順**で保存しており、画面で上流レグを追加すると配列の先頭に挿入します。各逆向きステップで時間が負になるため高度は増加し、ファイナル開始に必要な高度から上流の必要高度が求まります。別途、接続後の予想高度を`simulatePath(..., reverse=false)`で前進再計算します。検証：`engine.test.mjs::T09`、`v013-turn-target.test.mjs`。

## 5. 沈下モデルB：簡易抗力モデル

`flight/polar.ts::createSinkModel()`で`model='quadratic-drag'`を選んだ場合、対気速度\(V>0\ \mathrm{m/s}\)、バンク角\(\phi\)、入力した最大滑空比\(E>0\)とその速度\(V_g>0\ \mathrm{m/s}\)から

\[
\begin{aligned}
a&=\frac{1}{2EV_g^2},& b&=\frac{V_g^2}{2E},\\
w_B(V,\phi)&=aV^3+\frac{b}{V\cos^2\phi}.
\end{aligned}
\]

を評価します。\(w_B\)は**下向きの正の沈下率**（m/s）、\(a\)の単位はs²/m²、\(b\)はm²/s²です。無バンク・最良滑空速度では

\[
w_B(V_g,0)=\frac{V_g}{E},\qquad \frac{V_g}{w_B(V_g,0)}=E
\]

となります。**最大滑空比と最小沈下率を混同しない**でください。この式は簡易パラメータモデルであり、機体の重量・実測ポーラー・密度・失速制限に対する精度保証はありません。入力バンク角は保存形式では1〜70°（旋回レグ）、最大バンク制限は5〜70°で検証します。接続区間のバンクは形状から算出するため、制限超過は診断で扱います。対応：`flight/polar.ts::createSinkModel()`、`io/project.ts::validateProject()`。検証：`engine.test.mjs::T03`、`T12`、`source/tests/v012-polar.test.mjs`。

## 6. 沈下モデルC：速度別ポーラーモデル

モデルC（`model='polar-table'`）では、\((V_i,w_i)\)を**速度の厳密な昇順、4〜64点**で保存します。内部単位はm/s、画面入力は速度km/h・沈下率m/sです。各沈下率は正の有限値、速度は有効範囲内でなければなりません。実測データは同梱されていません。

### 6.1 PCHIP補間

`flight/polar.ts::pchip()`は\(h_i=V_{i+1}-V_i\)、\(\delta_i=(w_{i+1}-w_i)/h_i\)を計算します。内点では\(\delta_{i-1}\delta_i\le0\)なら傾き\(m_i=0\)、同符号なら

\[
m_i=\frac{q_1+q_2}{q_1/\delta_{i-1}+q_2/\delta_i},\quad
q_1=2h_i+h_{i-1},\quad q_2=h_i+2h_{i-1}
\]

という重み付き調和平均です。端点は片側3点から推定した導関数を、隣接差分の符号と大きさに応じて制限します。区間\([V_i,V_{i+1}]\)内の正規化座標\(t=(V-V_i)/h_i\)に対して、

\[
\begin{aligned}
P(V)={}&(2t^3-3t^2+1)w_i+(t^3-2t^2+t)h_im_i\\
&+(-2t^3+3t^2)w_{i+1}+(t^3-t^2)h_im_{i+1}
\end{aligned}
\]

を使います。**外挿はしません。** 範囲外（実装の許容差\(10^{-10}\) m/sを超える場合）は`POLAR_SPEED_OUT_OF_RANGE`、補間結果が正の有限値でなければ`POLAR_RESULT_INVALID`です。適用前の画面編集値は計算に使用されません。対応：`flight/polar.ts::pchip()`、`main.ts`の`polar-apply`、`io/project.ts`。検証：`source/tests/v012-polar.test.mjs`。

### 6.2 旋回補正と参照速度

旋回時の荷重倍数を\(n=1/\cos\phi\)と仮定し、コードは

\[
w_C(V,\phi)=n^{3/2}P\!\left(\frac{V}{\sqrt n}\right)
\]

で沈下率を算出します。**実際のレグ速度\(V\)と、補間テーブルを参照する速度\(V/\sqrt n\)は異なります。** レグ速度が範囲内でも参照速度が下限未満なら計算不能です。この式は荷重倍数の簡易相似補正であって、実測旋回ポーラーの再現保証ではありません。

入力ポーラーの重量・機体形態・IAS/TAS定義や大気密度をコードは検証しません。モデルC利用時は測定条件が揃った信頼できるデータを入力者が用意する必要があります。ポーラー参照に失敗したレグは`altitudeValid=false`を伝播し、位置・時間の計算を可能な限り継続します。接続区間でも参照に失敗すれば`feasible=false`かつ高度評価不能となります。対応：`flight/core.ts::safeSink()`、`planning/connector.ts::evaluateConnector()`、`planning/plan.ts::buildPlan()`。

## 7. 5次Bezier接続

### 7.1 境界条件と6制御点

離脱側終端\(A\)と、逆算した着陸側始端\(B\)を接続します。端点の**対地航跡**から単位接線\(\boldsymbol{t}_A,\boldsymbol{t}_B\)、直前・直後のバンクと風から**符号付き地上曲率**\(\kappa_A,\kappa_B\)を求めます。ハンドル長\(a,b>0\)を用いた6制御点は

\[
\begin{aligned}
P_0&=A, &P_1&=A+a\boldsymbol t_A,\\
P_2&=A+2a\boldsymbol t_A+\tfrac54\kappa_Aa^2\boldsymbol n_A,\\
P_3&=B-2b\boldsymbol t_B+\tfrac54\kappa_Bb^2\boldsymbol n_B,\\
P_4&=B-b\boldsymbol t_B, &P_5&=B,
\end{aligned}
\]

です。\(\boldsymbol n\)は進行方向ベクトルの平面左法線\((-t_N,t_E)\)。\(u\in[0,1]\)で

\[
\boldsymbol B(u)=\sum_{i=0}^{5}\binom{5}{i}(1-u)^{5-i}u^iP_i,\qquad
\kappa(u)=\frac{B'_E B''_N-B'_N B''_E}{\|\boldsymbol B'(u)\|^3}
\]

を評価します。この構築で端点位置・接線方向・符号付き曲率が一致します。**端点曲率をゼロへ強制していません。** 実際のロールイン／ロールアウト過程を自動生成するものでもありません。対応：`geometry/bezier.ts::{controlPoints,point,derivative,second,curvature}`。検証：`engine.test.mjs::T10`。

### 7.2 端点曲率と風の補正

`endpointCurvature()`は、\(\omega=g\tan\phi/V\)で変化する対気速度ベクトルと、一定風を足した対地速度\(\boldsymbol v_g\)を使い、

\[
\kappa_{\mathrm{end}}=\frac{v_{g,E}\dot v_{g,N}-v_{g,N}\dot v_{g,E}}{\|\boldsymbol v_g\|^3}
\]

を計算します（バンクがほぼ0°ならゼロ扱い）。曲線上では単位接線\(\boldsymbol t\)に対する風の接線・法線成分を\(W_t,W_n\)とし、

\[
q=\sqrt{V^2-W_n^2},\qquad G=W_t+q,\qquad
\omega=-\kappa\frac{G^2}{q},\qquad
\phi=\arctan\left(\frac{\omega V}{g}\right)
\]

で必要バンクを計算します。\(G>0\)、\(q>0\)、曲線の接線と曲率が有限であることを必要とします。バンクの絶対最大値が設定上限を超えれば`BANK_LIMIT_EXCEEDED`です。風・曲率条件不成立、自己交差（サンプル線分間の交差）、モデルCの範囲外も不成立判定に含みます。対応：`planning/connector.ts::endpointCurvature()`、`evaluateConnector()`。

### 7.3 曲線長・沈下積分・探索

`evaluateConnector()`はBezierを有限点でサンプリングします。通常評価の分割数は概ね**3 m当たり1分割**を目安に**160〜2000**の範囲（実装では直線距離と両ハンドル長の和から決定）とし、探索用軽量評価は**56分割**です。離散点の距離\(\Delta s\)を加算して曲線長を近似し、隣接2点の対地速度の逆数と沈下率の**台形則**で時間と高度損失を積分します。

`chooseConnector()`は端点距離\(d\)に対しハンドル長\(\max(10,rd)\)、\(r\in\{0.16,0.27,0.40,0.58,0.8,1.1\}\)の36候補を軽量評価し、成立性・最大バンク・曲線長等でスコアを比較して最良候補を通常精度で再評価します。**有限候補の探索なので、成立可能な解が存在しても発見する保証はありません。** 自動モードでは幾何・風・性能等の変更時に再探索し、ハンドルドラッグによる手動モードでは保存済み長さを直接評価します。対応：`planning/connector.ts::chooseConnector()`、`main.ts::recalculate()`。検証：`engine.test.mjs::T10`、`T11`、接続関連の`engine.test.mjs`と各版の回帰テスト。

## 8. ファイナル、高度プロファイル、高度余裕

### 8.1 ファイナルの幾何学的目標

ファイナルの地上長さを\(L_f\)、目標パス角を\(\gamma\)、着陸標高を\(H_L\)とすると、ファイナル開始位置は着陸目標点から進入方向の逆向きに\(L_f\)だけ戻し、

\[
H_{f,\mathrm{start}}=H_L+L_f\tan\gamma
\]

を必要高度の始点とします。ファイナルの目標線は距離\(x\)に対し\(H_L+(L_f-x)\tan\gamma\)です。**これはブレーキやフレアで実際に達成できる降下角の評価ではありません。** 対応：`planning/plan.ts::buildPlan()`、`FINAL_BRAKE_UNKNOWN`診断。

ファイナル開始の逆向き計算では、直前の着陸側レグの個別速度があればその速度を使用します。一方、表示用のファイナル予想線と沈下損失は**基準対気速度**で計算しています。これはv0.1.16の実装上の扱いであり、速度を全区間で厳密に連続させるモデルではありません。速度差の効果を実機精度で解釈しないでください。

### 8.2 二つの高度系列

- **離脱側予想高度（青）：** 離脱MSL高度から順方向に沈下積分し、接続曲線を含めて表示。接続後の着陸側・ファイナルについても別途順方向の予想高度を内部計算します。
- **着陸側逆算必要高度（緑）：** ファイナル開始の幾何学的目標高度から着陸側レグを逆積分。ファイナル目標線も緑の系列に含めます。

グラフの時間横軸は離脱からの**推定経過飛行時間**です。着陸側は逆向き積分で負の時刻を持つため、`ui/chart.ts::altitudeTimeSeries()`が接続後の時刻へ**表示上のみ**平行移動して配置します。未接続または接続制約不成立時の着陸側時刻は**仮配置**です。未設計（着陸側レグ0件・接続なし）ではグラフの代わりに案内を表示します。検証：`source/tests/chart.test.mjs`。

### 8.3 高度余裕の正確な定義

画面の主要指標「高度余裕」は接続境界の高さの差として

\[
\Delta H_{\mathrm{connector}}
=H_{\mathrm{departure\ end}}-\Delta H_{\mathrm{Bezier\ loss}}-H_{\mathrm{arrival\ start,required}}
\]

で計算します。\(\Delta H>0\)は**入力条件と計算モデルの下での**余剰、負なら不足です。表示条件は接続あり・`connector.feasible`・計算値有効・`PlanResult.valid`であることです。**正の値は安全余裕ではありません。** 対応：`main.ts::renderResults()`、`planning/plan.ts::buildPlan()`。

別途、ダウンウィンド開始チェックポイントにおける高度差は、接続後を順方向に計算した予想高度\(H_{d,\mathrm{pred}}\)と、同位置の逆算必要高度\(H_{d,\mathrm{req}}\)との差\(\Delta H_d=H_{d,\mathrm{pred}}-H_{d,\mathrm{req}}\)です。これは**ダウンウィンド印がある場合のみ**計算されます。チェックポイント一覧の高度余裕行と、画面上部の接続境界高度余裕は評価位置が異なります。対応：`planning/plan.ts::buildPlan()`の`marginAvailable`、`ui/checkpoints.ts::checkpointRows()`。検証：`engine.test.mjs::T14`。

「**—**」は初期未設計で評価指標を表示しない状態を表します。「**評価不可**」は接続が不成立・未接続、沈下率の参照失敗、または別の計算エラー等で必要条件を満たさない状態です。個別のエラーコードは`Notice`として保存され「計算上の診断」に表示されます。接続が成立しても、実際の操縦・飛行安全性を判定したことにはなりません。

## 9. 数値計算と検証の観点

| 対象 | 数値的な実装 | 検証・観察対象 |
|---|---|---|
| WGS84/ENU | ECEF変換、逆変換は緯度8回反復 | `engine.test.mjs::T02`の往復誤差 |
| 直線 | 長さ分割、一定風・一定速度・一定沈下率 | `T03`（距離／高度損失）、`T06`〜`T08`（風） |
| 旋回 | 方位積分とRK4相当の位置積分、角度分割 | `T04`（解析半径／所要時間）、`T05`（1周閉合）、`T24` |
| 順逆計算 | 着陸側の配列を逆順・負時間で計算 | `T09`、`v013-turn-target.test.mjs` |
| 沈下B | 2係数の簡易モデル | `T12`（最良滑空点）、`v012-polar.test.mjs` |
| 沈下C | 形状保存PCHIP、範囲外拒否、荷重倍数補正 | `v012-polar.test.mjs`の端点・補間・範囲外 |
| 5次Bezier | 解析的1・2階導関数、制御点、サンプル弧長 | `T10`（端点接線／曲率）、`T11`（弧長収束） |
| 結果・高度余裕 | 接続境界差、ダウンウィンド地点差、独立2系列 | `T14`、`chart.test.mjs`、`heights.test.mjs` |
| 保存互換性 | スキーマ検証・旧版移行 | `engine.test.mjs::T15`、`T16`、`v014-behavior.test.mjs`等 |

数値的な**検証**と、実飛行性能に対する**妥当性確認（validation）**は別です。同梱テストはコードの数値的一貫性・回帰を対象とし、実機データによる妥当性確認を提供しません。風、速度、旋回バンク、接続距離、ポーラー速度境界を変えた独立計算は追加検証として有用ですが、実施していない比較を「検証済」と呼ばないでください。

## 10. 診断と実装上の限界

### 10.1 代表的な診断（コード上の条件）

| コード | 主な意味 |
|---|---|
| `PLAN_NOT_STARTED` | 着陸側0件・接続なしの未設計（情報通知） |
| `CONNECTOR_MISSING` | 着陸側レグがあるのに接続曲線がない |
| `CONNECTOR_NO_SOLUTION` | 接続曲線の速度・風・曲率・探索条件が不成立 |
| `CONNECTOR_SELF_INTERSECTION` | 曲線サンプルが自己交差 |
| `BANK_LIMIT_EXCEEDED` | レグまたはBezierの必要バンクが設定上限超過 |
| `WIND_TRACK_INFEASIBLE` | 対地航跡やファイナルの指定航跡を風に対して維持できない |
| `TURN_TARGET_NOT_REACHED` | 対地航跡指定旋回で目標角度へ到達しない |
| `TURN_COMPLEXITY_LIMIT`／`PATH_COMPLEXITY_LIMIT` | 旋回または経路サンプル数の計算量上限 |
| `POLAR_SPEED_OUT_OF_RANGE`／`POLAR_RESULT_INVALID` | 参照速度範囲外または補間値不正で高度計算不能 |
| `SPEED_STEP_IDEALIZED`／`BOUNDARY_DISCONTINUITY` | 瞬時速度変更・機首方位不連続を理想化 |
| `ALTITUDE_DEFICIT`／`MARGIN_NOT_EVALUABLE` | ダウンウィンド地点の計算上の高度不足・評価不能 |
| `FINAL_BRAKE_UNKNOWN` | ファイナル角度は幾何学的目標でスポイラー／フレア未評価 |

コードは全診断を一覧表として表示するわけではなく、`Notice`のメッセージを最大12件表示しています。区別が必要な場合は`flight/core.ts`、`planning/connector.ts`、`planning/plan.ts`を確認してください。

### 10.2 モデル化していない主な現象

- 実機の失速速度・最大荷重倍数・操縦可能範囲、運用限界の判定。
- IAS/TASの厳密な換算、重量・重心・大気密度・高度による性能変化、実機ポーラーの自動補正。
- 旋回開始／終了時のロール速度・ロールイン／ロールアウト、および速度変更の加減速・エネルギー交換。
- 地形・障害物・空域・クリアランス、地表面の各点に対するAGL。
- 風の空間・時間・高度変動、鉛直気流、突風、乱流、サーマル。
- スポイラー・エアブレーキ・フレア・接地運動、滑走路運用の成立性。
- 有限候補のBezier探索で見つからない解、曲線の数値分割による幾何学・高度積分の近似誤差。

**結論：** この実装が表示するのは、選択した理想化モデルに対する経路と高度の**計算結果**です。現実の安全余裕や着陸成立を示していません。
