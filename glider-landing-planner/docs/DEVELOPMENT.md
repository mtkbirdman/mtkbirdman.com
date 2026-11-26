# 開発者ガイド — Glider Landing Planner v0.1.20

**対象：** ソースの変更・検証・配布を担当する開発者。**目的：** 変更箇所を特定し、ビルド／テスト／WordPress更新まで安全に実施できること。本文は同梱の**v0.1.20ソースコード**を基準にしています。利用者操作は[README](../README.md)、数式は[THEORY](THEORY.md)、変更経緯は[CHANGELOG](../CHANGELOG.md)を参照してください。

## 1. 全体像と動作環境

このアプリは、WordPressのショートコード`[glider_landing_planner]`が出力する`.glp-app-root`に、ブラウザでTypeScript由来のES Modulesを読み込んで動作します。計画の編集・計算・保存・印刷は閲覧者のブラウザ内で行い、独自のサーバーAPIやアカウントはありません。

- **実行時：** WordPress/PHP（プラグインヘッダーはPHP 7.4以降）、ブラウザのES Modules・IndexedDB・Leaflet 1.9.4、国土地理院の淡色地図タイル。Leafletは`unpkg.com`から、タイルは`cyberjapandata.gsi.go.jp`から読み込まれます。両配信先へのネットワーク接続が必要です。
- **ソース：** TypeScript 5.8.3（`source/package.json`）、`strict`、ES2022、ES2020モジュール、DOM API。Vanilla DOMでUIを生成し、地図のみLeafletを利用します。
- **配布：** `source/src/*.ts`を`tsc`で`assets/js/*.js`へ出力する**非バンドル**ES Modules構成。Vite設定・コマンドは存在しますが、現行配布経路の正本は`build:offline`です。Vitest/Playwrightを使うテスト構成は同梱されていません。
- **ローカル試用：** `preview-local.html`は簡易プレビュー入口。WordPress/Cocoon上での最終検証を代替しません。`file://`直開きではモジュール読込の制限を受ける場合があるため、ローカルHTTPサーバーを使用します。

## 2. 主要ディレクトリ

```text
glider-landing-planner/
├── README.md                  # 初めて利用する人の操作手順
├── CHANGELOG.md               # バージョン別変更履歴
├── docs/
│   ├── DEVELOPMENT.md         # この文書
│   ├── THEORY.md              # 数値モデル
│   └── archive/implementation-spec-v1.0.md  # 歴史的設計資料
├── glider-landing-planner.php # ショートコード・JS/CSS登録
├── preview-local.html         # 簡易ローカル起動
├── assets/
│   ├── css/app.css            # 本体UI
│   ├── css/print.css          # 印刷専用CSS
│   └── js/                    # 実際に読み込まれるコンパイル済みES Modules
└── source/
    ├── .gitignore            # node_modulesの除外
    ├── package.json          # build:offline / build:vite / test
    ├── tsconfig.json
    ├── version-assets.mjs     # 相対importへ版数クエリを付加
    ├── vite.config.ts         # 将来用の任意ビルド設定
    ├── src/
    │   ├── main.ts            # UI・状態・再計算・保存の制御
    │   ├── model/             # 型定義・初期値・プリセット
    │   ├── geometry/          # WGS84/ENU・Bezierの幾何
    │   ├── flight/            # 飛行運動・沈下率
    │   ├── planning/          # 計画統合・接続判定
    │   ├── ui/                # 地図・グラフ・表・印刷・表示変換
    │   └── io/                # JSON検証・IndexedDB
    └── tests/*.test.mjs       # Node.jsの組込みテストランナー
```

`assets/js`は**生成成果物ですがWordPress配布には必須**です。TypeScriptだけを編集して配布JSを更新し忘れると、公開画面は変わりません。`source/`・`docs/`・`preview-local.html`は実行自体には不要です。

## 3. 初めに読む型とモジュール

| 変更したいもの | 入口 | 主な参照先 |
|---|---|---|
| 保存対象となる計画項目 | `model/types.ts`の`ProjectV1`、`FlightLeg` | `model/defaults.ts`、`io/project.ts` |
| 初期位置・高度・速度・新規レグ・左右場周 | `model/defaults.ts`の`createDemo()`、`makeStraight()`、`makeTurn()`、`trafficPreset()` | `main.ts`の`new`、`preset-left/right` |
| 入力と編集・Undo/Redo | `main.ts`の`handleChange()`、`handleAction()` | `edit()`、`applyProject()` |
| 計算全体・高度余裕・チェックポイント | `planning/plan.ts`の`buildPlan()` | `flight/core.ts`、`planning/connector.ts` |
| 風・直線・旋回・実効旋回角 | `flight/core.ts`の`simulateLeg()` | `headingForTrack()`、`simulatePath()` |
| 沈下率計算 | `flight/sink.ts`の`createSinkModel()` | `PerformanceConfig`、`flight/core.ts` |
| 5次Bezier・接続成立性 | `planning/connector.ts`の`evaluateConnector()` | `chooseConnector()`、`geometry/bezier.ts` |
| 地図上の操作・描画 | `ui/map.ts`の`PlannerMap` | `main.ts`のコールバックと`render()` |
| 高度グラフ・結果表 | `ui/chart.ts`、`ui/checkpoints.ts` | `main.ts`の`renderResults()` |
| PDF用印刷プレビュー | `ui/print.ts`の`openPrintReport()` | `reportHtml()`、`assets/css/print.css` |
| JSON保存・読込 | `io/project.ts`の`serializeProject()`、`loadProject()` | `validateProject()`、`parseProject()` |
| ブラウザ自動保存 | `io/autosave.ts`の`saveSnapshot()`、`latestSnapshot()` | `main.ts`の`scheduleSave()`、起動処理 |
| WordPress統合・表示キャッシュ | `glider-landing-planner.php` | `source/version-assets.mjs` |

中心データの`ProjectV1`は、`schemaVersion: '1.3'`、計画名、速度・風・性能（`performance.bestLD`、`performance.bestLDAtMps`）、離脱／着陸条件、両側の`FlightLeg[]`、任意の`BezierConnector`、表示設定を保持します。`PlanResult`は計算の戻り値であり、`PathResult`・`ConnectionResult`・`Sample`・`Checkpoint`・`Notice`を含みます。高度は内部でm MSL、速度はm/s、水平位置は着陸地点を基準とする局所ENUのmです。

## 4. 起動、状態更新、計算の入口

### 4.1 起動と自動復元

`glider-landing-planner.php`がCSS・Leaflet・`assets/js/main.js`を登録し、`main.js`は各`.glp-app-root`に`new App(root)`を作成します。コンストラクターは`createDemo()`→`buildPlan()`で初期計画を生成し、ツールバー・サイドバー・結果・地図を描画します。その後、`latestSnapshot()`を非同期実行します。

保存済みデータがあり、読込中に利用者編集が発生していなければ、**`loadProject()`で検証後**に`buildPlan()`し、画面へ復元します。起動中の編集で`editRevision`が変われば遅れて到着した自動復元は適用しません。復元失敗は表示に通知し、保存データを暗黙に削除しません。

### 4.2 一つの編集が反映されるまで

1. DOMの`change`／`click`、または地図コールバックが`handleChange()`／`handleAction()`へ入る。
2. 変更を`edit(mutate)`へ渡し、既存の`ProjectV1`を`structuredClone`して操作する。
3. 着陸側が0レグから1レグ以上になった場合は接続曲線を初期生成し、逆に0レグになれば接続を解除する（**既存JSON読込には適用しない**）。
4. `validateProject()`がスキーマ・数値範囲・ID等を検証し、不正なら変更を適用しない。
5. 変更内容に応じて`buildPlan()`または既存`PlanResult`の再利用を選択する。`recalculate()`は`adjustmentMode: 'auto'`で幾何・速度・風などが変わった場合のみ`chooseConnector()`でハンドル長を再探索する。
6. `applyProject()`が実効旋回角を`rememberResolvedTurnAngles()`で保存し、状態・地図・グラフ・チェックポイント・診断を再描画。編集前状態は最大50件のUndo履歴に積む。
7. `scheduleSave()`が1秒後のIndexedDB書込を予約する。遅い保存結果が新しい編集を「保存済」と誤表示しないよう、`editRevision`で確認する。

計画名・表示設定の変更は計算結果を再利用します。チェックポイント名は再計算しますが、接続最適化の対象から除外します。Undo/RedoとJSON読込は**保存済みハンドルを勝手に最適化しない**ことが重要な互換性仕様です。

### 4.3 計算エンジンとの境界

`buildPlan(project)`はDOM・Leafletに依存しない計算入口です。着陸地点を原点とする`LocalFrame`を作成し、`createSinkModel()`を一度構築して両経路へ渡します。

- `simulateLeg()`は**1レグ**の速度指定・風補正・位置／時間／高度積分・診断・目標方位からの実効旋回角解決を担当。
- `simulatePath()`は**複数レグ**を指定順（着陸側は逆順）に繰り返し、サンプル・セグメント・診断を統合。
- `evaluateConnector()`は**指定済みの制御ハンドル長**で曲線を生成し、曲率・必要バンク・風・沈下・自己交差を評価。
- `chooseConnector()`はハンドル長の有限候補探索によって接続の候補を選ぶ。探索が成功しても実際の成立を保証しない。
- `buildPlan()`は離脱側を順方向、着陸側をファイナル開始点から逆方向に計算し、着陸側の逆算済み経路を正本として予想高度を構築。ファイナル速度・風補正は共通化し、全制約確認後に`status`と`metrics`を確定。予想高度系列と必要高度系列、ファイナル目標線、チェックポイント・診断をまとめる。

`SegmentResult.start/end`は実サンプル両端、`completed`は水平経路の完了を示します。`PlanResult.metrics`は成立時のみ数値、不成立・未設計時はnullです。`checkedSink()`は純粋な沈下式と適用範囲を分離します。

`main.ts`の`render()`が`PlannerMap.redraw()`、`chartSvg()`、`checkpointTable()`へ`PlanResult`を渡します。`printReport()`は`openPrintReport()`へ計画・計算・地図DOM・集計HTMLを渡します。数式の詳細は[THEORY](THEORY.md)を参照してください。

## 5. UIと出力の責務

- **`ui/map.ts`：** Leafletの淡色タイル、レグ線、距離リング、方位線、風向矢印、地点・「＋」・Bezierハンドルを描画。地図背景はパン可能。地点・距離・ファイナル長のドラッグ編集は**廃止済み**。新規計画の青い250 m方向ガイドは描画専用で、飛行計算に投入しない。
- **`ui/responsive.ts`：** `matchMedia`でコンパクトレイアウトを判定。PCでは左サイドバーの該当編集欄へスクロール、コンパクト画面では下部パネルを表示。分岐はCSS幅だけでなくJavaScriptで制御。
- **`ui/chart.ts`：** `altitudeTimeSeries()`は計算済みの正本時刻を使用し、無効高度をまたいで線を結ばない。`chartSvg()`は青・緑を別系列で描画。未設計時は案内メッセージを描く。
- **`ui/checkpoints.ts`：** MSL高度と「着陸標高との差」である滑空場基準高度を区別し、ダウンウィンド開始の高度余裕行を生成する。AGLではない。
- **`ui/print.ts`：** 表示中のLeaflet DOM（画像タイル・SVG重ね描き）を複製して、別タブに独立した印刷HTMLを出力。地図範囲は印刷ボタンを押した時点の表示範囲。印刷は`window.print()`／ブラウザ側のPDF保存であり、サーバーPDF生成ではない。

## 6. データの保存・検証

| 処理 | 実装 | 要点 |
|---|---|---|
| 現行値の検証 | `validateProject(input)` | 入力を変異させず、スキーマ1.3・数値範囲・配列・ID・必須値・列挙型を厳密に検査 |
| 計画読込 | `loadProject(input)` | `structuredClone()`後にスキーマ1.3を検証。旧形式の移行なし |
| ファイル読込 | `parseProject(text)` | JSON構文解析、2,000,000文字上限、`loadProject()` |
| ファイル書出し | `serializeProject()`→`downloadProject()` | 現行形式を整形JSONとしてダウンロード（`.glp.json`） |
| 自動保存 | `saveSnapshot()` | IndexedDB `glp-autosave-v1`の`snapshots`へ保存、直近10件を保持 |
| 自動復元 | `latestSnapshot()` | 最新スナップショットを返し、`App`側で検証・表示 |

`loadProject()`はスキーマ1.3のみを受け付けます。性能設定の旧構造、旧スキーマ1.0〜1.2の移行・補完処理はありません。自動復元も同じ検証を使用し、旧形式の場合は通知して新規計画を維持します。既存IndexedDBのストア構造は変更せず、旧データの強制削除は行いません。

## 7. ビルドとテスト

### 7.1 必要条件

Node.js 22.16以降（24系も対象）とnpmを用意します。`source/package.json`は開発依存に`typescript: 5.8.3`、`vite: 8.3.8`を指定しています。**配布ZIPに`node_modules`と`package-lock.json`は同梱されていません。** 開発環境では次のように依存関係を準備します（レジストリへの接続が必要です）。

```powershell
Set-Location .\glider-landing-planner\source
npm install --include=dev --no-package-lock
npm run build:offline
npm test
```

`npm run build:offline`は`tsc -p tsconfig.json && node version-assets.mjs`を実行し、`../assets/js`にJSを生成して相対importへ`?v=0.1.20`を付加します。`npm test`は`pretest`で同じビルドを行った後、`node --test tests/*.test.mjs`を実行します。`npm run build:vite`は任意・将来用であり、**現行配布方式との互換性を確認せず本番アセットへ置き換えない**でください。

配布前に `npm test` を実行し、すべてのテストの成功を確認してください。v0.1.20ではNode.js 22.16.0、TypeScript 5.8.3で175件のテストに成功しました（基準166件、追加9件）。`plan-consistency.test.mjs`、`schema.test.mjs`、`map-render.test.mjs` に正則性・未完了レグ・保存範囲・ID出力・編集の確定順序に関する回帰試験を統合しました。既存のID照合を検査する静的アサーション1件は可逆トークンの照合に合わせて更新し、数値精度基準は変更していません。Node.js 24実機、実DOM・IndexedDB（Chromiumでローカル検証ページへのアクセスが管理ポリシーにより拒否）、WordPress実サイト、FFFTP転送、PDF実印刷は別途確認が必要です。`source/.gitignore`で`node_modules/`を除外します。

### 7.2 テストの探し方

`source/tests/`はNode.jsの組込みテストです。`engine.test.mjs`で風・直線・旋回・順逆計算・Bezier・高度余裕・JSON、`heights.test.mjs`と`chart.test.mjs`で高度・グラフ、`map-render.test.mjs`で地図描画、`sink.test.mjs`で沈下率の解析式と有限値検証、`schema.test.mjs`で保存形式・UI・旧形式拒否、`v013-turn-target.test.mjs`で方位指定、`v016-initial-flow.test.mjs`で新規画面と初回操作を検証します。その他の`v0xx-*.test.mjs`は既存挙動の回帰テストです。名前に旧バージョンが含まれても、現行実装への回帰確認として保持されています。

テストの成功は**ブラウザE2Eの成功を意味しません**。特に公開画面のクリック／タップ、画面サイズ、CocoonのCSS競合、タイル取得、別タブ印刷を対象端末で確認してください。

## 8. WordPressへの導入・更新

1. 既存プラグインフォルダーをバックアップする。旧形式の計画は1.3では読み込めないため、必要なら旧版でJSONを保存し、旧版と併せて保管する。
2. `npm test`を実行してビルドと回帰テストを確認する。`glider-landing-planner.php`でプラグイン版数、`source/package.json`でパッケージ版数、`model/defaults.ts`の`ENGINE_VERSION`、`main.ts`の画面版数、`ui/print.ts`のアセットクエリを必要なリリース時に整合させる。
3. 配布する`glider-landing-planner/`フォルダー内に`glider-landing-planner.php`と`assets/css/`、`assets/js/`があることを確認する。配布用ZIPの内側を二重の`glider-landing-planner/`にしない。
4. FFFTPではZIPを**ローカルで展開**し、WordPressの`wp-content/plugins/glider-landing-planner/`にフォルダー内容をバイナリ転送する。WordPressサーバー上でTypeScriptをビルドする必要はない。v0.1.16からの更新では旧 `assets/js/flight/polar.js` を削除する（上書き転送では消えない）。`source/`や`docs/`を公開サーバーへ転送することは実行上必須ではない。
5. WordPressでプラグインを有効化し、固定ページに**`[glider_landing_planner]`**を1個配置。既存ページを使う場合はショートコードを維持する。
6. アプリ画面を再読込し、JS/CSSキャッシュを確認。ブラウザキャッシュ、WordPress／Cocoon／CDNのキャッシュが残れば必要に応じて消去する。
7. 新規計画、左右場周、離脱側レグ追加、直線・旋回編集、JSON保存と読込、旧JSONの拒否、印刷プレビューをPC／スマートフォンの対象環境でスモークテストする。

`glider-landing-planner.php`は入口の`assets/js/main.js`へWordPressの版数クエリを付けます。**相対ES Modulesのimportはブラウザが個別にキャッシュするため**、`source/version-assets.mjs`が内部importにも`?v=<package.jsonのversion>`を付与します。**ドキュメントだけの更新では不要な版数変更を行いません。** コード／アセットを変更して配布版を上げる際は、PHP・画面・JS相対import・印刷CSS・パッケージ等を整合させて再ビルドします。

## 9. Gitでのリリース運用

この配布ZIP自体には`.git`ディレクトリがないため、**実際のリポジトリのブランチ名や最新コミットは確認できません**。以下はリポジトリのルートを`glider-landing-planner`とした例です。作業前に`git status`で未コミット変更を確認してください。

```powershell
git status --short
# 変更ファイルを配置してから
npm --prefix source test
git add -u -- .
# 除外設定を確認し、必要な新規ファイルだけを登録
git check-ignore -v --no-index -- source/node_modules/typescript/package.json
git add -- source/tests/plan-consistency.test.mjs
git diff --cached --check
git diff --cached --stat
git commit -m "fix: stabilize connector validity, partial results, persistence and HTML output for v0.1.20"
```

スキーマ1.3は維持します。既存レグの`guidance`と新規レグの`ground-track`既定値は変更しません。正常な旧計画は読込可能ですが、計算適用範囲により不成立になる場合があります。プラグイン・パッケージ・配布アセットの版数を整合させます。

## 10. リリース前チェック

- [ ] 新規計画が**離脱側0件・着陸側0件・接続なし**で始まり、方向ガイドが計算に入らない。
- [ ] 着陸側最初のレグ追加・左右場周で自動接続し、着陸側0件に戻せば接続を外す。
- [ ] レグ編集、機首方位指定／進行方向指定、変化量指定／方向指定、多周旋回が既存テストで維持される。
- [ ] 簡易抗力式の沈下率、必要高度と予想高度、診断・評価不可が維持される。
- [ ] スキーマ1.3のJSON／IndexedDB復元、旧形式の拒否、Undo/Redoを確認する。
- [ ] TypeScriptソース変更時は`assets/js`を再生成し、相対importの版数が揃っている。
- [ ] `npm test`が成功し、`git diff --check`に問題がない。
- [ ] WordPress実機の地図／パネル／PDF・キャッシュを確認する（**自動テストで代替しない**）。

### v0.1.20のデータ契約と出力境界

- `model/limits.ts` はハンドル長の包含範囲1〜50,000 mを共有する。自動候補は範囲外を除外し、有限な手動ドラッグは範囲に制限する。非有限の操作値と不正JSONは拒否する。計算不成立は保存拒否の理由にしない。
- 編集・復元では候補を用意し、入力検証、計算、解決済み旋回角等の反映、最終検証を終えてから履歴と計画を確定する。自動保存でも保存前に検証する。Undo/Redo・JSON読込・自動保存復元で再最適化しない。
- `completed` のレグだけ正式チェックポイントを生成する。未完了の部分経路は破線とし、非正則Bezierの点列は幾何プレビューとしてのみ扱う。
- レグIDはHTML属性出力時に `htmlLegId()` でエスケープし、CR・NUL・サロゲート等を含むIDと予約接頭辞の衝突を可逆トークンで処理する。DOMからの操作と選択照合は `readLegId()` で元のIDを復元する。JSONのIDは変更しない。計画名、ラベル、診断、ツールチップ、PDFを含めて出力文脈を限定レビューした。WordPress全体のセキュリティ監査ではない。
