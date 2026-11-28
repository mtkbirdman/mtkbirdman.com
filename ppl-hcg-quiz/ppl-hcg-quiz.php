<?php
/**
 * Plugin Name: PPL HCG Quiz
 * Description: 自家用操縦士（滑空機）の学科試験問題をショートコードで出題します。
 * Version: 1.1.0
 * Requires at least: 6.0
 * Requires PHP: 7.4
 * Text Domain: ppl-hcg-quiz
 */

if (!defined('ABSPATH')) {
    exit;
}

add_shortcode('ppl_hcg_quiz', 'ppl_hcg_quiz_render');

function ppl_hcg_quiz_render(): string
{
    // 問題データを読み込む。CSVの列順は既存の question.csv に合わせる。
    $file = @fopen(__DIR__ . '/question.csv', 'r');
    if ($file === false) {
        return '<p>問題データを読み込めませんでした。</p>';
    }

    $header = fgetcsv($file, 0, ',', '"', '');
    $expected_header = ['Label', 'Year', 'Month', 'Subject ', 'No', 'QuestionState', 'Answer', 'Figure', 'Option1', 'Option2', 'Option3', 'Option4'];
    if ($header !== false) {
        $header[0] = preg_replace('/^\xEF\xBB\xBF/', '', (string) $header[0]);
    }
    if ($header !== $expected_header) {
        fclose($file);
        return '<p>問題データのヘッダーが正しくありません。</p>';
    }

    $questions = [];
    while (($row = fgetcsv($file, 0, ',', '"', '')) !== false) {
        if ($row === [null]) {
            continue;
        }
        if (count($row) !== 12 || $row[0] === '' || !in_array($row[6], ['1', '2', '3', '4'], true)) {
            fclose($file);
            return '<p>問題データに不正な行があります。</p>';
        }

        $questions[] = [
            'label' => $row[0],
            'subject' => $row[3],
            'number' => $row[4],
            'statement' => $row[5],
            'answer' => (int) $row[6],
            'figure' => $row[7],
            'options' => array_slice($row, 8, 4),
        ];
    }
    fclose($file);

    // HTMLのscript要素内に安全に埋め込む（CSVにタグが含まれても実行させない）。
    $json = wp_json_encode($questions, JSON_HEX_TAG | JSON_UNESCAPED_UNICODE);
    if (!$questions || $json === false) {
        return '<p>有効な問題データがありません。</p>';
    }

    // ショートコードを使用したページだけで読み込む。
    wp_enqueue_style('ppl-hcg-quiz', plugins_url('style.css', __FILE__), [], '1.1.0');
    wp_enqueue_script('ppl-hcg-quiz', plugins_url('question.js', __FILE__), [], '1.1.0', true);

    ob_start();
    ?>
    <div class="ppl-hcg-quiz" data-figure-url="<?php echo esc_url(plugins_url('fig/', __FILE__)); ?>">
        <label class="ppl-hcg-quiz-range">
            出題範囲を選択してください:
            <select class="ppl-hcg-quiz-select">
                <option value="" data-mode="all">全科目シャッフル</option>
                <option value="航空気象(P22)" data-mode="subject">航空気象(P22)シャッフル</option>
                <option value="航空工学(P26)" data-mode="subject">航空工学(P26)シャッフル</option>
                <option value="空中航法(P49)" data-mode="subject">空中航法(P49)シャッフル</option>
                <option value="航空法規(P28)" data-mode="subject">航空法規(P28)シャッフル</option>
            </select>
        </label>
        <noscript><p>クイズを利用するにはJavaScriptを有効にしてください。</p></noscript>
        <p class="ppl-hcg-quiz-message" role="status" hidden></p>
        <div class="ppl-hcg-quiz-panel" hidden>
            <p class="ppl-hcg-quiz-label"></p>
            <p class="ppl-hcg-quiz-number"></p>
            <p class="ppl-hcg-quiz-statement"></p>
            <p class="ppl-hcg-quiz-figure" hidden></p>
            <div class="ppl-hcg-quiz-options">
                <?php for ($choice = 1; $choice <= 4; $choice++): ?>
                    <button type="button" data-choice="<?php echo $choice; ?>"></button>
                <?php endfor; ?>
            </div>
            <div class="ppl-hcg-quiz-actions">
                <button type="button" class="ppl-hcg-quiz-next" disabled>次の問題へ</button>
                <button type="button" data-choice="0">分からない</button>
            </div>
            <p class="ppl-hcg-quiz-answer" role="status"></p>
        </div>
        <script type="application/json" class="ppl-hcg-quiz-data"><?php echo $json; ?></script>
    </div>
    <?php
    return ob_get_clean();
}
