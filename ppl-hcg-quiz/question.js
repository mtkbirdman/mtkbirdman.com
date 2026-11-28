// 各ショートコードのクイズを独立して初期化する。
// WordPressがフッターでこのスクリプトを読み込むため、HTMLは生成済み。
for (const root of document.querySelectorAll('.ppl-hcg-quiz')) {
    const questions = JSON.parse(root.querySelector('.ppl-hcg-quiz-data').textContent);
    const select = root.querySelector('.ppl-hcg-quiz-select');
    const message = root.querySelector('.ppl-hcg-quiz-message');
    const panel = root.querySelector('.ppl-hcg-quiz-panel');
    const label = root.querySelector('.ppl-hcg-quiz-label');
    const number = root.querySelector('.ppl-hcg-quiz-number');
    const statement = root.querySelector('.ppl-hcg-quiz-statement');
    const figure = root.querySelector('.ppl-hcg-quiz-figure');
    const optionButtons = [...root.querySelectorAll('.ppl-hcg-quiz-options button')];
    const answerButtons = [...root.querySelectorAll('button[data-choice]')];
    const nextButton = root.querySelector('.ppl-hcg-quiz-next');
    const answer = root.querySelector('.ppl-hcg-quiz-answer');

    let selectedQuestions = [];
    let index = 0;

    // 年月・科目別の出題範囲をCSVの登場順で追加する。
    for (const name of new Set(questions.map(question => question.label))) {
        const option = document.createElement('option');
        option.value = name;
        option.dataset.mode = 'label';
        option.textContent = name;
        select.append(option);
    }

    // 出題範囲の選択、重複除去、シャッフル、最初の問題の表示。
    function startQuiz() {
        const mode = select.selectedOptions[0].dataset.mode;
        const value = select.value;
        selectedQuestions = questions.filter(question =>
            mode === 'all' ||
            (mode === 'subject' && question.subject === value) ||
            (mode === 'label' && question.label === value)
        );

        // 年月・科目別はCSVの順番と全問題を維持する。
        if (mode !== 'label') {
            // 既存仕様と同じく、問題文・正解番号が一致する問題を重複とする。
            const seen = new Set();
            selectedQuestions = selectedQuestions.filter(question => {
                const key = JSON.stringify([question.statement, question.answer]);
                if (seen.has(key)) return false;
                seen.add(key);
                return true;
            });

            // Fisher–Yates法：残っている要素から均等に1つ選ぶ。
            for (let i = selectedQuestions.length - 1; i > 0; i--) {
                const j = Math.floor(Math.random() * (i + 1));
                [selectedQuestions[i], selectedQuestions[j]] = [selectedQuestions[j], selectedQuestions[i]];
            }
        }

        index = 0;
        panel.hidden = selectedQuestions.length === 0;
        message.hidden = !panel.hidden;
        message.textContent = panel.hidden ? '選択した範囲に問題がありません。' : '';
        if (!panel.hidden) showQuestion();
    }

    // 問題を表示し、前の問題の採点結果をすべてリセットする。
    function showQuestion() {
        const question = selectedQuestions[index];
        label.textContent = question.label;
        number.textContent = question.number;
        statement.textContent = question.statement;
        optionButtons.forEach((button, i) => {
            button.textContent = question.options[i];
        });

        figure.replaceChildren();
        figure.hidden = question.figure === '-';
        if (!figure.hidden) {
            const image = document.createElement('img');
            image.alt = '問題の図';
            image.addEventListener('error', () => {
                image.replaceWith('問題の図を読み込めませんでした。');
            }, { once: true });
            image.src = root.dataset.figureUrl + encodeURIComponent(question.figure);
            figure.append(image);
        }

        for (const button of answerButtons) {
            button.classList.remove('is-correct', 'is-incorrect');
            button.removeAttribute('aria-disabled');
        }
        answer.textContent = '';
        nextButton.disabled = true;
    }

    select.addEventListener('change', startQuiz);

    // 選択肢・「分からない」を共通の採点処理で扱う。
    for (const button of answerButtons) {
        button.addEventListener('click', () => {
            if (!nextButton.disabled) return; // 採点は1問につき1回だけ。

            const correctAnswer = selectedQuestions[index].answer;
            const isCorrect = Number(button.dataset.choice) === correctAnswer;
            if (!isCorrect) button.classList.add('is-incorrect');
            optionButtons[correctAnswer - 1].classList.add('is-correct');
            answer.textContent = isCorrect ? '正解' : '不正解';

            for (const answerButton of answerButtons) {
                answerButton.setAttribute('aria-disabled', 'true');
            }
            nextButton.disabled = false;
        });
    }

    nextButton.addEventListener('click', () => {
        index = (index + 1) % selectedQuestions.length;
        showQuestion();
    });

    startQuiz();
}
