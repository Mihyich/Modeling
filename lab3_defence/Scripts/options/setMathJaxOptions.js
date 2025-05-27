window.MathJax = {
    tex: {
      inlineMath: [['$', '$'], ['\\(', '\\)']],
      processEscapes: true
    },
    options: {
      enableMenu: false
    }
};

document.querySelectorAll('.math-container').forEach(el => {
    el.style.cursor = 'pointer';
    el.title = 'Нажмите, чтобы скопировать LaTeX';

    el.addEventListener('click', () => {
            const text = el.textContent.trim();
            navigator.clipboard.writeText(text).then(() => {
            alert('Скопировано: ' + text);
        }).catch(err => { console.error('Ошибка копирования: ', err); });
    });
});