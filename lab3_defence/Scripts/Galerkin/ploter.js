export function createPlot(xTitle, yTitle) {
    if (window.solutionChart instanceof Chart)
        return;

    const canvas = document.getElementById("solutionChart");
    if (!canvas) {
        console.warn("Canvas не найден");
        return;
    }

    const ctx = canvas.getContext("2d");
    if (!ctx) {
        console.error("Не удалось получить 2D-контекст canvas");
        return;
    }

    window.solutionChart = new Chart(
        ctx,
        {
            type: "line",
            options:
            {
                responsive: true,
                scales:
                {
                    x: { title: { display: true, text: xTitle } },
                    y: { title: { display: true, text: yTitle } }
                },
                plugins:
                {
                    legend: { display: true },
                    tooltip: { enabled: true }
                }
            }
        }
    );
}

export function clearPlot() {
    if (window.solutionChart instanceof Chart) {
        window.solutionChart.data.datasets = [];
        window.solutionChart.update();
    }
}

export function plotSolution(f, a, b, n, label, color) {
    const xs = [];
    const ys = [];

    for (let i = 0; i <= n; ++i) {
        const x = a + (b - a) * (i / n);
        const y = f(x);
        xs.push(x.toFixed(3));
        ys.push(isNaN(y) || !isFinite(y) ? null : y);
    }

    // Если график уже создан — обновляем его
    if (window.solutionChart instanceof Chart) {
        // Обновляем метки (x)
        window.solutionChart.data.labels = xs;

        // Добавляем новую линию как новый датасет
        window.solutionChart.data.datasets.push({
            label: label,
            data: ys,
            borderColor: color,
            fill: false,
            tension: 0.1
        });

        window.solutionChart.update();
    }
}