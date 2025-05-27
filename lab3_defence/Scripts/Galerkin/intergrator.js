/**
 * @param {number} n
 * N = 2n число элементарных отрезков одинаковой длины
 */
export function integrateSimpson(f, a, b, n) {
    const N = n % 2 ? n + 1 : n;
    const h = (b - a) / N;
    let sum = f(a) + f(b);
    let x;
    let factor;

    for (let i = 1; i < N; ++i) {
        x = a + i * h;
        factor = i % 2 ? 4 : 2;
        sum += factor * f(x);
    }

    return h * sum / 3;
}