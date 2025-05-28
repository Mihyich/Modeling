def integrateSimpson(
        a: float, b: float,
        f: callable,
        n: int = 100
) -> float:
    """Численное интегрирование методом Симпсона"""
    
    N = n + 1 if n % 2 else n
    h = (b - a) / N
    sum = f(a) + f(b)

    for i in range(1, N):
        x = a + i * h
        factor = 4 if i % 2 else 2
        sum += factor * f(x)

    return h * sum / 3