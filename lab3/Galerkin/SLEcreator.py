# Вычисление коэфициентов для базиса из экспоненциальных функций: x^i
expKoefFactor = lambda i, j: 2 * (1 - i) / (i + j + 1) - i * j / (i + j - 1)

# Вычисление результирующего вектора для базиса из экспоненциальных функций: x^i
expRes = lambda j: 1 / (j + 2) - 1

# Вычисление приближенного значения для коэффициентов для экспоненциальной функции: x^i
expBasis = lambda x, koef: sum(koef[i - 1] * x**i for i in range(1, len(koef) + 1))

def createSLE(N: int, koef_func: callable, res_func: callable) -> tuple[list[list[float]], list[float]]:
    """Создать систему линейных уравнений"""
    return (
        [[koef_func(i, j) for i in range(1, N + 1)] for j in range(1, N + 1)] + [[i for i in range(1, N + 1)]],
        [res_func(j) for j in range(1, N + 1)] + [1]
    )