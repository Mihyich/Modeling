from math import exp


# Температурное поле
TempFieldFunc = lambda x, T0, Tw, R, p: (Tw - T0) * (x / R)**p + T0

# Функция Планка
PlanckFunc = lambda x, T0, Tw, R, p: 0.0003084 / (exp(47990 / TempFieldFunc(x, T0, Tw, R, p)) - 1)


def solute(
        N: int,
        T0: float, Tw: float,
        R: float, p: float,
        T: list[float], k: list[float],
        upFunc: callable) -> tuple[list[float], list[float]]:
    """Решение ОДУ"""

    h = R / N
    r = [i / N * R for i in range(0, N + 1)]
    u_p = [upFunc(ri) for ri in r]

    A = [0] * (N+1)
    B = [0] * (N+1)
    C = [0] * (N+1)
    F = [0] * (N+1)
    u = [0] * (N+1)

    rmp = lambda i: r[i] + h / 2
    rmm = lambda i: r[i] - h / 2

    kmp = lambda i: (k[i] + k[i + 1]) / 2
    kmm = lambda i: (k[i - 1] + k[i]) / 2

    alphamp = lambda i: rmp(i) / (h * kmp(i))
    alphamm = lambda i: rmm(i) / (h * kmm(i))

    betamp = lambda i: 3 * h * rmp(i) * kmp(i) / 4
    betamm = lambda i: 3 * h * rmm(i) * kmm(i) / 4

    Af = lambda i: alphamm(i) - betamm(i)
    Bf = lambda i: -alphamp(i) - alphamm(i) - betamm(i) - betamp(i)
    Cf = lambda i: alphamp(i) - betamp(i)
    Ff = lambda i: -betamm(i) * u_p[i - 1] - (betamm(i) + betamp(i)) * u_p[i] - betamp(i) * u_p[i + 1]

    # Граничное условие в центре (i=0)
    A[0] = 0.0
    B[0] = 1.0
    C[0] = -1.0  # u[0] = u[1]
    F[0] = 0.0

    # Внутренние точки (i = 1 ... N-1)
    for i in range(1, N):
        A[i] = Af[i]
        B[i] = Bf[i]
        C[i] = Cf[i]
        F[i] = Ff[i]
    
    # Граничное условие при r=R (i=N)
    A[N] = 1.0
    B[N] = -4.0
    C[N] = 3 + 2.34 * h * k[N]
    F[N] = 0.0

    # Метод прогонки
    Alpha = [0] * (N+1)
    Beta = [0] * (N+1)

    # Прямой ход
    Alpha[1] = -C[1] / B[1]
    Beta[1] = F[1] / B[1]
    for i in range(2, N+1):
        denom = B[i] + A[i] * Alpha[i-1]
        Alpha[i] = -C[i] / denom if i < N else 0.0
        Beta[i] = (F[i] - A[i] * Beta[i-1]) / denom
    
    # Обратный ход
    u[N] = (4 * Beta[N-1] - Beta[N-2]) / (3 + 2.34 * h * k[N] + 4 * Alpha[N-1] - Alpha[N-2])
    for i in range(N-1, -1, -1):
        u[i] = Alpha[i] * u[i+1] + Beta[i] if i > 0 else u[1]
    
    return r, u


def main():
    pass


if __name__ == "__main__":
    main()