from math import exp
from RadiationTransfer.absoption import logInterpKT

# Температурное поле
TempFieldFunc = lambda x, T0, Tw, R, p: (Tw - T0) * (x / R)**p + T0

# Функция Планка
PlanckFunc = lambda x, T0, Tw, R, p: 0.0003084 / (exp(47990 / TempFieldFunc(x, T0, Tw, R, p)) - 1)

# Функция уточнения коэффициента поглощения
RectificateAbsorptionCoefFunc = lambda x, T, k, T0, Tw, R, p: logInterpKT(T, k, TempFieldFunc(x, T0, Tw, R, p))


def solute(
        N: int,
        T0: float, Tw: float,
        R: float, p: float,
        T: list[float], k: list[float]) -> tuple[list[float], list[float]]:
    """Решение ОДУ"""

    h = R / N
    r = [i / N * R for i in range(0, N + 1)]
    u_p = [PlanckFunc(ri, T0, Tw, R, p) for ri in r]
    kinterp = [RectificateAbsorptionCoefFunc(ri, T, k, T0, Tw, R, p) for ri in r]

    A = [0] * (N+1)
    B = [0] * (N+1)
    C = [0] * (N+1)
    F = [0] * (N+1)

    rmp = lambda i: r[i] + h / 2
    rmm = lambda i: r[i] - h / 2

    kmp = lambda i: (kinterp[i] + kinterp[i + 1]) / 2
    kmm = lambda i: (kinterp[i - 1] + kinterp[i]) / 2

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
        A[i] = Af(i)
        B[i] = Bf(i)
        C[i] = Cf(i)
        F[i] = Ff(i)
    
    # Граничное условие при r=R (i=N)
    A[N] = 1.0
    B[N] = -4.0
    C[N] = 3 + 2.34 * h * kinterp[N]
    F[N] = 0.0

    # Метод прогонки
    u = [0] * (N+1)
    zita = [0] * (N+1)
    ita = [0] * (N+1)

    # Граничное условие при i=0 для u[i], r=0 (центр)
    zita[0] = 1.0
    ita[0] = 0.0

    # Прямой ход
    for i in range(1, N+1):
        denom = B[i] + A[i] * zita[i-1]
        zita[i] = -C[i] / denom
        ita[i] = (F[i] - A[i] * ita[i-1]) / denom

    # Граничное условие при i=N для u[i], r=R (стенка)
    denom = 3 + 2.34 * h * kinterp[N] + 4 * zita[N] - zita[N-1] * zita[N]
    u[N] = (4 * ita[N] - zita[N-1] * ita[N] - ita[N-1]) / denom
    
    # Обратный ход
    for i in range(N-1, 0, -1):
        u[i] = zita[i] * u[i+1] + ita[i]
    
    # Граничное условие при i=N для u[i], r=R (стенка): u_0 = u_1
    u[0] = u[1]
    
    return r, u