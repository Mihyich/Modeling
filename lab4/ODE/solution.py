from Core.utility import interp, linearInterp, Plank, soluteTriSys


def soluteODE(
    N: int,
    T0: float, Tw: float,
    R: float, p: float,
    T: list[float],
    T_table: list[float],
    K_table: list[float]
) -> tuple[list[float], list[float], list[float], list[float]]:
    """Решение ОДУ"""

    h = R / N
    r = [i * h for i in range(N+1)]
    # u_p = [Plank(ri, T0, Tw, R, p) for ri in r]
    u_p = [Plank(ti, T0, Tw, R, p) for ti in T]
    kinterp = [interp(T_table, K_table, ti, linearInterp) for ti in T]

    # Коэффициенты системы
    A = [1 / k for k in kinterp]
    B = [3 * k for k in kinterp]
    C = [3 * k * up for k, up in zip(kinterp, u_p)]
    F = [2 * (A[i] * A[i+1]) / (A[i] + A[i+1]) for i in range(N)]

    rh = [(r[i] + r[i+1]) / 2 for i in range(N)]

    y = [0.0] * (N+1)
    zita = [0.0] * (N+1)
    ita = [0.0] * (N+1)

    # Левое граничное условие (при r=0)
    zita[0] = 1.0
    ita[0] = 0.0

    # Прямой ход - вычисление прогоночных коэффициентов
    for i in range(1, N):
        rhl = rh[i-1]
        rhr = rh[i]
        V = (rhr**2 - rhl**2) / 2
        
        a_i = rhl * F[i-1] / (R**2 * h)
        d_i = rhr * F[i] / (R**2 * h)
        b_i = a_i + d_i + B[i] * V
        f_i = C[i] * V
        
        denominator = b_i - a_i * zita[i-1]
        zita[i] = d_i / denominator
        ita[i] = (a_i * ita[i-1] + f_i) / denominator

    # Правое граничное условие (при z=R)
    Mn = (rh[N-1] * F[N-1]) / (R**2 * h)
    Kn = -(r[N] * 3 * 0.39 / R + (rh[N-1] * F[N-1]) / (R**2 * h) + B[N] * r[N] * h / 2)
    Qn = -C[N] * r[N] * h / 2
    
    y[N] = (Qn - Mn * ita[N-1]) / (Mn * zita[N-1] + Kn)

    # Обратный ход - вычисление решения
    for i in range(N-1, -1, -1):
        y[i] = zita[i] * y[i+1] + ita[i]

    # try:
    #     u = soluteTriSys(A, B, C, F)
    # except ZeroDivisionError as e:
    #     print(f"[Ошибка]: {e}")
    #     u = [0.0] * (N+1)

    return r, kinterp, y, u_p