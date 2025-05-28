from Core.defines import c
from Core.utility import TempField, ElectricFieldStrenght, logInterpXT
from Core.integrator import integrateSimpson
from ODE.solution import soluteODE

def crank_nicolson_step(
    T_prev: list[float],
    r: list[float],
    h: float,
    tau: float,
    q: list[float],
    E: float,
    T_tablePDE: list[float],
    sigma_table: list[float],
    lambda_table: list[float],
    cT_table: list[float],
    eps: float = 1e-5,
    max_iter: int = 100
) -> list[float]:
    """
    Обновляет температурное поле методом Кранка–Николсона с итерациями.
    
    Параметры:
        T_prev       - предыдущее распределение температуры
        r            - радиальная сетка
        h            - шаг по пространству
        tau          - шаг по времени
        q            - источник тепла (q(r))
        sigma        - электропроводность σ(T)
        lambd        - теплопроводность λ(T)
        c_T_list     - теплоёмкость c_T(T)
        eps          - точность итераций
        max_iter     - максимальное число итераций
    
    Возвращает:
        T_new - новое распределение температуры
    """
    N = len(r) - 1
    T_new = T_prev.copy()
    T_old = T_new.copy()

    # Вспомогательные функции
    def r_half(i, direction):
        return (r[i] + r[i + direction]) / 2
    
    rmp = lambda i: r[i] + h / 2
    rmm = lambda i: r[i] - h / 2

    for _ in range(max_iter):
        # обновление коэффициентов
        lambdas = [logInterpXT(T_tablePDE, lambda_table, T_new[i]) for i in range(N+1)]
        sigmas = [logInterpXT(T_tablePDE, sigma_table, T_new[i]) for i in range(N+1)]
        cTs = [logInterpXT(T_tablePDE, cT_table, T_new[i]) for i in range(N+1)]

        # Пересчёт коэффициентов на основе текущего приближения T_new
        A = [0.0] * (N + 1)
        B = [0.0] * (N + 1)
        C = [0.0] * (N + 1)
        F = [0.0] * (N + 1)

        # Внутренние точки
        for i in range(1, N):
            ri = r[i]
            rhp = rmp(i)
            rhm = rmm(i)

            lhp = (lambdas[i] + lambdas[i+1]) / 2
            lhm = (lambdas[i-1] + lambdas[i]) / 2

            ap = rhp * lhp
            am = rhm * lhm

            # Коэффициенты неявной части
            A[i] = (tau / (2 * cTs[i] * ri * h**2)) * am
            B[i] = 1 + (tau / (2 * cTs[i] * ri * h**2)) * (ap + am)
            C[i] = -(tau / (2 * cTs[i] * ri * h**2)) * ap

            # Коэффициенты явной части
            As = (tau / (2 * cTs[i] * ri * h**2)) * am
            Bs = 1 - (tau / (2 * cTs[i] * ri * h**2)) * (ap + am)
            Cs = (tau / (2 * cTs[i] * ri * h**2)) * ap

            # Правая часть
            F[i] = (
                As * T_prev[i - 1] +
                Bs * T_prev[i] +
                Cs * T_prev[i + 1] +
                (tau / cTs[i]) * (sigmas[i] * E**2 - q[i])
            )
        
        # Краевое условие в центре цилиндра
        A[0] = 0.0
        B[0] = 1.0
        C[0] = 0.0
        F[0] = T_prev[0]

        # Краевое условие на стенке
        rhm = rmm(N)
        lhm = (lambdas[N-1] + lambdas[N]) / 2
        am = rhm * lhm

        A[N] = (tau / (2 * cTs[N] * ri * h**2)) * am
        B[N] = 1.0
        C[N] = -(tau / (2 * cTs[N] * ri * h**2)) * am

        As = A[N]
        Bs = B[N]
        Cs = C[N]

        F[N] = (
            As * T_prev[N - 1] +
            Bs * T_prev[N] +
            Cs * T_prev[N] +
            (tau / cTs[N]) * (sigmas[i] * E**2 - q[N])
        )

        # Метод прогонки
        P = [0.0] * (N + 1)
        Q = [0.0] * (N + 1)

        P[0] = -C[0] / B[0]
        Q[0] = F[0] / B[0]

        for i in range(1, N + 1):
            denom = B[i] + A[i] * P[i - 1]
            P[i] = -C[i] / denom
            Q[i] = (F[i] - A[i] * Q[i - 1]) / denom

        # Обратный ход
        T_new[N] = Q[N]
        for i in range(N - 1, -1, -1):
            T_new[i] = P[i] * T_new[i + 1] + Q[i]

        # Условие симметрии в центре
        T_new[0] = T_new[1]

        # Проверка сходимости
        error = max(abs((T_new[i] - T_old[i]) / T_new[i]) if T_new[i] != 0 else float('inf') for i in range(N+1))

        if error < eps:
            break

        T_old = T_new.copy()

    return T_new


def solutePDE(
    N: int,
    R: float,
    T0: float,
    Tw: float,
    p: float,
    tau: float,
    I_max: float,
    t_max: float,
    S: int,
    T_tableODE: list[float],
    K_table: list[float],
    T_tablePDE: list[float],
    sigma_table: list[float],
    lambda_table: list[float],
    cT_table: list[float]
) -> tuple[list[float], list[list[float]]]:
    """Решение ДУЧП"""

    h = R / N
    r = [n * h for n in range(N+1)]
    Tprev = [TempField(ri, T0, Tw, R, p) for ri in r]
    Thistory = [Tprev.copy()]

    for s in range(S + 1):
        # print(f"Шаг по времени {s}/{S}")

        t_current = s * tau

        # 1. Решение u(r)
        _, kinterp, u, u_p = soluteODE(N, T0, Tw, R, p, Tprev, T_tableODE, K_table)

        # 2. Вычисление источника тепла q(r)
        q = [c * kinterp[i] * (u_p[i] - u[i]) for i in range(N+1)]

        # 3. Вычисление напряжённости электрического поля E(t)
        E = ElectricFieldStrenght(
            t_current,
            I_max,
            t_max,
            T0=T0,
            Tw=Tw,
            R=R,
            p=p,
            T_table=T_tablePDE,
            sigma_table=sigma_table,
            integrate=integrateSimpson
        )

        # 4. Метод итераций
        T_new = crank_nicolson_step(Tprev, r, h, tau, q, E, T_tablePDE, sigma_table, lambda_table, cT_table)

        # 5. Сохранить результат
        Thistory.append(T_new.copy())
        Tprev = T_new.copy()

    return r, Thistory
