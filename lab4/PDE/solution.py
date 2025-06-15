from Core.defines import c
from Core.utility import TempField, PulseCurrent, ElectricFieldStrenght, interp, logInterp, linearInterp, soluteTriSys
from Core.integrator import integrateSimpson
from Core.plotter import Plotter
from ODE.solution import soluteODE


def crank_nicolson_step(
    T_prev: list[float],
    r: list[float],
    h: float,
    tau: float,
    T0: float,
    Tw: float,
    R: float,
    p: float,
    I_max: float,
    t_max: float,
    T_tableODE: list[float],
    K_table: list[float],
    T_tablePDE: list[float],
    sigma_table: list[float],
    lambda_table: list[float],
    cT_table: list[float],
    plotter: Plotter,
    s: int,
    eps: float = 1e-7,
    max_iter: int = 100
) -> list[float]:
    """Обновляет температурное поле методом Кранка–Николсона с итерациями."""

    N = len(r) - 1

    # Скопировать предыдущий слой
    T_new = T_prev.copy()
    T_old = T_new.copy()

    # Текущий временной слой
    cur_tau = s * tau

    # Пересчет импулься тока
    I = PulseCurrent(cur_tau, t_max, I_max)

    for _ in range(max_iter):
        # 1. Обновление коэффициентов
        lambdas = [interp(T_tablePDE, lambda_table, T_new[i], logInterp) for i in range(N+1)]
        sigmas = [interp(T_tablePDE, sigma_table, T_new[i], linearInterp) for i in range(N+1)]
        cTs = [interp(T_tablePDE, cT_table, T_new[i], linearInterp) for i in range(N+1)]

        # 2. Решение ODE для u(r) с текущим T_new
        _, kinterp, u, u_p = soluteODE(N, T0, Tw, R, p, T_new, T_tableODE, K_table)

        # 3. Пересчёт источника тепла q(r)
        q = [c * kinterp[i] * (u_p[i] - u[i]) for i in range(N)]

        # 4. Расчёт напряжённости электрического поля E(t)
        E = ElectricFieldStrenght(
            r, T_new, R, I,
            T_tablePDE, sigma_table,
            integrateSimpson
        )

        kappas = [(lambdas[i] + lambdas[i+1]) / 2 for i in range(N)]
        rh = [(r[i] + r[i+1]) / 2 for i in range(N)]

        # 5. Пересчёт коэффициентов разностной схемы на основе текущего приближения T_new
        A = [0.0] * (N+1)
        B = [0.0] * (N+1)
        C = [0.0] * (N+1)
        F = [0.0] * (N+1)

        # Внутренние точки
        for i in range(1, N):
            A[i] = rh[i-1] * (kappas[i-1] / h) * tau

            C[i] = rh[i] * (kappas[i] / h) * tau

            B[i] = A[i] + C[i] + cTs[i] * r[i] * h

            F[i] = (cTs[i] * T_prev[i] * r[i] * h +
                    sigmas[i] * (E ** 2) * r[i] * h * tau -
                    q[i] * r[i] * h * tau)

        # 6. Метод прогонки
        T_new = soluteTriSys(A, B, C, F)

        # 7. Проверка сходимости
        error = max(abs((T_new[i] - T_old[i]) / T_new[i]) if T_new[i] != 0 else float('inf') for i in range(len(T_new)))

        if error < eps:
            break

        T_old = T_new.copy()

    # plotter.plot_q(r, q, s)
    # plotter.plot_u_and_up(r, u, u_p, s)

    return T_new, E


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
    cT_table: list[float],
    plotter: Plotter
) -> tuple[list[float], list[list[float]]]:
    """Решение ДУЧП"""

    h = R / N
    r = [i * h for i in range(N+1)]
    Tprev = [TempField(ri, T0, Tw, R, p) for ri in r]
    Thistory = [Tprev.copy()]

    Evalues = []

    # Цикл по времени
    for s in range(S + 1):
        print(f"Шаг по времени {s}/{S}")

        # 1. Метод итераций
        T_new, E = crank_nicolson_step(Tprev, r, h, tau, T0, Tw, R, p, I_max, t_max, T_tableODE, K_table, T_tablePDE, sigma_table, lambda_table, cT_table, plotter, s)
        Evalues += [E]

        # 2. Сохранить результат
        Thistory.append(T_new.copy())
        Tprev = T_new.copy()

    times = [s * tau for s in range(S+1)]
    plotter.plot_E_over_time(times, Evalues)

    # print(Evalues)

    return r, Thistory
