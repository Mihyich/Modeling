import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

# Константы
c = 3e10  # скорость света, см/с
R = 0.35  # радиус цилиндра, см
Tw = 2000  # температура стенки, K
T0 = 10000  # начальная температура, K
dr = 0.01  # шаг по r
p = 4

# Таблица значений k(T)
t_values = np.array([2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
k_values = np.array(
    [8.200e-3, 2.768e-2, 6.560e-2, 1.281e-1, 2.214e-1, 3.516e-1, 5.248e-1, 7.472e-1, 1.025e0]
)

# Интерполяция k(T) в логарифмическом масштабе
log_k_interp = interp1d(
    np.log(t_values), np.log(k_values), kind="linear", fill_value="extrapolate"
)


def k_interp(T):
    """Вычисляет интерполированное значение k(T) с использованием логарифмической интерполяции."""
    return np.exp(log_k_interp(np.log(T)))


# Функция температуры
# T = lambda r: T0 * (1 - (r / R)) + Tw * (r / R)

T = lambda r: (Tw - T0) * (r / R)**p + T0

# Функция Планка
up = lambda r: 3.084e-4 / (np.exp(4.799e4 / T(r)) - 1)

# Система уравнений
# def system(y, r, k, up_r):
#     u, F = y
#     # return [F - (3 * k / c) * (u - up_r), -u + (F * dr / (3 * k))]
#     # return [-F * (3 * k / c) * (u - up_r), -F / r - c * k * (u - up_r)]
#     # return [-F / r - c * k * (u - up_r), -F * (3 * k / c)]

# # Метод трапеций
# def trapezoidal_method(R, dr, tol=1e-6):
#     r = np.arange(0, R + dr, dr)
#     y = np.zeros((len(r), 2))
#     y[0] = [up(0), 0]  # начальные условия

#     for i in range(1, len(r)):
#         k_i = k_interp(T(r[i-1]))
#         k_i1 = k_interp(T(r[i]))
#         up_r_i = up(r[i])
#         up_r_i1 = up(r[i-1])

#         func = lambda y_i: system(y_i, r[i], (k_i + k_i1) / 2, (up_r_i + up_r_i1) / 2)
#         y[i] = fsolve(func, y[i-1])

#     return r, y


def trapezoidal_method(R, dr, u0):
    r = np.arange(0, R + dr, dr)
    y = np.zeros((len(r), 2))
    y[0] = [0, u0]
    epsilon = 1e-6
    
    for i in range(1, len(r)):
        k_i = k_interp(T(r[i - 1]))
        k_i1 = k_interp(T(r[i]))
        h = (k_i + k_i1) / 2
        up_r_i = up(r[i - 1])
        up_r_i1 = up(r[i])

        denominator = (
            4 * r[i - 1] * r[i]
            + 2 * h * r[i - 1]
            - 3 * h**2 * k_i * k_i1 * r[i - 1] * r[i] * up_r_i1
        )

        if denominator == 0.0:
            print(f"Warning: denominator is zero at r = {r[i]} with k_i = {k_i}, k_i1 = {k_i1}, up_r_i = {up_r_i}, up_r_i1 = {up_r_i1}")
            denominator = epsilon
            # y[i][0] = y[i - 1][0]
            # y[i][1] = y[i - 1][1]
            # continue

        y[i][1] = (
            4 * r[i - 1] * r[i] * y[i - 1][1]
            - 2 * h * r[i] * y[i - 1][1]
            - 2 * r[i - 1] * r[i] * h * c * k_i * y[i - 1][0] * up_r_i
            - 2 * c * y[i - 1][0] * r[i - 1] * r[i] * h * k_i1 * up_r_i1
            + 3 * h**2 * k_i * y[i - 1][1] * r[i - 1] * r[i] * k_i1 * up_r_i1
        ) / denominator

        y[i][0] = (
            2 * c * y[i - 1][0] - 3 * h * k_i * y[i - 1][1] - 3 * h * k_i1 * y[i][1]
        ) / (2 * c)

        # Отладочные выводы
        print(f"Step {i}: r = {r[i]}, k_i = {k_i}, k_i1 = {k_i1}, up_r_i = {up_r_i}, up_r_i1 = {up_r_i1}, y[i][0] = {y[i][0]}, y[i][1] = {y[i][1]}")
        print(f"Denominator: {denominator}")

    return r, y



# Метод стрельбы
def shooting_method(tol=1e-8, max_iter=200):  # Увеличение точности и числа итераций
    chi_min, chi_max = 0.00, 1.0  # Начальные условия
    for iteration in range(max_iter):
        chi_mid = (chi_min + chi_max) / 2
        u0 = chi_mid * up(0)
        r, y = trapezoidal_method(R, dr, u0)
        u_R = y[-1, 0]
        F_R = y[-1, 1]

        # Отладочная информация
        print(f"Iteration: {iteration}, chi_mid: {chi_mid:.10f}, F_R: {F_R:.10e}, u_R: {u_R:.10e}")

        # Условие для завершения
        if abs(-F_R + 0.39 * u_R) < tol:
            print(f"Оптимальное значение chi найдено: {chi_mid:.10f}")
            return chi_mid
        elif -F_R + 0.39 * u_R > 0:
            chi_max = chi_mid  # Уменьшаем верхнюю границу
        else:
            chi_min = chi_mid  # Увеличиваем нижнюю границу

    print("Не удалось найти оптимальное значение chi.")
    return None  # Возврат None, если решение не найдено


# Поиск оптимального значения chi
chi_opt = shooting_method()

if chi_opt is not None:
    u0_opt = chi_opt * up(0)
    r_vals, u_vals = trapezoidal_method(R, dr, u0_opt)

    # Построение графиков
    plt.figure(figsize=(12, 8))


    plt.subplot(2, 1, 1)
    plt.plot(r_vals, [v + 0.000000043 for v in u_vals[:, 0]], "b-", linewidth=2)
    plt.plot(r_vals, up(r_vals), label=f'up(r)', color='red', linestyle='--')
    plt.title(f"Решение u(r) при chi = {chi_opt:.6f}")
    plt.xlabel("r")
    plt.ylabel("u(r)")
    plt.grid(True)

    # plt.subplot(2, 1, 2)
    # plt.plot(r_vals, u_vals[:, 1], "r-", linewidth=2)
    # plt.title("Функция F(r)")
    # plt.xlabel("r")
    # plt.ylabel("F(r)")
    # plt.grid(True)

    plt.tight_layout()
    plt.show()

    # Вывод оптимального значения chi
    print(f"Оптимальное значение chi: {chi_opt}")
else:
    print("Не удалось найти оптимальное значение chi")
