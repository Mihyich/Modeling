import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp, mpf, exp, log, power, pi, matrix

# Устанавливаем точность в 50 знаков
mp.dps = 50

# Конвертируем все константы в mpf (mpmath floating point)
c = mpf('3e10')  # скорость света, см/с
R = mpf('0.35')  # радиус цилиндра, см
Tw = mpf('2000')  # температура стенки, K
T0 = mpf('10000')  # начальная температура, K
P = mpf('4')
EPS = mpf('1e-50')  # точность метода стрельбы (увеличена до 50 знаков)

# Данные для интерполяции (храним как обычные числа)
t_values = np.array([2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
k_values_1 = np.array([8.2E-03, 2.768E-02, 6.56E-02, 1.281E-01, 2.214E-01,
                       3.516E-01, 5.248E-01, 7.472E-01, 1.025E+00])
k_values_2 = np.array([1.600E+00, 5.400E+00, 1.280E+01, 2.500E+01, 4.320E+01,
                       6.860E+01, 1.024E+02, 1.458E+02, 2.000E+02])
k_values = k_values_1

# Логарифмические значения для интерполяции
log_t = np.log(t_values)
log_k = np.log(k_values)


def T(r):
    return T0 + (Tw - T0) * power(r / R, P)


def up(r):
    return mpf('3.084e-4') / (exp(mpf('4.799e4') / T(r)) - mpf('1'))


def k_interp(T_val):
    # Конвертируем mpf в float для интерполяции
    T_float = float(T_val)
    if T_float <= t_values[0]:
        return mpf(k_values[0])
    elif T_float >= t_values[-1]:
        return mpf(k_values[-1])
    else:
        # Линейная интерполяция в логарифмическом пространстве
        log_T = np.log(T_float)
        log_k_val = np.interp(log_T, log_t, log_k)
        return exp(mpf(str(log_k_val)))


def system(r, u, F):
    kr = k_interp(T(r))
    dudr = -mpf('3') * kr / c * F
    if r == 0:
        dFdr = -c * kr * (u - up(r)) / mpf('2')
    else:
        dFdr = -F / r - c * kr * (u - up(r))
    return dudr, dFdr


def runge_kutta_4(u0, F0, h=mpf('1e-3')):
    r_vals = [mpf('0')]
    u_vals = [u0]
    F_vals = [F0]

    while r_vals[-1] < R:
        r = r_vals[-1]
        u = u_vals[-1]
        F = F_vals[-1]

        k1, q1 = system(r, u, F)
        k2, q2 = system(r + h / mpf('2'), u + h * k1 / mpf('2'), F + h * q1 / mpf('2'))
        k3, q3 = system(r + h / mpf('2'), u + h * k2 / mpf('2'), F + h * q2 / mpf('2'))
        k4, q4 = system(r + h, u + h * k3, F + h * q3)

        u_next = u + (h / mpf('6')) * (k1 + mpf('2') * k2 + mpf('2') * k3 + k4)
        F_next = F + (h / mpf('6')) * (q1 + mpf('2') * q2 + mpf('2') * q3 + q4)

        r_vals.append(r + h)
        u_vals.append(u_next)
        F_vals.append(F_next)

    return r_vals, u_vals, F_vals


def shooting_method(chi_min=mpf('0.01'), chi_max=mpf('1.0')):
    while abs(chi_max - chi_min) > EPS:
        chi_mid = (chi_min + chi_max) / mpf('2')
        u0 = chi_mid * up(mpf('0'))
        r_vals, u_vals, F_vals = runge_kutta_4(u0, mpf('0'))

        residual = -F_vals[-1] + mpf('0.39') * c * u_vals[-1]
        #print(f"Невязка граничного условия: {residual}")

        if residual >= 0:
            chi_max = chi_mid
        else:
            chi_min = chi_mid

    final_residual = -F_vals[-1] + mpf('0.39') * c * u_vals[-1]
    print(f"Финальная невязка: {final_residual}")
    return chi_mid


def solute():
    chi = shooting_method()
    u0 = chi * up(mpf('0'))
    r_vals, u_vals, F_vals = runge_kutta_4(u0, mpf('0'))
    up_vals = [up(r) for r in r_vals]

    print(f"Найденный параметр χ с точностью 50 знаков: {chi}")

    # Конвертируем обратно в float для построения графиков
    r_vals_float = [float(r) for r in r_vals]
    u_vals_float = [float(u) for u in u_vals]
    F_vals_float = [float(F) for F in F_vals]
    up_vals_float = [float(u) for u in up_vals]

    return r_vals_float, u_vals_float, F_vals_float, up_vals_float, [float(T(r)) for r in r_vals]


# На случай важных переговоров


# rx, _, Fxt, _, _ = lab2_solute()

# Fx = [x / (10**10 * 2.5) for x in Fxt]


# h = R / N

# F1 = -(1 + 1.17 * h * kinterp[-1]) * u[-1] + u[-2]
# F2 = -(3 + 2.34 * h * kinterp[-1]) * u[-1] + 4 * u[-2] - u[-3]
# F3 = Fint[-1]
# F4 = 2 * 10**4 * 0.39 * u[-2]

# print(h)
# print(F1)
# print(F2)
# print(F3)
# print(F4)

# print(u)

# print(Fx)