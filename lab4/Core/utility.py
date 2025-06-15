from Core.integrator import integrateSimpson
from math import exp, log, pi
from numpy import interp

from Core.defines import Tw

# Температурное поле (Тут вроде нормально, как и должно быть)
TempField = lambda x, T0, Tw, R, p: (Tw - T0) * ((x / R)**p) + T0
# TempField = lambda x, T0, Tw, R, p: (Tw - T0) * x**p + T0

# Функция Планка (Почему-то именно так???)
# Plank = lambda x, T0, Tw, R, p: 3.084e-4 / (exp(4.799e4 / TempField(x, T0, Tw, R, p)) - 1)
Plank = lambda x, T0, Tw, R, p: 3.084e-4 / (exp(4.799e4 / max(x, 1e-10)) - 1)

# Электрический ток в момент времени t в течение импульса
PulseCurrent = lambda t, t_max, I_max: (I_max / t_max) * t * exp(-(t / t_max - 1))


def ElectricFieldStrenght(
    r: list[float],      # Сетка по радиусу
    T: list[float],      # Текущее распределение температуры
    R: float,
    I: float,
    T_table: list[float],
    sigma_table: list[float],
    integrate: callable = integrateSimpson
) -> float:
    """Напряженность электрического поля"""

    def integrand(x: float) -> float:
        """Подынтегральная функция: sigma(T(x)) * x"""
        t = interp(r, T, x, linearInterp)
        sigma_x = interp(T_table, sigma_table, t, linearInterp)
        return sigma_x * x

    integral = integrate(0, R, integrand)
    
    return I / (2 * pi * integral)


def linearInterp(
    a1: float,
    a2: float,
    b1: float,
    b2: float,
    x: float
) -> float:
    """Линейная интерполяция"""

    return b1 + (x - a1) / (a2 - a1) * (b2 - b1)


def logInterp(
    a1: float,
    a2: float,
    b1: float,
    b2: float,
    x: float
) -> float:
    """Логарифмическая интерполяция"""

    ln_a1, ln_a2 = log(a1), log(a2)
    ln_b1, ln_b2 = log(b1), log(b2)
    ln_x = log(x)

    return exp(linearInterp(ln_a1, ln_a2, ln_b1, ln_b2, ln_x))


def interp(
    A: list[float], 
    B: list[float], 
    x: float,
    interp_func: callable = logInterp
) -> float:
    """Интерполяция B(A), массив A должен быть упорядочен по возрастанию"""

    # Проверка на выход за границы
    if x <= A[0]: return B[0]
    if x >= A[-1]: return B[-1]

    # Проверка длин массивов
    if len(A) != len(B):
        raise ValueError("Массивы имеют разные длины")
    
    # Проверка упорядоченности массива A
    if not all(A[i] < A[i+1] for i in range(len(A)-1)):
        raise ValueError("Массив A должен быть строго возрастающим")
    
    # Поиск интервала, в который попадает x В массиве A
    for i in range(len(A) - 1):
        if A[i] <= x <= A[i+1]:
            a1, a2 = A[i], A[i+1]
            b1, b2 = B[i], B[i+1]
            return interp_func(a1, a2, b1, b2, x)

    # Этого по логике никогда не должно случиться при проверках выше 
    raise ValueError("Не найден промежуток в массиве A, в который попадает x")


def soluteTriSys(
    A: list[float],
    B: list[float],
    C: list[float],
    F: list[float]
) -> list[float]:
    """Прогонка"""

    if len(A) != len(B) or len(B) != len(C) or len(C) != len(F):
        raise ValueError("Длины массивов коэффициентов разностной схемы различны")

    N = len(A)

    # Инициализация
    zita = [0.0] * (N)
    ita = [0.0] * (N)
    
    zita[0] = 1
    ita[0] = 0

    # Прямой ход
    for i in range(1, N-1):
        denominator = B[i] - A[i] * zita[i - 1]

        if abs(denominator) < 1e-12:
            raise ZeroDivisionError(f"Деление на ноль при i={i}")
        
        zita[i] = C[i] / denominator
        ita[i] = (A[i] * ita[i - 1] + F[i]) / denominator
    
    y = [0.0] * N
    y[N - 1] = Tw

    # Обратный ход
    for i in range(N-2, -1, -1):
        y[i] = zita[i] * y[i+1] + ita[i]

    return y