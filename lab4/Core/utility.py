from math import exp, log, pi

# Температурное поле
TempField = lambda x, T0, Tw, R, p: (Tw - T0) * (x / R)**p + T0

# Функция Планка
Plank = lambda x, T0, Tw, R, p: 0.0003084 / (exp(47990 / TempField(x, T0, Tw, R, p)) - 1)

# Электрический ток в момент времени t в течение импульса
PulseCurrent = lambda t, I_max, t_max: (I_max / t_max) * t * exp(-(t / t_max - 1))


def ElectricFieldStrenght(
    t: float,
    I_max: float,
    t_max: float,
    T0: float,
    Tw: float,
    R: float,
    p: float,
    T_table: list[float],
    sigma_table: list[float],
    integrate: callable
) -> float:
    """Напряженность электрического поля"""

    def integrand(x: float) -> float:
        """Подынтегральная функция: sigma(T(x)) * x"""
        T_x = TempField(x, T0, Tw, R, p)
        sigma_x = logInterpXT(T_table, sigma_table, T_x)
        return sigma_x * x

    integral = integrate(0, R, integrand)
    current = PulseCurrent(t, I_max, t_max)
    
    return current / (2 * pi * integral)


def logInterp(
    k1: float,
    k2: float,
    T1: float,
    T2: float,
    t: float
) -> float:
    """Логарифмическая интерполяция"""
    ln_k1 = log(k1)
    ln_k2 = log(k2)
    ln_k = ln_k1 + (ln_k2 - ln_k1) * (t - T1) / (T2 - T1)
    return exp(ln_k)


def logInterpXT(
    T: list[float], 
    X: list[float], 
    t: float
) -> float:
    """Логарифмическая интерполяция X(T)"""

    # Проверка на выход за границы таблицы
    if t <= T[0]: return X[0]
    if t >= T[-1]: return X[-1]

    res = X[-1]
    
    # Поиск интервала, в который попадает T
    for i in range(len(T) - 1):
        if T[i] <= t <= T[i+1]:
            T1, T2 = T[i], T[i+1]
            k1, k2 = X[i], X[i+1]
            res = logInterp(k1, k2, T1, T2, t)
            break
    
    return res