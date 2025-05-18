from pathlib import Path
import csv
from math import log, exp


def initAbsorptionCoefTable(filePath: Path, variant: int) -> tuple[list[float], list[float]]:
    """Создание таблицы коэффициентов поглощения"""
    v = min(2, max(1, variant))
    T = []
    k = []

    try:
        with open(filePath, mode='r', newline='', encoding='utf-8') as file:
            table = [(r[0], r[v]) for r in csv.reader(file)]
            T = [float(t[0]) for t in table]
            k = [float(t[1]) for t in table]
    except FileNotFoundError:
        print(f"Ошибка: Файл '{filePath}' не найден.")
    except IOError:
        print(f"Ошибка: Проблема с чтением файла '{filePath}'.")
    except Exception as e:
        print(f"Произошла ошибка: {e}")
    finally:
        return (T, k)


def logInterp(k1: float, k2: float, T1: float, T2: float, t: float) -> float:
    """Логарифмическая интерполяция"""
    ln_k1 = log(k1)
    ln_k2 = log(k2)
    ln_k = ln_k1 + (ln_k2 - ln_k1) * (t - T1) / (T2 - T1)
    return exp(ln_k)


def logInterpKT(T: list[float], k: list[float], t: float) -> float:
    """Логарифмическая интерполяция k(T) по таблице из T и k"""

    # Проверка на выход за границы таблицы
    if t <= T[0]: return k[0]
    if t >= T[-1]: return k[-1]

    res = k[-1]
    
    # Поиск интервала, в который попадает T
    for i in range(len(T) - 1):
        if T[i] <= t <= T[i+1]:
            T1, T2 = T[i], T[i+1]
            k1, k2 = k[i], k[i+1]
            res = logInterp(k1, k2, T1, T2, t)
            break
    
    return res


def main():
    filePath = Path('data.csv')
    variant = 1
    T, k = initAbsorptionCoefTable(filePath, variant)
    t = 3000
    print(t, logInterpKT(T, k, t))


if __name__ == "__main__":
    main()