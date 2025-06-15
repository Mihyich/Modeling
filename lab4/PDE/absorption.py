from pathlib import Path
import csv

def initAbsorptionCoefTable(filePath: Path) -> tuple[list[float], list[float], list[float], list[float]]:
    """Создание таблицы коэффициентов T, sigma, λ, c"""
    T_t = []
    sigma_t = []
    lambda_t = []
    c_t = []

    try:
        with open(filePath, mode='r', newline='', encoding='utf-8') as file:
            table = [(r[0], r[1], r[2], r[3]) for r in csv.reader(file)]
            T_t = [float(t[0]) for t in table]
            sigma_t = [float(t[1]) for t in table]
            lambda_t = [float(t[2]) for t in table]
            c_t = [float(t[3]) for t in table]
    except FileNotFoundError:
        print(f"Ошибка: Файл '{filePath}' не найден.")
        exit(1)
    except IOError:
        print(f"Ошибка: Проблема с чтением файла '{filePath}'.")
        exit(1)
    except Exception as e:
        print(f"Произошла ошибка: {e}")
        exit(1)
    finally:
        return (T_t, sigma_t, lambda_t, c_t)