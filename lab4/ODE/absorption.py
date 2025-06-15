from pathlib import Path
import csv

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
        exit(1)
    except IOError:
        print(f"Ошибка: Проблема с чтением файла '{filePath}'.")
        exit(1)
    except Exception as e:
        print(f"Произошла ошибка: {e}")
        exit(1)
    finally:
        return (T, k)