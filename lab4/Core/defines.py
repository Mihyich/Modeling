from pathlib import Path

# Путь к таблицам
filePathODE = Path.cwd() / Path("lab4/dataODE.csv")
filePathPDE = Path.cwd() / Path("lab4/dataPDE.csv")

# Папка сохранения графиков
folderPathPlot = Path.cwd() / Path("lab4/plots")

# Вариант
variant = 1

#
N = 100

T0 = 8000

Tw = 1800

p = 2

# Радиус цилиндра
R = 0.35

# Скорость света
c = 3*10**10

# Амплитуда импульса тока
I_max = 300

# Время достижения амплитуды импульса тока
t_max = 80e-6

# Шаг времени
tau = 1e-7

# Количество слоев
S = int(t_max / tau)