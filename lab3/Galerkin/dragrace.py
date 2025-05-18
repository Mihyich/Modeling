import matplotlib.pyplot as plt

def dragrace(n: int):
    '''Прогонка'''
    h = 1.0/n
    x = [i / n for i in range(0, n + 1)]
    
    # Инициализация массивов
    a = [0] * (n + 1)
    b = [0] * (n + 1)
    c = [0] * (n + 1)
    d = [0] * (n + 1)
    
    # Граничное условие слева
    b[0] = 1
    c[0] = 0
    d[0] = 0
    
    # Внутренние точки
    for i in range(1, n):
        a[i] = 1/h**2 + x[i]/h
        b[i] = -2/h**2 + 2
        c[i] = 1/h**2 - x[i]/h
        d[i] = x[i]
    
    # Граничное условие справа (u'(1)=1)
    a[n] = -1/h
    b[n] = 1/h
    d[n] = 1
    
    # Метод прогонки
    u = [0] * (n + 1)
    
    # Прямой ход
    for i in range(1, n+1):
        m = a[i]/b[i-1]
        b[i] = b[i] - m*c[i-1]
        d[i] = d[i] - m*d[i-1]
    
    # Обратный ход
    u[n] = d[n]/b[n]
    for i in range(n-1, -1, -1):
        u[i] = (d[i] - c[i]*u[i+1])/b[i]
    
    return x, u

def solute(NexusCount: int) -> tuple[list[float], list[float]]:
    return dragrace(NexusCount)