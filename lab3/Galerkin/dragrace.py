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

def main():
    NexusCount = 100
    x, u = dragrace(NexusCount)

    plt.plot(x, u, "b-", linewidth=2, label='u(x)')
    plt.title(f"Решение u'' - 2xu' + 2u = x при кол-ве точек дискретизации = {NexusCount}")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.grid(True)
    plt.legend()
    plt.gcf().canvas.manager.set_window_title("График решения дифференциального уравнения (Прогонка)")
    plt.show()
    

# if __name__ == "__main__":
#     main()