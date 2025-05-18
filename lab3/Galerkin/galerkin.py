import matplotlib.pyplot as plt
from Galerkin.SLEcreator import createSLE, expKoefFactor, expRes, expBasis

# Транспонирование матрицы
transpose = lambda matrix: [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

# Умножение матриц A (m*n) и B (n*k)
matmul = lambda A, B: [[sum(A[i][l] * B[l][j] for l in range(len(A[0]))) for j in range(len(B[0]))] for i in range(len(A))]

# Вычисление невязки Ax - b
residual = lambda A, x, b: [sum(A[i][j] * x[j] for j in range(len(x))) - b[i] for i in range(len(b))]

# Вычисление сумму квадратов невязок:
residual_ssq = lambda A, x, b: sum(r**2 for r in residual(A, x, b))

def gauss_elimination(A, b):
    """Решение СЛАУ Ax = b методом Гаусса."""
    n = len(A)
    # Прямой ход
    for col in range(n):
        # Выбор главного элемента
        max_row = max(range(col, n), key=lambda i: abs(A[i][col]))
        A[col], A[max_row] = A[max_row], A[col]
        b[col], b[max_row] = b[max_row], b[col]
        # Обнуление элементов ниже диагонали
        for row in range(col + 1, n):
            factor = A[row][col] / A[col][col]
            b[row] -= factor * b[col]
            for c in range(col, n):
                A[row][c] -= factor * A[col][c]
    # Обратный ход
    x = [0] * n
    for row in reversed(range(n)):
        x[row] = b[row] / A[row][row]
        for col in range(row + 1, n):
            x[row] -= A[row][col] * x[col] / A[row][row]
    return x

def least_squares(A, b):
    """МНК-решение системы Ax = b."""
    At = transpose(A)
    AtA = matmul(At, A)
    Atb = matmul(At, [[bi] for bi in b])
    return gauss_elimination(AtA, [row[0] for row in Atb])

def solute(koefCount: int, h: float) -> tuple[list[float], list[float], float]:
    A, b = createSLE(koefCount, expKoefFactor, expRes)
    K = least_squares(A, b)

    Xaxis = [x * h for x in range(0, int(1 / h) + 1)]
    Yaxis = [expBasis(x, K) for x in Xaxis]

    return (Xaxis, Yaxis, residual_ssq(A, K, b))

def main():
    koefCount = 6
    h = 0.001
    A, b = createSLE(koefCount, expKoefFactor, expRes)
    K = least_squares(A, b)
    print("Коэффициенты: ", K)
    print("\nНевязка: ", residual(A, K, b))
    print("\nСумма квадратов невязок: ", residual_ssq(A, K, b))

    Xaxis = [x * h for x in range(0, int(1 / h) + 1)]
    Yaxis = [expBasis(x, K) for x in Xaxis]

    plt.plot(Xaxis, Yaxis, "b-", linewidth=2, label='u(x)')
    plt.title(f"Решение u'' - 2xu' + 2u = x при кол-ве коэффициентов: {koefCount}")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.grid(True)
    plt.legend()
    plt.gcf().canvas.manager.set_window_title("График решения дифференциального уравнения (МНК)")
    plt.show()
    

# if __name__ == "__main__":
#     main()