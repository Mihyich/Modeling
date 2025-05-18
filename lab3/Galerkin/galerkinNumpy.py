import matplotlib.pyplot as plt
import numpy as np
from Galerkin.SLEcreator import createSLE, expKoefFactor, expRes, expBasis

def solute(koefCount: int, h: float) -> tuple[list[float], list[float], any]:
    A, b = createSLE(koefCount, expKoefFactor, expRes)
    Anp, bnp = np.array(A), np.array(b)

    Aaug = np.column_stack((Anp, bnp))
    rankA = np.linalg.matrix_rank(A)
    rankAaug = np.linalg.matrix_rank(Aaug)

    if rankA == rankAaug:
        K, rssq = np.linalg.solve(Anp[:koefCount], bnp[:koefCount]), "Нет"
    else:
        K, residuals, _, _ = np.linalg.lstsq(Anp, bnp, rcond=None)
        rssq = sum(r**2 for r in residuals)

    Xaxis = [x * h for x in range(0, int(1 / h) + 1)]
    Yaxis = [expBasis(x, K) for x in Xaxis]

    return (Xaxis, Yaxis, rssq)

def main():
    koefCount = 6
    h = 0.001
    A, b = createSLE(koefCount, expKoefFactor, expRes)
    Anp, bnp = np.array(A), np.array(b)

    Aaug = np.column_stack((Anp, bnp))
    rankA = np.linalg.matrix_rank(A)
    rankAaug = np.linalg.matrix_rank(Aaug)

    print(f'rank(A) = {rankA}\nrank(A|b) = {rankAaug}\n')

    if rankA == rankAaug:
        K = np.linalg.solve(Anp[:koefCount], bnp[:koefCount])

        print("Решение существует")
        print("\nРешение через solve:", K)
        print("\nНевязка: нет")
    else:
        K, residuals, _, _ = np.linalg.lstsq(Anp, bnp, rcond=None)

        print("Точного решения нет")
        print("\nПриближённое решение:", K)
        print("\nНевязка:", residuals)

    Xaxis = [x * h for x in range(0, int(1 / h) + 1)]
    Yaxis = [expBasis(x, K) for x in Xaxis]
    
    plt.plot(Xaxis, Yaxis, "b-", linewidth=2, label='u(x)')
    plt.title(f"Решение u'' - 2xu' + 2u = x при кол-ве коэффициентов: {koefCount}")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.grid(True)
    plt.legend()
    plt.gcf().canvas.manager.set_window_title("График решения дифференциального уравнения (Numpy)")
    plt.show()

# if __name__ == "__main__":
#     main()