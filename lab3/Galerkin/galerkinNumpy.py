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