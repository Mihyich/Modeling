import { clearPlot, createPlot, plotSolution } from './ploter.js';
import { psiPolyBasis, psiPolyBasisDeriv, psiTrigBasis, psiTrigBasisDeriv} from './basisSamples.js';
import { integrateSimpson } from './intergrator.js';

function buildG(a, b, alpha1, beta1, gamma1, alpha2, beta2, gamma2) {
    const D = (alpha1 + beta1 * a) * beta2 - (alpha2 + beta2 * b) * beta1;

    if (D === 0.0)
        throw new Error("Матрицы условий Робина вырождена");

    const DA = gamma1 * beta2 - gamma2 * beta1;
    const DB = (alpha1 + beta1 * a) * gamma2 - (alpha2 + beta2 * b) * gamma1;

    const A = DA / D;
    const B = DB / D;

    const g = function(x) {
        return A * x + B;
    };

    const dg = function(x) {
        return A;
    }

    return {
        g: g,
        dg: dg
    };
}

function buildPhi(i, a, b, alpha1, beta1, alpha2, beta2, psi, dpsi) {
    const D = alpha1 * alpha2 - (alpha2 + beta2 * (b - a)) * (alpha1 + beta1 * (a - b));

    if (D === 0.0)
        throw new Error("Матрицы условий Робина вырождена");

    const psi_a = psi(a, i, a, b);
    const psi_b = psi(b, i, a, b);

    const dpsi_a = dpsi(a, i, a, b);
    const dpsi_b = dpsi(b, i, a, b);

    const DAi = alpha2 * (alpha1 * dpsi_a + beta1 * psi_a) - (alpha2 * dpsi_b + beta2 * psi_b) * (alpha1 + beta1 * (a - b));
    const DBi = alpha1 * (alpha2 * dpsi_b + beta2 * psi_b) - (alpha2 + beta2 * (b - a)) * (alpha1 * dpsi_a + beta1 * psi_a);

    const Ai = DAi / D;
    const Bi = DBi / D;

    const phi = function(x) {
        return psi(x, i, a, b) - Ai * (x - a) - Bi * (x - b);
    };

    const dphi = function(x) {
        return dpsi(x, i, a, b) - Ai - Bi;
    }

    return {
        phi: phi,
        dphi: dphi
    };
}

function buildBasisFunctions(N, a, b, alpha1, beta1, alpha2, beta2, psi, dpsi) {
    const basis = [];

    for (let i = 1; i <= N; i++) {
        const phi = buildPhi(i, a, b, alpha1, beta1, alpha2, beta2, psi, dpsi);
        basis.push(phi);
    }

    return basis;
}

function buildStiffnessMatrix(basis, p, q, a, b, integrateMethod) {
    const N = basis.length;
    const K = Array.from({ length: N }, () => Array(N).fill(0));

    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
            const I1Func = function(x) {
                return basis[j].dphi(x) * basis[i].dphi(x);
            };

            const I2Func = function(x) {
                return p(x) * basis[j].dphi(x) * basis[i].phi(x);
            }

            const I3Func = function(x) {
                return q(x) * basis[j].phi(x) * basis[i].phi(x);
            }

            const I1 = basis[j].dphi(b) * basis[i].phi(b) - basis[j].dphi(a) * basis[i].phi(a) - integrateMethod(I1Func, a, b, 100);
            
            const I2 = integrateMethod(I2Func, a, b, 100);

            const I3 = integrateMethod(I3Func, a, b, 100);


            K[i][j] = I1 + I2 + I3;
        }
    }

    return K;
}

function buildForceVector(basis, G, p, q, f, a, b, integrateMethod) {
    const N = basis.length;
    const F = [];

    for (let i = 0; i < N; i++) {
        const I2Func = function(x) {
            return p(x) * basis[i].phi(x);
        };

        const I3Func = function(x) {
            return q(x) * G.g(x) * basis[i].phi(x);
        };

        const I4Func = function(x) {
            return f(x) * basis[i].phi(x);   
        }

        const I1 = 0;
        const I2 = G.dg(1) * integrateMethod(I2Func, a, b, 100); // Неважен аргумент D.dg, там все равно производная = константа
        const I3 = integrateMethod(I3Func, a, b, 100);
        const I4 = integrateMethod(I4Func, a, b, 100);
        
        F[i] = -I1 - I2 - I3 + I4;
    }

    return F;
}

function parseFunction(exprStr) {
    try {
        const node = math.parse(exprStr);
        return function(x) {
            try {
                const result = node.evaluate({ x });
                return isNaN(result) || !isFinite(result) ? 0 : result;
            } catch (e) {
                return 0;
            }
        };
    } catch (e) {
        alert("Ошибка в выражении: " + e.message);
        return () => 0;
    }
}

function solveGalerkinSystem(K, F) {
    try {
        return numeric.solve(K, F);
    } catch (e) {
        console.error("Ошибка при решении системы", e);
        return null;
    }
}

function buildSolution(G, c, basis) {
    return function uN(x) {
        let sum = 0;
        for (let i = 0; i < basis.length; ++i)
            sum += c[i] * basis[i].phi(x);
        return G.g(x) + sum;
    };
}

function solveEquation() {
    const pFunc = parseFunction(document.getElementById("function-p").value);
    const qFunc = parseFunction(document.getElementById("function-q").value);
    const fFunc = parseFunction(document.getElementById("function-f").value);

    const a = parseFloat(document.getElementById("interval-a").value);
    const b = parseFloat(document.getElementById("interval-b").value);
    const alpha1 = parseFloat(document.getElementById("coef-alpha1").value);
    const beta1 = parseFloat(document.getElementById("coef-beta1").value);
    const gamma1 = parseFloat(document.getElementById("coef-gamma1").value);
    const alpha2 = parseFloat(document.getElementById("coef-alpha2").value);
    const beta2 = parseFloat(document.getElementById("coef-beta2").value);
    const gamma2 = parseFloat(document.getElementById("coef-gamma2").value);
    const N = parseInt(document.getElementById("num-points").value);

    createPlot("x", "y");
    clearPlot();

    // Шаг 1: g(x) = Ax + B
    const G = buildG(a, b, alpha1, beta1, gamma1, alpha2, beta2, gamma2);

    // Шаг 2: φi(x)
    const basis = buildBasisFunctions(N, a, b, alpha1, beta1, alpha2, beta2, psiTrigBasis, psiTrigBasisDeriv);

    // Шаг 3: Создать матрицу жесткости K
    const K = buildStiffnessMatrix(basis, pFunc, qFunc, a, b, integrateSimpson);

    // Шаг 4: Создать вектор правой части F
    const F = buildForceVector(basis, G, pFunc, qFunc, fFunc, a, b, integrateSimpson);

    // Шаг 5: Решить систему - найти коэффициенты
    const c = solveGalerkinSystem(K, F);

    // Шаг 6: Создать приближенное решение
    const uN = buildSolution(G, c, basis);

    plotSolution(G.g, a, b, N, "g(x)", "green");
    plotSolution(G.dg, a, b, N, "g'(x)", "red");
    plotSolution(uN, a, b, N, "uN(x)", "blue");

    console.log(K);
    console.log(F);
}

document.getElementById("solveBtn").addEventListener("click", function () {
    solveEquation();
});