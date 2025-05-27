export function psiPolyBasis(x, i) {
    return Math.pow(x, i);
}

export function psiPolyBasisDeriv(x, i) {
    return i * Math.pow(x, i - 1);
}


export function psiTrigBasis(x, i, a, b) {
    return Math.sin((i * Math.PI * (x - a)) / (b - a));
}

export function psiTrigBasisDeriv(x, i, a, b) {
    const k = i * Math.PI / (b - a);
    return k * Math.cos(k * (x - a));
}