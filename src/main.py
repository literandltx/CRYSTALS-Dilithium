import sympy 
from sympy.abc import x
from sympy import div
from sympy import symbols, GF, div, Poly
import random

# todo set ALPHA and Q
ALPHA = 0
Q = 0

X = symbols('X')

def random_polynomial_Zq(n: int, q: int, eta=None):
    Fq = GF(q)
    
    # Create the polynomial with random coefficients in Z_q
    coefficients = None
    if eta is None:
        coefficients = [random.randint(0, q-1) for _ in range(n)]
    else: # if coeficients should be small
        coefficients = [random.randint(0, eta) for _ in range(n)]

    # Generate the polynomial
    polynomial = sum(c * X**i for i, c in enumerate(coefficients))
    
    return Poly(polynomial, X, domain=Fq)

def matrix_mul(A :list[list], B :list[list], modulus):
    n = len(A)
    m = len(A[0])
    if len(B) != m:
        raise ValueError("Inappropriate matrices")
    p = len(B[0])

    C = [[None for j in range(p)] for i in range(n)]
    for i in range(n):
        for j in range(p):
            sum = 0
            for k in range(m):
                sum = div(sum + div(A[i][k] * B[k][j], modulus)[1], modulus)[1]
            C[i][j] = sum
    
    return C

def matrix_sum(A :list[list], B :list[list], modulus):
    n = len(A)
    m = len(A[0])
    if len(B) != n or len(B[0]) != m:
        raise ValueError("Inappropriate matrices")

    C = [[None for j in range(m)] for i in range(n)]
    for i in range(n):
        for j in range(m):
            C[i][j] = div(A[i][j] + B[i][j], modulus)[1]
    
    return C


def decompose(a):
    # Centralized remainder mod ALPHA
    t = a & 0x7FFFF
    t += (a >> 19) << 9
    t -= ALPHA // 2 + 1
    t += (t >> 31) & ALPHA
    t -= ALPHA // 2 - 1
    a -= t

    # Divide by ALPHA (possible to avoid)
    u = a - 1
    u >>= 31
    a = (a >> 19) + 1
    a -= u & 1

    # Border case
    a0 = Q + t - (a >> 4)
    a &= 0xF
    return a, a0


def Gen(k: int, l: int, n: int, q: int, eta: int):
    modulus = X**n + 1

    A = [[random_polynomial_Zq(n=n, q=q) for j in range(l)] for i in range(k)]

    s1 = [[random_polynomial_Zq(n=n, q=q, eta=eta)] for j in range(l)]
    s2 = [[random_polynomial_Zq(n=n, q=q, eta=eta)] for j in range(k)]

    t = matrix_sum(matrix_mul(A, s1, modulus=modulus), s2, modulus=modulus)

    return ((A, t), (A, t, s1, s2))

def Sign(sk, M):
    pass

def Verify():
    pass

def example():
    n = 256

    # Define the finite field Z_q
    q = 2**23 - 2**13 + 1
    Fq = GF(q)

    modulus = X**n + 1

    polynomial1 = random_polynomial_Zq(n=n, q=q)
    polynomial2 = random_polynomial_Zq(n=n, q=q)

    print(f"polynomial1: {polynomial1}")
    print(f"polynomial2: {polynomial2}")

    
    result_add = div(polynomial1 + polynomial2, modulus)[1] # Take remainder after division
    print(f"Addition result modulo {modulus}: {result_add}")

    result_mul = div(polynomial1 * polynomial2, modulus)[1] # Take remainder after division
    print(f"Multiplication result modulo {modulus}: {result_mul}")

def main():
    n = 256
    q = 2**23 - 2**13 + 1 

    k=4
    l=4
    eta = 5

    pk, sk = Gen(k=k, l=l, n=n, q=q, eta=eta)
    print(pk)
    print("------------------------------------------")
    print(sk)

if __name__ == "__main__":
    main()