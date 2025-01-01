import sympy 
from sympy.abc import x
from sympy import div
from sympy import symbols, GF, div, Poly
import random

N = 256
Q = 2**23 - 2**13 + 1 
GAMMA1 = (Q - 1) // 16
GAMMA2 = GAMMA1 // 2
BETA = 375
ALPHA = 2 * GAMMA2

K = 3
L = 2
ETA = 7
SETABITS = 4
BETA = 375
OMEGA = 64

# x = symbols('X')
Fq = GF(Q)

MODULUS = x**N + 1

def random_polynomial_Zq(use_eta: bool=False):
    # Create the polynomial with random coefficients in Z_q
    coefficients = None
    if use_eta: # if coeficients should be small
        coefficients = [random.randint(0, ETA) for _ in range(N)]
    else:
        coefficients = [random.randint(0, Q-1) for _ in range(N)]

    # Generate the polynomial
    polynomial = sum(c * x**i for i, c in enumerate(coefficients))
    
    return Poly(polynomial, x, domain=Fq)

def matrix_mul(A :list[list], B :list[list]):
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
                sum = div(sum + div(A[i][k] * B[k][j], MODULUS)[1], MODULUS)[1]
            C[i][j] = sum
    
    return C

def matrix_sum(A :list[list], B :list[list]):
    n = len(A)
    m = len(A[0])
    if len(B) != n or len(B[0]) != m:
        raise ValueError("Inappropriate matrices")

    C = [[None for j in range(m)] for i in range(n)]
    for i in range(n):
        for j in range(m):
            C[i][j] = div(A[i][j] + B[i][j], MODULUS)[1]
    
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


def Gen():
    A = [[random_polynomial_Zq() for j in range(L)] for i in range(K)]

    s1 = [[random_polynomial_Zq(use_eta=True)] for j in range(L)]
    s2 = [[random_polynomial_Zq(use_eta=True)] for j in range(K)]

    t = matrix_sum(matrix_mul(A, s1), s2)

    return ((A, t), (A, t, s1, s2))

def Sign(sk, M):
    pass

def Verify():
    pass

def example():
    polynomial1 = random_polynomial_Zq()
    polynomial2 = random_polynomial_Zq()

    print(f"polynomial1: {polynomial1}")
    print(f"polynomial2: {polynomial2}")

    
    result_add = div(polynomial1 + polynomial2, MODULUS)[1] # Take remainder after division
    print(f"Addition result modulo {MODULUS}: {result_add}")

    result_mul = div(polynomial1 * polynomial2, MODULUS)[1] # Take remainder after division
    print(f"Multiplication result modulo {MODULUS}: {result_mul}")

def main():
    pk, sk = Gen()
    print(pk)
    print("------------------------------------------")
    print(sk)

if __name__ == "__main__":
    main()
    # example()