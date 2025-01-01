import sympy 
from sympy.abc import x
from sympy import div
from sympy import symbols, GF, div, Poly
import random
import hashlib
import math

N = 256
Q = 2**23 - 2**13 + 1 
GAMMA1 = (Q - 1) // 16
GAMMA2 = GAMMA1 // 2
K = 3
L = 2
BETA = 175
ETA = 3

ALPHA = 2 * GAMMA2

Fq = GF(Q)

MODULUS = x**N + 1

SHAKE128_LENGTH = (256 // 4)
# SHAKE128_LENGTH = 4

def SHAKE128(s: str) -> list[int]:
    number = int(hashlib.shake_128(s.encode('utf-8')).hexdigest(SHAKE128_LENGTH), base=16)
    
    c = [0 for i in range(256)]

    for i in range(60):
        c[i] = pow(-1, number & 1)
        number = number >> 1
    
    return Poly(c, x, domain=Fq)

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

def matrix_const_mul(A :list[list], c :int):
    n = len(A)
    m = len(A[0])

    C = [[None for j in range(m)] for i in range(n)]
    for i in range(n):
        for j in range(m):
            C[i][j] = div(c * A[i][j], MODULUS)[1]
    
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

def matrix_sub(A :list[list], B :list[list]):
    n = len(A)
    m = len(A[0])
    if len(B) != n or len(B[0]) != m:
        raise ValueError(f"Inappropriate matrices:: (nxm) = {n}x{m}, (nxm) = {len(B)}x{len(B[0])}")

    C = [[None for j in range(m)] for i in range(n)]
    for i in range(n):
        for j in range(m):
            C[i][j] = div(A[i][j] - B[i][j], MODULUS)[1]
    
    return C

def decompose(a):
    r = a % Q
    r0 = modpm(r, modulus=ALPHA)
    if r - r0 == Q - 1:
        r1 = 0
        r0 = r0 - 1
    else:
        r1 = (r - r0) // ALPHA

    return (r1, r0)

def Gen():
    A = [[random_polynomial_Zq() for j in range(L)] for i in range(K)]

    s1 = [[random_polynomial_Zq(use_eta=True)] for j in range(L)]
    s2 = [[random_polynomial_Zq(use_eta=True)] for j in range(K)]

    t = matrix_sum(matrix_mul(A, s1), s2)

    return ((A, t), (A, t, s1, s2))

def HighBits(a):
    r1, r0 = decompose(a)
    return r1

def LowBits(a):
    r1, r0 = decompose(a)
    return r0

def modpm(a, modulus: int=Q):
    res = a % modulus
    if res > (modulus // 2):
        res = abs(res - modulus)
    
    return res

def Norm_polyvec(polyveck):
    res = -math.inf;

    for poly in polyveck:
        for coeff in poly[0].all_coeffs():
            res = max(res, modpm(coeff))

    return res;

def Norm_polyvec_check(polyveck, b):
    for poly in polyveck:
        for coeff in poly[0].all_coeffs():
            if modpm(coeff) >= b:
                return True

    return False;

def Norm_matr(matr):
    res = -math.inf;

    for row in matr:
        for cell in row:
            res = max(res, modpm(cell))

    return res;

def Norm_matr_check(matr, b):
    for row in matr:
        for cell in row:
            if modpm(cell) >= b:
                return True

    return False

def Sign(sk, M):
    A, t, s1, s2 = sk

    z = None

    while z is None:
        y = [[random_polynomial_Zq(use_eta=True)] for i in range(L)]
        w1 = [[HighBits(coef) for coef in poly[0].all_coeffs()] for poly in matrix_mul(A, y)]
        c = SHAKE128(M + str(w1))
        z = matrix_sum(y, matrix_const_mul(c=c, A=s1))
        temp_value = matrix_sub(matrix_mul(A, y), matrix_const_mul(s2, c=c))
        temp_value = [[LowBits(coef) for coef in poly[0].all_coeffs()] for poly in temp_value]
        
        ch1 = Norm_polyvec_check(z, b=(GAMMA1 - BETA))
        ch2 = Norm_matr_check(temp_value, b=(GAMMA2 - BETA))
        if ch1 or ch2:          
            z = None

    return (z, c)

def Verify(pk, M, sigma):
    A, t = pk
    z, c = sigma

    temp_value = matrix_sub(matrix_mul(A, z), matrix_const_mul(t, c=c))
    w1 = [[HighBits(coef) for coef in poly[0].all_coeffs()] for poly in temp_value]

    return Norm_polyvec(z) < (GAMMA1 - BETA) and c == SHAKE128(M + str(w1))

def example(): # polinomial multipication modulo MODULUS (which is polinom)
    polynomial1 = random_polynomial_Zq()
    polynomial2 = random_polynomial_Zq()

    print(f"polynomial1: {polynomial1}")
    print(f"polynomial2: {polynomial2}")

    
    result_add = div(polynomial1 + polynomial2, MODULUS)[1] # Take remainder after division
    print(f"Addition result modulo {MODULUS}: {result_add}")

    result_mul = div(polynomial1 * polynomial2, MODULUS)[1] # Take remainder after division
    print(f"Multiplication result modulo {MODULUS}: {result_mul}")

def main():
    print("Key generation...")
    pk, sk = Gen()
    print("Key is generated")
    print("")

    M1 = "aaaaaaaaaaaaa"
    M2 = "aaaaaaaaaaaaaaa"
    print("Sign first message...")
    sigma1 = Sign(sk=sk, M=M1)
    print("Sign second message...")
    sigma2 = Sign(sk=sk, M=M2)
    print("Messages are signed...")
    print("")

    print("Verify signature of the messages")
    print(f"First signature correctness:", Verify(pk=pk, M=M1, sigma=sigma1))
    print(f"Seccond signature correctness:", Verify(pk=pk, M=M2, sigma=sigma2))
    print("")

    print("Sanity check...")
    print("Check Verify(pk, M1, sigma2):", Verify(pk=pk, M=M1, sigma=sigma2))
    print("Check Verify(pk, M2, sigma1):", Verify(pk=pk, M=M2, sigma=sigma1))

if __name__ == "__main__":
    main()

