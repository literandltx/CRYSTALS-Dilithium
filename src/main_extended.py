import sympy 
from sympy.abc import x
from sympy import div
from sympy import symbols, GF, div, Poly
import random
import hashlib
import math

N = 256
Q = 2**23 - 2**13 + 1
D = 14
GAMMA1 = (Q - 1) // 16
GAMMA2 = GAMMA1 // 2
K = 3
L = 2
BETA = 175
ETA = 3
OMEGA = 64

ALPHA = 2 * GAMMA2

Fq = GF(Q)

MODULUS = x**N + 1

SHAKE128_LENGTH = 256 // 4
CRH_OUTPUT_LENGTH = 384 // 4
# SHAKE128_LENGTH = 4

def SHAKE128(s, length: int=SHAKE128_LENGTH) -> list[int]:
    number = s

    if isinstance(s, str):
        number = int(hashlib.shake_128(s.encode('utf-8')).hexdigest(length), base=16)
    
    c = [0 for i in range(256)]

    for i in range(60):
        c[i] = pow(-1, number & 1)
        number = number >> 1
    
    return Poly(c, x, domain=Fq)

def random_polynomial_Zq(use_eta: bool=False, eta: int=ETA):
    # Create the polynomial with random coefficients in Z_q
    coefficients = None
    if use_eta: # if coeficients should be small
        coefficients = [random.randint(0, eta) for _ in range(N)]
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

def Power2Round(r):
    r = r % Q
    twopd = pow(2, D)
    r0 = modpm(a=r, modulus=twopd)

    return(r - r0) // twopd, r0

def Power2Round_polyvec(polyvec):
    p2r = [[Power2Round(coeff) for coeff in poly[0].all_coeffs()] for poly in polyvec]

    # there is a hard enumeration becasue of the main paper
    polyvec1 = []
    polyvec1_coefs = [[cell[0] for cell in line] for line in p2r]
    for i, line in enumerate(polyvec1_coefs):
        poly = sum(c * x**i for i, c in enumerate(line))
        poly = Poly(poly, x, domain=Fq)
        polyvec1.append([poly])

    polyvec0 = []
    polyvec0_coefs = [[cell[0] for cell in line] for line in p2r]
    for i, line in enumerate(polyvec0_coefs):
        poly = sum(c * x**i for i, c in enumerate(line))
        poly = Poly(poly, x, domain=Fq)
        polyvec0.append([poly])

    return polyvec1, polyvec0

def MakeHint(z, r):
    r1 = HighBits(r)
    v1 = HighBits(r + z)

    return r1 != v1

def UseHint(h, r):
    m = (Q - 1) // ALPHA
    r1, r0 = decompose(r)

    if h == 1 and r0 > 0:
        return (r1 + 1) % m

    if h == 1 and r0 <= 0:
        return (r1 - 1) % m
    
    return r1

def ExpandA(rho):
    A = [[None for j in range(L)] for i in range(K)]
    
    for i in range(K):
        for j in range(L):
            seed = (rho << 8) + pow(2, 4) * j + i
            random.seed(int(hashlib.shake_128(repr(seed).encode()).hexdigest(SHAKE128_LENGTH), base=16))
            A[i][j] = random_polynomial_Zq()
    
    return A

def CRH(s: str) -> list[int]:
    number = int(hashlib.shake_128(s.encode('utf-8')).hexdigest(CRH_OUTPUT_LENGTH), base=16)

    return number

def Gen():
    rho = random.randint(0, pow(2, 256) - 1)
    K_ = random.randint(0, pow(2, 256) - 1)
    

    s1 = [[random_polynomial_Zq(use_eta=True)] for j in range(L)]
    s2 = [[random_polynomial_Zq(use_eta=True)] for j in range(K)]

    # A = [[random_polynomial_Zq() for j in range(L)] for i in range(K)]
    A = ExpandA(rho=rho)

    t = matrix_sum(matrix_mul(A, s1), s2)

    t1, t0 = Power2Round_polyvec(polyvec=t)

    tr = CRH(str(rho) + str(t1))

    return ((rho, t1), (rho, K_, tr, s1, s2, t0))

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

def ExpandMask(K_, mu, k_):
    y = [None for i in range(L)]
    for i in range(L):
        seed = (K_ << (48 * 8 + 2 * 8)) + (mu << (2 * 8)) + k_ + i
        random.seed(int(hashlib.shake_128(repr(seed).encode()).hexdigest(SHAKE128_LENGTH), base=16))
        
        y[i] = [random_polynomial_Zq(use_eta=True, eta=GAMMA1 - 1)]

    return y


def decompose_polyvec(polyvec):
    p2r = [[decompose(coeff) for coeff in poly[0].all_coeffs()] for poly in polyvec]

    # there is a hard enumeration becasue of the main paper
    # polyvec1 = []
    polyvec1_coefs = [[cell[0] for cell in line] for line in p2r]
    # for i, line in enumerate(polyvec1_coefs):
    #     poly = sum(c * x**i for i, c in enumerate(line))
    #     poly = Poly(poly, x, domain=Fq)
    #     polyvec1.append([poly])

    polyvec0 = []
    polyvec0_coefs = [[cell[0] for cell in line] for line in p2r]
    for i, line in enumerate(polyvec0_coefs):
        poly = sum(c * x**i for i, c in enumerate(line))
        poly = Poly(poly, x, domain=Fq)
        polyvec0.append([poly])

    # return polyvec1, polyvec0
    return polyvec1_coefs, polyvec0

def Sign(sk, M):
    rho, K_, tr, s1, s2, t0 = sk

    A = ExpandA(rho=rho)
    mu = CRH(str(tr) + M)
    k_ = 0
    z = None
    h = None

    while z is None and h is None:
        # y = [[random_polynomial_Zq(use_eta=True)] for i in range(L)]
        y = ExpandMask(K_=K_, mu=mu, k_=k_)
        
        w = matrix_mul(A, y)
        w1 = [[HighBits(coef) for coef in poly[0].all_coeffs()] for poly in w]
        
        c = SHAKE128(str(mu) + str(w1))
        z = matrix_sum(y, matrix_const_mul(c=c, A=s1))
        
        temp_value = matrix_sub(w, matrix_const_mul(s2, c=c))
        # temp_value = [[LowBits(coef) for coef in poly[0].all_coeffs()] for poly in temp_value]
        r1, r0 = decompose_polyvec(temp_value)
        
        ch1 = Norm_polyvec_check(z, b=(GAMMA1 - BETA))
        ch2 = Norm_polyvec_check(r0, b=(GAMMA2 - BETA))
        ch3 = (r1 != w1)
        if ch1 or ch2 or ch3:
            z = None
            h = None
            print("...next iteration... (start)")
        else:
            ct0 = matrix_const_mul(t0, c=c)
            cs0 = matrix_const_mul(s2, c=c)
            h = [[MakeHint(a, b) for a, b in zip(poly_a[0].all_coeffs(), poly_b[0].all_coeffs())] for poly_a, poly_b in zip(matrix_const_mul(t0, c=(-1 * c)), matrix_sub(w, matrix_sum(cs0, ct0)))]

            n_of_ones = 0
            for poly in h:
                for coeff in poly:
                    if coeff == 1:
                        n_of_ones += 1

            if Norm_polyvec_check(ct0, b=GAMMA2) or n_of_ones > OMEGA:
                z = None
                h = None
                print("...next iteration... (final)")
        k_ = k_ + 1

    return (z, h, c)

def Verify(pk, M, sigma):
    rho, t1 = pk
    z, h, c = sigma

    A = ExpandA(rho=rho)
    mu = CRH(str(CRH(str(rho) + str(t1))) + M)
    
    temp_value = matrix_sub(matrix_mul(A, z), matrix_const_mul(t1, c=c  * pow(2, D)))
    w1 = [[UseHint(a, b) for a, b in zip(poly_a, poly_b[0].all_coeffs())] for poly_a, poly_b in zip(h, temp_value)]

    n_of_ones = 0
    for poly in h:
        for coeff in poly:
            if coeff == 1:
                n_of_ones += 1

    print(Norm_polyvec(z) < (GAMMA1 - BETA), c == SHAKE128(str(mu) + str(w1)), n_of_ones <= OMEGA)
    return Norm_polyvec(z) < (GAMMA1 - BETA) and c == SHAKE128(str(mu) + str(w1)) and n_of_ones <= OMEGA

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
    # print("Sign second message...")
    # sigma2 = Sign(sk=sk, M=M2)
    print("Messages are signed...")
    print("")

    print("Verify signature of the messages")
    print(f"First signature correctness:", Verify(pk=pk, M=M1, sigma=sigma1))
    # print(f"Seccond signature correctness:", Verify(pk=pk, M=M2, sigma=sigma2))
    print("")

    # print("Sanity check...")
    # print("Check Verify(pk, M1, sigma2):", Verify(pk=pk, M=M1, sigma=sigma2))
    # print("Check Verify(pk, M2, sigma1):", Verify(pk=pk, M=M2, sigma=sigma1))

if __name__ == "__main__":
    main()
