"""
verify_formula_BD.py
====================
Verify that the Bavier-Doleshal formula (Proposition in paper)
reproduces J(t) computed directly by SageMath's Knot class.

Usage: sage verify_formula_BD.py
"""

B = BraidGroup(3)

def compute_k_ell(q):
    p = 3
    q_hat = q % p
    for inv in range(1, p):
        if (q_hat * inv) % p == 1:
            q_hat_inv = inv
            break
    k = min(q_hat_inv, p - q_hat_inv)
    for i_val in [0, 1, 2]:
        val = q * k - i_val + 1
        if val % p == 0 and val // p >= 0:
            return k, val // p
    raise ValueError(f"Cannot compute ell for q={q}")


def jones_from_BD(q, n):
    """
    Compute J(t) from the Bavier-Doleshal formula.
    Uses the Lemma 3.12 auxiliary polynomial and converts A^4 -> t.
    """
    R_A = LaurentPolynomialRing(ZZ, 'A')
    A = R_A.gen()
    p = 3
    k, ell = compute_k_ell(q)
    s = 2 * n

    def Xprime(pp, qq):
        return 1 - A**(4*(pp+1)) - A**(4*(qq+1)) + A**(4*(pp+qq))

    exp_fixed = 4*s - 2 + 4*(k + ell - k*ell) + 2*(k*q + ell*p)
    S = sum((-1)**i * A**(-4*i + exp_fixed) * Xprime(p-2*k, q-2*ell)
            for i in range(s))
    numerator_A = ((-1)**s) * A**(2*(p-1)*(q-1)+2*s) * (Xprime(p,q) - (-1)**s * S)

    # Convert A^(4m) -> t^m, checking all exponents are multiples of 4
    coeffs_A = numerator_A.dict()
    R_t = LaurentPolynomialRing(QQ, 't')
    t = R_t.gen()

    num_t = R_t(0)
    for exp_A, c in coeffs_A.items():
        if exp_A % 4 != 0:
            raise ValueError(f"Non-multiple-of-4 exponent A^{exp_A} for q={q}, n={n}")
        num_t += c * t**(exp_A // 4)

    # Divide by (1 - t^2) using polynomial long division
    R_poly = ZZ['x']
    x = R_poly.gen()

    num_coeffs = num_t.dict()
    if not num_coeffs:
        return R_t(0)

    min_t = min(num_coeffs.keys())
    shifted = sum(c * x**(e - min_t) for e, c in num_coeffs.items())
    denom = 1 - x**2
    q_poly, rem = shifted.quo_rem(denom)

    if rem != 0:
        raise ValueError(f"(1-t^2) does not divide N(t) for q={q}, n={n}")

    J_formula = R_t(0)
    for exp_shifted, c in q_poly.dict().items():
        J_formula += c * t**(exp_shifted + min_t)

    return J_formula


def jones_direct(q, n):
    braid_word = [1, 2] * q + [1] * (2 * n)
    K = Knot(B(braid_word))
    return K.jones_polynomial()


def verify_BD_formula(q_list, n_list):
    R_t = LaurentPolynomialRing(QQ, 't')
    all_ok = True
    print("Verifying Bavier-Doleshal formula against direct SageMath computation\n")

    for q in q_list:
        if gcd(3, q) != 1:
            continue
        for n in n_list:
            try:
                J_bd = jones_from_BD(q, n)
                J_direct = R_t(jones_direct(q, n))
                match = (J_bd == J_direct)
                if not match:
                    print(f"  FAIL q={q}, n={n}")
                    print(f"    BD:     {J_bd}")
                    print(f"    Direct: {J_direct}")
                    all_ok = False
            except Exception as e:
                print(f"  ERROR q={q}, n={n}: {e}")
                all_ok = False

    if all_ok:
        print("All cases PASS.")
        print(f"Verified for q in {q_list} and n in {n_list}.")

    return all_ok


if __name__ == "__main__":
    Q_LIST = [4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20]
    N_LIST = list(range(1, 9))
    verify_BD_formula(Q_LIST, N_LIST)
