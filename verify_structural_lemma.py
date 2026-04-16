"""
verify_structural_lemma.py
==========================
Verify Lemma (Structural lemma) of the paper:
  N(A) has exactly 6 nonzero monomials for n >= 2, and 4 for n = 1,
  at the explicit exponents g1,...,g6 with coefficients +1,-1,-1,+1,-1,+1.

Usage: sage verify_structural_lemma.py
"""

B = BraidGroup(3)

def compute_k_ell(q):
    """Compute k and ell from Lemmas 4.3 and Corollary 3.5 of Bavier-Doleshal."""
    p = 3
    q_hat = q % p
    q_hat_inv = None
    for inv in range(1, p):
        if (q_hat * inv) % p == 1:
            q_hat_inv = inv
            break
    k = min(q_hat_inv, p - q_hat_inv)
    for i_val in [0, 1, 2]:
        val = q * k - i_val + 1
        if val % p == 0 and val // p >= 0:
            return k, val // p
    raise ValueError(f"Could not compute ell for q={q}")


def numerator_explicit(q, n):
    """
    Compute N(A) from the Bavier-Doleshal formula (Lemma 3.12 / Proposition in paper).
    Returns a LaurentPolynomial in A.
    """
    R = LaurentPolynomialRing(ZZ, 'A')
    A = R.gen()
    p = 3
    k, ell = compute_k_ell(q)
    s = 2 * n

    def Xprime(pp, qq):
        return 1 - A**(4*(pp+1)) - A**(4*(qq+1)) + A**(4*(pp+qq))

    exp_fixed = 4*s - 2 + 4*(k + ell - k*ell) + 2*(k*q + ell*p)
    S = sum((-1)**i * A**(-4*i + exp_fixed) * Xprime(p-2*k, q-2*ell)
            for i in range(s))
    numerator = ((-1)**s) * A**(2*(p-1)*(q-1) + 2*s) * (Xprime(p, q) - (-1)**s * S)
    return numerator


def expected_exponents(q, n):
    """
    Return the expected nonzero exponents of N(A) from the Structural Lemma.
    """
    delta = 0 if q % 3 == 1 else 1
    uz = n + 2*q - delta
    uz_A = 4 * (uz + 1)

    g1 = 4*(n + q - 1)
    g2 = 4*(n + q + 3)
    g3 = uz_A
    g4 = uz_A + 4
    if n >= 2:
        g5 = uz_A + 8*n - 4 + 8*delta
        g6 = uz_A + 8*n + 8*delta
        return sorted([g1, g2, g3, g4, g5, g6])
    else:
        return sorted([g1, g2, g3, g4])


def verify_structural_lemma(q_list, n_list):
    all_ok = True
    print("Verifying Structural Lemma: N(A) has at most 6 nonzero monomials")
    print("at explicit exponents with coefficients +1,-1,-1,+1,-1,+1\n")

    for q in q_list:
        if gcd(3, q) != 1:
            continue
        for n in n_list:
            N = numerator_explicit(q, n)
            coeffs = {e: c for e, c in N.dict().items() if c != 0}
            actual_exps = sorted(coeffs.keys())
            expected_exps = expected_exponents(q, n)

            exp_match = (actual_exps == expected_exps)

            # Check coefficients: should be +1,-1,-1,+1,-1,+1 (or +1,-1,-1,+1 for n=1)
            expected_coeffs = [1, -1, -1, 1, -1, 1] if n >= 2 else [1, -1, -1, 1]
            actual_coeffs = [coeffs[e] for e in actual_exps] if exp_match else []
            coeff_match = (actual_coeffs == expected_coeffs) if exp_match else False

            ok = exp_match and coeff_match
            if not ok:
                print(f"  FAIL q={q}, n={n}:")
                print(f"    actual   exps: {actual_exps}")
                print(f"    expected exps: {expected_exps}")
                if exp_match:
                    print(f"    actual   coeffs: {actual_coeffs}")
                    print(f"    expected coeffs: {expected_coeffs}")
                all_ok = False

    if all_ok:
        print("All cases PASS.")
        print(f"Verified for q in {q_list} and n in {n_list}.")
    return all_ok


if __name__ == "__main__":
    Q_LIST = [4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20]
    N_LIST = list(range(1, 9))
    verify_structural_lemma(Q_LIST, N_LIST)
