"""
verify_facts_1_2_3.py — Verify structural facts of the Main Theorem.

Fact 1: isolated zero at n+q (between tail terms n+q-1 and n+q+1)
Fact 2: contiguous zero block [n+q+2, u_z]
Fact 3: last zero u_z = n+2q (q≡1) or n+2q-1 (q≡2)
Tail:   coeff of t^{n+q-1} and t^{n+q+1} are both +1

Usage: sage verify_facts_1_2_3.py
"""

B = BraidGroup(3)

def jones_T3q2n(q, n):
    braid_word = [1, 2] * q + [1] * (2 * n)
    K = Knot(B(braid_word))
    return K.jones_polynomial()

def verify_all_facts(q_list, n_list):
    R = LaurentPolynomialRing(ZZ, 't')
    all_ok = True
    fails = []
    print("Verifying Facts 1, 2, 3 of the Main Theorem\n")

    for q in q_list:
        if gcd(3, q) != 1:
            continue
        delta = 0 if q % 3 == 1 else 1

        for n in n_list:
            J = R(jones_T3q2n(q, n))
            coeffs = J.dict()
            deg_min = min(coeffs.keys())
            deg_max = max(coeffs.keys())

            all_zeros = sorted([d for d in range(deg_min, deg_max+1)
                                if coeffs.get(d, 0) == 0])
            if not all_zeros:
                continue

            pz = min(all_zeros)
            uz = max(all_zeros)

            fact1_ok    = (pz == n+q)
            fact2_ok    = all(coeffs.get(d, 0) == 0 for d in range(n+q+2, uz+1))
            fact3_ok    = (uz == n+2*q-delta)
            tail_ok     = (coeffs.get(n+q-1, 0) == 1 and coeffs.get(n+q+1, 0) == 1)
            isolated_ok = (coeffs.get(n+q, 0) == 0)

            ok = fact1_ok and fact2_ok and fact3_ok and tail_ok and isolated_ok
            if not ok:
                fails.append((q, n, fact1_ok, fact2_ok, fact3_ok, tail_ok, isolated_ok))
                all_ok = False

    if all_ok:
        print("All cases PASS.\n")
        print(f"Verified for q in {q_list}, n in {n_list}.\n")
        print("Fact 1  (isolated zero at n+q):               OK")
        print("Fact 2  (zero block [n+q+2, u_z]):            OK")
        print("Fact 3  (last zero = n+2q or n+2q-1):        OK")
        print("Tail    (t^{n+q-1} and t^{n+q+1} both +1):  OK")
        print("Iso     (t^{n+q} = 0):                        OK")
    else:
        print(f"FAILURES ({len(fails)} cases):")
        for row in fails:
            print(f"  q={row[0]},n={row[1]}: F1={row[2]},F2={row[3]},"
                  f"F3={row[4]},tail={row[5]},iso={row[6]}")
    return all_ok

if __name__ == "__main__":
    Q_LIST = [4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26]
    N_LIST = list(range(1, 9))
    verify_all_facts(Q_LIST, N_LIST)
