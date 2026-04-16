"""
jones_T3q2n.py
==============
Compute the Jones polynomial of the twisted torus knot T(3,q,2,n)
using SageMath's Knot class.

Usage (from within SageMath or sage -python):
    sage jones_T3q2n.py

Requirements: SageMath >= 9.0
"""

B = BraidGroup(3)

def jones_T3q2n(q, n):
    """
    Return the Jones polynomial of T(3,q,2,n).

    Parameters
    ----------
    q : int   (gcd(3,q) must equal 1)
    n : int   (positive)

    Returns
    -------
    Element of LaurentPolynomialRing(ZZ, 't')
    """
    if gcd(3, q) != 1:
        raise ValueError(f"gcd(3,{q}) = {gcd(3,q)} != 1")
    if n < 1:
        raise ValueError("n must be >= 1")

    # Braid word: (sigma_1 sigma_2)^q * sigma_1^(2n)
    braid_word = [1, 2] * q + [1] * (2 * n)
    K = Knot(B(braid_word))
    return K.jones_polynomial()


if __name__ == "__main__":
    print("Jones polynomial of T(3,q,2,n) for small values:\n")
    for q in [4, 5, 7]:
        for n in [1, 2, 3]:
            J = jones_T3q2n(q, n)
            print(f"  T(3,{q},2,{n}): J(t) = {J}")
        print()
