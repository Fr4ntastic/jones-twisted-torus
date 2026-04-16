# Jones Polynomial of Twisted Torus Knots T(3,q,2,n)

Computational verification code for the paper:

> **Fine Structure of the Jones Polynomial of Twisted Torus Knots T(3,q,2,n)**
> Francesco Cancello, Independent Researcher, Naples, Italy

## Paper Summary

We prove that the Jones polynomial J(t) of T(3,q,2,n), where gcd(3,q)=1
and n≥1, decomposes into three explicit blocks:

1. **Tail**: two positive terms at degrees n+q-1 and n+q+1
2. **Zero block**: contiguous zeros from degree n+q+3 to u_z
3. **Head block**: alternating block of length L from u_z+1 to deg_max

where u_z = n+2q (if q≡1 mod 3) or n+2q-1 (if q≡2 mod 3), and
deg_max = 3n+2q-1 (if q≡1 mod 3) or 3n+2q (if q≡2 mod 3).

The key result is that the auxiliary polynomial numerator N(A) has at most
**6 nonzero monomials**, independent of n.

## Requirements

- [SageMath](https://www.sagemath.org/) >= 9.0

No other dependencies are needed. All scripts use SageMath's built-in
`Knot`, `BraidGroup`, and `LaurentPolynomialRing` classes.

## Files

| File | Description |
|------|-------------|
| `jones_T3q2n.py` | Compute J(t) for T(3,q,2,n) via SageMath |
| `verify_structural_lemma.py` | Verify the 6-monomial structure of N(A) |
| `verify_facts_1_2_3.py` | Verify Facts 1, 2, 3 of the Main Theorem |
| `verify_formula_BD.py` | Compare Bavier-Doleshal formula with direct computation |

## Usage

Run any script with SageMath:

```bash
sage jones_T3q2n.py
sage verify_structural_lemma.py
sage verify_facts_1_2_3.py
sage verify_formula_BD.py
```

Or load interactively in a SageMath session:

```python
sage: load('verify_structural_lemma.py')
```

## Expected Output

### verify_structural_lemma.py
```
Verifying Structural Lemma: N(A) has at most 6 nonzero monomials
at explicit exponents with coefficients +1,-1,-1,+1,-1,+1

All cases PASS.
Verified for q in [4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20] and n in [1, ..., 8].
```

### verify_facts_1_2_3.py
```
Verifying Facts 1, 2, 3 of the Main Theorem

All cases PASS.
Fact 1 (first zero = n+q+3):                    OK
Fact 2 (contiguous zero block):                  OK
Fact 3 (last zero = n+2q or n+2q-1 by q mod 3): OK
```

### verify_formula_BD.py
```
Verifying Bavier-Doleshal formula against direct SageMath computation

All cases PASS.
```

## Quick Example

```python
sage: load('jones_T3q2n.py')
sage: J = jones_T3q2n(4, 2)
sage: print(J)
-t^13 + t^12 - t^11 + t^7 + t^5
```

## Reference

B. Bavier and B. Doleshal,
*The Jones polynomial for a torus knot with twists*,
Topology and its Applications, 2024. arXiv:2308.00502.

## License

MIT License. See `LICENSE` for details.
