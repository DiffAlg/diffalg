# DiffAlg — Differential Algebra for SageMath

A port of the [Macaulay2 package DiffAlg](https://doi.org/10.2140/jsag.2019.9.11) to SageMath.

DiffAlg computes the standard operations with polynomial differential forms and vector fields. Its main purpose is to associate algebraic objects to differential operators in the exterior algebra of differential forms.

**Original Macaulay2 package by** M. Dubinsky, C. Massri, A. Molinuevo, F. Quallbrunn — [DiffAlg — Macaulay2 package documentation](https://macaulay2.com/doc/Macaulay2-1.21/share/doc/Macaulay2/DiffAlg/html/index.html)

---

## Requirements

- [SageMath](https://www.sagemath.org/) ≥ 9.0
- Python ≥ 3.8

## Installation

### Standard

Install directly from GitHub into your SageMath environment:

```bash
sage -pip install git+https://github.com/DiffAlg/diffalg.git
```

Or clone and install in editable mode (changes to the source take effect immediately):

```bash
git clone https://github.com/DiffAlg/diffalg.git
cd diffalg
sage -pip install -e .
```

### Arch Linux (and other PEP 668 distributions)

On Arch Linux, SageMath uses the system Python, which is governed by PEP 668 and
blocks `sage -pip` from writing to system site-packages. Use `pip` directly with
`--user` and `--break-system-packages`:

```bash
git clone https://github.com/DiffAlg/diffalg.git
cd diffalg
pip install --user --break-system-packages -e .
```

`--user` installs to `~/.local/lib/python*/site-packages/`, which Sage loads
automatically — no system files are modified.

After installation, `from diff_alg import *` works in any Sage session.

## Quick Start

```python
from diff_alg import *

# Generic linear 1-form in A^3 with parameters a0..a8
w = new_form(2, 1, 1, "a")
print(w)
# (a0*x0 + a3*x1 + a6*x2)*dx0 + (a1*x0 + a4*x1 + a7*x2)*dx1 + (a2*x0 + a5*x1 + a8*x2)*dx2

# Degree: (n, r, d) = (space dimension, form degree, polynomial degree)
print(w.degree)  # (2, 1, 1)

# Radial (Euler) vector field
R = radial(2)
print(R)
# x0*d/dx0 + x1*d/dx1 + x2*d/dx2

# Exterior differential
dw = w.exterior_diff()

# Lie derivative
Lw = w.lie_derivative(R)

# Contraction (interior product)
iRw = w.contract(R)

# Kernel of contraction with R: forms w such that i_R(w) = 0
h = new_form(2, 1, 2, "b")           # generic 1-form of degree 2 in A^3
ker = gen_ker(h.contract(R), h)       # returns basis of the kernel
```

## Features

### Classes

| Class | Description |
|-------|-------------|
| `DiffAlgForm` | Differential $r$-form with polynomial coefficients |
| `DiffAlgField` | Polynomial vector field |
| `DiffAlgDistribution` | Distribution of vector fields |

### Operations on `DiffAlgForm`

| Method | Description |
|--------|-------------|
| `degree` | `(n, r, d)` — space dimension, form degree, polynomial degree |
| `wedge(other)` | Exterior product $\omega \wedge \eta$ |
| `exterior_diff()` | Exterior differential $d\omega$ |
| `contract(field)` | Interior product $\iota_X \omega$ |
| `lie_derivative(field)` | Lie derivative $\mathcal{L}_X \omega$ |
| `pullback(map)` | Pullback $\phi^* \omega$ |
| `singular_ideal()` | Ideal of tangency/singular locus |
| `moduli_ideal()` | Moduli ideal |
| `is_homogeneous()` | Test homogeneity of coefficients |
| `homogenize()` | Homogenize the form |
| `projectivize()` | Project to $\mathbb{P}^n$ |
| `randomize()` | Replace parameters with random integer coefficients |

### Operations on `DiffAlgField`

| Method | Description |
|--------|-------------|
| `degree` | `(n, d)` — space dimension, polynomial degree |
| `bracket(other)` | Lie bracket $[X, Y]$ |
| `is_homogeneous()` | Test homogeneity of coefficients |
| `randomize()` | Replace parameters with random integer coefficients |

### Operations on `DiffAlgDistribution`

| Method | Description |
|--------|-------------|
| `rank()` | Rank of the distribution at the generic point |
| `is_involutive()` | Test involutivity (Frobenius) |

### Constructors

| Function | Description |
|----------|-------------|
| `new_form(n, r, d, var)` | Generic $r$-form of degree $d$ in $n+1$ variables |
| `new_form(expr)` | Parse a form from a string expression |
| `new_field(n, d, var)` | Generic vector field of degree $d$ |
| `new_field(expr)` | Parse a vector field from a string |
| `radial(n)` | Radial (Euler) field $\sum x_i \partial/\partial x_i$ |
| `logarithmic_form(n, degrees, var)` | Logarithmic form for hypersurfaces of degrees `[d0,d1,...]` |
| `linear_comb(forms, var)` | Parametric linear combination |
| `dist(*fields)` | Distribution from a list of fields |

### Kernel and Image

`gen_ker(expr, var)` and `gen_im(expr, var)` work on a linear expression `expr`
in the free parameters of `var`.

```python
R = radial(2)

# Kernel: find all 1-forms h of degree 2 such that i_R(h) = 0
h = new_form(2, 1, 2, "b")
ker = gen_ker(h.contract(R), h)   # list of basis forms

# Non-homogeneous system: find h such that i_R(h) = w0 (fixed form)
w0 = new_form(2, 0, 1, "c").randomize()
basis, particular = gen_ker(h.contract(R) - w0, h)

# Image: span of d(h) as h ranges over all linear 1-forms
h = new_form(2, 1, 1, "b")
im = gen_im(h.exterior_diff(), h)  # list of independent 2-forms
```

## String Parsing

Forms and fields can be defined from string expressions.  
Supports both `x0`/`dx0` and `x_0`/`dx_0` notation, and `^` or `**` for exponentiation.

```python
w = new_form("x_0*dx_1 - x_1*dx_0")
X = new_field("x0^2 * D0 - x1*x2 * D2")
```

## Running the Tests

```bash
# Standard (most systems)
sage -python -m pytest tests/

# Arch Linux (Sage uses system Python)
python -m pytest tests/
```

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

## References

- M\. Dubinsky, C\. Massri, A\. Molinuevo, F\. Quallbrunn. *DiffAlg: a Macaulay2 package for differential algebra.* J. Softw. Algebra Geom. 9 (2019), 11–16. [DOI: 10.2140/jsag.2019.9.11](https://doi.org/10.2140/jsag.2019.9.11)
- [DiffAlg — Macaulay2 package documentation](https://macaulay2.com/doc/Macaulay2-1.21/share/doc/Macaulay2/DiffAlg/html/index.html)
