# DiffAlg — Differential Algebra for SageMath

A port of the [Macaulay2 package DiffAlg](https://doi.org/10.2140/jsag.2019.9.11) to SageMath.

DiffAlg computes the standard operations with polynomial differential forms and vector fields. Its main purpose is to associate algebraic objects to differential operators in the exterior algebra of differential forms.

**Original Macaulay2 package by** M. Dubinsky, C. Massri, A. Molinuevo, F. Quallbrunn.  
*The Journal of Software for Algebra and Geometry*, vol. 9, 2019. [DOI: 10.2140/jsag.2019.9.11](https://doi.org/10.2140/jsag.2019.9.11)

---

## Requirements

- [SageMath](https://www.sagemath.org/) ≥ 9.0
- Python ≥ 3.8

## Installation

Clone the repository and run SageMath from the `DiffAlg/` directory:

```bash
git clone https://github.com/amoli79/diffalg.git
cd diffalg
sage -python -c "from diff_alg import *"
```

Or run interactively:

```bash
sage
```
```python
sage: import sys; sys.path.insert(0, '/path/to/diffalg')
sage: from diff_alg import *
```

## Quick Start

```python
from diff_alg import *

# Generic linear 1-form in A^3 with parameters a0..a8
w = new_form(2, 1, 1, "a")
print(w)
# (a0*x0 + a3*x1 + a6*x2)*dx0 + (a1*x0 + a4*x1 + a7*x2)*dx1 + (a2*x0 + a5*x1 + a8*x2)*dx2

# Radial (Euler) vector field
R = radial(2)
print(R)
# x0*∂/∂x0 + x1*∂/∂x1 + x2*∂/∂x2

# Exterior differential
dw = w.exterior_diff()

# Lie derivative
Lw = w.lie_derivative(R)

# Contraction
iRw = w.contraction(R)

# Kernel of contraction with R (logarithmic forms)
log_forms = gen_ker(lambda f: f.contraction(R), 2, 1, 2, "a")
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
| `wedge(other)` | Exterior product $\omega \wedge \eta$ |
| `exterior_diff()` | Exterior differential $d\omega$ |
| `contraction(field)` | Interior product $\iota_X \omega$ |
| `lie_derivative(field)` | Lie derivative $\mathcal{L}_X \omega$ |
| `pullback(map)` | Pullback $\phi^* \omega$ |
| `singular_ideal()` | Ideal of tangency/singular locus |
| `moduli_ideal()` | Moduli ideal |
| `homogenize()` | Homogenize the form |
| `projectivize()` | Project to $\mathbb{P}^n$ |

### Operations on `DiffAlgField`

| Method | Description |
|--------|-------------|
| `bracket(other)` | Lie bracket $[X, Y]$ |
| `is_homogeneous()` | Test homogeneity of coefficients |

### Constructors

| Function | Description |
|----------|-------------|
| `new_form(n, r, d, var)` | Generic $r$-form of degree $d$ in $n+1$ variables |
| `new_form(expr)` | Parse a form from a string expression |
| `new_field(n, d, var)` | Generic vector field of degree $d$ |
| `new_field(expr)` | Parse a vector field from a string |
| `radial(n)` | Radial (Euler) field $\sum x_i \partial/\partial x_i$ |
| `logarithmic_form(n, r, d, var)` | Generic logarithmic $r$-form |
| `linear_comb(forms, var)` | Parametric linear combination |
| `dist(*fields)` | Distribution from a list of fields |

### Kernel and Image

```python
# Kernel of a linear operator on forms
# Returns a list of forms spanning the kernel
ker = gen_ker(lambda f: f.contraction(R), n=2, r=1, d=2, var="a")

# Image of a linear operator on forms
im = gen_im(lambda f: f.exterior_diff(), n=2, r=1, d=1, var="a")
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
sage -python -m pytest tests/
```

## References

- M. Dubinsky, C. Massri, A. Molinuevo, F. Quallbrunn. *DiffAlg: a Macaulay2 package for differential algebra.* J. Softw. Algebra Geom. 9 (2019), 11–16. [DOI: 10.2140/jsag.2019.9.11](https://doi.org/10.2140/jsag.2019.9.11)
