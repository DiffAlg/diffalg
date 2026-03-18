# DiffAlg — Differential Algebra for SageMath

**Port of the Macaulay2 package DiffAlg to SageMath.**

DiffAlg is a differential algebra package. It computes the usual operations with polynomial differential forms and vector fields. Its main purpose is to associate algebraic objects to differential operators in the exterior algebra of differential forms.

Based on the Macaulay2 package by M. Dubinsky, C. Massri, A. Molinuevo, F. Quallbrunn.
Published in *The Journal of Software for Algebra and Geometry*, vol. 9, 2019. [DOI: 10.2140/jsag.2019.9.11](https://doi.org/10.2140/jsag.2019.9.11)

---

## Getting Started

```python
from diff_alg import *
```

Define a linear differential 1-form $\omega$ and the radial vector field $R$ in 3-dimensional space:

```python
>>> w = new_form(2, 1, 1, "a")
>>> w
(a0*x0 + a3*x1 + a6*x2)*dx0 + (a1*x0 + a4*x1 + a7*x2)*dx1 + (a2*x0 + a5*x1 + a8*x2)*dx2

>>> R = radial(2)
>>> R
x0*∂/∂x0 + x1*∂/∂x1 + x2*∂/∂x2
```

The form $\omega$ is a generic homogeneous 1-form in $\mathbb{A}^3$ with linear polynomial coefficients. The scalar parameters `a0, ..., a8` are free. The radial field is the Euler vector field $R = \sum x_i \partial/\partial x_i$.

---

## Types

### `DiffAlgForm`

The class of all differential forms. A `DiffAlgForm` of type $(n, r, d)$ is an $r$-form in $(n+1)$ affine variables with homogeneous polynomial coefficients of degree $d$.

### `DiffAlgField`

The class of all polynomial vector fields. A `DiffAlgField` of type $(n, d)$ is a vector field in $(n+1)$ variables with homogeneous polynomial coefficients of degree $d$.

### `DiffAlgDistribution`

A distribution of vector fields, used to test involutivity and compute rank.

---

## Constructors

### `new_form(n, r, d, var_name)` — constructor of a differential form

Creates a generic homogeneous $r$-form in $(n+1)$-dimensional affine space with polynomial coefficients of degree $d$.

| Parameter | Type | Description |
|-----------|------|-------------|
| `n` | int | Number of variables minus one |
| `r` | int | Degree of the differential form (exterior degree) |
| `d` | int | Degree of the polynomial coefficients |
| `var_name` | str | Name prefix for the scalar coefficients |

```python
>>> w = new_form(2, 1, 1, "a")
>>> w
(a0*x0 + a3*x1 + a6*x2)*dx0 + (a1*x0 + a4*x1 + a7*x2)*dx1 + (a2*x0 + a5*x1 + a8*x2)*dx2

>>> w.degree
(2, 1, 1)
```

The degree is a triple $(n, r, d)$ — the `n` in front denotes the ambient dimension parameter.

A 2-form with quadratic coefficients in $\mathbb{A}^4$:

```python
>>> eta = new_form(3, 2, 2, "b")
>>> eta.degree
(3, 2, 2)
```

---

### `new_form(expression)` — parse a differential form from a string

Defines a particular differential form from a string expression, analogous to `newForm(String)` in the Macaulay2 package.

Coordinates are written as `x0, x1, ...` (or `x_0, x_1, ...`) and differentials as `dx0, dx1, ...` (or `dx_0, dx_1, ...`). The exterior product is written as ordinary multiplication `*`. Exponentiation can use `^` or `**`. Any unrecognized symbols become scalar parameters in the coefficient ring.

```python
>>> w = new_form("x_0*dx_1 - x_1*dx_0")
>>> w
-x1*dx0 + x0*dx1

>>> w = new_form("a*x_1*dx_0*dx_1")
>>> w
a*x1*dx0*dx1
>>> w.param_names
{'a'}

# The dimension n is inferred from the highest variable index
>>> v = new_form("dx_5")
>>> v.n
5

# 0-forms (polynomial functions) are also supported
>>> f = new_form("x0^2 + x1^2 + x2^2")
>>> f.exterior_diff()
2*x0*dx0 + 2*x1*dx1 + 2*x2*dx2
```

---

### `new_field(n, d, var_name)` — constructor of a vector field

Creates a generic homogeneous vector field in $(n+1)$-dimensional affine space with polynomial coefficients of degree $d$.

| Parameter | Type | Description |
|-----------|------|-------------|
| `n` | int | Number of variables minus one |
| `d` | int | Degree of the polynomial coefficients |
| `var_name` | str | Name prefix for the scalar coefficients |

```python
>>> X = new_field(2, 2, "a")
>>> X.degree
(2, 2)
```

---

### `new_field(expression)` — parse a vector field from a string

Defines a particular vector field from a string expression, analogous to `newField(String)` in the Macaulay2 package.

Coordinates are written as `x0, x1, ...` (or `x_0, x_1, ...`) and partial derivatives $\partial/\partial x_i$ as `ax0, ax1, ...` (or `ax_0, ax_1, ...`). Any unrecognized symbols become scalar parameters.

```python
>>> X = new_field("x_0^2*ax_1 + x_1*ax_0")
>>> X
x1*∂/∂x0 + x0^2*∂/∂x1

>>> X = new_field("a*x0*ax0 + b*x1*ax1")
>>> X.param_names
{'a', 'b'}
```

---

### `radial(n)` — the radial (Euler) vector field

Defines the radial vector field $R = \sum_{i=0}^{n} x_i \frac{\partial}{\partial x_i}$.

```python
>>> R = radial(2)
>>> R
x0*∂/∂x0 + x1*∂/∂x1 + x2*∂/∂x2
```

---

### `logarithmic_form(n, degrees, var_name, projective=False)` — logarithmic forms

A logarithmic form of type $(d_0, \ldots, d_k)$ is a differential 1-form $\omega$ that can be written as:

$$\omega = \left(\prod f_i\right) \sum \frac{df_i}{f_i}$$

where $f_i$ is a homogeneous polynomial of degree $d_i$.

| Parameter | Type | Description |
|-----------|------|-------------|
| `n` | int | Number of variables minus one |
| `degrees` | list | List of degrees $[d_0, \ldots, d_k]$ |
| `var_name` | str | Name prefix for the scalar coefficients |
| `projective` | bool | If `True`, create a form that descends to $\mathbb{P}^n$ |

When the list has length two, the form is called *rational*.

```python
# Generic logarithmic form in A^3 of type (1,1,2)
>>> w = logarithmic_form(2, [1, 1, 2], "a")
>>> w.degree
(2, 1, 4)

# Random projective rational form in P^2
>>> l = logarithmic_form(2, [1, 1], "a", projective=True).randomize()
>>> R = radial(2)
>>> l.contract(R).is_zero   # projective: i_R(l) = 0
True
```

---

### `linear_comb(elements, var_name)` — generic linear combination

Produces a generic linear combination $\sum a_i \cdot e_i$ where $e_i$ are the given elements and $a_i$ are new free scalar parameters.

```python
>>> w1 = new_form(2, 1, 1, "a").randomize()
>>> w2 = new_form(2, 1, 1, "b").randomize()
>>> h = linear_comb([w1, w2], "c")
# h = c0*w1 + c1*w2, with c0, c1 free
```

---

### `dist(fields)` — distribution from vector fields

Creates a `DiffAlgDistribution` from a list of vector fields.

```python
>>> X = new_field(2, 0, "a").randomize()
>>> Y = new_field(2, 0, "b").randomize()
>>> D = dist([X, Y])
>>> D.rank()
2
>>> D.is_involutive()
True
```

Constant-coefficient fields always span an involutive distribution (the bracket of constant fields is zero).

---

## Operations on Forms

### `w.wedge(h)` or `w ^ h` — exterior product

```python
>>> w = new_form(2, 1, 1, "a").randomize()
>>> h = new_form(2, 1, 1, "b").randomize()
>>> omega = w ^ h
>>> omega.degree
(2, 2, 2)

# Anticommutativity: w^h = -h^w for 1-forms
>>> (w ^ h) == -(h ^ w)
True

# w^w = 0 for any 1-form
>>> (w ^ w).is_zero
True
```

---

### `w.exterior_diff()` — exterior derivative $d\omega$

Computes the exterior derivative. Satisfies $d \circ d = 0$.

```python
>>> w = new_form(2, 1, 1, "a")
>>> dw = w.exterior_diff()
>>> dw.degree
(2, 2, 0)

# d(dw) = 0 always
>>> dw.exterior_diff().is_zero
True
```

---

### `w.contract(X)` — contraction $\iota_X(\omega)$

Contraction (interior product) of a differential form with a vector field.

```python
>>> w = new_form(2, 1, 1, "a")
>>> R = radial(2)
>>> w.contract(R)   # i_R(w)
```

Euler's formula: for a homogeneous $r$-form $\omega$ with polynomial coefficients of degree $d$,

$$\iota_R(\omega) \wedge d\omega + d(\iota_R \omega) \wedge \omega = (r + d) \cdot \omega \wedge d\omega$$

```python
>>> w = new_form(2, 1, 1, "a")
>>> R = radial(2)
>>> dw = w.exterior_diff()
>>> Rw = w.contract(R)
>>> lhs = Rw.wedge(dw) + Rw.exterior_diff().wedge(w)
>>> rhs = 2 * w.wedge(dw)   # r + d = 1 + 1 = 2
>>> lhs == rhs
True
```

---

### `w.lie_derivative(X)` — Lie derivative $\mathcal{L}_X(\omega)$

Cartan's formula: $\mathcal{L}_X(\omega) = \iota_X(d\omega) + d(\iota_X\omega)$.

```python
>>> w = new_form(2, 1, 2, "a")
>>> X = new_field(2, 1, "b")
>>> Lx = w.lie_derivative(X)
```

---

### `w.pullback(morphism)` — pull-back by a polynomial map

Pull-back of a differential form by a polynomial map $\phi: \mathbb{A}^m \to \mathbb{A}^n$.

```python
>>> w = new_form(2, 1, 1, "a")
>>> # morphism is a list of polynomials [phi_0, ..., phi_n]
```

---

### `w.singular_ideal()` — singular ideal

The ideal generated by the polynomial coefficients of the form. Its vanishing locus is the singular set of the associated foliation.

```python
>>> w = new_form(2, 1, 1, "a").randomize()
>>> I = w.singular_ideal()
```

---

### `w.moduli_ideal()` — moduli ideal

The ideal generated by the scalar coefficients of the form, viewed as polynomials in the parameters. Used internally by `gen_ker`.

```python
>>> w = new_form(2, 1, 1, "a")
>>> J = w.moduli_ideal()
```

---

### `w.randomize(ring=ZZ, density=1.0, height=10)` — substitute random values

Replaces all free parameters with random values from the specified ring.

```python
>>> w = new_form(2, 1, 1, "a")
>>> w_rand = w.randomize()
>>> len(w_rand.param_names)   # no more free parameters
0
```

---

### `w.homogenize()` / `w.projectivize()` — extend to projective space

`homogenize()` adds a new variable $x_{n+1}$ to make all coefficients have the same total degree.

`projectivize()` applies Jouanolou's projectivization: extends to $\mathbb{P}^{n+1}$ and adds the contraction with the new radial term.

```python
>>> w = new_form(2, 1, 1, "a").randomize()
>>> wp = w.projectivize()
>>> wp.degree
(3, 1, 1)
```

---

## Operations on Fields

### `X.bracket(Y)` or `X | Y` — Lie bracket $[X, Y]$

```python
>>> X = new_field(2, 1, "a").randomize()
>>> Y = new_field(2, 1, "b").randomize()

# Antisymmetry
>>> (X | Y) == -(Y | X)
True

# [X, X] = 0
>>> (X | X).is_zero
True
```

The Jacobi identity holds: $[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$.

---

## Kernel and Image of Linear Operators

These are the main functions for computing algebraic objects associated to differential operators.

### `gen_ker(expr, var)` — basis of the kernel of a linear expression

Given `expr`, an expression that is linear in the parameters of `var`, computes a basis of the kernel — that is, all values of `var` such that `expr = 0`.

For non-homogeneous expressions, returns `[basis, particular_solution]`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr` | DiffAlgForm or DiffAlgField | Expression linear in `var`'s parameters |
| `var` | DiffAlgForm or DiffAlgField | Variable with free scalar coefficients |

**Example 1: Projective 1-forms.** A projective form is a 1-form $h$ satisfying $\iota_R(h) = 0$. We compute a basis of projective 1-forms in $\mathbb{P}^4$ with linear coefficients:

```python
>>> h = new_form(4, 1, 1, "a")
>>> R = radial(4)
>>> T = gen_ker(h.contract(R), h)
>>> len(T)
10

# Verify that all basis elements are projective
>>> all(t.contract(R).is_zero for t in T)
True
```

The 10 basis elements are the forms $x_i dx_j - x_j dx_i$ for $i < j$.

**Example 2: Tangent directions.** Given a foliation $\omega$, the tangent directions are the 1-forms $h$ satisfying $\omega \wedge dh + h \wedge d\omega = 0$. We compute them for a random rational foliation in $\mathbb{P}^4$:

```python
>>> w = logarithmic_form(4, [1, 1], "a", projective=True).randomize()
>>> H = linear_comb(gen_ker(w.contract(radial(4)), new_form(4, 1, 1, "b")), "b")
>>> dw = w.exterior_diff()
>>> dH = H.exterior_diff()
>>> L = gen_ker(w.wedge(dH) + dw.wedge(H), H)
>>> len(L)   # number of tangent directions
4
```

**Example 3: Particular solution (non-homogeneous).** Solve $\omega_1 \wedge h = \omega_3$ for $h$:

```python
>>> w1 = new_form(4, 1, 1, "a").randomize()
>>> w2 = new_form(4, 1, 1, "a").randomize()
>>> w3 = w1 ^ w2
>>> h = new_form(4, 1, 1, "a")
>>> result = gen_ker(w1.wedge(h) - w3, h)
# result = [kernel_basis, particular_solution]
>>> particular = result[1]
>>> (w1 ^ particular) == w3
True
```

**Example 4: Field annihilator.** Find vector fields $X$ such that $\iota_X(\omega) = 0$:

```python
>>> w = new_form(2, 2, 2, "a").randomize()
>>> X = new_field(2, 1, "b")
>>> L = gen_ker(w.contract(X), X)
>>> all(w.contract(fld).is_zero for fld in L)
True
```

---

### `gen_im(expr, var)` — basis of the image of a linear expression

Computes a basis of the image: $\{\ \mathrm{expr}(v) : v \text{ ranges over all values of var}\ \}$.

```python
>>> w = new_form(3, 2, 2, "a").randomize()
>>> h = new_form(3, 1, 1, "b")
>>> dw = w.exterior_diff()
>>> dh = h.exterior_diff()
>>> expr = w.wedge(dh) + h.wedge(dw)
>>> im = gen_im(expr, h)
>>> ker = gen_ker(expr, h)

# Rank-nullity: dim(ker) + dim(im) = number of parameters
>>> len(ker) + len(im) == len(h.param_names)
True
```

---

## Properties

### `w.degree` / `X.degree`

Returns the degree tuple: $(n, r, d)$ for forms, $(n, d)$ for fields.

### `w.is_zero` / `X.is_zero`

Tests if the element is identically zero.

### `w.is_homogeneous()` / `X.is_homogeneous()`

Tests if all polynomial coefficients have the same degree.

### `w.param_names` / `X.param_names`

Set of the names of free scalar parameters.

---

## Distributions

```python
>>> X = new_field(3, 1, "a").randomize()
>>> Y = new_field(3, 1, "b").randomize()
>>> Z = X | Y

>>> D = dist([X, Y, Z])
>>> D.is_involutive()
True

>>> D.rank()
2   # Z is in the span of X and Y (generically)
```

---

## Ring Architecture

- **Polynomial ring**: $\mathbb{Q}[a_0, \ldots, a_k, x_0, \ldots, x_n]$ — parameters first, then coordinates
- **Exterior algebra**: $\bigwedge(R)\langle dx_0, \ldots, dx_n \rangle$ over $R$
- **Vector fields**: stored as dict $\{i: P_i\}$ where $X = \sum P_i \frac{\partial}{\partial x_i}$

When two elements with different rings interact (e.g., wedge, contraction), their rings are automatically merged to a common ring containing all variables and parameters.

---

## Operator Overloading

| Syntax | Operation |
|--------|-----------|
| `w ^ h` | Exterior product (`w.wedge(h)`) |
| `w + h` | Addition |
| `w - h` | Subtraction |
| `-w` | Negation |
| `3 * w` | Scalar multiplication |
| `w / 2` | Scalar division |
| `X \| Y` | Lie bracket (`X.bracket(Y)`) |

---

## Comparison with Macaulay2 DiffAlg

| Macaulay2 | SageMath (this port) |
|-----------|---------------------|
| `newForm(n,r,d,"a")` | `new_form(n, r, d, "a")` |
| `newField(n,d,"a")` | `new_field(n, d, "a")` |
| `radial n` | `radial(n)` |
| `diff w` | `w.exterior_diff()` |
| `w * h` or `w ^ h` | `w.wedge(h)` or `w ^ h` |
| `X _ w` or `w _ X` | `w.contract(X)` |
| `X \| Y` | `X.bracket(Y)` or `X \| Y` |
| `singularIdeal w` | `w.singular_ideal()` |
| `moduliIdeal w` | `w.moduli_ideal()` |
| `genKer(expr, var)` | `gen_ker(expr, var)` |
| `genIm(expr, var)` | `gen_im(expr, var)` |
| `linearComb(L, "a")` | `linear_comb(L, "a")` |
| `logarithmicForm(n,L,"a")` | `logarithmic_form(n, L, "a")` |
| `homogenize w` | `w.homogenize()` |
| `projectivize w` | `w.projectivize()` |
| `random w` | `w.randomize()` |
| `isHomogeneous w` | `w.is_homogeneous()` |
| `isInvolutive D` | `D.is_involutive()` |
| `dist L` | `dist(L)` |
| `rank D` | `D.rank()` |

---

## Caveat

It is recommended to operate in low degrees and dimensions because of the computational time needed to handle the number of variables generated in every degree.

---

## Running Tests

```bash
cd DiffAlg
python -m pytest tests/test_diffalg.py -v
```

47 tests covering all operations, ported from the original M2 test suite plus additional coverage.

---

## Citation

```bibtex
@article{DiffAlgArticle,
  title   = {{DiffAlg: a Differential algebra package}},
  author  = {Manuel Dubinsky and Cesar Massri and Ariel Molinuevo and Federico Quallbrunn},
  journal = {The Journal of Software for Algebra and Geometry},
  volume  = {9},
  year    = {2019},
  doi     = {10.2140/jsag.2019.9.11}
}
```
