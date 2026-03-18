"""
Core implementation of DiffAlg for SageMath.

Ring architecture:
    - Coordinates:  x_0, ..., x_n  in a polynomial ring over QQ
    - Parameters:   a_0, a_1, ...  (scalar coefficients) in the same polynomial ring
    - Differentials: dx_0, ..., dx_n  generators of an ExteriorAlgebra over the
      polynomial ring of coordinates+parameters

A DiffAlgForm of type (n, r, d) is an r-form in (n+1) variables with
homogeneous polynomial coefficients of degree d.

A DiffAlgField of type (n, d) is a vector field in (n+1) variables with
homogeneous polynomial coefficients of degree d.
"""

from itertools import combinations
from sage.all import (
    QQ,
    PolynomialRing,
    ExteriorAlgebra,
    binomial,
    Combinations,
    matrix,
    vector,
    ZZ,
    SR,
)


# ─────────────────── Ring factory ───────────────────────────────

def _var_names(prefix, n):
    """Generate variable names: prefix_0, ..., prefix_n."""
    return [f"{prefix}{i}" for i in range(n + 1)]


def _make_rings(n, param_names=None):
    """
    Build the polynomial ring and exterior algebra for n+1 variables.

    Returns (poly_ring, ext_alg, x_list, dx_list)
    where x_list are the coordinate generators and dx_list the exterior gens.
    If param_names is given, the polynomial ring includes those as generators
    (they come BEFORE the x_i in the ring).
    """
    x_names = _var_names("x", n)
    if param_names:
        all_names = list(param_names) + x_names
    else:
        all_names = x_names
    R = PolynomialRing(QQ, all_names)
    dx_names = _var_names("dx", n)
    E = ExteriorAlgebra(R, dx_names)
    nparams = len(param_names) if param_names else 0
    x_list = [R.gen(nparams + i) for i in range(n + 1)]
    dx_list = list(E.gens())
    param_list = [R.gen(i) for i in range(nparams)] if param_names else []
    return R, E, x_list, dx_list, param_list


def _merge_rings(*elements):
    """
    Given DiffAlgForm / DiffAlgField objects, build a common ring that
    contains all their variables and parameters. Returns (R, E, x, dx, params).
    Reembeds all elements into this common ring.
    """
    max_n = max(e.n for e in elements)
    all_params = set()
    for e in elements:
        all_params.update(e.param_names)
    param_names = sorted(all_params)
    R, E, x, dx, params = _make_rings(max_n, param_names if param_names else None)
    return R, E, x, dx, params


def _reembed_form(form, R, E, x, dx):
    """Reembed a DiffAlgForm into a new exterior algebra E over R."""
    mc = form.element.monomial_coefficients()
    result = E.zero()
    old_R = form.poly_ring
    # Build substitution map: old ring gens -> new ring gens by name
    name_map = {str(g): g for g in R.gens()}
    for basis_key, coeff in mc.items():
        # coeff is in the old polynomial ring
        new_coeff = _sub_poly(coeff, old_R, R, name_map)
        # basis_key encodes which dx_i appear — reconstruct in new E
        new_basis = _reconstruct_basis(basis_key, E)
        result += new_coeff * new_basis
    return result


def _reembed_field_coeffs(field, R):
    """Reembed a DiffAlgField's polynomial coefficients into ring R."""
    old_R = field.poly_ring
    name_map = {str(g): g for g in R.gens()}
    new_coeffs = {}
    for i, p in field.coeffs.items():
        new_coeffs[i] = _sub_poly(p, old_R, R, name_map)
    return new_coeffs


def _sub_poly(poly, old_R, new_R, name_map=None):
    """Substitute a polynomial from old_R into new_R by matching variable names."""
    if name_map is None:
        name_map = {str(g): g for g in new_R.gens()}
    if poly in QQ:
        return new_R(poly)
    result = new_R.zero()
    for coeff, monom in poly:
        term = new_R(coeff)
        exp = monom.exponents()[0] if hasattr(monom, 'exponents') else None
        # Use the dict-based approach
        md = poly.dict() if hasattr(poly, 'dict') else None
    # Safer approach: use string-based substitution through the polynomial
    # Actually, let's use the hom approach
    old_gens = old_R.gens()
    images = [name_map.get(str(g), new_R.zero()) for g in old_gens]
    phi = old_R.hom(images, new_R)
    return phi(poly)


def _reconstruct_basis(basis_key, E):
    """
    Given a basis_key from monomial_coefficients() (a FrozenBitset-like thing),
    reconstruct the corresponding wedge product in E.
    """
    # basis_key is a FrozenBitset; iterating gives the set bit positions
    indices = sorted(basis_key)
    if not indices:
        return E.one()
    result = E.one()
    for i in indices:
        result = result * E.gen(i)
    return result


def _basis_key_to_indices(basis_key):
    """Extract sorted list of generator indices from a basis key."""
    return sorted(basis_key)


def _degree_of_basis_key(basis_key):
    """Form degree (number of dx_i) encoded in a basis key."""
    return len(list(basis_key))


# ─────────────────── Generic element constructors ───────────────

def _generic_param_names(var_name, count):
    """Generate parameter names: var_name0, var_name1, ..."""
    return [f"{var_name}{i}" for i in range(count)]


def _homogeneous_monomials(ring, variables, degree):
    """
    Return list of monomials of given degree in the specified variables
    (which must be generators of ring).
    """
    if degree == 0:
        return [ring.one()]
    n = len(variables)
    monoms = []
    for exp in _partitions_with_length(degree, n):
        m = ring.one()
        for i, e in enumerate(exp):
            m *= variables[i] ** e
        monoms.append(m)
    return monoms


def _partitions_with_length(d, n):
    """
    All n-tuples of non-negative integers summing to d.
    (compositions into n bins)
    """
    if n == 1:
        yield (d,)
        return
    for i in range(d + 1):
        for rest in _partitions_with_length(d - i, n - 1):
            yield (i,) + rest


# ─────────────────── DiffAlgForm ────────────────────────────────

class DiffAlgForm:
    """
    A differential r-form with homogeneous polynomial coefficients.

    Internally wraps an element of ExteriorAlgebra(R, [dx_0,...,dx_n])
    where R = QQ[params, x_0,...,x_n].
    """

    def __init__(self, element, poly_ring, ext_alg, n, param_names=None):
        self.element = element
        self.poly_ring = poly_ring
        self.ext_alg = ext_alg
        self.n = n  # space has n+1 coordinates x_0,...,x_n
        self.param_names = set(param_names) if param_names else set()

    @property
    def degree(self):
        """Return (n, r, d) where r = form degree, d = polynomial degree."""
        mc = self.element.monomial_coefficients()
        if not mc:
            return (self.n, 0, 0)
        r_max = 0
        d_max = 0
        for bk, coeff in mc.items():
            r = _degree_of_basis_key(bk)
            if r > r_max:
                r_max = r
            # polynomial degree in x_i only (not parameters)
            d = _poly_degree_in_x(coeff, self.poly_ring, self.n)
            if d > d_max:
                d_max = d
        return (self.n, r_max, d_max)

    @property
    def is_zero(self):
        return self.element == self.ext_alg.zero()

    def __repr__(self):
        if self.is_zero:
            return "0"
        return str(self.element)

    def __eq__(self, other):
        if isinstance(other, DiffAlgForm):
            R, E, x, dx, params = _merge_rings(self, other)
            e1 = _reembed_form(self, R, E, x, dx)
            e2 = _reembed_form(other, R, E, x, dx)
            return e1 == e2
        if other == 0:
            return self.is_zero
        return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        if isinstance(other, DiffAlgForm):
            R, E, x, dx, params = _merge_rings(self, other)
            e1 = _reembed_form(self, R, E, x, dx)
            e2 = _reembed_form(other, R, E, x, dx)
            pnames = self.param_names | other.param_names
            n = max(self.n, other.n)
            return DiffAlgForm(e1 + e2, R, E, n, pnames)
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, DiffAlgForm):
            return self + (-other)
        return NotImplemented

    def __neg__(self):
        return DiffAlgForm(-self.element, self.poly_ring, self.ext_alg,
                           self.n, self.param_names)

    def __mul__(self, other):
        """Scalar multiplication (right): form * scalar."""
        if isinstance(other, DiffAlgForm):
            return self.wedge(other)
        # Scalar
        return DiffAlgForm(self.element * self.poly_ring(other),
                           self.poly_ring, self.ext_alg, self.n, self.param_names)

    def __rmul__(self, other):
        """Scalar multiplication (left): scalar * form."""
        if isinstance(other, DiffAlgForm):
            return other.wedge(self)
        return DiffAlgForm(self.poly_ring(other) * self.element,
                           self.poly_ring, self.ext_alg, self.n, self.param_names)

    def __xor__(self, other):
        """Wedge product: w ^ h. Python's ^ operator."""
        if isinstance(other, DiffAlgForm):
            return self.wedge(other)
        return NotImplemented

    def __truediv__(self, scalar):
        """Division by scalar."""
        return self * QQ(1) / QQ(scalar)

    def wedge(self, other):
        """Exterior product of two forms."""
        R, E, x, dx, params = _merge_rings(self, other)
        e1 = _reembed_form(self, R, E, x, dx)
        e2 = _reembed_form(other, R, E, x, dx)
        pnames = self.param_names | other.param_names
        n = max(self.n, other.n)
        return DiffAlgForm(e1 * e2, R, E, n, pnames)

    def exterior_diff(self):
        """Exterior derivative d(self)."""
        mc = self.element.monomial_coefficients()
        E = self.ext_alg
        R = self.poly_ring
        n = self.n
        # x variables are the last n+1 generators of R
        nparams = len(R.gens()) - (n + 1)
        x_gens = [R.gen(nparams + i) for i in range(n + 1)]
        dx_gens = list(E.gens())

        result = E.zero()
        for bk, coeff in mc.items():
            # d(coeff * dx_I) = sum_j (d coeff/d x_j) dx_j ^ dx_I
            for j in range(n + 1):
                dc = coeff.derivative(x_gens[j])
                if dc != 0:
                    result += dc * dx_gens[j] * _reconstruct_basis(bk, E)
        return DiffAlgForm(result, R, E, n, self.param_names)

    def contract(self, field):
        """Contraction i_X(self) of this form with a vector field X."""
        if not isinstance(field, DiffAlgField):
            raise TypeError("Can only contract with a DiffAlgField")
        R, E, x, dx, params = _merge_rings(self, field)
        w = _reembed_form(self, R, E, x, dx)
        X_coeffs = _reembed_field_coeffs(field, R)
        n = max(self.n, field.n)

        mc = w.monomial_coefficients()
        result = E.zero()
        for bk, coeff in mc.items():
            indices = _basis_key_to_indices(bk)
            # i_X(dx_{i1} ^ ... ^ dx_{ir}) = sum_k (-1)^k X_{ik} dx_{i1}^...^hat{dx_{ik}}^...^dx_{ir}
            for k_pos, ik in enumerate(indices):
                X_ik = X_coeffs.get(ik, R.zero())
                if X_ik == 0:
                    continue
                sign = (-1) ** k_pos
                remaining = indices[:k_pos] + indices[k_pos + 1:]
                basis_rem = E.one()
                for idx in remaining:
                    basis_rem = basis_rem * E.gen(idx)
                result += sign * coeff * X_ik * basis_rem
        pnames = self.param_names | field.param_names
        return DiffAlgForm(result, R, E, n, pnames)

    def lie_derivative(self, field):
        """Lie derivative L_X(self) = i_X(d self) + d(i_X self) (Cartan formula)."""
        return self.exterior_diff().contract(field) + self.contract(field).exterior_diff()

    def pullback(self, morphism):
        """
        Pull-back by a morphism (list of DiffAlgForm 0-forms).
        morphism = [F_0, ..., F_n] where F_i are 0-forms (polynomials as forms).
        """
        R, E, x, dx, params = _merge_rings(self, *morphism)
        n = max(e.n for e in [self] + morphism)
        pnames = set()
        for e in [self] + morphism:
            pnames |= e.param_names

        # Build substitution: x_i -> F_i, dx_i -> dF_i
        nparams = len(R.gens()) - (n + 1)
        x_gens = [R.gen(nparams + i) for i in range(n + 1)]
        dx_gens = list(E.gens())

        # F_i as polynomials in R
        F_polys = []
        for fi in morphism:
            fi_embed = _reembed_form(fi, R, E, x, dx)
            mc = fi_embed.monomial_coefficients()
            # 0-form: key should be empty set
            for bk, coeff in mc.items():
                if _degree_of_basis_key(bk) == 0:
                    F_polys.append(coeff)
                    break
            else:
                F_polys.append(R.zero())

        # dF_i = sum_j (dF_i/dx_j) dx_j
        dF = []
        for fi_poly in F_polys:
            dfi = E.zero()
            for j in range(n + 1):
                dc = fi_poly.derivative(x_gens[j])
                if dc != 0:
                    dfi += dc * dx_gens[j]
            dF.append(dfi)

        # Now substitute in self
        w = _reembed_form(self, R, E, x, dx)
        mc = w.monomial_coefficients()
        result = E.zero()
        for bk, coeff in mc.items():
            indices = _basis_key_to_indices(bk)
            # coeff(F_0,...,F_n) * dF_{i1} ^ ... ^ dF_{ir}
            # First substitute x_i -> F_i in coeff
            sub_dict = {x_gens[i]: F_polys[i] for i in range(min(len(F_polys), n + 1))}
            new_coeff = coeff.subs(sub_dict)
            # Then wedge the dF's
            wedge_part = E.one()
            for idx in indices:
                if idx < len(dF):
                    wedge_part = wedge_part * dF[idx]
            result += new_coeff * wedge_part
        return DiffAlgForm(result, R, E, n, pnames)

    def singular_ideal(self):
        """
        Ideal generated by the polynomial coefficients (in x_i) of the form.
        Returns an ideal in QQ[params, x_0,...,x_n].
        """
        mc = self.element.monomial_coefficients()
        polys = [c for c in mc.values() if c != 0]
        if not polys:
            return self.poly_ring.ideal([self.poly_ring.zero()])
        return self.poly_ring.ideal(polys)

    def moduli_ideal(self):
        """
        Ideal generated by the scalar coefficients (parameters only).
        For a form with polynomial coefficients that are linear in parameters,
        extracts the scalar (parameter) coefficients.
        Returns an ideal in QQ[params].
        """
        mc = self.element.monomial_coefficients()
        n = self.n
        R = self.poly_ring
        nparams = len(R.gens()) - (n + 1)
        if nparams == 0:
            # No parameters — all coefficients are the scalars
            all_coeffs = []
            for poly in mc.values():
                all_coeffs.extend(_extract_all_coefficients(poly, R, n))
            if not all_coeffs:
                return R.ideal([R.zero()])
            return R.ideal(all_coeffs)

        # Extract coefficients of each x-monomial, which are polynomials in params
        param_ring = PolynomialRing(QQ, [str(R.gen(i)) for i in range(nparams)])
        all_scalars = []
        for poly in mc.values():
            scalars = _extract_param_coefficients(poly, R, n, param_ring)
            all_scalars.extend(scalars)
        if not all_scalars:
            return param_ring.ideal([param_ring.zero()])
        return param_ring.ideal(all_scalars)

    def is_homogeneous(self):
        """Test if the form has all terms of the same total degree (poly + form degree)."""
        mc = self.element.monomial_coefficients()
        if not mc:
            return True
        degrees = set()
        for bk, coeff in mc.items():
            r = _degree_of_basis_key(bk)
            d = _poly_degree_in_x(coeff, self.poly_ring, self.n)
            degrees.add(r + d)
        return len(degrees) <= 1

    def homogenize(self, new_var_index=None):
        """
        Homogenize the form by adding a new variable x_{n+1}.
        Returns a form in n+2 variables.
        """
        new_n = self.n + 1
        if new_var_index is None:
            new_var_index = new_n

        # Find max total degree in x
        mc = self.element.monomial_coefficients()
        if not mc:
            return self  # zero form

        max_deg = 0
        for bk, coeff in mc.items():
            d = _poly_degree_in_x(coeff, self.poly_ring, self.n)
            if d > max_deg:
                max_deg = d

        # Build new rings with n+2 variables
        param_names = sorted(self.param_names) if self.param_names else None
        R, E, x, dx, params = _make_rings(new_n, param_names)

        result = E.zero()
        old_R = self.poly_ring
        name_map = {str(g): g for g in R.gens()}
        for bk, coeff in mc.items():
            new_coeff = _sub_poly(coeff, old_R, R, name_map)
            d = _poly_degree_in_x(coeff, old_R, self.n)
            # Multiply by x_{new_n}^(max_deg - d)
            if max_deg - d > 0:
                new_coeff = new_coeff * x[new_var_index] ** (max_deg - d)
            indices = _basis_key_to_indices(bk)
            basis = E.one()
            for i in indices:
                basis = basis * E.gen(i)
            result += new_coeff * basis

        return DiffAlgForm(result, R, E, new_n, self.param_names)

    def projectivize(self):
        """
        Projectivize: homogenize and then ensure i_R(w) = 0 where R is the
        radial field. If the form already descends to projective space, just
        homogenize. Otherwise, apply the correction.
        """
        wh = self.homogenize()
        new_n = self.n + 1
        R_field = radial(new_n)
        iR_wh = wh.contract(R_field)
        if iR_wh.is_zero:
            return wh
        # Correction: x_{new_n} * wh - dx_{new_n} ^ iR_wh
        R, E, x, dx, params = _merge_rings(wh, iR_wh)
        e_wh = _reembed_form(wh, R, E, x, dx)
        e_iR = _reembed_form(iR_wh, R, E, x, dx)
        nparams = len(R.gens()) - (new_n + 1)
        x_new = R.gen(nparams + new_n)
        dx_new = E.gen(new_n)
        result = x_new * e_wh - dx_new * e_iR
        return DiffAlgForm(result, R, E, new_n, wh.param_names)

    def randomize(self, ring=ZZ, density=1.0, height=10):
        """
        Replace all parameter variables with random values.
        Returns a form with no parameters.
        """
        import random as pyrandom
        R = self.poly_ring
        n = self.n
        nparams = len(R.gens()) - (n + 1)
        if nparams == 0:
            return self

        # Build new ring without parameters
        x_names = _var_names("x", n)
        new_R = PolynomialRing(QQ, x_names)
        dx_names = _var_names("dx", n)
        new_E = ExteriorAlgebra(new_R, dx_names)

        # Random substitution for parameters
        max_attempts = 10
        for _ in range(max_attempts):
            sub_dict = {}
            for i in range(nparams):
                if pyrandom.random() <= density:
                    sub_dict[R.gen(i)] = R(pyrandom.randint(-height, height))
                else:
                    sub_dict[R.gen(i)] = R.zero()
            # Apply substitution
            mc = self.element.monomial_coefficients()
            result = new_E.zero()
            name_map = {str(g): g for g in new_R.gens()}
            for bk, coeff in mc.items():
                new_coeff = coeff.subs(sub_dict)
                # Now move to new_R
                nc = _sub_poly(new_coeff, R, new_R, name_map)
                basis = _reconstruct_basis(bk, new_E)
                result += nc * basis
            if result != new_E.zero():
                return DiffAlgForm(result, new_R, new_E, n)
        raise ValueError("Random substitution produced zero form after multiple attempts")


# ─────────────────── DiffAlgField ───────────────────────────────

class DiffAlgField:
    """
    A polynomial vector field X = sum P_i(x) * d/dx_i.

    Stored as a dict {i: P_i} over a polynomial ring QQ[params, x_0,...,x_n].
    """

    def __init__(self, coeffs, poly_ring, n, param_names=None):
        """
        coeffs: dict {i: polynomial} for X = sum coeffs[i] * d/dx_i
        poly_ring: the polynomial ring
        n: space dimension - 1
        """
        self.coeffs = {i: c for i, c in coeffs.items() if c != 0}
        self.poly_ring = poly_ring
        self.n = n
        self.param_names = set(param_names) if param_names else set()

    @property
    def degree(self):
        """Return (n, d) where d = polynomial degree of coefficients."""
        if not self.coeffs:
            return (self.n, 0)
        d_max = max(_poly_degree_in_x(c, self.poly_ring, self.n)
                     for c in self.coeffs.values())
        return (self.n, d_max)

    @property
    def is_zero(self):
        return all(c == 0 for c in self.coeffs.values())

    def __repr__(self):
        if self.is_zero:
            return "0"
        terms = []
        for i in sorted(self.coeffs.keys()):
            c = self.coeffs[i]
            if c == 0:
                continue
            c_str = str(c)
            if c == 1:
                terms.append(f"d/dx{i}")
            elif c == -1:
                terms.append(f"-d/dx{i}")
            else:
                if '+' in c_str or '-' in c_str[1:]:
                    terms.append(f"({c_str})*d/dx{i}")
                else:
                    terms.append(f"{c_str}*d/dx{i}")
        return " + ".join(terms).replace(" + -", " - ")

    def __eq__(self, other):
        if isinstance(other, DiffAlgField):
            R, E, x, dx, params = _merge_rings(self, other)
            c1 = _reembed_field_coeffs(self, R)
            c2 = _reembed_field_coeffs(other, R)
            n = max(self.n, other.n)
            for i in range(n + 1):
                if c1.get(i, R.zero()) != c2.get(i, R.zero()):
                    return False
            return True
        if other == 0:
            return self.is_zero
        return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        if isinstance(other, DiffAlgField):
            R, E, x, dx, params = _merge_rings(self, other)
            c1 = _reembed_field_coeffs(self, R)
            c2 = _reembed_field_coeffs(other, R)
            n = max(self.n, other.n)
            coeffs = {}
            for i in range(n + 1):
                s = c1.get(i, R.zero()) + c2.get(i, R.zero())
                if s != 0:
                    coeffs[i] = s
            pnames = self.param_names | other.param_names
            return DiffAlgField(coeffs, R, n, pnames)
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, DiffAlgField):
            return self + (-other)
        return NotImplemented

    def __neg__(self):
        return DiffAlgField({i: -c for i, c in self.coeffs.items()},
                            self.poly_ring, self.n, self.param_names)

    def __mul__(self, scalar):
        R = self.poly_ring
        s = R(scalar)
        return DiffAlgField({i: s * c for i, c in self.coeffs.items()},
                            R, self.n, self.param_names)

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __truediv__(self, scalar):
        return self * (QQ(1) / QQ(scalar))

    def __or__(self, other):
        """Lie bracket: X | Y."""
        if isinstance(other, DiffAlgField):
            return self.bracket(other)
        return NotImplemented

    def bracket(self, other):
        """Lie bracket [self, other]."""
        R, E, x, dx, params = _merge_rings(self, other)
        n = max(self.n, other.n)
        nparams = len(R.gens()) - (n + 1)
        x_gens = [R.gen(nparams + i) for i in range(n + 1)]

        X = _reembed_field_coeffs(self, R)
        Y = _reembed_field_coeffs(other, R)

        coeffs = {}
        for i in range(n + 1):
            val = R.zero()
            for j in range(n + 1):
                Xj = X.get(j, R.zero())
                Yi = Y.get(i, R.zero())
                Yj = Y.get(j, R.zero())
                Xi = X.get(i, R.zero())
                val += Xj * Yi.derivative(x_gens[j]) - Yj * Xi.derivative(x_gens[j])
            if val != 0:
                coeffs[i] = val
        pnames = self.param_names | other.param_names
        return DiffAlgField(coeffs, R, n, pnames)

    def is_homogeneous(self):
        """Test if all coefficients have the same polynomial degree in x."""
        if not self.coeffs:
            return True
        degrees = set()
        for c in self.coeffs.values():
            if c != 0:
                degrees.add(_poly_degree_in_x(c, self.poly_ring, self.n))
        return len(degrees) <= 1

    def randomize(self, ring=ZZ, density=1.0, height=10):
        """Replace all parameter variables with random values."""
        import random as pyrandom
        R = self.poly_ring
        n = self.n
        nparams = len(R.gens()) - (n + 1)
        if nparams == 0:
            return self

        x_names = _var_names("x", n)
        new_R = PolynomialRing(QQ, x_names)

        max_attempts = 10
        for _ in range(max_attempts):
            sub_dict = {}
            for i in range(nparams):
                if pyrandom.random() <= density:
                    sub_dict[R.gen(i)] = R(pyrandom.randint(-height, height))
                else:
                    sub_dict[R.gen(i)] = R.zero()
            new_coeffs = {}
            name_map = {str(g): g for g in new_R.gens()}
            for idx, c in self.coeffs.items():
                nc = c.subs(sub_dict)
                nc = _sub_poly(nc, R, new_R, name_map)
                if nc != 0:
                    new_coeffs[idx] = nc
            if new_coeffs:
                return DiffAlgField(new_coeffs, new_R, n)
        raise ValueError("Random substitution produced zero vector field")


# ─────────────────── DiffAlgDistribution ────────────────────────

class DiffAlgDistribution:
    """A distribution: a list of vector fields."""

    def __init__(self, fields):
        if not all(isinstance(f, DiffAlgField) for f in fields):
            raise TypeError("All elements must be DiffAlgField instances")
        self.fields = list(fields)

    def __repr__(self):
        return f"Distribution({', '.join(str(f) for f in self.fields)})"

    def is_involutive(self):
        """Test if the distribution is involutive (closed under Lie bracket)."""
        if len(self.fields) <= 1:
            return True
        # Merge all into common ring
        R, E, x, dx, params = _merge_rings(*self.fields)
        n = max(f.n for f in self.fields)
        nparams = len(R.gens()) - (n + 1)
        x_gens = [R.gen(nparams + i) for i in range(n + 1)]

        # Reembed all fields
        fields_re = [_reembed_field_coeffs(f, R) for f in self.fields]

        # Check each bracket [X_i, X_j] is in span of the fields
        for i in range(len(self.fields)):
            for j in range(len(self.fields)):
                br = self.fields[i].bracket(self.fields[j])
                br_coeffs = _reembed_field_coeffs(br, R)
                # Check if br is in span of fields — build matrix system
                if not _is_in_span(br_coeffs, fields_re, R, n):
                    return False
        return True

    def rank(self):
        """Compute the rank of the distribution."""
        if not self.fields:
            return 0
        R, E, x, dx, params = _merge_rings(*self.fields)
        n = max(f.n for f in self.fields)
        fields_re = [_reembed_field_coeffs(f, R) for f in self.fields]
        # Build matrix with coefficients
        from sage.all import matrix as sage_matrix
        frac = R.fraction_field()
        mat_rows = []
        for fc in fields_re:
            row = [frac(fc.get(i, R.zero())) for i in range(n + 1)]
            mat_rows.append(row)
        M = sage_matrix(frac, mat_rows)
        return M.rank()


def _is_in_span(target_coeffs, basis_coeffs_list, R, n):
    """Check if target vector field is in span of basis fields modulo ideal."""
    # Use fraction field for the check
    frac = R.fraction_field()
    m = len(basis_coeffs_list)
    # Build system: sum lambda_k * basis_k[i] = target[i] for all i
    # This is a linear system over the fraction field
    mat_rows = []
    rhs = []
    for i in range(n + 1):
        row = [frac(bc.get(i, R.zero())) for bc in basis_coeffs_list]
        mat_rows.append(row)
        rhs.append(frac(target_coeffs.get(i, R.zero())))
    from sage.all import matrix as sage_matrix, vector as sage_vector
    try:
        M = sage_matrix(frac, mat_rows)
        b = sage_vector(frac, rhs)
        # Check if b is in column space of M
        # M * lambda = b where lambda is m-vector
        # Transpose: M^T has rows that are the columns of M
        Mt = M.transpose()
        try:
            Mt.solve_left(b)
            return True
        except ValueError:
            return False
    except (TypeError, ValueError):
        return False


# ─────────────────── Constructor functions ──────────────────────

def new_form(n_or_expr, r=None, d=None, var_name=None):
    """
    Create a differential form.

    new_form(n, r, d, "a")  — generic r-form in (n+1) vars with degree-d coefficients
    new_form("x0*dx1 - x1*dx0") — parse a form from a string expression

    String syntax:
        Coordinates: x0, x1, ... (or x_0, x_1, ...)
        Differentials: dx0, dx1, ... (or dx_0, dx_1, ...)
        Any other symbols become parameters in the coefficient ring.
        Use * for multiplication, ^ for powers.
    """
    if isinstance(n_or_expr, (int, Integer_type())):
        n = int(n_or_expr)
        assert r is not None and d is not None and var_name is not None
        return _new_generic_form(n, r, d, var_name)
    elif isinstance(n_or_expr, str):
        return _parse_form_string(n_or_expr)
    else:
        raise TypeError(
            f"new_form expects int or string, got {type(n_or_expr)}"
        )


def Integer_type():
    from sage.rings.integer import Integer
    return Integer


import re as _re

def _normalize_var_names(expr):
    """Normalize x_0 -> x0, dx_0 -> dx0, ax_0 -> ax0, ^ -> ** in the expression."""
    expr = _re.sub(r'\bdx_(\d+)', r'dx\1', expr)
    expr = _re.sub(r'\bax_(\d+)', r'ax\1', expr)
    expr = _re.sub(r'\bx_(\d+)', r'x\1', expr)
    expr = expr.replace('^', '**')
    return expr


def _parse_expr_tokens(expr):
    """
    Parse a DiffAlg expression string. Returns (n, param_names, coord_names, diff_names).
    coord_names are like x0, x1, ... ; diff_names are dx0/ax0 etc.
    """
    # Find all identifiers
    tokens = set(_re.findall(r'\b([a-zA-Z_]\w*)\b', expr))
    # Classify tokens
    dx_indices = []
    ax_indices = []
    x_indices = []
    param_names = []
    for tok in tokens:
        m_dx = _re.fullmatch(r'dx(\d+)', tok)
        m_ax = _re.fullmatch(r'ax(\d+)', tok)
        m_x = _re.fullmatch(r'x(\d+)', tok)
        if m_dx:
            dx_indices.append(int(m_dx.group(1)))
        elif m_ax:
            ax_indices.append(int(m_ax.group(1)))
        elif m_x:
            x_indices.append(int(m_x.group(1)))
        else:
            param_names.append(tok)
    all_indices = x_indices + dx_indices + ax_indices
    n = max(all_indices) if all_indices else 0
    return n, sorted(param_names), dx_indices, ax_indices


def _parse_form_string(expr):
    """Parse a string into a DiffAlgForm."""
    expr = _normalize_var_names(expr)
    n, param_names, dx_indices, ax_indices = _parse_expr_tokens(expr)
    if ax_indices:
        raise ValueError(
            "Found ax_i tokens in a form expression. "
            "Use dx0, dx1, ... for differentials, or use new_field() for vector fields."
        )
    R, E, x, dx, params = _make_rings(n, param_names if param_names else None)
    # Build a local namespace for eval
    ns = {}
    for i in range(n + 1):
        ns[f'x{i}'] = R.gen((len(param_names) if param_names else 0) + i)
        ns[f'dx{i}'] = E.gen(i)
    for j, pname in enumerate(param_names):
        ns[pname] = R.gen(j)
    result = eval(expr, {"__builtins__": {}}, ns)  # noqa: S307
    if result.parent() is R:
        # It's a 0-form (function), wrap it into the exterior algebra
        result = E(result)
    return DiffAlgForm(result, R, E, n, set(param_names))


def _parse_field_string(expr):
    """Parse a string into a DiffAlgField."""
    expr = _normalize_var_names(expr)
    n, param_names, dx_indices, ax_indices = _parse_expr_tokens(expr)
    if dx_indices:
        raise ValueError(
            "Found dx_i tokens in a field expression. "
            "Use ax0, ax1, ... for partials, or use new_form() for differential forms."
        )
    R, E, x, dx, params = _make_rings(n, param_names if param_names else None)
    # Build a local namespace: x_i are polynomial gens, ax_i are symbolic markers
    # We use a temporary polynomial ring with ax_i as extra variables to parse,
    # then extract coefficients.
    ax_names = [f'ax{i}' for i in range(n + 1)]
    all_tmp_names = (param_names if param_names else []) + \
                    [f'x{i}' for i in range(n + 1)] + ax_names
    R_tmp = PolynomialRing(QQ, all_tmp_names)
    ns = {}
    nparams = len(param_names) if param_names else 0
    for j, pname in enumerate(param_names):
        ns[pname] = R_tmp.gen(j)
    for i in range(n + 1):
        ns[f'x{i}'] = R_tmp.gen(nparams + i)
        ns[f'ax{i}'] = R_tmp.gen(nparams + n + 1 + i)
    parsed = eval(expr, {"__builtins__": {}}, ns)  # noqa: S307
    # Extract coefficient of each ax_i
    coeffs = {}
    nparams_r = len(param_names) if param_names else 0
    name_map = {str(g): g for g in R.gens()}
    for i in range(n + 1):
        ax_var = R_tmp.gen(nparams + n + 1 + i)
        # Coefficient of ax_i in the parsed polynomial
        c = parsed.coefficient(ax_var)
        if c != 0:
            # Re-embed c into R (it's in R_tmp but only uses x_j and param vars)
            coeffs[i] = _sub_poly(c, R_tmp, R, name_map)
    return DiffAlgField(coeffs, R, n, set(param_names))


def _new_generic_form(n, r, d, var_name):
    """
    Create a generic homogeneous r-form in (n+1) variables with degree-d
    polynomial coefficients. Parameters are named var_name0, var_name1, ...
    """
    num_x_monomials = int(binomial(n + d, d))
    num_dx_monomials = int(binomial(n + 1, r))
    num_params = num_x_monomials * num_dx_monomials
    param_names = _generic_param_names(var_name, num_params)

    R, E, x, dx, params = _make_rings(n, param_names)

    # Build the generic form
    x_monoms = _homogeneous_monomials(R, x, d)
    dx_combos = list(combinations(range(n + 1), r))

    result = E.zero()
    param_idx = 0
    for dx_combo in dx_combos:
        dx_prod = E.one()
        for k in dx_combo:
            dx_prod = dx_prod * dx[k]
        for x_mon in x_monoms:
            result += params[param_idx] * x_mon * dx_prod
            param_idx += 1

    return DiffAlgForm(result, R, E, n, set(param_names))


def new_field(n_or_expr, d=None, var_name=None):
    """
    Create a vector field.

    new_field(n, d, "a")  — generic field in (n+1) vars with degree-d coefficients
    new_field("x0^2*ax1 + x1*ax0") — parse a field from a string expression

    String syntax:
        Coordinates: x0, x1, ... (or x_0, x_1, ...)
        Partials: ax0, ax1, ... (or ax_0, ax_1, ...) for ∂/∂x_i
        Any other symbols become parameters in the coefficient ring.
        Use * for multiplication, ^ for powers.
    """
    if isinstance(n_or_expr, (int, Integer_type())):
        n = int(n_or_expr)
        assert d is not None and var_name is not None
        return _new_generic_field(n, d, var_name)
    elif isinstance(n_or_expr, str):
        return _parse_field_string(n_or_expr)
    else:
        raise TypeError(
            f"new_field expects int or string, got {type(n_or_expr)}"
        )


def _new_generic_field(n, d, var_name):
    """
    Create a generic vector field in (n+1) variables with degree-d
    polynomial coefficients.
    """
    num_x_monomials = int(binomial(n + d, d))
    num_params = num_x_monomials * (n + 1)
    param_names = _generic_param_names(var_name, num_params)

    R, E, x, dx, params = _make_rings(n, param_names)

    coeffs = {}
    param_idx = 0
    x_monoms = _homogeneous_monomials(R, x, d)
    for i in range(n + 1):
        val = R.zero()
        for x_mon in x_monoms:
            val += params[param_idx] * x_mon
            param_idx += 1
        if val != 0:
            coeffs[i] = val
    return DiffAlgField(coeffs, R, n, set(param_names))


def radial(n):
    """
    The radial (Euler) vector field R = sum x_i * d/dx_i in (n+1) variables.
    """
    R, E, x, dx, params = _make_rings(n)
    coeffs = {i: x[i] for i in range(n + 1)}
    return DiffAlgField(coeffs, R, n)


def linear_comb(elements, var_name):
    """
    Generic linear combination of a list of forms or fields.
    Returns sum(var_name_i * elements[i]).
    """
    if not elements:
        raise ValueError("Empty list")

    k = len(elements)
    param_names = _generic_param_names(var_name, k)

    if isinstance(elements[0], DiffAlgForm):
        all_elems = elements
        max_n = max(e.n for e in elements)
        all_params = set()
        for e in elements:
            all_params.update(e.param_names)
        all_params.update(param_names)

        R, E, x, dx, params = _make_rings(max_n, sorted(all_params))
        result = E.zero()
        # Find the new parameters by name
        name_to_gen = {str(g): g for g in R.gens()}
        for i, elem in enumerate(elements):
            lam = name_to_gen[param_names[i]]
            e_re = _reembed_form(elem, R, E, x, dx)
            result += lam * e_re
        return DiffAlgForm(result, R, E, max_n, all_params)

    elif isinstance(elements[0], DiffAlgField):
        max_n = max(e.n for e in elements)
        all_params = set()
        for e in elements:
            all_params.update(e.param_names)
        all_params.update(param_names)

        R, E, x, dx, params = _make_rings(max_n, sorted(all_params))
        name_to_gen = {str(g): g for g in R.gens()}
        coeffs = {}
        for i, elem in enumerate(elements):
            lam = name_to_gen[param_names[i]]
            ec = _reembed_field_coeffs(elem, R)
            for j, c in ec.items():
                coeffs[j] = coeffs.get(j, R.zero()) + lam * c
        return DiffAlgField(coeffs, R, max_n, all_params)
    else:
        raise TypeError("Elements must be DiffAlgForm or DiffAlgField")


def logarithmic_form(n, degrees, var_name, projective=False):
    """
    Create a logarithmic form of type (d_0,...,d_k) in (n+1) variables.

    A logarithmic form w = (prod f_i) * sum (a_i * df_i / f_i)
    where f_i is a generic polynomial of degree degrees[i].

    If projective=True, impose that the form descends to projective space.
    """
    k = len(degrees)

    # Create generic polynomials f_0, ..., f_{k-1}
    F_forms = []
    for i in range(k):
        F_forms.append(new_form(n, 0, degrees[i], f"{var_name}{i + 1}"))

    # Build the logarithmic form:
    # w = sum_i a_i * (prod_{j≠i} f_j) * df_i
    coeff_param_names = _generic_param_names(f"{var_name}0", k)
    # We need the a_i parameters

    # Actually, let's build it step by step merging rings
    # First create the a_i as 0-forms
    all_forms = list(F_forms)
    # Build products and differentials
    terms = []
    for i in range(k):
        # Product of all f_j except f_i
        prod_form = None
        for j in range(k):
            if j != i:
                if prod_form is None:
                    prod_form = F_forms[j]
                else:
                    prod_form = prod_form.wedge(F_forms[j])
        # df_i
        dfi = F_forms[i].exterior_diff()
        # a_i * prod * df_i
        if prod_form is not None:
            term = prod_form.wedge(dfi)
        else:
            term = dfi
        terms.append(term)

    # Sum with coefficients a_0, ..., a_{k-1}
    # Build common ring including the a_i
    result = terms[0]
    for i in range(1, k):
        result = result + terms[i]

    # Actually, we want a_i coefficients. Let's use linear_comb
    result = linear_comb(terms, f"{var_name}0")

    if projective:
        # Impose i_R(w) = 0: substitute a_0 = -(sum a_i * d_i/d_0 for i>0)
        R_field = radial(n)
        # For logarithmic forms, the projective condition gives
        # a_0 = -(a_1*d_1/d_0 + a_2*d_2/d_0 + ... + a_{k-1}*d_{k-1}/d_0)
        # We substitute in the parameter ring
        R = result.poly_ring
        name_map = {str(g): g for g in R.gens()}
        a0_name = f"{var_name}00"
        if a0_name in name_map:
            a0_gen = name_map[a0_name]
            sub_expr = R.zero()
            for i in range(1, k):
                ai_name = f"{var_name}0{i}"
                if ai_name in name_map:
                    ai_gen = name_map[ai_name]
                    sub_expr -= ai_gen * QQ(degrees[i]) / QQ(degrees[0])
            # Apply substitution
            mc = result.element.monomial_coefficients()
            E = result.ext_alg
            new_elem = E.zero()
            for bk, coeff in mc.items():
                new_coeff = coeff.subs({a0_gen: sub_expr})
                new_elem += new_coeff * _reconstruct_basis(bk, E)
            # Remove a0 from param names
            new_params = result.param_names - {a0_name}
            result = DiffAlgForm(new_elem, R, E, result.n, new_params)

    return result


def dist(fields):
    """Create a DiffAlgDistribution from a list of vector fields."""
    return DiffAlgDistribution(fields)


# ─────────────────── genKer / genIm ─────────────────────────────

def gen_ker(expr, var):
    """
    Compute a basis of the kernel of a linear expression.

    Given `expr` linear in the parameters of `var`, find all values of `var`
    such that `expr = 0`.

    Parameters
    ----------
    expr : DiffAlgForm or DiffAlgField
        Expression that is linear in the parameters of var.
    var : DiffAlgForm or DiffAlgField
        The variable with free parameters (b_0, ..., b_{m-1}).

    Returns
    -------
    list : Basis of the kernel (matching var's type).
           If expr has a constant part (non-homogeneous in var's params),
           returns [basis_list, particular_solution].
    """
    from sage.all import matrix as sage_matrix, vector as sage_vector
    from collections import defaultdict

    R, E, x, dx, params = _merge_rings(expr, var)
    n = max(expr.n, var.n)
    nparams_total = len(R.gens()) - (n + 1)

    var_params = sorted(var.param_names)
    m = len(var_params)
    if m == 0:
        is_zero = (isinstance(expr, DiffAlgForm) and expr.is_zero) or \
                  (isinstance(expr, DiffAlgField) and expr.is_zero)
        return [var] if is_zero else []

    name_map = {str(g): g for g in R.gens()}
    var_gens = [name_map[p] for p in var_params]
    var_gen_idx_set = set()
    for vg in var_gens:
        for i, g in enumerate(R.gens()):
            if g == vg:
                var_gen_idx_set.add(i)
                break

    other_param_names = sorted(
        (expr.param_names | var.param_names) - set(var_params)
    )
    has_other = len(other_param_names) > 0

    # ── Step 1: extract polynomial coefficients from expr ──
    if isinstance(expr, DiffAlgForm):
        expr_elem = _reembed_form(expr, R, E, x, dx)
        poly_coeffs = list(expr_elem.monomial_coefficients().values())
    else:
        expr_coeffs = _reembed_field_coeffs(expr, R)
        poly_coeffs = list(expr_coeffs.values())

    # ── Step 2: decompose linearly in var_gens ──
    # For each polynomial P: P = sum_j b_j * dP/db_j + P(b=0)
    # dP/db_j gives the "A column" (doesn't involve b's since P is linear in b)
    # P(b=0) gives the "constant" part C

    sub_zero = {vg: R.zero() for vg in var_gens}

    # For each P_k and each b_j: derivative = coefficient of b_j
    A_polys = []   # A_polys[j] = [dP_0/db_j, dP_1/db_j, ...]
    for j in range(m):
        A_polys.append([P.derivative(var_gens[j]) for P in poly_coeffs])
    C_polys = [P.subs(sub_zero) for P in poly_coeffs]

    # ── Step 3: flatten to scalar equations ──
    # Group each polynomial by x-exponent to extract scalar equations.
    # Scalar coefficients live in QQ (no other params) or QQ[other_params].

    def group_by_x(poly):
        """Group poly by x-exponent. Returns {x_exp: param_poly_in_R}."""
        if poly == 0:
            return {}
        result = defaultdict(lambda: R.zero())
        for exp, coeff in poly.dict().items():
            x_exp = tuple(exp[nparams_total:])
            param_monom = R.one()
            for i in range(nparams_total):
                if exp[i] > 0:
                    param_monom *= R.gen(i) ** exp[i]
            result[x_exp] += QQ(coeff) * param_monom
        return dict(result)

    # Flatten A columns and C into {(poly_idx, x_exp) → value}
    A_flat = []
    for j in range(m):
        flat_j = {}
        for k, P in enumerate(A_polys[j]):
            for x_exp, val in group_by_x(P).items():
                key = (k, x_exp)
                flat_j[key] = flat_j.get(key, R.zero()) + val
        A_flat.append(flat_j)

    C_flat = {}
    for k, P in enumerate(C_polys):
        for x_exp, val in group_by_x(P).items():
            key = (k, x_exp)
            C_flat[key] = C_flat.get(key, R.zero()) + val

    # Collect all equation keys
    all_keys = set()
    for af in A_flat:
        all_keys |= set(af.keys())
    all_keys |= set(C_flat.keys())
    sorted_keys = sorted(all_keys, key=str)

    # ── Step 4: build the coefficient field and matrix ──
    if has_other:
        S = PolynomialRing(QQ, other_param_names)
        K = S.fraction_field()
        # Hom from R → S (maps other_params to S gens, everything else to 0)
        S_name_map = {str(g): g for g in S.gens()}
        def to_K(val):
            if val == 0:
                return K.zero()
            if val in QQ:
                return K(val)
            s_val = S.zero()
            for exp, c in val.dict().items():
                term = S(QQ(c))
                for i in range(nparams_total):
                    if exp[i] > 0:
                        gn = str(R.gen(i))
                        if gn in S_name_map:
                            term *= S_name_map[gn] ** exp[i]
                s_val += term
            return K(s_val)
    else:
        K = QQ
        def to_K(val):
            if val == 0:
                return K.zero()
            try:
                return K(val)
            except (TypeError, ValueError):
                return K.zero()

    eq_rows = []
    eq_rhs = []
    for key in sorted_keys:
        row = [to_K(A_flat[j].get(key, R.zero())) for j in range(m)]
        c_val = to_K(C_flat.get(key, R.zero()))
        if any(r != K.zero() for r in row) or c_val != K.zero():
            eq_rows.append(row)
            eq_rhs.append(-c_val)

    # ── Step 5: solve the linear system ──
    if not eq_rows:
        # No constraints: full space is kernel
        return _basis_from_params(var, var_params, var_gens, R, E, n)

    M = sage_matrix(K, eq_rows)
    b_vec = sage_vector(K, eq_rhs)

    ker = M.right_kernel()
    has_particular = any(r != K.zero() for r in eq_rhs)

    # ── Step 6: reconstruct basis elements ──
    # Reembed var in merged ring
    if isinstance(var, DiffAlgForm):
        var_mc = _reembed_form(var, R, E, x, dx).monomial_coefficients()
    else:
        var_coeffs_re = _reembed_field_coeffs(var, R)

    def _k_to_R(val):
        """Convert K element to R polynomial."""
        if K is QQ:
            return R(val)
        # val in Frac(S): clear denominator
        num = val.numerator()
        den = val.denominator()
        num_R = _S_to_R(num, S, R, name_map)
        den_R = _S_to_R(den, S, R, name_map)
        if den_R == R.one():
            return num_R
        if den_R in QQ and den_R != 0:
            return num_R * R(QQ(1) / QQ(den_R))
        raise ValueError(f"Non-polynomial kernel element: {val}")

    def reconstruct(vec):
        """Reconstruct form/field from a solution vector in K^m."""
        # Clear denominators if working over Frac(S)
        if has_other:
            from sage.all import lcm
            denoms = [v.denominator() for v in vec if v != 0]
            if denoms:
                common_den = denoms[0]
                for d in denoms[1:]:
                    common_den = lcm(common_den, d)
                vec = [v * K(common_den) for v in vec]

        sub = {}
        for j in range(m):
            sub[var_gens[j]] = _k_to_R(vec[j])

        pnames = set(other_param_names) if other_param_names else set()
        if isinstance(var, DiffAlgForm):
            new_elem = E.zero()
            for bk, coeff in var_mc.items():
                nc = coeff.subs(sub)
                if nc != 0:
                    new_elem += nc * _reconstruct_basis(bk, E)
            return DiffAlgForm(new_elem, R, E, n, pnames)
        else:
            new_coeffs = {}
            for i, c in var_coeffs_re.items():
                nc = c.subs(sub)
                if nc != 0:
                    new_coeffs[i] = nc
            return DiffAlgField(new_coeffs, R, n, pnames)

    basis = [reconstruct(v) for v in ker.basis()]

    if has_particular:
        zero_elem = _make_zero(var, R, E, n)
        try:
            x_part = M.solve_right(b_vec)
            particular = reconstruct(x_part)
            return [basis if basis else [zero_elem], particular]
        except ValueError:
            # Non-homogeneous system with no solution
            return [basis if basis else [zero_elem], None]

    if not basis:
        return [_make_zero(var, R, E, n)]
    return basis


def _S_to_R(poly_S, S, R, name_map):
    """Convert polynomial in S = QQ[other_params] to R."""
    if S is QQ or poly_S in QQ:
        return R(poly_S)
    images = [name_map.get(str(g), R.zero()) for g in S.gens()]
    phi = S.hom(images, R)
    return phi(poly_S)


def _basis_from_params(var, var_params, var_gens, R, E, n):
    """Full basis: one element per free parameter."""
    name_map = {str(g): g for g in R.gens()}
    if isinstance(var, DiffAlgForm):
        var_mc = _reembed_form(var, R, E, [R.gen(len(R.gens()) - n - 1 + i) for i in range(n+1)],
                               list(E.gens())).monomial_coefficients()
    else:
        var_mc = None
        var_coeffs_re = _reembed_field_coeffs(var, R)

    basis = []
    for fp_idx, fp in enumerate(var_params):
        sub = {var_gens[j]: (R.one() if j == fp_idx else R.zero())
               for j in range(len(var_gens))}
        if isinstance(var, DiffAlgForm):
            new_elem = E.zero()
            for bk, coeff in var_mc.items():
                nc = coeff.subs(sub)
                if nc != 0:
                    new_elem += nc * _reconstruct_basis(bk, E)
            if new_elem != E.zero():
                basis.append(DiffAlgForm(new_elem, R, E, n))
        else:
            new_coeffs = {}
            for i, c in var_coeffs_re.items():
                nc = c.subs(sub)
                if nc != 0:
                    new_coeffs[i] = nc
            if new_coeffs:
                basis.append(DiffAlgField(new_coeffs, R, n))
    return basis


def _make_zero(var, R, E, n):
    """Create a zero element matching var's type."""
    if isinstance(var, DiffAlgForm):
        return DiffAlgForm(E.zero(), R, E, n)
    else:
        return DiffAlgField({}, R, n)


def gen_im(expr, var):
    """
    Compute a basis of the image of a linear expression.

    Given `expr` which is linear and homogeneous in the parameters of `var`,
    compute the image {expr(v) : v ranges over all values of var}.

    Parameters
    ----------
    expr : DiffAlgForm or DiffAlgField
        Homogeneous expression linear in var's parameters.
    var : DiffAlgForm or DiffAlgField
        The variable with free parameters.

    Returns
    -------
    list : Basis of the image (matching expr's type).
    """
    from sage.all import matrix as sage_matrix

    R, E, x, dx, params = _merge_rings(expr, var)
    n = max(expr.n, var.n)
    var_params = sorted(var.param_names)
    m = len(var_params)
    name_map = {str(g): g for g in R.gens()}
    var_gens = [name_map[p] for p in var_params]
    pnames = (expr.param_names | var.param_names) - set(var_params)

    # Image = span of {expr(b_j = delta_{j,k}) : k = 0..m-1}
    # i.e., evaluate expr at each "unit vector" in parameter space
    if isinstance(expr, DiffAlgForm):
        expr_re = _reembed_form(expr, R, E, x, dx)
        mc = expr_re.monomial_coefficients()
        images = []
        for k in range(m):
            sub = {var_gens[j]: (R.one() if j == k else R.zero())
                   for j in range(m)}
            img = E.zero()
            for bk, coeff in mc.items():
                nc = coeff.subs(sub)
                if nc != 0:
                    img += nc * _reconstruct_basis(bk, E)
            images.append(img)

        # Find linearly independent subset using matrix rank
        # Vectorize each image
        all_keys = set()
        img_dicts = []
        for img in images:
            d = {}
            if img != E.zero():
                for bk, coeff in img.monomial_coefficients().items():
                    for exp, c in coeff.dict().items():
                        key = (tuple(sorted(bk)), exp)
                        d[key] = d.get(key, QQ.zero()) + QQ(c)
            img_dicts.append(d)
            all_keys |= set(d.keys())

        if not all_keys:
            return [DiffAlgForm(E.zero(), R, E, n, pnames)]

        sorted_keys = sorted(all_keys, key=str)
        mat_rows = []
        for d in img_dicts:
            mat_rows.append([QQ(d.get(k, 0)) for k in sorted_keys])

        M = sage_matrix(QQ, mat_rows)
        pivots = M.pivot_rows()
        basis = []
        for idx in pivots:
            basis.append(DiffAlgForm(images[idx], R, E, n, pnames))
        return basis if basis else [DiffAlgForm(E.zero(), R, E, n, pnames)]

    elif isinstance(expr, DiffAlgField):
        expr_coeffs = _reembed_field_coeffs(expr, R)
        images = []
        for k in range(m):
            sub = {var_gens[j]: (R.one() if j == k else R.zero())
                   for j in range(m)}
            new_coeffs = {}
            for i, c in expr_coeffs.items():
                nc = c.subs(sub)
                if nc != 0:
                    new_coeffs[i] = nc
            images.append(new_coeffs)

        # Remove zero images and find independent subset
        all_keys = set()
        for img in images:
            for i, c in img.items():
                if c != 0:
                    for exp, coeff in c.dict().items():
                        all_keys.add((i, exp))

        if not all_keys:
            return [DiffAlgField({}, R, n, pnames)]

        sorted_keys = sorted(all_keys, key=str)
        mat_rows = []
        for img in images:
            row = []
            for (i, exp) in sorted_keys:
                c = img.get(i, R.zero())
                if c == 0:
                    row.append(QQ.zero())
                else:
                    d = c.dict()
                    row.append(QQ(d.get(exp, 0)))
            mat_rows.append(row)

        M = sage_matrix(QQ, mat_rows)
        pivots = M.pivot_rows()
        basis = []
        for idx in pivots:
            basis.append(DiffAlgField(images[idx], R, n, pnames))
        return basis if basis else [DiffAlgField({}, R, n, pnames)]


# ─────────────────── Helper: build forms/fields from ring elts ──

def form_from_dict(n, terms, param_names=None):
    """
    Build a DiffAlgForm from a dict of {(i1,...,ir): polynomial}.
    Example: form_from_dict(2, {(0,): x0, (1,): -x1})  for x0*dx0 - x1*dx1.

    The polynomials should be elements of the same polynomial ring.
    """
    if not terms:
        R, E, x, dx, params = _make_rings(n, sorted(param_names) if param_names else None)
        return DiffAlgForm(E.zero(), R, E, n, param_names)

    # Infer the polynomial ring from coefficients
    first_val = next(iter(terms.values()))
    R_orig = first_val.parent()

    R, E, x, dx, params = _make_rings(n, sorted(param_names) if param_names else None)
    name_map = {str(g): g for g in R.gens()}

    result = E.zero()
    for dx_indices, coeff in terms.items():
        if isinstance(dx_indices, int):
            dx_indices = (dx_indices,)
        basis = E.one()
        for i in dx_indices:
            basis = basis * E.gen(i)
        if coeff.parent() is R:
            result += coeff * basis
        else:
            result += _sub_poly(coeff, R_orig, R, name_map) * basis

    return DiffAlgForm(result, R, E, n, param_names)


def field_from_dict(n, coeffs, param_names=None):
    """
    Build a DiffAlgField from a dict {i: polynomial}.
    Example: field_from_dict(2, {0: x0, 1: x1, 2: x2}) for the radial field.
    """
    if not coeffs:
        R, E, x, dx, params = _make_rings(n, sorted(param_names) if param_names else None)
        return DiffAlgField({}, R, n, param_names)

    first_val = next(iter(coeffs.values()))
    R_orig = first_val.parent()
    R, E, x, dx, params = _make_rings(n, sorted(param_names) if param_names else None)
    name_map = {str(g): g for g in R.gens()}

    new_coeffs = {}
    for i, c in coeffs.items():
        if c.parent() is R:
            new_coeffs[i] = c
        else:
            new_coeffs[i] = _sub_poly(c, R_orig, R, name_map)

    return DiffAlgField(new_coeffs, R, n, param_names)


# ─────────────────── Helper utilities ───────────────────────────

def _poly_degree_in_x(poly, R, n):
    """
    Degree of a polynomial only counting the x_0,...,x_n variables.
    Parameters (first gens of R) are ignored.
    """
    if poly == 0:
        return 0
    nparams = len(R.gens()) - (n + 1)
    # Get exponent vectors
    max_deg = 0
    for exp in poly.exponents():
        # exp is a tuple; entries nparams..nparams+n are x_0..x_n
        x_deg = sum(exp[nparams:nparams + n + 1])
        if x_deg > max_deg:
            max_deg = x_deg
    return max_deg


def _extract_all_coefficients(poly, R, n):
    """Extract all scalar (QQ) coefficients from a polynomial."""
    return list(poly.coefficients())


def _extract_param_coefficients(poly, R, n, param_ring):
    """
    Given a polynomial in R = QQ[params, x_0,...,x_n],
    extract the coefficient of each x-monomial (which is a polynomial in params).
    """
    nparams = len(R.gens()) - (n + 1)
    x_vars = [R.gen(nparams + i) for i in range(n + 1)]

    # Collect coefficients of all x-monomials
    result = []
    # Use the .coefficient() method for each monomial
    # Or iterate over terms and group by x-part
    if poly == 0:
        return []
    for c, m in poly:
        # c is in QQ, m is a monomial
        exp = m.exponents()[0]
        # The param part
        param_part = param_ring.one()
        for i in range(nparams):
            if exp[i] > 0:
                param_part *= param_ring.gen(i) ** exp[i]
        result.append(QQ(c) * param_part)
    # Group by x-monomial → sum the param contributions
    from collections import defaultdict
    grouped = defaultdict(lambda: param_ring.zero())
    for term_c, term_m in poly:
        exp = term_m.exponents()[0]
        x_key = tuple(exp[nparams:nparams + n + 1])
        param_part = param_ring.one()
        for i in range(nparams):
            if exp[i] > 0:
                param_part *= param_ring.gen(i) ** exp[i]
        grouped[x_key] += QQ(term_c) * param_part
    return list(grouped.values())


def _extract_monomial_coefficients(poly, R, n):
    """
    Extract coefficients of each x-monomial from a polynomial in R.
    Returns a list of polynomials in R (the coefficients grouped by x-part).
    """
    if poly == 0:
        return []
    nparams = len(R.gens()) - (n + 1)

    from collections import defaultdict
    grouped = defaultdict(lambda: R.zero())
    for c, m in poly:
        exp = m.exponents()[0]
        x_key = tuple(exp[nparams:nparams + n + 1])
        # Reconstruct the parameter part as element of R
        param_monomial = R.one()
        for i in range(nparams):
            if exp[i] > 0:
                param_monomial *= R.gen(i) ** exp[i]
        grouped[x_key] += QQ(c) * param_monomial

    return list(grouped.values())
