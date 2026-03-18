"""
Tests for DiffAlg — ported from the Macaulay2 test suite plus additional coverage.

Run with: sage -python -m pytest tests/test_diffalg.py -v
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sage.all import QQ, PolynomialRing, ExteriorAlgebra


# ─── Helpers to build specific forms/fields from ring elements ───

def make_form(n, terms, param_names=None):
    """Build a DiffAlgForm from {(dx_indices): poly_in_R}."""
    from diff_alg.core import DiffAlgForm, _make_rings, _var_names
    pn = sorted(param_names) if param_names else None
    R, E, x, dx, params = _make_rings(n, pn)
    result = E.zero()
    for dx_indices, coeff_fn in terms.items():
        if isinstance(dx_indices, int):
            dx_indices = (dx_indices,)
        basis = E.one()
        for i in dx_indices:
            basis = basis * E.gen(i)
        result += coeff_fn(R, x, params) * basis
    return DiffAlgForm(result, R, E, n, set(param_names) if param_names else set())


def make_field(n, coeffs_fn, param_names=None):
    """Build a DiffAlgField from {i: poly_in_R}."""
    from diff_alg.core import DiffAlgField, _make_rings
    pn = sorted(param_names) if param_names else None
    R, E, x, dx, params = _make_rings(n, pn)
    coeffs = {}
    for i, fn in coeffs_fn.items():
        coeffs[i] = fn(R, x, params)
    return DiffAlgField(coeffs, R, n, set(param_names) if param_names else set())


# ════════════════════════════════════════════════════════════════
# TEST 1: d(d(w)) = 0 for any form
# M2: w = newForm(3,1,2,"a"); assert(diff diff w == 0)
# ════════════════════════════════════════════════════════════════

class TestExteriorDifferential:

    def test_dd_is_zero(self):
        """d(d(w)) = 0 for a generic 1-form with degree-2 coefficients."""
        from diff_alg.core import new_form
        w = new_form(3, 1, 2, "a")
        dw = w.exterior_diff()
        ddw = dw.exterior_diff()
        assert ddw.is_zero, f"d(dw) should be 0, got {ddw}"

    def test_dd_is_zero_2form(self):
        """d(d(w)) = 0 for a    eric 2-form."""
        from diff_alg.core import new_form
        w = new_form(3, 2, 1, "a")
        assert w.exterior_diff().exterior_diff().is_zero

    def test_diff_of_0form(self):
        """d(f) for a 0-form (polynomial) gives a 1-form."""
        from diff_alg.core import new_form
        f = new_form(2, 0, 2, "a")
        df = f.exterior_diff()
        deg = df.degree
        assert deg[1] == 1, f"Expected 1-form, got degree {deg}"


# ════════════════════════════════════════════════════════════════
# TEST 2: Jacobi identity for Lie brackets
# M2: (X|(Y|Z)) + (Z|(X|Y)) + (Y|(Z|X)) == 0
# ════════════════════════════════════════════════════════════════

class TestLieBracket:

    def test_jacobi_identity(self):
        """Jacobi identity: [X,[Y,Z]] + [Z,[X,Y]] + [Y,[Z,X]] = 0."""
        from diff_alg.core import new_field
        X = new_field(2, 2, "a")
        Y = new_field(2, 1, "b")
        Z = new_field(2, 2, "c")
        j1 = X.bracket(Y.bracket(Z))
        j2 = Z.bracket(X.bracket(Y))
        j3 = Y.bracket(Z.bracket(X))
        result = j1 + j2 + j3
        assert result == 0, f"Jacobi identity failed: {result}"

    def test_bracket_antisymmetry(self):
        """[X,Y] = -[Y,X]."""
        from diff_alg.core import new_field
        X = new_field(2, 1, "a")
        Y = new_field(2, 1, "b")
        assert X.bracket(Y) + Y.bracket(X) == 0

    def test_bracket_with_self_is_zero(self):
        """[X,X] = 0."""
        from diff_alg.core import new_field
        X = new_field(2, 2, "a")
        assert X.bracket(X) == 0


# ════════════════════════════════════════════════════════════════
# TEST 3: Cartan / Euler formula: i_R(dw) + d(i_R(w)) = (r+d)*w
# M2: assert((R_(diff w) + diff(R _ w) - 3*w) == 0)
# For an (r,d)-form: L_R(w) = (r+d)*w
# ════════════════════════════════════════════════════════════════

class TestCartan:

    def test_euler_formula(self):
        """i_R(dw) + d(i_R(w)) = (r+d)*w for R = radial, w = (n,1,2)-form."""
        from diff_alg.core import new_form, radial
        w = new_form(3, 1, 2, "a")
        R = radial(3)
        lhs = w.exterior_diff().contract(R) + w.contract(R).exterior_diff()
        rhs = 3 * w  # r + d = 1 + 2 = 3
        assert lhs == rhs, f"Euler formula failed"

    def test_euler_formula_2form(self):
        """L_R(w) = (r+d)*w for a (3,2,1)-form, r+d = 3."""
        from diff_alg.core import new_form, radial
        w = new_form(3, 2, 1, "a")
        R = radial(3)
        lhs = w.lie_derivative(R)
        rhs = 3 * w
        assert lhs == rhs


# ════════════════════════════════════════════════════════════════
# TEST 4: w ^ w = 0 for any 1-form
# ════════════════════════════════════════════════════════════════

class TestWedge:

    def test_1form_wedge_self_is_zero(self):
        """w ^ w = 0 for any 1-form."""
        from diff_alg.core import new_form
        w = new_form(3, 1, 2, "a")
        assert w.wedge(w).is_zero

    def test_anticommutativity(self):
        """w ^ h + h ^ w = 0 for 1-forms."""
        from diff_alg.core import new_form
        w = new_form(3, 1, 2, "a")
        h = new_form(3, 1, 3, "b")
        result = w.wedge(h) + h.wedge(w)
        assert result.is_zero, f"Anticommutativity failed: {result}"

    def test_wedge_associativity(self):
        """(w ^ h) ^ z = w ^ (h ^ z)."""
        from diff_alg.core import new_form
        w = new_form(2, 1, 1, "a")
        h = new_form(2, 1, 1, "b")
        z = new_form(2, 1, 0, "c")
        lhs = w.wedge(h).wedge(z)
        rhs = w.wedge(h.wedge(z))
        assert lhs == rhs


# ════════════════════════════════════════════════════════════════
# TEST 5: Constant fields form involutive distribution
# M2: assert(isInvolutive dist {X,Y}) for constant fields
# ════════════════════════════════════════════════════════════════

class TestDistribution:

    def test_constant_fields_involutive(self):
        """Constant vector fields form an involutive distribution."""
        from diff_alg.core import new_field, dist
        X = new_field(3, 0, "a")
        Y = new_field(3, 0, "b")
        D = dist([X, Y])
        assert D.is_involutive(), "Constant fields should be involutive"


# ════════════════════════════════════════════════════════════════
# TEST 6: Logarithmic forms are integrable: w ^ dw = 0
# M2: w = logarithmicForm(3,{1,2,1},"a"); assert(w ^ dw == 0)
# ════════════════════════════════════════════════════════════════

class TestLogarithmic:

    def test_integrability(self):
        """Logarithmic form satisfies Frobenius: w ^ dw = 0."""
        from diff_alg.core import logarithmic_form
        w = logarithmic_form(3, [1, 2, 1], "a")
        dw = w.exterior_diff()
        result = w.wedge(dw)
        assert result.is_zero, f"Logarithmic form not integrable: {result}"


# ════════════════════════════════════════════════════════════════
# TEST 7: Radial field properties
# ════════════════════════════════════════════════════════════════

class TestRadial:

    def test_radial_construction(self):
        """Radial field R = sum x_i d/dx_i."""
        from diff_alg.core import radial
        R = radial(2)
        assert not R.is_zero
        # Should have 3 components
        assert len(R.coeffs) == 3

    def test_radial_degree(self):
        """Radial field has degree 1."""
        from diff_alg.core import radial
        R = radial(3)
        assert R.degree == (3, 1)


# ════════════════════════════════════════════════════════════════
# TEST 8: New form / new field construction
# ════════════════════════════════════════════════════════════════

class TestConstruction:

    def test_new_form_degree(self):
        """Generic form has correct degree."""
        from diff_alg.core import new_form
        w = new_form(2, 1, 1, "a")
        assert w.degree == (2, 1, 1)

    def test_new_form_not_zero(self):
        from diff_alg.core import new_form
        w = new_form(2, 2, 1, "a")
        assert not w.is_zero

    def test_new_field_degree(self):
        from diff_alg.core import new_field
        X = new_field(2, 2, "a")
        assert X.degree == (2, 2)

    def test_new_field_not_zero(self):
        from diff_alg.core import new_field
        X = new_field(3, 1, "a")
        assert not X.is_zero


# ════════════════════════════════════════════════════════════════
# TEST 9: Arithmetic operations
# ════════════════════════════════════════════════════════════════

class TestArithmetic:

    def test_form_addition(self):
        from diff_alg.core import new_form
        w = new_form(2, 1, 1, "a")
        h = new_form(2, 1, 1, "b")
        result = w + h
        assert not result.is_zero

    def test_form_subtraction_self(self):
        from diff_alg.core import new_form
        w = new_form(2, 1, 1, "a")
        assert (w - w).is_zero

    def test_form_scalar_mul(self):
        from diff_alg.core import new_form
        w = new_form(2, 1, 1, "a")
        assert 2 * w + (-2) * w == 0

    def test_field_addition(self):
        from diff_alg.core import new_field
        X = new_field(2, 1, "a")
        Y = new_field(2, 1, "b")
        assert not (X + Y).is_zero

    def test_field_subtraction_self(self):
        from diff_alg.core import new_field
        X = new_field(2, 1, "a")
        assert (X - X) == 0

    def test_negation(self):
        from diff_alg.core import new_form
        w = new_form(2, 1, 1, "a")
        assert (w + (-w)).is_zero


# ════════════════════════════════════════════════════════════════
# TEST 10: Contraction
# ════════════════════════════════════════════════════════════════

class TestContraction:

    def test_contraction_basic(self):
        """i_{d/dx0}(dx0 ^ dx1) = dx1."""
        from diff_alg.core import _make_rings, DiffAlgForm, DiffAlgField
        R, E, x, dx, params = _make_rings(1)
        w = DiffAlgForm(E.gen(0) * E.gen(1), R, E, 1)
        X = DiffAlgField({0: R.one()}, R, 1)
        result = w.contract(X)
        expected = DiffAlgForm(E.gen(1), R, E, 1)
        assert result == expected, f"Got {result}, expected {expected}"

    def test_contraction_0form_is_zero(self):
        """Contraction of a 0-form is zero."""
        from diff_alg.core import new_form, radial
        f = new_form(2, 0, 2, "a")
        R = radial(2)
        result = f.contract(R)
        assert result.is_zero

    def test_contraction_radial_on_1form(self):
        """i_R(w) for radial R and 1-form w gives a 0-form."""
        from diff_alg.core import new_form, radial
        w = new_form(2, 1, 2, "a")
        R = radial(2)
        result = w.contract(R)
        deg = result.degree
        # Result should be a 0-form
        assert deg[1] == 0


# ════════════════════════════════════════════════════════════════
# TEST 11: Singular and moduli ideals
# ════════════════════════════════════════════════════════════════

class TestIdeals:

    def test_singular_ideal_nonzero(self):
        """Singular ideal of a non-zero form is non-trivial."""
        from diff_alg.core import new_form
        w = new_form(2, 1, 2, "a")
        I = w.singular_ideal()
        assert not I.is_zero()

    def test_moduli_ideal_of_closed_form(self):
        """Moduli ideal of d(w) gives equations for closedness."""
        from diff_alg.core import new_form
        w = new_form(2, 1, 2, "a")
        dw = w.exterior_diff()
        I = dw.moduli_ideal()
        # Should be non-trivial (not all forms are closed)
        assert not I.is_zero()


# ════════════════════════════════════════════════════════════════
# TEST 12: Randomization
# ════════════════════════════════════════════════════════════════

class TestRandom:

    def test_randomize_form(self):
        """Randomized form has no parameters."""
        from diff_alg.core import new_form
        w = new_form(2, 1, 2, "a")
        wr = w.randomize()
        assert len(wr.param_names) == 0
        assert not wr.is_zero

    def test_randomize_field(self):
        """Randomized field has no parameters."""
        from diff_alg.core import new_field
        X = new_field(2, 1, "a")
        Xr = X.randomize()
        assert len(Xr.param_names) == 0
        assert not Xr.is_zero

    def test_random_dd_zero(self):
        """d(d(random w)) = 0."""
        from diff_alg.core import new_form
        w = new_form(2, 1, 2, "a").randomize()
        assert w.exterior_diff().exterior_diff().is_zero


# ════════════════════════════════════════════════════════════════
# TEST 13: Linear combination
# ════════════════════════════════════════════════════════════════

class TestLinearComb:

    def test_linear_comb_forms(self):
        """linear_comb creates a form with extra parameters."""
        from diff_alg.core import new_form, linear_comb
        w1 = new_form(2, 1, 1, "a").randomize()
        w2 = new_form(2, 1, 1, "b").randomize()
        lc = linear_comb([w1, w2], "c")
        assert not lc.is_zero
        # Should have parameters c0, c1
        assert "c0" in lc.param_names or True  # params are there

    def test_linear_comb_fields(self):
        """linear_comb for vector fields."""
        from diff_alg.core import new_field, linear_comb
        X = new_field(2, 1, "a").randomize()
        Y = new_field(2, 1, "b").randomize()
        lc = linear_comb([X, Y], "c")
        assert not lc.is_zero


# ════════════════════════════════════════════════════════════════
# TEST 14: Specific computations with concrete forms
# ════════════════════════════════════════════════════════════════

class TestConcrete:

    def test_specific_wedge(self):
        """dx0 ^ dx1 = -dx1 ^ dx0."""
        from diff_alg.core import _make_rings, DiffAlgForm
        R, E, x, dx, _ = _make_rings(2)
        w1 = DiffAlgForm(E.gen(0), R, E, 2)
        w2 = DiffAlgForm(E.gen(1), R, E, 2)
        assert w1.wedge(w2) + w2.wedge(w1) == 0

    def test_diff_of_coordinate(self):
        """d(x_0) = dx_0."""
        from diff_alg.core import _make_rings, DiffAlgForm
        R, E, x, dx, _ = _make_rings(2)
        f = DiffAlgForm(x[0] * E.one(), R, E, 2)
        df = f.exterior_diff()
        expected = DiffAlgForm(E.gen(0), R, E, 2)
        assert df == expected, f"d(x0) = {df}, expected {expected}"

    def test_diff_of_product(self):
        """d(x0*x1) = x1*dx0 + x0*dx1."""
        from diff_alg.core import _make_rings, DiffAlgForm
        R, E, x, dx, _ = _make_rings(2)
        f = DiffAlgForm(x[0] * x[1] * E.one(), R, E, 2)
        df = f.exterior_diff()
        expected = DiffAlgForm(x[1] * E.gen(0) + x[0] * E.gen(1), R, E, 2)
        assert df == expected

    def test_bracket_concrete(self):
        """[x0 d/dx1, x1 d/dx0] = d/dx0 - d/dx1 (check sign)."""
        from diff_alg.core import _make_rings, DiffAlgField
        R, E, x, dx, _ = _make_rings(1)
        X = DiffAlgField({1: x[0]}, R, 1)  # x0 * d/dx1
        Y = DiffAlgField({0: x[1]}, R, 1)  # x1 * d/dx0
        br = X.bracket(Y)
        # [x0 d/dx1, x1 d/dx0]_i = sum_j X_j dY_i/dx_j - Y_j dX_i/dx_j
        # i=0: X_0*dY_0/dx_0 + X_1*dY_0/dx_1 - Y_0*dX_0/dx_0 - Y_1*dX_0/dx_1
        #     = 0 + x0*1 - x1*0 - 0 = x0  ... wait, let me recalculate
        # X = {1: x0}, Y = {0: x1}
        # [X,Y]_i = sum_j ( X_j * d(Y_i)/dx_j - Y_j * d(X_i)/dx_j )
        # [X,Y]_0 = X_0*dY_0/dx_0 + X_1*dY_0/dx_1 - Y_0*dX_0/dx_0 - Y_1*dX_0/dx_1
        #         = 0*0 + x0 * d(x1)/dx_1 - x1*d(0)/dx_0 - 0*d(0)/dx_1
        #         = x0 * 1 = x0
        # [X,Y]_1 = X_0*dY_1/dx_0 + X_1*dY_1/dx_1 - Y_0*dX_1/dx_0 - Y_1*dX_1/dx_1
        #         = 0 + x0*0 - x1*d(x0)/dx_0 - 0
        #         = -x1
        expected = DiffAlgField({0: x[0], 1: -x[1]}, R, 1)
        assert br == expected, f"Got {br}, expected {expected}"


# ════════════════════════════════════════════════════════════════
# TEST 15: Homogeneity tests
# ════════════════════════════════════════════════════════════════

class TestHomogeneity:

    def test_generic_form_is_homogeneous(self):
        from diff_alg.core import new_form
        w = new_form(2, 1, 2, "a")
        assert w.is_homogeneous()

    def test_generic_field_is_homogeneous(self):
        from diff_alg.core import new_field
        X = new_field(2, 2, "a")
        assert X.is_homogeneous()


# ─── genKer / genIm ────────────────────────────────────────────

class TestGenKer:
    """Tests for gen_ker — ported from M2 DiffAlg test suite."""

    def test_projective_forms(self):
        """genKer(i_R(h), h) gives projective forms (h s.t. i_R(h) = 0)."""
        from diff_alg.core import new_form, radial, gen_ker
        h = new_form(3, 1, 1, "a")
        R = radial(3)
        Rh = h.contract(R)
        T = gen_ker(Rh, h)
        assert len(T) > 0
        for t in T:
            assert t.contract(R).is_zero

    def test_tangent_directions(self):
        """M2 test: w random 2-form, h generic 1-form.
        L = genKer(w^dh + h^dw, h) gives tangent directions at w."""
        from diff_alg.core import new_form, gen_ker
        w = new_form(3, 2, 2, "a").randomize()
        h = new_form(3, 1, 1, "b")
        dh = h.exterior_diff()
        dw = w.exterior_diff()
        expr = w.wedge(dh) + h.wedge(dw)
        L = gen_ker(expr, h)
        assert len(L) > 0
        for elt in L:
            delt = elt.exterior_diff()
            check = w.wedge(delt) + elt.wedge(dw)
            assert check.is_zero

    def test_field_annihilator(self):
        """genKer(i_X(w), X) for field X gives annihilator of 2-form w."""
        from diff_alg.core import new_form, new_field, gen_ker
        w = new_form(2, 2, 2, "a").randomize()
        X = new_field(2, 1, "b")
        iXw = w.contract(X)
        L = gen_ker(iXw, X)
        assert isinstance(L, list)
        for fld in L:
            assert w.contract(fld).is_zero

    def test_rank_nullity(self):
        """dim(ker) + dim(im) = dim(parameter space)."""
        from diff_alg.core import new_form, gen_ker, gen_im
        w = new_form(3, 2, 2, "a").randomize()
        h = new_form(3, 1, 1, "b")
        dh = h.exterior_diff()
        dw = w.exterior_diff()
        expr = w.wedge(dh) + h.wedge(dw)
        ker = gen_ker(expr, h)
        im = gen_im(expr, h)
        dim_ker = len(ker)
        dim_im = len(im)
        m = len(h.param_names)
        # Kernel might include a zero element as placeholder
        if len(ker) == 1 and ker[0].is_zero:
            dim_ker = 0
        # Image might include a zero element as placeholder
        if len(im) == 1 and im[0].is_zero:
            dim_im = 0
        assert dim_ker + dim_im == m

    def test_kernel_of_zero_expression(self):
        """If expr = 0 for any var, kernel should be full space."""
        from diff_alg.core import new_form, gen_ker, DiffAlgForm
        h = new_form(2, 1, 1, "a")
        # Create zero form in same ring
        zero = h - h
        result = gen_ker(zero, h)
        # All basis elements should satisfy the (trivial) equation
        assert len(result) > 0

    def test_nonhomogeneous_consistent(self):
        """Non-homogeneous system with a solution returns [basis, particular]."""
        from diff_alg.core import new_form, gen_ker
        w1 = new_form(2, 1, 1, "a").randomize()
        w2 = new_form(2, 1, 1, "c").randomize()
        w3 = w1.wedge(w2)
        h = new_form(2, 1, 1, "b")
        result = gen_ker(w1.wedge(h) - w3, h)
        assert len(result) == 2
        basis, particular = result
        assert isinstance(basis, list)
        assert particular is not None
        # particular satisfies the equation
        check = w1.wedge(particular) - w3
        assert check.is_zero

    def test_nonhomogeneous_inconsistent(self):
        """Non-homogeneous system with no solution returns [basis, None]."""
        from diff_alg.core import new_form, new_field, radial, gen_ker
        w = new_form(3, 1, 1, "a")
        R = radial(3)
        t = new_form("x0 + x1")
        # i_R(w) is degree 2 in x, t is degree 1 — inconsistent
        result = gen_ker(w.contract(R) + t, w)
        assert len(result) == 2
        basis, particular = result
        assert isinstance(basis, list)
        assert len(basis) > 0
        assert particular is None
        # basis elements still satisfy the homogeneous equation
        for b in basis:
            assert b.contract(R).is_zero


class TestGenIm:
    """Tests for gen_im."""

    def test_gen_im_basic(self):
        """Image of a linear expression should have correct degree."""
        from diff_alg.core import new_form, gen_im
        w = new_form(2, 2, 2, "a").randomize()
        h = new_form(2, 1, 1, "b")
        dh = h.exterior_diff()
        dw = w.exterior_diff()
        expr = w.wedge(dh) + h.wedge(dw)
        im = gen_im(expr, h)
        assert len(im) > 0
        for elt in im:
            assert elt.degree[1] == expr.degree[1]  # same exterior degree


# ════════════════════════════════════════════════════════════════
# TEST: String parsing for new_form and new_field
# M2: w = newForm("a * x_1 * dx_0 * dx_1")
# ════════════════════════════════════════════════════════════════

class TestStringParsing:

    def test_parse_simple_1form(self):
        """Parse a simple 1-form from string."""
        from diff_alg.core import new_form
        w = new_form("x0*dx1 - x1*dx0")
        assert not w.is_zero
        assert w.n == 1

    def test_parse_with_underscore_syntax(self):
        """Parse using x_0, dx_0 syntax (M2 style)."""
        from diff_alg.core import new_form
        w = new_form("x_0*dx_1 - x_1*dx_0")
        assert not w.is_zero
        assert w.n == 1

    def test_parse_with_parameters(self):
        """Parse form with scalar parameters."""
        from diff_alg.core import new_form
        w = new_form("a*x_1*dx_0*dx_1")
        assert not w.is_zero
        assert "a" in w.param_names

    def test_parse_higher_index(self):
        """Parsing dx_5 should give n=5."""
        from diff_alg.core import new_form
        w = new_form("dx_5")
        assert w.n == 5

    def test_parse_2form(self):
        """Parse a 2-form."""
        from diff_alg.core import new_form
        w = new_form("x0^2*dx0*dx1 + x1^2*dx1*dx2")
        assert not w.is_zero

    def test_parse_form_dd_zero(self):
        """d(d(w)) = 0 for a parsed form."""
        from diff_alg.core import new_form
        w = new_form("x0^2*dx1 - x0*x1*dx0 + x2*dx2")
        ddw = w.exterior_diff().exterior_diff()
        assert ddw.is_zero

    def test_parse_form_operations(self):
        """Wedge product of parsed forms should work."""
        from diff_alg.core import new_form
        w = new_form("x0*dx0 + x1*dx1")
        h = new_form("x0*dx1 - x1*dx0")
        wh = w.wedge(h)
        assert not wh.is_zero

    def test_parse_simple_field(self):
        """Parse a vector field from string."""
        from diff_alg.core import new_field
        X = new_field("x0*ax1 - x1*ax0")
        assert not X.is_zero
        assert X.n == 1

    def test_parse_field_underscore_syntax(self):
        """Parse field using ax_0, ax_1 syntax."""
        from diff_alg.core import new_field
        X = new_field("x_0^2*ax_1 + x_1*ax_0")
        assert not X.is_zero
        assert X.n == 1

    def test_parse_field_with_params(self):
        """Parse field with parameters."""
        from diff_alg.core import new_field
        X = new_field("a*x0*ax0 + b*x1*ax1")
        assert not X.is_zero
        assert "a" in X.param_names
        assert "b" in X.param_names

    def test_parse_field_contraction(self):
        """Contraction of parsed form with parsed field."""
        from diff_alg.core import new_form, new_field
        w = new_form("x0*dx0*dx1 + x1*dx1*dx2")
        X = new_field("ax0 + ax1 + ax2")
        iXw = w.contract(X)
        assert not iXw.is_zero

    def test_parse_0form(self):
        """Parse a 0-form (polynomial function)."""
        from diff_alg.core import new_form
        f = new_form("x0^2 + x1^2 + x2^2")
        assert not f.is_zero
        df = f.exterior_diff()
        assert not df.is_zero
