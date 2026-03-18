"""
DiffAlg — Differential Algebra package for SageMath.

Port of the Macaulay2 package DiffAlg by Dubinsky, Massri, Molinuevo, Quallbrunn.
Specialized routines for homogeneous differential forms and vector fields.

Operations:
    - Wedge products and exterior differential of differential forms
    - Contraction and Lie derivative with respect to a vector field
    - Lie brackets between vector fields
    - Kernel and image of degree-one differential operators
    - Logarithmic forms, distributions, involutivity
"""

from .core import (
    DiffAlgForm,
    DiffAlgField,
    DiffAlgDistribution,
    new_form,
    new_field,
    radial,
    linear_comb,
    logarithmic_form,
    dist,
    gen_ker,
    gen_im,
)
