Solid Models
======================

Solid mechanics constitutive models govern the relationship between measures of 
deformation (and its history) and material state (e.g. stress) in a solid 
material.
The measures of deformation and history used as kinematic input to a 
constitutive relation typically depend on the assumptions of
`infinitesimal <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_ or 
`finite <https://en.wikipedia.org/wiki/Finite_strain_theory>`_ 
strain, as well as the rate dependence of the constitutive
relation itself.

Infinitesimal/Small Strain/Deformation
-------------------------------------------------
An incremental constitutive update for stress at some time :math:`t^{n+1}` 
under infinitesimal strain assumptions can be expressed as 

.. math::
   \bm{\sigma^{n+1}} = \bm{\sigma}(\bm{\epsilon^{n}}, \bm{ \dot{{\epsilon}}}, Q^n, \Delta t),
where :math:`\bm{\epsilon}` is the small strain tensor, 
:math:`\bm{\dot{{\epsilon}}}` rate of the small strain tensor, and :math:`Q^n`
is the collection of material state variables (which includes stress), and 
:math:`\Delta t = t^{n+1} - t^n`.
Note that for path and rate independent models, such as linear elasticity,
the constitutive update may also be formulated in terms of total strain:

.. math::
   \bm{\sigma^{n+1}} = \bm{\sigma}(\bm{\epsilon^{n+1}}).


Finite/Large Strain/Deformation
-------------------------------------------------
Similarly, an incremental constitutive update for under finite strain 
assumptions may be expressed as 

.. math::
   \bm{\sigma^{n+1}} = \bm{\sigma}(\mathbf{F^n}, \mathbf{ \dot{{F}}}, Q^n, \Delta t),
where :math:`\mathbf{F}` and :math:`\mathbf{\dot{{F}}}` are the 
finite strain deformation gradient and its rate. 
Note that the finite deformation stress measures such as one of the
Piola-Kirchhoff stress tensors, and one of the finite deformation strain tensors
are commonly used in finite deformation constitutive relations.
Of course, any of these finite deformation tensors can be recovered with 
knowledge of the Cauchy stress (:math:`\bm{\sigma}`) and the deformation 
gradient :math:`\mathbf{F}`.

When modeling finite deformations, a constitutive update may expressed in 
terms of a hypoelastic/plastic update.
A hypoelastic/plastic update is typically a relation similar to the infinitesimal 
strain update:

.. math::
   \bm{\sigma^{n+1}} = \bm{\sigma}(\mathbf{D}, \mathbf{R}, Q^n, \Delta t),

where the rate/increment in the symmetric small strain tensor is replaced by an 
objective incremental rate of deformation tensor (:math:`\hat{\mathbf{D}}`), and an objective 
rotation tensor (:math:`\mathbf{R}`) which is applied to the state variables in :math:`Q^n`.
A commonly used method for attaining a weakly objective kinematics for a hypo-update
is given by the `Hughes-Winget <https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620151210>`__ algorithm:

.. math::
      \hat{\tensor{D}} &= \frac{1}{2}( \hat{\tensor{L}} + \hat{\tensor{L}}^T ), \\
      \hat{\tensor{W}} &= \frac{1}{2}( \hat{\tensor{L}} - \hat{\tensor{L}}^T ), \\
      \hat{\tensor{R}} &= ( \tensor{I} - \frac{1}{2} \hat{\tensor{W}} )^{-1} ( \tensor{I} + \frac{1}{2} \hat{\tensor{W}} ), \\
      \tensor{\bar{\sigma}}^{n+1} &= \hat{\tensor{R}} \tensor{\sigma}^{n} \hat{\tensor{R}}^T, \\
      \tensor{\sigma}^{n+1} &= \tensor{\bar{\sigma}}^{n+1} + \Delta \tensor{\sigma}

where :math:`\hat{\tensor{L}}` is the `velocity gradient <https://en.wikipedia.org/wiki/Strain-rate_tensor>`_ multiplied by :math:`\Delta t`.

In a hyperelastic/plastic update, the elastic portion of the relation is 
expressed in terms of a stored/strain energy function that serves as the
potential for the stress:

.. math::
   \mathbf{S} = \frac{\partial \psi (\tensor{C})}{ \partial \tensor{C} },
where :math:`\mathbf{S}` is the second Piola-Kirchoff stress, :math:`\psi` is 
the stored energy potential, and :math:`\tensor{C}` is the right Cauchy-Green 
deformation tensor.



Voigt Representation of Rank-2 Symmetric Tensors
------------------------------------------------
In GEOSX we express rank-2 symmetric tensors using 
`Voigt notation <https://en.wikipedia.org/wiki/Voigt_notation>`_.
Note that in Voigt notation, the stress tensors are represented by the direct 
tensor values, while deformation measures are expressed in engineering strain 
s.t. the shear components of strain are stored with a factor of 2.
This dual representation, the string energy density is simply the inner product
of the Voigt representation of stress and strain.
Thus , the Voigt representation of the common rank-2 symmetric tensors are
expressed as:

.. math::
   \bm{\sigma} = 
      \begin{bmatrix} \
         \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12}
      \end{bmatrix},
   \mathbf{S} = 
      \begin{bmatrix} 
         S_{11} \\ S_{22} \\ S_{33} \\ S_{23} \\ S_{13} \\ S_{12}
      \end{bmatrix},
   \bm{\epsilon} = 
      \begin{bmatrix} 
         \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\\ 
         2 \epsilon_{23} \\ 2 \epsilon_{13} \\ 2 \epsilon_{12}
      \end{bmatrix}, 
   \mathbf{D} = 
      \begin{bmatrix}
         D_{11} \\ D_{22} \\ D_{33} \\ 2 D_{23} \\ 2 D_{13} \\ 2 D_{12}
      \end{bmatrix},
   \mathbf{E} = 
      \begin{bmatrix}
         E_{11} \\ E_{22} \\   E_{33} \\ 2 E_{23} \\ 2 E_{13} \\ 2 E_{12}
      \end{bmatrix},
   \mathbf{B} = 
      \begin{bmatrix}
         B_{11} \\ B_{22} \\ B_{33} \\ 2 B_{23} \\ 2 B_{13} \\ 2 B_{12}
   \end{bmatrix},
   \mathbf{C} = 
      \begin{bmatrix} 
        C_{11} \\ C_{22} \\ C_{33} \\ 2 C_{23} \\ 2 C_{13} \\ 2 C_{12}
      \end{bmatrix},

where 
:math:`\bm{\sigma}` is the Cauchy stress, 
:math:`\mathbf{S}` is the second Piola-Kirchhoff stress,
:math:`\bm{\epsilon}` is the small strain tensor,
:math:`\mathbf{D}` is the rate of deformation tensor,
:math:`\mathbf{E}` is the Lagrangian finite strain tensor,
:math:`\mathbf{B}` is the left Cauchy-Green tensor,
:math:`\mathbf{C}` is the right Cauchy-Green deformation tensor.
Therefore, when examining symmetric tensorial deformation measures produced from 
GEOSX, please note that the values of the shear components will be 2x the actual
shear components.



.. toctree::
   :maxdepth: 2

   LinearElastic
   
   TwoInvariantPlasticity

   PoroElastic
   
   DamageModel
