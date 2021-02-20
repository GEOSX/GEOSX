.. _PoroElastic:

############################################
Porous solids
############################################

Linear isotropic poroelasticity
=========================

Overview
----------------------------------------

This model may be used to represents a porous material with a linear poroelastic isotropic response to coupled deformation-diffusion processes.
The relationship between stress and strain is typically formulated within the framework of the `Biot theory of poroelasticity <https://doi.org/10.1016/B978-0-08-040615-2.50011-3>`__,
which for the case of linear poroelasticity, may be expressed as:

.. math::
   \sigma_{ij} = \sigma\prime_{ij}  - b p \delta_{ij} = \lambda \epsilon_{kk} + 2 \mu \epsilon_{ij} - b p \delta_{ij},

where :math:`\sigma_{ij}` is the :math:`ij` component of the total stress tensor,
:math:`\sigma\prime_{ij}` is the :math:`ij` component of the effective (Cauchy) stress tensor,
:math:`\epsilon_{ij}` is the :math:`ij` component of the the strain tensor,
:math:`\lambda` is the Lames elastic constant,
:math:`\mu` is the elastic shear modulus,
:math:`b` is Biot's coefficient,
:math:`p` is fluid pressure,
and :math:`\delta` is Kronecker delta.

Hooke's Law may also be expressed using `Voigt notation <https://en.wikipedia.org/wiki/Voigt_notation>`__ for the effective stress and strain tensors as:

.. math::
   \tensor{\sigma\prime} = \tensor{C} \cdot \tensor{\epsilon},

or,

.. math::
    \begin{bmatrix}
      \sigma\prime_{11} \\
      \sigma\prime_{22} \\
      \sigma\prime_{33} \\
      \sigma\prime_{23} \\
      \sigma\prime_{13} \\
      \sigma\prime_{12}
    \end{bmatrix}
    =
    \begin{bmatrix}
      2\mu+\lambda  &   \lambda     &   \lambda   & 0   & 0 & 0 \\
          \lambda     &  2\mu+\lambda   &   \lambda   & 0   & 0 & 0 \\
          \lambda     &    \lambda    & 2\mu+\lambda & 0  & 0 & 0 \\
          0         &       0     &       0 &\mu  & 0 & 0 \\
      0         &           0     & 0       & 0   & \mu & 0 \\
      0         &       0     & 0       & 0   & 0 & \mu
    \end{bmatrix}
    \begin{bmatrix}
      \epsilon_{11} \\
      \epsilon_{22} \\
      \epsilon_{33} \\
      2\epsilon_{23} \\
      2\epsilon_{13} \\
      2\epsilon_{12}
    \end{bmatrix}.

Variations
-----------------------
The application of linear poroelasticity as presented above is typically restricted to the case of
`infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__.
For the case of `infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__, the
above relations are applied directly.
For the case of `finite strain <https://en.wikipedia.org/wiki/Finite_strain_theory>`__,
the above relations may be slightly modified to an incremental update and rotation:

.. math::
   \Delta \tensor{\sigma\prime} &= \tensor{C} \cdot \hat{\tensor{D}},\\
   \tensor{\sigma\prime}^{n+1} &= \hat{\tensor{R}}( \tensor{\sigma\prime}^{n} + \Delta \tensor{\sigma\prime} ) \hat{\tensor{R}}^T,

where :math:`\hat{\tensor{D}}` is the "incremental rate of deformation tensor" and :math:`\hat{\tensor{R}}` is the incremental rotation tensor, which are
typically calculated from the `velocity gradient <https://en.wikipedia.org/wiki/Strain-rate_tensor>`__.
This extension into finite strain constitutes a `hypo-elastic update <https://en.wikipedia.org/wiki/Hypoelastic_material>`__,
and the choice of method to calculate :math:`\hat{\tensor{D}}`, and :math:`\hat{\tensor{R}}` determines if the update is objective.
One commonly used method is the `Hughes-Winget <https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620151210>`__ algorithm.


Parameters
-----------------------

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/PoroLinearElasticIsotropic.rst

Example
-----------------------

.. code-block:: xml

  <Constitutive>
    <PoroLinearElasticIsotropic name="shale"
                                defaultDensity="2700"
                                defaultBulkModulus="61.9e6"
                                defaultShearModulus="28.57e6"
                                BiotCoefficient="1.0"/>
  </Constitutive>
  
Linear anisotropic poroelasticity
=============================

  Overview
----------------------------------

This model may be used to represents a porous material with a linear poroelastic anisotropic response to coupled deformation-diffusion processes.
The relationship between stress and strain is typically formulated within the framework of the `Biot theory of poroelasticity <https://doi.org/10.1016/B978-0-08-040615-2.50011-3>`__,
which for the case of linear poroelasticity, may be expressed as:

.. math::
   \sigma_{ij} = \sigma\prime_{ij}  - b p \delta_{ij} = C_{ijkl} \epsilon_{kl}  - b p \delta_{ij},

where :math:`\sigma_{ij}` is the :math:`ij` component of the total stress tensor,
:math:`\sigma\prime_{ij}` is the :math:`ij` component of the effective (Cauchy) stress tensor,
:math:`\epsilon_{ij}` is the :math:`ij` component of the the strain tensor,
:math:`\lambda` is the Lames elastic constant,
:math:`\mu` is the elastic shear modulus,
:math:`b` is Biot's coefficient,
:math:`p` is fluid pressure,
and :math:`\delta` is Kronecker delta.

Hooke's Law may also be expressed using `Voigt notation <https://en.wikipedia.org/wiki/Voigt_notation>`__ for the effective stress and strain tensors as:

.. math::
   \tensor{\sigma\prime} = \tensor{C} \cdot \tensor{\epsilon},

or,

.. math::
    \begin{bmatrix}
      \sigma\prime_{11} \\
      \sigma\prime_{22} \\
      \sigma\prime_{33} \\
      \sigma\prime_{23} \\
      \sigma\prime_{13} \\
      \sigma\prime_{12}
    \end{bmatrix}
    =
    \begin{bmatrix}
      c_{1111} & c_{1122} & c_{1133} & c_{1123} & c_{1113} & c_{1112} \\
      c_{2211} & c_{2222} & c_{2233} & c_{2223} & c_{2213} & c_{2212} \\
      c_{3311} & c_{3322} & c_{3333} & c_{3323} & c_{3313} & c_{3312} \\
      c_{2311} & c_{2322} & c_{2333} & c_{2323} & c_{2313} & c_{2312} \\
      c_{1311} & c_{1322} & c_{1333} & c_{1323} & c_{1313} & c_{1312} \\
      c_{2311} & c_{2322} & c_{2333} & c_{2323} & c_{2313} & c_{2312}
    \end{bmatrix}
    \begin{bmatrix}
      \epsilon_{11} \\
      \epsilon_{22} \\
      \epsilon_{33} \\
      2\epsilon_{23} \\
      2\epsilon_{13} \\
      2\epsilon_{12}
    \end{bmatrix}.


Variations
-----------------------
The application of linear poroelasticity as presented above is typically restricted to the case of
`infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__.
For the case of `infinitesimal strain <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__, the
above relations are applied directly.
For the case of `finite strain <https://en.wikipedia.org/wiki/Finite_strain_theory>`__,
the above relations may be slightly modified to an incremental update and rotation:

.. math::
   \Delta \tensor{\sigma\prime} &= \tensor{C} \cdot \hat{\tensor{D}},\\
   \tensor{\sigma\prime}^{n+1} &= \hat{\tensor{R}}( \tensor{\sigma\prime}^{n} + \Delta \tensor{\sigma\prime} ) \hat{\tensor{R}}^T,

where :math:`\hat{\tensor{D}}` is the "incremental rate of deformation tensor" and :math:`\hat{\tensor{R}}` is the incremental rotation tensor, which are
typically calculated from the `velocity gradient <https://en.wikipedia.org/wiki/Strain-rate_tensor>`__.
This extension into finite strain constitutes a `hypo-elastic update <https://en.wikipedia.org/wiki/Hypoelastic_material>`__,
and the choice of method to calculate :math:`\hat{\tensor{D}}`, and :math:`\hat{\tensor{R}}` determines if the update is objective.
One commonly used method is the `Hughes-Winget <https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620151210>`__ algorithm.


Parameters
-----------------------

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/PoroLinearElasticAnisotropic.rst

Example
-----------------------

.. code-block:: xml

  <Constitutive>
    <PoroLinearElasticAnisotropic name="shale"
                                  defaultDensity="2700"
                                  c11=1.0e10 c12=1.1e9  c13=1.2e9  c14=1.3e9 c15=1.4e9 c16=1.5e9
                                  c21=2.0e9  c22=2.1e10 c23=2.2e9  c24=2.3e9 c25=2.4e9 c26=2.5e9
                                  c31=3.0e9  c32=3.1e9  c33=3.2e10 c34=3.3e9 c35=3.4e9 c36=3.5e9
                                  c41=4.0e9  c42=4.1e9  c43=4.2e9  c44=4.3e9 c45=4.4e9 c46=4.5e9
                                  c51=5.0e9  c52=5.1e9  c53=5.2e9  c54=5.3e9 c55=5.4e9 c56=5.5e9
                                  c61=6.0e9  c62=6.1e9  c63=6.2e9  c64=6.3e9 c65=6.4e9 c66=6.5e9
                                  BiotCoefficient="1.0"/>
  </Constitutive>

