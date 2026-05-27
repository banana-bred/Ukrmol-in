! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
module const
   use precisn
      
   !> unit for standard output
   integer, parameter :: stdout = 6
   !> unit for standard input
   integer, parameter :: stdin  = 5

   !> length of one line
   integer, parameter :: line_len = 132

   integer, parameter :: len_ufmat = 11
   character(len_ufmat-2), parameter :: fmat = 'formatted'
   character(len_ufmat), parameter :: ufmat = 'unformatted'

   !> Default character string used for the int_input_output%header variable.
   character(len=19), parameter :: no_header = 'No header specified'

   !> Character parameters used by the method l_mol_basis%molecular_integrals to identify the type of molecular integral requested. These headers are used by the user to request the specific integral.
   character(len=*), parameter :: cgto_overlap_kinetic_ints = 'Overlap and kinetic energy integrals for contracted spherical GTOs'
   character(len=*), parameter :: cgto_cgto_two_el_ints = 'Two-electron integrals for contracted spherical GTOs'
   character(len=*), parameter :: bto_overlap_kinetic_ints = 'Overlap and kinetic energy integrals for BTOs'
   character(len=*), parameter :: cgto_bto_overlap_kinetic_ints = 'Overlap and kinetic energy integrals between contracted spherical GTOs and BTOs'
   character(len=*), parameter :: cgto_cgto_bloch_terms = 'Matrix elements of the Bloch operator between contracted spherical CMS GTOs'
   character(len=*), parameter :: cgto_bto_bloch_terms = 'Matrix elements of the Bloch operator between BTOs and contracted spherical CMS GTOs'
   character(len=*), parameter :: bto_bto_bloch_terms = 'Matrix elements of the Bloch operator between BTOs'

   !> translation_formula_obj: maximum principal quantum number for which we have the translation formulas. This is constant is also used as a maximum GTO L value in all integral calculations.
   !> \warning This parameter is related to the parameters of the Boys function (see boys_function). Therefore any change in this parameter must be followed by thorough tests of accuracy of evaluation of
   !> the Boys function and changes to the related parameters (also see below). 
   integer, parameter :: tf_max_l = 6
   !> translation_formula_obj: maximum azimuthal quantum number (corresponding to tf_max_l) for which we have the translation formulas.
   integer, parameter :: tf_max_m = 6
   !> translation_formula_obj: number of terms in the longest translation formula from the ones that are included in translation_formula_obj
   integer, parameter :: tf_max_n = 36
   
   !> Convergence parameter for the calculation of the Boys function using boys_function. The convergence is chosen as the machine epsilon for double precision reals.
   real(kind=wp), parameter :: boys_tol_dprec = epsilon(1.0_wp)

   !> Empirically found value for the largest power of T required to calculate the Boys function \f$F_{m}(T)\f$ for \f$T=0,\dots,60\f$ with the accuracy epsilon(1.0_wp).
   integer, parameter :: boys_maxi_dprec = 140

   !> Convergence parameter for the calculation of the Boys function using boys_function_quad. The convergence is chosen as the machine epsilon for quad precision reals.
   real(kind=wp), parameter :: boys_tol_qprec = epsilon(1.0_ep1)

   !> The Boys function calculated using the routine boys_function will be calculated using the asymptotic expansion if the input parameter T is .ge. this value. Note that this value
   !> was determined using Mathematica assuming tf_max_l = 6, i.e. mmax = 24 (see boys_function for details). For T=85 and mmax=24 the error in using the asymptotic expansion is not greater than ~3.10^-15.
   !> For T=60 the values of the Boys function for m .ge. 12 are .le. 10^{-16}. For m = 12 the relative precision of the asymptotic formula is ~10^{-14}. Therefore for T > 60 and double precision we can 
   !> safely calculate the Boys function using the asymptotic expansion. It is also important to note that for T=60 the full numerical evaluation involving series fails in double precision for large mmax.
   !> Therefore increasing boys_f_dprec_asym_thr beyond 60 cannot go without modifying boys_function to branch into boys_function_quad for T>60 .and. T<boys_f_dprec_asym_thr in order to avoid numerical 
   !> errors in this range of T. In general the larger the tf_max_l the larger the value of boys_f_dprec_asym_thr. Clearly, the numerical stability of boys_function effectively determines the limit on the
   !> largest GTO L in the basis: if large L are needed then boys_function_quad must be used (for certain range of T values) which is also much slower than boys_function.
   real(kind=wp), parameter :: boys_f_dprec_asym_thr = 60.0_wp

   !> See comment for boys_f_dprec_asym_thr. This value is the threshold for the use of the asymptotic expansion in the quadruple precision function boys_function_quad (assuming mmax=24). For T=140 the 
   !> asymptotic formula gives results accurate to full quadruple precision.
   real(kind=wp), parameter :: boys_f_qprec_asym_thr = 140.0_wp

   !> Molecular integrals with relative precision last than or equal to this value will trigger an error. This parameter is used only as a default in integral_options%prec, so this
   !> parameter can be effectively adjusted on run-time. 
   !> \warning Not all integral calculation routines are necessarily using this parameter.
   real(kind=wp), parameter :: int_rel_prec = 10d-10

   !> Molecular integrals (contracted) smaller than this value will be neglected. Similarly to int_rel_prec this value is used only as a sensible default in integral_options%tol.
   real(kind=wp), parameter :: int_del_thr = 10d-11

   !> Maximum allowed value for a cross overlap of two orbitals after the symmetric orthogonalization has been performed.
   !> \warning This value should be actually equal to the threshold value for the integrals.
   real(kind=wp), parameter :: thrs_symm_ortho = 10d-6

   !> Maximum allowed value for a cross overlap of two orbitals after the Gramm-Schmidt orthogonalization has been performed.
   real(kind=wp), parameter :: thrs_gs_ortho = 10d-10

   !> Coefficients in the transformation matrix for the symmetric orthogonalization smaller than this value will be neglected.
   real(kind=wp), parameter :: thrs_cf_sym_ortho_trans = 10d-10

   !> Threshold value for self-overlap of an orbital (before normalization) using during the Gramm-Schmidt orthogonalization. 
   !> If an orbital has a self-overlap (before normalization) smaller than this value then we assume that linear dependency in the orbital basis is present. In other words we decide that self-overlaps
   !> smaller than this value would lead to numerical problems.
   real(kind=wp), parameter :: thrs_lin_dep_gs_ortho = 10d-7

   !> Absolute precision for the numerical quadrature routine dqags.
   real(kind=wp), parameter :: EPSABS = 10d-10
   !> Relative precision for the numerical quadrature routine dqags.
   real(kind=wp), parameter :: EPSREL = 10d-10
   !> Determines the maximum number of subintervals in the partition of the given integration interval in dqags.
   integer, parameter :: LIMIT = 1000
   !> Dimensioning parameter for dqags.
   integer, parameter :: LENW = 4*LIMIT

   !> Numerical identifier of the C1 point group-symmetry
   integer, parameter :: C1_id = 1
   !> Numerical identifier of the Cs point group-symmetry
   integer, parameter :: Cs_id = 2
   !> Numerical identifier of the C2 point group-symmetry
   integer, parameter :: C2_id = 3
   !> Numerical identifier of the Ci point group-symmetry
   integer, parameter :: Ci_id = 4
   !> Numerical identifier of the C2v point group-symmetry
   integer, parameter :: C2v_id = 5
   !> Numerical identifier of the C2h point group-symmetry
   integer, parameter :: C2h_id = 6
   !> Numerical identifier of the D2 point group-symmetry
   integer, parameter :: D2_id = 7
   !> Numerical identifier of the D2h point group-symmetry
   integer, parameter :: D2h_id = 8

   !> Length of the character variable 'name' in the nucleus_type object.
   integer, parameter :: nuc_nam_len = 2

   !> Length of the character variable sym_op specifying the symmetry operation in the geometry_obj object.
   integer, parameter :: sym_op_nam_len = 3

end module
