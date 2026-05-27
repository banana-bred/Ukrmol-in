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
program t
  use precisn
  use cgto_pw_expansions_mod
  use basis_data_generic_mod
  use bspline_grid_mod

  type(CGTO_shell_data_obj) :: cgto_shell_A, cgto_shell_B
  type(CGTO_shell_pair_pw_expansion_obj) :: CGTO_pair_pw

  real(kind=cfp) :: a = 10.0_cfp
  integer :: first_bspline_index,max_bspline_l,max_prop_l,n_rng_knot,max_l_legendre
  type(bspline_grid_obj) :: bspline_grid

     call cgto_shell_A%make_space(2)
     cgto_shell_A%l = 2
     cgto_shell_A%center(1:3) = -0.50_cfp !(/0.0_cfp,0.0_cfp,-0.7005214769179526_cfp/)
     cgto_shell_A%number_of_primitives = 2
     cgto_shell_A%exponents(1:cgto_shell_A%number_of_primitives) = (/2.5_cfp,0.75_cfp/)
     cgto_shell_A%contractions(1:cgto_shell_A%number_of_primitives) = (/1.0_cfp,0.25_cfp/)
     cgto_shell_A%non_zero_at_boundary = .false.
     cgto_shell_A%number_of_functions = 2*cgto_shell_A%l+1
     call cgto_shell_A%normalize
  
     call cgto_shell_B%make_space(3)
     cgto_shell_B%l = 1 !0
     cgto_shell_B%center(1:3) = 5.0_cfp !(/0.0_cfp,0.0_cfp,0.7005214769179526_cfp/)
     cgto_shell_B%number_of_primitives = 3
     cgto_shell_B%exponents(1:cgto_shell_B%number_of_primitives) = (/1.5_cfp,0.15_cfp,3.0_cfp/)
     cgto_shell_B%contractions(1:cgto_shell_B%number_of_primitives) = (/1.10_cfp,2.0_cfp,1.5_cfp/)
     cgto_shell_B%non_zero_at_boundary = .false.
     cgto_shell_B%number_of_functions = 2*cgto_shell_B%l+1
     call cgto_shell_B%normalize

     first_bspline_index = 2
     max_bspline_l = 0
     max_prop_l = 2
     n_rng_knot = 1
     max_l_legendre = 40
     call bspline_grid%init_grid(2.0_cfp,a,10,12)
     call init_CGTO_pw_expansions_mod(max_l_legendre,max(cgto_shell_B%l,cgto_shell_A%l))

     call CGTO_pair_pw%init_CGTO_shell_pair_pw_expansion(cgto_shell_A,1,cgto_shell_B,2)
     call CGTO_pair_pw%eval_regular_grid(0.0_cfp,a,1.0_cfp)
     call CGTO_pair_pw%eval_CGTO_shell_pair_pw_expansion

end program t
