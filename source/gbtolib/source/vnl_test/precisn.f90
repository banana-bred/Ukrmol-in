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
module precisn
! Definition of precision-related parameters
! from Martin Plummer, Nov 2009
  implicit none
  public

! Uncomment if needed
!  logical, parameter    :: real64 = .TRUE.
!  integer, parameter    :: ibyte = 4
!  integer, parameter    :: rbyte = 8
!  integer, parameter    :: zbyte = 16

!  Set kind type parameters.
!  The argument to the function `selected_real_kind' is the 
!  requested minimum number of decimal digits of accuracy.
  integer, parameter    :: sp=selected_real_kind(6)   !> `single' precision 
  integer, parameter    :: wp=selected_real_kind(12)  !> `double' precision
  integer, parameter    :: ep1=selected_real_kind(19) !> `quad', extended, precision
  integer, parameter    :: ep = MAX(ep1, wp) !> extended precision if possible; double precision if not

! Uncomment if needed
!  real(wp), parameter   :: acc8 = 2.0e-16_wp
!  real(wp), parameter   :: acc16 = 3.0e-33_ep
!  real(wp), parameter   :: fpmax = 1.0e60_wp
!  real(wp), parameter   :: fpmin = 1.0e-60_wp
!  real(wp), parameter   :: fplmax = 140.0_wp
!  real(wp), parameter   :: fplmin = -140.0_wp

end module precisn
