! Copyright 2025
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

!> \brief   Unit test for `symmetry_gbl::Xlm_mgvn`
!> \authors J Benda
!> \date    2025
!>
!> Calculates irreducible representations of the first few real spherical harmonics and compares them to reference numbers.
!> These were obtained in Maxima 5.47.0 with the following program:
!>
!> \code{.mac}
!>   declare([x, y, z], real)$
!>   declare([l, m], integer)$
!>
!>   Xlm(l, m, r, x, y, z) :=
!>       if     m < 0 then (-1)^m*sqrt(2)*sqrt((2*l+1)/(4*%pi)*(l+m)!/(l-m)!)*assoc_legendre_p(l,-m,z/r)*imagpart((x/r+%i*y/r)**-m)
!>       elseif m > 0 then (-1)^m*sqrt(2)*sqrt((2*l+1)/(4*%pi)*(l-m)!/(l+m)!)*assoc_legendre_p(l,m,z/r)*realpart((x/r+%i*y/r)**m)
!>       else              sqrt((2*l+1)/(4*%pi))*legendre_p(l,z/r)$
!>
!>   for l: 0 thru 6 do
!>   (
!>       for m: -l thru +l do
!>       (
!>           X0: radcan(Xlm(l, m, r,  x,  y,  z)),
!>           Xx: radcan(Xlm(l, m, r, -x,  y,  z)),
!>           Xy: radcan(Xlm(l, m, r,  x, -y,  z)),
!>           Xz: radcan(Xlm(l, m, r,  x,  y, -z)),
!>
!>           sign_x: Xx/X0,
!>           sign_y: Xy/X0,
!>           sign_z: Xz/X0,
!>
!>           ax: (1 - sign_x)/2,
!>           ay: (1 - sign_y)/2,
!>           az: (1 - sign_z)/2,
!>           axy: (1 - sign_x*sign_y)/2,
!>           axz: (1 - sign_x*sign_z)/2,
!>           ayz: (1 - sign_y*sign_z)/2,
!>           axyz: (1 - sign_x*sign_y*sign_z)/2,
!>
!>           C1: 0,
!>           Ci: axyz,
!>           Cs: ax,
!>           C2: axy,
!>           D2: axz + 2*ayz,
!>           C2h: az + 2*axy,
!>           C2v: ax + 2*ay,
!>           D2h: ax + 2*ay + 4*az,
!>
!>           printf(true, "reference(:,~3d) = [~2d,~2d,~2d,~2d,~2d,~2d,~2d,~2d ]~%", l*l+l+m+1, C1, Ci, Cs, C2, D2, C2h, C2v, D2h)
!>       )
!>   )$
!> \endcode
!>
program test_xlm_mgvn

   use const_gbl,    only: sym_op_nam_len
   use symmetry_gbl, only: Xlm_mgvn

   implicit none

   integer, parameter :: lmax = 6
   integer :: irr, irr0, l, m, pg, reference(8, (lmax + 1)**2), nsymel(8)
   character(len=sym_op_nam_len) :: symel(3, 8)

   ! symmetry operations
   nsymel(1) = 0;  symel(:, 1) = [ '   ', '   ', '   ' ]  ! C1
   nsymel(2) = 1;  symel(:, 2) = [ 'xyz', '   ', '   ' ]  ! Ci
   nsymel(3) = 1;  symel(:, 3) = [ 'x  ', '   ', '   ' ]  ! Cs
   nsymel(4) = 1;  symel(:, 4) = [ 'xy ', '   ', '   ' ]  ! C2
   nsymel(5) = 2;  symel(:, 5) = [ 'xz ', 'yz ', '   ' ]  ! D2
   nsymel(6) = 2;  symel(:, 6) = [ 'z  ', 'xy ', '   ' ]  ! C2h
   nsymel(7) = 2;  symel(:, 7) = [ 'x  ', 'y  ', '   ' ]  ! C2v
   nsymel(8) = 3;  symel(:, 8) = [ 'x  ', 'y  ', 'z  ' ]  ! D2h

   ! irreducible representation per point group (C1, Ci, Cs, C2, D2, C2h, C2v, D2h) obtained from Maxima
   reference(:,  1) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{0,+0}
   reference(:,  2) = [ 0, 1, 0, 1, 2, 2, 2, 2 ]  ! X_{1,-1}
   reference(:,  3) = [ 0, 1, 0, 0, 3, 1, 0, 4 ]  ! X_{1,+0}
   reference(:,  4) = [ 0, 1, 1, 1, 1, 2, 1, 1 ]  ! X_{1,+1}
   reference(:,  5) = [ 0, 0, 1, 0, 3, 0, 3, 3 ]  ! X_{2,-2}
   reference(:,  6) = [ 0, 0, 0, 1, 1, 3, 2, 6 ]  ! X_{2,-1}
   reference(:,  7) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{2,+0}
   reference(:,  8) = [ 0, 0, 1, 1, 2, 3, 1, 5 ]  ! X_{2,+1}
   reference(:,  9) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{2,+2}
   reference(:, 10) = [ 0, 1, 0, 1, 2, 2, 2, 2 ]  ! X_{3,-3}
   reference(:, 11) = [ 0, 1, 1, 0, 0, 1, 3, 7 ]  ! X_{3,-2}
   reference(:, 12) = [ 0, 1, 0, 1, 2, 2, 2, 2 ]  ! X_{3,-1}
   reference(:, 13) = [ 0, 1, 0, 0, 3, 1, 0, 4 ]  ! X_{3,+0}
   reference(:, 14) = [ 0, 1, 1, 1, 1, 2, 1, 1 ]  ! X_{3,+1}
   reference(:, 15) = [ 0, 1, 0, 0, 3, 1, 0, 4 ]  ! X_{3,+2}
   reference(:, 16) = [ 0, 1, 1, 1, 1, 2, 1, 1 ]  ! X_{3,+3}
   reference(:, 17) = [ 0, 0, 1, 0, 3, 0, 3, 3 ]  ! X_{4,-4}
   reference(:, 18) = [ 0, 0, 0, 1, 1, 3, 2, 6 ]  ! X_{4,-3}
   reference(:, 19) = [ 0, 0, 1, 0, 3, 0, 3, 3 ]  ! X_{4,-2}
   reference(:, 20) = [ 0, 0, 0, 1, 1, 3, 2, 6 ]  ! X_{4,-1}
   reference(:, 21) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{4,+0}
   reference(:, 22) = [ 0, 0, 1, 1, 2, 3, 1, 5 ]  ! X_{4,+1}
   reference(:, 23) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{4,+2}
   reference(:, 24) = [ 0, 0, 1, 1, 2, 3, 1, 5 ]  ! X_{4,+3}
   reference(:, 25) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{4,+4}
   reference(:, 26) = [ 0, 1, 0, 1, 2, 2, 2, 2 ]  ! X_{5,-5}
   reference(:, 27) = [ 0, 1, 1, 0, 0, 1, 3, 7 ]  ! X_{5,-4}
   reference(:, 28) = [ 0, 1, 0, 1, 2, 2, 2, 2 ]  ! X_{5,-3}
   reference(:, 29) = [ 0, 1, 1, 0, 0, 1, 3, 7 ]  ! X_{5,-2}
   reference(:, 30) = [ 0, 1, 0, 1, 2, 2, 2, 2 ]  ! X_{5,-1}
   reference(:, 31) = [ 0, 1, 0, 0, 3, 1, 0, 4 ]  ! X_{5,+0}
   reference(:, 32) = [ 0, 1, 1, 1, 1, 2, 1, 1 ]  ! X_{5,+1}
   reference(:, 33) = [ 0, 1, 0, 0, 3, 1, 0, 4 ]  ! X_{5,+2}
   reference(:, 34) = [ 0, 1, 1, 1, 1, 2, 1, 1 ]  ! X_{5,+3}
   reference(:, 35) = [ 0, 1, 0, 0, 3, 1, 0, 4 ]  ! X_{5,+4}
   reference(:, 36) = [ 0, 1, 1, 1, 1, 2, 1, 1 ]  ! X_{5,+5}
   reference(:, 37) = [ 0, 0, 1, 0, 3, 0, 3, 3 ]  ! X_{6,-6}
   reference(:, 38) = [ 0, 0, 0, 1, 1, 3, 2, 6 ]  ! X_{6,-5}
   reference(:, 39) = [ 0, 0, 1, 0, 3, 0, 3, 3 ]  ! X_{6,-4}
   reference(:, 40) = [ 0, 0, 0, 1, 1, 3, 2, 6 ]  ! X_{6,-3}
   reference(:, 41) = [ 0, 0, 1, 0, 3, 0, 3, 3 ]  ! X_{6,-2}
   reference(:, 42) = [ 0, 0, 0, 1, 1, 3, 2, 6 ]  ! X_{6,-1}
   reference(:, 43) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{6,+0}
   reference(:, 44) = [ 0, 0, 1, 1, 2, 3, 1, 5 ]  ! X_{6,+1}
   reference(:, 45) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{6,+2}
   reference(:, 46) = [ 0, 0, 1, 1, 2, 3, 1, 5 ]  ! X_{6,+3}
   reference(:, 47) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{6,+4}
   reference(:, 48) = [ 0, 0, 1, 1, 2, 3, 1, 5 ]  ! X_{6,+5}
   reference(:, 49) = [ 0, 0, 0, 0, 0, 0, 0, 0 ]  ! X_{6,+6}

   do l = 0, lmax
      do m = -l, l
         do pg = 1, 8
            irr0 = reference(pg, l*l + l + m + 1)
            irr = Xlm_mgvn(nsymel(pg), symel(:, pg), l, m)
            if (irr /= irr0) then
               print '(*(a,i0))', 'Wrong MGVN for l = ', l, ', m = ', m, ', pg = ', pg, ': expected ', irr0, ', got ', irr
               stop 1
            end if
         end do
      end do
   end do

   print '(a)', 'All tests OK!'

end program
