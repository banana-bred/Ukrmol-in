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

!> \brief   Unit test for `symmetry_gbl::is_centre_symmetric`
!> \authors J Benda
!> \date    2025
!>
!> Verifies if the specified centre is symmetric with respect to the provided symmetry operations.
!>
program test_is_centre_symmetric

   use const_gbl,    only: sym_op_nam_len
   use precisn_gbl,  only: cfp
   use symmetry_gbl, only: is_centre_symmetric

   implicit none

   real(cfp) :: c0(3) = [0._cfp, 0._cfp, 0._cfp]
   real(cfp) :: cx(3) = [1._cfp, 0._cfp, 0._cfp]
   real(cfp) :: cy(3) = [0._cfp, 1._cfp, 0._cfp]
   real(cfp) :: cz(3) = [0._cfp, 0._cfp, 1._cfp]

   real(cfp) :: cxy(3) = [1._cfp, 1._cfp, 0._cfp]
   real(cfp) :: cxz(3) = [1._cfp, 0._cfp, 1._cfp]
   real(cfp) :: cyz(3) = [0._cfp, 1._cfp, 1._cfp]

   real(cfp) :: cxyz(3) = [1._cfp, 1._cfp, 1._cfp]

   character(len=sym_op_nam_len) :: c1(0)
   character(len=sym_op_nam_len) :: ci(1) = ['xyz']
   character(len=sym_op_nam_len) :: cs(1) = ['x  ']
   character(len=sym_op_nam_len) :: c2(1) = ['xy ']
   character(len=sym_op_nam_len) :: d2(2) = ['xz ', 'yz ']
   character(len=sym_op_nam_len) :: c2h(2) = ['z  ', 'xy ']
   character(len=sym_op_nam_len) :: c2v(2) = ['x  ', 'y  ']
   character(len=sym_op_nam_len) :: d2h(3) = ['x  ', 'y  ', 'z  ']

   call test(c0,   c1, .true.)
   call test(cx,   c1, .true.)
   call test(cy,   c1, .true.)
   call test(cz,   c1, .true.)
   call test(cxy,  c1, .true.)
   call test(cyz,  c1, .true.)
   call test(cxz,  c1, .true.)
   call test(cxyz, c1, .true.)

   call test(c0,   ci, .true.)
   call test(cx,   ci, .false.)
   call test(cy,   ci, .false.)
   call test(cz,   ci, .false.)
   call test(cxy,  ci, .false.)
   call test(cyz,  ci, .false.)
   call test(cxz,  ci, .false.)
   call test(cxyz, ci, .false.)

   call test(c0,   cs, .true.)
   call test(cx,   cs, .false.)
   call test(cy,   cs, .true.)
   call test(cz,   cs, .true.)
   call test(cxy,  cs, .false.)
   call test(cyz,  cs, .true.)
   call test(cxz,  cs, .false.)
   call test(cxyz, cs, .false.)

   call test(c0,   c2, .true.)
   call test(cx,   c2, .false.)
   call test(cy,   c2, .false.)
   call test(cz,   c2, .true.)
   call test(cxy,  c2, .false.)
   call test(cyz,  c2, .false.)
   call test(cxz,  c2, .false.)
   call test(cxyz, c2, .false.)

   call test(c0,   d2, .true.)
   call test(cx,   d2, .false.)
   call test(cy,   d2, .false.)
   call test(cz,   d2, .false.)
   call test(cxy,  d2, .false.)
   call test(cyz,  d2, .false.)
   call test(cxz,  d2, .false.)
   call test(cxyz, d2, .false.)

   call test(c0,   c2h, .true.)
   call test(cx,   c2h, .false.)
   call test(cy,   c2h, .false.)
   call test(cz,   c2h, .false.)
   call test(cxy,  c2h, .false.)
   call test(cyz,  c2h, .false.)
   call test(cxz,  c2h, .false.)
   call test(cxyz, c2h, .false.)

   call test(c0,   c2v, .true.)
   call test(cx,   c2v, .false.)
   call test(cy,   c2v, .false.)
   call test(cz,   c2v, .true.)
   call test(cxy,  c2v, .false.)
   call test(cyz,  c2v, .false.)
   call test(cxz,  c2v, .false.)
   call test(cxyz, c2v, .false.)

   call test(c0,   d2h, .true.)
   call test(cx,   d2h, .false.)
   call test(cy,   d2h, .false.)
   call test(cz,   d2h, .false.)
   call test(cxy,  d2h, .false.)
   call test(cyz,  d2h, .false.)
   call test(cxz,  d2h, .false.)
   call test(cxyz, d2h, .false.)

   print '(a)', 'All tests OK!'

contains

   subroutine test(centre, ops, expected)

      real(cfp), intent(in) :: centre(3)
      character(len=sym_op_nam_len), intent(in) :: ops(:)
      logical, intent(in) :: expected

      if (is_centre_symmetric(centre, size(ops), ops) .neqv. expected) then
         print '(a,3f4.1)', 'Centre ', centre, ' symmetry check failed.'
         print '(a,*(a))', 'Symmetry operations: ', ops
         print '(a,l1,a,l1)', 'Expected ', expected, ', got ', .not. expected
         stop 1
      end if

   end subroutine test

end program
