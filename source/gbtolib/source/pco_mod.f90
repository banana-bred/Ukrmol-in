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

!> \brief   Pseudo-continuum orbitals
!> \authors D Darby-Lewis
!> \date    2019
!>
module pco_gbl
   use precisn_gbl
   use const_gbl, only: level2

   private

   public rm_cont_gt_pco, generate_PCO_exponents

contains

   subroutine generate_PCO_exponents(min_l,max_l,max_num,num_PCOs,alpha0,beta,exponents,thrs)
      implicit none
      integer, intent(in) :: min_l, max_l, num_PCOs(min_l:max_l),max_num
      real(kind=cfp), intent(in)    :: alpha0(min_l:max_l), beta(min_l:max_l)
      real(kind=cfp), intent(inout) :: exponents(1:max_num,min_l:max_l), thrs(min_l:max_l)
      integer :: i, j

      ! Add the single gto functions as single primitive CGTOs
      do i=min_l,max_l
         if(thrs(i) .lt. 0.0_cfp) thrs(i) = abs( alpha0(i)*(beta(i)-1.0_cfp) )
         do j=1,num_PCOs(i)
            ! The PCOs are being added in here using the formula alpha = alpha0 * beta ** (j - 1)
            exponents(j,i) = alpha0(i)*beta(i)**(real(j,kind=cfp)-1.0_cfp)
         enddo !j
      enddo !i

   end subroutine

   subroutine rm_cont_gt_pco(thrs,exponents,min_l,max_l,no_cont_exps,max_num,n_c_shells, &
                             PCO_exponents,min_PCO_l,max_PCO_l,num_PCOs,max_num_PCOs)
      implicit none
      integer, intent(in) :: min_l, max_l,no_cont_exps(min_l:max_l),max_num
      integer, intent(in) :: min_PCO_l,max_PCO_l, num_PCOs(min_PCO_l:max_PCO_l), max_num_PCOs
      real(kind=cfp), intent(inout) :: exponents(1:max_num,min_l:max_l)
      real(kind=cfp), intent(in)    :: PCO_exponents(1:max_num_PCOs,min_PCO_l:max_PCO_l), thrs(min_PCO_l:max_PCO_l)
      integer :: i, j, n_c_shells,x
      real(kind=cfp) :: min_exp

      write(level2,'(/,"Checking Continum exponents and removing those greater than the smallest PCO components for the same l.")')
      do i=max(min_l,min_PCO_l),min(max_l,max_PCO_l)
         x=0
         min_exp=minval(PCO_exponents(1:num_PCOs(i),i))
         do j=1,no_cont_exps(i)
            if (exponents(j,i) .ge. min_exp - thrs(i) .and. min_exp .gt. 0  ) then
               write(level2,'(/,"Removing the continuum shell = ",i3," for l = ",i3,3F12.5/)') j,i,exponents(j,i),min_exp, thrs(i)
               exponents(j,i) = -1.0_cfp
               n_c_shells= n_c_shells - 1
               x=x+1
            end if
         end do
         write(level2,'(/,"Finished checking continum exponents for l = ",i3," removed ", i3)') i, x
      end do

   end subroutine

end module pco_gbl
