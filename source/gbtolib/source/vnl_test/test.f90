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
program vnl_test
   use precisn
   use vnl_module, only: vnl_potential, imu, vnl_direct

   implicit none

   integer, parameter :: L = 20, N = 20, points = 200
   real(kind=wp) :: r, vnl(0:L,0:N), v, diff, prec
   integer :: i, j, k

      prec = 0.0_wp

      !small r values
      do k=1,points
         r = 1.0_wp/(k*1.0_wp)
         call vnl_potential(vnl,r,N,L)
         do i=0,N
            do j=0,L
              v = vnl_direct(r,i,j)
              diff = abs((vnl(j,i)-v)/v)
              if (diff > prec) prec = diff
              write (*,'(2i3,4e25.15)') i,j,r,vnl(j,i),v,diff
            enddo
         enddo
      enddo

      !large r values
      do k=1,points
         r = k*1.0_wp
         call vnl_potential(vnl,r,N,L)
         do i=0,N
            do j=0,L
              v = vnl_direct(r,i,j)
              diff = abs((vnl(j,i)-v)/v)
              if (diff > prec) prec = diff
              write (*,'(2i3,4e25.15)') i,j,r,vnl(j,i),v,diff
            enddo
         enddo
      enddo

      !the worst precision of the V_{nl} potential value
      print *,prec

end program vnl_test
