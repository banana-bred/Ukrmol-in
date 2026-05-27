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

   integer, parameter :: L = 1, N = 10, points = 200
   real(kind=wp), parameter :: pi = 3.14_wp
   real(kind=wp) :: r, vnl(0:L,0:N), vnl1(0:L,0:N), pot, theta
   integer :: i, k

      r = 1.0_wp
      call vnl_potential(vnl1,r,N,L)

      do k=1,points
         theta = pi*k/(2.0_wp*points)
         r = tan(theta)
         call vnl_potential(vnl,r,N,L)
         pot = 0.0_wp
         do i=0,N
            pot = pot + vnl1(1,i)*vnl(1,i)
         enddo
         pot = 6.0_wp/pi*pot
         write(*,'(2e25.15)') theta, pot
      enddo

end program vnl_test
