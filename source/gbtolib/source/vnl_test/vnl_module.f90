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
module vnl_module
   use precisn
   use general_quadrature

   implicit none

   complex(kind=wp), parameter :: imu = (0.0_wp,1.0_wp)

   real(kind=wp), parameter, private :: smallest = tiny(1.0_wp)*1.0e+04_wp
   real(kind=wp), parameter, private :: eps = 1.0e-15_wp !accuracy setting for the continued fractions computation

   !> \todo implement a Chebyshev approximation of the Vnl in the range where the recurrence is inaccurate.
   !> \todo calculate the Vnl solving a system of equations with starting values determined using the Chebyshev approx.
   type, extends(bound_user_function) :: vnl_integrand
      integer :: n
      integer :: l
      real(kind=wp) :: r
      logical :: eval_re = .true. !evaluate the real part of vnl
   contains
      procedure :: eval => eval_vnl
   end type

contains

   subroutine lpn(n,x,pn,pd)
!
!  ===============================================
!  Purpose: Compute Legendre polynomials Pn(x)
!           and their derivatives Pn'(x)
!  Input :  x --- Argument of Pn(x)
!           n --- Degree of Pn(x) ( n = 0,1,...)
!  Output:  PN(n) --- Pn(x)
!           PD(n) --- Pn'(x)
!  ===============================================
!
   integer, intent(in) :: n
   real(kind=wp), intent(in) :: x
   real(kind=wp), intent(out) :: pn(0:n),pd(0:n)
   real(kind=wp) :: p0,p1,pf
   integer k

   pn(0)=1.0_wp
   pn(1)=x
   pd(0)=0.0_wp
   pd(1)=1.0_wp
   p0=1.0_wp
   p1=x
   do k=2,n
      pf=(2.0_wp*k-1.0_wp)/k*x*p1-(k-1.0_wp)/k*p0
      pn(k)=pf
      if (dabs(x).eq.1.0_wp) then
         pd(k)=0.5_wp*x**(k+1)*k*(k+1.0_wp)
      else
         pd(k)=k*(p1-x*pf)/(1.0_wp-x*x)
      endif
      p0=p1
      p1=pf
   enddo

   end subroutine

   real(wp) function eval_vnl(data,x) result(r)
      class(vnl_integrand) :: data
      real(wp), intent(in) :: x
      real(wp) :: p(0:data%l), pd(0:data%l)
      complex(wp) :: z

         !todo: note that the number of evaluations of z can be halved if the values n,l,r are saved after each eval and an inquiry is made at the start about the current vs. last parameters to see whether 
         !they have changed or not 
         call lpn(data%l,x,p,pd)
         z = -imu/data%r
   
         z = (-z-x)**data%n/(z-x)**(data%n+1)

         if (data%eval_re) then
            z = real(z)
         else
            z = aimag(z)
         endif

         r = p(data%l)*z

   end function eval_vnl

   function vnl_direct(r,n,l)
      use const, only: epsabs, epsrel, limit, lenw
      implicit none
      real(kind=wp), intent(in) :: r
      integer, intent(in) :: n, l
      real(kind=wp) :: vnl_direct

      integer :: iwork(limit)
      real(kind=wp) :: work(lenw)
      real(kind=wp), parameter :: rt2 = sqrt(2.0_wp)

      real(kind=wp) :: re, im, A, B, abserr
      complex(kind=wp) :: s
      integer :: neval, ier, last, j
      type(vnl_integrand) :: integrand

         if (n < 0 .or. l < 0 .or. r .le. 0.0_wp) stop "vnl_direct: (n < 0 .or. l < 0 .or. r .le. 0.0_wp)"

         integrand%n = n
         integrand%l = l
         integrand%r = r
         A = -1.0_wp
         B = 1.0_wp

         integrand%eval_re = .true.
         call dqags(integrand,A,B,epsabs,epsrel,re,abserr,neval,ier,limit,lenw,last,iwork,work)

         integrand%eval_re = .false.
         call dqags(integrand,A,B,epsabs,epsrel,im,abserr,neval,ier,limit,lenw,last,iwork,work)

         !s = (-imu)**(l)*imu
         if (mod(l+1,2) .eq. 0) then !l+1 is even
            s = 1.0_wp
         else !l+1 is odd
            s = imu
         endif
         j = l/2
         if (mod(j,2) .ne. 0) s = -s
 
         vnl_direct = real(-rt2*s/(r*2.0_wp)*(re + imu*im))

   end function vnl_direct

   !Radial potentials V_{nl}(r) for the coulomb resolution based on the Laguerre generator.
   !For given values of r,n_max,l_max this routine generates the potentials V_{nl}(r) for n=0,...,n_max; l=0,...,l_max
   !For r .ge. 1 and n .le 10 and l .le. 10 the potentials are obtained using the recurrent relations given in P.M.W. Gill, et al., Chem. Phys. 356 (2009) 86-90
   !For parameters outside this range the potentials are obtained from their definition  - equation (18) ibid. - using an adaptive numerical quadrature.
   !In this way the precision of the generated V_{nl}(r) is approximately 10d-10 for any combination of the input parameters.
   !todo vectorize over a list of r values
   subroutine vnl_potential(vnl,r,n_max,l_max)
      implicit none
      real(kind=wp), intent(out) :: vnl(0:l_max,0:n_max)
      real(kind=wp), intent(in) :: r
      integer, intent(in) :: n_max, l_max

      integer, parameter :: n_lim = 10, l_lim = 10
      real(kind=wp), parameter :: r_lim = 1.0_wp

      integer :: k, l, n, s, l_m, n_m
      real(kind=wp) :: theta, hlp(1:n_lim), vnl_back(0:l_lim,0:n_lim) !, delta(0:l_max,0:n_max), dp(0:l_max,0:n_max)
      complex(kind=wp) :: z, q(0:l_lim)

      !compute the vnl potentials for n,l,r values lying outside of the stability domain of the recurrent relations
      if (r < r_lim .or. n_max > n_lim .or. l_max > l_lim) then

         do n=0,n_max
            do l=0,l_max
               if (r < r_lim .or. n > n_lim .or. l > l_lim) vnl(l,n) = vnl_direct(r,n,l)
            enddo
         enddo

         !finish if everything has been calculated numerically
         if (r < r_lim) return

      endif

      !we can use the recurrences only for a range of n,l values: n .le. n_lim .and. l .le. l_lim
      n_m = min(n_lim,n_max)
      l_m = min(l_lim,l_max)

      !1: generate the Vn0 potentials
      theta = atan(r)

      do k=1,n_m !vectorized
         hlp(k) = sin(2*k*theta)/(k*theta)
      enddo

      vnl(0,0) = 0.0_wp
      s = 1
      do k=1,n_m
         s = -s
         vnl(0,k) = vnl(0,k-1) + s*hlp(k) !hlp(k)=sin(2*k*theta)/(k*theta)
      enddo

      vnl(0,0:n_m) = (vnl(0,0:n_m) + 1.0_wp)*sqrt(2.0_wp)*theta/r

      !2:
      !generate the V0l potentials
      z = imu/r
      call legendre_ql(q,z,l_m) !calculate Q_{l}(imu/r) for l=0,...,l_m
   
      z = sqrt(2.0_wp)/r
      do l=0,l_m
         z = z*imu
         vnl(l,0) = real(z*q(l))
      enddo

      !3: recurrence to generate Vn1 potentials
      !there is a mild loss of precision at points n1 and n1+1, where vnl(1,n1) and vnl(n1+1) changes sign
      if (l_m .ge. 1) then
         do n=1,n_m
            vnl(1,n)=-(vnl(0,n)+vnl(0,n-1))/r+vnl(1,n-1)
         enddo
      endif

      !4: recurrence to generate the Vnl potentials for the rest of the n,l values
      !unstable for r smaller than ~ 20a.u.!!!
      do n=0,n_m-1
         do l=1,l_m-1
            vnl(l+1,n+1) = l/(l+1.0_wp)*(vnl(l-1,n+1)-vnl(l-1,n)) - (2*l+1.0_wp)/((l+1.0_wp)*r)*(vnl(l,n+1)+vnl(l,n))+vnl(l+1,n)
         enddo
      enddo

      !formulation using the differences, this is a tiny bit more accurate
      !calculate the starting differences and the column with l=1,n=0,...,n_m
!      do n=0,n_m-1
!         dp(0,n) = vnl(0,n+1)+vnl(0,n)
!         delta(0,n) = vnl(0,n+1) - vnl(0,n)
!         delta(1,n) = -dp(0,n)/r
!         vnl(1,n+1) = delta(1,n) + vnl(1,n)
!         dp(1,n) = vnl(1,n+1) + vnl(1,n)
!      enddo

      !4: use the difference formulation of the recurrence
!      do n=0,n_m-1
!         do l=1,l_m-1
!            delta(l+1,n) = l/(l+1.0_wp)*delta(l-1,n) - dp(l,n)*(2*l+1)/(r*(l+1.0_wp))
!            vnl(l+1,n+1) = delta(l+1,n) + vnl(l+1,n)
!            dp(l+1,n) = vnl(l+1,n+1) + vnl(l+1,n)
!         enddo
!      enddo

      !starting from the corner n_m,l_m; this is OK for larger n,l but gets worse for smaller n,l; the optimal solution is to combine the results obtained using the forward and the backward recursion (see 5:)
      do l=0,l_m
         vnl_back(l,n_m) = vnl_direct(r,n_m,l)
      enddo

      do n=0,n_m
         vnl_back(l_m,n) = vnl_direct(r,n,l_m)
         if (l_m .ge. 1) vnl_back(l_m-1,n) = vnl_direct(r,n,l_m-1)
      enddo

      !4: recurrence to generate the Vnl potentials for the rest of the n,l values
      do n=n_m-1,0,-1
         do l=l_m-1,1,-1
            vnl_back(l-1,n) = -(2*l+1)/(r*l)*(vnl_back(l,n+1) + vnl_back(l,n)) - (l+1)/(l*1.0_wp)*(vnl_back(l+1,n+1) - vnl_back(l+1,n)) + vnl_back(l-1,n+1)
         enddo
      enddo

     !5: output the best values only - forward recurrent results for small n,l; backward recurrent results for large n,l
     do n=1,n_m
        do l=1,l_m
           if (n > n_m/2 .or. l > l_m/2) vnl(l,n) = vnl_back(l,n)
        enddo
     enddo

   end subroutine vnl_potential

   !computes a string of associated legendre functions of the second kind Q_{l}(z) for one complex z and l=0,...,L
   !for details of the method see: B.I.Schneider et al., CPC, 181 (2010), 2091-2097
   !These functions correspond to the Mathematica function: LegendreQ[l, 0, 3, z], i.e. to the Associated Legendre functions of the second kind of type 3; 
   !type 3 functions have a single branch cut going from -\infty to 1. Types 1,2 have a different branch structure.
   !todo assume z is an array of input values instead of just one value
   !todo rewrite the routines fr_q and br_q to vectorize over a list of different z values.
   subroutine legendre_ql(q,z,L)
      implicit none
      complex(kind=wp), intent(out) :: q(0:L)
      complex(kind=wp), intent(in) :: z
      integer, intent(in) :: L

         if (L < 0) stop "legendre_ql: L<0"

         !make sure we guarantee ~15 accurate digits of the result
         if (abs(z-1) .le. 10d-2) stop "legendre_ql: precision of Q_{l}(z) may be smaller than 14 digits."

         if (L .le. 1 .or. abs(z) < 1) then !the forward recurrence contains exact forms for L=0,1; for abs(z) < 1 the forward recurrence is stable
            call fr_q_quad(q,z,L)
         else !L > 1 .and. abs(z) .ge. 1: the backward recurrence is stable for abs(z) > 1
            if (z .eq. 1.0_wp) stop "legendre_ql: z .eq. 1" !Q_{l}(z) has a pole for z == 1
            call br_q(q,z,L)
         endif

   end subroutine legendre_ql

   !evaluation of Q_{l}(z) (complex z) using the forward recurrence for Q_{l}(z), l=0,...,L
   subroutine fr_q(q,z,L)
      implicit none
      complex(kind=wp), intent(out) :: q(0:L)
      complex(kind=wp), intent(in) :: z
      integer, intent(in) :: L

      integer :: l_it, n_1, n_2, n_3

         !if (L < 0) stop "fr_q: (L < 0)"
         !if (abs(z) > 1) stop "fr_q: (abs(z) < 1); backward recursions are unstable for this z"
         !if (z .eq. 1.0d0) stop "fr_q: (z .eq. 1.0d0)"

         q(0) = 0.5_wp*log((z+1)/(z-1))
         q(1) = z*q(0) - 1.0_wp

         !forward recurrence: Abramowitz and Stegun, eq. 8.5.3
         n_1 = 3
         n_2 = 1
         n_3 = 2
         do l_it = 2, L
            q(l_it) = (n_1*z*q(l_it-1) - n_2*q(l_it-2))/n_3
            n_1 = n_1 + 2
            n_2 = n_2 + 1
            n_3 = n_3 + 1
         enddo

   end subroutine fr_q

   !evaluation of Q_{l}(z) (complex z) using the forward recurrence for Q_{l}(z), l=0,...,L and quadruple precision; the result is output in the double precision array.
   !for large L values the forward recurrence is unstable (even for abs(z) < 1)
   !todo try to resolve this without using quad prec
   subroutine fr_q_quad(q,z,L)
      implicit none
      complex(kind=wp), intent(out) :: q(0:L)
      complex(kind=wp), intent(in) :: z
      integer, intent(in) :: L

      integer :: l_it, n_1, n_2, n_3
      complex(kind=ep1) :: z_quad, q_quad(0:L)

         !if (L < 0) stop "fr_q: (L < 0)"
         !if (abs(z) > 1) stop "fr_q: (abs(z) < 1); backward recursions are unstable for this z"
         !if (z .eq. 1.0d0) stop "fr_q: (z .eq. 1.0d0)"

         z_quad = real(z,ep1) + cmplx(0.0_ep1,1.0_ep1)*real(aimag(z),ep1)
         q_quad(0) = 0.5_ep1*log((z_quad+1.0_ep1)/(z_quad-1.0_ep1))
         q_quad(1) = z_quad*q_quad(0) - 1.0_ep1

         !forward recurrence: Abramowitz and Stegun, eq. 8.5.3
         n_1 = 3
         n_2 = 1
         n_3 = 2
         do l_it = 2, L
            q_quad(l_it) = (n_1*z_quad*q_quad(l_it-1) - n_2*q_quad(l_it-2))/n_3
            n_1 = n_1 + 2
            n_2 = n_2 + 1
            n_3 = n_3 + 1
         enddo

         q = q_quad

   end subroutine fr_q_quad

   !evaluation of Q_{l}(z) (complex z) using the backwards recurrence for Q_{l}(z), l=0,...,L
   subroutine br_q(q,z,L)
      implicit none
      complex(kind=wp), intent(out) :: q(0:L)
      complex(kind=wp), intent(in) :: z
      integer, intent(in) :: L

      integer :: l_it, n_1, n_2, n_3
      complex(kind=wp) :: c

         !if (L .le. 1) stop "br_q: (L .le. 1); use the explicit expressions"
         !if (abs(z) < 1) stop "br_q: (abs(z) < 1); backward recursions are unstable for this z"
         !if (z .eq. 1.0d0) stop "br_q: (z .eq. 1.0d0)"

         !calculate Q_{L}(z)/Q_{L-1}(z) and use it to set up Q_L and Q_{L-1} (unnormalized) starting values for the recursion
         q(L) = cf_q(z,L)*smallest
         q(L-1) = smallest

         !backward recurrent relations for the A-L functions: Abramowitz and Stegun, eq. 8.5.3
         n_1 = -L
         n_2 = 2*L-1
         n_3 = L-1
         do l_it = L-2,0,-1
            q(l_it) = (n_2*z*q(l_it+1) + n_1*q(l_it+2))/n_3
            n_1 = n_1 + 1
            n_2 = n_2 - 2
            n_3 = n_3 - 1
         enddo

         !normalization using the exact value of Q_{0}(z) = 0.5d0*log((z+1)/(z-1))
         c = 0.5_wp*log((z+1)/(z-1))/q(0)
         q = q*c !this does not vectorize by default using ifort v13.1.3 claiming ineffectivity

   end subroutine br_q

   !returns the value of Q_{L}(z)/Q_{L-1}(z) (for a complex value of z) evaluated using the continued fractions and the modified Lentz algorithm.
   function cf_q(z,L)
      implicit none
      complex(kind=wp) :: cf_q
      complex(kind=wp), intent(in) :: z
      integer, intent(in) :: L

      integer :: L0
      real(kind=wp) :: a, test
      complex(kind=wp) :: b, D, C, Del

         !if (L .le. 0) stop "cf_q: (L .le. 0)"

         cf_q = smallest
         C = cf_q
         D = 0.0_wp
         test = 1.0_wp
         a = 1.0_wp
         L0 = L
         b = 0.0_wp
         do
            b = ((2*L0 + 1.0_wp)*z)/L0

            D = b + a*D
            if (abs(D) == 0.0_wp) then
               D = smallest
            end if

            C = b + a/C
            if (abs(C) == 0.0_wp) then
               C = smallest
            end if

            D = 1.0_wp/D
            Del = C*D
            cf_q = cf_q*Del

            test = abs(Del - 1.0_wp)
            if (test .le. eps) exit

            a = -(L0 + 1.0_wp)/L0
            L0 = L0 + 1
         end do

   end function cf_q

end module vnl_module
