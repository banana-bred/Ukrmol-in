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
!> This module contains routines to calculate the 1-electron integrals in the B-spline basis.
!> \todo Add checking that the radial quadratures have high enough order to get exact results.
!> Integrals involving radial B-splines with a non-zero derivative at r=A (i.e. at the origin of the B-spline grid) are not evaluated.
module bto_integrals_mod
   use utils
   use precisn
   use bto_data_obj
   use coupling_obj
   implicit none

   !> Used to get various coupling coefficients, mostly the real Gaunt coefficients.
   type(couplings_type), private :: cpl

contains

   !> Constructs a quadrature grid for the B-spline basis described by a given knot sequence. The final quadrature grid comprises 
   !> effectively a series of quadratures over subintervals generated in between each pair of distinct knots.
   subroutine construct_bspline_quadrature_grid(knots,x,w,n,r_points,weights,n_total_points)
      use general_quadrature, only: gl_expand_A_B
      implicit none
      !IN/OUT:
      integer, intent(in) :: n
      real(kind=cfp), intent(in) :: knots(:), x(2*n+1), w(2*n+1)
      real(kind=cfp), allocatable :: r_points(:), weights(:)
      integer, intent(out) :: n_total_points

      integer :: i, n_points, n_knots, err, cnt

         n_knots = size(knots)

         if (n_knots .le. 1 .or. n .eq. 0) call xermsg('bto_integrals_mod','construct_quadrature_grid','Invalid knot grid or bad Gaussian quadrature rule.',1,1)

         !Calculate how many quadrature points there will be in total.
         n_points = 2*n+1 !number of points in the elementary Gaussian quadrature.
         n_total_points = 0
         do i=1,n_knots-1
            if (knots(i+1) > knots(i)) n_total_points = n_total_points + n_points
         enddo !i

         if (allocated(r_points)) deallocate(r_points)
         if (allocated(weights)) deallocate(weights)

         allocate(r_points(n_total_points),weights(n_total_points),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','construct_quadrature_grid','Memory allocation failed.',err,1)
         r_points = 0.0_cfp; weights = 0.0_cfp

         !Construct the full grid expanding the elementary Gaussian rule on each distinct interval of knots.
         cnt = 0
         do i=1,n_knots-1
            if (knots(i+1) > knots(i)) then
               call gl_expand_A_B(x,w,n,r_points(cnt+1:cnt+n_points),weights(cnt+1:cnt+n_points),knots(i),knots(i+1))
               cnt = cnt + n_points
            endif
         enddo !i

   end subroutine construct_bspline_quadrature_grid

   !> Calculates overlap and kinetic energy integrals over a set of BTO orbitals with angular momentum up to l=max_l. The tail integrals
   !> for KEI are automatically subtracted.
   subroutine olap_kei_bto(cms_bto,max_l,first_index,olap,kei,int_index)
      use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
      use bspline_base
      implicit none
      !IN/OUT:
      type(bto_data), intent(inout) :: cms_bto
      integer :: max_l, first_index
      real(kind=cfp), allocatable :: olap(:), kei(:)
      integer, allocatable :: int_index(:,:)

      integer :: i, j, n_total_points, err, ind, A_ind, B_ind, l, m, lm, n_angular, n_radial, bloch_ind, last_index
      real(kind=cfp), allocatable :: r_points(:), weights(:), B_vals(:,:,:), temp_r(:), bspline_boundary_val(:,:)
      integer, allocatable :: bspline_start_end(:,:), bto_indices(:,:)
      real(kind=cfp) :: bto_norm, radial_olap, radial_kei_1, radial_kei_2, der_val, test, val, bloch_el
      logical :: bloch

         if (max_l < 0) call xermsg('bto_integrals_mod','olap_kei_bto','max_l < 0 but it must be .ge. 0.',1,1)

         err = cms_bto%check()
         if (err .ne. 0) call xermsg('bto_integrals_mod','olap_kei_bto','bto_data%check failed with an error message.',err,1)

         call generate_bto_indices(cms_bto,max_l,first_index,bto_indices,last_index)

         !Construct the quadrature grid that will be used to integrate over all pairs of B-splines
         call construct_bspline_quadrature_grid(cms_bto%knots,x_10,w_10,n_10,r_points,weights,n_total_points)

         !Find the start and end of each B-spline in the quadrature grid.
         call map_knots_to_grid(cms_bto%knots,cms_bto%order,cms_bto%n,r_points,bspline_start_end)
         
         allocate(B_vals(2,n_total_points,cms_bto%n),temp_r(n_total_points),bspline_boundary_val(2,cms_bto%n),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','olap_kei_bto','Memory allocation 1 failed.',err,1)
         B_vals = 0.0_cfp
         cms_bto%bcoef = 0.0_cfp
         bloch = .false.
         bspline_boundary_val = 0.0_cfp
         do ind=cms_bto%ind_0_der,cms_bto%n
            cms_bto%bcoef(ind) = 1.0_cfp
 
            !Calculate norm of the radial B-spline
            cms_bto%ind = ind
            cms_bto%r1 = cms_bto%knots(ind)
            cms_bto%r2 = cms_bto%knots(ind+cms_bto%order)
            call cms_bto%normalize
            bto_norm = cms_bto%norm
 
            !Evaluate it on the grid
            do i=bspline_start_end(1,ind),bspline_start_end(2,ind)
               !no derivative:
               B_vals(1,i,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,0,r_points(i),cms_bto%inbv,cms_bto%work)
               !2nd derivative:
               B_vals(2,i,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,2,r_points(i),cms_bto%inbv,cms_bto%work)
            enddo !i

            !Calculate the value and the first derivative of the B-spline at the boundary:
            bspline_boundary_val(1,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,0,cms_bto%B,cms_bto%inbv,cms_bto%work)
            bspline_boundary_val(2,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,1,cms_bto%B,cms_bto%inbv,cms_bto%work)

            cms_bto%bcoef(ind) = 0.0_cfp !clean-up for the next
         enddo !ind

         !Precalculate 1/r**2 on the quadrature grid.
         !The weights could be included here but to keep the code transparent we don't do it here.
         temp_r(1:n_total_points) = 1.0_cfp/(r_points(1:n_total_points))**2
        
         !Calculate the overlap integrals using the integrals over all unique pairs of the radial B-splines:
         n_angular = (max_l+1)**2
         n_radial = cms_bto%n*(cms_bto%n+1)/2
         i = n_radial*n_angular
!
         allocate(olap(i),kei(i),int_index(2,i),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','olap_kei_bto','Memory allocation 2 failed.',err,1)
         olap = 0.0_cfp; kei = 0.0_cfp; int_index = 0
!
         do i=cms_bto%ind_0_der,cms_bto%n
            do j=cms_bto%ind_0_der,i

               !Integrate over the radial interval where the B-splines overlap:
               A_ind = max(bspline_start_end(1,i),bspline_start_end(1,j))
               B_ind = min(bspline_start_end(2,i),bspline_start_end(2,j))

               if (B_ind .le. A_ind) cycle

               !Radial overlap integral:
               radial_olap = sum(B_vals(1,A_ind:B_ind,i)*B_vals(1,A_ind:B_ind,j)*weights(A_ind:B_ind))

               !Two integrals that enter the evaluation of the KEI for a pair of BTOs:
               radial_kei_1 = sum(B_vals(1,A_ind:B_ind,i)*B_vals(2,A_ind:B_ind,j)*weights(A_ind:B_ind))
               radial_kei_2 = sum(B_vals(1,A_ind:B_ind,i)*B_vals(1,A_ind:B_ind,j)*temp_r(A_ind:B_ind)*weights(A_ind:B_ind))

               !Save the integrals to the appropriate place in the final array of integrals:
               !we use the fact that the overlaps are non-zero only for (l1,m1) .eq. (l2,m2) = (l,m)
               do l=0,max_l
                  do m=-l,l
                     lm = l*l+l+m+1
                     !save the integrals as a 2D array (lm) X (i,j)
                     ind = lm + (i*(i-1)/2+j-1)*n_angular
                     olap(ind) = radial_olap
                     kei(ind) = -0.5_cfp*(radial_kei_1 - l*(l+1)*radial_kei_2)
                     !Calculate the Bloch element (most of them are zero but their calculation is so fast it doesn't matter):
                     bloch_el = 0.5_cfp*bspline_boundary_val(1,i)*bspline_boundary_val(2,j)
                     kei(ind) = kei(ind) +  bloch_el
                     !generate the int_index values:
                     int_index(1,ind) = max(bto_indices(lm,i),bto_indices(lm,j))
                     int_index(2,ind) = min(bto_indices(lm,i),bto_indices(lm,j))
                  enddo
               enddo !l
            enddo !j
         enddo !i

   end subroutine olap_kei_bto

   !> Generates property integrals centered on the CMS for the whole set of BTOs.
   subroutine property_ints_bto(cms_bto,max_l,first_index,max_prop_l,prop,int_index)
      use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
      use bspline_base
      use phys_const, only: fourpi
      implicit none
      !IN/OUT:
      type(bto_data), intent(inout) :: cms_bto
      integer :: max_l, first_index, max_prop_l
      real(kind=cfp), allocatable :: prop(:)
      integer, allocatable :: int_index(:,:)

      integer :: i, j, n_total_points, err, ind, A_ind, B_ind, l, m, lm, n_angular, n_radial, n_ang_prop, l1, m1, l2, m2, l1m1, l2m2, last_index
      real(kind=cfp), allocatable :: r_points(:), weights(:), B_vals(:,:), temp_r(:,:)
      real(kind=cfp) :: fac
      integer, allocatable :: bspline_start_end(:,:), bto_indices(:,:)
      real(kind=cfp) :: bto_norm, radial_prop 

         if (max_l < 0 .or. max_prop_l < 0) call xermsg('bto_integrals_mod','property_ints_bto','max_l < 0 or max_prop_l < 0 but they must be .ge. 0.',1,1)

         err = cms_bto%check()
         if (err .ne. 0) call xermsg('bto_integrals_mod','property_ints_bto','bto_data%check failed with an error message.',err,1)

         call generate_bto_indices(cms_bto,max_l,first_index,bto_indices,last_index)

         !Precalculate the coupling coefficients: real Gaunt cfs.
         call cpl%prec_cgaunt(max(max_l,max_prop_l))

         !Construct the quadrature grid that will be used to integrate over all pairs of B-splines
         call construct_bspline_quadrature_grid(cms_bto%knots,x_10,w_10,n_10,r_points,weights,n_total_points)

         !Find the start and end of each B-spline in the quadrature grid.
         call map_knots_to_grid(cms_bto%knots,cms_bto%order,cms_bto%n,r_points,bspline_start_end)
         
         allocate(B_vals(n_total_points,cms_bto%n),temp_r(n_total_points,0:max(1,max_prop_l)),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','property_ints_bto','Memory allocation 1 failed.',err,1)
         B_vals = 0.0_cfp; temp_r = 0.0_cfp
         cms_bto%bcoef = 0.0_cfp
         do ind=cms_bto%ind_0_der,cms_bto%n
            cms_bto%bcoef(ind) = 1.0_cfp
 
            !Calculate norm of the radial B-spline
            cms_bto%ind = ind
            cms_bto%r1 = cms_bto%knots(ind)
            cms_bto%r2 = cms_bto%knots(ind+cms_bto%order)
            call cms_bto%normalize
            bto_norm = cms_bto%norm
 
            !Evaluate it on the grid
            do i=bspline_start_end(1,ind),bspline_start_end(2,ind)
               !no derivative:
               B_vals(i,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,0,r_points(i),cms_bto%inbv,cms_bto%work)
            enddo !i

            cms_bto%bcoef(ind) = 0.0_cfp !clean-up for the next
         enddo !ind

         !Precalculate sqrt(fourpi/(2*l+1.0_cfp))*r**l on the quadrature grid
         do l=0,max_prop_l
            fac = sqrt(fourpi/(2*l+1.0_cfp))
            do i=1,n_total_points
               temp_r(i,l) = fac*r_points(i)**l
            enddo !i
         enddo !l
        
         !Calculate the overlap integrals using the integrals over all unique pairs of the radial B-splines:
         n_angular = (max_l+1)**2
         n_angular = n_angular*(n_angular+1)/2 * (max_prop_l+1)**2
         n_radial = cms_bto%n*(cms_bto%n+1)/2
         i = n_radial*n_angular
!
         allocate(prop(i),int_index(3,i),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','property_ints_bto','Memory allocation 2 failed.',err,1)
         prop = 0.0_cfp; int_index = 0

         lm = (max_l+1)**2
         n_angular = lm*(lm+1)/2
         n_ang_prop = (max_prop_l+1)**2

         do i=cms_bto%ind_0_der,cms_bto%n
            do j=cms_bto%ind_0_der,i

               !Integrate over the radial interval where the B-splines overlap:
               A_ind = max(bspline_start_end(1,i),bspline_start_end(1,j))
               B_ind = min(bspline_start_end(2,i),bspline_start_end(2,j))

               if (B_ind .le. A_ind) cycle

               do l=0,max_prop_l

                  !Radial property integral:
                  radial_prop = sum(B_vals(A_ind:B_ind,i)*B_vals(A_ind:B_ind,j)*temp_r(A_ind:B_ind,l)*weights(A_ind:B_ind))

                  do m=-l,l
                     lm = l*l+l+m+1
                     fac =(-1)**m

                     do l1=0,max_l
                        do m1=-l1,l1
                           l1m1 = l1*l1+l1+m1+1
                           do l2=abs(l-l1),l1 !0,l1
                              do m2=-l2,l2
                                 l2m2 = l2*l2+l2+m2+1
                                 if (l1m1 .ge. l2m2) then
                                    !save in the order: (l1m1,l2m2,lm,ij)
                                    ind = l1m1*(l1m1-1)/2+l2m2 + (lm-1)*n_angular + (i*(i-1)/2+j-1)*n_angular*n_ang_prop
                                    prop(ind) = fac*radial_prop*cpl%rgaunt(l,l1,l2,m,m1,m2)
                                    !generate the int_index values:
                                    int_index(1,ind) = max(bto_indices(l1m1,i),bto_indices(l2m2,j))
                                    int_index(2,ind) = min(bto_indices(l1m1,i),bto_indices(l2m2,j))
                                    int_index(3,ind) = lm
                                 endif
                              enddo !m2
                           enddo !l2
                        enddo !m1
                     enddo !l1

                  enddo !m

               enddo !l

            enddo !j
         enddo !i

   end subroutine property_ints_bto

   !> Generates nuclear attraction integrals for the whole basis of BTOs and the given set of centers and their charges.
   subroutine nuclear_attraction_ints_bto(cms_bto,max_l,first_index,centers,charges,n_centers,nai,int_index)
      use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
      use bspline_base
      use phys_const, only: fourpi
      use special_functions, only: cfp_resh
      implicit none
      !IN/OUT:
      type(bto_data), intent(inout) :: cms_bto
      integer :: max_l, first_index, n_centers
      real(kind=cfp), intent(in) :: centers(3,n_centers), charges(n_centers)
      real(kind=cfp), allocatable :: nai(:)
      integer, allocatable :: int_index(:,:)

      integer :: i, j, k, n_total_points, err, ind, A_ind, B_ind, l, m, lm, n_angular, n_radial, l1, m1, l2, m2, l1m1, l2m2, last_index
      real(kind=cfp), allocatable :: r_points(:), weights(:), B_vals(:,:), temp_r(:,:,:), radial_int(:,:), xlm(:,:), re_sph_harm_center(:,:)
      real(kind=cfp) :: fac, R
      integer, allocatable :: bspline_start_end(:,:), bto_indices(:,:)
      real(kind=cfp) :: bto_norm, s, t(n_centers)

         if (max_l < 0) call xermsg('bto_integrals_mod','nuclear_attraction_ints_bto','max_l < 0 but it must be .ge. 0.',1,1)

         err = cms_bto%check()
         if (err .ne. 0) call xermsg('bto_integrals_mod','nuclear_attraction_ints_bto','bto_data%check failed with an error message.',err,1)

         do i=1,n_centers
            R = sqrt(dot_product(centers(1:3,i),centers(1:3,i)))
            !This assumption is easy to overcome. It is only the first
            !implementation that is not completely general.
            if (R > cms_bto%A) call xermsg('bto_integrals_mod','nuclear_attraction_ints_bto','In the current implementation the B-splines&
                                           &are allowed to start only beyond the nucleus furtherst away from CMS.',1,1)
         enddo !i

         call generate_bto_indices(cms_bto,max_l,first_index,bto_indices,last_index)

         !Precalculate the coupling coefficients: real Gaunt cfs.
         call cpl%prec_cgaunt(2*max_l)

         !Construct the quadrature grid that will be used to integrate over all pairs of B-splines
         call construct_bspline_quadrature_grid(cms_bto%knots,x_10,w_10,n_10,r_points,weights,n_total_points)

         !Find the start and end of each B-spline in the quadrature grid.
         call map_knots_to_grid(cms_bto%knots,cms_bto%order,cms_bto%n,r_points,bspline_start_end)
         
         allocate(B_vals(n_total_points,cms_bto%n),temp_r(n_total_points,0:max(1,2*max_l),n_centers),&
                 &xlm(-2*max_l:max(1,2*max_l),0:max(1,2*max_l)),re_sph_harm_center(n_centers,(2*max_l+1)**2),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','nuclear_attraction_ints_bto','Memory allocation 1 failed.',err,1)
         B_vals = 0.0_cfp; temp_r = 0.0_cfp
         cms_bto%bcoef = 0.0_cfp
         do ind=cms_bto%ind_0_der,cms_bto%n
            cms_bto%bcoef(ind) = 1.0_cfp
 
            !Calculate norm of the radial B-spline
            cms_bto%ind = ind
            cms_bto%r1 = cms_bto%knots(ind)
            cms_bto%r2 = cms_bto%knots(ind+cms_bto%order)
            call cms_bto%normalize
            bto_norm = cms_bto%norm
 
            !Evaluate it on the grid
            do i=bspline_start_end(1,ind),bspline_start_end(2,ind)
               !no derivative:
               B_vals(i,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,0,r_points(i),cms_bto%inbv,cms_bto%work)
            enddo !i

            cms_bto%bcoef(ind) = 0.0_cfp !clean-up for the next
         enddo !ind

         !Precalculate fourpi/(2*l+1.0_cfp)*r**2*R**l/r**(l+1) on the quadrature grid
         do i=1,n_centers
            R = sqrt(dot_product(centers(1:3,i),centers(1:3,i)))
            do l=0,2*max_l
               fac = fourpi/(2*l+1.0_cfp)
               do j=1,n_total_points
                  temp_r(j,l,i) = fac*r_points(j)*(R/r_points(j))**l
               enddo !j
            enddo !l
         enddo !i

         !Precalculate the values of the real spherical harmonics for each nucleus:
         if (max_l > 0) then
            do k=1,n_centers
               call cfp_resh(xlm,centers(1,k),centers(2,k),centers(3,k),2*max_l)
               do l=0,2*max_l
                  do m=-l,l
                     lm = l*l+l+m+1
                     re_sph_harm_center(k,lm) = xlm(m,l)
                  enddo !m
               enddo !l
            enddo !k
         else !l .eq. 0
            re_sph_harm_center = 1.0_cfp/sqrt(fourpi)
         endif
        
         !Calculate the overlap integrals using the integrals over all unique pairs of the radial B-splines:
         n_angular = (max_l+1)**2
         n_angular = n_angular*(n_angular+1)/2
         n_radial = cms_bto%n*(cms_bto%n+1)/2
         i = n_radial*n_angular
!
         allocate(nai(i),int_index(2,i),radial_int(n_centers,0:max(1,2*max_l)),stat=err)
         if (err .ne. 0) call xermsg('bto_integrals_mod','property_ints_bto','Memory allocation 2 failed.',err,1)
         nai = 0.0_cfp; int_index = 0

         lm = (max_l+1)**2
         n_angular = lm*(lm+1)/2

         do i=cms_bto%ind_0_der,cms_bto%n
            do j=cms_bto%ind_0_der,i

               !Integrate over the radial interval where the B-splines overlap:
               A_ind = max(bspline_start_end(1,i),bspline_start_end(1,j))
               B_ind = min(bspline_start_end(2,i),bspline_start_end(2,j))

               if (B_ind .le. A_ind) cycle

               do k=1,n_centers
                  do l=0,2*max_l
                     !Radial integral from the Legendre expansion assuming the particle described by BTOs is the electron:
                     radial_int(k,l) = -charges(k)*sum(B_vals(A_ind:B_ind,i)*B_vals(A_ind:B_ind,j)*temp_r(A_ind:B_ind,l,k)*weights(A_ind:B_ind))
                  enddo !l
               enddo !k

               do l1=0,max_l
                  do m1=-l1,l1
                     l1m1 = l1*l1+l1+m1+1
                     do l2=0,l1
                        do m2=-l2,l2
                           l2m2 = l2*l2+l2+m2+1
                           if (l1m1 .ge. l2m2) then
                              !save in the order: (l1m1,l2m2,lm,ij)
                              ind = l1m1*(l1m1-1)/2+l2m2 + (i*(i-1)/2+j-1)*n_angular
                              do l=0,2*max_l
                                 s = 0.0_cfp
                                 do m=-l,l
                                    fac = cpl%rgaunt(l,l1,l2,m,m1,m2)
                                    lm = l*l+l+m+1
                                    do k=1,n_centers
                                       s = s +fac*radial_int(k,l)*re_sph_harm_center(k,lm)
                                    enddo !k
                                 enddo !m
                                 !Accumulate the final integral from the contributions of the terms from the Leg. expansion:
                                 nai(ind) = nai(ind) + s
                                 !generate the int_index values:
                                 int_index(1,ind) = max(bto_indices(l1m1,i),bto_indices(l2m2,j))
                                 int_index(2,ind) = min(bto_indices(l1m1,i),bto_indices(l2m2,j))
                                 !write(stdout,'(6i4,e25.15)') i,j, l1,m1,l2,m2,nai(ind)
                              enddo !l
                           endif
                        enddo !m2
                     enddo !l2
                  enddo !m1
               enddo !l1

            enddo !j
         enddo !i

   end subroutine nuclear_attraction_ints_bto

end module bto_integrals_mod
