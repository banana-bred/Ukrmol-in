! Copyright 2024
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
module grid_gbl
  use basis_data_generic_gbl, only: CGTO_shell_data_obj
  use bspline_grid_gbl
  use common_obj_gbl, only: nucleus_type
  use precisn_gbl
  use const_gbl, only: stdout, level2, level3
  use utils_gbl, only: xermsg
  implicit none
  
  private

  public grid_r1_r2_obj, radial_grid, radial_grid_CGTO_pair, grid_r1_r2!, grid_r1_r2_test

  type :: grid_r1_r2_obj
     !> Quadrature points and weights for r1 and r2 radial grids for the Legendre expansion. 
     !> The r1 grid spans the radial distance from the start of the B-spline grid up to the larger of: the end of the B-spline grid and the (adjusted) R-matrix radius.
     !> The value of the adjusted R-matrix radius is obtained as the end point of the last B-spline that starts in the inner region and extends to the outer region.
     !> The r2 grid is restricted to the (adjusted) inner region only: it is not needed in the outer region.
     real(kind=cfp), allocatable :: r1(:), w1(:), r2(:), w2(:)
     integer :: n1_total_points,n2_total_points,max_bspline_l,max_prop_l,first_bspline_index,n_unique_pairs,max_l_legendre
     !> Grid of radial B-splines which will be used to calculate the auxiliary functions.
     type(bspline_grid_obj) :: bspline_grid
     !> R-matrix radius: used to determine which knots lie inside the sphere and which lie outside.
     !> The adjusted value corresponds to the end-point of the last B-spline which starts in the inner region and (perhaps) extends into the outer region.
     real(kind=cfp), private :: rmat_radius = -1.0_cfp, rmat_radius_adjusted = -1.0_cfp
     !> Index of the last B-spline which lies in the inner region (r <= R-matrix sphere).
     integer :: last_bspline_inner = 0
     !> Normalized B-splines (and pairs of B-splines) evaluated on the r1 and r2 grids. We don't need BB pairs on the r2 grid.
     real(kind=cfp), allocatable :: B_vals_r1(:,:), B_vals_r2(:,:), BB_vals_r1(:,:)
     !> The evaluation of the 1-electron mixed integrals involves quadrature over r1 grid involving a product of GTO-related part and a BTO-related part. 
     !> The part related to the BTOs is stored in the arrays below. The radial part of the BTO is defined as: B(r)/r. 
     !> The values stored include multiplication by the radial part of the Jacobian (r**2) and the quadrature weights for the r1 grid.
     real(kind=cfp), allocatable :: bto_radial_olap(:,:), bto_radial_kei(:,:,:), bto_radial_prop(:,:,:), bto_end_points(:,:)
     !> Starting and ending points of the B-splines (and pairs of B-splines) on the r1 and r2 grids. 
     !> The value in bspline_start_end_*(3,ind) contains the index of the last end-point of the B-spline with index ind such that it is never larger than the last point in the inner region.
     integer, allocatable :: bspline_start_end_r1(:,:), bspline_start_end_r2(:,:), BB_start_end_r1(:,:)
     !> For a given pair of B-splines it contains the index of the corresponding unique pair of B-splines.
     integer, allocatable :: pairs_to_unique_pairs(:,:)
     !> For a given index of the unique pair of B-splines it contains the indices of the two corresponding B-splines.
     integer, allocatable :: unique_pairs_to_pairs(:,:)
     !> Y_l function on the r1 grid for all pairs of B-splines.
     real(kind=cfp), allocatable :: lambda_BB_r1_l_ij(:,:,:)
     !> Lebedev grid and the real spherical harmonics evaluated at the Lebedev angular points.
     real(kind=cfp), allocatable :: Xlm_Lebedev(:), leb_r1(:), leb_r2(:), leb_r3(:), leb_w(:)
     !> Normalization factors for the radial B-splines.
     real(kind=cfp), allocatable :: bto_norm(:)
     !> Number of points of the chosen Lebedev rule.
     integer :: lebedev_order = -1
     !> Set to .true. following a call to init.
     logical, private :: initialized = .false.
  contains
     !> Initializes the grid given the B-spline grid, the elementary quadrature rules, etc.
     procedure :: init
     !> Constructs the quadrature grids appropriate for the Legendre expansion and evaluates other auxiliary arrays.
     procedure :: construct_r1_r2_grids
     !> Calculates the Y function for each unique pair of radial B-splines.
     procedure :: eval_lambda_l_BB
     !> Construct the Lebedev grid of a given (or nearest larger) order.
     procedure :: construct_lebedev_grid
     !> Constructs angular grid quadrature grid based on the Gauss-Kronrod product rule.
     procedure :: construct_GK_angular_grid
     !> Construct the Lebedev grid and evaluate the real spherical harmonics on the angular Lebedev grid.
     procedure :: eval_Xlm_on_lebedev_grid
     !> Deallocates everything.
     procedure :: final
  end type grid_r1_r2_obj

  !> Grid that is used for various calculations.
  type(grid_r1_r2_obj) :: grid_r1_r2

  ! Grid that is used in testing of the B-spline two-p-continuum integrals
  !type(grid_r1_r2_obj) :: grid_r1_r2_test

contains

  subroutine init(this,bspline_grid,first_bspline_index,max_bspline_l,max_prop_l,nuclei,delta_r1,n1,n2,only_on_bto_grid,&
                  rmat_radius,dipole_damp_factor,max_l_legendre)
     use general_quadrature_gbl, only: get_gaussrule !,n_10,w_10,x_10, n_7,w_7,x_7
     implicit none
     class(grid_r1_r2_obj) :: this
     type(nucleus_type) :: nuclei(:)
     type(bspline_grid_obj) :: bspline_grid
     integer, intent(in) :: n1, n2, first_bspline_index
     integer, intent(in) :: max_bspline_l, max_prop_l, max_l_legendre
     logical, intent(in) :: only_on_bto_grid
     real(kind=cfp), intent(in) :: delta_r1, rmat_radius, dipole_damp_factor

     real(kind=cfp) :: x1(n1), w1(n1), x2(n2), w2(n2), quad_start, quad_end
     integer :: n

        write(level3,'("--------->","grid_r1_r2_obj:init")')

        if (this%initialized) call this%final

        call get_gaussrule(n1,x1,w1)

        call get_gaussrule(n2,x2,w2)

        call this%construct_r1_r2_grids(bspline_grid,first_bspline_index,max_bspline_l,max_prop_l,dipole_damp_factor,&
                                        max_l_legendre,nuclei,delta_r1,x1,w1,n1,x2,w2,n2,only_on_bto_grid,rmat_radius)

        this%initialized = .true.

        write(level3,'("<---------","grid_r1_r2_obj:init")')

  end subroutine init

  subroutine construct_r1_r2_grids(this,inp_bspline_grid,first_bspline_index,max_bspline_l,max_prop_l,dipole_damp_factor,&
                                   max_l_legendre,nuclei,delta_r1,x1,w1,n1,x2,w2,n2,only_on_bto_grid,rmat_radius)
     use general_quadrature_gbl, only: gl_expand_A_B
     use sort_gbl, only: sort_float
     use bspline_base_gbl, only: bvalu, map_knots_to_grid
     use phys_const_gbl, only: fourpi
     implicit none
     class(grid_r1_r2_obj) :: this
     type(nucleus_type) :: nuclei(:)
     type(bspline_grid_obj) :: inp_bspline_grid
     integer, intent(in) :: n1, n2, first_bspline_index
     real(kind=cfp), intent(in) :: x1(n1), w1(n1), x2(n2), w2(n2), delta_r1, rmat_radius, dipole_damp_factor
     integer, intent(in) :: max_bspline_l, max_prop_l, max_l_legendre
     logical, intent(in) :: only_on_bto_grid

     real(kind=cfp), allocatable :: centers(:)
     real(kind=cfp) :: R_min, R_max
     integer :: n_points, err, i, j, k, l, ind, max_bspline_l_adj ,max_prop_l_adj, ind1, ind2, s_point, e_point
     real(kind=cfp) :: quad_start, quad_end, bto_val_2
     real(kind=cfp), allocatable :: list(:,:)
     real(kind=cfp) :: fac, r_start, r_end, test
     real(kind=cfp), allocatable :: tmp_r(:), tmp_w(:)

     real(kind=cfp), allocatable :: r(:,:), eval_CGTO_shell(:,:)

        write(level3,'("--------->","grid_r1_r2_obj:construct_r1_r2_grids")')

        this%bspline_grid = inp_bspline_grid
        this%max_bspline_l = max_bspline_l
        this%max_prop_l = max_prop_l
        this%max_l_legendre = max_l_legendre
        this%first_bspline_index = first_bspline_index
        this%rmat_radius = rmat_radius

        if (only_on_bto_grid) then
           R_min = this%bspline_grid%A
        else
           R_min = 0.0_cfp
        endif
        R_max = this%bspline_grid%B

        if (allocated(this%bto_norm)) deallocate(this%bto_norm)
        allocate(centers(size(nuclei)),this%bto_norm(this%bspline_grid%n),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation 1 failed.',err,1)
        do i=1,size(nuclei)
           centers(i) = sqrt(dot_product(nuclei(i)%center,nuclei(i)%center))
        enddo

        !todo replace with a more flexible grid determination, i.e. delta_r1
        write(level2,'("Maximum length of the r1 quadrature-grids: ",e25.15)') delta_r1
        call radial_grid(centers,R_min,R_max,this%bspline_grid,delta_r1,x1,w1,n1,this%r1,this%w1,this%n1_total_points)
!
!------ Determine index of the last B-spline which lies at least partially inside the R-matrix sphere.
        this%last_bspline_inner = this%bspline_grid%get_last_inner_bspline(this%rmat_radius)
!
!------ Determine where the last B-spline straddling the R-matrix sphere ends (this%rmat_radius_adjusted) and use this value in the routine below which maps the knots to the grid.
!       The mixed integrals over the inner-region B-splines and the continuum GTOs (which we assume always extend to the outer region) are evaluated using
!       the inner region routines. Setting the end-point of the inner-region quadrature to the end point of the last B-spline which straddles the
!       R-matrix radius ensures that these mixed integrals will be evaluated integrating over the full radial domain of the B-spline so the overlap with the continuum GTO is taken
!       into account correctly also in the outer region.
        call this%bspline_grid%bspline_range(this%last_bspline_inner,r_start,this%rmat_radius_adjusted)

        if (this%rmat_radius_adjusted .ne. this%rmat_radius) then
           write(stdout,'("End-point for the inner-region quadratures will be adjusted to: ",e25.15)') this%rmat_radius_adjusted
        endif
!
!------ Determine mapping of the endpoints of each B-spline with indices in the r1 array.
        call map_knots_to_grid(this%bspline_grid%knots,&
                               this%bspline_grid%order,&
                               this%bspline_grid%n,&
                               this%r1,&
                               this%rmat_radius_adjusted,&
                               this%bspline_start_end_r1)
!
!------ Evaluate the B-splines on the grid
        if (allocated(this%B_vals_r1)) deallocate(this%B_vals_r1)
        allocate(this%B_vals_r1(this%n1_total_points,this%bspline_grid%n),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation failed.',err,1)

        if (allocated(this%bto_radial_olap)) deallocate(this%bto_radial_olap)
        if (allocated(this%bto_radial_kei)) deallocate(this%bto_radial_kei)
        if (allocated(this%bto_radial_prop)) deallocate(this%bto_radial_prop)
        if (allocated(this%bto_end_points)) deallocate(this%bto_end_points)

        max_prop_l_adj = max(1,max_prop_l)
        max_bspline_l_adj = max(1,max_bspline_l)

        allocate(this%bto_radial_olap(this%n1_total_points,this%bspline_grid%n),&
                 this%bto_radial_kei(this%n1_total_points,this%bspline_grid%n,0:max_bspline_l_adj),&
                 this%bto_radial_prop(this%n1_total_points,this%bspline_grid%n,0:max_prop_l_adj),&
                 this%bto_end_points(2,this%bspline_grid%n),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation 2 failed.',err,1)

        this%B_vals_r1 = 0.0_cfp
        this%bspline_grid%bcoef = 0.0_cfp
        do ind=this%first_bspline_index,this%bspline_grid%n
           this%bspline_grid%bcoef(ind) = 1.0_cfp

           this%bto_norm(ind) = this%bspline_grid%normalize(ind)

           !Evaluate the normalized B-spline on the grid
           do i=this%bspline_start_end_r1(1,ind),this%bspline_start_end_r1(3,ind)
              this%B_vals_r1(i,ind) = this%bto_norm(ind)&
                                      *bvalu(this%bspline_grid%knots,&
                                             this%bspline_grid%bcoef,&
                                             this%bspline_grid%n,&
                                             this%bspline_grid%order,&
                                             0,&
                                             this%r1(i),&
                                             this%bspline_grid%inbv,&
                                             this%bspline_grid%work)

              !Evaluate the 2nd derivative of the B-spline times norm
              bto_val_2 = this%bto_norm(ind)&
                                      *bvalu(this%bspline_grid%knots,&
                                             this%bspline_grid%bcoef,&
                                             this%bspline_grid%n,&
                                             this%bspline_grid%order,&
                                             2,&
                                             this%r1(i),&
                                             this%bspline_grid%inbv,&
                                             this%bspline_grid%work)

              this%bto_radial_olap(i,ind) = this%w1(i)*this%B_vals_r1(i,ind)*this%r1(i) !B-spline part times Jac of the radial integrand times the quadrature weight
              do l=0,max_bspline_l_adj
                 this%bto_radial_kei(i,ind,l) = this%w1(i)*(bto_val_2*this%r1(i)-l*(l+1)/this%r1(i)*this%B_vals_r1(i,ind)) !precision loss is possible for small values of r
              enddo !l
              do l=0,max_prop_l_adj
                 this%bto_radial_prop(i,ind,l) = this%w1(i)*this%B_vals_r1(i,ind)*this%r1(i)**(l+1)
              enddo !l
              !include the exponential damping for the dipole operator if requested
              if (dipole_damp_factor /= 0.0_cfp) then
                 l=1
                 this%bto_radial_prop(i,ind,l) = this%bto_radial_prop(i,ind,l) * exp(-dipole_damp_factor*this%r1(i))
              endif
           enddo !i

           !Indices in the angular_integrals_at_knots of the r1, r2 points   
           r_start = this%bspline_grid%knots(ind)
           r_end = this%bspline_grid%knots(ind + this%bspline_grid%order)

           !first derivative of the B-spline at its end points: only the B-splines at the end points of the interval may give rise to Bloch terms
           this%bto_end_points(1,ind) = this%bto_norm(ind)&
                                        *bvalu(this%bspline_grid%knots,&
                                               this%bspline_grid%bcoef,&
                                               this%bspline_grid%n,&
                                               this%bspline_grid%order,&
                                               1,&
                                               r_start,&
                                               this%bspline_grid%inbv,&
                                               this%bspline_grid%work)*r_start !1/r*Jac*BTO_derivative at r=r_start
           this%bto_end_points(2,ind) = this%bto_norm(ind)&
                                        *bvalu(this%bspline_grid%knots,&
                                               this%bspline_grid%bcoef,&
                                               this%bspline_grid%n,&
                                               this%bspline_grid%order,&
                                               1,&
                                               r_end,&
                                               this%bspline_grid%inbv,&
                                               this%bspline_grid%work)*r_end !1/r*Jac*BTO_derivative at r=r_end

           this%bspline_grid%bcoef(ind) = 0.0_cfp !clean-up for the next B-spline
        enddo !ind

        k = this%n1_total_points+size(this%bspline_grid%knots)+size(centers)
        allocate(list(k,1),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation 0 failed.',err,1)

        n_points = n2

        k = 0
        do i=1,this%n1_total_points
           if ((this%r1(i) .ge. this%bspline_grid%A) .and. (this%r1(i) .le. this%rmat_radius_adjusted)) then
              k = k + 1
              list(k,1) = this%r1(i)
           endif
        enddo

        do i=1,size(this%bspline_grid%knots)
           k = k + 1
           list(k,1) = this%bspline_grid%knots(i)
        enddo

        do i=1,size(centers)
           if (centers(i) >= this%bspline_grid%A) then
              k = k + 1
              list(k,1) = centers(i)
           endif
        enddo

        call sort_float(k,1,list)

        !Construct the r2 quadrature grid making sure the r1 points and knots are taken as the end points of the quadrature subintervals.
        !This is crucial to integrate accurately around the Coulomb cusp.
        !Dry run to determine array sizes.
        test = 10*F1MACH(4,cfp_dummy)
        this%n2_total_points = 0
        do i=2,k
           if (min(list(i,1)-test,this%bspline_grid%B) > list(i-1,1)) then
              this%n2_total_points = this%n2_total_points + n_points
           endif
        enddo !i

        if (allocated(this%r2)) deallocate(this%r2)
        if (allocated(this%w2)) deallocate(this%w2)
        allocate(this%r2(this%n2_total_points),this%w2(this%n2_total_points),tmp_r(n_points),tmp_w(n_points),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation failed.',err,1)

        this%n2_total_points = 0
        do i=2,k
           if (min(list(i,1)-test,this%bspline_grid%B) > list(i-1,1)) then
              quad_start = list(i-1,1)
              quad_end =   list(i  ,1)
              !print *,'w2',w2
              !print *,'range 2',quad_start,quad_end,tmp_w
              call gl_expand_A_B(x2,w2,(n_points-1)/2,tmp_r,tmp_w,quad_start,quad_end)
              this%r2(this%n2_total_points+1:this%n2_total_points+n_points) = tmp_r(1:n_points)
              this%w2(this%n2_total_points+1:this%n2_total_points+n_points) = tmp_w(1:n_points)
              this%n2_total_points = this%n2_total_points + n_points
           endif
        enddo !i

        write(level2,'("Total number of r1 and r2 points: ",2i15)') this%n1_total_points, this%n2_total_points
!
!------ Determine mapping of the endpoints of each B-spline with indices in the r2 array
        call map_knots_to_grid(this%bspline_grid%knots,&
                               this%bspline_grid%order,&
                               this%bspline_grid%n,&
                               this%r2,&
                               this%rmat_radius_adjusted,&
                               this%bspline_start_end_r2)

!
!------ Evaluate the B-splines on the grid
        if (allocated(this%B_vals_r2)) deallocate(this%B_vals_r2)
        allocate(this%B_vals_r2(this%n2_total_points,this%bspline_grid%n),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation failed.',err,1)

        this%B_vals_r2 = 0.0_cfp
        this%bspline_grid%bcoef = 0.0_cfp
        do ind=1,this%bspline_grid%n
           this%bspline_grid%bcoef(ind) = 1.0_cfp

           !Evaluate the normalized B-spline on the grid
           do i=this%bspline_start_end_r2(1,ind),this%bspline_start_end_r2(3,ind)

              this%B_vals_r2(i,ind) = this%bto_norm(ind)&
                                      *bvalu(this%bspline_grid%knots,&
                                             this%bspline_grid%bcoef,&
                                             this%bspline_grid%n,&
                                             this%bspline_grid%order,&
                                             0,&
                                             this%r2(i),&
                                             this%bspline_grid%inbv,&
                                             this%bspline_grid%work)
           enddo !i

           this%bspline_grid%bcoef(ind) = 0.0_cfp !clean-up for the next
        enddo !ind
!
!------ Evaluate pairs of B-splines on the r1 and r2 grids:
!
        !Get the mapping between each pair of B-splines and the unique set of overlapping B-splines
        call this%bspline_grid%get_unique_pairs(this%pairs_to_unique_pairs,this%n_unique_pairs)

        write(level2,'("Number of unique overlapping pairs of radial B-splines: ",i15)') this%n_unique_pairs

        !
        !Find the mapping between the unique pairs and the B-spline indices
        !
        if (allocated(this%unique_pairs_to_pairs)) deallocate(this%unique_pairs_to_pairs)
        allocate(this%unique_pairs_to_pairs(2,this%n_unique_pairs),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation 9 failed.',err,1)

        do ind1=1,this%bspline_grid%n
           do ind2=1,ind1

              j = this%pairs_to_unique_pairs(ind1,ind2)
              if (j .eq. 0) cycle !the b-splines do not overlap

              this%unique_pairs_to_pairs(1:2,j) = (/ind1,ind2/)

           enddo !ind2
        enddo !ind1
        !
        !
        !
        if (allocated(this%BB_vals_r1)) deallocate(this%BB_vals_r1)
        if (allocated(this%BB_start_end_r1)) deallocate(this%BB_start_end_r1)
        allocate(this%BB_vals_r1(this%n1_total_points,this%n_unique_pairs),&
        &this%BB_start_end_r1(2,this%n_unique_pairs),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_r1_r2_grids','Memory allocation 10 failed.',err,1)
        this%BB_vals_r1 = 0.0_cfp
        this%BB_start_end_r1(1,:) = 1
        this%BB_start_end_r1(2,:) = 0

        !Evaluate
        do ind1=1,this%bspline_grid%n
           do ind2=1,ind1
              j = this%pairs_to_unique_pairs(ind1,ind2)
              if (j .eq. 0) cycle !the b-splines do not overlap
              !Note that here the start and end points are chosen regardless of
              !whether the B-splines span the R-matrix radius (or are completely outside of it): we use the values from the 3rd row of this%bspline_start_end_r*
              !This is necessary since when integrating over a product of B-splines we always want to integrate over their full domain - not only over their inner region parts.
              s_point = max(this%bspline_start_end_r1(1,ind1),this%bspline_start_end_r1(1,ind2))
              e_point = min(this%bspline_start_end_r1(3,ind1),this%bspline_start_end_r1(3,ind2))
              if (s_point .le. e_point) then
                  !values on the r1 grid:
                  this%BB_start_end_r1(1,j) = s_point
                  this%BB_start_end_r1(2,j) = e_point
                  do i=s_point,e_point
                     this%BB_vals_r1(i,j) = this%B_vals_r1(i,ind1)*this%B_vals_r1(i,ind2)
                  enddo
              else
                  call xermsg('grid_mod','construct_r1_r2_grids','Error in r2 grid or in bspline_start_end_r2.',1,1)
              endif
           enddo !ind2
        enddo !ind1

        !B-spline testing
!        deallocate(r,eval_CGTO_shell)

        write(level3,'("<---------","grid_r1_r2_obj:construct_r1_r2_grids")')

  end subroutine construct_r1_r2_grids

  subroutine final(this)
     implicit none
     class(grid_r1_r2_obj) :: this

        if (allocated(this%r1)) deallocate(this%r1)
        if (allocated(this%w1)) deallocate(this%w1)
        if (allocated(this%r2)) deallocate(this%r2)
        if (allocated(this%w2)) deallocate(this%w2)
        if (allocated(this%B_vals_r1)) deallocate(this%B_vals_r1)
        if (allocated(this%B_vals_r2)) deallocate(this%B_vals_r2)
        if (allocated(this%pairs_to_unique_pairs)) deallocate(this%pairs_to_unique_pairs)
        if (allocated(this%BB_vals_r1)) deallocate(this%BB_vals_r1)
        if (allocated(this%bto_radial_olap)) deallocate(this%bto_radial_olap)
        if (allocated(this%bto_radial_kei)) deallocate(this%bto_radial_kei)
        if (allocated(this%bto_radial_prop)) deallocate(this%bto_radial_prop)
        if (allocated(this%bto_end_points)) deallocate(this%bto_end_points)
        if (allocated(this%bspline_start_end_r1)) deallocate(this%bspline_start_end_r1)
        if (allocated(this%bspline_start_end_r2)) deallocate(this%bspline_start_end_r2)
        if (allocated(this%BB_start_end_r1)) deallocate(this%BB_start_end_r1)
        if (allocated(this%lambda_BB_r1_l_ij)) deallocate(this%lambda_BB_r1_l_ij)
        if (allocated(this%Xlm_Lebedev)) deallocate(this%Xlm_Lebedev)
        if (allocated(this%leb_r1)) deallocate(this%leb_r1)
        if (allocated(this%leb_r2)) deallocate(this%leb_r2)
        if (allocated(this%leb_r3)) deallocate(this%leb_r3)
        if (allocated(this%leb_w)) deallocate(this%leb_w)
        if (allocated(this%bto_norm)) deallocate(this%bto_norm)

        this%n_unique_pairs = 0
        this%lebedev_order = -1
        this%rmat_radius = -1.0_cfp
        this%rmat_radius_adjusted = -1.0_cfp
        this%last_bspline_inner = 0
        this%initialized = .false.

  end subroutine final

  subroutine construct_lebedev_grid(this,n)
     use phys_const_gbl, only: fourpi
     use lebedev_gbl
     implicit none
     class(grid_r1_r2_obj) :: this
     integer, intent(in) :: n

     integer :: i, err, n_rule, available

        if (n .le. 0 .or. this%max_bspline_l < 0) stop "error in input"

        !Loop over the radial r1 points and calculate the nuclear attraction integrals at the Lebedev points
        !todo instead of +2 there should be the cgto max l
        n_rule = min(max(n,this%max_bspline_l + 2),rule_max)
        do i=n_rule,rule_max
           available = available_table(i)
           if (available == 1) then
              !In case of CMS-only functions the order must be high enough to integrate the product of the spherical harmonics otherwise the rule can fall on the nodes of Xlm
              this%lebedev_order = order_table(i)
              print *,'lebedev order',this%lebedev_order

              if (allocated(this%leb_r1)) deallocate(this%leb_r1)
              if (allocated(this%leb_r2)) deallocate(this%leb_r2)
              if (allocated(this%leb_r3)) deallocate(this%leb_r3)
              if (allocated(this%leb_w)) deallocate(this%leb_w)
              allocate(this%leb_r1(this%lebedev_order),&
                       this%leb_r2(this%lebedev_order),&
                       this%leb_r3(this%lebedev_order),&
                       this%leb_w(this%lebedev_order),stat=err)
              if (err /= 0) call xermsg('grid_mod','construct_lebedev_grid','Memory allocation failed.',err,1)

              call ld_by_order(this%lebedev_order,this%leb_r1,this%leb_r2,this%leb_r3,this%leb_w)
              print *,'got leb',this%lebedev_order
              exit
           endif
        enddo !i

  end subroutine construct_lebedev_grid

  subroutine construct_GK_angular_grid(this,n_div,r,w1,n_points)
     use phys_const_gbl, only: twopi, pi
     use lebedev_gbl
     implicit none
     class(grid_r1_r2_obj) :: this
     integer, intent(in) :: n_div
     integer, intent(out) :: n_points
     real(kind=cfp), allocatable :: r(:,:), w1(:)

     integer, parameter :: n = 8
     real(kind=cfp), parameter :: x(2*n+1) = (/0.003310062059141922032055965490164602_cfp, &
                                               0.019855071751231884158219565715263505_cfp, &
                                               0.052939546576271789025819491230874338_cfp, &
                                               0.101666761293186630204223031762084782_cfp, &
                                               0.163822964527420661421844630953584475_cfp, &
                                               0.237233795041835507091130475405376825_cfp, &
                                               0.319649451035934021403725688851554280_cfp, &
                                               0.408282678752175097530261928819908010_cfp, &
                                               0.500000000000000000000000000000000000_cfp, &
                                               0.591717321247824902469738071180091990_cfp, &
                                               0.680350548964065978596274311148445720_cfp, &
                                               0.762766204958164492908869524594623175_cfp, &
                                               0.836177035472579338578155369046415525_cfp, &
                                               0.898333238706813369795776968237915218_cfp, &
                                               0.947060453423728210974180508769125662_cfp, &
                                               0.980144928248768115841780434284736495_cfp, &
                                               0.996689937940858077967944034509835398_cfp/)
     real(kind=cfp), parameter :: w(2*n+1) = (/0.00891119166035517757639348060137489490_cfp, &
                                               0.024719697501069654250181984723498447_cfp, &
                                               0.0412411494656791653443125967228039477_cfp, &
                                               0.055823185413419806611054079466970742_cfp, &
                                               0.0681315546275861076311693726272531016_cfp, &
                                               0.078326303084094200245124044243484369_cfp, &
                                               0.0860353042776056559286474401019285433_cfp, &
                                               0.09070001253401732153087426258627522_cfp, &
                                               0.0922232028723458217644854778528214649_cfp, &
                                               0.09070001253401732153087426258627522_cfp, &
                                               0.0860353042776056559286474401019285433_cfp, &
                                               0.078326303084094200245124044243484369_cfp, &
                                               0.0681315546275861076311693726272531016_cfp, &
                                               0.055823185413419806611054079466970742_cfp, &
                                               0.0412411494656791653443125967228039477_cfp, &
                                               0.024719697501069654250181984723498447_cfp, &
                                               0.00891119166035517757639348060137489490_cfp/)

     integer :: i, j, k, err, p, q, ni, n_phi_div
     real(kind=cfp) :: th_max, phi_max, th, phi, cth, sth, cphi, sphi, th_min, phi_min, dth, dphi, &
                       intervals(2*n_div+1), alp, A, B

        if (n .le. 0) stop "error in input"

        n_phi_div = n_div !max(1,n_div/2)
        n_points = (2*n+1)*2*n_div*(2*n+1)*n_phi_div

        if (allocated(r)) deallocate(r)
        if (allocated(w1)) deallocate(w1)
        allocate(r(3,n_points),w1(n_points),stat=err)
        if (err /= 0) call xermsg('grid_mod','construct_GK_angular_grid','Memory allocation failed.',err,1)

!        !Interval [0;pi/2] divided as ~1/theta where points accumulate towards theta = 0.0

        alp = 1.0_cfp/10.0_cfp
        ni = n_div
        A = pi/2.0_cfp/(1.0_cfp/alp - 1.0_cfp/abs(alp+ni))
        B = pi/2.0_cfp*(1.0_cfp - 1.0_cfp/(1.0_cfp-alp/abs(alp+ni)))
        do i=0,n_div
           intervals(i+1) = A/abs(i-alp-ni) + B
        enddo
        do i=1,n_div
           intervals(n_div+1+i) = (pi/2.0_cfp-intervals(n_div+1-i))+pi/2.0_cfp
        enddo
        print *,'ints',intervals(1:2*n_div+1)

        k = 0
        do p=1,n_div*2
           !th_max = intervals(p) !pi/(n_div*1.0_cfp)*p
           !th_min = intervals(p+1) !pi/(n_div*1.0_cfp)*(p-1)
           th_max = intervals(p+1) !pi/(n_div*1.0_cfp)*p
           th_min = intervals(p) !pi/(n_div*1.0_cfp)*(p-1)
           dth = th_max-th_min
           do i=1,2*n+1
              th = x(i)*dth + th_min
              cth = cos(th)
              sth = sin(th)
              do q=1,n_phi_div
                 phi_max = twopi/(n_phi_div*1.0_cfp)*q
                 phi_min = twopi/(n_phi_div*1.0_cfp)*(q-1)
                 dphi = phi_max-phi_min
                 do j=1,2*n+1
                    k = k + 1
                    phi = x(j)*dphi + phi_min
                    cphi = cos(phi)
                    sphi = sin(phi)
                    r(1,k) = sth*cphi
                    r(2,k) = sth*sphi
                    r(3,k) = cth
                    w1(k) = w(i)*w(j)*sth*dphi*dth
                 enddo !j
              enddo !q
           enddo !i
        enddo !p

  end subroutine construct_GK_angular_grid

  subroutine eval_Xlm_on_lebedev_grid(this,n)
     use phys_const_gbl, only: fourpi
     use lebedev_gbl
     use special_functions_gbl, only: real_harmonics_obj
     implicit none
     class(grid_r1_r2_obj) :: this
     integer, intent(in) :: n

     integer :: point, l, m, err, lm, stride
     real(kind=cfp), allocatable :: RH(:,:)
     real(kind=cfp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     type(real_harmonics_obj) :: rh_obj

        call this%construct_lebedev_grid(n)

        if (allocated(this%Xlm_Lebedev)) deallocate(this%Xlm_Lebedev)
        allocate(this%Xlm_Lebedev(this%lebedev_order*(this%max_bspline_l+1)**2),&
                 RH(-this%max_bspline_l:this%max_bspline_l,0:this%max_bspline_l),stat=err)
        if (err /= 0) call xermsg('grid_mod','eval_Xlm_on_lebedev_grid','Memory allocation failed.',err,1)

        do point=1,this%lebedev_order
           if (this%max_bspline_l > 0) then
              call rh_obj%resh(RH,this%leb_r1(point),this%leb_r2(point),this%leb_r3(point),this%max_bspline_l)
           else
              RH(0,0) = norm
           endif
   
           lm = 0
           do l=0,this%max_bspline_l
              do m=-l,l
                 lm = lm + 1
                 this%Xlm_Lebedev(point+this%lebedev_order*(lm-1)) = RH(m,l)
              enddo !m
           enddo !l
        enddo !point

  end subroutine eval_Xlm_on_lebedev_grid

  subroutine eval_lambda_l_BB(this,max_l)
     use phys_const_gbl, only: fourpi
     use function_integration_gbl
     use general_quadrature_gbl
     use omp_lib
!     use const_gbl, only: epsabs, epsrel
     implicit none
     class(grid_r1_r2_obj) :: this
     integer, intent(in) :: max_l

     integer :: pair_index, err, i, si, ei, l, k, k2, B_start, B_end
     real(kind=cfp) :: fac, r12
     real(kind=cfp), allocatable :: f_l_w1_w2(:)

     integer :: ind1,ind2, neval, ier, last, n, n_intervals, thread_id, n_threads
     integer, parameter :: limit = 1000, lenw = 4*limit
     integer, allocatable :: iwork(:)
     real(kind=cfp), allocatable :: work(:)
     real(kind=cfp) :: res, abserr, A, B, res1, res2, epsabs, epsrel, res_interval
     real(kind=cfp), allocatable :: interval(:,:)
     type(BB_legendre_integrand) :: integrand

     ! FOR OLD VERSION
     integer :: max_l_d, j, B_i_start, B_i_end, B_j_start, B_j_end

        write(level3,'("--------->","grid_r1_r2_obj:eval_lambda_l_BB")')

        !Calculate the Y_l function for all unique pairs of radial B-splines
        max_l_d = min(2*this%max_bspline_l,this%max_l_legendre)
        pair_index = this%bspline_grid%n*(this%bspline_grid%n+1)/2
        allocate(this%lambda_BB_r1_l_ij(this%n1_total_points,0:max(max_l_d,1),pair_index),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','new_BG_mixed_2el_initialize','Memory allocation 2 failed.',err,1)

        this%lambda_BB_r1_l_ij(:,:,:) = 0.0_cfp

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP & PRIVATE(l,fac,k,k2,r12,f_l_w1_w2,i,j,pair_index,B_i_start,B_i_end,B_j_start,B_j_end,B_start,B_end,err) &
        !$OMP & SHARED(this,max_l_d)
        allocate(f_l_w1_w2(this%n2_total_points),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','new_BG_mixed_2el_initialize','Memory allocation 3 failed.',err,1)
        !$OMP DO SCHEDULE(DYNAMIC)
        do l=0,max_l_d
           fac = fourpi/(2*l+1.0_cfp)
           do k=1,this%n1_total_points !r1
              do k2=1,this%n2_total_points !r2
                 !Calculate w1*w2*fac*r1*r<**l/r>**(l+1) from the Legendre expansion
                 if (this%r2(k2) .le. this%r1(k)) then
                    r12 = this%r2(k2)/this%r1(k) !to make sure ifort 16.0.1 produces correct results
                    f_l_w1_w2(k2) = this%w2(k2)*this%w1(k)*fac*(r12)**l
                 else
                    r12 = this%r1(k)/this%r2(k2) !to make sure ifort 16.0.1 produces correct results
                    f_l_w1_w2(k2) = this%w2(k2)*this%w1(k)*fac*(r12)**(l+1)
                 endif
              enddo !k2

              !Now loop over all unique pairs of B-splines
              do i=this%first_bspline_index,this%bspline_grid%n
                 do j=this%first_bspline_index,i
                    pair_index = i*(i-1)/2+j
                    B_i_start = this%bspline_start_end_r2(1,i)
                    B_i_end = this%bspline_start_end_r2(2,i)
                    B_j_start = this%bspline_start_end_r2(1,j)
                    B_j_end = this%bspline_start_end_r2(2,j)
                    B_start = max(B_i_start,B_j_start)
                    B_end = min(B_i_end,B_j_end)
                    if (B_start > B_end) cycle !the B-splines do not overlap
                    !this%r1(k) = Jacobian(r1)/r1: 1/r1 comes from 1/r1*Y(r1).
                    this%lambda_BB_r1_l_ij(k,l,pair_index) = sum(f_l_w1_w2(B_start:B_end)*this%B_vals_r2(B_start:B_end,i) &
                                                            *this%B_vals_r2(B_start:B_end,j))
                 enddo !j
              enddo !i

           enddo !k
        enddo !l
        !$OMP END DO
        !$OMP END PARALLEL

        write(level3,'("<---------","grid_r1_r2_obj:eval_lambda_l_BB")')

  end subroutine eval_lambda_l_BB

  !> Determines the radial grid needed to describe, integrals involving the CGTO and the B-splines on the given grid.
  subroutine radial_grid(centers,R_min,R_max,bspline_grid,delta_cgto_grid,x,w,n_points,r_points,weights,n_total_points)
     use sort_gbl, only: sort_float
     use general_quadrature_gbl, only: gl_expand_A_B
     implicit none
     type(bspline_grid_obj), intent(in) :: bspline_grid
     integer, intent(in) :: n_points
     real(kind=cfp), intent(in) :: R_min, R_max, x(n_points), w(n_points), centers(:), delta_cgto_grid
     !OUTPUT variables:
     integer, intent(out) :: n_total_points
     real(kind=cfp), allocatable :: r_points(:), weights(:)

     integer :: err, interval, cnt, i, n_intervals, n_ranges, j, k, rng, n_int
     real(kind=cfp) :: d, A, B, R_lim, x_AB(n_points), w_AB(n_points), range_end, range_start, R_start, R_end, center, test
     real(kind=cfp), allocatable :: list(:,:), range_start_end(:,:), tmp_r(:), tmp_w(:)

         if (R_max .le. R_min) call xermsg('grid_mod','radial_grid','On input R_max .le. R_min: wrong input.',1,1)

         R_start = R_min
         R_end = R_max
         !R_end = min(R_max,bspline_grid%B)

         if (R_start > bspline_grid%A .and. R_start < bspline_grid%B) then
            !Find the nearest knot lying below the start point:
             do i=2,bspline_grid%no_knots
                if (bspline_grid%knots(i) > R_start) then
                   R_start = bspline_grid%knots(i-1)
                   exit
                endif
             enddo !i
         endif

         if (R_end > bspline_grid%A .and. R_end < bspline_grid%B) then
            !Find the nearest knot lying beyond the end point:
             do i=1,bspline_grid%no_knots
                if (bspline_grid%knots(i) > R_end) then
                   R_end = bspline_grid%knots(i)
                   exit
                endif
             enddo !i
         endif

         if (delta_cgto_grid .le. 0.0_cfp) then
            call xermsg('grid_mod','radial_grid','On input delta_cgto_grid was .le. 0.0_cfp.',3,1)
         endif

         j = ceiling((R_end-R_start)/delta_cgto_grid)
         allocate(list(bspline_grid%no_knots+size(centers)+j+2,1),range_start_end(2,bspline_grid%no_knots+size(centers)+j),stat=err)
         if (err /= 0) call xermsg('grid_mod','radial_grid','Memory allocation error.',err,1)
         list = 0.0_cfp
         range_start_end = 0.0_cfp

         k = 0
         do i=1,bspline_grid%no_knots
            if (bspline_grid%knots(i) .ge. R_start .and. bspline_grid%knots(i) .le. R_end) then
               k = k + 1
               list(k,1) = bspline_grid%knots(i)
            endif
         enddo

         do i=1,size(centers)
            if (centers(i) .ge. R_start) then
               k = k + 1
               list(k,1) = centers(i)
            endif
         enddo

         do i=1,j
            k = k + 1
            list(k,1) = min(R_start + i*delta_cgto_grid,R_end)
         enddo

         k = k + 1
         list(k,1) = R_start
         k = k + 1
         list(k,1) = R_end

         call sort_float(k,1,list)

         n_ranges = 0
         test = 10*F1MACH(4,cfp_dummy)
         do i=2,k
            if (min(list(i,1)-test,R_end) > list(i-1,1)) then
               n_ranges = n_ranges + 1
               range_start_end(1,n_ranges) = list(i-1,1)
               range_start_end(2,n_ranges) = list(i  ,1)
               !print *,'range 1',range_start_end(1:2,n_ranges)
            endif
         enddo !i

         !Total number of quadrature points in the interval [R_start,R_end].
         n_total_points = n_ranges*n_points

         if (allocated(r_points)) deallocate(r_points)
         if (allocated(weights)) deallocate(weights)
         allocate(r_points(n_total_points),weights(n_total_points),tmp_r(n_total_points),tmp_w(n_total_points),stat=err)
         if (err /= 0) call xermsg('grid_mod','radial_grid','Memory allocation failed.',err,1)

         !Construct quadratures for the individual ranges:
         cnt = 0
         do j=1,n_ranges
            A = range_start_end(1,j)
            B = range_start_end(2,j)
            !Prepare the quadrature grid for the interval [A,B] expanding the canonical G-L quadrature.
            !print *,'quad',A,B,cnt
            call gl_expand_A_B(x,w,(n_points - 1)/2,tmp_r,tmp_w,A,B)
            r_points(cnt+1:cnt+n_points) = tmp_r(1:n_points)
            weights(cnt+1:cnt+n_points) = tmp_w(1:n_points)
            cnt = cnt + n_points
         enddo
!
         if (cnt .ne. n_total_points) call xermsg('grid_mod','radial_grid','Error constructing the radial grid.',3,1)

  end subroutine radial_grid

  subroutine radial_grid_CGTO_pair(cgto_A,cgto_B,R_min,R_max,x,w,n_points,delta_cgto_grid,r_points,weights,n_total_points)
     use sort_gbl, only: sort_float
     use general_quadrature_gbl, only: gl_expand_A_B
     implicit none
     integer, intent(in) :: n_points
     real(kind=cfp), intent(in) :: R_min, R_max, x(n_points), w(n_points)
     type(CGTO_shell_data_obj), intent(in) :: cgto_A, cgto_B
     !OUTPUT variables:
     integer, intent(out) :: n_total_points
     real(kind=cfp), allocatable :: r_points(:), weights(:)
     real(kind=cfp), intent(in) :: delta_cgto_grid !todo this should be determined adaptively so that in between consecutive end points the CGTO falls down only by a given amount, etc.

     integer :: err, cnt, i, n_intervals, n_ranges, j, k, p, q
     real(kind=cfp) :: range_end, range_start, R_start, R_end, prod_alp, prod_P(3), RP
     real(kind=cfp), allocatable :: list(:,:), range_start_end(:,:), tmp_r(:), tmp_w(:)

         !write(stdout,'("CGTO negligible beyond radius [a.u.]: ",e25.15)') R_max

         !todo determine properly
         R_start = R_min
         R_end = R_max

         j = ceiling((R_max-R_min)/delta_cgto_grid)
         n_ranges = j+cgto_A%number_of_primitives*cgto_B%number_of_primitives+1
         allocate(list(n_ranges,1),range_start_end(2,n_ranges),stat=err)
         if (err /= 0) call xermsg('grid_mod','radial_grid_CGTO_pair','Memory allocation 1 failed.',err,1)

         !print *,'input',R_min,R_max
         !print *,'adj',R_start,R_end
 
         k = 0
         do i=0,j
            k = k + 1
            list(k,1) = min(R_start + i*delta_cgto_grid,R_end)
         enddo

         i = 0
         do p=1,cgto_A%number_of_primitives
            do q=1,cgto_B%number_of_primitives
               i = i + 1
               prod_alp = cgto_A%exponents(p)+cgto_B%exponents(q) !exponent of the product GTO
               prod_P(1:3) = (cgto_A%exponents(p)*cgto_A%center(1:3)+cgto_B%exponents(q)*cgto_B%center(1:3))/prod_alp !center of the product GTO
               RP = sqrt(dot_product(prod_P(1:3),prod_P(1:3)))
               if (RP .le. F1MACH(4,cfp_dummy)) RP = 0.0_cfp !regard this pair as sitting on CMS
               k = k + 1
               list(k,1) = RP
            enddo !q
         enddo !p

         call sort_float(k,1,list)

         n_ranges = 0
         do i=2,k
            if (list(i,1) > list(i-1,1)) then
               n_ranges = n_ranges + 1
               range_start_end(1,n_ranges) = list(i-1,1)
               range_start_end(2,n_ranges) = list(i  ,1)
               !print *,'range pair',range_start_end(1:2,n_ranges)
            endif
         enddo !i

         !Total number of quadrature points in the interval [R_start,R_end].
         n_total_points = n_ranges*n_points

         if (allocated(r_points)) deallocate(r_points)
         if (allocated(weights)) deallocate(weights)
         allocate(r_points(n_total_points),weights(n_total_points),tmp_r(n_total_points),tmp_w(n_total_points),stat=err)
         if (err /= 0) call xermsg('grid_mod','radial_grid_CGTO_pair','Memory allocation 2 failed.',err,1)

         !Construct quadratures for the individual ranges:
         cnt = 0
         do i=1,n_ranges
           !Prepare the quadrature grid for the interval [range_start,range_end] expanding the canonical G-L quadrature.
           call gl_expand_A_B(x,w,(n_points - 1)/2,tmp_r,tmp_w,range_start_end(1,i),range_start_end(2,i))
           r_points(cnt+1:cnt+n_points) = tmp_r(1:n_points)
           weights(cnt+1:cnt+n_points) = tmp_w(1:n_points)
           cnt = cnt + n_points
           !print *,'gto',range_start,range_end
         enddo
!
         if (cnt .ne. n_total_points) then
            call xermsg('grid_mod','radial_grid_CGTO_pair','Error constructing the radial grid.',3,1)
         endif

  end subroutine radial_grid_CGTO_pair

end module grid_gbl
