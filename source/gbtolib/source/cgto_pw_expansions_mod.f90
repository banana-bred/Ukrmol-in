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
module cgto_pw_expansions_gbl
  use general_quadrature_gbl
  use basis_data_generic_gbl
  use bspline_grid_gbl
  use common_obj_gbl, only: nucleus_type, resize_array, resize_copy_2d_array
  use coupling_obj_gbl
  use precisn_gbl
  use const_gbl, only: mib, line_len, stdout, level1, level2, level3
  use special_functions_gbl, only: real_harmonics_obj
  use utils_gbl, only: xermsg
  implicit none
  
  private

  public CGTO_shell_pw_expansion_obj, CGTO_shell_pair_pw_expansion_obj, init_CGTO_pw_expansions_mod, cpl
  public precalculate_Xlm_for_nuclei, pw_expansion_obj, calculate_lambda_couplings!, legendre_grid_r1_r2_obj

!ONLY FOR DEBUGING
  integer, parameter :: n_alp = 10 !11
  integer, parameter :: dbg_l = 0
  real(kind=cfp), parameter :: dbg_alp(1:n_alp)= (/0.127610_cfp,0.101448_cfp,0.181731_cfp,0.166052_cfp,0.153288_cfp,&
                0.142773_cfp,0.134058_cfp,0.126815_cfp,0.120778_cfp,0.115690_cfp/) !d-type no tails
!10 d-type
!(/0.127610,0.101448,0.081731,0.066052,0.053288,0.042773,0.034058,0.026815,0.020778,0.015690/) !10 d-type
!(/0.160912_cfp,0.121121_cfp,0.092329_cfp,0.070542_cfp,0.053752_cfp,0.040714_cfp,0.030575_cfp,0.022709_cfp,0.016634_cfp,0.011967_cfp,0.008375_cfp/) !11 s-type
!(/2.160912_cfp,2.121121_cfp,2.092329_cfp,2.070542_cfp,2.053752_cfp,2.040714_cfp,2.030575_cfp,2.022709_cfp,2.016634_cfp,2.011967_cfp,2.008375_cfp/) !no tail subtr. necessary
  type(CGTO_shell_data_obj) :: dbg_cgto(n_alp)
  public dbg_cgto, init_dbg

  type :: pw_expansion_obj
     real(kind=cfp), allocatable :: r_points(:), weights(:)
     integer :: n_total_points = -1
     logical :: initialized = .false.
  contains
     procedure :: eval_regular_grid
     procedure :: eval_exponential_grid
     procedure :: assign_grid
  end type

  type, extends(pw_expansion_obj) :: CGTO_shell_pw_expansion_obj
     type(CGTO_shell_data_obj) :: cgto_shell
     real(kind=cfp), allocatable :: angular_integrals(:,:), angular_integrals_at_knots(:,:) !(radial points,CGTO m,lm)
     real(kind=cfp), allocatable :: gaunt_angular_integrals(:,:), gaunt_angular_integrals_at_knots(:,:), Y_lm_mixed(:,:,:)
     integer, allocatable :: Y_lm_non_neg_indices(:,:,:) !Used for Y_lm_mixed to map the m_cgto,lm,lbmb indices onto one linear index.
     integer, allocatable :: non_neg_indices_l(:,:), non_neg_indices_l_lp(:,:,:) !Used to map the m_cgto,lm,lbmb indices for gaunt_angular_integrals into one linear index.
     integer, allocatable :: non_neg_indices_l_at_knots(:,:), non_neg_indices_l_lp_at_knots(:,:,:) !Used to map the m_cgto,lm,lbmb indices for gaunt_angular_integrals into one linear index.
     real(kind=cfp), allocatable :: at_lebedev_points(:,:,:)
     character(len=line_len) :: Y_lm_path = ''
     logical :: Y_lm_on_disk = .false.
     integer :: Y_lm_unit = -1, Y_lm_dim_on_disk(1:4) = -1
     integer, allocatable :: Y_lm_disk_offset(:,:)
     real(kind=cfp) :: Y_lm_size_mib = 0.0_cfp
     real(kind=cfp), allocatable :: NAI_X_lm_projections(:,:,:), Xlm(:), transl_cfs(:,:,:)
     integer :: cgto_shell_index = -1
  contains
     procedure :: init_CGTO_shell_pw_expansion
     !procedure :: eval_CGTO_radial_grid
     procedure :: eval_CGTO_shell_pw_expansion
     procedure :: eval_CGTO_single_projection_shell_pw_expansion
     procedure :: eval_BTO_CGTO_Y_lm
     procedure :: eval_BTO_CGTO_Y_lm_full
     procedure :: eval_BTO_CGTO_Y_lm_inside
     procedure :: eval_BTO_CGTO_Y_lm_endpts
     procedure :: eval_BTO_CGTO_Y_lm_outside
     procedure :: eval_NAI_X_lm_projections
     !procedure :: expand_pw_in_bsplines
     procedure :: eval_at_lebedev_points
     procedure :: final => final_CGTO_shell_pw_expansion_obj
     procedure :: write => write_CGTO_shell_pw_expansion_obj
     procedure :: write_Y_lm_to_file
     procedure :: read_Y_lm_from_file
  end type CGTO_shell_pw_expansion_obj

  type, extends(pw_expansion_obj) :: CGTO_shell_pair_pw_expansion_obj
     type(CGTO_shell_data_obj) :: cgto_shell_A, cgto_shell_B
     integer :: cgto_shell_A_index = -1, cgto_shell_B_index = -1
     real(kind=cfp), allocatable :: angular_integrals(:,:,:) !(radial points,CGTO m,lm)
     real(kind=cfp), allocatable :: radial_lm_BB_GG(:,:,:)
     logical, allocatable :: neglect_m_lm(:,:)
     !Needed for the Lebedev method:
     real(kind=cfp), allocatable :: lebedev_points(:,:), coulomb_integrals(:,:,:)
     integer :: n_leb_points = -1
     logical :: order_ab
  contains
     procedure :: init_CGTO_shell_pair_pw_expansion
     procedure :: eval_CGTO_pair_radial_grid
     procedure :: eval_CGTO_shell_pair_pw_expansion
     procedure :: eval_radial_GG_BB
     procedure :: eval_coulomb_integrals
     procedure :: eval_damped_dipole_integrals
     procedure :: eval_GGG_integrals
     procedure :: eval_GBG_integrals
     procedure :: storage_occupied => storage_occupied_CGTO_shell_pair_pw_expansion
     procedure :: final => final_CGTO_shell_pair_pw_expansion
  end type CGTO_shell_pair_pw_expansion_obj

  !> Set to .true. following a call to init_CGTO_pw_expansions_mod.
  logical :: module_initialized = .false.

  !> Largest L in the pw expansions of the CGTOs that is required. This is set by init_CGTO_pw_expansions_mod and must not be changed by any other routine.
  integer :: max_l_pw = -1

  !> Used to evaluate various coupling coefficients.
  type(couplings_type) :: cpl

  !> Used to evaluate the real spherical and solid harmonics. PRIVATE FOR ALL THREADS.
  type(real_harmonics_obj) :: real_harmonics

  !> These are used to keep the largest L for which the coefficients a,b,c and as,bs,cs for the real spherical and real solid harmonic have been precalculated.
  integer :: max_l = -1, max_ls = -1
  !> Coefficients used for calculation of real spherical and real solid harmonics using the routines resh, solh.
  real(kind=cfp), allocatable :: a(:), b(:), c(:), as(:), bs(:), cs(:)

  !> Values of the real spherical harmonics at the positions of the nuclei
!  real(kind=cfp), allocatable :: Xlm_nuclei(:)
!  !> Stride in array Xlm_nuclei.
!  integer :: n_Xlm_nuclei = 0

  !> Object that is used to evaluate the product of a pair of CGTOs and a real spherical harmonic at arbitrary points in space. This is used when integrating this function adaptively over sphere.
  !> Used only for debugging of the pw expansion routines.
  type, extends(function_2d) :: Xlm_x_pair_cgto_surface
     !> Contains the data on the CGTO.
     type(CGTO_shell_data_obj) :: cgto_A, cgto_B
     !>
     integer :: cgto_A_m = 0, cgto_B_m = 0
     !> radial distance for which to evaluate.
     real(kind=cfp) :: r = 0.0_cfp
     !> Angular numbers of the real spherical harmonic centered on CSM multiplying the GTO.
     integer :: l = 0, m = 0
  contains
     procedure :: eval => eval_Xlm_x_pair_cgto_surface
  end type Xlm_x_pair_cgto_surface

contains

  subroutine init_CGTO_pw_expansions_mod(inp_max_l_pw,max_l_cgto)
     implicit none
     integer, intent(in) :: inp_max_l_pw, max_l_cgto

     integer :: err, ind, max_l_aux

        max_l_pw = inp_max_l_pw
        max_l_aux = max_l_pw+2*max_l_cgto

        !The two routines below must be called before the OpenMP region to ensure thread safety of resh, solh.
        call calc_resh_coefficients(max_l_aux)
        call calc_solh_coefficients(max_l_aux)

        !Precalculate all couplings needed to evaluate the overlap-type integrals and the nuclear attraction integrals.
        call cpl%prec_G_cf(max_l_aux)

        write(level2,'(/,10X,"cgto_pw_expansions_mod initialized with max_l_pw = ",i3)') max_l_pw

        module_initialized = .true.

  end subroutine init_CGTO_pw_expansions_mod

  subroutine init_CGTO_shell_pw_expansion(this,cgto_shell,cgto_shell_index)
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this
     type(CGTO_shell_data_obj) :: cgto_shell
     integer, intent(in) :: cgto_shell_index

     real(kind=cfp) :: RA

        if (.not. module_initialized) then
            call xermsg ('cgto_pw_expansions_mod', 'init_CGTO_shell_pw_expansion', &
                         'The module has not been initialized: run init_CGTO_pw_expansions_mod first.', 1, 1)
        end if

        this%cgto_shell = cgto_shell
        this%cgto_shell_index = cgto_shell_index

        this%Y_lm_on_disk = .false.
        this%Y_lm_size_mib = 0.0_cfp

        RA = sqrt(dot_product(this%cgto_shell%center,this%cgto_shell%center))

        !Real spherical harmonics for the CGTO center: result in the module array this%Xlm
        call real_harmonics%precalculate_Xlm_for_CGTO_center(this%cgto_shell%center,this%cgto_shell%l,this%Xlm)

        !Precalculate the coefficients in the translation formula for the solid harmonics: this requires this%Xlm
        call precalculate_solh_translation_coeffs(this%cgto_shell%l,RA,this%Xlm,this%transl_cfs)

        this%initialized = .true.

  end subroutine init_CGTO_shell_pw_expansion

  subroutine final_CGTO_shell_pw_expansion_obj(this)
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this

        if (allocated(this%angular_integrals)) deallocate(this%angular_integrals)
        if (allocated(this%r_points)) deallocate(this%r_points)
        if (allocated(this%weights)) deallocate(this%weights)
        if (allocated(this%angular_integrals)) deallocate(this%angular_integrals)
        if (allocated(this%angular_integrals_at_knots)) deallocate(this%angular_integrals_at_knots)
        if (allocated(this%gaunt_angular_integrals)) deallocate(this%gaunt_angular_integrals)
        if (allocated(this%gaunt_angular_integrals_at_knots)) deallocate(this%gaunt_angular_integrals_at_knots)
        if (allocated(this%Y_lm_mixed)) deallocate(this%Y_lm_mixed)
        if (allocated(this%Y_lm_non_neg_indices)) deallocate(this%Y_lm_non_neg_indices)
        if (allocated(this%non_neg_indices_l)) deallocate(this%non_neg_indices_l)
        if (allocated(this%non_neg_indices_l_lp)) deallocate(this%non_neg_indices_l_lp)
        if (allocated(this%non_neg_indices_l_at_knots)) deallocate(this%non_neg_indices_l_at_knots)
        if (allocated(this%non_neg_indices_l_lp_at_knots)) deallocate(this%non_neg_indices_l_lp_at_knots)
        if (allocated(this%at_lebedev_points)) deallocate(this%at_lebedev_points)
        if (allocated(this%Y_lm_disk_offset)) deallocate(this%Y_lm_disk_offset)

        this%initialized = .false.
        this%n_total_points = -1
        this%Y_lm_on_disk = .false.
        this%Y_lm_size_mib = 0.0_cfp

  end subroutine final_CGTO_shell_pw_expansion_obj

  subroutine init_CGTO_shell_pair_pw_expansion(this,cgto_shell_A,cgto_shell_A_index,cgto_shell_B,cgto_shell_B_index)
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     type(CGTO_shell_data_obj) :: cgto_shell_A, cgto_shell_B
     integer, intent(in) :: cgto_shell_A_index, cgto_shell_B_index

     integer :: max_l_aux

        if (.not. module_initialized) then
            call xermsg ('cgto_pw_expansions_mod', 'init_CGTO_shell_pair_pw_expansion', &
                         'The module has not been initialized: run init_CGTO_pw_expansions_mod first.', 1, 1)
        end if

        this%cgto_shell_A = cgto_shell_A
        this%cgto_shell_A_index = cgto_shell_A_index

        this%cgto_shell_B = cgto_shell_B
        this%cgto_shell_B_index = cgto_shell_B_index

        max_l_aux = max_l_pw+max(cgto_shell_B%l,cgto_shell_A%l)

        if (max_l_aux .gt. cpl%L) then
            call xermsg ('cgto_pw_expansions_mod', 'init_CGTO_shell_pair_pw_expansion', &
                         'Requesting higher max_l than precalculated.', 1, 1)
        end if

        this%initialized = .true.

  end subroutine init_CGTO_shell_pair_pw_expansion

  subroutine final_CGTO_shell_pair_pw_expansion(this)
     class(CGTO_shell_pair_pw_expansion_obj) :: this

        this%cgto_shell_A_index = -1
        this%cgto_shell_B_index = -1
        this%n_leb_points = -1

        if (allocated(this%angular_integrals)) deallocate(this%angular_integrals)
        if (allocated(this%radial_lm_BB_GG)) deallocate(this%radial_lm_BB_GG)
        if (allocated(this%neglect_m_lm)) deallocate(this%neglect_m_lm)
        if (allocated(this%lebedev_points)) deallocate(this%lebedev_points)
        if (allocated(this%coulomb_integrals)) deallocate(this%coulomb_integrals)

  end subroutine final_CGTO_shell_pair_pw_expansion

  function storage_occupied_CGTO_shell_pair_pw_expansion(this)
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     real(kind=cfp) :: storage_occupied_CGTO_shell_pair_pw_expansion

        storage_occupied_CGTO_shell_pair_pw_expansion = 0.0_cfp

        if (allocated(this%angular_integrals)) then
           storage_occupied_CGTO_shell_pair_pw_expansion = storage_occupied_CGTO_shell_pair_pw_expansion &
           + size(this%angular_integrals,1)*size(this%angular_integrals,2)*size(this%angular_integrals,3) &
           * cfp_bytes / (Mib * 1.0_cfp)
        endif

        if (allocated(this%radial_lm_BB_GG)) then
           storage_occupied_CGTO_shell_pair_pw_expansion = storage_occupied_CGTO_shell_pair_pw_expansion &
           + size(this%radial_lm_BB_GG,1)*size(this%radial_lm_BB_GG,2)*size(this%radial_lm_BB_GG,3) &
           * cfp_bytes / (Mib * 1.0_cfp)
        endif

        if (allocated(this%lebedev_points)) then
           storage_occupied_CGTO_shell_pair_pw_expansion = storage_occupied_CGTO_shell_pair_pw_expansion &
           + size(this%lebedev_points,1)*size(this%lebedev_points,2) &
           * cfp_bytes / (Mib * 1.0_cfp)
        endif

        if (allocated(this%coulomb_integrals)) then
           storage_occupied_CGTO_shell_pair_pw_expansion = storage_occupied_CGTO_shell_pair_pw_expansion &
           + size(this%coulomb_integrals,1)*size(this%coulomb_integrals,2)*size(this%coulomb_integrals,3) &
           * cfp_bytes / (Mib * 1.0_cfp)
        endif

        if (allocated(this%neglect_m_lm)) then
           storage_occupied_CGTO_shell_pair_pw_expansion = storage_occupied_CGTO_shell_pair_pw_expansion &
           + size(this%neglect_m_lm,1)*size(this%neglect_m_lm,2) &
           * storage_unit_int / (Mib * 1.0_cfp)
        endif

  end function storage_occupied_CGTO_shell_pair_pw_expansion
 
  !> USED ONLY FOR DEBUGGING. 
  subroutine init_dbg(bspline_grid)
     use gto_routines_gbl, only: cms_gto_norm
     implicit none
     type(bspline_grid_obj) :: bspline_grid
     real(kind=cfp) :: dbg_norm
     integer :: i

        do i=1,n_alp
           call dbg_cgto(i)%make_space(1)
           dbg_cgto(i)%l = dbg_l
           !todo DBG
           !if (i .eq. n_alp) dbg_cgto(i)%l = dbg_l+1
           dbg_cgto(i)%center = 0.0_cfp
           dbg_cgto(i)%number_of_primitives = 1
           dbg_cgto(i)%contractions(1) = 1.0_cfp
           dbg_cgto(i)%non_zero_at_boundary = .false.
           dbg_cgto(i)%number_of_functions = 2*dbg_cgto(i)%l+1
           dbg_cgto(i)%exponents(1) = dbg_alp(i)
           call dbg_cgto(i)%normalize
           dbg_norm = cms_gto_norm(bspline_grid%B,dbg_cgto(i)%l,dbg_cgto(i)%number_of_primitives,&
                            dbg_cgto(i)%exponents,dbg_cgto(i)%contractions,dbg_cgto(i)%norm,dbg_cgto(i)%norms)
           dbg_cgto(i)%norm = dbg_cgto(i)%norm*dbg_norm
        enddo !i

  end subroutine init_dbg

  subroutine eval_regular_grid(this,A,B,delta_r)
     use general_quadrature_gbl, only: n_7, x_7, w_7, gl_expand_A_B
     use const_gbl, only: epsabs
     implicit none
     class(pw_expansion_obj) :: this
     real(kind=cfp), intent(in) :: A, B, delta_r

     integer :: i, n, err, n_points, cnt
     real(kind=cfp) :: r, r_prev

        if (allocated(this%r_points)) deallocate(this%r_points)
        if (allocated(this%weights)) deallocate(this%weights)

        if (A < 0.0_cfp .or. B .le. 0.0_cfp .or. B .le. A .or. delta_r .le. 0.0_cfp .or. delta_r > B-A) then
           print *,A,B,delta_r
           call xermsg('pw_expansion_obj','eval_regular_grid','On input at least one of A,B,delta_r were invalid.',1,1)
        endif

        !Count the number of radial points
        n_points = 2*n_7+1
        r_prev = A
        n = 0
        do 
           r = min(r_prev+delta_r,B)
           if (r-r_prev .ge. epsabs) n = n + n_points
           if (r .eq. B) exit
           r_prev = r
        enddo
        this%n_total_points = n

        allocate(this%r_points(n),this%weights(n),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_regular_grid','Memory allocation failed.',err,1)

        r_prev = A
        cnt = 0
        do
           r = min(r_prev+delta_r,B)
           if (r-r_prev .ge. epsabs) then
              call gl_expand_A_B(x_7,w_7,n_7,this%r_points(cnt+1:cnt+n_points),this%weights(cnt+1:cnt+n_points),r_prev,r)
              cnt = cnt + n_points
           endif
           if (r .eq. B) exit
           r_prev = r
        enddo

  end subroutine eval_regular_grid

  !> Exponential grid on [0,1] with points accumulating towards 1.
  subroutine eval_exponential_grid(this,alp,n,order_01_or_10)
     use general_quadrature_gbl, only: n_7, x_7, w_7, gl_expand_A_B
     implicit none
     class(pw_expansion_obj) :: this
     real(kind=cfp), intent(in) :: alp
     integer, intent(in) :: n
     logical, intent(in) :: order_01_or_10

     real(kind=cfp), allocatable :: intervals(:)
     real(kind=cfp) :: A, B
     integer :: i, n_points, err, cnt

        if (allocated(this%r_points)) deallocate(this%r_points)
        if (allocated(this%weights)) deallocate(this%weights)

        n_points = 2*n_7+1
        this%n_total_points = (n-1)*n_points

        allocate(this%r_points(this%n_total_points),stat=err)
        if (err .ne. 0) call xermsg('pw_expansion_obj','eval_exponential_grid','Memory allocation 1 failed.',err,1)

        allocate(this%weights(this%n_total_points),stat=err)
        if (err .ne. 0) call xermsg('pw_expansion_obj','eval_exponential_grid','Memory allocation 2 failed.',err,1)
        
        allocate(intervals(n),stat=err)
        if (err .ne. 0) call xermsg('pw_expansion_obj','eval_exponential_grid','Memory allocation 3 failed.',err,1)

        if (order_01_or_10) then
           do i=1,n
              intervals(i) = 1.0_cfp - (exp(alp*((i-1.0_cfp)/(n-1.0_cfp))) - 1.0_cfp)/(exp(alp) - 1.0_cfp)
           enddo !i
        else
           do i=1,n
              intervals(i) = (exp(alp*((i-1.0_cfp)/(n-1.0_cfp))) - 1.0_cfp)/(exp(alp) - 1.0_cfp)
           enddo !i
        endif

        cnt = 0
        do i=2,n
           A = intervals(i-1)
           B = intervals(i)
           call gl_expand_A_B(x_7,w_7,n_7,this%r_points(cnt+1:cnt+n_points),this%weights(cnt+1:cnt+n_points),A,B)
           cnt = cnt + n_points
        enddo !i

  end subroutine eval_exponential_grid

  subroutine assign_grid(this,r_points,weights)
     implicit none
     class(pw_expansion_obj) :: this
     real(kind=cfp), intent(in) :: r_points(:), weights(:)

     integer :: err, n

        if (.not.(this%initialized)) call xermsg('pw_expansion_obj','assign_grid','Object not initialized.',1,1)

        n = size(r_points)
        if (n .ne. size(weights)) call xermsg('pw_expansion_obj','assign_grid','r1 grid incompatible with weights.',2,1)
        this%n_total_points = n

        if (allocated(this%r_points)) deallocate(this%r_points)
        if (allocated(this%weights)) deallocate(this%weights)

        allocate(this%r_points,source=r_points,stat=err)
        if (err .ne. 0) call xermsg('pw_expansion_obj','assign_grid','Memory allocation 1 failed.',err,1)

        allocate(this%weights,source=weights,stat=err)
        if (err .ne. 0) call xermsg('pw_expansion_obj','assign_grid','Memory allocation 2 failed.',err,1)

  end subroutine assign_grid

  !> Don't compute the double angular projections and the projections for the knot points.
  subroutine eval_CGTO_single_projection_shell_pw_expansion(this,threshold,lpp)
     use const_gbl, only: epsrel, epsabs
     use omp_lib
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this
     integer, intent(in) :: lpp
     real(kind=cfp), intent(in) :: threshold

     real(kind=wp) :: start_t, end_t
     real(kind=cfp), allocatable :: dummy_proj(:,:)
     integer, allocatable :: dummy_ind(:,:,:)

        write(level3,'("--------->","CGTO_shell_pw_expansion_obj:eval_CGTO_single_projection_shell_pw_expansion")')

        start_t = omp_get_wtime()

        if (.not. this % initialized) then
           call xermsg('cgto_pw_expansions_mod','eval_CGTO_single_projection_shell_pw_expansion','Object not initialized.',1,1)
        endif
!
        if (.not. allocated(this%r_points)) then
           call xermsg('cgto_pw_expansions_mod','eval_CGTO_single_projection_shell_pw_expansion','r1 grid not allocated.',2,1)
        endif
        if (.not. allocated(this%weights)) then
           call xermsg('cgto_pw_expansions_mod','eval_CGTO_single_projection_shell_pw_expansion', &
                       'weights for r1 grid not allocated.',3,1)
        endif
        if (size(this%r_points) /= size(this%weights)) then
           call xermsg('cgto_pw_expansions_mod','eval_CGTO_single_projection_shell_pw_expansion', &
                       'r1 grid incompatible with weights.',4,1)
        endif
!
!------ Calculate the projections of the CGTOs on the real spherical harmonics (i.e. PW-expansion) for all radial grid points and maximum L = lpp
        call omp_calculate_CGTO_pw_coefficients_analytic(threshold, &
                                                         lpp, &
                                                         0,&!-1
                                                         this%cgto_shell, &
                                                         this%r_points, &
                                                         this%angular_integrals, &
                                                         this%non_neg_indices_l, &
                                                         dummy_proj, &
                                                         dummy_ind)

        end_t = omp_get_wtime()

        write(level3,'("<---------",&
                      &"CGTO_shell_pw_expansion_obj:eval_CGTO_single_projection_shell_pw_expansion and took [s]: ",f8.3)') &
                      end_t-start_t

  end subroutine eval_CGTO_single_projection_shell_pw_expansion

  subroutine eval_CGTO_shell_pw_expansion(this,knots,max_bspline_l,max_prop_l,max_l_legendre)
     use const_gbl, only: epsrel, epsabs
     use omp_lib
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this
     integer, intent(in) :: max_bspline_l,max_prop_l,max_l_legendre
     real(kind=cfp), allocatable :: knots(:)

     real(kind=cfp) :: threshold
     real(kind=wp) :: start_t, end_t

        write(level3,'("--------->","CGTO_shell_pw_expansion_obj:eval_CGTO_shell_pw_expansion")')

        start_t = omp_get_wtime()

        if (.not. this % initialized) then
            call xermsg ('cgto_pw_expansions_mod', 'eval_CGTO_shell_pw_expansion', 'Object not initialized.', 1, 1)
        end if
!
        if (.not. allocated(this % r_points)) then
            call xermsg ('cgto_pw_expansions_mod', 'eval_CGTO_shell_pw_expansion', 'r1 grid not allocated.', 2, 1)
        end if
        if (.not. allocated(this % weights)) then
            call xermsg ('cgto_pw_expansions_mod', 'eval_CGTO_shell_pw_expansion', 'weights for r1 grid not allocated.', 3, 1)
        end if
        if (size(this % r_points) /= size(this % weights)) then
            call xermsg ('cgto_pw_expansions_mod', 'eval_CGTO_shell_pw_expansion', 'r1 grid incompatible with weights.', 4, 1)
        end if

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs
!
!------ Calculate the projections of the CGTOs on the real spherical harmonics (i.e. PW-expansion) for all radial grid points.
        call omp_calculate_CGTO_pw_coefficients_analytic (threshold, max(max_l_legendre, max_prop_l), max_bspline_l, &
                this % cgto_shell, this % r_points, this % angular_integrals, this % non_neg_indices_l, &
                this % gaunt_angular_integrals, this % non_neg_indices_l_lp)
!
!------ Calculate the projections of the CGTOs on the real spherical harmonics (i.e. PW-expansion) for the knot points: this is for the Bloch operator
        call omp_calculate_CGTO_pw_coefficients_analytic (threshold, max_bspline_l, 0, this % cgto_shell, &
                knots, this % angular_integrals_at_knots, this % non_neg_indices_l_at_knots, &
                this % gaunt_angular_integrals_at_knots, this % non_neg_indices_l_lp_at_knots)

        end_t = omp_get_wtime()

        write(level3,'("<---------","CGTO_shell_pw_expansion_obj:eval_CGTO_shell_pw_expansion and took [s]: ",f8.3)') end_t-start_t

  end subroutine eval_CGTO_shell_pw_expansion

  subroutine write_CGTO_shell_pw_expansion_obj(this)
     use const_gbl, only: line_len
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this

     integer :: lu, p, lm, cgto_m, i
     character(len=line_len) :: path, a

        if (.not. this % initialized) then
            call xermsg ('cgto_pw_expansions_mod', 'write_CGTO_shell_pw_expansion_obj', 'Object not initialized.', 1, 1)
        end if

        write(a,'(i2)') this%cgto_shell_index
        path='cgto_pw.'//trim(adjustl(a))
        open(file=path,newunit=lu,form='formatted')

        do lm=1,size(this%non_neg_indices_l,2)
           do cgto_m=1,size(this%non_neg_indices_l,1)
              p = this%non_neg_indices_l(cgto_m,lm)
              if (p .eq. 0) cycle
              write(lu,'("angular part: ",i0,1x,i0)') cgto_m,lm
              do i=1,size(this%angular_integrals,1)
                 write(lu,'(e25.15)') this%angular_integrals(i,p)
              enddo !i
           enddo !cgto_m
        enddo !lm

        close(lu)

  end subroutine write_CGTO_shell_pw_expansion_obj


  !> \brief   Integrate two-electron mixed integral over one variable
  !> \authors Z Masin, J Benda
  !> \date    2017 - 2025
  !>
  !> This subroutine evaluates the integral
  !> \f[
  !>     \mathcal{Y}_{lm}^{ab}(r_1) = \frac{4\pi}{2l + 1} r_1 \int_{t_b}^{t_{b+o}} B_{i_b}(r_2) \frac{r_<^l}{r_>^{l+1}} \int_{4\pi}
  !>     X_{l_b m_b}(\hat{\mathbf{r}}_2) G_a(\mathbf{r}_2) X_{lm}(\hat{\mathbf{r}}_2) \mathrm{d}\Omega_2 r_2 \mathrm{d}r_2
  !> \f]
  !> for one CGTO \f$ a \f$, all B-splines \f$ b \equiv (i_b, l_b, m_b) \f$, and all Legendre expansion quantum numbers
  !> \f$ l, m \f$. The angular integral in performed first by expanding the CGTO in partial waves and performing the angular
  !> integration at all radii \f$ r_2 \f$. The remaining radial integral over \f$ r_2 \f$ is performed in the following way:
  !>
  !>  - Phase I:   If \f$ r_1 \f$ is within the support of the radial B-spline \f$ b \f$, straightforward numerical quadrature
  !>               is performed. This phase is implemented in subroutine eval_BTO_CGTO_Y_lm_phase1.
  !>  - Phase II:  Auxiliary integrals at the starting and ending knots of individual B-splines are calculated. These are
  !>               not needed by the Gauss-Legendre quadrature, but they are used later in Phase III. This phase is implemented
  !>               in subroutine eval_BTO_CGTO_Y_lm_phase2.
  !>  - Phase III: The remaining integrals are calculated for cases where \f$ r_1 \f$ lies outside of the support of the given
  !>               radial B-spline. Then the \f$ r_1 \f$-dependence factorizes out of the integral and so the values can be obtained
  !>               simply by scaling the auxiliary values precomputed in Phase II:
  !>               \f[
  !>                     \mathcal{Y}_{lm}^{ab}(r_1 < t_{i_b}) = (r_1/t_{i_b})^l \mathcal{Y}_{lm}^{ab}(t_{i_b})
  !>               \f]
  !>               \f[
  !>                     \mathcal{Y}_{lm}^{ab}(r_1 > t_{i_b+o}) = (t_{i_b+o}/r_1)^{l+1} \mathcal{Y}_{lm}^{ab}(t_{i_b+o})
  !>               \f]
  !>               This phase is implemented in subroutine eval_BTO_CGTO_Y_lm_phase3.
  !>
  subroutine eval_BTO_CGTO_Y_lm(this,grid_r1_r2)
     use const_gbl, only: epsrel, epsabs
     use omp_lib,   only: omp_get_wtime
     use grid_gbl,  only: grid_r1_r2_obj
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this
     class(grid_r1_r2_obj) :: grid_r1_r2

     logical, parameter :: factor_Y_lm = .true.  ! take advantage of factorization for non-overlapping B-splines

     integer :: err, k, p, q, ind
     real(kind=cfp) :: threshold,distance,R_min,R_max,r1,r2
     real(kind=cfp), allocatable :: bto_leg_part(:), f_l_w1_w2(:), inv_r1(:), inv_r2(:), Y_lm_begin(:,:), Y_lm_end(:,:)
     real(kind=wp) :: start_t, end_t

        write(level3,'("--------->","CGTO_shell_pw_expansion_obj:eval_BTO_CGTO_Y_lm")')

        start_t = omp_get_wtime()

        if (.not.(this%initialized)) call xermsg('cgto_pw_expansions_mod','eval_BTO_CGTO_Y_lm','Object not initialized.',1,1)

        call this%assign_grid(grid_r1_r2%r2,grid_r1_r2%w2)

        !Calculate the partial wave expansion of the CGTO on the r2 grid
        call this % eval_CGTO_shell_pw_expansion(grid_r1_r2 % bspline_grid % knots, &
                                                 grid_r1_r2 % max_bspline_l, &
                                                 grid_r1_r2 % max_prop_l, &
                                                 grid_r1_r2 % max_l_legendre)

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs

        call this%cgto_shell%estimate_shell_radius(1,this%cgto_shell%number_of_primitives,threshold,distance,R_min,R_max)
        do ind=grid_r1_r2%first_bspline_index,grid_r1_r2%bspline_grid%n
           p = grid_r1_r2%bspline_start_end_r2(1,ind)
           q = grid_r1_r2%bspline_start_end_r2(2,ind)
           r1 = grid_r1_r2%r2(p)
           r2 = grid_r1_r2%r2(q)
           !if (r1 < R_min .and. r2 < R_max .or. r1 > R_max .and. r2 > R_max) then
           !   print *,'bspline',ind,'out of range'
           !endif
        enddo

        write(level1,'("CGTO is non-negligible in the radial range: ",2e25.15)') R_min,R_max

        !Assign the non-negligible indices from this%gaunt_angular_integrals valid for the r2 grid into Y_lm_non_neg_indices
        if (allocated(this%Y_lm_non_neg_indices)) deallocate(this%Y_lm_non_neg_indices)

        allocate(this % Y_lm_non_neg_indices(size(this % non_neg_indices_l_lp, 1), &
                                             size(this % non_neg_indices_l_lp, 2), &
                                             size(this % non_neg_indices_l_lp, 3)), stat = err)

        this%Y_lm_non_neg_indices = this%non_neg_indices_l_lp

        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_BTO_CGTO_Y_lm','Memory allocation 1 failed.',err,1)

        if (allocated(this%Y_lm_mixed)) deallocate(this%Y_lm_mixed)
        p = size(this%gaunt_angular_integrals,2) !maxval(this%Y_lm_non_neg_indices)
        this % Y_lm_size_mib = grid_r1_r2 % n1_total_points * p * (grid_r1_r2 % bspline_grid % n &
                                            - grid_r1_r2 % first_bspline_index + 1) * cfp_bytes / (Mib * 1.0_cfp)
        write(level1,'("Memory (MiB) required for the Y_lm function: ",f25.10)') this%Y_lm_size_mib
        allocate (this%Y_lm_mixed(grid_r1_r2%n1_total_points, p, grid_r1_r2%first_bspline_index:grid_r1_r2%bspline_grid%n), &
                  Y_lm_begin(p, grid_r1_r2%bspline_grid%n), &
                  Y_lm_end(p, grid_r1_r2%bspline_grid%n), &
                  inv_r1(grid_r1_r2%n1_total_points), &
                  inv_r2(grid_r1_r2%n2_total_points), &
                  stat = err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_BTO_CGTO_Y_lm','Memory allocation 2 failed.',err,1)
        this%Y_lm_mixed = 0.0_cfp

        !Precalculate inverses of the r1 and r2 grid coordinates to replace the
        !division in the calculation of the radial part of the Legendre
        !resolution below.
        do k=1,grid_r1_r2%n1_total_points !r1
           inv_r1(k) = 1.0_cfp/grid_r1_r2%r1(k)
        enddo !k

        do k=1,grid_r1_r2%n2_total_points !r2
           inv_r2(k) = 1.0_cfp/grid_r1_r2%r2(k)
        enddo !k

        !$OMP PARALLEL DEFAULT(NONE) PRIVATE(err,bto_leg_part,f_l_w1_w2) SHARED(this,grid_r1_r2,inv_r1,inv_r2,Y_lm_begin,Y_lm_end)
        allocate(bto_leg_part(grid_r1_r2%n2_total_points),f_l_w1_w2(grid_r1_r2%n2_total_points),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_BTO_CGTO_Y_lm','Memory allocation 3 failed.',err,1)
        !$OMP BARRIER
        if (factor_Y_lm) then
           call this % eval_BTO_CGTO_Y_lm_inside(grid_r1_r2, bto_leg_part, f_l_w1_w2, inv_r1, inv_r2)
           call this % eval_BTO_CGTO_Y_lm_endpts(grid_r1_r2, bto_leg_part, inv_r2, Y_lm_begin, Y_lm_end)
           call this % eval_BTO_CGTO_Y_lm_outside(grid_r1_r2, Y_lm_begin, Y_lm_end)
        else
           call this % eval_BTO_CGTO_Y_lm_full(grid_r1_r2, bto_leg_part, f_l_w1_w2, inv_r1, inv_r2)
        end if
        !$OMP END PARALLEL

        this%Y_lm_on_disk = .false.

        end_t = omp_get_wtime()

        write(level1,'("Calculation of (GB|Ylm) functions took [s]: ",F0.1)') end_t - start_t
        write(level3,'("<---------","CGTO_shell_pw_expansion_obj:eval_BTO_CGTO_Y_lm")')

  end subroutine eval_BTO_CGTO_Y_lm


  !> \brief   Original calculation of Y_lm_mixed
  !> \authors Z Masin, J Benda
  !> \date    2017 - 2025
  !>
  !> This subroutine does not take advantage of factorization of Y_lm_mixed for non-overlapping B-splines.
  !> It is kept as an option for debugging.
  !>
  subroutine eval_BTO_CGTO_Y_lm_full(this, grid_r1_r2, bto_leg_part, f_l_w1_w2, inv_r1, inv_r2)

     use grid_gbl,       only: grid_r1_r2_obj
     use phys_const_gbl, only: fourpi

     class(CGTO_shell_pw_expansion_obj)   :: this
     class(grid_r1_r2_obj), intent(in)    :: grid_r1_r2
     real(kind=cfp),        intent(in)    :: inv_r1(:), inv_r2(:)
     real(kind=cfp),        intent(inout) :: bto_leg_part(:), f_l_w1_w2(:)

     real(kind=cfp) :: fac, r12, Y_lm
     integer        :: B_start, B_end, j, k, k2, l, m, lm, lbmb, cgto_m, p, n_cgto_m

        n_cgto_m = 2*this%cgto_shell%l + 1

        !$OMP DO SCHEDULE(DYNAMIC)
        do k = 1, grid_r1_r2 % n1_total_points !r1
           do l = 0, grid_r1_r2 % max_l_legendre
              fac = fourpi/(2*l+1.0_cfp)
              do k2=1,grid_r1_r2%n2_total_points !r2
                 !Calculate w1*w2*fac*r1*r<**l/r>**(l+1) from the Legendre expansion
                 if (grid_r1_r2%r2(k2) .le. grid_r1_r2%r1(k)) then
                    r12 = grid_r1_r2%r2(k2)*inv_r1(k) !/grid_r1_r2%r1(k) !to make sure ifort 16.0.1 produces correct results
                    f_l_w1_w2(k2) = grid_r1_r2%w2(k2)*grid_r1_r2%w1(k)*fac*(r12)**l
                 else
                    r12 = grid_r1_r2%r1(k)*inv_r2(k2) !/grid_r1_r2%r2(k2) !to make sure ifort 16.0.1 produces correct results
                    f_l_w1_w2(k2) = grid_r1_r2%w2(k2)*grid_r1_r2%w1(k)*fac*(r12)**(l+1)
                 endif
              enddo !k2
              do j=grid_r1_r2%first_bspline_index,grid_r1_r2%bspline_grid%n
                 B_start = grid_r1_r2%bspline_start_end_r2(1,j)
                 B_end = grid_r1_r2%bspline_start_end_r2(2,j)
                 bto_leg_part(B_start:B_end) = f_l_w1_w2(B_start:B_end) &
                                *this%r_points(B_start:B_end)*grid_r1_r2%B_vals_r2(B_start:B_end,j)
                 do m=-l,l
                    lm = l*l+l+m+1
                    do lbmb=1,(grid_r1_r2%max_bspline_l+1)**2
                       do cgto_m=1,n_cgto_m
                          p = this%Y_lm_non_neg_indices(cgto_m,lbmb,lm)
                          if (p .eq. 0) cycle
                          this%Y_lm_mixed(k,p,j) = sum(bto_leg_part(B_start:B_end)*this%gaunt_angular_integrals(B_start:B_end,p))
                       enddo !cgto_m
                    enddo !lbmb
                 enddo !m
              enddo !j
           end do !l
        end do !k
        !$OMP END DO

  end subroutine eval_BTO_CGTO_Y_lm_full


  !> \brief   Phase I of evaluation of partially integrated mixed two-electron integrals (Ylm)
  !> \authors Z Masin, J Benda
  !> \date    2017 - 2025
  !>
  !> Evaluates Y_lm_mixed for the non-factorizable case, where r1 lies inside the given B-spline's support.
  !> This subroutine is essentially a clone of `eval_BTO_CGTO_Y_lm_full`, except that evaluation of `Y_lm_mixed(k,p,j)`
  !> is not performed for radial positions `k` that are outside of the j-th B-spline support. That is left
  !> for a later call to `eval_BTO_CGTO_Y_lm_outside`.
  !>
  subroutine eval_BTO_CGTO_Y_lm_inside(this, grid_r1_r2, bto_leg_part, f_l_w1_w2, inv_r1, inv_r2)

     use grid_gbl,       only: grid_r1_r2_obj
     use phys_const_gbl, only: fourpi

     class(CGTO_shell_pw_expansion_obj)   :: this
     class(grid_r1_r2_obj), intent(in)    :: grid_r1_r2
     real(kind=cfp),        intent(in)    :: inv_r1(:), inv_r2(:)
     real(kind=cfp),        intent(inout) :: bto_leg_part(:), f_l_w1_w2(:)

     real(kind=cfp) :: fac, r12
     integer        :: B_start, B_end, j, k, k2, l, m, lm, lbmb, cgto_m, p, n_cgto_m

        n_cgto_m = 2*this%cgto_shell%l + 1

        !$OMP DO SCHEDULE(DYNAMIC)
        do k = 1, grid_r1_r2 % n1_total_points !r1
           do l = 0, grid_r1_r2 % max_l_legendre
              fac = fourpi/(2*l+1.0_cfp)
              do k2=1,grid_r1_r2%n2_total_points !r2
                 !Calculate w1*w2*fac*r1*r<**l/r>**(l+1) from the Legendre expansion
                 if (grid_r1_r2%r2(k2) .le. grid_r1_r2%r1(k)) then
                    r12 = grid_r1_r2%r2(k2)*inv_r1(k) !/grid_r1_r2%r1(k) !to make sure ifort 16.0.1 produces correct results
                    f_l_w1_w2(k2) = grid_r1_r2%w2(k2)*grid_r1_r2%w1(k)*fac*(r12)**l
                 else
                    r12 = grid_r1_r2%r1(k)*inv_r2(k2) !/grid_r1_r2%r2(k2) !to make sure ifort 16.0.1 produces correct results
                    f_l_w1_w2(k2) = grid_r1_r2%w2(k2)*grid_r1_r2%w1(k)*fac*(r12)**(l+1)
                 endif
              enddo !k2
              do j=grid_r1_r2%first_bspline_index,grid_r1_r2%bspline_grid%n
                 B_start = grid_r1_r2%bspline_start_end_r1(1,j)
                 B_end = grid_r1_r2%bspline_start_end_r1(2,j)
                 if (k < B_start) exit
                 if (B_end < k) cycle
                 B_start = grid_r1_r2%bspline_start_end_r2(1,j)
                 B_end = grid_r1_r2%bspline_start_end_r2(2,j)
                 bto_leg_part(B_start:B_end) = f_l_w1_w2(B_start:B_end) &
                                *this%r_points(B_start:B_end)*grid_r1_r2%B_vals_r2(B_start:B_end,j)
                 do m=-l,l
                    lm = l*l+l+m+1
                    do lbmb=1,(grid_r1_r2%max_bspline_l+1)**2
                       do cgto_m=1,n_cgto_m
                          p = this%Y_lm_non_neg_indices(cgto_m,lbmb,lm)
                          if (p .eq. 0) cycle
                          this%Y_lm_mixed(k,p,j) = sum(bto_leg_part(B_start:B_end)*this%gaunt_angular_integrals(B_start:B_end,p))
                       enddo !cgto_m
                    enddo !lbmb
                 enddo !m
              enddo !j
           end do !l
        end do !k
        !$OMP END DO

  end subroutine eval_BTO_CGTO_Y_lm_inside


  !> \brief   Phase II of evaluation of partially integrated mixed two-electron integrals (Ylm)
  !> \authors J Benda
  !> \date    2025
  !>
  !> Evaluates Y_lm_mixed at B-spline endpoints for the factorizable case. These are used later in Phase III.
  !>
  subroutine eval_BTO_CGTO_Y_lm_endpts(this, grid_r1_r2, bto_leg_part, inv_r2, Y_lm_begin, Y_lm_end)

     use grid_gbl,       only: grid_r1_r2_obj
     use phys_const_gbl, only: fourpi

     class(CGTO_shell_pw_expansion_obj)   :: this
     class(grid_r1_r2_obj), intent(in)    :: grid_r1_r2
     real(kind=cfp),        intent(in)    :: inv_r2(:)
     real(kind=cfp),        intent(out)   :: Y_lm_begin(:, :), Y_lm_end(:, :)
     real(kind=cfp),        intent(inout) :: bto_leg_part(:)

     real(kind=cfp) :: fac
     integer        :: B_start, B_end, j, l, m, lm, lbmb, cgto_m, p, n_cgto_m

        n_cgto_m = 2*this%cgto_shell%l + 1

        !$OMP DO SCHEDULE(DYNAMIC)
        do j = grid_r1_r2%first_bspline_index, grid_r1_r2%bspline_grid%n
            B_start = grid_r1_r2%bspline_start_end_r2(1, j)
            B_end = grid_r1_r2%bspline_start_end_r2(2, j)
            do l = 0, grid_r1_r2%max_l_legendre
               fac = fourpi/(2*l+1.0_cfp)
               ! II.1. starting knot
               bto_leg_part(B_start:B_end) = fac * grid_r1_r2%w2(B_start:B_end) * inv_r2(B_start:B_end)**l &
                           * grid_r1_r2%B_vals_r2(B_start:B_end, j)
               do m = -l, l
                   lm = l*l + l + m + 1
                   do lbmb = 1, (grid_r1_r2%max_bspline_l + 1)**2
                      do cgto_m = 1, n_cgto_m
                         p = this%Y_lm_non_neg_indices(cgto_m, lbmb, lm)
                         if (p == 0) cycle
                         Y_lm_begin(p, j) = sum(bto_leg_part(B_start:B_end) * this%gaunt_angular_integrals(B_start:B_end, p))
                      end do !cgto_m
                   end do !lbmb
               end do !m
               ! II.2. ending knot
               bto_leg_part(B_start:B_end) = fac * grid_r1_r2%w2(B_start:B_end) * this%r_points(B_start:B_end)**(l + 1) &
                           * grid_r1_r2%B_vals_r2(B_start:B_end, j)
               do m = -l, l
                   lm = l*l + l + m + 1
                   do lbmb = 1, (grid_r1_r2%max_bspline_l + 1)**2
                      do cgto_m = 1, n_cgto_m
                         p = this%Y_lm_non_neg_indices(cgto_m, lbmb, lm)
                         if (p == 0) cycle
                         Y_lm_end(p, j) = sum(bto_leg_part(B_start:B_end) * this%gaunt_angular_integrals(B_start:B_end, p))
                      end do !cgto_m
                   end do !lbmb
               end do !m
            end do
        end do
        !$OMP END DO

  end subroutine eval_BTO_CGTO_Y_lm_endpts


  !> \brief   Phase III of evaluation of partially integrated mixed two-electron integrals (Ylm)
  !> \authors J Benda
  !> \date    2025
  !>
  !> Evaluates Y_lm_mixed for the factorizable case, where r1 lies outside of the given B-spline's support.
  !>
  subroutine eval_BTO_CGTO_Y_lm_outside(this, grid_r1_r2, Y_lm_begin, Y_lm_end)

     use grid_gbl, only: grid_r1_r2_obj

     class(CGTO_shell_pw_expansion_obj) :: this
     class(grid_r1_r2_obj), intent(in)  :: grid_r1_r2
     real(kind=cfp),        intent(in)  :: Y_lm_begin(:, :), Y_lm_end(:, :)

     integer :: B_start, B_end, j, k, l, m, lm, lbmb, cgto_m, p, n_cgto_m

        n_cgto_m = 2*this%cgto_shell%l + 1

        !$OMP DO SCHEDULE(DYNAMIC)
        do k = 1, grid_r1_r2 % n1_total_points !r1
           do j = grid_r1_r2%first_bspline_index, grid_r1_r2%bspline_grid%n
              B_start = grid_r1_r2%bspline_start_end_r1(1, j)
              B_end = grid_r1_r2%bspline_start_end_r1(2, j)
              if (B_start <= k .and. k <= B_end) cycle
              do l = 0, grid_r1_r2 % max_l_legendre
                 do m = -l, l
                    lm = l*l + l + m + 1
                    do lbmb = 1, (grid_r1_r2%max_bspline_l + 1)**2
                       do cgto_m = 1, n_cgto_m
                          p = this%Y_lm_non_neg_indices(cgto_m, lbmb, lm)
                          if (p == 0) cycle
                          if (k < B_start) then
                             this%Y_lm_mixed(k, p, j) = Y_lm_begin(p, j) * grid_r1_r2%r1(k)**(l + 1) * grid_r1_r2%w1(k)
                          else
                             this%Y_lm_mixed(k, p, j) = Y_lm_end(p, j) * grid_r1_r2%r1(k)**(-l) * grid_r1_r2%w1(k)
                          end if
                       end do !cgto_m
                    end do !lbmb
                 end do !m
              end do !l
           end do !j
        end do !k
        !$OMP END DO

  end subroutine eval_BTO_CGTO_Y_lm_outside


  subroutine write_Y_lm_to_file(this,scratch_directory)
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this

     character(len=line_len), intent(in) :: scratch_directory
     character(len=line_len) :: path
     integer :: j, k, p, err, l3, u3

        write(level3,'("--------->","CGTO_shell_pw_expansion_obj:write_Y_lm_to_file")')

        if (.not.(this%initialized)) call xermsg('cgto_pw_expansions_mod','write_Y_lm_to_file','Object not initialized.',1,1)
        if (.not. allocated(this % Y_lm_mixed)) then
            call xermsg ('cgto_pw_expansions_mod', 'write_Y_lm_to_file', 'The Y_lm has not been allocated.', 2, 1)
        end if

        write(path,'(i10)') this%cgto_shell_index
        this%Y_lm_path = trim(scratch_directory)//'Y_lm_'//trim(adjustl(path))
        this%Y_lm_on_disk = .true.
        l3 = lbound(this%Y_lm_mixed,3) !grid_r1_r2%first_bspline_index
        u3 = ubound(this%Y_lm_mixed,3) !grid_r1_r2%bspline_grid%n
        this%Y_lm_dim_on_disk(1:4) = (/size(this%Y_lm_mixed,1),l3,u3,size(this%Y_lm_mixed,2)/)
        write(level1,'("Writing Y_lm to the file: ",a)') trim(this%Y_lm_path)

        open(newunit=this%Y_lm_unit,file=this%Y_lm_path,status='replace',form='unformatted',access='stream',iostat=err)
        if (err .ne. 0)  call xermsg('cgto_pw_expansions_mod','write_Y_lm_to_file','Error opening Y_lm file.',err,1)
       
        !write(this%Y_lm_unit) size(this%Y_lm_mixed,1), grid_r1_r2%first_bspline_index, grid_r1_r2%bspline_grid%n, size(this%Y_lm_mixed,3), cfp_bytes
        if (allocated(this%Y_lm_disk_offset)) deallocate(this%Y_lm_disk_offset)
        allocate(this%Y_lm_disk_offset(size(this%Y_lm_mixed,2),l3:u3),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_BTO_CGTO_Y_lm','Memory allocation 4 failed.',err,1)

        k = size(this%Y_lm_mixed,1)
        do j=l3,u3
           do p=1,size(this%Y_lm_mixed,2)
              inquire(unit=this%Y_lm_unit,pos=this%Y_lm_disk_offset(p,j)) !save position of the next byte 
              write(this%Y_lm_unit) this%Y_lm_mixed(1:k,p,j)
           enddo !p
        enddo !j

        write(level2,'("Memory required for the offset array (MiB): ",f25.10)') &
            (u3-l3+1)*size(this%Y_lm_disk_offset,1)*(bit_size(err)/8)/(Mib*1.0)

        write(level3,'("<---------","CGTO_shell_pw_expansion_obj:write_Y_lm_to_file")')

  end subroutine write_Y_lm_to_file

  subroutine read_Y_lm_from_file(this,p,j,Y_lm_mixed)
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this
     integer, intent(in) :: p, j
     real(kind=cfp), allocatable :: Y_lm_mixed(:) !must be allocated before

     integer :: pos_p_j, j_offset, n_points, n_p

        if (.not. this % initialized) then
            call xermsg ('cgto_pw_expansions_mod', 'read_Y_lm_from_file', 'Object not initialized.', 1, 1)
        end if

        if (.not. this % Y_lm_on_disk) then
            call xermsg ('cgto_pw_expansions_mod', 'read_Y_lm_from_file', 'Y_lm has not been saved to disk.', 2, 1)
        end if

        j_offset = j-this%Y_lm_dim_on_disk(2)+1
        if (j_offset < 0) call xermsg('cgto_pw_expansions_mod','read_Y_lm_from_file','On input the value of j was incorrect.',3,1)

        n_points = this%Y_lm_dim_on_disk(1)
        n_p = this%Y_lm_dim_on_disk(4)
        !pos_p_j = ((j_offset-1)*n_p + p-1)*n_points*cfp_bytes + 1
        pos_p_j = this%Y_lm_disk_offset(p,j)

        read(this%Y_lm_unit,pos=pos_p_j) Y_lm_mixed(1:n_points)

  end subroutine read_Y_lm_from_file

  subroutine eval_NAI_X_lm_projections(this,grid_r1_r2,nuclei)
     use special_functions_gbl, only: cfp_besi, cfp_eval_poly_horner_many
     use phys_const_gbl, only: pi, fourpi
     use grid_gbl, only: grid_r1_r2_obj
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this
     class(grid_r1_r2_obj) :: grid_r1_r2
     type(nucleus_type) :: nuclei(:)

     type(pw_expansion_obj) :: u_grid
     integer :: i,u,n,p,lm,l,m,n_Xlm,l_max,lpmp,ma,lp,mp,lap,map,m_ind
     real(kind=cfp) :: val, S_vec(3), R_AC_vec(3), R_AC_sq, fac, arg, tol, exp_part, one_m_u_sq, u_sq, RA, cf
     integer, parameter :: kode = 2
     real(kind=cfp), parameter :: half = 0.5_wp
     integer :: nz, err, besi_dim, max_l_aux
     real(kind=cfp), allocatable :: y(:), S(:), exp_fac(:), S_X_lm(:,:), RH(:,:), bessel_fac(:), int_tmp(:,:), &
                                    transl_cfs(:,:,:), Xlm_CGTO_center(:), d_cfs(:,:), NAI_s_X_lm_projections(:,:)
     logical :: non_zero, cgto_center_is_cms, R_AC_sq_is_zero

        write(level3,'("--------->","CGTO_shell_pw_expansion_obj:eval_NAI_X_lm_projections")')

        if (.not.(this%initialized)) call xermsg('cgto_pw_expansions_mod','eval_NAI_X_lm_projections','Object not initialized.',1,1)

        tol = F1MACH(4,cfp_dummy)

        !Maximum L for which the auxiliary angular projections have to be calculated:
        max_l_aux = grid_r1_r2%max_bspline_l + this%cgto_shell%l

        besi_dim = max_l_aux+1
        allocate(y(besi_dim),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_NAI_X_lm_projections','Memory allocation 1 failed.',err,1)

        val = 20.0_cfp !min(maxval(this%cgto_shell%exponents),20.0_cfp) !exponents larger than approx. 20 lead to a too small spacing between the first couple of integration intervals (for 50 intervals total).
        call u_grid%eval_exponential_grid(val,55,.true.)
        !call u_grid%eval_regular_grid(0.0_cfp,1.0_cfp,0.005_cfp)

        n_Xlm = (max_l_aux+1)**2

        allocate(S(u_grid%n_total_points),exp_fac(u_grid%n_total_points),S_X_lm(u_grid%n_total_points,n_Xlm),&
                 RH(-max_l_aux:max_l_aux,0:max_l_aux+1),bessel_fac(0:max_l_aux+1),&
                 NAI_s_X_lm_projections(grid_r1_r2%n1_total_points,n_Xlm),int_tmp(n_Xlm,grid_r1_r2%n1_total_points),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_NAI_X_lm_projections','Memory allocation 2 failed.',err,1)

        cgto_center_is_cms = .true.
        do i=1,3
           if (this%cgto_shell%center(i) .ne. 0.0_cfp) cgto_center_is_cms = .false.
        enddo

        NAI_s_X_lm_projections = 0.0_cfp !projection(r1,lm)
        do n=1,size(nuclei)
           if (abs(nuclei(n)%charge - nuclei(n)%ecp%n_core) .le. tol) cycle

           R_AC_vec = -(this%cgto_shell%center - nuclei(n)%center)
           R_AC_sq = dot_product(R_AC_vec,R_AC_vec)
           R_AC_sq_is_zero = (abs(R_AC_sq) .le. tol)

           do p=1,this%cgto_shell%number_of_primitives
              !todo why should there not be -nuclei(n)%charge? It appears that
              !to get agreement with BG_nai_integrals the minus sign must be
              !omitted here...
              fac = 2.0_cfp/sqrt(pi) * (nuclei(n)%charge - nuclei(n)%ecp%n_core) &
                 * this%cgto_shell%contractions(p)*this%cgto_shell%norm*this%cgto_shell%norms(p)
              do u=1,u_grid%n_total_points
                 u_sq = u_grid%r_points(u)**2
                 one_m_u_sq = (1.0_cfp-u_grid%r_points(u))*(1.0_cfp+u_grid%r_points(u))

                 if (R_AC_sq_is_zero .and. cgto_center_is_cms) then !Prevent division by S(u) = 0 that would happen in this case, see the else branch.
                    S(u) = 0.0_cfp
                    S_vec = 0.0_cfp
                 else
                    S_vec(1:3) = u_sq*R_AC_vec + this%cgto_shell%center
                    S(u) = sqrt(dot_product(S_vec,S_vec))
                    S_vec(1:3) = S_vec(1:3)/S(u) !normalize the S-vector as required by resh
                 endif
                 call resh(RH,S_vec(1),S_vec(2),S_vec(3),max_l_aux)

                 lm = 0
                 do l=0,max_l_aux
                    do m=-l,l
                       lm = lm + 1
                       S_X_lm(u,lm) = RH(m,l)
                    enddo !m
                 enddo !l

                 exp_fac(u) = exp(-this%cgto_shell%exponents(p)*u_sq*R_AC_sq) &
                                *sqrt(this%cgto_shell%exponents(p))*(one_m_u_sq)**(-1.5_cfp)

              enddo !u
              
              do i=1,grid_r1_r2%n1_total_points

                 !Integrate over u for each X_lm projection and each r1 point:
                 int_tmp = 0.0_cfp
                 do u=u_grid%n_total_points,1,-1 !sum from the smallest values
                    u_sq = u_grid%r_points(u)**2
                    one_m_u_sq = (1.0_cfp-u_grid%r_points(u))*(1.0_cfp+u_grid%r_points(u))

                    if (grid_r1_r2%r1(i)*S(u) .le. tol) then !Evaluate limit for S*r -> 0
                       bessel_fac(0) = 1.0_cfp
                       bessel_fac(1:max_l_aux) = 0.0_cfp
                    else
                       arg = 2*this%cgto_shell%exponents(p)/one_m_u_sq*grid_r1_r2%r1(i)*S(u)
                       call cfp_besi(arg, half, kode, besi_dim, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
                       val = sqrt(pi/(2.0_cfp*arg))
                       do l=0,max_l_aux
                          bessel_fac(l) = y(l+1)*val
                       enddo !l
                    endif

                    exp_part = fac*fourpi*exp_fac(u)*exp(-this%cgto_shell%exponents(p)/one_m_u_sq*(grid_r1_r2%r1(i)-S(u))**2) !*u_grid%weights(t)
                   
                    lm = 0
                    do l=max_l_aux,0,-1
                       do m=-l,l
                          lm = l*l+l+m+1 !lm + 1
                          val = exp_part*bessel_fac(l)*S_X_lm(u,lm)
                          !NAI_s_X_lm_projections(i,lm) = NAI_s_X_lm_projections(i,lm) + val*u_grid%weights(u)
                          int_tmp(lm,i) = int_tmp(lm,i) + val*u_grid%weights(u)
                          !if (lm .eq. 1) write(100+lm+((max_l_aux+1)**2)*(n-1),'(3e25.15)') grid_r1_r2%r1(i),u_grid%r_points(u),val
                       enddo !m
                    enddo !l
                 enddo !u
                 do lm=1,(max_l_aux+1)**2
                    NAI_s_X_lm_projections(i,lm) = NAI_s_X_lm_projections(i,lm) + int_tmp(lm,i)
                    !if (lm .eq. 1) write(100+lm+((max_l_aux+1)**2)*(n-1),'("")')
                 enddo

              enddo !i

           enddo !p
        enddo !n

        if (allocated(this%NAI_X_lm_projections)) deallocate(this%NAI_X_lm_projections)
        allocate(this%NAI_X_lm_projections(grid_r1_r2%n1_total_points,2*this%cgto_shell%l+1,(grid_r1_r2%max_bspline_l+1)**2),&
                 stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_NAI_X_lm_projections','Memory allocation 3 failed.',err,1)

        this%NAI_X_lm_projections = 0.0_cfp

        if (this%cgto_shell%l > 0) then

           !Real spherical harmonics for the nuclei: result in the module array Xlm_CGTO_center
           call precalculate_Xlm_for_CGTO_center(this%cgto_shell%center,this%cgto_shell%l,Xlm_CGTO_center)
      
           RA = sqrt(dot_product(this%cgto_shell%center,this%cgto_shell%center))
   
           !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
           call precalculate_solh_translation_coeffs(this%cgto_shell%l,RA,Xlm_CGTO_center,transl_cfs)
   
           allocate(d_cfs(grid_r1_r2%n1_total_points,1:this%cgto_shell%l+1),stat=err)
           if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_NAI_X_lm_projections','Memory allocation 4 failed.',err,1)
   
           do l=0,grid_r1_r2%max_bspline_l
              do m=-l,l
                 lm = l*l+l+m+1
   
                 do ma=-this%cgto_shell%l,this%cgto_shell%l
                    m_ind = ma+this%cgto_shell%l+1

                    d_cfs = 0.0_cfp
                    do lap=0,this%cgto_shell%l

                       non_zero = .false.
   
                       do lp=abs(l-lap),l+lap
                          do mp=-lp,lp
                             lpmp = lp*lp+lp+mp+1
   
                             cf = 0.0_cfp
                             do map=-lap,lap
                                cf = cf + transl_cfs(map+lap+1,lap,m_ind)*cpl%rgaunt(l,lap,lp,m,map,mp)
                             enddo !map

                             if (cf .ne. 0.0_cfp) then
                                non_zero = .true.
                                d_cfs(1:grid_r1_r2%n1_total_points,lap+1) = d_cfs(1:grid_r1_r2%n1_total_points,lap+1) &
                                                    + cf*NAI_s_X_lm_projections(1:grid_r1_r2%n1_total_points,lpmp)
                             endif
                          enddo !mp
                       enddo !lp
   
                    enddo !lap

                    if (non_zero) then
                       call cfp_eval_poly_horner_many (this % cgto_shell % l, &
                                                       grid_r1_r2 % r1, grid_r1_r2 % n1_total_points, d_cfs, &
                                                       this % NAI_X_lm_projections(1 : grid_r1_r2 % n1_total_points, m_ind, lm))
                    end if
   
                 enddo !ma
   
              enddo !m
           enddo !l

        else
           lm = (grid_r1_r2%max_bspline_l+1)**2
           this%NAI_X_lm_projections(1:grid_r1_r2%n1_total_points,1,1:lm) &
                = NAI_s_X_lm_projections(1:grid_r1_r2%n1_total_points,1:lm)
        endif

        deallocate(NAI_s_X_lm_projections)

        write(level3,'("<---------","CGTO_shell_pw_expansion_obj:eval_NAI_X_lm_projections")')

  end subroutine eval_NAI_X_lm_projections

  subroutine eval_at_lebedev_points(this,grid_r1_r2)
     use grid_gbl, only: grid_r1_r2_obj
     implicit none
     class(CGTO_shell_pw_expansion_obj) :: this
     class(grid_r1_r2_obj) :: grid_r1_r2
     
     integer :: err, i, j, k, m
     real(kind=cfp) :: SH(-this%cgto_shell%l:this%cgto_shell%l,0:this%cgto_shell%l+1)
     real(kind=cfp) :: r(3), r_square, sum_exp, x

        if (.not.(this%initialized)) call xermsg('cgto_pw_expansions_mod','eval_at_lebedev_points','Object not initialized.',1,1)

        if (allocated(this%at_lebedev_points)) deallocate(this%at_lebedev_points)
        allocate(this%at_lebedev_points(grid_r1_r2%lebedev_order,grid_r1_r2%n1_total_points,2*this%cgto_shell%l+1),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_at_lebedev_points','Memory allocation failed.',err,1 )
        this%at_lebedev_points = 0.0_cfp

        do i=1,grid_r1_r2%n1_total_points
           x = grid_r1_r2%r1(i)
           do j=1,grid_r1_r2%lebedev_order
              r(1:3) = (/x*grid_r1_r2%leb_r1(j)-this%cgto_shell%center(1),x*grid_r1_r2%leb_r2(j) &
                            -this%cgto_shell%center(2),x*grid_r1_r2%leb_r3(j)-this%cgto_shell%center(3)/)
      
              if (this%cgto_shell%l > 0) then
                 call solh(SH,r(1),r(2),r(3),this%cgto_shell%l)
              else
                 SH(0,0) = 1.0_cfp
              endif
      
              r_square = dot_product(r,r)
              sum_exp = 0.0_wp
              do k=1,this%cgto_shell%number_of_primitives
                 sum_exp = sum_exp + this%cgto_shell%contractions(k)*this%cgto_shell%norms(k)&
                            *exp(-this%cgto_shell%exponents(k)*r_square)
              enddo
              sum_exp = this%cgto_shell%norm*sum_exp
      
               do m=-this%cgto_shell%l,this%cgto_shell%l
                  this%at_lebedev_points(j,i,m+this%cgto_shell%l+1) = SH(m,this%cgto_shell%l)*sum_exp
               enddo !m_b
            enddo !j
         enddo !i

  end subroutine eval_at_lebedev_points

  subroutine eval_CGTO_shell_pair_pw_expansion(this, neglected)
     use general_quadrature_gbl, only: n_10,w_10,x_10
     use const_gbl, only: epsabs
     use omp_lib
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     integer, optional, intent(out) :: neglected

     integer :: err, lm, cgto_m, i
     real(kind=cfp) :: threshold, max_value, max_value_global

        if (.not. this % initialized) then
            call xermsg ('cgto_pw_expansions_mod', 'eval_CGTO_shell_pair_pw_expansion', 'Object not initialized.', 1, 1)
        end if
!
        !todo is x_10 necessary? Can we use x_7 only?
!        call this%eval_CGTO_pair_radial_grid(a,x_10,w_10,n_10)
!
!------ Calculate the projections of the CGTOs on the real spherical harmonics (i.e. PW-expansion) for all radial grid points.
        call omp_calculate_CGTO_pair_pw_coefficients_analytic (max_l_pw, this % cgto_shell_A, this % cgto_shell_B, &
                                                                this % r_points, this % angular_integrals)

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs

        if (allocated(this%neglect_m_lm)) deallocate(this%neglect_m_lm)
        allocate(this%neglect_m_lm(size(this%angular_integrals,2),size(this%angular_integrals,3)),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_CGTO_shell_pair_pw_expansion','Memory allocation failed.',err,1)

        if (present(neglected)) neglected = 0
        this%neglect_m_lm = .false.
        max_value_global = 0.0_cfp
        do lm=1,size(this%angular_integrals,3)
           do cgto_m=1,size(this%angular_integrals,2)
              max_value = 0.0_cfp
              do i=1,size(this%angular_integrals,1)
                 max_value = max(max_value,abs(this%angular_integrals(i,cgto_m,lm)))
              enddo !i
              if (max_value < threshold) then
                 if (present(neglected)) neglected = neglected + 1
                 this%neglect_m_lm(cgto_m,lm) = .true.
              else
                 max_value_global = max(max_value,max_value_global)
              endif
           enddo !cgto_m
        enddo !lm

        !write(stdout,'("Neglected vs total m,lm indices: ",i0,1x,i0)') &
        !    neglected, size(this%angular_integrals,2)*size(this%angular_integrals,3)
        !if (neglected .eq. 0) write(stdout,'("Pw expansion probably too short, max_value is: ",e25.15)') max_value_global

  end subroutine eval_CGTO_shell_pair_pw_expansion

  subroutine eval_coulomb_integrals(this,grid_r1_r2)
     use const_gbl, only: epsrel, epsabs
     use phys_const_gbl, only: fourpi
     use cgto_hgp_gbl, only: sph_nari
     use grid_gbl, only: grid_r1_r2_obj
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     class(grid_r1_r2_obj) :: grid_r1_r2

     integer :: err, i, j, k, ind_a, ind_b, n_integrals
     real(kind=cfp), allocatable :: nari(:)
     real(kind=cfp) :: point(3)
     integer, allocatable :: int_index(:,:)

        if (.not.(this%initialized)) call xermsg('cgto_pw_expansions_mod','eval_coulomb_integrals','Object not initialized.',1,1)
!
        n_integrals = (2*this%cgto_shell_A%l+1)*(2*this%cgto_shell_B%l+1)

        if (allocated(this%coulomb_integrals)) deallocate(this%coulomb_integrals)
        allocate(nari(n_integrals),int_index(2,n_integrals),this%coulomb_integrals(grid_r1_r2%lebedev_order,&
                    grid_r1_r2%n1_total_points,n_integrals),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_coulomb_integrals','Memory allocation failed.',err,1)

        if (this%cgto_shell_A%l .ge. this%cgto_shell_B%l) then
           this%order_AB = .true.
        else
           this%order_AB = .false.
        endif

        ind_a = 1
        ind_b = 1+2*this%cgto_shell_A%l+1
        do i=1,grid_r1_r2%n1_total_points
           do j=1,grid_r1_r2%lebedev_order
              point(1:3) = (/grid_r1_r2%r1(i)*grid_r1_r2%leb_r1(j),grid_r1_r2%r1(i) &
                                    *grid_r1_r2%leb_r2(j),grid_r1_r2%r1(i)*grid_r1_r2%leb_r3(j)/)

              call sph_nari (this % cgto_shell_A % number_of_primitives, &
                             this % cgto_shell_A % center(1), &
                             this % cgto_shell_A % center(2), &
                             this % cgto_shell_A % center(3), &
                             this % cgto_shell_A % norm, &
                             this % cgto_shell_A % norms, &
                             this % cgto_shell_A % l, &
                             this % cgto_shell_A % exponents, &
                             this % cgto_shell_A % contractions, ind_a, &
                             this % cgto_shell_B % number_of_primitives, &
                             this % cgto_shell_B % center(1), &
                             this % cgto_shell_B % center(2), &
                             this % cgto_shell_B % center(3), &
                             this % cgto_shell_B % norm, &
                             this % cgto_shell_B % norms, &
                             this % cgto_shell_B % l, &
                             this % cgto_shell_B % exponents, &
                             this % cgto_shell_B % contractions, ind_b, &
                             point(1), point(2), point(3), nari, int_index)

              this%coulomb_integrals(j,i,1:n_integrals) = nari(1:n_integrals)*grid_r1_r2%leb_w(j)*fourpi
              
           enddo !j
        enddo !i

  end subroutine eval_coulomb_integrals

  !> Calculates the radial integral 2-electron integrals between all pairs of BTOs and all partial waves of the GG pair.
  subroutine eval_radial_GG_BB(this,grid_r1_r2)
     use grid_gbl, only: grid_r1_r2_obj
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     class(grid_r1_r2_obj) :: grid_r1_r2

     integer :: i, j, pair_index, cgto_m, n_cgto_m, l, m, lm, max_l, err, n, n1_points

        if ((.not.allocated(this%angular_integrals)) .or. (.not.allocated(grid_r1_r2%lambda_BB_r1_l_ij))) then
           call xermsg('cgto_pw_expansions_mod','eval_radial_GG_BB',&
                       'At least one of angular_integrals, lambda_BB_r1_l_ij is not initialized.',1,1)
        endif

        n_cgto_m = size(this%angular_integrals,2)
        max_l = ubound(grid_r1_r2%lambda_BB_r1_l_ij,2)
        n = grid_r1_r2%bspline_grid%n

        err = resize_array(this%radial_lm_BB_GG, (max_l+1)**2, n_cgto_m, n*(n+1)/2)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','eval_radial_GG_BB','Memory allocation failed.',err,1)

        n1_points = grid_r1_r2%n1_total_points

        do i=grid_r1_r2%first_bspline_index,grid_r1_r2%bspline_grid%n
           do j=grid_r1_r2%first_bspline_index,i
              pair_index = i*(i-1)/2+j
              do cgto_m=1,n_cgto_m
                 do l=0,max_l
                    do m=-l,l
                       lm = l*l+l+m+1
                       this%radial_lm_BB_GG(lm,cgto_m,pair_index) &
                            = sum(this%angular_integrals(1:n1_points,cgto_m,lm)*&
                                  grid_r1_r2%lambda_BB_r1_l_ij(1:n1_points,l,pair_index)*grid_r1_r2%r1(1:n1_points))
                    enddo !m
                 enddo !l
              enddo !cgto_m
           enddo !j
        enddo !i

  end subroutine eval_radial_GG_BB

  subroutine eval_damped_dipole_integrals(this,shell_A,shell_B,A,B,starting_index_A,starting_index_B,dipole_damp_factor,&
                                          prop_column,integrals,int_index)
     use phys_const_gbl, only : fourpi
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     type(CGTO_shell_data_obj), intent(in) :: shell_A, shell_B
     integer, intent(in) :: A,B, starting_index_A, starting_index_B, prop_column
     real(kind=cfp), intent(in) :: dipole_damp_factor
     !We assume that these two arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)

     
     integer :: err, lp, mp, lp_mp, cgto_m, n_cgto_m, i, Mg_A_ind, Mg_B_ind, Mg_A, Mg_B, n_cgto_B_m
     real(kind=cfp) :: fac, prop
     real(kind=cfp), allocatable :: radial_prop(:)

       if (this%n_total_points <= 0) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_damped_dipole_integrals', 'The radial grid has not been initialized.', 1, 1)
       endif

       call this%init_CGTO_shell_pair_pw_expansion(shell_A,A,shell_B,B)
       call this%eval_CGTO_shell_pair_pw_expansion

       allocate(radial_prop(this%n_total_points),stat=err)
       if (err /= 0) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_damped_dipole_integrals', 'Memory allocation failed.', err, 1)
       endif

       lp = 1
       n_cgto_m = (2*shell_B%l+1)*(2*shell_A%l+1)

       fac = sqrt(fourpi/(2*lp+1.0_cfp)) !this is for the dipolar solid harmonic
       do i = 1, this%n_total_points
          radial_prop(i) = fac*this%r_points(i)*exp(-dipole_damp_factor*this%r_points(i))*this%weights(i) * this%r_points(i)**2 !r^2 is Jacobian
       enddo

       !For all property components (lp,mp)
       do mp=-lp,lp
          lp_mp = lp*lp + lp+mp+1
          !For all m-values of the CGTO-CGTO pair
          do cgto_m=1,n_cgto_m
             if (.not.this%neglect_m_lm(cgto_m,lp_mp)) then
                !note that radial_prop contains the quadrature weights
                prop = sum(radial_prop(1:this%n_total_points)*this%angular_integrals(1:this%n_total_points,cgto_m,lp_mp))
             else
                prop = 0.0_cfp
             endif
             !todo before saving I can print out the integrals(cgto_m,prop_column+lp_mp-1) which should be the dipoles calculated over all space (cf. GG_shell_integrals) so I can compare them here.
             !     However, remember that the analytic results in integrals(:,:) are over all space, i.e. they don't have the tails subtracted.
             !if (abs((integrals(cgto_m,prop_column+lp_mp-1)-prop)/prop) > 10e-10_cfp .and. prop /= 0.0_cfp) then
             !   print *,'damped dipole',integrals(cgto_m,prop_column+lp_mp-1),prop
             !endif
             !The integrals are stored in the order: (CGTO_B_m,CGTO_A_m)
             integrals(cgto_m,prop_column+lp_mp-1) = prop
          enddo !cgto_m
       enddo !mp
!
!------ Generate the basis function indices
!
       n_cgto_B_m = 2*shell_B%l+1
       i = 0
       do Mg_A=-shell_A%l,shell_A%l
          Mg_A_ind = Mg_A+shell_A%l+1
          do Mg_B=-shell_B%l,shell_B%l
             Mg_B_ind = Mg_B+shell_B%l+1
             i = i + 1
             int_index(1,i) = max(Mg_A_ind-1+starting_index_A, Mg_B_ind-1+starting_index_B)
             int_index(2,i) = min(Mg_A_ind-1+starting_index_A, Mg_B_ind-1+starting_index_B)
          enddo !Mg_B
       enddo !Mg_A

  end subroutine eval_damped_dipole_integrals

  subroutine eval_GGG_integrals(this,cgto_shell,starting_index_A,starting_index_B,ggg_integrals,int_index)
     use phys_const_gbl, only : fourpi
     use gto_routines_gbl, only: index_1el
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     type(CGTO_shell_data_obj), intent(in) :: cgto_shell
     integer, intent(in) :: starting_index_A, starting_index_B
     !We assume that these two arrays have been allocated to the appropriate dimensions:
     real(kind=cfp), allocatable :: ggg_integrals(:,:)
     integer, allocatable :: int_index(:,:)
     
     integer :: err, lp, mp, lp_mp, cgto_m, n_cgto_m, i, j
     real(kind=cfp) :: ggg, r_square
     real(kind=cfp), allocatable :: contraction(:)
     logical, parameter :: calculate_indices = .false.

       if (this%n_total_points <= 0) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_GGG_integrals', 'The radial grid has not been initialized.', 1, 1)
       endif

       allocate(contraction(this%n_total_points),stat=err)
       if (err /= 0) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_GGG_integrals', 'Memory allocation failed.', err, 1)
       endif

       if (.not.allocated(this%angular_integrals)) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_GGG_integrals', 'The PW must be precalculated.', 2, 1)
       endif

       n_cgto_m = (2*this%cgto_shell_B%l+1)*(2*this%cgto_shell_A%l+1)

       do i=1,this%n_total_points
          !Radial part of the CGTO assumed to be in the center of coordinates
          r_square = this%r_points(i)*this%r_points(i)
          contraction(i) = 0.0_cfp
          do j=1,cgto_shell%number_of_primitives
             contraction(i) = contraction(i) &
                               + cgto_shell%contractions(j)*cgto_shell%norms(j)*exp(-cgto_shell%exponents(j)*r_square)
          enddo
          contraction(i) = contraction(i)*cgto_shell%norm*sqrt(fourpi/(2*cgto_shell%l+1.0_cfp))*this%r_points(i)**cgto_shell%l&
                          &*this%weights(i)*r_square !r_square is the Jacobian
       enddo !i

       !For all cgto_shell angular components (lp,mp)
       lp = cgto_shell%l
       do mp=-lp,lp
          lp_mp = lp*lp + lp+mp+1
          !For all m-values of the CGTO-CGTO pair
          do cgto_m=1,n_cgto_m
             if (.not.this%neglect_m_lm(cgto_m,lp_mp)) then
                !note that contraction contains the quadrature weights
                ggg = sum(contraction(1:this%n_total_points)*this%angular_integrals(1:this%n_total_points,cgto_m,lp_mp))
             else
                ggg = 0.0_cfp
             endif
             !The integrals are stored in the order: (CGTO_B_m,CGTO_A_m)
             ggg_integrals(cgto_m,lp+mp+1) = ggg
          enddo !cgto_m
       enddo !mp
!
!------ Check that the integrals have been generated in the expected order given by int_index
!
       call index_1el(this%cgto_shell_B%l,this%cgto_shell_A%l,starting_index_B,starting_index_A,1,int_index,calculate_indices)

  end subroutine eval_GGG_integrals

  subroutine eval_GBG_integrals(this,bto_shell,starting_index_A,starting_index_B,gbg_integrals,int_index)
     use phys_const_gbl, only : fourpi
     use gto_routines_gbl, only: index_1el
     use bspline_base_gbl, only: bvalu
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     type(BTO_shell_data_obj), target, intent(in) :: bto_shell
     integer, intent(in) :: starting_index_A, starting_index_B
     !We assume that these two arrays have been allocated to the appropriate dimensions:
     real(kind=cfp), allocatable :: gbg_integrals(:,:)
     integer, allocatable :: int_index(:,:)
     
     type(BTO_shell_data_obj) :: bto_shell_dummy
     integer :: err, lp, mp, lp_mp, cgto_m, n_cgto_m, i, i_min, i_max
     real(kind=cfp) :: gbg, r_jac, bspline_val, r1, r2
     real(kind=cfp), allocatable :: bto(:)
     logical, parameter :: calculate_indices = .false.

       if (this%n_total_points <= 0) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_GBG_integrals', 'The radial grid has not been initialized.', 1, 1)
       endif

       bto_shell_dummy = bto_shell

       allocate(bto(this%n_total_points),stat=err)
       if (err /= 0) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_GBG_integrals', 'Memory allocation failed.', err, 1)
       endif

       if (.not.allocated(this%angular_integrals)) then
          call xermsg ('cgto_pw_expansions_mod', 'eval_GBG_integrals', 'The PW must be precalculated.', 2, 1)
       endif

       n_cgto_m = (2*this%cgto_shell_B%l+1)*(2*this%cgto_shell_A%l+1)

       i_min = this%n_total_points + 1
       i_max = 0

       do i=1,this%n_total_points
          !Radial part of the BTO assumed to be in the center of coordinates
          r_jac = this%r_points(i) !This is the Jacobian r**2 multiplied by 1/r coming from the B-spline

          call bto_shell_dummy%bspline_grid%bspline_range(bto_shell_dummy%bspline_index,r1,r2)

          bspline_val = 0.0_cfp
          if (this%r_points(i) >= r1 .and. this%r_points(i) <= r2) then

             i_min = min(i,i_min)
             i_max = max(i,i_max)

             bto_shell_dummy%bspline_grid%bcoef = 0.0_cfp
             bto_shell_dummy%bspline_grid%bcoef(bto_shell_dummy%bspline_index) = 1.0_cfp
             bspline_val = bto_shell_dummy % norm * bvalu(bto_shell_dummy % bspline_grid % knots, &
                                               bto_shell_dummy % bspline_grid % bcoef, &
                                               bto_shell_dummy % bspline_grid % n,&
                                               bto_shell_dummy % bspline_grid % order, 0, this%r_points(i), &
                                               bto_shell_dummy % bspline_grid % inbv, &
                                               bto_shell_dummy % bspline_grid % work)

          endif

          bto(i) = bspline_val*this%weights(i)*r_jac
       enddo !i

       !For all bto_shell_dummy angular components (lp,mp)
       lp = bto_shell_dummy%l
       do mp=-lp,lp
          lp_mp = lp*lp + lp+mp+1
          !For all m-values of the CGTO-CGTO pair
          do cgto_m=1,n_cgto_m
             if (.not.this%neglect_m_lm(cgto_m,lp_mp)) then
                !note that contraction contains the quadrature weights
                gbg = sum(bto(i_min:i_max)*this%angular_integrals(i_min:i_max,cgto_m,lp_mp))
             else
                gbg = 0.0_cfp
             endif
             !The integrals are stored in the order: (CGTO_B_m,CGTO_A_m)
             gbg_integrals(cgto_m,lp+mp+1) = gbg
          enddo !cgto_m
       enddo !mp
!
!------ Check that the integrals have been generated in the expected order given by int_index
!
       call index_1el(this%cgto_shell_B%l,this%cgto_shell_A%l,starting_index_B,starting_index_A,1,int_index,calculate_indices)

  end subroutine eval_GBG_integrals

  subroutine eval_CGTO_pair_radial_grid(this,a,x,w,n)
     use const_gbl, only: epsrel, epsabs
     implicit none
     class(CGTO_shell_pair_pw_expansion_obj) :: this
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: x(2*n+1), w(2*n+1), a

     real(kind=cfp) :: threshold,d,R_min,R_max, delta = 1.0_cfp !todo determine adaptively depending on the compactness of the product pair

        if (.not. this % initialized) then
            call xermsg ('cgto_pw_expansions_mod', 'eval_CGTO_pair_radial_grid', 'Object not initialized.', 1, 1)
        end if

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs
!
!------ For the present CGTO shrink the radial grid defined for the whole B-spline range [bspline_grid%A,bspline_grid%B]
!       into the range which corresponds to the extent of the CGTO [bspline_grid%A,R].
!       The resulting quadrature points and weights are in the arrays r_points,weights.
!
!------ Estimate radius of the CGTO
!        call this%cgto_shell%estimate_shell_radius(threshold,d,R_min,R_max)
        !todo test
        R_min = 0.0_cfp
        R_max = a
!
        call radial_grid_CGTO_pair (this % cgto_shell_A, this % cgto_shell_B, R_min, R_max, x, w, n, delta, &
                                    this % r_points, this % weights, this % n_total_points)

  end subroutine eval_CGTO_pair_radial_grid

  !> Radial grid used when integrating over the molecular orbitals purely numerically.
  subroutine radial_grid_mo(centers,R_min,R_max,bspline_grid,delta_cgto_grid,x,w,n,r_points,weights,n_total_points)
     use sort_gbl, only: sort_float
     implicit none
     type(bspline_grid_obj), intent(in) :: bspline_grid
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: R_min, R_max, x(2*n+1), w(2*n+1), centers(:), delta_cgto_grid
     !OUTPUT variables:
     integer, intent(out) :: n_total_points
     real(kind=cfp), allocatable :: r_points(:), weights(:)

     integer :: n_points, err, interval, cnt, i, j, k, n_ranges
     real(kind=cfp) :: d, A, B, delta, R_start, R_end, test
     real(kind=cfp), allocatable :: list(:,:), range_start_end(:,:)

         n_points = 2*n+1 !Number of quadrature points within each quadrature interval
         if (R_min >= bspline_grid%B) then
            call xermsg ('cgto_pw_expansions_mod', 'radial_grid_mo', &
                         'All CGTOs must start below the end point of the B-spline grid.', 2, 1)
         end if

         !write(stdout,'("CGTO negligible beyond radius [a.u.]: ",e25.15)') R_max

         R_start = R_min
         R_end = min(R_max,bspline_grid%B)

         if (delta_cgto_grid <= 0.0_cfp) then
            call xermsg ('cgto_pw_expansions_mod', 'radial_grid_mo', 'On input delta_cgto_grid was .le. 0.0_cfp.', 3, 1)
         end if

         j = ceiling((R_end-R_start)/delta_cgto_grid)
         allocate(list(bspline_grid%no_knots+size(centers)+j+2,1),range_start_end(2,bspline_grid%no_knots+size(centers)+j),stat=err)
         if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','radial_grid_mo','Memory allocation error.',err,1)
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
            if (centers(i) .ge. R_start .and. centers(i) .le. R_end) then
               k = k + 1
               list(k,1) = centers(i)
            endif
         enddo

         do i=1,j
            test = R_start + i*delta_cgto_grid
            if (test .ge. R_start .and. test .le. R_end) then
               k = k + 1
               list(k,1) = test
            endif
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
         allocate(r_points(n_total_points),weights(n_total_points),stat=err)
         if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','radial_grid_mo','Memory allocation failed.',err,1)

         !Construct quadratures for the individual ranges:
         cnt = 0
         do j=1,n_ranges
            A = range_start_end(1,j)
            delta = (range_start_end(2,j)-range_start_end(1,j))
            B = A+delta
            B = range_start_end(2,j)
            !Prepare the quadrature grid for the interval [A,B] expanding the canonical G-L quadrature.
            !print *,'quad',A,B,cnt
            call gl_expand_A_B(x,w,n,r_points(cnt+1:cnt+n_points),weights(cnt+1:cnt+n_points),A,B)
            cnt = cnt + n_points
            A = A + delta
         enddo
!
         if (cnt /= n_total_points) then
            call xermsg ('cgto_pw_expansions_mod', 'radial_grid_mo', 'Error constructing the radial grid.', 3, 1)
         end if

  end subroutine radial_grid_mo

  !> Determines the radial grid needed to describe, integrals involving the CGTO and the B-splines on the given grid. We assume that if the CGTO spanned the whole B-spline grid
  !> it would be sufficient to use the given quadrature rule (x,w,n) within each knot interval. For CGTOs whose radial extent is smaller than the range of the B-spline grid the number of quadrature points
  !> within each knot interval is expanded by the factor ceiling(n_int/R), where R is the extent of the CGTO (determined to coincide with the nearest larger knot) and n_int is the number of distinct 
  !> intervals of knots in the B-spline basis.
  subroutine radial_grid (centers,only_on_bto_grid,R_min,R_max,bspline_grid,&
                            delta_cgto_grid,x,w,n,n_rng_knot,r_points,weights,n_total_points)
     use sort_gbl, only: sort_float
     implicit none
     type(bspline_grid_obj), intent(in) :: bspline_grid
     integer, intent(in) :: n, n_rng_knot
     real(kind=cfp), intent(in) :: R_min, R_max, x(2*n+1), w(2*n+1), centers(:), delta_cgto_grid
     logical, intent(in) :: only_on_bto_grid
     !OUTPUT variables:
     integer, intent(out) :: n_total_points
     real(kind=cfp), allocatable :: r_points(:), weights(:)

     integer :: n_points, err, interval, cnt, i, n_intervals, n_ranges, j, k, rng, n_int
     real(kind=cfp) :: d, A, B, R_lim, x_AB(2*n+1), w_AB(2*n+1), delta, range_end, range_start, R_start, R_end, center, test
!     real(kind=cfp), parameter :: delta_cgto_grid = 0.20_cfp !1.0_cfp !todo this should be determined adaptively so that in between consecutive end points the CGTO falls down only by a given amount, etc.
     real(kind=cfp), allocatable :: list(:,:), range_start_end(:,:)

         n_points = 2*n+1 !Number of quadrature points within each quadrature interval
         if (R_max <= bspline_grid % A .and. only_on_bto_grid) then
            call xermsg ('cgto_pw_expansions_mod', 'radial_grid', &
                         'CGTO is negligible on the BTO grid. This routine would not work properly.', 1, 1)
         end if
         if (R_min >= bspline_grid % B) then
            call xermsg ('cgto_pw_expansions_mod', 'radial_grid', &
                         'All CGTOs must start below the end point of the B-spline grid.', 2, 1)
         end if

         !write(stdout,'("CGTO negligible beyond radius [a.u.]: ",e25.15)') R_max

         if (only_on_bto_grid) then
            R_start = max(R_min,bspline_grid%A)
         else
            R_start = min(R_min,bspline_grid%A)
         endif

         R_end = min(R_max,bspline_grid%B)

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

         !print *,'input',R_min,R_max
         !print *,'adj',R_start,R_end

         if (delta_cgto_grid <= 0.0_cfp) then
            call xermsg ('cgto_pw_expansions_mod', 'radial_grid', 'On input delta_cgto_grid was .le. 0.0_cfp.', 3, 1)
         end if

         j = ceiling((R_end-R_start)/delta_cgto_grid)
         allocate(list(bspline_grid%no_knots+size(centers)+j+2,1),range_start_end(2,bspline_grid%no_knots+size(centers)+j),stat=err)
         if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','radial_grid','Memory allocation error.',err,1)
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
            k = k + 1
            list(k,1) = centers(i)
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
         n_total_points = n_ranges*n_points*n_rng_knot

         if (allocated(r_points)) deallocate(r_points)
         if (allocated(weights)) deallocate(weights)
         allocate(r_points(n_total_points),weights(n_total_points),stat=err)
         if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','radial_grid','Memory allocation failed.',err,1)

         !Construct quadratures for the individual ranges:
         cnt = 0
         do j=1,n_ranges
            A = range_start_end(1,j)
            delta = (range_start_end(2,j)-range_start_end(1,j))/(n_rng_knot*1.0_cfp)
            do i=1,n_rng_knot
               B = A+delta
               if (i .eq. n_rng_knot) B = range_start_end(2,j)
               !Prepare the quadrature grid for the interval [A,B] expanding the canonical G-L quadrature.
               !print *,'quad',A,B,cnt
               call gl_expand_A_B(x,w,n,r_points(cnt+1:cnt+n_points),weights(cnt+1:cnt+n_points),A,B)
               cnt = cnt + n_points
               A = A + delta
            enddo !i
         enddo
!
         if (cnt .ne. n_total_points) call xermsg('cgto_pw_expansions_mod','radial_grid','Error constructing the radial grid.',3,1)

  end subroutine radial_grid

  subroutine radial_grid_CGTO_pair(cgto_A,cgto_B,R_min,R_max,x,w,n,delta_cgto_grid,r_points,weights,n_total_points)
     use sort_gbl, only: sort_float
     implicit none
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: R_min, R_max, x(2*n+1), w(2*n+1)
     type(CGTO_shell_data_obj), intent(in) :: cgto_A, cgto_B
     !OUTPUT variables:
     integer, intent(out) :: n_total_points
     real(kind=cfp), allocatable :: r_points(:), weights(:)
     real(kind=cfp), intent(in) :: delta_cgto_grid !todo this should be determined adaptively so that in between consecutive end points the CGTO falls down only by a given amount, etc.

     integer :: n_points, err, cnt, i, n_intervals, n_ranges, j, k, p, q
     real(kind=cfp) :: range_end, range_start, R_start, R_end, prod_alp, prod_P(3), RP
     real(kind=cfp), allocatable :: list(:,:), range_start_end(:,:)

         n_points = 2*n+1 !Number of quadrature points within each quadrature interval

         !write(stdout,'("CGTO negligible beyond radius [a.u.]: ",e25.15)') R_max

         !todo determine properly
         R_start = R_min
         R_end = R_max

         j = ceiling((R_max-R_min)/delta_cgto_grid)
         n_ranges = j+cgto_A%number_of_primitives*cgto_B%number_of_primitives+1
         allocate(list(n_ranges,1),range_start_end(2,n_ranges),stat=err)
         if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','radial_grid_CGTO_pair','Memory allocation 1 failed.',err,1)

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
         allocate(r_points(n_total_points),weights(n_total_points),stat=err)
         if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','radial_grid_CGTO_pair','Memory allocation 2 failed.',err,1)

         !Construct quadratures for the individual ranges:
         cnt = 0
         do i=1,n_ranges
           !Prepare the quadrature grid for the interval [range_start,range_end] expanding the canonical G-L quadrature.
           call gl_expand_A_B (x, w, n, r_points(cnt+1:cnt+n_points), weights(cnt+1:cnt+n_points), &
                                range_start_end(1,i), range_start_end(2,i))
           cnt = cnt + n_points
           !print *,'gto',range_start,range_end
         enddo
!
         if (cnt /= n_total_points) then
            call xermsg ('cgto_pw_expansions_mod', 'radial_grid_CGTO_pair', 'Error constructing the radial grid.', 3, 1)
         end if

  end subroutine radial_grid_CGTO_pair

  !> Calculates the partial wave expansion of CGTO up to partial wave L=max_l_inp and on the grid of radial points r_points. From this it generates the double angular projections for X_lm up to l=max_lp.
  !> It uses auxiliary arrays declared on top of this routine.
  subroutine omp_calculate_CGTO_pw_coefficients_analytic (threshold,max_l_inp,max_lp,cgto_shell,r_points,angular_integrals,&
                                                non_neg_indices_l,gaunt_angular_integrals,non_neg_indices_l_lp)
      use phys_const_gbl, only: fourpi, twopi, pi
      use special_functions_gbl, only: cfp_besi, cfp_eval_poly_horner
      use omp_lib
      implicit none
      integer, intent(in) :: max_l_inp, max_lp
      type(CGTO_shell_data_obj), intent(in) :: cgto_shell
      real(kind=cfp), intent(in) :: threshold
      real(kind=cfp), allocatable :: angular_integrals(:,:), gaunt_angular_integrals(:,:)
      integer, allocatable :: non_neg_indices_l_lp(:,:,:), non_neg_indices_l(:,:)
      real(kind=cfp), allocatable :: r_points(:)

      integer, parameter :: kode = 2
      real(kind=cfp), parameter :: half = 0.5_cfp

      integer :: i, j, err, l, m, m_ind, lm, Mg, l_besi, nz, Lg_Mg, max_l, max_l_aux, n_threads, iam, n_cgto_m, n_Xlm, &
                    l_lim, lp, mp, lpp, mpp, p_l, p, p_lp, lpp_mpp, lp_mp, i_thread
      real(kind=cfp) :: r_square, arg, fnu, asym, RA, prec, coupling, max_val, frac !, start_t, end_t
      real(kind=cfp), allocatable :: besi_vals(:,:), contraction_besi(:,:), contraction(:), prim_fac(:,:), besi_args(:,:), &
                                    transl_cfs(:,:,:), c_lambda(:,:,:), Xlm_CGTO_center(:), tmp(:)
      integer :: n_total_points, lambda, CGTO_L
      logical :: is_zero, in_parallel

      real(kind=cfp), allocatable :: angular_integrals_tmp(:,:,:), tmp_lp(:,:), gaunt_angular_integrals_tmp(:,:)
      integer, allocatable :: non_neg_indices_l_tmp(:,:,:), non_neg_indices_l_lp_tmp(:,:,:), p_lp_last(:)
      integer, allocatable :: non_neg_indices_l_lp_owners(:,:,:)

         if (omp_in_parallel()) then
            in_parallel = .true.
         else
            in_parallel = .false.
         endif

         if (max_lp < 0 .or. max_l_inp < 0) then
            call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                         'On input max_lp or max_l_inp were < 0.', 1, 1)
         end if

         prec = F1MACH(4,cfp_dummy)

         max_l = max_l_inp+max_lp !Maximum L that will be needed in the <CGTO|X_{lm}> angular projections to then calculate the double angular projections.
         n_total_points = size(r_points)
         n_cgto_m = 2*cgto_shell%l+1
         n_Xlm = (max_l+1)**2

         if (allocated(angular_integrals)) deallocate(angular_integrals)
         if (allocated(gaunt_angular_integrals)) deallocate(gaunt_angular_integrals)
         if (allocated(non_neg_indices_l)) deallocate(non_neg_indices_l)
         if (allocated(non_neg_indices_l_lp)) deallocate(non_neg_indices_l_lp)

         allocate(non_neg_indices_l(n_cgto_m,n_Xlm),non_neg_indices_l_lp(n_cgto_m,(max_lp+1)**2,(max_l_inp+1)**2),stat=err)
         if (err /= 0) then
            call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                         'Memory allocation 1 failed.', err, 1)
         end if
         non_neg_indices_l = 0
         non_neg_indices_l_lp = 0

         RA = sqrt(dot_product(cgto_shell%center,cgto_shell%center))

         !Calculate the angular integrals trivially if CGTO is at the CMS.
         if (RA .le. prec) then

            !The only projections that are non-zero are those where the real spherical harmonic has the same (l,m) values as the CGTO:
            p_l = n_cgto_m
            if (max_l < cgto_shell%l) p_l = 1 !in this case we just allocate the smallest possible number of elements and exit

            allocate(angular_integrals(n_total_points,p_l),stat=err)
            if (err /= 0) then
                call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                             'Memory allocation 1 failed.', err, 1)
            end if
            angular_integrals = 0.0_cfp

            if (max_l < cgto_shell%l) return !see below

            allocate(contraction(n_total_points),stat=err)
            if (err /= 0) then
                call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                             'Memory allocation 2 failed.', err, 1)
            end if

            do i=1,n_total_points
               !Radial part of the CGTO
               r_square = r_points(i)*r_points(i)
               contraction(i) = 0.0_cfp
               do j=1,cgto_shell%number_of_primitives
                  contraction(i) = contraction(i) &
                                    + cgto_shell%contractions(j)*cgto_shell%norms(j)*exp(-cgto_shell%exponents(j)*r_square)
               enddo
               contraction(i) = contraction(i)*cgto_shell%norm*sqrt(fourpi/(2*cgto_shell%l+1.0_cfp))*r_points(i)**cgto_shell%l
            enddo !i
   
            p_l = 0
            do m=-cgto_shell%l,cgto_shell%l
               lm = cgto_shell%l*cgto_shell%l+cgto_shell%l+m+1
               p_l = p_l + 1
               non_neg_indices_l(cgto_shell%l+m+1,lm) = p_l
               do i=1,n_total_points
                  angular_integrals(i,p_l) = contraction(i)
               enddo !i
            enddo
            !return
       
         else
!
!-------    Use single centre expansion of the CGTO to evaluate the partial wave coefficients.
!           Angular integrals for s-type CGTOs can be always evaluated accurately using a single centre expansion; tests show that for higher-L CGTOs accuracy is preserved as well.
            l_besi = max_l+cgto_shell%l+1
   
            FNU = l_besi-1+half 
            !Compute the smallest value of the Bessel I argument for which the asymptotic expansion can be used
            !todo the asym value is not used anywhere at the moment
            if (cfp .eq. wp) then
               asym = MAX(17.00_cfp,0.550_cfp*FNU*FNU)
            else !quad precision case
               asym = MAX(76.00_cfp,0.550_cfp*FNU*FNU)
            endif 
   
            !Buffer for exponentially scaled modified Bessel values for each CGTO primitive exponent.
            allocate(besi_vals(1:l_besi,1:cgto_shell%number_of_primitives),contraction_besi(0:l_besi-1,n_total_points),stat=err)
            if (err /= 0) then
                call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                             'Memory allocation 3 failed.', err, 1)
            end if
   
            do i=1,n_total_points
               if (r_points(i) .le. prec) then
                  l_lim = 0
                  contraction_besi(:,i) = 0
                  do j=1,cgto_shell%number_of_primitives
                     !take the limit of the modified (exponentially scalled) spherical Bessel function for r->0: only the l=0 term is non-zero
                     besi_vals(1:l_besi,j) = 0.0_cfp
                     besi_vals(1,j) = fourpi
                  enddo !j
               else
                  l_lim = l_besi-1
                  do j=1,cgto_shell%number_of_primitives
                     arg = 2.0_cfp*cgto_shell%exponents(j)*r_points(i)*RA
                     call cfp_besi(arg,half,kode,l_besi,besi_vals(1:l_besi,j),nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=l+1 is needed to get y_{alpha+l}.
                     besi_vals(1:l_besi,j) = twopi*sqrt(pi)*besi_vals(1:l_besi,j)/sqrt(cgto_shell%exponents(j)*r_points(i)*RA)
                  enddo !j
               endif
      
               r_square = (r_points(i)-RA)**2
               do l=0,l_lim
                  !Radial part of the CGTO multiplied by the modified Bessel function
                  contraction_besi(l,i) = 0.0_cfp
                  do j=1,cgto_shell%number_of_primitives
                     contraction_besi(l,i) = contraction_besi(l,i) &
                            + cgto_shell%contractions(j)*cgto_shell%norms(j)*exp(-cgto_shell%exponents(j)*r_square)*besi_vals(l+1,j)
                  enddo !j
                  contraction_besi(l,i) = contraction_besi(l,i)*cgto_shell%norm
               enddo !l
            enddo !i
   
            if (cgto_shell%l .eq. 0) then
   
               !Real spherical harmonics for the nuclei: result in the array Xlm_CGTO_center
               call precalculate_Xlm_for_CGTO_center(cgto_shell%center,max_l,Xlm_CGTO_center)

               !$OMP PARALLEL IF(in_parallel .eqv. .false.) DEFAULT(NONE) PRIVATE(l,m,lm,m_ind,i,iam,tmp,max_val,p_l,err) &
               !$OMP & SHARED(max_l,cgto_shell,n_total_points,contraction_besi,Xlm_CGTO_center,&
               !$OMP & n_threads,angular_integrals_tmp,non_neg_indices_l_tmp,n_Xlm,n_cgto_m,threshold,in_parallel)

               !$OMP SINGLE
               n_threads = omp_get_num_threads()
               p_l = max((max_l + 1)**2 / n_threads + 1, n_threads)

               allocate (angular_integrals_tmp(n_total_points,p_l,n_threads), stat=err)
               if (err /= 0) then
                  call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                               'Memory allocation 41 failed.', err, 1)
               end if

               allocate (non_neg_indices_l_tmp(n_cgto_m,n_Xlm,n_threads), stat=err)
               if (err /= 0) then
                 call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                              'Memory allocation 42 failed.', err, 1)
               end if
               angular_integrals_tmp = 0.0_cfp
               non_neg_indices_l_tmp = 0
               !$OMP END SINGLE

               allocate(tmp(n_total_points),stat=err)
               if (err /= 0) then
                  call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                               'Memory allocation 4a failed.', err, 1)
               end if
               tmp = 0.0_cfp
               !$OMP BARRIER
   
               iam = omp_get_thread_num()
               m_ind = 1
               p_l = 0
               do l=0,max_l
                  do m=-l,l
                     lm = l*l+l+m+1
                     if (mod(lm,n_threads) .ne. iam) cycle !work distribution

                     do i=1,n_total_points
                        tmp(i) = Xlm_CGTO_center(lm)*contraction_besi(l,i)
                     enddo !i

                     max_val = 0.0_cfp
                     max_val = max(max_val,maxval(abs(tmp)))
   
                     !Save only those projections that are non-negligible
                     if (max_val .ge. threshold) then
                        p_l = p_l + 1
                        non_neg_indices_l_tmp(m_ind,lm,iam+1) = p_l
                        angular_integrals_tmp(1:n_total_points,p_l,iam+1) = tmp(1:n_total_points)
                     endif
   
                  enddo !m
               enddo !l
               !$OMP END PARALLEL
   
            else !cgto_shell%l .ne. 0
               max_l_aux = max_l+cgto_shell%l
               !Real spherical harmonics for the nuclei: result in the module array Xlm_CGTO_center
               call precalculate_Xlm_for_CGTO_center(cgto_shell%center,max_l_aux,Xlm_CGTO_center)
   
               !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
               call precalculate_solh_translation_coeffs(cgto_shell%l,RA,Xlm_CGTO_center,transl_cfs)
   
               !start_t = omp_get_wtime()
               !Evaluate the SCE of the CGTO at the radial points
               !$OMP PARALLEL IF(in_parallel .eqv. .false.) DEFAULT(NONE) &
               !$OMP & PRIVATE(l,m,lm,m_ind,Mg,Lg_Mg,i,iam,c_lambda,is_zero,lambda,CGTO_L,tmp,max_val,p_l,err) &
               !$OMP & SHARED(max_l,cgto_shell,n_total_points,l_besi,angular_integrals,r_points,contraction_besi,asym, &
               !$OMP &        prim_fac,besi_args,transl_cfs,Xlm_CGTO_center,n_threads,angular_integrals_tmp, &
               !$OMP &        non_neg_indices_l_tmp,n_Xlm,n_cgto_m,threshold)
   
               !$OMP SINGLE
               n_threads = omp_get_num_threads()
               p_l = ((max_l + 1)**2 / n_threads + 1) * (2*cgto_shell%l + 1)

               allocate (angular_integrals_tmp(n_total_points,p_l,n_threads), stat=err)
               if (err /= 0) then
                  call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                               'Memory allocation 41 failed.', err, 1)
               end if

               allocate (non_neg_indices_l_tmp(n_cgto_m,n_Xlm,n_threads), stat=err)
               if (err /= 0) then
                  call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                               'Memory allocation 42 failed.', err, 1)
               end if

               angular_integrals_tmp = 0.0_cfp
               non_neg_indices_l_tmp = 0
               !$OMP END SINGLE
   
               allocate(tmp(n_total_points),stat=err)
               if (err /= 0) then
                  call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                               'Memory allocation 4a failed.', err, 1)
               end if
               tmp = 0.0_cfp
               !$OMP BARRIER
   
               iam = omp_get_thread_num()
               p_l = 0
               do l=0,max_l
                  do m=-l,l
                     lm = l*l+l+m+1
                     if (mod(lm,n_threads) .ne. iam) cycle !work distribution
                     !Calculate the coupling coefficients needed to evaluate the projection using SCE of GTO.
                     call calculate_lambda_couplings(cgto_shell%l,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
                     m_ind = 0
                     do Mg=-cgto_shell%l,cgto_shell%l
                        m_ind=m_ind+1
                        Lg_Mg = cgto_shell%l*cgto_shell%l+cgto_shell%l+mg+1
   
                        !Skip the symmetry-forbidden projections
                        is_zero = .true.
                        do lambda=0,l+cgto_shell%l
                           do CGTO_L=0,cgto_shell%l
                              if (c_lambda(lambda,CGTO_L,m_ind) .ne. 0.0_cfp) then
                                 is_zero = .false.
                                 exit
                              endif
                           enddo !CGTO_L
                        enddo !lambda
                        if (is_zero) cycle
   
                        do i=1,n_total_points
                           !angular_integrals(i,m_ind,lm) = CGTO_pw_coefficient(r_points(i),l,cgto_shell%l,Mg,c_lambda,contraction_besi(0:l_besi-1,i))
                           tmp(i) = CGTO_pw_coefficient(r_points(i),l,cgto_shell%l,Mg,c_lambda,contraction_besi(0:l_besi-1,i))
                        enddo !i
   
                        max_val = 0.0_cfp
                        max_val = max(max_val,maxval(abs(tmp)))
   
                        !Save only those projections that are non-negligible
                        if (max_val .ge. threshold) then
                           p_l = p_l + 1
                           non_neg_indices_l_tmp(m_ind,lm,iam+1) = p_l
                           angular_integrals_tmp(1:n_total_points,p_l,iam+1) = tmp(1:n_total_points)
                        endif
   
                     enddo !Mg
                  enddo !m
               enddo !l
               !$OMP END PARALLEL
               !end_t = omp_get_wtime()

            endif !cgto_shell%l .eq. 0

            !Merge the data obtained by each thread into one final array.
   
            !Count the total number of non-negligible projections:
            p_l = 0
            do iam=0,n_threads-1
               p_l = p_l + maxval(non_neg_indices_l_tmp(1:n_cgto_m,1:n_Xlm,iam+1))
            enddo !iam
   
            write(level3,'("Non-negligible part of the angular projections: ",f8.3,"%")') p_l/(n_cgto_m*n_Xlm*1.0_cfp)*100.0_cfp
            write(level3,'("Memory (MiB) required: ",f25.10)') n_total_points*max(p_l,1)*cfp_bytes/(Mib*1.0_cfp)
   
            allocate(angular_integrals(n_total_points,max(p_l,1)),stat=err)
            if (err /= 0) then
               call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pw_coefficients_analytic', &
                            'Memory allocation 5 failed.', err, 1)
            end if
            angular_integrals = 0.0_cfp
   
            p_l = 0
            do l=0,max_l
               do m=-l,l
                  lm = l*l+l+m+1
                  iam = mod(lm,n_threads)
                  do m_ind=1,n_cgto_m
                     p = non_neg_indices_l_tmp(m_ind,lm,iam+1)
                     if (p .ne. 0) then
                        p_l = p_l + 1
                        non_neg_indices_l(m_ind,lm) = p_l
                        angular_integrals(1:n_total_points,p_l) = angular_integrals_tmp(1:n_total_points,p,iam+1)
                     endif
                  enddo !m_ind
               enddo !m
            enddo !l

         endif !RA .le. prec

         allocate (non_neg_indices_l_lp_owners, mold = non_neg_indices_l_lp)
         non_neg_indices_l_lp_owners = -1

         !$OMP PARALLEL IF(.not. in_parallel) DEFAULT(none) &
         !$OMP& FIRSTPRIVATE(gaunt_angular_integrals_tmp, tmp_lp) &
         !$OMP& PRIVATE(m, lm, mp, lp_mp, lpp, mpp, coupling, lpp_mpp, m_ind, p_l, i, l, lp, max_val, i_thread, n_threads, p_lp) &
         !$OMP& SHARED(gaunt_angular_integrals, max_l_inp, max_lp, non_neg_indices_l, n_total_points, threshold, p_lp_last, &
         !$OMP&        n_cgto_m, angular_integrals, non_neg_indices_l_lp, cpl, non_neg_indices_l_lp_owners)
         p_lp = 0
         n_threads = omp_get_num_threads()
         i_thread = omp_get_thread_num()
         allocate (gaunt_angular_integrals_tmp(n_total_points, 20))  ! will grow on demand
         allocate (tmp_lp(n_total_points, n_cgto_m))
         !$OMP DO COLLAPSE(2) SCHEDULE(dynamic, 1)
         do l=0,max_l_inp
            do lp=0,max_lp
               do m=-l,l
                  lm = l*l+l+m+1
                  do mp=-lp,lp
                     lp_mp = lp*lp+lp+mp+1

                     tmp_lp = 0.0_cfp
                     do lpp=abs(l-lp),l+lp
                        do mpp=-lpp,lpp
                           coupling = cpl%rgaunt(lp,l,lpp,mp,m,mpp)
                           if (coupling .ne. 0.0_cfp) then
                              lpp_mpp = lpp*lpp+lpp+mpp+1
                              do m_ind=1,n_cgto_m

                                 p_l = non_neg_indices_l(m_ind,lpp_mpp)
                                 if (p_l .eq. 0) cycle

                                 do i=1,n_total_points
                                    tmp_lp(i,m_ind) = tmp_lp(i,m_ind) + coupling*angular_integrals(i,p_l)
                                 enddo !i

                              enddo !m_ind
                           endif
                        enddo !mpp
                     enddo !lpp

                     do m_ind=1,n_cgto_m
                        max_val = 0.0_cfp
                        max_val = max(max_val,maxval(abs(tmp_lp(1:n_total_points,m_ind))))
                        if (max_val .ge. threshold) then
                           p_lp = p_lp + 1
                           if (p_lp > size(gaunt_angular_integrals_tmp,2)) then
                              i = nint(p_lp*1.5_cfp) !allocate a bit more space
                              !print *,'resize needed',size(gaunt_angular_integrals_tmp,2),i
                              call resize_copy_2d_array(gaunt_angular_integrals_tmp,n_total_points,i)
                           endif
                           non_neg_indices_l_lp_owners(m_ind,lp_mp,lm) = i_thread
                           non_neg_indices_l_lp(m_ind,lp_mp,lm) = p_lp  ! position within the thread-private buffer
                           gaunt_angular_integrals_tmp(1:n_total_points,p_lp) = tmp_lp(1:n_total_points,m_ind)
                        endif
                     enddo !m_ind
      
                  enddo !mp
               enddo !m
            enddo !lp
         enddo !l
         !$OMP SINGLE
         allocate (p_lp_last(n_threads))
         !$OMP END SINGLE
         p_lp_last(i_thread + 1) = p_lp
         !$OMP BARRIER
         !$OMP SINGLE
         do i = 2, n_threads
            p_lp_last(i) = p_lp_last(i) + p_lp_last(i - 1)   ! running sum; global index of the last column owned by a thread
         end do
         allocate (gaunt_angular_integrals(n_total_points, p_lp_last(n_threads)))
         !$OMP END SINGLE
         where (non_neg_indices_l_lp_owners == i_thread)
            non_neg_indices_l_lp = non_neg_indices_l_lp + (p_lp_last(i_thread + 1) - p_lp)   ! shift local to global indices
         end where
         if (p_lp > 0) then
            gaunt_angular_integrals(:, p_lp_last(i_thread + 1) - p_lp + 1 : p_lp_last(i_thread + 1)) &
               = gaunt_angular_integrals_tmp(:, 1:p_lp)
         end if
         !$OMP END PARALLEL
         p_lp = size(gaunt_angular_integrals, 2)

         write(level3,'("Non-negligible part of the double angular projections: ",f8.3,"%")') &
                                        p_lp/(n_cgto_m*(max_lp+1)**2*(max_l_inp+1)**2*1.0_cfp)*100.0_cfp
         write(level3,'("Memory (MiB) required: ",f25.10)') n_total_points*p_lp*cfp_bytes/(Mib*1.0_cfp)
         flush(stdout)   

  end subroutine omp_calculate_CGTO_pw_coefficients_analytic

  !> Calculates the partial wave expansion of a product of two CGTOs up to partial wave L=max_l and on the grid of radial points r_points.
  subroutine omp_calculate_CGTO_pair_pw_coefficients_analytic(max_l,cgto_A,cgto_B,r_points,angular_integrals)
      use phys_const_gbl, only: fourpi, twopi, pi
      use special_functions_gbl, only: cfp_besi, cfp_eval_poly_horner
      use omp_lib
      implicit none
      integer, intent(in) :: max_l
      type(CGTO_shell_data_obj), intent(in) :: cgto_A, cgto_B
      real(kind=cfp), allocatable :: angular_integrals(:,:,:)
      real(kind=cfp), allocatable :: r_points(:)

      integer, parameter :: kode = 2
      real(kind=cfp), parameter :: half = 0.5_cfp

      integer :: i, j, err, l, m, lm, l_besi, nz, max_l_aux, n_threads, iam, n_cgto_A_m, n_cgto_B_m, n_Xlm, p, q, ij, &
                 min_l, max_l_adj
      real(kind=cfp) :: r_square, arg, fnu, asym, term, R_A, R_B, R_AB(3), R_AB_square, fac, K_pref, prec !, start_t, end_t
      real(kind=cfp), allocatable :: besi_vals(:,:), contraction_besi(:,:,:), contraction(:), prim_fac(:,:), &
                                     besi_args(:,:), transl_cfs(:,:,:), c_lambda(:,:,:), prod_alp(:)
      real(kind=cfp), allocatable :: prod_P(:,:), RP(:), exp_fac(:), prod_contr(:), Xlm_CGTO_A_center(:), Xlm_CGTO_B_center(:), &
                                     Xlm_product_CGTO_center(:), transl_cfs_AB(:,:,:,:,:), c_pair_lambda(:,:,:,:,:)
      integer :: n_total_points, n_contr_pairs, Mg_A, Mg_B, Mg_A_ind, BA_ind
      logical :: in_parallel
      !UNCOMMENT BELOW TO DEBUG
!      type(Xlm_x_pair_cgto_surface) :: xlm_pair_cgto_r
!      real(kind=cfp), parameter :: th_min = 0.0_cfp, th_max = pi, phi_min = 0.0_cfp, phi_max = twopi
!      real(kind=cfp) :: SH_A(-cgto_A%l:cgto_A%l,0:cgto_A%l+1),SH_B(-cgto_B%l:cgto_B%l,0:cgto_B%l+1)

         if (omp_in_parallel()) then
            in_parallel = .true.
         else
            in_parallel = .false.
         endif

         prec = F1MACH(4,cfp_dummy)

         if (allocated(angular_integrals)) deallocate(angular_integrals)
         n_total_points = size(r_points)
         n_cgto_A_m = 2*cgto_A%l+1
         n_cgto_B_m = 2*cgto_B%l+1
         n_Xlm = (max_l+1)**2
         allocate(angular_integrals(n_total_points,n_cgto_B_m*n_cgto_A_m,n_Xlm),stat=err)
         if (err /= 0) then
            call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pair_pw_coefficients_analytic', &
                         'Memory allocation 1 failed.', err, 1)
         end if

         angular_integrals = 0.0_cfp

         R_A = sqrt(dot_product(cgto_A%center,cgto_A%center))
         R_B = sqrt(dot_product(cgto_B%center,cgto_B%center))

         !Calculate the angular integrals trivially if both CGTOs sit at the CMS.
         !todo OpenMP parallelize this section
         if (R_A .le. prec .and. R_B .le. prec) then

            allocate(contraction(n_total_points),stat=err)
            if (err /= 0) then
                call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pair_pw_coefficients_analytic', &
                             'Memory allocation 2 failed.', err, 1)
            end if

            do i=1,n_total_points
               !Radial part of the CGTO
               r_square = r_points(i)*r_points(i)
               contraction(i) = 0.0_cfp
               do p=1,cgto_A%number_of_primitives
                  fac = cgto_A%contractions(p)*cgto_A%norms(p)
                  do q=1,cgto_B%number_of_primitives
                     contraction(i) = contraction(i) + fac*cgto_B%contractions(q)*cgto_B%norms(q) &
                                        *exp(-(cgto_A%exponents(p)+cgto_B%exponents(q))*r_square)
                  enddo !k
               enddo !j
               contraction(i) = contraction(i)*cgto_A%norm*cgto_B%norm*fourpi &
                                    /sqrt((2*cgto_A%l+1.0_cfp)*(2*cgto_B%l+1.0_cfp))*r_points(i)**(cgto_A%l+cgto_B%l)
            enddo !i
   
            !The only projections that are non-zero are those for which the real Gaunt coefficient <l,m|L_A,M_A|L_B,M_B> is non-zero:
            min_l = abs(cgto_A%l-cgto_B%l)
            max_l_adj = cgto_A%l+cgto_B%l
            !if (max_l_adj > max_l) call xermsg('cgto_pw_expansions_mod','omp_calculate_CGTO_pair_pw_coefficients_analytic',&
            !&'On input max_l was too small to allow correct evaluation of pw expansion for a pair of CGTOs sitting on CMS.',1,1)
            max_l_adj = min(max_l_adj,max_l) !ensure we don't overflow if max_l is small
            do l=min_l,max_l_adj
               if (mod(l+cgto_A%l+cgto_B%l,2) .ne. 0) cycle !selection rule for Gaunt coefficients
               do m=-l,l
                  lm = l*l+l+m+1
                  do Mg_A=-cgto_A%l,cgto_A%l
                     Mg_A_ind = Mg_A+cgto_A%l+1
                     do Mg_B=-cgto_B%l,cgto_B%l
                        term = cpl%rgaunt(l,cgto_A%l,cgto_B%l,m,Mg_A,Mg_B)
                        if (term .eq. 0.0_cfp) cycle
                        BA_ind = Mg_B+cgto_B%l+1 + n_cgto_B_m*(Mg_A_ind-1)
                        do i=1,n_total_points
                           angular_integrals(i,BA_ind,lm) = term*contraction(i)
                        enddo !i
                     enddo !Mg_B
                  enddo !Mg_A
               enddo !m
            enddo !l
            !goto 1 !UNCOMMENT THIS LINE TO DEBUG
            return

         endif
!
!------- Use single centre expansion of the CGTO to evaluate the partial wave coefficients.
!        Angular integrals for s-type CGTOs can be always evaluated accurately using a single centre expansion; tests show that for higher-L CGTOs accuracy is preserved too.
         min_l = 0
         max_l_adj = max_l
         l_besi = max_l_adj+cgto_A%l+cgto_B%l+1

         FNU = l_besi-1+half 
         !Compute the smallest value of the Bessel I argument for which the asymptotic expansion can be used
         !todo the asym value is not used anywhere at the moment
         if (cfp .eq. wp) then
            asym = MAX(17.00_cfp,0.550_cfp*FNU*FNU)
         else !quad precision case
            asym = MAX(76.00_cfp,0.550_cfp*FNU*FNU)
         endif 

         n_contr_pairs = cgto_A%number_of_primitives*cgto_B%number_of_primitives

         !Buffer for exponentially scaled modified Bessel values for each CGTO primitive exponent.
         allocate(besi_vals(1:l_besi,1:n_contr_pairs), &
                  contraction_besi(n_contr_pairs,0:l_besi-1,n_total_points), &
                  prod_P(3,n_contr_pairs),prod_alp(n_contr_pairs),&
                  RP(n_contr_pairs),exp_fac(n_contr_pairs),prod_contr(n_contr_pairs),stat=err)
         if (err /= 0) then
            call xermsg ('cgto_pw_expansions_mod', 'omp_calculate_CGTO_pair_pw_coefficients_analytic', &
                         'Memory allocation 3 failed.', err, 1)
         end if

         R_AB(1:3) = cgto_A%center(1:3)-cgto_B%center(1:3)
         R_AB_square = dot_product(R_AB,R_AB)
         fac = cgto_A%norm*cgto_B%norm
         ij = 0
         do p=1,cgto_A%number_of_primitives
            do q=1,cgto_B%number_of_primitives
               ij = ij + 1
               prod_alp(ij) = cgto_A%exponents(p)+cgto_B%exponents(q) !exponent of the product GTO
               prod_P(1:3,ij) = (cgto_A%exponents(p)*cgto_A%center(1:3)+cgto_B%exponents(q)*cgto_B%center(1:3))/prod_alp(ij) !center of the product GTO
               K_pref = exp(-cgto_A%exponents(p)*cgto_B%exponents(q)/prod_alp(ij)*R_AB_square) !exponential prefactor for the product GTO
               RP(ij) = sqrt(dot_product(prod_P(1:3,ij),prod_P(1:3,ij)))
               if (RP(ij) .le. prec) RP(ij) = 0.0_cfp !regard this pair as sitting on CMS
               prod_contr(ij) = K_pref*fac*cgto_A%contractions(p)*cgto_A%norms(p)*cgto_B%contractions(q)*cgto_B%norms(q)
            enddo !q
         enddo !p

         do i=1,n_total_points
            do ij=1,n_contr_pairs
               if (RP(ij) .le. prec) then !the product GTO is centered on CMS
                  !take the limit of the modified (exponentially scalled) spherical Bessel function for r->0: only the l=0 term is non-zero
                  besi_vals(1:l_besi,ij) = 0.0_cfp
                  besi_vals(1,ij) = fourpi
               else
                  arg = 2.0_cfp*prod_alp(ij)*r_points(i)*RP(ij)
                  call cfp_besi(arg,half,kode,l_besi,besi_vals(1:l_besi,ij),nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=l+1 is needed to get y_{alpha+l}.
                  besi_vals(1:l_besi,ij) = twopi*sqrt(pi)*besi_vals(1:l_besi,ij)/sqrt(prod_alp(ij)*r_points(i)*RP(ij))
               endif
            enddo !ij

            ij = 0
            do p=1,cgto_A%number_of_primitives
               do q=1,cgto_B%number_of_primitives
                  ij = ij + 1
                  r_square = (r_points(i)-RP(ij))**2
                  exp_fac(ij) = prod_contr(ij)*exp(-prod_alp(ij)*r_square)
               enddo !q
            enddo !p

            do l=0,l_besi-1
               !Radial parts of the product GTOs multiplied by the modified Bessel functions
               ij = 0
               do p=1,cgto_A%number_of_primitives
                  do q=1,cgto_B%number_of_primitives
                     ij = ij + 1
                     contraction_besi(ij,l,i) = exp_fac(ij)*besi_vals(l+1,ij)
                  enddo !q
               enddo !p
            enddo !l
         enddo !i

         !Real spherical harmonics for the CGTO A center: result in the module array Xlm_CGTO_A_center
         call precalculate_Xlm_for_CGTO_center(cgto_A%center,max_l_adj+cgto_A%l,Xlm_CGTO_A_center)

         !Real spherical harmonics for the CGTO B center: result in the module array Xlm_CGTO_B_center
         call precalculate_Xlm_for_CGTO_center(cgto_B%center,max_l_adj+cgto_B%l,Xlm_CGTO_B_center)

         max_l_aux = max_l_adj+cgto_A%l+cgto_B%l
         call precalculate_Xlm_for_CGTO_product_center(n_contr_pairs,prod_P,max_l_aux,Xlm_product_CGTO_center)

         !Precalculate the coefficients in the translation formula for the pair of solid harmonics sitting on centers A,B: this requires Xlm_CGTO_center
         call precalculate_pair_solh_translation_coeffs(cgto_A%l,R_A,Xlm_CGTO_A_center,cgto_B%l,R_B,Xlm_CGTO_B_center,transl_cfs_AB)

         !start_t = omp_get_wtime()
         !Evaluate the SCE of the CGTO at the radial points
         !$OMP PARALLEL IF(in_parallel .eqv. .false.) DEFAULT(NONE) &
         !$OMP & PRIVATE(l,m,lm,Mg_A,Mg_A_ind,Mg_B,BA_ind,i,n_threads,iam,c_pair_lambda) &
         !$OMP & SHARED(n_cgto_B_m,min_l,max_l_adj,cgto_A,cgto_B,n_total_points,l_besi,angular_integrals,r_points, &
         !$OMP &        contraction_besi,n_contr_pairs,Xlm_product_CGTO_center,transl_cfs_AB,in_parallel)
         n_threads = omp_get_num_threads()
         iam = omp_get_thread_num()
         do l=min_l,max_l_adj
            do m=-l,l
               lm = l*l+l+m+1
               if (mod(lm,n_threads) .ne. iam) cycle !work distribution
               !Calculate the coupling coefficients needed to evaluate the projection using SCE of GTO.
               call calculate_pair_lambda_couplings(l,m,cgto_A%l,cgto_B%l,n_contr_pairs,&
                                                    Xlm_product_CGTO_center,transl_cfs_AB,c_pair_lambda)
               do Mg_A=-cgto_A%l,cgto_A%l
                  Mg_A_ind = Mg_A+cgto_A%l+1
                  do Mg_B=-cgto_B%l,cgto_B%l
                     BA_ind = Mg_B+cgto_B%l+1 + n_cgto_B_m*(Mg_A_ind-1)
                     do i=1,n_total_points
                        angular_integrals(i,BA_ind,lm) = &
                        &CGTO_pair_pw_coefficient(r_points(i),l,cgto_A%l,Mg_A,cgto_B%l,Mg_B,c_pair_lambda,&
                                                  n_contr_pairs,contraction_besi(1:n_contr_pairs,0:l_besi-1,i))
                        !write(111,'(e25.15,3i,e25.15)') r_points(i),Mg_B,Mg_A,lm,angular_integrals(i,BA_ind,lm)
                     enddo !i
                  enddo !Mg_B
               enddo !Mg_A
            enddo !m
         enddo !l
         !$OMP END PARALLEL
         !end_t = omp_get_wtime()

         !UNCOMMENT BELOW TO DEBUG
!   1     xlm_pair_cgto_r%cgto_A = cgto_A
!         xlm_pair_cgto_r%cgto_B = cgto_B
!         print *,'r=',r_points(10)
!         do l=min_l,max_l_adj
!         xlm_pair_cgto_r%l=l
!         do m=-l,l
!         xlm_pair_cgto_r%m = m
!         do Mg_A=-cgto_A%l,cgto_A%l
!         do Mg_B=-cgto_B%l,cgto_B%l
!         xlm_pair_cgto_r%neval = 0; xlm_pair_cgto_r%ndiv = 0
!         xlm_pair_cgto_r%cgto_A_m = Mg_A !m value of the CGTO A
!         xlm_pair_cgto_r%cgto_B_m = Mg_B !m value of the CGTO B
!         xlm_pair_cgto_r%r = r_points(10)
!         lm  = xlm_pair_cgto_r%l*xlm_pair_cgto_r%l+xlm_pair_cgto_r%l+xlm_pair_cgto_r%m+1
!         BA_ind = Mg_B+cgto_B%l+1 + n_cgto_B_m*(Mg_A+cgto_A%l+1-1)
!         term = quad2d(xlm_pair_cgto_r,th_min,th_max,phi_min,phi_max,10e-10_cfp)
!         write(*,'("numeric",2i4,2i4,3e25.15)') l,m,Mg_B,Mg_A,term,angular_integrals(10,BA_ind,lm),term/angular_integrals(10,BA_ind,lm)
!         enddo
!         enddo
!         enddo
!         enddo

  end subroutine omp_calculate_CGTO_pair_pw_coefficients_analytic

  !> Evaluates the product of CGTO and a real spherical harmonic at a given point in space. x,y are the spherical polar coordinates on the sphere.
  function eval_Xlm_x_pair_cgto_surface(this,x,y)
      use phys_const_gbl, only: fourpi
      implicit none
      class(Xlm_x_pair_cgto_surface) :: this
      real(kind=cfp), intent(in) :: x, y
      real(kind=cfp) :: eval_Xlm_x_pair_cgto_surface

      real(kind=cfp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
      real(kind=cfp) :: xc,yc,zc,s, r_A(3), RH(-this%l:this%l,0:this%l+1)

         s = sin(x) !sin(theta)
         xc = s*cos(y)
         yc = s*sin(y)
         zc = cos(x)

         r_A(1:3) = (/xc,yc,zc/)
         if (this%l > 0) then
            call resh(RH,xc,yc,zc,this%l)
         else
            RH(0,0) = norm
         endif

         eval_Xlm_x_pair_cgto_surface = gto_pair_eval_R(this%cgto_A,this%cgto_A_m,this%cgto_B,&
                                                        this%cgto_B_m,this%r,r_A)*RH(this%m,this%l)*s

         this%neval = this%neval + 1

  end function eval_Xlm_x_pair_cgto_surface

  !> Evaluates the product of a pair of CGTOs for the given m values at radial distance x and for angular direction given by vector R.
  function gto_pair_eval_R(cgto_A,m_A,cgto_B,m_B,x,R)
      implicit none
      type(CGTO_shell_data_obj), intent(in) :: cgto_A, cgto_B
      real(kind=cfp), intent(in) :: R(3), x
      integer, intent(in) :: m_A, m_B
      real(kind=cfp) :: gto_pair_eval_R

      real(kind=cfp) :: sum_exp_A, sum_exp_B, r_A(3), r_B(3), r_A_square, r_B_square
      real(kind=cfp) :: SH_A(-cgto_A%l:cgto_A%l,0:cgto_A%l+1), SH_B(-cgto_B%l:cgto_B%l,0:cgto_B%l+1)
      integer :: j

         r_A(1:3) = (/x*R(1)-cgto_A%center(1),x*R(2)-cgto_A%center(2),x*R(3)-cgto_A%center(3)/)
         r_B(1:3) = (/x*R(1)-cgto_B%center(1),x*R(2)-cgto_B%center(2),x*R(3)-cgto_B%center(3)/)

         if (cgto_A%l > 0) then
            call solh(SH_A,r_A(1),r_A(2),r_A(3),cgto_A%l)
         else
            SH_A(0,0) = 1.0_cfp
         endif

         if (cgto_B%l > 0) then
            call solh(SH_B,r_B(1),r_B(2),r_B(3),cgto_B%l)
         else
            SH_B(0,0) = 1.0_cfp
         endif

         r_A_square = dot_product(r_A,r_A)
         sum_exp_A = 0.0_cfp
         do j=1,cgto_A%number_of_primitives
            sum_exp_A = sum_exp_A + cgto_A%contractions(j)*cgto_A%norms(j)*exp(-cgto_A%exponents(j)*r_A_square)
         enddo
         sum_exp_A = cgto_A%norm*sum_exp_A

         r_B_square = dot_product(r_B,r_B)
         sum_exp_B = 0.0_cfp
         do j=1,cgto_B%number_of_primitives
            sum_exp_B = sum_exp_B + cgto_B%contractions(j)*cgto_B%norms(j)*exp(-cgto_B%exponents(j)*r_B_square)
         enddo
         sum_exp_B = cgto_B%norm*sum_exp_B

         gto_pair_eval_R = SH_A(m_A,cgto_A%l)*sum_exp_A*SH_B(m_B,cgto_B%l)*sum_exp_B

  end function gto_pair_eval_R

  function CGTO_pair_pw_coefficient(r,l,Lg_A,Mg_A,Lg_B,Mg_B,c_pair_lambda,n_contr_pairs,contraction_besi)
      use special_functions_gbl, only: cfp_eval_poly_horner
     implicit none
     integer, intent(in) :: l,Lg_A,Mg_A,Lg_B,Mg_B,n_contr_pairs
     real(kind=cfp), intent(in) :: r, c_pair_lambda(1:,0:,0:,:,:), contraction_besi(1:,0:)
     real(kind=cfp) :: CGTO_pair_pw_coefficient

     integer :: lp, Mg_A_ind, Mg_B_ind, Lg, lambda, ij
     real(kind=cfp) :: lambda_cf(0:Lg_A+Lg_B+1), contr

        Mg_A_ind = Lg_A+Mg_A+1
        Mg_B_ind = Lg_B+Mg_B+1
        Lg = Lg_A+Lg_B

        do lp=0,Lg
           lambda_cf(lp) = 0.0_cfp
           do lambda=max(0,l-lp),l+lp !here lp = lap+lbp, but c_pair_lambda for each lp is obtained from various combinations of (lap,lbp) for which lap+lbp=lp
              contr = 0.0_cfp
              do ij=1,n_contr_pairs
                 contr = contr + c_pair_lambda(ij,lambda,lp,Mg_B_ind,Mg_A_ind)*contraction_besi(ij,lambda)
              enddo !ij
              lambda_cf(lp) = lambda_cf(lp) + contr
           enddo !lambda
        enddo !lp

        !Use the Horner scheme to evaluate the polynomial in r: this shouldn't be necessary since cgto%l is typically low but it doesn't hurt to do it accurately.
        CGTO_pair_pw_coefficient = cfp_eval_poly_horner(Lg,r,lambda_cf(0:Lg+1))

  end function CGTO_pair_pw_coefficient

  !> \warning Requires precalculated values of the real spherical harmonics at the position of the CGTO nucleus. The coupling coefficients should also be precalculated for performance reasons.
  subroutine precalculate_pair_solh_translation_coeffs(CGTO_A_L,RA_A,Xlm_CGTO_A_center,&
                                                       CGTO_B_L,RA_B,Xlm_CGTO_B_center,transl_cfs_AB)
     use phys_const_gbl, only: fourpi
     implicit none
     real(kind=cfp), allocatable :: Xlm_CGTO_A_center(:), Xlm_CGTO_B_center(:), transl_cfs_AB(:,:,:,:,:)
     integer, intent(in) :: CGTO_A_L, CGTO_B_L
     real(kind=cfp), intent(in) :: RA_A, RA_B

     integer :: n_mp, err, CGTO_M, lp, mp, max_lp, lp_min, lp_max, CGTO_A_M, CGTO_A_M_ind, CGTO_B_M, CGTO_B_M_ind
     integer :: la_p, lb_p, la_p_lb_p, ma_p, mb_p
     real(kind=cfp), allocatable :: transl_cfs_A(:,:,:), transl_cfs_B(:,:,:)
   
        !Translation coefficients for the individual CGTOs
        call precalculate_solh_translation_coeffs(CGTO_A_L,RA_A,Xlm_CGTO_A_center,transl_cfs_A)
        call precalculate_solh_translation_coeffs(CGTO_B_L,RA_B,Xlm_CGTO_B_center,transl_cfs_B)

        max_lp = CGTO_A_L+CGTO_B_L
        n_mp = 2*max_lp+1

        if (allocated(transl_cfs_AB)) deallocate(transl_cfs_AB)
        allocate(transl_cfs_AB(n_mp,0:max(max_lp,1),0:max(max_lp,1),2*CGTO_B_L+1,2*CGTO_A_L+1),stat=err)
        if (err /= 0) then
            call xermsg ('cgto_pw_expansions_mod', 'precalculate_pair_solh_translation_coeffs', 'Memory allocation failed.', err, 1)
        end if

        transl_cfs_AB = 0.0_cfp

        do CGTO_A_M = -CGTO_A_L,CGTO_A_L
           CGTO_A_M_ind = CGTO_A_M+CGTO_A_L+1  
           do CGTO_B_M = -CGTO_B_L,CGTO_B_L
              CGTO_B_M_ind = CGTO_B_M+CGTO_B_L+1
              do la_p=0,CGTO_A_L
                 do lb_p=0,CGTO_B_l
                    la_p_lb_p = la_p+lb_p
                    lp_min = abs(la_p-lb_p)
                    lp_max = la_p+lb_p
                    do ma_p=-la_p,la_p
                       do mb_p=-lb_p,lb_p
                          do lp=lp_min,lp_max
                             if (mod(lp+la_p+lb_p,2) .ne. 0) cycle !selection rule for Gaunt coefficients
                             do mp=-lp,lp
                                transl_cfs_AB(mp+lp+1,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) = &
                                    transl_cfs_AB(mp+lp+1,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) + &
                                    transl_cfs_A(ma_p+la_p+1,la_p,CGTO_A_M_ind) &
                                    *transl_cfs_B(mb_p+lb_p+1,lb_p,CGTO_B_M_ind)*cpl%rgaunt(lp,la_p,lb_p,mp,ma_p,mb_p)
                             enddo !mp
                          enddo !lp
                       enddo !mb_p
                    enddo !ma_p
                 enddo !lb_p
              enddo !la_p
           enddo !CGTO_B_M
        enddo !CGTO_A_M

  end subroutine precalculate_pair_solh_translation_coeffs

  subroutine calculate_pair_lambda_couplings(l,m,CGTO_A_L,CGTO_B_L,n_contr_pairs,Xlm_product_CGTO_center,&
                                                transl_cfs_AB,c_pair_lambda)
     implicit none
     real(kind=cfp), allocatable :: transl_cfs_AB(:,:,:,:,:), Xlm_product_CGTO_center(:)
     integer, intent(in) :: CGTO_A_L, CGTO_B_L, l, m, n_contr_pairs
     !OUTPUT:
     real(kind=cfp), allocatable :: c_pair_lambda(:,:,:,:,:)

     real(kind=cfp) :: coupling 
     integer :: CGTO_M, CGTO_M_ind, lm, lp, mp, lambda, mu, base, lambda_max, err, d1, d2, d3, d4, d5, mp_ind
     integer :: CGTO_A_M, CGTO_B_M, la_p_lb_p, ij, CGTO_A_M_ind, CGTO_B_M_ind

        lambda_max = l+CGTO_A_L+CGTO_B_L
        d1 = n_contr_pairs
        d2 = max(lambda_max,1)
        d3 = max(CGTO_A_L+CGTO_B_L,1)
        d4 = 2*CGTO_B_L+1
        d5 = 2*CGTO_A_L+1

        if (.not. allocated(c_pair_lambda) .or. &
            ubound(c_pair_lambda,1) < d1 .or. ubound(c_pair_lambda,2) < d2 .or. &
            ubound(c_pair_lambda,3) < d3 .or. ubound(c_pair_lambda,4) < d4 .or. &
            ubound(c_pair_lambda,5) < d5) then
           if (allocated(c_pair_lambda)) deallocate(c_pair_lambda)
           allocate(c_pair_lambda(1:d1,0:d2,0:d3,1:d4,1:d5),stat=err)
           if (err .ne. 0) then
              print *,d1,d2,d3,d4,d5
              call xermsg('cgto_pw_expansions_mod','calculate_pair_lambda_couplings','Memory allocation failed.',err,1)
           endif
        endif

        c_pair_lambda = 0.0_cfp
        do lp=0,CGTO_A_L+CGTO_B_L
           do lambda=abs(l-lp),l+lp
              if (mod(lambda+lp+l,2) .ne. 0) cycle !selection rule for Gaunt coefficients
              do mp=-lp,lp
                 mp_ind = lp+mp+1
                 do mu=-lambda,lambda
                    coupling = cpl%rgaunt(l,lp,lambda,m,mp,mu)
                    if (coupling .eq. 0.0_cfp) cycle
                    base = n_contr_pairs*(lambda*lambda+lambda+mu)
                    do CGTO_A_M = -CGTO_A_L,CGTO_A_L
                       CGTO_A_M_ind = CGTO_A_M+CGTO_A_L+1
                       do CGTO_B_M = -CGTO_B_L,CGTO_B_L
                          CGTO_B_M_ind = CGTO_B_M+CGTO_B_L+1
                          do la_p_lb_p=0,CGTO_A_L+CGTO_B_L
                             do ij=1,n_contr_pairs
                                c_pair_lambda(ij,lambda,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) = &
                                    c_pair_lambda(ij,lambda,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind)&
                                    & + coupling*transl_cfs_AB(mp_ind,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind)&
                                        * Xlm_product_CGTO_center(ij+base)
                             enddo !ij
                          enddo !la_p_lb_p
                       enddo !CGTO_B_M
                    enddo !CGTO_A_M
                 enddo !mu
              enddo !mp
           enddo !lambda
        enddo !lp

  end subroutine calculate_pair_lambda_couplings

  function CGTO_pw_coefficient(r,l,Lg,Mg,c_lambda,contraction_besi)
      use special_functions_gbl, only: cfp_eval_poly_horner
     implicit none
     integer, intent(in) :: l,Lg,Mg
     real(kind=cfp), intent(in) :: r, c_lambda(0:,0:,:), contraction_besi(0:)
     real(kind=cfp) :: CGTO_pw_coefficient

     integer :: lp, Mg_ind
     real(kind=cfp) :: lambda_cf(0:Lg+1)

        Mg_ind = Lg+Mg+1
        do lp=0,Lg
           lambda_cf(lp) = sum(c_lambda(abs(l-lp):(l+lp),lp,Mg_ind)*contraction_besi(abs(l-lp):(l+lp)))
        enddo !lp

        !Use the Horner scheme to evaluate the polynomial in r: this shouldn't be necessary since cgto%l is typically low but it doesn't hurt to do it accurately.
        CGTO_pw_coefficient = cfp_eval_poly_horner(Lg,r,lambda_cf(0:Lg+1))

  end function CGTO_pw_coefficient

  function CGTO_pw_coefficient_stable(asym,r,l,lm,Lg,Lg_Mg,c_lambda,contraction_besi,prim_fac,n_prim,besi_args)
      use special_functions_gbl, only: cfp_eval_poly_horner
     implicit none
     integer, intent(in) :: l,Lg,lm,Lg_Mg,n_prim
     real(kind=cfp), intent(in) :: asym, r, c_lambda(0:,0:,:,:), contraction_besi(0:), prim_fac(n_prim), besi_args(n_prim)
     real(kind=cfp) :: CGTO_pw_coefficient_stable

     integer :: lp, j
     real(kind=cfp) :: lambda_cf(0:Lg+1), contraction(0:max(l+Lg,1))

        if (all(besi_args(:) .ge. asym)) then
           do lp=0,Lg
              lambda_cf(lp) = 0.0_cfp
              do j=1,n_prim
                 lambda_cf(lp) = lambda_cf(lp) + prim_fac(j) &
                                * sum_besi_half_asym_lambda_cf(besi_args(j),c_lambda(0:l+lp,lp,lm,Lg_Mg),l+lp,abs(l-lp),l+lp)
              enddo !j
           enddo !lp
        else
           do lp=0,Lg
              lambda_cf(lp) = sum(c_lambda(abs(l-lp):(l+lp),lp,lm,Lg_Mg)*contraction_besi(abs(l-lp):(l+lp)))
           enddo !lp
        endif

        !Use the Horner scheme to evaluate the polynomial in r: this shouldn't be necessary since cgto%l is typically low but it doesn't hurt to do it accurately.
        CGTO_pw_coefficient_stable = cfp_eval_poly_horner(Lg,r,lambda_cf(0:Lg+1))

  end function CGTO_pw_coefficient_stable

  !> \warning Requires precalculated values of the real spherical harmonics at the position of the CGTO nucleus. The coupling coefficients should also be precalculated for performance reasons.
  subroutine precalculate_solh_translation_coeffs(CGTO_L,RA,Xlm_CGTO_center,transl_cfs)
     use phys_const_gbl, only: fourpi
     implicit none
     real(kind=cfp), allocatable :: transl_cfs(:,:,:), Xlm_CGTO_center(:)
     integer, intent(in) :: CGTO_L
     real(kind=cfp), intent(in) :: RA

     integer :: n_mp, err, CGTO_M, lp, lpp, mp, mpp, ind
     real(kind=cfp) :: fac, sum_mpp
   
        n_mp = 2*CGTO_L+1

        if (allocated(transl_cfs)) deallocate(transl_cfs)
        allocate(transl_cfs(n_mp,0:max(CGTO_L,1),n_mp),stat=err)
        if (err /= 0) then
            call xermsg ('cgto_pw_expansions_mod', 'precalculate_solh_translation_coeffs', 'Memory allocation failed.', err, 1)
        end if

        fac = sqrt(fourpi/(2*CGTO_L+1.0_cfp))
        transl_cfs = 0.0_cfp
        do CGTO_M=-CGTO_L,CGTO_L
           do lp=0,CGTO_L
              lpp = CGTO_L-lp
              do mp=-lp,lp
                 sum_mpp = 0.0_cfp
                 do mpp=-lpp,lpp
                    ind = lpp*lpp+lpp+mpp+1
                    sum_mpp = sum_mpp + cpl%G_real_cf(CGTO_L,lp,CGTO_M,mp,mpp)*Xlm_CGTO_center(ind)
                 enddo !mpp
                 transl_cfs(mp+lp+1,lp,CGTO_M+CGTO_L+1) = sum_mpp*fac*(-1)**(lpp)*RA**lpp
              enddo !mp
           enddo !lp
        enddo !CGTO_M

  end subroutine precalculate_solh_translation_coeffs

  !> \warning Requires precalculated values of the real spherical harmonics at the position of the CGTO nucleus. The coupling coefficients should also be precalculated for performance reasons.
  subroutine calculate_lambda_couplings(CGTO_L,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
     implicit none
     real(kind=cfp), allocatable :: transl_cfs(:,:,:), Xlm_CGTO_center(:)
     integer, intent(in) :: CGTO_L, l, m
     !OUTPUT:
     real(kind=cfp), allocatable :: c_lambda(:,:,:)

     real(kind=cfp) :: cf, sum_mu
     integer :: CGTO_M, CGTO_M_ind, lm, lp, lpp, mp, mu, lambda, ind, lambda_max, err, d1, d2, d3

        lambda_max = l+CGTO_L

        d1 = max(lambda_max,1)
        d2 = max(CGTO_L,1)
        d3 = 2*CGTO_L+1
        if (allocated(c_lambda)) then
           if (ubound(c_lambda,1) < d1 .or. ubound(c_lambda,2) < d2 .or. ubound(c_lambda,3) < d3) then
              deallocate (c_lambda)
           end if
        end if
        if (.not. allocated(c_lambda)) then
           allocate(c_lambda(0:d1,0:d2,d3),stat=err)
           if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','calculate_lambda_couplings','Memory allocation failed.',err,1)
        endif

        do CGTO_M=-CGTO_L,CGTO_L
           CGTO_M_ind = CGTO_L+CGTO_M+1
           lm = l*l+l+m+1
           do lp=0,CGTO_L
              lpp = CGTO_L-lp

              c_lambda(0:lambda_max,lp,CGTO_M_ind) = 0.0_cfp
              do lambda=abs(l-lp),l+lp
                 if (mod(l+lp+lambda,2) .ne. 0) cycle !selection rule for the Gaunt coefficients
                 do mp=-lp,lp
                    cf = transl_cfs(mp+lp+1,lp,CGTO_M_ind)
                    if (cf .eq. 0.0_cfp) cycle
                    sum_mu = 0.0_cfp
                    do mu=-lambda,lambda
                       ind = lambda*lambda+lambda+mu+1
                       sum_mu = sum_mu + cf*cpl%rgaunt(l,lp,lambda,m,mp,mu)*Xlm_CGTO_center(ind)
                    enddo !mu
                    c_lambda(lambda,lp,CGTO_M_ind) = c_lambda(lambda,lp,CGTO_M_ind) + sum_mu
                 enddo !mp
              enddo !lambda

           enddo !lp
        enddo !CGTO_M
     
  end subroutine calculate_lambda_couplings

  function sum_besi_half_asym_lambda_cf(z,couplings,l,lmin,lmax)
     use phys_const_gbl, only: twopi
     implicit none
     integer, intent(in) :: l, lmin,lmax
     real(kind=cfp), intent(in) :: z, couplings(0:l)
     real(kind=cfp) :: sum_besi_half_asym_lambda_cf

     integer :: k, klim, lambda
     real(kind=cfp) :: nu, z8, s, t(lmin:lmax+1)

        if (cfp .eq. wp) then
           klim = 25
        else
           klim = 30
        endif
   
        z8 = 8.0_cfp*z
        sum_besi_half_asym_lambda_cf = sum(couplings(lmin:lmax))
        t(:) = 1.0_cfp
        do k=1,klim
           s = 0.0_cfp
           do lambda=lmin,lmax
              nu = 4*(lambda+0.5_cfp)**2
              t(lambda) = t(lambda)*(nu-(2*k-1)**2)/(k*z8)
              s = s + (-1)**k*t(lambda)*couplings(lambda)
           enddo
           sum_besi_half_asym_lambda_cf = sum_besi_half_asym_lambda_cf + s
        enddo !k
        sum_besi_half_asym_lambda_cf = sum_besi_half_asym_lambda_cf/sqrt(twopi*z)

  end function sum_besi_half_asym_lambda_cf

  function besi_half_asym(z,l)
     use phys_const_gbl, only: twopi
     implicit none
     real(kind=cfp), intent(in) :: z
     integer, intent(in) :: l
     real(kind=cfp) :: besi_half_asym
 
     integer :: k, klim
     real(kind=cfp) :: nu, z8, t

        if (cfp .eq. wp) then 
           klim = 25
        else
           klim = 30
        endif
   
        z8 = 8.0_cfp*z
        nu = 4*(l+0.5_cfp)**2
        besi_half_asym = 1.0_cfp
        t = 1.0_cfp
        do k=1,klim
           t = t*(nu-(2*k-1)**2)/(k*z8)
           besi_half_asym = besi_half_asym + (-1)**k*t
        enddo !k
        besi_half_asym = besi_half_asym/sqrt(twopi*z)

  end function besi_half_asym

  !> Calculates the real spherical harmonics for all l,m up to l=max_l for directions corresponding to the nuclear positions.
  subroutine precalculate_Xlm_for_nuclei(nuclei,max_l,Xlm_nuclei,n_Xlm_nuclei)
     use phys_const_gbl, only: fourpi
     use common_obj_gbl, only: nucleus_type
     use const_gbl, only: epsabs
     implicit none
     type(nucleus_type), intent(in) :: nuclei(:)
     integer, intent(in) :: max_l
     real(kind=cfp), allocatable :: Xlm_nuclei(:)
     integer, intent(out) :: n_Xlm_nuclei

     real(kind=cfp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     integer :: n_nuc, ind, err, l, m, base
     real(kind=cfp) :: RH(-max_l:max_l,0:max_l+1), inv_distance, R(3), distance

        if (allocated(Xlm_nuclei)) deallocate(Xlm_nuclei)

        n_nuc = size(nuclei)
        n_Xlm_nuclei = (max_l+1)**2

        allocate(Xlm_nuclei(n_nuc*n_Xlm_nuclei),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','precalculate_Xlm_for_nuclei','Memory allocation 1 failed.',err,1)

        Xlm_nuclei = 0.0_cfp
        RH = 0.0_cfp
        do ind=1,n_nuc
           distance = sqrt(dot_product(nuclei(ind)%center,nuclei(ind)%center))
           if (max_l > 0 .and. distance > epsabs) then
              inv_distance = 1.0_cfp/distance
              R(1:3) = nuclei(ind)%center(1:3)*inv_distance !normalize the position vectors to 1
              call resh(RH,R(1),R(2),R(3),max_l)
           else
              RH = 0.0_cfp
              RH(0,0) = norm
           endif

           base = (ind-1)*n_Xlm_nuclei
           do l=0,max_l
              do m=-l,l
                 Xlm_nuclei(base+l*l+l+m+1) = RH(m,l) !save Xlm in the order: m,l,nucleus
              enddo !m
           enddo !l
        enddo !ind

  end subroutine precalculate_Xlm_for_nuclei

  !> Calculates the real spherical harmonics for all l,m up to l=max_l for directions corresponding to the nuclear positions.
  subroutine precalculate_Xlm_for_CGTO_center(RA,max_l,Xlm_CGTO_center)
     use phys_const_gbl, only: fourpi
     implicit none
     real(kind=cfp), intent(in) :: RA(3)
     integer, intent(in) :: max_l
     !OUTPUT:
     real(kind=cfp), allocatable :: Xlm_CGTO_center(:)

     real(kind=cfp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     integer :: err, l, m, n_Xlm
     real(kind=cfp) :: RH(-max_l:max_l,0:max_l+1), R(3), d

        if (allocated(Xlm_CGTO_center)) deallocate(Xlm_CGTO_center)

        n_Xlm = (max_l+1)**2

        allocate(Xlm_CGTO_center(n_Xlm),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','precalculate_Xlm_for_CGTO_center','Memory allocation 1 failed.',err,1)

        !Normalize the direction vector
        d = sqrt(dot_product(RA,RA))
        if (d .ne. 0.0_cfp) then
           R(1:3) = RA(1:3)/sqrt(dot_product(RA,RA))
        else
           Xlm_CGTO_center(1) = norm
           Xlm_CGTO_center(2:(max_l+1)**2) = 0.0_cfp
           return
        endif

        RH = 0.0_cfp
        if (max_l > 0) then
           call resh(RH,R(1),R(2),R(3),max_l)
        else
           RH(0,0) = norm
        endif

        do l=0,max_l
           do m=-l,l
              Xlm_CGTO_center(l*l+l+m+1) = RH(m,l) !save Xlm in the order: m,l
           enddo !m
        enddo !l

  end subroutine precalculate_Xlm_for_CGTO_center

  subroutine precalculate_Xlm_for_CGTO_product_center(n_contr_pairs,RA,max_l,Xlm_product_CGTO_center)
     use phys_const_gbl, only: fourpi
     implicit none
     integer, intent(in) :: max_l, n_contr_pairs
     real(kind=cfp), intent(in) :: RA(3,n_contr_pairs)
     !OUTPUT:
     real(kind=cfp), allocatable :: Xlm_product_CGTO_center(:)

     real(kind=cfp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     integer :: err, l, m, n_Xlm, i, base, ind
     real(kind=cfp) :: RH(-max_l:max_l,0:max_l+1), R(3), d

        if (allocated(Xlm_product_CGTO_center)) deallocate(Xlm_product_CGTO_center)

        n_Xlm = (max_l+1)**2

        allocate(Xlm_product_CGTO_center(n_Xlm*n_contr_pairs),stat=err)
        if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','precalculate_Xlm_for_CGTO_center','Memory allocation 1 failed.',err,1)

        do i=1,n_contr_pairs
           !Normalize the direction vector
           d = sqrt(dot_product(RA(1:3,i),RA(1:3,i)))
           if (d .ne. 0.0_cfp) then
              R(1:3) = RA(1:3,i)/d
           else
              R = 0.0_cfp
           endif

           RH = 0.0_cfp
           if (max_l > 0) then
              call resh(RH,R(1),R(2),R(3),max_l)
           else
              RH(0,0) = norm
           endif

           do l=0,max_l
              do m=-l,l
                 ind = l*l+l+m+1
                 base = n_contr_pairs*(ind-1)
                 Xlm_product_CGTO_center(i+base) = RH(m,l) !save Xlm in the order: (i),(m,l) i.e. ~2D array
              enddo !m
           enddo !l
        enddo !i

  end subroutine precalculate_Xlm_for_CGTO_product_center

  subroutine calc_resh_coefficients(L)
    implicit none
    integer, intent(in) :: L

    integer :: l_it, m_it, ind

       !$OMP CRITICAL
       if (L > max_l) then
          max_l = L
          if (allocated(a)) deallocate(a,b,c)
          allocate(a((L+1)**2),b((L+1)**2),c(L+1))
          ind = 0
          do l_it=1,L-1
             c(l_it) = sqrt((2.0_cfp*l_it+3.0_cfp)/(2.0_cfp*l_it+2.0_cfp))
             do m_it=-l_it,l_it
                a(ind+m_it+l_it+1) = sqrt((4*l_it*l_it+8*l_it+3.0_cfp)/((l_it+1)**2-m_it*m_it))
                b(ind+m_it+l_it+1) = sqrt((l_it*l_it-m_it*m_it)/(4*l_it*l_it-1.0_cfp))
             enddo
             ind = ind + 2*l_it+1
          enddo
       endif
       !$OMP END CRITICAL

  end subroutine calc_resh_coefficients

  !> Calculates the real spherical harmonics assuming that the input values X,Y,Z lie on the unit sphere.
  !> \warning This routine is not threadsafe unless calc_resh_coefficients was called before the parallel region for high enough L to ensure the coefficients are never recalculated while using this routine.
  subroutine resh(SH,X,Y,Z,L)
    use precisn_gbl
    use phys_const_gbl, only: fourpi
    implicit none
    integer, intent(in) :: L
    real(kind=cfp), intent(out) :: SH(-L:L,0:L)
    real(kind=cfp), intent(in) :: x, y, z
  
    integer :: l_it, m_it, ind
    real(kind=cfp), parameter :: norm1 = sqrt(3.0_cfp/fourpi), norm0 = 1.0_cfp/sqrt(fourpi)
  
       SH = 0.0_cfp !vectorized
 
       if (L .eq. 0) then
          SH(0,0) = norm0
          return
       endif
     
       !initialize the starting values
       SH(0,0) = norm0
       SH(-1,0) = 0.0_cfp
     
       !recursion
       !l_it = 0 case: L=1
       SH(1,1) = norm1*x
       SH(0,1) = norm1*z
       SH(-1,1)= norm1*y
   
       !check that we have the precalculated coefficients a,b,c
       if (L > max_l) then
          call xermsg ('cgto_pw_expansions_mod', 'resh', &
                       'Requesting resh for L larger than for which it has been initialized', L, 1)
       endif
   
       !diagonal recursions
       ind = 0
       do l_it = 1, L-1
          SH(l_it+1, l_it+1) = c(l_it)*(x*SH(l_it,l_it)-y*SH(-l_it,l_it))
          SH(-l_it-1,l_it+1) = c(l_it)*(y*SH(l_it,l_it)+x*SH(-l_it,l_it))
          !vertical recursions (vectorized)
          !$OMP SIMD
          do m_it = -l_it,l_it
             SH(m_it,l_it+1) = a(ind+m_it+l_it+1)*(z*SH(m_it,l_it)-b(ind+m_it+l_it+1)*SH(m_it,l_it-1))
          end do
          ind = ind + 2*l_it+1
       end do

  end subroutine resh

  subroutine calc_solh_coefficients(L)
    implicit none
    integer, intent(in) :: L
    integer :: l_it, m_it, l2p1, ind
    real(kind=cfp) :: lp1

       !$OMP CRITICAL
       if (L > max_ls) then
          max_ls = L
          if (allocated(as)) deallocate(as,bs,cs)
          allocate(as((L+1)**2),bs((L+1)**2),cs((L+1)**2))
          ind = 0
          do l_it = 1, L-1
             cs(l_it) = sqrt((2.0_cfp*l_it+1.0_cfp)/(2.0_cfp*l_it+2.0_cfp))
             l2p1 = 2*l_it + 1
             lp1 = l_it + 1.0_cfp
             do m_it = -l_it,l_it
                as(ind+m_it+l_it+1) = l2p1/sqrt((lp1-m_it)*(lp1+m_it))
                bs(ind+m_it+l_it+1) = sqrt((l_it+m_it)*(l_it-m_it)/((lp1-m_it)*(lp1+m_it)))
             end do
             ind = ind + 2*l_it+1
          enddo
        endif
       !$OMP END CRITICAL

  end subroutine calc_solh_coefficients

  !> Calculates the real solid harmonic.
  !> \warning This routine is not threadsafe unless calc_solh_coefficients was called before the parallel region for high enough L to ensure the coefficients are never recalculated while using this routine.
  subroutine solh(SH,x,y,z,L)
    use precisn_gbl
    implicit none
    integer, intent(in) :: L
    real(kind=cfp), intent(out) :: SH(-L:L,0:L+1)
    real(kind=cfp), intent(in) :: x, y, z
    
    integer :: l_it, m_it, ind
    real(kind=cfp) :: rsq
    
       SH = 0.0_cfp !vectorized
   
       !initialize the starting values
       SH(0,0) = 1.0_cfp
       SH(-1,0) = 0.0_cfp
       rsq = x*x + y*y + z*z 
   
       !recursion
       !l_it = 0 case: L=1
       SH(1,1) = x
       SH(0,1) = z
       SH(-1,1)= y
   
       !check that we have the precalculated coefficients a,b,c
       if (L > max_ls) then
          call xermsg ('cgto_pw_expansions_mod', 'solh', &
                       'Requesting solh for L larger than for which it has been initialized', L, 1)
       endif
   
       !diagonal recursions
       ind = 0
       do l_it = 1, L-1
          SH(l_it+1, l_it+1) = cs(l_it)*(x*SH(l_it,l_it)-y*SH(-l_it,l_it))
          SH(-l_it-1,l_it+1) = cs(l_it)*(y*SH(l_it,l_it)+x*SH(-l_it,l_it))
          !vertical recursions (vectorized)
          !$OMP SIMD
          do m_it = -l_it,l_it
             SH(m_it,l_it+1) = z*SH(m_it,l_it)*as(ind+m_it+l_it+1)-bs(ind+m_it+l_it+1)*rsq*SH(m_it,l_it-1)
          end do
          ind = ind + 2*l_it+1
       end do
  
  end subroutine solh

end module cgto_pw_expansions_gbl
