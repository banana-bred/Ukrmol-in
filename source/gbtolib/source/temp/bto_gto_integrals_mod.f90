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
module bto_gto_integrals_mod
  use lebedev, only: available_table, order_table, ld_by_order, rule_max, mmax
  use general_quadrature
  use gto_data_obj
  use bto_data_obj
  use common_obj, only: nucleus_type
  use coupling_obj
  implicit none

  private

  public BTO_GTO_one_electron_integrals, calculate_Tl_matrix, solve_poisson_equation
  public BTO_GTO_two_el_BB_GG, BTO_GTO_two_el_BG_BG, BTO_GTO_two_el_BG_GG

  !> These routines must be made accessible for the module: bto_gto_test_routines_mod
  public calc_resh_coefficients, calc_solh_coefficients, radial_grid_cgto_pair, omp_calculate_cgto_pair_pw_coefficients_analytic, omp_calculate_CGTO_pw_coefficients_analytic

  !> These are used to keep the largest L for which the coefficients a,b,c and as,bs,cs for the real spherical and real solid harmonic have been precalculated.
  integer :: max_l = -1, max_ls = -1
  !> Coefficients used for calculation of real spherical and real solid harmonics using the routines resh, solh.
  real(kind=cfp), allocatable :: a(:), b(:), c(:), as(:), bs(:), cs(:)

  !> Number of points in each of the Lebedev quadrature rules (n_points) and the offset in Xlm_Lebedev for each quadrature rule.
  integer, allocatable :: n_points(:), offset(:)

  !> Maximum L of the Xlm real spherical harmonics for which the values at the Lebedev quadrature points have been precalculated.
  integer :: max_l_precalculated = -1

  !> Values of the real spherical harmonics evaluated at all Lebedev quadrature points (Xlm_Lebedev) and values of the real spherical harmonics at the positions of the nuclei (Xlm_nuclei) and
  !> values of the real spherical harmonics evaluated at the center of the CGTO.
  real(kind=cfp), allocatable :: Xlm_Lebedev(:), Xlm_nuclei(:)
  !> Stride in array Xlm_nuclei.
  integer :: n_Xlm_nuclei = 0

  !> Used to evaluate various coupling coefficients.
  type(couplings_type) :: cpl

  !> Object that is used to evaluate the product of a CGTO and a real spherical harmonic at arbitrary points in space. This is used when integrating this function adaptively over sphere.
  type, extends(function_2d) :: Xlm_x_cgto_surface
     !> Contains the data on the CGTO.
     type(cgto_data) :: cgto
     !> radial distance for which to evaluate.
     real(kind=cfp) :: r = 0.0_cfp
     !> Angular numbers of the real spherical harmonic centered on CSM multiplying the GTO.
     integer :: l = 0, m = 0
  contains
     procedure :: eval => eval_Xlm_x_cgto_surface
  end type Xlm_x_cgto_surface

  !> Object that is used to evaluate the product of a pair of CGTOs and a real spherical harmonic at arbitrary points in space. This is used when integrating this function adaptively over sphere.
  type, extends(function_2d) :: Xlm_x_pair_cgto_surface
     !> Contains the data on the CGTO.
     type(cgto_data) :: cgto_A, cgto_B
     !> radial distance for which to evaluate.
     real(kind=cfp) :: r = 0.0_cfp
     !> Angular numbers of the real spherical harmonic centered on CSM multiplying the GTO.
     integer :: l = 0, m = 0
  contains
     procedure :: eval => eval_Xlm_x_pair_cgto_surface
  end type Xlm_x_pair_cgto_surface

  !> Used to evaluate matrix elements of the radial free Schroedinger equation in the radial B-spline basis.
  type, extends(bound_user_function) :: radial_free_bspline
     !> Contains the radial grid for evaluation of the B-splines.
     type(bto_data) :: bspline_radial_grid
  contains
     procedure :: wp_eval => wp_eval_radial_free_bspline
     procedure :: ep_eval => ep_eval_radial_free_bspline
  end type radial_free_bspline

  public wp_eval_radial_free_bspline, ep_eval_radial_free_bspline, eval_Xlm_x_cgto_surface

contains

  !> Calculates the real spherical harmonics for all (l,m) combinations up to l=L on all available Lebedev grids. The values are stored in the module private array Xlm_Lebedev.
  subroutine precalculate_Xlm_on_Lebedev_grid(L)
     use const, only: Mib
     use phys_const, only: fourpi
     implicit none
     integer, intent(in) :: L

     integer :: i, j, available, n, err, m, point
     real(kind=wp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     real(kind=cfp) :: r1(mmax), r2(mmax), r3(mmax), w(mmax)
     real(kind=cfp) :: RH(-L:L,0:L+1)

        if (L .le. max_l_precalculated) return

        if (allocated(Xlm_Lebedev)) deallocate(Xlm_Lebedev)
        if (allocated(n_points)) deallocate(n_points)
        if (allocated(offset)) deallocate(offset)

        allocate(offset(rule_max),n_points(rule_max),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','precalculate_Xlm_on_Lebedev_grid','Memory allocation 1 failed.',err,1)

        offset = 0
        n_points = 0

        !Loop over all quadratures and count the total number of points for which Xlm will be required
        i = 0
        do n = 1,rule_max
           available = available_table(n)

           if (available == 1) then
              n_points(n) = order_table(n)
              offset(n) = i*(L+1)**2 !total number of preceeding points
              i = i + n_points(n)
           endif

        enddo !n

        i = i*(L+1)**2
        allocate(Xlm_Lebedev(i),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','precalculate_Xlm_on_Lebedev_grid','Memory allocation 2 failed.',err,1)

        !Precalculate the Xlm at the quadrature points and store them in the order (m,l,point,n).
        RH = 0.0_cfp
        do n = 1,rule_max
           available = available_table(n)

           if (available == 1) then
              call ld_by_order(n_points(n),r1,r2,r3,w)

              j = offset(n)
              do point=1,n_points(n)
                 if (L > 0) then
                    call resh(RH,r1(point),r2(point),r3(point),L)
                 else
                    RH(0,0) = norm
                 endif
   
                 do i=0,L
                    do m=-i,i
                       j = j + 1
                       Xlm_Lebedev(j) = RH(m,i)
                    enddo !m
                 enddo !i
              enddo !point

           endif
        enddo !n

        max_l_precalculated = L

  end subroutine precalculate_Xlm_on_Lebedev_grid

  !> Calculates the real spherical harmonics for all l,m up to l=max_l for directions corresponding to the nuclear positions.
  subroutine precalculate_Xlm_for_nuclei(nuclei,max_l)
     use phys_const, only: fourpi
     use common_obj, only: nucleus_type
     implicit none
     type(nucleus_type), intent(in) :: nuclei(:)
     integer, intent(in) :: max_l

     real(kind=wp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     integer :: n_nuc, ind, err, l, m, base
     real(kind=cfp) :: RH(-max_l:max_l,0:max_l+1), inv_distance, R(3)

        if (allocated(Xlm_nuclei)) deallocate(Xlm_nuclei)

        n_nuc = size(nuclei)
        n_Xlm_nuclei = (max_l+1)**2

        allocate(Xlm_nuclei(n_nuc*n_Xlm_nuclei),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','precalculate_Xlm_for_nuclei','Memory allocation 1 failed.',err,1)

        RH = 0.0_cfp
        do ind=1,n_nuc
           inv_distance = 1.0_cfp/sqrt(dot_product(nuclei(ind)%center,nuclei(ind)%center))
           R(1:3) = nuclei(ind)%center(1:3)*inv_distance !normalize the position vectors to 1
           if (max_l > 0) then
              call resh(RH,R(1),R(2),R(3),max_l)
           else
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
     use phys_const, only: fourpi
     implicit none
     real(kind=cfp), intent(in) :: RA(3)
     integer, intent(in) :: max_l
     !OUTPUT:
     real(kind=cfp), allocatable :: Xlm_CGTO_center(:)

     real(kind=wp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     integer :: err, l, m, n_Xlm
     real(kind=cfp) :: RH(-max_l:max_l,0:max_l+1), R(3)

        if (allocated(Xlm_CGTO_center)) deallocate(Xlm_CGTO_center)

        n_Xlm = (max_l+1)**2

        allocate(Xlm_CGTO_center(n_Xlm),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','precalculate_Xlm_for_CGTO_center','Memory allocation 1 failed.',err,1)

        !Normalize the direction vector
        R(1:3) = RA(1:3)/sqrt(dot_product(RA,RA))

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
     use phys_const, only: fourpi
     implicit none
     real(kind=cfp), intent(in) :: RA(3,n_contr_pairs)
     integer, intent(in) :: max_l, n_contr_pairs
     !OUTPUT:
     real(kind=cfp), allocatable :: Xlm_product_CGTO_center(:)

     real(kind=wp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
     integer :: err, l, m, n_Xlm, i, base, ind
     real(kind=cfp) :: RH(-max_l:max_l,0:max_l+1), R(3)

        if (allocated(Xlm_product_CGTO_center)) deallocate(Xlm_product_CGTO_center)

        n_Xlm = (max_l+1)**2

        allocate(Xlm_product_CGTO_center(n_Xlm*n_contr_pairs),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','precalculate_Xlm_for_CGTO_center','Memory allocation 1 failed.',err,1)

        do i=1,n_contr_pairs
           !Normalize the direction vector
           R(1:3) = RA(1:3,i)/sqrt(dot_product(RA(1:3,i),RA(1:3,i)))

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

  !> Calculates the one electron mixed integrals between the given shell of CGTOs and all BTOs. It is required that the first break point of the BTO grid lies radially behind the center of the CGTO.
  !> The nuclei that form the nuclear potential is also required to lie radially before the BTOs. The calculated integrals are output in the array integrals in the order: 
  !> (BTO m,BTO l,BTO ind),(CGTO m),(integral type). The types of integrals calculated here are the following: overlap=1,kinetic energy=2,nuclear attraction=3,property (up to l=2)=4,...,12.
  !> \warning The integrals are calculated for the whole basis of radial B-splines but the Bloch term is only added for r=cms_bto%B.
  !> \todo make sure the integrals are calculated only for rad. b-splines starting from cms_bto%ind_0_der!!!!
  !> \todo this routine can be called only if the CGTO is indeed non-negligible on the BTO grid otherwise the routine radial_grid_CGTO will not work
  subroutine BTO_GTO_one_electron_integrals(cms_bto,cgto,nuclei,integrals)
     use const, only: epsrel, epsabs
     use phys_const, only: fourpi
     use general_quadrature, only: n_10,w_10,x_10
     use bspline_base, only: bvalu
     use gto_function_obj
     use common_obj, only: nucleus_type
     use omp_lib
     implicit none
     type(bto_data), intent(inout) :: cms_bto
     type(cgto_data), intent(inout) :: cgto
     real(kind=cfp), allocatable :: integrals(:,:,:)
     type(nucleus_type), intent(in) :: nuclei(:)

     integer, parameter :: n_rng_knot = 1 !Number of intervals in between each pair of knots within which the quadrature will be applied
     integer, parameter :: max_prop_l = 2
     !Maximum L to use in the Legendre expansion of the Coulomb potential.
     integer, parameter :: max_l_legendre = 30

     real(kind=cfp) :: A, B, R, RA, threshold, olap, kei, prop, fac, R_true
     real(kind=cfp), allocatable :: r_points(:), weights(:), bto_norms(:), nuclei_RA(:)
     integer :: n_total_points, err, max_l_pw, n_Xlm, n_cgto_m, i, ind, m_ind, lm, base, l, m, n_threads, rank, iam, n_properties, n_Xlm_bto, lp_mp, lp, mp
     integer :: n_nuc, max_l_aux
     real(kind=cfp), allocatable :: angular_integrals(:,:,:) !(radial points,CGTO m,lm)
     real(kind=cfp), allocatable :: gaunt_prop_integrals(:,:,:,:), bto_val_0(:,:), bto_val_2(:,:)
     real(kind=cfp), allocatable :: gaunt_angular_integrals(:,:,:,:), lp_integrals(:,:,:,:), bto_radial(:,:)

        call cms_bto%init_grid
        call cms_bto%normalize

        call cgto%normalize

        RA = sqrt(dot_product(cgto%RA,cgto%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','The CGTO center must lie radially before the start of the BTO grid.',1,1)

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs
!
!------ For the present CGTO shrink the radial grid defined for the whole B-spline range [cms_bto%A,cms_bto%B] into the range which corresponds to the extent of the CGTO [cms_bto%A,R].
!       The resulting quadrature points and weights are in the arrays r_points,weights.
        call radial_grid_CGTO(cgto,cms_bto,threshold,x_10,w_10,n_10,n_rng_knot,r_points,weights,n_total_points,R,R_true)
!
        n_properties = (max_prop_l+1)**2

        !cms_bto%l is assumed to be set to the largest angular momentum in the whole BTO basis.
        max_l_pw = cms_bto%l+max(max_prop_l,max_l_legendre)
        n_Xlm = (max_l_pw+1)**2 !total number of partial wave projections, i.e. Xlm
        n_Xlm_bto = (cms_bto%l+1)**2 !total number of BTO (l,m) combinations
        n_cgto_m = 2*cgto%l+1
        n_nuc = size(nuclei)

        max_l_aux = max_l_pw+cgto%l

        if (n_nuc .le. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','On input the array of nuclei has not been allocated.',2,1)

        if (allocated(integrals)) deallocate(integrals)

        allocate(integrals(n_Xlm_bto*cms_bto%n,n_cgto_m,3+n_properties),bto_norms(cms_bto%n),nuclei_RA(n_nuc),bto_val_0(cms_bto%n,n_total_points),bto_val_2(cms_bto%n,n_total_points),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','Memory allocation failed.',err,1)

        integrals = 0.0_cfp

        do i=1,n_nuc
           nuclei_RA(i) = sqrt(dot_product(nuclei(i)%center,nuclei(i)%center))
           if (nuclei_RA(i) .ge. cms_bto%A) &
          &call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','All nuclei must lie radially before the start of the BTO grid.',3,1)
        enddo
!
!------ Precalculate all coupling coefficients needed for the calculations.

        !The two routines below must be called before the OpenMP region to ensure thread safety of resh, solh.
        call calc_resh_coefficients(max_l_aux)
        call calc_solh_coefficients(max_l_aux)

        !Precalculate all couplings needed to evaluate the overlap-type integrals and the nuclear attraction integrals.
        print *,'gaunt'
        call cpl%prec_cgaunt(max_l_aux)

        !Real spherical harmonics at the Lebedev quadrature points: result in the module array Xlm_Lebedev.
        !print *,'leb'
        !call precalculate_Xlm_on_Lebedev_grid(max_l_pw)
!
!------ Calculate the projections of the CGTOs on the real spherical harmonics (i.e. PW-expansion) for all radial grid points.
        print *,'pw'
        call omp_calculate_CGTO_pw_coefficients_analytic(max_l_pw,cgto,r_points,angular_integrals)
        print *,'done'
!
!------ Calculate and store norms of all BTOs
        cms_bto%bcoef = 0.0_cfp
        do ind=1,cms_bto%n
           cms_bto%bcoef(ind) = 1.0_cfp
           cms_bto%ind = ind
           cms_bto%r1 = cms_bto%knots(ind)
           cms_bto%r2 = cms_bto%knots(ind+cms_bto%order)
           call cms_bto%normalize
           cms_bto%bcoef(ind) = 0.0_cfp !clean-up for the next
           bto_norms(ind) = cms_bto%norm
        enddo
!
!------ Evaluate the B-splines at the quadrature points. This is done in serial: it seems like bvalu is not threadsafe for some reason!
        cms_bto%bcoef = 0.0_cfp
        do i=1,n_total_points
           do ind=1,cms_bto%n
              if (r_is_inside_A_B(cms_bto%knots(ind),cms_bto%knots(ind+cms_bto%order),r_points(i))) then !this B-spline is non-zero at this quadrature point
                 !Evaluate the B-spline part times Jac of the radial integrand; todo precalculate the B-spline values so I can parallelize the loop over the radial points
                 cms_bto%bcoef(ind) = 1.0_cfp
                 bto_val_0(ind,i) = bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,0,r_points(i),cms_bto%inbv,cms_bto%work)
                 bto_val_2(ind,i) = bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,2,r_points(i),cms_bto%inbv,cms_bto%work)
                 cms_bto%bcoef(ind) = 0.0_cfp !clean-up for the next
              endif
           enddo !ind
        enddo !i
!
!------ Multiply-in and sum over the coupling coefficients: needed for the property integrals
        call omp_mult_gaunt_into_angular_integrals(angular_integrals,cms_bto%l,max_prop_l,gaunt_prop_integrals)
!
!------ Calculate the radial integrals for all BTOs: at each quadrature point calculate its contribution to integrals over all BTOs non-zero at that point.
        allocate(bto_radial(n_total_points,cms_bto%n),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','Memory allocation 2 failed.',err,1)
        !$OMP PARALLEL DEFAULT(NONE) SHARED(n_total_points,n_Xlm,n_Xlm_bto,n_cgto_m,r_points,weights,angular_integrals,cms_bto,integrals,bto_val_0,bto_val_2,bto_norms,gaunt_prop_integrals,bto_radial,cgto) &
        !$OMP & PRIVATE(i,iam,n_threads,ind,olap,m_ind,base,lm,l,kei,m,lp,fac,mp,prop,lp_mp)
        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()
!
!------ Overlap integrals
        !$OMP SINGLE
        bto_radial = 0.0_cfp
        !$OMP END SINGLE
        !$OMP DO
        do i=1,n_total_points
          do ind=1,cms_bto%n
             if (r_is_inside_A_B(cms_bto%knots(ind),cms_bto%knots(ind+cms_bto%order),r_points(i))) then !this B-spline is non-zero at this quadrature point
                bto_radial(i,ind) = weights(i)*bto_val_0(ind,i)*bto_norms(ind)*r_points(i) !B-spline part times Jac of the radial integrand times the quadrature weight
             endif
          enddo !ind
        enddo !i
        !$OMP END DO

        !For all m-values of the CGTO function
        do m_ind=1,n_cgto_m
           !For all BTO angular parts
           do lm=1,n_Xlm_bto
              if (mod(lm,n_threads) .ne. iam) cycle !work distribution
              !For all radial B-splines
              do ind=1,cms_bto%n
                 olap = 0.0_cfp
                 do i=1,n_total_points
                    olap = olap + bto_radial(i,ind)*angular_integrals(i,m_ind,lm)
                 enddo
                 !The integrals are stored in the order: (BTO m,BTO l,BTO radial ind),(CGTO m),(olap=1)
                 base = n_Xlm_bto*(ind-1)
                 integrals(lm+base,m_ind,1) = olap
              enddo !ind
           enddo !lm
        enddo !m_ind
        !$OMP BARRIER
!
!------ Kinetic energy integrals
        !For all m-values of the CGTO function
        do m_ind=1,n_cgto_m
           !For all BTO angular parts
           do l=0,cms_bto%l
               !$OMP SINGLE
               bto_radial = 0.0_cfp
               !$OMP END SINGLE
               !$OMP DO
               do i=1,n_total_points
                 do ind=1,cms_bto%n
                    if (r_is_inside_A_B(cms_bto%knots(ind),cms_bto%knots(ind+cms_bto%order),r_points(i))) then !this B-spline is non-zero at this quadrature point
                       bto_radial(i,ind) = weights(i)*bto_norms(ind)*(bto_val_2(ind,i)*r_points(i)-l*(l+1)/r_points(i)*bto_val_0(ind,i)) !precision loss is possible for small values of r
                    endif
                 enddo !ind
              enddo !i
              !$OMP END DO
              do m=-l,l
                 lm = l*l+l+m+1
                 if (mod(lm,n_threads) .ne. iam) cycle !work distribution
                 !For all radial B-splines
                 do ind=1,cms_bto%n
                    kei = 0.0_cfp
                    do i=1,n_total_points
                       kei = kei + bto_radial(i,ind)*angular_integrals(i,m_ind,lm)
                    enddo
                    !The integrals are stored in the order: (BTO m,BTO l,BTO radial ind),(CGTO m),(kei=2)
                    base = n_Xlm_bto*(ind-1)
                    integrals(lm+base,m_ind,2) = -0.5_cfp*kei
                 enddo !ind
              enddo !m
              !$OMP BARRIER
           enddo !l
        enddo !m_ind
        !Add the Bloch terms at the end if needed (CGTO is a continuum basis function centered at CMS).
        !$OMP SINGLE
        if (cgto%nuc .eq. 0) call add_cgto_bto_bloch_elements(bto_norms,cms_bto,cgto,integrals(1:n_Xlm_bto*cms_bto%n,1:n_cgto_m,2))
        !$OMP END SINGLE
!
!------ Property integrals: (lp,mp) = property l,m values
        !For all property (l,m)
        do lp=0,max_prop_l
           do mp=-lp,lp
              lp_mp = lp*lp + lp+mp+1

              fac = (-1)**mp*sqrt(fourpi/(2*lp+1.0_cfp))
              !$OMP SINGLE
              bto_radial = 0.0_cfp
              !$OMP END SINGLE
              !$OMP DO
              do i=1,n_total_points
                 prop = fac*r_points(i)**lp
                 do ind=1,cms_bto%n
                    if (r_is_inside_A_B(cms_bto%knots(ind),cms_bto%knots(ind+cms_bto%order),r_points(i))) then !this B-spline is non-zero at this quadrature point
                       bto_radial(i,ind) = weights(i)*bto_val_0(ind,i)*bto_norms(ind)*r_points(i)*prop
                    endif
                 enddo !ind
              enddo !i
              !$OMP END DO
              !For all BTO angular parts
              do l=0,cms_bto%l
                 do m=-l,l
                    lm = l*l+l+m+1
                    if (mod(lm,n_threads) .ne. iam) cycle !work distribution
                    !For all m-values of the CGTO function
                    do m_ind=1,n_cgto_m
                       !For all radial B-splines
                       do ind=1,cms_bto%n
                          prop = 0.0_cfp
                          do i=1,n_total_points
                             prop = prop + bto_radial(i,ind)*gaunt_prop_integrals(i,m_ind,lp_mp,lm)
                          enddo
                          !The integrals are stored in the order: (BTO m,BTO l,BTO radial ind),(CGTO m),(kei=2)
                          base = n_Xlm_bto*(ind-1)
                          integrals(lm+base,m_ind,3+lp_mp) = prop
                       enddo !ind
                    enddo !m_ind
                 enddo !m
              enddo !l
              !$OMP BARRIER
           enddo !mp
        enddo !lp
        !$OMP END PARALLEL
!
!------ Nuclear attraction integrals
        print *,'step 1'
        call omp_mult_gaunt_into_angular_integrals(angular_integrals,cms_bto%l,max_l_legendre,gaunt_angular_integrals)
        deallocate(angular_integrals)
        print *,'step 2'
        call sum_mp_legendre(gaunt_angular_integrals,r_points,weights,cgto%l,max_l_legendre,cms_bto%l,nuclei,nuclei_RA,lp_integrals)
        deallocate(gaunt_angular_integrals)
        print *,'step 3'
        call integrate_legendre(lp_integrals,bto_norms,bto_val_0,cms_bto%knots,cms_bto%order,cms_bto%n,cgto%l,max_l_legendre,cms_bto%l,integrals(1:n_Xlm_bto*cms_bto%n,1:n_cgto_m,3))

        !get rid of all arrays except 'integrals' which is the final output
        deallocate(r_points,weights,gaunt_prop_integrals,bto_norms,bto_val_0,bto_val_2)

  end subroutine BTO_GTO_one_electron_integrals
 
  !> Determines the radial grid needed to describe, within the accuracy threshold, integrals involving the CGTO and the B-splines on the given grid. We assume that if the CGTO spanned the whole B-spline grid
  !> it would be sufficient to use the given quadrature rule (x,w,n) within each knot interval. For CGTOs whose radial extent is smaller than the range of the B-spline grid the number of quadrature points
  !> within each knot interval is expanded by the factor ceiling(n_int/R), where R is the extent of the CGTO (determined to coincide with the nearest larger knot) and n_int is the number of distinct 
  !> intervals of knots in the B-spline basis.
  subroutine radial_grid_CGTO(cgto,cms_bto,threshold,x,w,n,n_rng_knot,r_points,weights,n_total_points,R_grid,R_true)
     use gto_function_obj
     implicit none
     type(cgto_data), intent(inout) :: cgto
     type(bto_data), intent(inout) :: cms_bto
     integer, intent(in) :: n, n_rng_knot
     real(kind=cfp), intent(in) :: threshold, x(2*n+1), w(2*n+1)
     !OUTPUT variables:
     integer, intent(out) :: n_total_points
     real(kind=cfp), allocatable :: r_points(:), weights(:)
     real(kind=cfp), intent(out) :: R_grid, R_true

     integer :: n_points, err, interval, cnt, i, n_intervals, n_ranges, j, k, rng, n_int
     real(kind=cfp) :: d, A, B, R_lim, x_AB(2*n+1), w_AB(2*n+1), range_start_end(2,cms_bto%no_knots), delta
     type(contracted_gto_function) :: cgto_fn

         n_points = 2*n+1 !Number of quadrature points within each quadrature interval
!
!------- Estimate radius of the CGTO
         err = cgto_fn%init(cgto)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','cgto_fn init failed.',err,1)

         call cgto_fn%estimate_shell_radius(threshold,d,R_true)
         print *,'radius=',d,R_true
         if (R_true .le. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','CGTO is negligible on the BTO grid. This routine would not work properly.',1,1)

         write(stdout,'("CGTO negligible beyond radius [a.u.]: ",e25.15)') R_true
         write(stdout,'("CGTO threshold: ",e25.15)') threshold
!
!------- Determine the range of the quadrature grid so that the end point coincides with the nearest knot larger than R_true.
         range_start_end = 0.0_cfp
         n_ranges = 0
         R_grid = R_true
         do i=1,cms_bto%no_knots-1
            if (cms_bto%knots(i+1) .eq. cms_bto%knots(i)) cycle
            n_ranges = n_ranges + 1
            range_start_end(1:2,n_ranges) = (/cms_bto%knots(i),cms_bto%knots(i+1)/)
            if (cms_bto%knots(i+1) .ge. R_grid) then
               R_grid = cms_bto%knots(i+1)
               exit
            endif
         enddo !i
         R_grid = min(R_grid,cms_bto%B)

         write(stdout,'("Radial integration will be performed in the range [a.u.]: ",2e25.15)') cms_bto%A,R_grid
!
!------- Determine the number of distinct intervals of knots in the B-spline grid.
         n_int = 0
         do i=1,cms_bto%no_knots-1
            if (cms_bto%knots(i+1) .eq. cms_bto%knots(i)) cycle
            n_int = n_int+1
         enddo !i

         !Number of quadrature intervals dividing the radial interval between two consecutive knots.
         n_intervals = n_rng_knot*ceiling(n_int/R_grid)
         !print *,'intervals',n_intervals,n_int

         !Total number of quadrature points in the interval [cms_bto%A,R_grid].
         n_total_points = n_ranges*n_intervals*n_points
         !print *,'n_ranges',n_ranges

         if (allocated(r_points)) deallocate(r_points)
         if (allocated(weights)) deallocate(weights)
         allocate(r_points(n_total_points),weights(n_total_points),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','Memory allocation failed.',err,1)
!
!------- Prepare the quadrature grid for the interval [cms_bto%A,R] expanding the G-L quadrature within each range and each interval within the range.
         cnt = 0
         do rng=1,n_ranges
            delta = (range_start_end(2,rng)-range_start_end(1,rng))/n_intervals
            A = range_start_end(1,rng)
            do interval=1,n_intervals
               !Quadrature points and weights on the interval [A,B] with B = A + delta:
               B = A + delta
               if (interval .eq. n_intervals) B = range_start_end(2,rng)
               call gl_expand_A_B(x,w,n,r_points(cnt+1:cnt+n_points),weights(cnt+1:cnt+n_points),A,B)
               cnt = cnt + n_points
               A = A + delta
            enddo !interval
         enddo !rng
!
         if (cnt .ne. n_total_points) call xermsg('bto_gto_integrals_mod','BTO_GTO_one_electron_integrals','Error constructing the radial grid.',2,1)

  end subroutine radial_grid_CGTO

  !> Determines the radial grid needed to describe, within the accuracy threshold, integrals involving the CGTO and the B-splines on the given grid. We assume that if the CGTO spanned the whole B-spline grid
  !> it would be sufficient to use the given quadrature rule (x,w,n) within each knot interval. For CGTOs whose radial extent is smaller than the range of the B-spline grid the number of quadrature points
  !> within each knot interval is expanded by the factor ceiling(n_int/R), where R is the extent of the CGTO (determined to coincide with the nearest larger knot) and n_int is the number of distinct 
  !> intervals of knots in the B-spline basis.
  subroutine radial_grid_CGTO_pair(cgto_1,cgto_2,threshold,rmat_radius,x,w,n,n_rng_knot,r_points,weights,n_total_points,R_min,R_max)
     use gto_function_obj
     implicit none
     type(cgto_data), intent(inout) :: cgto_1, cgto_2
     integer, intent(in) :: n, n_rng_knot
     real(kind=cfp), intent(in) :: threshold, x(2*n+1), w(2*n+1), rmat_radius
     !OUTPUT variables:
     integer, intent(out) :: n_total_points
     real(kind=cfp), allocatable :: r_points(:), weights(:)
     real(kind=cfp), intent(out) :: R_min, R_max

     integer :: n_points, err, interval, cnt, i, n_intervals, n_ranges, rng, n_int
     real(kind=cfp) :: A, B, x_AB(2*n+1), w_AB(2*n+1), delta, val, radius_1, radius_2, R_1, R_2
     type(contracted_gto_function) :: cgto_1_fn, cgto_2_fn
     logical :: found

         n_points = 2*n+1 !Number of quadrature points within each quadrature interval
!
!------- Estimate the domain of the product of the two CGTOs
         err = cgto_1_fn%init(cgto_1)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','radial_grid_CGTO_pair','cgto_1_fn init failed.',err,1)
!
         err = cgto_2_fn%init(cgto_2)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','radial_grid_CGTO_pair','cgto_2_fn init failed.',err,1)
!
         R_1 = sqrt(dot_product(cgto_1%RA,cgto_1%RA))      !radial distance from CMS to the CGTO 1 center
         R_2 = sqrt(dot_product(cgto_2%RA,cgto_2%RA))      !radial distance from CMS to the CGTO 2 center

         !Start the estimate of R_max from the radial distance corresponding to the center of the CGTO further away from CMS
         radius_1 = max(R_1,R_2)-R_1
         radius_2 = max(R_1,R_2)-R_2
         do
            radius_1 = radius_1 + 0.1_cfp
            radius_2 = radius_2 + 0.1_cfp
            !Product of the radial parts of the CGTOs
            val = abs(cgto_1_fn%eval_rad(radius_1,0) * cgto_2_fn%eval_rad(radius_2,0))
            if (val .le. threshold) exit
         enddo

         !Distance from CMS beyond which this pair of CGTOs can be neglected.
         R_max = radius_1+R_1

         !Now estimate R_min
         radius_1 = abs(min(R_1,R_2)-R_1)
         radius_2 = abs(min(R_1,R_2)-R_2)
         do
            radius_1 = radius_1 + 0.1_cfp
            radius_2 = radius_2 + 0.1_cfp
            if (radius_1 > R_1) then
               radius_1 = R_1
               exit
            endif
            !Product of the radial parts of the CGTOs
            val = abs(cgto_1_fn%eval_rad(radius_1,0) * cgto_2_fn%eval_rad(radius_2,0))
            if (val .le. threshold) exit
         enddo 

         !Distance from CMS below which this pair of CGTOs can be neglected.
         R_min = R_1-radius_1

         write(stdout,'("Radial domain of the CGTO pair [a.u.]: ",2e25.15)') R_min,R_max
         write(stdout,'("CGTO pair threshold: ",e25.15)') threshold
         R_min = max(floor(R_min)*1.0_cfp,0.0_cfp)
         R_max = min(ceiling(R_max)*1.0_cfp,rmat_radius) !don't integrate further than up to R-matrix radius
         write(stdout,'("Radial integration will be performed in the range [a.u.]: ",2e25.15)') R_min,R_max
!
         n_int = nint(R_max-R_min) !number of intervals of length 1 a.u.

         !Number of quadrature intervals dividing the radial interval of length 1 a.u.
         n_intervals = n_rng_knot*n_int
         !print *,'intervals',n_intervals,n_int

         !Total number of quadrature points in the interval [R_min,R_max].
         n_total_points = n_intervals*n_points
         !print *,'n_ranges',n_ranges

         if (allocated(r_points)) deallocate(r_points)
         if (allocated(weights)) deallocate(weights)
         allocate(r_points(n_total_points),weights(n_total_points),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','radial_grid_CGTO_pair','Memory allocation failed.',err,1)
!
!------- Prepare the quadrature grid for the interval [R_min,R_max] expanding the G-L quadrature within each range and each interval within the range.
         cnt = 0
         delta = 1.0_cfp/(n_intervals*1.0_cfp)
         do rng=1,n_int
            A = R_min + rng-1.0_cfp
            do interval=1,n_rng_knot
               B = A + delta
               if (interval .eq. n_rng_knot) B = R_min + rng
               call gl_expand_A_B(x,w,n,r_points(cnt+1:cnt+n_points),weights(cnt+1:cnt+n_points),A,B)
               cnt = cnt + n_points
               A = A + delta
            enddo !interval
         enddo !rng
!
         if (cnt .ne. n_total_points) call xermsg('bto_gto_integrals_mod','radial_grid_CGTO_pair','Error constructing the radial grid.',2,1)

  end subroutine radial_grid_CGTO_pair

  subroutine BTO_GTO_two_el_BB_GG(cgto_shell_1,cgto_shell_2,cms_bto,integrals)
     use const, only: epsrel, epsabs
     use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
     use gto_function_obj
     use bspline_base, only: map_knots_to_grid
     use omp_lib
     implicit none
     real(kind=cfp), allocatable :: integrals(:)
     type(bto_data), intent(inout) :: cms_bto
     type(cgto_data), intent(inout) :: cgto_shell_1,cgto_shell_2

     integer, parameter :: n_rng_knot = 1 !Number of intervals in between each pair of knots within which the quadrature will be applied
     !Maximum L to use in the Legendre expansion of the Coulomb potential.
     integer :: max_l_legendre

     real(kind=cfp), allocatable :: cgto_1_x_2_pw(:,:,:,:) !(radial points,CGTO 2 m,CGTO 1 m,lm)
     real(kind=cfp), allocatable :: r1(:), w1(:), r2(:), w2(:), R_l_ints(:,:,:)
     real(kind=cfp), allocatable :: legendre_rad(:,:,:), R_lm(:,:,:), GG_r1(:,:,:,:), cpl_buf(:)
     real(kind=cfp) :: RA, threshold, R1_min, R1_max, integral
     integer :: n1_total_points, n2_total_points, cnt, i, ind, err, max_l_pw, max_l_aux, j, k, n_cgto_1_m, n_cgto_2_m, l2, m2, last_bspline, m_1_ind, m_2_ind, lm, l, m
     integer :: nres, n, block_total, iam, n_threads, stride, my_element, n_gq, n_pairs, m1m2_ind, l_min, l_max, ij, l1, m1
     integer, allocatable :: bspline_start_end_r2(:,:), thread_cnt(:)

        !For the (BB|GG) class the Legendre expansion is finite. The maximum L=2*BTO_L
        max_l_legendre = 2*cms_bto%l

        call cms_bto%init_grid
        call cms_bto%normalize

        call cgto_shell_1%normalize
        call cgto_shell_2%normalize

        RA = sqrt(dot_product(cgto_shell_1%RA,cgto_shell_1%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BB_GG','The CGTO 1 center must lie radially before the start of the BTO grid.',1,1)

        RA = sqrt(dot_product(cgto_shell_2%RA,cgto_shell_2%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BB_GG','The CGTO 2 center must lie radially before the start of the BTO grid.',2,1)

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs
!
!------ Determine the radial grid appropriate for integrations involving this CGTO shell. The resulting quadrature points and weights are in the arrays r1,w1.
        call radial_grid_CGTO_pair(cgto_shell_1,cgto_shell_2,threshold,cms_bto%B,x_10,w_10,n_10,n_rng_knot,r1,w1,n1_total_points,R1_min,R1_max)
        print *,'pair domain=',R1_min,R1_max
        !todo branch into multipole expansion since the domains for GG and BB don't overlap!!!
        if (R1_max .le. cms_bto%A) then
           print *, "multipole expansion should be used here; this routine is too expensive"
           stop "error"
        endif

        !todo calculate the r2 grid for integrations between pairs of B-splines and the radial Legendre function.
        n_gq = 2*n_10+1 !number of points of the G-L quadrature rule
        if (2*cms_bto%order > 2*n_gq-1)  call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BB_GG','The B-spline order is too large for the n_10 Gauss-Legendre quadrature rule. Decrease the B-spline order or increase the quadrature rule order.',3,1)
!
!------ Evaluate the B-splines on the r2-grid.
        call bspline_quadrature_grid(cms_bto,x_10,w_10,n_10,r2,w2,n2_total_points)
!
!------ Determine mapping of the endpoints of each B-spline with indices in the r2 array.
        last_bspline = cms_bto%n
        call map_knots_to_grid(cms_bto%knots,cms_bto%order,last_bspline,r2,bspline_start_end_r2)
!
!------ Precalculate the radial terms from the Legendre expansion needed to assemble the Y function on the r1 X r2 (2D) grid.
        call omp_calc_legendre_radial(max_l_legendre,r1,r2,legendre_rad)
!
!------ Calculate the radial integrals between each unique pair of B-splines and the radial Legendre resolution. The final integrals
!       over r2 are multiplied by the Jacobian in the r1 coordinate (Jac=r1**2) and the 1/r1 factor completing with Y the radial 
!       terms of the Legendre expansion. The values are also multiplied by the weights w1 for integrals over r1.
        call omp_radial_integrals_legendre_bsplines(cms_bto,bspline_start_end_r2,legendre_rad,max_l_legendre,r2,w2,&
                                               &n2_total_points,r1,w1,n1_total_points,R_l_ints)
!
!------ Precalculate all coupling coefficients needed for the calculations.
        max_l_pw = cms_bto%l+max_l_legendre
        max_l_aux = max_l_pw+max(cgto_shell_1%l,cgto_shell_2%l)

        !The two routines below must be called before any OpenMP region to
        !ensure thread safety of resh, solh.
        call calc_resh_coefficients(max_l_aux)
        call calc_solh_coefficients(max_l_aux)

        !Precalculate all couplings needed
        print *,'gaunt'
        call cpl%prec_cgaunt(max_l_aux)
!
!------ Calculate the projections of the CGTO shell 1 on the real spherical harmonics (i.e. PW-expansion) for all radial points of the grid 1.
        print *,'pw 1'
        call omp_calculate_CGTO_pair_pw_coefficients_analytic(max_l_pw,cgto_shell_1,cgto_shell_2,r1,cgto_1_x_2_pw)
!        call calculate_CGTO_pair_pw_coefficients_numerical(max_l_pw,cgto_shell_1,cgto_shell_2,r1,cgto_1_x_2_pw,epsrel,epsabs)
        print *,'done'
!
!------ Assemble the 2-electron integrals:
        n_pairs = cms_bto%n*(cms_bto%n+1)/2 !number of unique pairs of B-splines.
        n_cgto_1_m = 2*cgto_shell_1%l+1
        n_cgto_2_m = 2*cgto_shell_2%l+1
!
        l_max = 2*cms_bto%l
        allocate(GG_r1(n1_total_points,0:(2*cms_bto%l),n_cgto_2_m,n_cgto_1_m),R_lm(n_cgto_2_m,n_cgto_1_m,n_pairs),&
                &cpl_buf((l_max+1)**2),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BB_GG','Memory allocation failed.',err,1)
!
        do l2=0,cms_bto%l
           do m2=-l2,l2
              print *,'l2,m2',l2,m2
              do l1=0,cms_bto%l
                 do m1=-l1,l1

                    !Preload the coupling coefficients:
                    l_min = abs(l2-l1)
                    l_max = l2+l1
                    do l=l_min,l_max
                       lm = l*l+l+1
                       do m=-l,l
                          cpl_buf(lm+m) = cpl%rgaunt(l2,l1,l,m2,m1,m)
                       enddo !m
                    enddo !l

                    !Multiply-in the coupling coefficients with the PW expansion of GG and sum over the m-values:
                    !m-values of the CGTO of the first electron
                    do m_1_ind=1,n_cgto_1_m
                       !m-values of the CGTO of the first electron
                       do m_2_ind=1,n_cgto_2_m
                          m1m2_ind = m_2_ind + n_cgto_2_m*(m_1_ind-1)
                          do l=l_min,l_max
                             GG_r1(1:n1_total_points,l,m_2_ind,m_1_ind) = 0.0_cfp
                             lm = l*l+l+1
                             do m=-l,l
                                if (cpl_buf(lm+m) .ne. 0.0_cfp) then
                                   do i=1,n1_total_points
                                      GG_r1(i,l,m_2_ind,m_1_ind) = GG_r1(i,l,m_2_ind,m_1_ind) + cpl_buf(lm+m)*cgto_1_x_2_pw(i,m_2_ind,m_1_ind,lm+m)
                                   enddo !i
                                endif
                             enddo !m
                          enddo !l
                       enddo !m_2_ind
                    enddo !m_1_ind

                    !Perform the integrals over r1: obtain the final 2-electron integrals.
                    R_lm = 0.0_cfp
                    !unique pairs of radial B-splines
                    do ij=1,n_pairs
                       do m_1_ind=1,n_cgto_1_m
                          !m-values of the CGTO of the first electron
                          do m_2_ind=1,n_cgto_2_m
                             !Sum over the Legendre terms:
                             do l=l_min,l_max
                                !Integral over r1:
                                R_lm(m_2_ind,m_1_ind,ij) = R_lm(m_2_ind,m_1_ind,ij) + &
                                &sum(GG_r1(1:n1_total_points,l,m_2_ind,m_1_ind)*R_l_ints(1:n1_total_points,l,ij))
                             enddo !l
                             if (R_lm(m_2_ind,m_1_ind,ij) .ne. 0.0_cfp) &
                             &write(stdout,'(4i3,1i4,2i3,e25.15)') l2,m2,l1,m1,ij,m_2_ind,m_1_ind, R_lm(m_2_ind,m_1_ind,ij)
                          enddo !m_2_ind
                       enddo !m_1_ind
                    enddo !ij

                 enddo !m1
              enddo !l1
           enddo !m2
        enddo !l2

  end subroutine BTO_GTO_two_el_BB_GG

  !> \todo Parts of this routine can be ran only once such as the evaluation of the B-spline pairs on the r2 quadrature grid.
  subroutine omp_radial_integrals_legendre_bsplines(cms_bto,bspline_start_end_r2,legendre_rad,max_l_legendre,r2,w2,&
                                               &n2_total_points,r1,w1,n1_total_points,R_l_ints)
     use bspline_base, only: bvalu
     use phys_const, only: fourpi
     use omp_lib
     implicit none
     type(bto_data), intent(inout) :: cms_bto
     integer, intent(in) :: max_l_legendre, n2_total_points, n1_total_points, bspline_start_end_r2(:,:)
     real(kind=cfp), allocatable :: legendre_rad(:,:,:), r2(:), w2(:), r1(:), w1(:)
     !OUTPUT:
     real(kind=cfp), allocatable :: R_l_ints(:,:,:) !r1,l,(i,j); (i,j) = unique B-spline pair i*(i-1)/2+j,i.ge.j

     real(kind=cfp) :: bto_norm
     real(kind=cfp), allocatable :: B_vals(:,:), bspline_prod(:)
     integer :: ind, i, j, err, j_start, j_end, l, cnt, ind_i, ind_j, iam, n_threads

        if (allocated(R_l_ints)) deallocate(R_l_ints)

        allocate(B_vals(n2_total_points,cms_bto%n),R_l_ints(n1_total_points,0:max(max_l_legendre,1),cms_bto%n*(cms_bto%n+1)/2),&
        &bspline_prod(n2_total_points),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','radial_integrals_legendre_bsplines','Memory allocation failed.',err,1)

        B_vals = 0.0_cfp
        cms_bto%bcoef = 0.0_cfp
        do ind=1,cms_bto%n
           cms_bto%bcoef(ind) = 1.0_cfp

           !Calculate norm of the radial B-spline
           cms_bto%ind = ind
           cms_bto%r1 = cms_bto%knots(ind)
           cms_bto%r2 = cms_bto%knots(ind+cms_bto%order)
           call cms_bto%normalize
           bto_norm = cms_bto%norm

           !Evaluate it on the grid
           do i=bspline_start_end_r2(1,ind),bspline_start_end_r2(2,ind)
              B_vals(i,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,0,r2(i),cms_bto%inbv,cms_bto%work)
           enddo !i

           cms_bto%bcoef(ind) = 0.0_cfp !clean-up for the next
        enddo !ind

        !Integrate over the radial resolution and the product of the B-splines:
        !$OMP PARALLEL DEFAULT(NONE) PRIVATE(cnt,ind_i,ind_j,j_start,j_end,bspline_prod,l,i,iam,n_threads) &
        !$OMP & SHARED(cms_bto,bspline_start_end_r2,max_l_legendre,n1_total_points,R_l_ints,legendre_rad,B_vals,r1,w1,w2,stdout)
        iam = omp_get_thread_num()
        n_threads = omp_get_num_threads()
        cnt = 0
        do ind_i=1,cms_bto%n
           do ind_j=1,ind_i
              cnt = cnt + 1
              if (mod(cnt,n_threads) .ne. iam) cycle !work redistribution
              !Select the r2 radial range where both B-splines overlap (if they do):
              j_start = max(bspline_start_end_r2(1,ind_i),bspline_start_end_r2(1,ind_j))
              j_end   = min(bspline_start_end_r2(2,ind_i),bspline_start_end_r2(2,ind_j))
              if (j_start > j_end) cycle
              bspline_prod(j_start:j_end) = B_vals(j_start:j_end,ind_i)*B_vals(j_start:j_end,ind_j)*w2(j_start:j_end)
              do l=0,max_l_legendre
                 do i=1,n1_total_points
                    !For each r1 we compute the integral over r2 and multiply by the corresponding r1 and w1:
                    R_l_ints(i,l,cnt) = sum(bspline_prod(j_start:j_end)*legendre_rad(j_start:j_end,i,l))*r1(i)*w1(i)
                    !write(stdout,*) R_l_ints(i,l,cnt),i,l,cnt
                 enddo !i
              enddo !l
           enddo !ind_j
        enddo !ind_i
        !$OMP END PARALLEL

  end subroutine omp_radial_integrals_legendre_bsplines

  !> Construct the G-L quadrature rule of the given order between each consecutive pair of knots.
  subroutine bspline_quadrature_grid(cms_bto,x,w,n,r_points,weights,n_total_points)
     implicit none
     type(bto_data), intent(inout) :: cms_bto
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: w(2*n+1), x(2*n+1)
     real(kind=cfp), allocatable :: r_points(:), weights(:)
     integer, intent(out) :: n_total_points

     integer :: ind, cnt, n_points, err
     real(kind=cfp) :: A, B

        n_points = 2*n+1

        if (allocated(r_points)) deallocate(r_points)
        if (allocated(weights)) deallocate(weights)

        n_total_points = 0
        do ind=1,cms_bto%no_knots-1
           A = cms_bto%knots(ind)
           B = cms_bto%knots(ind+1)
           if (B > A) n_total_points = n_total_points + n_points
        enddo !ind

        allocate(r_points(n_total_points),weights(n_total_points),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','bspline_quadrature_grid','Memory allocation failed.',err,1)

        cnt = 0
        do ind=1,cms_bto%no_knots-1
           A = cms_bto%knots(ind)
           B = cms_bto%knots(ind+1)
           if (B > A) then
              call gl_expand_A_B(x,w,n,r_points(cnt+1:cnt+n_points),weights(cnt+1:cnt+n_points),A,B)
              cnt = cnt + n_points
           endif
        enddo !ind
        
  end subroutine bspline_quadrature_grid

  subroutine BTO_GTO_two_el_BG_GG(cgto_shell_1,cgto_shell_2,cms_bto,cgto_shell_3,integrals)
     use const, only: epsrel, epsabs
     use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
     use bspline_base, only: bvalu, map_knots_to_grid
     use gto_function_obj
     use omp_lib
     implicit none
     real(kind=cfp), allocatable :: integrals(:)
     type(bto_data), intent(inout) :: cms_bto
     type(cgto_data), intent(inout) :: cgto_shell_1,cgto_shell_2,cgto_shell_3

     integer, parameter :: n_rng_knot = 1 !Number of intervals in between each pair of knots within which the quadrature will be applied
     !Maximum L to use in the Legendre expansion of the Coulomb potential.
     integer, parameter :: max_l_legendre = 40

     real(kind=cfp), allocatable :: cgto_1_x_2_pw(:,:,:,:), cgto_pw_shell_3(:,:,:) !(radial points,CGTO 2 m,CGTO 1 m,lm)
     real(kind=cfp), allocatable :: r1(:), w1(:), r2(:), w2(:)
     real(kind=cfp), allocatable :: legendre_rad(:,:,:), Y(:,:,:,:), Y_mp(:,:,:,:), B_w_r2(:,:), R_lm(:,:,:,:), epstab(:), res3la(:)
     real(kind=cfp) :: RA, threshold, R1_min, R1_max, R2_max, r, val, res, abserr, relative_precision, overflow, prec_goal, extrapolated, integral, abserr_lowest, R2_true
     integer :: n1_total_points, n2_total_points, cnt, i, ind, err, max_l_pw, max_l_aux, j, k, n_cgto_1_m, n_cgto_2_m, l2, m2, last_bspline, m_1_ind, m_2_ind, lm, l, m
     integer :: nres, n, n_wynn, to_converge, block_total, iam, n_threads, stride, my_element, n_cgto_3_m, m_3_ind
     integer, allocatable :: bspline_start_end_r2(:,:), thread_cnt(:)
     logical :: no_knot

        call cms_bto%init_grid
        call cms_bto%normalize

        call cgto_shell_1%normalize
        call cgto_shell_2%normalize
        call cgto_shell_3%normalize

        RA = sqrt(dot_product(cgto_shell_1%RA,cgto_shell_1%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','The CGTO 1 center must lie radially before the start of the BTO grid.',1,1)

        RA = sqrt(dot_product(cgto_shell_2%RA,cgto_shell_2%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','The CGTO 2 center must lie radially before the start of the BTO grid.',2,1)

        RA = sqrt(dot_product(cgto_shell_3%RA,cgto_shell_3%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','The CGTO 3 center must lie radially before the start of the BTO grid.',3,1)

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs
!
!------ Determine the radial grid appropriate for integrations involving this CGTO shell. The resulting quadrature points and weights are in the arrays r1,w1.
        call radial_grid_CGTO_pair(cgto_shell_1,cgto_shell_2,threshold,cms_bto%B,x_10,w_10,n_10,n_rng_knot,r1,w1,n1_total_points,R1_min,R1_max)
        print *,'pair domain=',R1_min,R1_max
        !todo branch into multipole expansion since the domains for GG and GB don't overlap!!!
        if (R1_max .le. cms_bto%A) then
           print *, "multipole expansion should be used here; this routine won't work"
           stop "error"
        endif
!
!------ Radial grid for CGTO 3 although we only need the extent of the grid: R2_max.
        call radial_grid_CGTO(cgto_shell_3,cms_bto,threshold,x_10,w_10,n_10,n_rng_knot,r2,w2,n2_total_points,R2_max,R2_true)
        if (R2_true < cms_bto%A) then
           print *,"the integrals are zero due to shell 3"
           stop "not implemented fully"
        endif
!
!------ The end point for the radial integration of the Y function of the second electron with the radial density of the first electron.
        R = R2_max
!
!------ Determine index of the last B-spline for which the integrals will be computed. This is done including all B-splines which start before the significance radius R.
        last_bspline = last_bspline_before_R(R2_true,cms_bto)
        print *,'last_bspline',last_bspline
        if (last_bspline .eq. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Error either in the routine radial_grid_CGTO or in B-spline grid.',3,1)
!
!------ Construct the dense grid r2,w2 needed to calculate the Y function. This grid is constructed by constructing the G-L grid of order 7 in between each pair of r1 quadrature points.
!       If there is a knot lying in between the two r1 quadrature points then the interval is divided in two: [r1(i),knot] and [knot,r1(i+1)]. This ensures that when calculating the Y function
!       corresponding to the radial B-splines that the intervals of the integrand close to the end points of the B-splines are always accurately integrated.
        call construct_dense_grid_from_sparse_grid(cms_bto,R,x_10,w_10,n_10,n_rng_knot,r1,n1_total_points,x_7,w_7,n_7,n2_total_points,r2,w2)
!
!------ Precalculate all coupling coefficients needed for the calculations.
        max_l_pw = cms_bto%l+max_l_legendre
        max_l_aux = max_l_pw+max(cgto_shell_1%l,cgto_shell_2%l,cgto_shell_3%l)

        !The two routines below must be called before any OpenMP region to ensure thread safety of resh, solh.
        call calc_resh_coefficients(max_l_aux)
        call calc_solh_coefficients(max_l_aux)

        !Precalculate all couplings needed
        print *,'gaunt'
        call cpl%prec_cgaunt(max_l_aux)
!
!------ Calculate the projections of the CGTO shell 1 on the real spherical harmonics (i.e. PW-expansion) for all radial points of the grid 1.
        print *,'pw 1'
        call omp_calculate_CGTO_pair_pw_coefficients_analytic(max_l_pw,cgto_shell_1,cgto_shell_2,r1,cgto_1_x_2_pw)
        print *,'done'
!
!------ Calculate the projections of the CGTO shell 2 on the real spherical harmonics (i.e. PW-expansion) for all radial points of the grid 2. This is the dense grid.
        print *,'pw 2'
        call omp_calculate_CGTO_pw_coefficients_analytic(max_l_pw,cgto_shell_3,r2,cgto_pw_shell_3)
        print *,'done',size(r2)
!
!------ Determine mapping of the endpoints of each B-spline with indices in the r2 array.
        call map_knots_to_grid(cms_bto%knots,cms_bto%order,last_bspline,r2,bspline_start_end_r2)
!
!------ Multiply w2 by Jacobian*1/r = r: this is the cancellation of the 1/r factor coming from the radial part of the BTO (B(r)/r).
        do i=1,n2_total_points
           w2(i) = w2(i)*r2(i)
        enddo !i
!
!------ Precalculate the B-splines on the r2 quadrature grid and multiply them by their norms and quadrature weights w2
        call eval_bspline_on_grid_times_weights(cms_bto,last_bspline,bspline_start_end_r2,r2,w2,B_w_r2)
        deallocate(w2)
!
!------ Multiply w1 by Jacobian*1/r = r: this is the cancellation of the 1/r factor in front of the Y function in the integral over r1 with the Jacobian (= r1**2).
        do i=1,n1_total_points
           w1(i) = w1(i)*r1(i)
        enddo !i
!
!------ Precalculate the radial terms from the Legendre expansion (4*pi/(2*l+1)*r<^l/r>^(l+1)) on the r1 X r2 (2D) grid.
        print *,'rad legendre'
        call omp_calc_legendre_radial(max_l_legendre,r1,r2,legendre_rad)
        print *,'done'
!
!------ Big loop over angular numbers of the B-splines of the second electron
        n_cgto_1_m = 2*cgto_shell_1%l+1
        n_cgto_2_m = 2*cgto_shell_2%l+1
        n_cgto_3_m = 2*cgto_shell_3%l+1
!
        allocate(R_lm(last_bspline,n_cgto_3_m,n_cgto_2_m,n_cgto_1_m),epstab(52),res3la(3),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Memory allocation 5 failed.',err,1)
        prec_goal = F1MACH(4,cfp_dummy)*100
        overflow = F1MACH(2,cfp_dummy)
        block_total = n_cgto_1_m*n_cgto_2_m*n_cgto_3_m*last_bspline !number of integrals in the block corresponding to the BTO angular numbers (l2,m2).
        stride = (2*cache_line_size)*8/bit_size(i)
        print *,'prec_goal',prec_goal
        !$OMP PARALLEL DEFAULT(NONE) PRIVATE(n_threads,iam,err) SHARED(thread_cnt,stride)
        !$OMP SINGLE
        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()
        !The array thread_cnt will be used below to store the number of integrals, within each block, converged by each thread.
        allocate(thread_cnt(0:max(n_threads-1,1)*stride),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Memory allocation 6 failed.',err,1)
        !$OMP END SINGLE
        !$OMP END PARALLEL
!
        !Loop over angular numbers of the BTO of the second electron
        do l2=0,cms_bto%l
           do m2=-l2,l2
              print *,'Y',l2,m2
              call omp_calculate_Y_BG(cms_bto,last_bspline,l2,m2,cgto_pw_shell_3,bspline_start_end_r2,max_l_legendre,n1_total_points,n2_total_points,n_cgto_3_m,legendre_rad,B_w_r2,Y_mp,Y,r1)
              !$OMP PARALLEL DEFAULT(NONE) PRIVATE(to_converge,m_1_ind,ind,i,m_2_ind,m_3_ind,l,m,lm,n,n_wynn,res,abserr,epstab,res3la,nres,n_threads,iam,relative_precision,&
              !$OMP & my_element,err,extrapolated,integral,abserr_lowest) &
              !$OMP & SHARED(cms_bto,R_lm,n_cgto_1_m,last_bspline,n_cgto_2_m,n_cgto_3_m,B_w_r2,cgto_1_x_2_pw,Y,overflow,prec_goal,block_total,stdout,l2,m2,stride,thread_cnt,n1_total_points,w1)
              n_threads = omp_get_num_threads()
              iam = omp_get_thread_num()
              my_element = stride*iam !there is a separation of lenght 'stride' elements between the values set by each element: this prevents false sharing between threads
              if (size(thread_cnt) < max(n_threads-1,1)*stride+1) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Number of threads used for allocation of thread_cnt &
                 &is smaller than the number of threads used during computation.',6,1)
              !$OMP SINGLE
              R_lm = 0.0_cfp
              thread_cnt = 0
              !$OMP END SINGLE
              !m-values of the CGTO of the first electron
              do m_1_ind=1,n_cgto_1_m
                !m-values of the CGTO of the first electron
                 do m_2_ind=1,n_cgto_2_m
                    i = m_2_ind + n_cgto_2_m*(m_1_ind-1) 
                    if (mod(i,n_threads) .ne. iam) cycle !work redistribution
                    !m-values of the CGTO of the second electron
                    do m_3_ind=1,n_cgto_3_m
                       !radial B-splines of the second electron
                       do ind=1,last_bspline
                          !Prepare the epsilon table
                          epstab = 0.0_cfp; res3la = 0.0_cfp; nres = 0
                          !Loop over l,m-angular numbers in the Legendre expansion
                          !write(111,'(2i4,2i4,4i3)') l2,m2, l1,m1, ind_2,m_2_ind,ind,m_1_ind
                          integral = 0.0_cfp !This value is the integral as obtained by straight summation of the Legendre terms
                          do l=0,min(max_l_legendre,49)

                             lm = l*l+l+1
                             do m=-l,l
                                !Integrate over r1 the product of the Y function of the second electron with the radial density of the first electron.
                                integral = integral + sum(cgto_1_x_2_pw(1:n1_total_points,m_2_ind,m_1_ind,lm+m)*Y(1:n1_total_points,lm+m,ind,m_3_ind)*w1(1:n1_total_points))
                             enddo !m
   
                             !write(111,'(i,e25.15)') l,integral
   
                             !Wynn's extrapolation: we need at least 3 points to start the algorithm
                             n = l+1
                             if (n .ge. 3) then
                                !The value of n_wynn can get changed by dqelg if numerical difficulties are encountered:
                                !typically this will happen for very small (negligible) values of the final integrals.
                                n_wynn = n
                                epstab(n) = integral
                                call dqelg(n_wynn, epstab(1:52), res, abserr, res3la(1:3), nres)
                                if (res .eq. overflow) then !The extrapolated value overflows - we cannot improve on the result by using more terms so terminate the expansion.
                                   thread_cnt(my_element) = thread_cnt(my_element)+1
                                   exit !converged
                                elseif (abserr < abserr_lowest) then !Save the extrapolated value as the estimate of the final integral only if it is a more precise estimate.
                                   R_lm(ind,m_3_ind,m_2_ind,m_1_ind) = res
                                   abserr_lowest = abserr
                                   !Estimate relative error in the extrapolated integral using the previous value res3la(2).
                                   relative_precision = abs((res-res3la(2))/res)
                                   !write(111,'("extrapolation",2e25.15)') R_lm(ind,m_3_ind,m_2_ind,m_1_ind), abserr
                                endif
                             elseif (n .eq. 2) then !First estimate of the final integrals
                                R_lm(ind,m_3_ind,m_2_ind,m_1_ind) = integral
                                abserr_lowest = overflow
                                relative_precision = 1.0_cfp
                             endif

                          enddo !l
                          !Ignore convergence of the integrals which are too small.
                          if (relative_precision .le. prec_goal .or. abs(R_lm(ind,m_3_ind,m_2_ind,m_1_ind)) < epsabs) thread_cnt(my_element) = thread_cnt(my_element)+1
                       enddo !ind
                    enddo !m_3_ind
                 enddo !m_2_ind
              enddo !m_1_ind
              !$OMP BARRIER

              !$OMP SINGLE
              do m_1_ind=1,n_cgto_1_m
                 do m_2_ind=1,n_cgto_2_m
                    do m_3_ind=1,n_cgto_3_m
                       do ind=1,last_bspline 
                          write(stdout,'(2i4,4i3,e25.15)') l2,m2, ind,m_3_ind,m_2_ind,m_1_ind, R_lm(ind,m_3_ind,m_2_ind,m_1_ind)
                       enddo !ind_2
                    enddo !m_2_ind
                 enddo !ind
              enddo !m_1_ind
              to_converge = block_total-sum(thread_cnt)
              write(stdout,'("not converged",2i)') to_converge,block_total
              !$OMP END SINGLE

              !$OMP END PARALLEL
           enddo !m2
        enddo !l2

  end subroutine BTO_GTO_two_el_BG_GG

  !> Takes the quadrature grid r_sparse and constructs from it another grid which contains 2*n+1 quadrature points of the rule x,w in between the quadrature points of the r_sparse grid. The resulting
  !> quadrature grid is in r_dense,w_dense and has n_dense_total_points points. The dense grid extends only up to the r_dense point not greater than R_max.
  subroutine construct_dense_grid_from_sparse_grid(cms_bto,R_max,x_t,w_t,n_t,n_rng_knot,r_sparse,n_sparse_total_points,x,w,n,n_dense_total_points,r_dense,w_dense)
     implicit none
     integer, intent(in) :: n_sparse_total_points, n, n_rng_knot, n_t
     real(kind=cfp), intent(in) :: R_max, r_sparse(n_sparse_total_points), x(2*n+1), w(2*n+1), x_t(2*n_t+1), w_t(2*n_t+1)
     type(bto_data), intent(in) :: cms_bto
     !OUTPUT:
     integer, intent(out) :: n_dense_total_points
     real(kind=cfp), allocatable :: r_dense(:), w_dense(:)

     integer :: n_points, i, err, ind, i_start
     logical :: no_knot
     real(kind=cfp) :: r, delta, B

        !We have to make sure that the value R_max lies beyond the starting point of the r_sparse grid otherwise no points for the r_dense grid would be generated.
        if (R_max < r_sparse(1)) call xermsg('bto_gto_integrals_mod','construct_dense_grid_from_sparse_grid','Error in input data: R_max must be larger than r_sparse(1).',1,1)

        !Find the first point r_sparse(i_start) such that r_sparse(i_start) > cms_bto%A. This guarantees that the r_dense quadrature is only over the radial range r .ge. cms_bto%A.
        i_start = 0
        do
           i_start = i_start + 1
           if (i_start > n_sparse_total_points) call xermsg('bto_gto_integrals_mod','construct_dense_grid_from_sparse_grid','Error in input data: r_sparse doesnt extent into B-spline knots.',2,1)
           if (r_sparse(i_start) > cms_bto%A) exit !we've found the starting point which can be used to build the quadrature r_dense.
        enddo !i

        n_points = 2*n+1
        !Initial sweep to determine the total number of points
        n_dense_total_points = 0
        r = cms_bto%A
        do i=i_start,n_sparse_total_points
           if (r > R_max) exit !We integrate only up to r_sparse = R_max
           no_knot = .true.
           do ind=1,cms_bto%no_knots
              !Test if there is a knot in between r and r_sparse(i):
              if (cms_bto%knots(ind) > r .and. cms_bto%knots(ind) < r_sparse(i)) then
                 !Rule for the intervals [r,knots(ind)] and [knots(ind),r_sparse(i)]. 
                 !This ensures that the Y function calculated for each B-spline will include properly the behavior of the B-splines around their end points.
                 n_dense_total_points = n_dense_total_points + n_points
                 n_dense_total_points = n_dense_total_points + n_points
                 no_knot = .false.
                 exit
              endif
           enddo
           !Rule for the interval [r,r_sparse(i)]
           if (no_knot) then
              n_dense_total_points = n_dense_total_points + n_points
           endif
           r = r_sparse(i)
        enddo !i
        !Add the points in the interval [r,R_max]:
        n_points = 2*n_t+1
        do
           if (r > R_max) exit !We integrate only up to r = R_max
           
           !Rules for the interval [r,knots(ind)]
           no_knot = .true.
           do ind=1,cms_bto%no_knots
              !Find the nearest larger knot:
              if (cms_bto%knots(ind) > r) then
                 n_dense_total_points = n_dense_total_points + n_points*n_rng_knot
                 i = ind
                 no_knot = .false.
                 exit
              endif
           enddo
           if (no_knot) exit !there is no knot larger than the current r
           r = cms_bto%knots(i)
           if (r .ge. R_max) exit !We integrate only up to r = R_max
        enddo

        !Now construct the actual rules
        if (allocated(r_dense)) deallocate(r_dense)
        if (allocated(w_dense)) deallocate(w_dense)
        allocate(r_dense(n_dense_total_points),w_dense(n_dense_total_points),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','construct_dense_grid_from_sparse_grid','Memory allocation failed.',err,1)
!
        n_points = 2*n+1
        n_dense_total_points = 0
        r = cms_bto%A
        do i=i_start,n_sparse_total_points
           if (r > R_max) exit !We integrate only up to r_sparse = R_max
           no_knot = .true.
           do ind=1,cms_bto%no_knots
              !Test if there is a knot in between r and r_sparse(i):
              if (cms_bto%knots(ind) > r .and. cms_bto%knots(ind) < r_sparse(i)) then
                 !Rule for the intervals [r,knots(ind)] and [knots(ind),r_sparse(i)]. 
                 !This ensures that the Y function calculated for each B-spline will include properly the behavior of the B-splines around their end points.
                 call gl_expand_A_B(x,w,n,r_dense(n_dense_total_points+1:n_dense_total_points+n_points),w_dense(n_dense_total_points+1:n_dense_total_points+n_points),r,cms_bto%knots(ind))
                 n_dense_total_points = n_dense_total_points + n_points
                 call gl_expand_A_B(x,w,n,r_dense(n_dense_total_points+1:n_dense_total_points+n_points),w_dense(n_dense_total_points+1:n_dense_total_points+n_points),cms_bto%knots(ind),r_sparse(i))
                 n_dense_total_points = n_dense_total_points + n_points
                 no_knot = .false.
                 exit
              endif
           enddo !ind
           !Rule for the interval [r,r_sparse(i)]
           if (no_knot) then
              call gl_expand_A_B(x,w,n,r_dense(n_dense_total_points+1:n_dense_total_points+n_points),w_dense(n_dense_total_points+1:n_dense_total_points+n_points),r,r_sparse(i))
              n_dense_total_points = n_dense_total_points + n_points
           endif
           r = r_sparse(i)
        enddo !i
        !Add the points in the interval [r,R_max]:
        n_points = 2*n_t+1
        do
           if (r > R_max) exit !We integrate only up to r = R_max

           !Rules for the interval [r,knots(ind)]
           no_knot = .true.
           do ind=1,cms_bto%no_knots
              !Find the nearest larger knot and divide the interval between the knots into n_ngr_knot ranges:
              if (cms_bto%knots(ind) > r) then
                 delta = (cms_bto%knots(ind)-r)/n_rng_knot
                 do i=1,n_rng_knot
                    B = r + delta*i
                    call gl_expand_A_B(x_t,w_t,n_t,r_dense(n_dense_total_points+1:n_dense_total_points+n_points),w_dense(n_dense_total_points+1:n_dense_total_points+n_points),r,B)
                    n_dense_total_points = n_dense_total_points + n_points
                    r = B
                 enddo
                 no_knot = .false.
                 i = ind
                 exit
              endif
           enddo
           if (no_knot) exit !there is no knot larger than the current r
           r = cms_bto%knots(i)
           if (r .ge. R_max) exit !We integrate only up to r = R_max
        enddo

  end subroutine construct_dense_grid_from_sparse_grid

  !> Returns the index of the last B-spline in the basis which starts at a point not greater than R.
  function last_bspline_before_R(R,cms_bto)
     implicit none
     real(kind=cfp), intent(in) :: R
     type(bto_data), intent(in) :: cms_bto
     integer :: last_bspline_before_R

     integer :: ind

        ind = 0
        do
           ind = ind + 1
           if (cms_bto%knots(ind) .ge. R .or. ind .eq. cms_bto%n) then
              last_bspline_before_R = max(1,ind-1)
              exit
           endif
        enddo !ind

  end function last_bspline_before_R

  !> Evaluates the B-splines on the given quadrature grid multiplying them by their norms and the quadrature weights.
  subroutine eval_bspline_on_grid_times_weights(cms_bto,last_bspline,bspline_start_end_r,r,w,B_w)
     use bspline_base, only: bvalu
     implicit none
     integer, intent(in) :: last_bspline, bspline_start_end_r(2,last_bspline)
     type(bto_data), intent(inout) :: cms_bto
     real(kind=cfp), intent(in) :: r(:), w(:)
     real(kind=cfp), allocatable :: B_w(:,:)

     integer :: err, ind, i, n_total_points
     real(kind=cfp) :: bto_norm

        n_total_points = size(r)

        if (allocated(B_w)) deallocate(B_w)
        allocate(B_w(n_total_points,last_bspline),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','eval_bspline_on_grid_times_weights','Memory allocation failed.',err,1)

        B_w = 0.0_cfp
        cms_bto%bcoef = 0.0_cfp
        do ind=1,last_bspline
           cms_bto%bcoef(ind) = 1.0_cfp

           !Calculate norm of the radial B-spline
           cms_bto%ind = ind
           cms_bto%r1 = cms_bto%knots(ind)
           cms_bto%r2 = cms_bto%knots(ind+cms_bto%order)
           call cms_bto%normalize
           bto_norm = cms_bto%norm

           !Evaluate it on the grid and multiply by the quadrature weight
           do i=bspline_start_end_r(1,ind),bspline_start_end_r(2,ind)
              B_w(i,ind) = bto_norm*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,0,r(i),cms_bto%inbv,cms_bto%work)*w(i)
           enddo !i

           cms_bto%bcoef(ind) = 0.0_cfp !clean-up for the next
        enddo !ind

  end subroutine eval_bspline_on_grid_times_weights

  !todo cgto_shell_1%l should never be smaller than cgto_shell_2%l for computational efficiency
  subroutine BTO_GTO_two_el_BG_BG(cgto_shell_1,cgto_shell_2,cms_bto,integrals)
     use const, only: epsrel, epsabs
     use general_quadrature, only: n_10,w_10,x_10, n_7,w_7,x_7
     use bspline_base, only: bvalu, map_knots_to_grid
     use gto_function_obj
     use omp_lib
     implicit none
     real(kind=cfp), allocatable :: integrals(:)
     type(bto_data), intent(inout) :: cms_bto
     type(cgto_data), intent(inout) :: cgto_shell_1,cgto_shell_2

     integer, parameter :: n_rng_knot = 1 !Number of intervals in between each pair of knots within which the quadrature will be applied
     !Maximum L to use in the Legendre expansion of the Coulomb potential.
     integer, parameter :: max_l_legendre = 40

     real(kind=cfp), allocatable :: cgto_pw_shell_1(:,:,:), cgto_pw_shell_2(:,:,:) !(radial points,CGTO m,lm)
     real(kind=cfp), allocatable :: gaunt_cgto_pw_shell_1(:,:,:,:), r1(:), w1(:), r2(:), w2(:)
     real(kind=cfp), allocatable :: legendre_rad(:,:,:), B_w_r2(:,:), Y(:,:,:,:), Y_mp(:,:,:,:), B_w_r1(:,:), R_lm(:,:,:,:), epstab(:), res3la(:)
     real(kind=cfp) :: RA, threshold, R1_max, R2_max, r, val, res, abserr, relative_precision, overflow, prec_goal, extrapolated, integral, abserr_lowest, R1_true, R2_true
     integer :: n1_total_points, n2_total_points, i, err, max_l_pw, max_l_aux, n_cgto_1_m, n_cgto_2_m, l2, m2, last_bspline_1, last_bspline_2, l1, m1, ind_1, ind_2, m_1_ind, m_2_ind, lm, l1_m1, l, m
     integer :: nres, n, n_wynn, to_converge, block_total, iam, n_threads, stride, my_element
     integer, allocatable :: bspline_start_end_r2(:,:), bspline_start_end_r1(:,:), thread_cnt(:)
     logical :: no_knot

        call cms_bto%init_grid
        call cms_bto%normalize

        call cgto_shell_1%normalize
        call cgto_shell_2%normalize

        RA = sqrt(dot_product(cgto_shell_1%RA,cgto_shell_1%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','The CGTO center must lie radially before the start of the BTO grid.',1,1)

        RA = sqrt(dot_product(cgto_shell_2%RA,cgto_shell_2%RA))      !radial distance from CMS to the GTO center
        if (RA .ge. cms_bto%A) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','The CGTO center must lie radially before the start of the BTO grid.',2,1)

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs
!
!------ Determine the radial grid appropriate for integrations involving this CGTO shell. The resulting quadrature points and weights are in the arrays r1,w1.
        call radial_grid_CGTO(cgto_shell_1,cms_bto,threshold,x_10,w_10,n_10,n_rng_knot,r1,w1,n1_total_points,R1_max,R1_true)

!------ Do the same for CGTO 2 although we only need the extent of the grid: R2_max.
        call radial_grid_CGTO(cgto_shell_2,cms_bto,threshold,x_10,w_10,n_10,n_rng_knot,r2,w2,n2_total_points,R2_max,R2_true)
!
!------ The end point for the radial integration of the Y function of the second electron with the radial density of the first electron.
!       This point coincides, by construction, with one of the B-spline knots.
        R = min(R1_max,R2_max)
!
!------ Determine index of the last B-spline for which the integrals will be computed. This is done including all B-splines which start before the significance radius R.
        !The last bspline of the first electron is given by the extent of first CGTO and the last bspline of the second electron is determined by the extent of the second CGTO 
        !although this is an overestimation of the actual integrals most of which are actually smaller.
        last_bspline_1 = last_bspline_before_R(R1_max,cms_bto)
        print *,'last_bspline_1',last_bspline_1
        if (last_bspline_1 .eq. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Error either in the routine radial_grid_CGTO or in B-spline grid.',3,1)
!
        last_bspline_2 = last_bspline_before_R(R2_max,cms_bto)
        print *,'last_bspline_2',last_bspline_2
        if (last_bspline_2 .eq. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Error either in the routine radial_grid_CGTO or in B-spline grid.',4,1)
!
!------ Construct the dense grid r2,w2 needed to calculate the Y function. This grid is constructed by constructing the G-L grid of order 7 in between each pair of r1 quadrature points.
!       If there is a knot lying in between the two r1 quadrature points then the interval is divided in two: [r1(i),knot] and [knot,r1(i+1)]. This ensures that when calculating the Y function
!       corresponding to the radial B-splines that the intervals of the integrand close to the end points of the B-splines are always accurately integrated.
        call construct_dense_grid_from_sparse_grid(cms_bto,R,x_10,w_10,n_10,n_rng_knot,r1,n1_total_points,x_7,w_7,n_7,n2_total_points,r2,w2)
!
!------ Precalculate all coupling coefficients needed for the calculations.
        max_l_pw = cms_bto%l+max_l_legendre
        max_l_aux = max_l_pw+max(cgto_shell_1%l,cgto_shell_2%l)

        !The two routines below must be called before any OpenMP region to ensure thread safety of resh, solh.
        call calc_resh_coefficients(max_l_aux)
        call calc_solh_coefficients(max_l_aux)

        !Precalculate all couplings needed
        print *,'gaunt'
        call cpl%prec_cgaunt(max_l_aux)
!
!------ Calculate the projections of the CGTO shell 1 on the real spherical harmonics (i.e. PW-expansion) for all radial points of the grid 1.
        print *,'pw 1'
        call omp_calculate_CGTO_pw_coefficients_analytic(max_l_pw,cgto_shell_1,r1,cgto_pw_shell_1)
        print *,'done'
        print *,'step 1'
        !Multiply in the Gaunt coefficients for L,M values of the Legendre expansion.
        call omp_mult_gaunt_into_angular_integrals(cgto_pw_shell_1,cms_bto%l,max_l_legendre,gaunt_cgto_pw_shell_1)
        deallocate(cgto_pw_shell_1)
!
!------ Calculate the projections of the CGTO shell 2 on the real spherical harmonics (i.e. PW-expansion) for all radial points of the grid 2. This is the dense grid.
        print *,'pw 2'
        call omp_calculate_CGTO_pw_coefficients_analytic(max_l_pw,cgto_shell_2,r2,cgto_pw_shell_2)
        print *,'done',size(r2)
!
!------ Determine mapping of the endpoints of each B-spline with indices in the r1 and r2 arrays.
        call map_knots_to_grid(cms_bto%knots,cms_bto%order,last_bspline_1,r1,bspline_start_end_r1)
        call map_knots_to_grid(cms_bto%knots,cms_bto%order,last_bspline_2,r2,bspline_start_end_r2)
!
!------ Precalculate the radial B-splines on the r1 quadrature grid and multiply them by their norms and quadrature weights w1
        call eval_bspline_on_grid_times_weights(cms_bto,last_bspline_1,bspline_start_end_r1,r1,w1,B_w_r1)
!
!------ Multiply w2 by Jacobian*1/r = r so that the resulting values in B_w_r2 will be: B-spline(r2(i))*r2(i)*w2(i)
        do i=1,n2_total_points
           w2(i) = w2(i)*r2(i)
        enddo !i
!
!------ Precalculate the B-splines on the r2 quadrature grid and multiply them by their norms and quadrature weights w2
        call eval_bspline_on_grid_times_weights(cms_bto,last_bspline_2,bspline_start_end_r2,r2,w2,B_w_r2)
        deallocate(w2)
!
!------ Precalculate the radial terms from the Legendre expansion (4*pi/(2*l+1)*r<^l/r>^(l+1)) on the r1 X r2 (2D) grid.
        print *,'rad legendre'
        call omp_calc_legendre_radial(max_l_legendre,r1,r2,legendre_rad)
        print *,'done'
!
!------ Big loop over angular numbers of the B-splines of the second electron
        n_cgto_1_m = 2*cgto_shell_1%l+1
        n_cgto_2_m = 2*cgto_shell_2%l+1
!
        allocate(R_lm(last_bspline_2,n_cgto_2_m,last_bspline_1,n_cgto_1_m),epstab(52),res3la(3),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Memory allocation 5 failed.',err,1)
        prec_goal = F1MACH(4,cfp_dummy)*100
        overflow = F1MACH(2,cfp_dummy)
        block_total = n_cgto_1_m*last_bspline_1*n_cgto_2_m*last_bspline_2 !number of integrals in the block corresponding to the BTO angular combinations (l1,m1),(l2,m2).
        stride = (2*cache_line_size)/(bit_size(i)/8)
        print *,'prec_goal',prec_goal
        !$OMP PARALLEL DEFAULT(NONE) PRIVATE(n_threads,iam,err) SHARED(thread_cnt,stride)
        !$OMP SINGLE
        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()
        !The array thread_cnt will be used below to store the number of integrals, within each block, converged by each thread.
        allocate(thread_cnt(0:max(n_threads-1,1)*stride),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Memory allocation 6 failed.',err,1)
        !$OMP END SINGLE
        !$OMP END PARALLEL
!
        do l2=0,cms_bto%l
           do m2=-l2,l2
              print *,'Y',l2,m2
              call omp_calculate_Y_BG(cms_bto,last_bspline_2,l2,m2,cgto_pw_shell_2,bspline_start_end_r2,max_l_legendre,n1_total_points,n2_total_points,n_cgto_2_m,legendre_rad,B_w_r2,Y_mp,Y,r1)
              !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l1,m1,l1_m1,to_converge,m_1_ind,ind_1,i,m_2_ind,ind_2,l,m,lm,n,n_wynn,res,abserr,epstab,res3la,nres,n_threads,iam,relative_precision,&
              !$OMP & my_element,err,extrapolated,integral,abserr_lowest,val) &
              !$OMP & SHARED(cms_bto,R_lm,n_cgto_1_m,last_bspline_1,last_bspline_2,n_cgto_2_m,bspline_start_end_r1,B_w_r1,&
              !$OMP & gaunt_cgto_pw_shell_1,Y,overflow,prec_goal,block_total,stdout,l2,m2,stride,thread_cnt,cgto_shell_1,cgto_shell_2)
              n_threads = omp_get_num_threads()
              iam = omp_get_thread_num()
              my_element = stride*iam !there is a separation of lenght 'stride' elements between the values set by each element: this prevents false sharing between threads
              if (size(thread_cnt) < max(n_threads-1,1)*stride+1) call xermsg('bto_gto_integrals_mod','BTO_GTO_two_el_BG_BG','Number of threads used for allocation of thread_cnt &
                 &is smaller than the number of threads used during computation.',6,1)
              !Loop over angular numbers of the B-splines of the first electron
              do l1=0,cms_bto%l
                 do m1=-l1,l1
                    l1_m1 = l1*l1+l1+m1+1

                    !$OMP SINGLE
                    R_lm = 0.0_cfp
                    thread_cnt = 0
                    !$OMP END SINGLE
                    !m-values of the CGTO of the first electron
                    do m_1_ind=1,n_cgto_1_m
                       !radial B-splines of the first electron
                       do ind_1=1,last_bspline_1
                          i = ind_1 + last_bspline_1*(m_1_ind-1)
                          if (mod(i,n_threads) .ne. iam) cycle !work redistribution
                          do m_2_ind=1,n_cgto_2_m
                             do ind_2=1,last_bspline_2
                                !Prepare the epsilon table
                                epstab = 0.0_cfp; res3la = 0.0_cfp; nres = 0
                                !Loop over l,m-angular numbers in the Legendre expansion
                                !write(111,'(2i4,2i4,4i3)') l2,m2, l1,m1, ind_2,m_2_ind,ind_1,m_1_ind
                                integral = 0.0_cfp !This value is the integral as obtained by straight summation of the Legendre terms
                                do l=0,min(max_l_legendre,49)

                                   do m=-l,l
                                      lm = l*l+l+m+1
                                      !Integrate over r1 the product of the Y function of the second electron with the radial density of the first electron.
                                      do i=bspline_start_end_r1(1,ind_1),bspline_start_end_r1(2,ind_1)
                                         integral = integral + B_w_r1(i,ind_1)*gaunt_cgto_pw_shell_1(i,m_1_ind,lm,l1_m1)*Y(i,lm,ind_2,m_2_ind)
                                      enddo !i
                                   enddo !m

                                   !write(111,'(i,e25.15)') l,integral

                                   !Wynn's extrapolation: we need at least 3 points to start the algorithm
                                   n = l+1
                                   if (n .ge. 3) then
                                      !The value of n_wynn can get changed by dqelg if numerical difficulties are encountered:
                                      !typically this will happen for very small (negligible) values of the final integrals.
                                      n_wynn = n
                                      epstab(n) = integral
                                      call dqelg(n_wynn, epstab(1:52), res, abserr, res3la(1:3), nres)
                                      if (res .eq. overflow) then !The extrapolated value overflows - we cannot improve on the result by using more terms so terminate the expansion.
                                         thread_cnt(my_element) = thread_cnt(my_element)+1
                                         exit !converged
                                      elseif (abserr < abserr_lowest) then !Save the extrapolated value as the estimate of the final integral only if it is a more precise estimate.
                                         R_lm(ind_2,m_2_ind,ind_1,m_1_ind) = res
                                         abserr_lowest = abserr
                                         !Estimate relative error in the extrapolated integral using the previous value res3la(2).
                                         relative_precision = abs((res-res3la(2))/res)
                                         !write(111,'("extrapolation",2e25.15)') R_lm(ind_2,m_2_ind,ind_1,m_1_ind), abserr
                                      endif
                                   elseif (n .eq. 2) then !First estimate of the final integrals
                                      R_lm(ind_2,m_2_ind,ind_1,m_1_ind) = integral
                                      abserr_lowest = overflow
                                      relative_precision = 1.0_cfp
                                   endif

                                enddo !l
                                !Ignore convergence of the integrals which are too small.
                                if (relative_precision .le. prec_goal .or. abs(R_lm(ind_2,m_2_ind,ind_1,m_1_ind)) < epsabs) thread_cnt(my_element) = thread_cnt(my_element)+1
                             enddo !ind_2
                          enddo !m_2_ind
                       enddo !ind_1
                    enddo !m_1_ind
                    !$OMP BARRIER

                    !$OMP SINGLE
                    do m_1_ind=1,n_cgto_1_m
                       do ind_1=1,last_bspline_1
                          do m_2_ind=1,n_cgto_2_m
                             do ind_2=1,last_bspline_2
                                write(stdout,'(2i4,2i4,4i3,e25.15)') l2,m2, l1,m1, ind_2,m_2_ind,ind_1,m_1_ind, R_lm(ind_2,m_2_ind,ind_1,m_1_ind)
                             enddo !ind_2
                          enddo !m_2_ind
                       enddo !ind_1
                    enddo !m_1_ind
                    to_converge = block_total-sum(thread_cnt)
                    write(stdout,'("not converged",2i)') to_converge,block_total
                    !$OMP END SINGLE

                 enddo !m1
              enddo !l1
              !$OMP END PARALLEL
           enddo !m2
        enddo !l2

  end subroutine BTO_GTO_two_el_BG_BG

  !> Calculates the part of the Legendre expansion which depends on the B-spline and the radial part of the Legendre expansion. This is calculated on the radial grid r2. The quadrature weights are 
  !> multiplied in as well. These terms are needed to compute the Y function of the second electron.
  subroutine omp_calc_legendre_radial(max_l_legendre,r1,r2,legendre_rad)
     use phys_const, only: fourpi
     implicit none
     !INPUT:
     integer, intent(in) :: max_l_legendre
     real(kind=cfp), allocatable :: r1(:), r2(:)
     !OUTPUT:
     real(kind=cfp), allocatable :: legendre_rad(:,:,:)

     integer :: n2_total_points, n1_total_points, err, i, j, l
     real(kind=cfp) :: fac
     real(kind=cfp), allocatable :: ratio(:,:,:)

       n2_total_points = size(r2)
       n1_total_points = size(r1)

       if (allocated(legendre_rad)) deallocate(legendre_rad)

       allocate(legendre_rad(n2_total_points,n1_total_points,0:max(1,max_l_legendre)),ratio(2,n2_total_points,n1_total_points),stat=err)
       if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calc_legendre_radial','Memory allocation failed.',err,1)

       !Evaluate the ratios of r1 to r2 and vice versa that enter the Leg. expansion
       do j=1,n1_total_points
          do i=1,n2_total_points
             ratio(1,i,j) = r2(i)/r1(j)
             ratio(2,i,j) = r1(j)/r2(i)
          enddo !i
       enddo !j

       !Compute the radial part of the Leg. expansion depending on r< and r> as it enters the expression for the Y function.
       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l,fac,i,j) SHARED(max_l_legendre,n1_total_points,n2_total_points,r1,r2,ratio,legendre_rad)
       !$OMP DO
       do l=0,max_l_legendre
          fac = fourpi/(2*l+1.0_cfp)
          do j=1,n1_total_points
             do i=1,n2_total_points
                if (r2(i) .le. r1(j)) then
                   legendre_rad(i,j,l) = fac*ratio(1,i,j)**l
                else
                   legendre_rad(i,j,l) = fac*ratio(2,i,j)**(l+1)
                endif
             enddo !i
          enddo !j
       enddo !l
       !$OMP END DO
       !$OMP END PARALLEL

       deallocate(ratio)

  end subroutine omp_calc_legendre_radial

  !> Calculates the Y function on the radial grid r1 for a given radial B-spline and all its angular parts and all angular parts of the CGTO. The grid at which the CGTO pw expansion has been evaluated
  !> is given in r2.
  subroutine omp_calculate_Y_BG(cms_bto,last_bspline,lb,mb,cgto_pw_expansion,bspline_start_end_r2,max_l_legendre,n1_total_points,n2_total_points,n_cgto_m,legendre_rad,B_w_r2,Y_mp,Y,r1)
     use omp_lib
     implicit none
     !INPUT:
     integer, intent(in) :: lb, mb, max_l_legendre, n1_total_points, n2_total_points, n_cgto_m, last_bspline
     real(kind=cfp), allocatable :: cgto_pw_expansion(:,:,:), legendre_rad(:,:,:), B_w_r2(:,:), r1(:)
     type(bto_data), intent(inout) :: cms_bto
     integer, allocatable :: bspline_start_end_r2(:,:), r2_eq_r1(:)
     real(kind=cfp), allocatable :: Y_mp(:,:,:,:)
     !OUTPUT:
     real(kind=cfp), allocatable :: Y(:,:,:,:)

     integer :: l,m,lp,mp,lp_mp,m_ind,mp_ind,lm,err,ind,i,j,n_threads,iam,lm_ind
     real(kind=cfp) :: coupling, start_t, end_t, integral
     real(kind=cfp), allocatable :: Y_tmp(:,:,:), cgto_lm(:,:,:)

        if (size(Y_mp,1) .ne. n1_total_points .or. size(Y_mp,2) .ne. last_bspline .or. size(Y_mp,3) .ne. n_cgto_m .or. size(Y_mp,3) .ne. 2*(max_l_legendre+lb)+1) then
           if (allocated(Y_mp)) deallocate(Y_mp)
           allocate(Y_mp(n1_total_points,last_bspline,n_cgto_m,2*(max_l_legendre+lb)+1),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_Y','Memory allocation 1 failed.',err,1)
        endif

        if (size(Y,1) .ne. n1_total_points .or. size(Y,2) .ne. (max_l_legendre+1)**2 .or. size(Y,3) .ne. last_bspline .or. size(Y,4) .ne. n_cgto_m) then
           if (allocated(Y)) deallocate(Y)
           allocate(Y(n1_total_points,(max_l_legendre+1)**2,last_bspline,n_cgto_m),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_Y','Memory allocation 2 failed.',err,1)
        endif

        allocate(cgto_lm(n2_total_points,n_cgto_m,2*max_l_legendre+1),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_Y','Memory allocation 3 failed.',err,1)

        start_t = omp_get_wtime()

        Y = 0.0_cfp
        !$OMP PARALLEL DEFAULT(NONE) PRIVATE(iam,n_threads,l,lp,mp,lp_mp,mp_ind,m_ind,ind,i,j,m,lm,coupling,Y_tmp,err,integral,lm_ind) &
        !$OMP & SHARED(cpl,max_l_legendre,lb,mb,n_cgto_m,last_bspline,n1_total_points,bspline_start_end_r2,Y_mp,cgto_pw_expansion,legendre_rad,B_w_r2,Y,stdout,cgto_lm,n2_total_points)
        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()
        do l=0,max_l_legendre

           lm = l*l
           !$OMP SINGLE
           cgto_lm = 0.0_cfp
           !$OMP END SINGLE

           do lp=abs(l-lb),l+lb
              if (mod(l+lb+lp,2) .ne. 0) cycle
              lp_mp = lp*lp
              !Precalculate all terms that don't depend on the m-value
              do mp=-lp,lp
                 mp_ind = lp_mp+lp+mp+1
                 do m=-l,l
                    lm_ind = l+m+1
                    if (mod(lm_ind,n_threads) .ne. iam) cycle !work redistribution
                    coupling = cpl%rgaunt(l,lp,lb,m,mp,mb)
                    if (coupling .eq. 0.0_cfp) cycle
                    do m_ind=1,n_cgto_m
                       do j=1,n2_total_points
                          cgto_lm(j,m_ind,lm_ind) = cgto_lm(j,m_ind,lm_ind) + coupling*cgto_pw_expansion(j,m_ind,mp_ind)
                       enddo !j
                    enddo !m_ind
                 enddo !m
              enddo !mp
           enddo !lp
           !$OMP BARRIER

           do m_ind=1,n_cgto_m
              do ind=1,last_bspline
                 if (mod(ind+last_bspline*(m_ind-1),n_threads) .ne. iam) cycle !work redistribution
                 do m=-l,l 
                    lm_ind = l+m+1
                    lp_mp = lm+lm_ind 
                    do i=1,n1_total_points
                       integral = 0.0_cfp
                       do j=bspline_start_end_r2(1,ind),bspline_start_end_r2(2,ind)
                          integral = integral + cgto_lm(j,m_ind,lm_ind)*legendre_rad(j,i,l)*B_w_r2(j,ind)
                       enddo !j
                       Y(i,lp_mp,ind,m_ind) = integral
                    enddo !i
                 enddo !m
              enddo !ind
           enddo !m_ind
           !$OMP BARRIER

        enddo !l
        !$OMP END PARALLEL
        end_t = omp_get_wtime()
        print *,'took',end_t-start_t

!        do lp_mp=1,(max_l_legendre+lb+1)**2
!           do m_ind=1,n_cgto_m
!              do j=1,n2_total_points
!                 write(stdout,'(3i,e25.15)') j,m_ind,lp_mp, cgto_pw_expansion(j,m_ind,lp_mp)
!              enddo
!           enddo !m_ind
!        enddo !lp_mp

!        write(stdout,'("Y=")')
!        do lm=1,(max_l_legendre+1)**2
!           do m_ind=1,n_cgto_m
!              do ind=1,last_bspline
!                 do i=1,n1_total_points
!                    write(stdout,'(e25.15,3i,e25.15)') r1(i),ind,m_ind,lm,Y(i,ind,m_ind,lm)
!                 enddo !i
!              enddo !ind
!           enddo !m_ind
!        enddo !lm

  end subroutine omp_calculate_Y_BG

  !> Multiplies the angular integrals by the Gaunt coefficients and sums over the pp indices.
  subroutine omp_mult_gaunt_into_angular_integrals(angular_integrals,l_max,lp_max,gaunt_angular_integrals)
     use omp_lib
     implicit none
     real(kind=cfp), allocatable :: angular_integrals(:,:,:), gaunt_angular_integrals(:,:,:,:)
     integer, intent(in) :: l_max, lp_max

     integer :: Mg, l,m, lp,mp, lpp,mpp, m_ind, n_cgto_m, lm, lp_mp, lpp_mpp, n_total_points, n_Xlm, n_Xlm_p, err, i, n_threads, iam
     real(kind=cfp) :: start_t, end_t, coupling

        n_total_points = size(angular_integrals,1)
        n_cgto_m = size(angular_integrals,2)
        n_Xlm = (l_max+1)**2
        n_Xlm_p = (lp_max+1)**2

        if (allocated(gaunt_angular_integrals)) deallocate(gaunt_angular_integrals)
        allocate(gaunt_angular_integrals(n_total_points,n_cgto_m,n_Xlm_p,n_Xlm),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_mult_gaunt_into_angular_integrals','Memory allocation failed.',err,1)

        gaunt_angular_integrals = 0.0_cfp
        start_t = omp_get_wtime()
        !$OMP PARALLEL DEFAULT(NONE) SHARED(l_max,lp_max,cpl,n_cgto_m,n_total_points,gaunt_angular_integrals,angular_integrals) &
        !$OMP & PRIVATE(l,m,lm,lp,mp,lp_mp,lpp,mpp,lpp_mpp,coupling,m_ind,i,n_threads,iam)
        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()
        do l=0,l_max
           do m=-l,l
              lm = l*l+l+m+1
              if (mod(lm,n_threads) .ne. iam) cycle !work distribution
              if (iam .eq. 0) print *,l,m
              do lp=0,lp_max
                 do mp=-lp,lp
                    lp_mp = lp*lp+lp+mp+1
                    do lpp=abs(l-lp),l+lp
                       do mpp=-lpp,lpp
                          lpp_mpp = lpp*lpp+lpp+mpp+1
                          coupling = cpl%rgaunt(lp,l,lpp,mp,m,mpp)
                          if (coupling .ne. 0.0_cfp) then
                             do m_ind=1,n_cgto_m
                                do i=1,n_total_points
                                   gaunt_angular_integrals(i,m_ind,lp_mp,lm) = gaunt_angular_integrals(i,m_ind,lp_mp,lm) + coupling*angular_integrals(i,m_ind,lpp_mpp)
                                enddo !i
                             enddo !m_ind
                          endif
                       enddo !mpp
                    enddo !lpp
                 enddo !mp
              enddo !lp
           enddo !m
        enddo !l
        !$OMP END PARALLEL
        end_t = omp_get_wtime()
        print *,'took',end_t-start_t

  end subroutine omp_mult_gaunt_into_angular_integrals

  subroutine sum_mp_legendre(gaunt_angular_integrals,r_points,weights,Lg,lp_max,l_max,nuclei,nuclei_RA,lp_integrals)
     use phys_const, only: fourpi
     use common_obj, only: nucleus_type
     implicit none
     real(kind=cfp), allocatable :: gaunt_angular_integrals(:,:,:,:), lp_integrals(:,:,:,:), nuclei_RA(:), r_points(:), weights(:)
     integer, intent(in) :: Lg, l_max, lp_max
     type(nucleus_type), intent(in) :: nuclei(:)

     integer :: l,m,lp,mp,nuc,i,m_ind,n_cgto_m,n_total_points,lm,lp_mp,n_Xlm,n_nuc,err
     real(kind=cfp) :: fac
     real(kind=cfp), allocatable :: rad_legendre(:,:), rad_sum(:)

        n_total_points = size(gaunt_angular_integrals,1)
        n_cgto_m = 2*Lg+1
        n_Xlm = (l_max+1)**2
        n_nuc = size(nuclei)

        if (allocated(lp_integrals)) deallocate(lp_integrals)
        allocate(lp_integrals(n_total_points,n_cgto_m,0:max(lp_max,1),n_Xlm),rad_legendre(n_total_points,n_nuc),rad_sum(n_total_points),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','sum_mp_legendre','Memory allocation failed.',err,1)

        !Real spherical harmonics for the nuclei: result in the module array Xlm_nuclei.
        print *,'nuc'
        call precalculate_Xlm_for_nuclei(nuclei,lp_max)

        lp_integrals = 0.0_cfp
        do l=0,l_max
           print *,l
           do m=-l,l
              lm = l*l+l+m+1
              do lp=0,lp_max

                 !Precalculate the radial terms from the Legendre expansion multiplied by the quadrature weights and nuclear charges and Jacobian of the radial integral and 1/r from BTO
                 fac = fourpi/(2*lp+1.0_cfp)
                 do nuc=1,n_nuc
                    do i=1,n_total_points
                       rad_legendre(i,nuc) = -weights(i)*fac*nuclei(nuc)%charge*(nuclei_RA(nuc)/r_points(i))**lp
                    enddo !i
                 enddo !nuc

                 do mp=-lp,lp
                    lp_mp = lp*lp+lp+mp+1

                    rad_sum = 0.0_cfp
                    do nuc=1,n_nuc
                       fac = Xlm_nuclei(lp_mp+n_Xlm_nuclei*(nuc-1))
                       do i=1,n_total_points
                          rad_sum(i) = rad_sum(i) + fac*rad_legendre(i,nuc)
                       enddo !i
                    enddo !nuc

                    do m_ind=1,n_cgto_m
                       do i=1,n_total_points
                          lp_integrals(i,m_ind,lp,lm) = lp_integrals(i,m_ind,lp,lm) + rad_sum(i)*gaunt_angular_integrals(i,m_ind,lp_mp,lm)
                       enddo !i
                    enddo !m_ind

                 enddo !mp
              enddo !lp
           enddo !m
        enddo !l

  end subroutine sum_mp_legendre

  subroutine integrate_legendre(lp_integrals,bto_norms,bto_val_0,bto_knots,bto_order,bto_n,Lg,lp_max,l_max,integrals)
     implicit none
     real(kind=cfp), allocatable :: lp_integrals(:,:,:,:), bto_norms(:), bto_val_0(:,:), bto_knots(:)
     real(kind=cfp), intent(inout) :: integrals(:,:)
     integer, intent(in) :: bto_order,bto_n,Lg,lp_max,l_max

     logical :: started
     real(kind=cfp) :: res, abserr, relative_precision, prec_goal, overflow
     real(kind=cfp), allocatable :: epstab(:,:,:), res3la(:,:,:), bto_vals(:,:)
     integer, allocatable :: nres(:,:), start_end(:,:)
     logical, allocatable :: nai_converged(:,:)
     integer :: n_Xlm_bto,ind,i,l,m,lm,lp,mp,n,m_ind,base,n_cgto_m,err,n_total_points,to_converge,n_wynn,lp_limit
     
        n_cgto_m = 2*Lg+1
        n_Xlm_bto = (l_max+1)**2
        n_total_points = size(lp_integrals,1)
        prec_goal = F1MACH(4,cfp_dummy)*10
        overflow = F1MACH(2,cfp_dummy)

        allocate(epstab(52,bto_n,n_cgto_m),res3la(3,bto_n,n_cgto_m),nres(bto_n,n_cgto_m),nai_converged(bto_n,n_cgto_m),start_end(2,bto_n),bto_vals(n_total_points,bto_n),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','integrate_legendre','Memory allocation failed.',err,1)

        lp_limit = min(lp_max,49) !49=52-3 as required by the extrapolation routine dqelg
        if (lp_limit < 3) call xermsg('bto_gto_integrals_mod','integrate_legendre','lp_max must be at least 3 for the extrapolation to work.',1,1)

        !Determine the range of indices of points for which each B-spline is non-zero and store the normalized values of the B-splines in bto_vals.
        !The algorithm below assumes that the functions BTO are nodeless as is indeed the case.
        bto_vals = 0.0_cfp
        do ind=1,bto_n
           started = .false.
           start_end(1,ind) = 1
           start_end(2,ind) = n_total_points
           do i=1,n_total_points
              if (bto_val_0(ind,i) .ne. 0.0_cfp) then
                 bto_vals(i,ind) = bto_val_0(ind,i)*bto_norms(ind)
                 if (.not.(started)) then
                    start_end(1,ind) = i
                    started = .true.
                 endif
              else
                 if (started) then
                    start_end(2,ind) = i
                    exit
                 endif
              endif
           enddo !i
        enddo !ind

        integrals = 0.0_cfp
        do l=0,l_max
           print *,l,n_cgto_m*bto_n
           do m=-l,l
              lm=l*l+l+m+1

              epstab = 0.0_cfp; res3la = 0.0_cfp; nres = 0
              to_converge = n_cgto_m*bto_n
              nai_converged = .false.
              do lp=0,lp_limit

                 n = lp+1
                 if (lp > 0) then
                    epstab(n,1:bto_n,1:n_cgto_m) = epstab(n-1,1:bto_n,1:n_cgto_m)
                 endif

                 do m_ind=1,n_cgto_m
                    do ind=1,bto_n
                       if (nai_converged(ind,m_ind)) cycle
                       base = n_Xlm_bto*(ind-1)
                       do i=start_end(1,ind),start_end(2,ind) !1,n_total_points
                          epstab(n,ind,m_ind) = epstab(n,ind,m_ind) + bto_vals(i,ind)*lp_integrals(i,m_ind,lp,lm)
                          !integrals(lm+base,m_ind) = integrals(lm+base,m_ind) + bto_val_0(ind,i)*bto_norms(ind)*lp_integrals(i,m_ind,lp,lm) !=direct summation of Leg. terms with no extrapolation
                       enddo !i
                    enddo !ind
                 enddo !m_ind

                 !Wynn's extrapolation: we need at least 3 points for each integral to start the algorithm
                 !todo correct this so that only the extrapolated values with lower abserr are saved into 'integrals'.
                 if (n .ge. 3) then
                    do m_ind=1,n_cgto_m
                       do ind=1,bto_n
                          if (nai_converged(ind,m_ind)) cycle
                          base = n_Xlm_bto*(ind-1)
                          !The value of n_wynn can get changed by dqelg if numerical difficulties are encountered: typically this will happen for very small (negligible) values of the final integrals.
                          n_wynn = n
                          call dqelg(n_wynn, epstab(1:52,ind,m_ind), res, abserr, res3la(1:3,ind,m_ind), nres(ind,m_ind))
                          if (res .eq. overflow) then
                             to_converge = to_converge-1
                             nai_converged(ind,m_ind) = .true.
                          else
                             relative_precision = abs((res-integrals(lm+base,m_ind))/res)
                             integrals(lm+base,m_ind) = res
                             if (relative_precision .le. prec_goal) then
                                to_converge = to_converge-1
                                nai_converged(ind,m_ind) = .true.
                             endif
                          endif
                       enddo !ind
                    enddo !m_ind
                 elseif (n .eq. 2) then !First estimate of the final integrals
                    do m_ind=1,n_cgto_m
                       do ind=1,bto_n
                          base = n_Xlm_bto*(ind-1)
                          integrals(lm+base,m_ind) = epstab(n,ind,m_ind)
                       enddo !ind
                    enddo !m_ind
                 endif

!                 write(100,'(3i)') lm,lp,to_converge
                 if (to_converge .eq. 0) exit !All integrals have been converged with respect to summation over lp

              enddo !lp
              print *,'not converged',to_converge,n
!              if (to_converge > 0) then
!                 do m_ind=1,n_cgto_m
!                    do ind=1,bto_n
!                       if (nai_converged(ind,m_ind)) cycle
!                       base = n_Xlm_bto*(ind-1)
!                       res = epstab(n,ind,m_ind)
!                       relative_precision = abs((res-integrals(lm+base,m_ind))/res)
!                       write(100,'("convergence",3i,2e25.15)') ind,m_ind,lm,integrals(lm+base,m_ind),relative_precision
!!                       do n=1,lp_limit+1
!!                          write(100,'("wynn",i,e25.15)') n,epstab(n,ind,m_ind)
!!                       enddo
!                    enddo !ind
!                 enddo !m_ind
!              endif

           enddo !m
        enddo !l

  end subroutine integrate_legendre

  !> Adds the Bloch terms to the kinetic energy integrals.
  subroutine add_cgto_bto_bloch_elements(bto_norms,cms_bto,cgto,integrals)
     use phys_const, only: fourpi
     use bspline_base, only: bvalu
     implicit none
     type(bto_data), intent(inout) :: cms_bto
     type(cgto_data), intent(in) :: cgto
     real(kind=cfp), intent(inout) :: integrals(:,:), bto_norms(:)

     real(kind=cfp) :: rmat_radius, bto_val_1, cgto_val, a0_square, fac, cgto_val_m, bloch
     integer :: base, n_Xlm_bto, m, lm, i, j

        rmat_radius = cms_bto%B

        n_Xlm_bto = (cms_bto%l+1)**2

        a0_square = rmat_radius*rmat_radius
        fac = cgto%norm*sqrt(fourpi/(2*cgto%l+1.0_cfp))*rmat_radius**(cgto%l+1) !1/r*CGTO*Jac at r=a0

        cms_bto%bcoef = 0.0_cfp
        do i=1,cms_bto%n
           if (cms_bto%knots(i+cms_bto%order) .eq. rmat_radius) then !this BTO ends at r=a0, i.e. there might be a surface term to get rid of using the Bloch op.
              base = (i-1)*n_Xlm_bto

              !first derivative of the B-spline at the end point of the radial grid: only the B-splines at the end point of the interval may give rise to Bloch terms
              cms_bto%bcoef(i) = 1.0_cfp
              bto_val_1 = bto_norms(i)*bvalu(cms_bto%knots,cms_bto%bcoef,cms_bto%n,cms_bto%order,1,rmat_radius,cms_bto%inbv,cms_bto%work)
              cms_bto%bcoef(i) = 0.0_cfp !clean-up for the next

              if (bto_val_1 .eq. 0.0_cfp) cycle

              !Radial part of the CGTO at r=a0
              cgto_val = 0.0_cfp
              do j=1,cgto%n_prim
                 cgto_val = cgto_val + cgto%ccf(j)*cgto%p_norm(j)*exp(-cgto%alp(j)*a0_square)
              enddo
              cgto_val = cgto_val*fac

              lm = cgto%l*cgto%l+cgto%l+1
              do m=-cgto%l,cgto%l
                 cgto_val_m = (-1)**m*cgto_val

                 bloch = 0.5_cfp*cgto_val_m*bto_val_1 !Bloch term for the R-matrix radius a0
                 print *,'bloch',bloch

                 !BTO (l,m) must be the same as CGTO (l,m) for the Bloch term to be non-zero
                 if (lm+m .le. n_Xlm_bto) integrals(base+lm+m,cgto%l+m+1) = integrals(base+lm+m,cgto%l+m+1) + bloch
              enddo !m

           endif
        enddo !i

  end subroutine add_cgto_bto_bloch_elements

  !> Naive cyclic swap: dimensions 1,2,3 -> 2,3,1. No cache blocking.
  subroutine swap_dim_1_with_3(array_in,array_out)
     implicit none
     real(kind=cfp), intent(inout) :: array_in(:,:,:), array_out(:,:,:)

     integer :: d1,d2,d3, i,j,k

        d1 = size(array_in,1)
        d2 = size(array_in,2)
        d3 = size(array_in,3)

        do k=1,d3
           do j=1,d2
              do i=1,d1
                 array_out(j,k,i) = array_in(i,j,k)
              enddo !i
           enddo !j
        enddo !k

  end subroutine swap_dim_1_with_3

  !> Assumes B .ge. A.
  function r_is_inside_A_B(A,B,r)
     implicit none
     real(kind=cfp), intent(in) :: A, B, r
     logical :: r_is_inside_A_B

        r_is_inside_A_B = .true.
        if (r > B .or. r < A) r_is_inside_A_B = .false.

  end function r_is_inside_A_B

  !> Calculates the partial wave expansion of CGTO up to partial wave L=max_l and on the grid of radial points r_points.
  subroutine omp_calculate_CGTO_pw_coefficients_analytic(max_l,cgto,r_points,angular_integrals)
      use phys_const, only: fourpi, twopi, pi
      use special_functions, only: cfp_besi, cfp_eval_poly_horner
      use omp_lib
      implicit none
      integer, intent(in) :: max_l
      type(cgto_data), intent(in) :: cgto
      real(kind=cfp), allocatable :: angular_integrals(:,:,:)
      real(kind=cfp), allocatable :: r_points(:)

      integer, parameter :: kode = 2
      real(kind=cfp), parameter :: half = 0.5_cfp

      integer :: i, j, err, l, m, m_ind, lm, Mg, l_besi, nz, Lg_Mg, max_l_aux, n_threads, iam, n_cgto_m, n_Xlm
      real(kind=cfp) :: r_square, arg, fnu, asym, term, RA, start_t, end_t
      real(kind=cfp), allocatable :: besi_vals(:,:), contraction_besi(:,:), contraction(:), prim_fac(:,:), besi_args(:,:), transl_cfs(:,:,:), c_lambda(:,:,:), Xlm_CGTO_center(:)
      integer :: n_total_points

         if (allocated(angular_integrals)) deallocate(angular_integrals)
         n_total_points = size(r_points)
         n_cgto_m = 2*cgto%l+1
         n_Xlm = (max_l+1)**2
         allocate(angular_integrals(n_total_points,n_cgto_m,n_Xlm),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_CGTO_pw_coefficients_analytic','Memory allocation 1 failed.',err,1)

         RA = sqrt(dot_product(cgto%RA,cgto%RA))

         !Calculate the angular integrals trivially if CGTO is at the CMS.
         if (RA .eq. 0.0_cfp) then

            allocate(contraction(n_total_points),stat=err)
            if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_CGTO_pw_coefficients_analytic','Memory allocation 2 failed.',err,1)

            do i=1,n_total_points
               !Radial part of the CGTO
               r_square = r_points(i)*r_points(i)
               contraction(i) = 0.0_cfp
               do j=1,cgto%n_prim
                  contraction(i) = contraction(i) + cgto%ccf(j)*cgto%p_norm(j)*exp(-cgto%alp(j)*r_square)
               enddo
               contraction(i) = contraction(i)*cgto%norm*sqrt(fourpi/(2*cgto%l+1.0_cfp))*r_points(i)**cgto%l
            enddo !i
   
            !The only projections that are non-zero are those where the real spherical harmonic has the same (l,m) values as the CGTO:
            lm = cgto%l*cgto%l+cgto%l+1
            do m=-cgto%l,cgto%l
               term = (-1)**m
               do i=1,n_total_points
                  angular_integrals(i,cgto%l+m+1,lm+m) = term*contraction(i)
               enddo !i
            enddo
            return

         endif
!
!------- Use single centre expansion of the CGTO to evaluate the partial wave coefficients.
!        Angular integrals for s-type CGTOs can be always evaluated accurately using a single centre expansion; tests show that for higher-L CGTOs accuracy is preserved as well.
         l_besi = max_l+cgto%l+1

         FNU = l_besi-1+half 
         !Compute the smallest value of the Bessel I argument for which the asymptotic expansion can be used
         !todo the asym value is not used anywhere at the moment
         if (cfp .eq. wp) then
            asym = MAX(17.00_cfp,0.550_cfp*FNU*FNU)
         else !quad precision case
            asym = MAX(76.00_cfp,0.550_cfp*FNU*FNU)
         endif 

         !Buffer for exponentially scaled modified Bessel values for each CGTO primitive exponent.
         allocate(besi_vals(1:l_besi,1:cgto%n_prim),contraction_besi(0:l_besi,n_total_points),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_CGTO_pw_coefficients_analytic','Memory allocation 3 failed.',err,1)

         do i=1,n_total_points
            do j=1,cgto%n_prim
               arg = 2.0_cfp*cgto%alp(j)*r_points(i)*RA
               call cfp_besi(arg,half,kode,l_besi,besi_vals(1:l_besi,j),nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=l+1 is needed to get y_{alpha+l}.
            enddo !j
   
            r_square = (r_points(i)-RA)**2
            do l=0,l_besi-1
               !Radial part of the CGTO multiplied by the modified Bessel function
               contraction_besi(l,i) = 0.0_cfp
               do j=1,cgto%n_prim
                  contraction_besi(l,i) = contraction_besi(l,i) + cgto%ccf(j)*cgto%p_norm(j)*exp(-cgto%alp(j)*r_square)*besi_vals(l+1,j)/sqrt(cgto%alp(j)*r_points(i)*RA)
               enddo !j
               contraction_besi(l,i) = contraction_besi(l,i)*cgto%norm*twopi*sqrt(pi)
            enddo !l
         enddo !i

         if (cgto%l .eq. 0) then
            !Real spherical harmonics for the nuclei: result in the module array Xlm_CGTO_center
            print *,'cent'
            call precalculate_Xlm_for_CGTO_center(cgto%RA,max_l,Xlm_CGTO_center)

            m_ind = 1
            do l=0,max_l
               lm = l*l+l+1
               do m=-l,l
                  do i=1,n_total_points
                     angular_integrals(i,m_ind,lm+m) = Xlm_CGTO_center(lm+m)*contraction_besi(l,i)
                  enddo !i
               enddo !m
            enddo !l
         else
!            allocate(prim_fac(cgto%n_prim,n_total_points),besi_args(cgto%n_prim,n_total_points),stat=err)
!            if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_CGTO_pw_coefficients_analytic','Memory allocation 4 failed.',err,1)
!
!            do i=1,n_total_points
!               do j=1,cgto%n_prim
!                  besi_args(j,i) = 2.0_cfp*cgto%alp(j)*r_points(i)*RA
!               enddo !j
!   
!               r_square = (r_points(i)-RA)**2
!               !Radial part of the CGTO multiplied by some factors from the SCE
!               do j=1,cgto%n_prim
!                  prim_fac(j,i) = cgto%norm*twopi*sqrt(pi)*cgto%ccf(j)*cgto%p_norm(j)*exp(-cgto%alp(j)*r_square)/sqrt(cgto%alp(j)*r_points(i)*RA)
!               enddo !j
!            enddo !i

            max_l_aux = max_l+cgto%l
            !Real spherical harmonics for the nuclei: result in the module array Xlm_CGTO_center
            call precalculate_Xlm_for_CGTO_center(cgto%RA,max_l_aux+cgto%l,Xlm_CGTO_center)

            !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
            call precalculate_solh_translation_coeffs(cgto%l,RA,Xlm_CGTO_center,transl_cfs)
!            print *,'prec done'

            start_t = omp_get_wtime()
            !Evaluate the SCE of the CGTO at the radial points
            !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l,m,lm,m_ind,Mg,Lg_Mg,i,n_threads,iam,c_lambda) &
            !$OMP & SHARED(max_l,cgto,n_total_points,l_besi,angular_integrals,r_points,contraction_besi,asym,prim_fac,besi_args,transl_cfs,Xlm_CGTO_center)
            n_threads = omp_get_num_threads()
            iam = omp_get_thread_num()
            do l=0,max_l
               do m=-l,l
                  lm = l*l+l+m+1
                  if (mod(lm,n_threads) .ne. iam) cycle !work distribution
                  !Calculate the coupling coefficients needed to evaluate the projection using SCE of GTO.
                  call calculate_lambda_couplings(cgto%l,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
                  m_ind = 0
                  do Mg=-cgto%l,cgto%l
                     m_ind=m_ind+1
                     Lg_Mg = cgto%l*cgto%l+cgto%l+mg+1
                     do i=1,n_total_points
                        angular_integrals(i,m_ind,lm) = CGTO_pw_coefficient(r_points(i),l,cgto%l,Mg,c_lambda,contraction_besi(0:l_besi,i))
!                        angular_integrals(i,m_ind,lm) = CGTO_pw_coefficient_stable(asym,r_points(i),l,lm,cgto%l,Lg_Mg,c_lambda,contraction_besi(0:l_besi,i),prim_fac,cgto%n_prim,besi_args(1:cgto%n_prim,i))
!                        write(111,'(e25.15,2i,e25.15)') r_points(i),m_ind,lm,angular_integrals(i,m_ind,lm)
                     enddo !i
                  enddo !Mg
               enddo !m
            enddo !l
            !$OMP END PARALLEL
            end_t = omp_get_wtime()
!            print *,'pw took',end_t-start_t
         endif

  end subroutine omp_calculate_CGTO_pw_coefficients_analytic

  function CGTO_pw_coefficient(r,l,Lg,Mg,c_lambda,contraction_besi)
      use special_functions, only: cfp_eval_poly_horner
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
      use special_functions, only: cfp_eval_poly_horner
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
                 lambda_cf(lp) = lambda_cf(lp) + prim_fac(j)*sum_besi_half_asym_lambda_cf(besi_args(j),c_lambda(0:l+lp,lp,lm,Lg_Mg),l+lp,abs(l-lp),l+lp)
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

  !> Calculates the partial wave expansion of a product of two CGTOs up to partial wave L=max_l and on the grid of radial points r_points.
  subroutine omp_calculate_CGTO_pair_pw_coefficients_analytic(max_l,cgto_A,cgto_B,r_points,angular_integrals)
      use phys_const, only: fourpi, twopi, pi
      use special_functions, only: cfp_besi, cfp_eval_poly_horner
      use omp_lib
      implicit none
      integer, intent(in) :: max_l
      type(cgto_data), intent(in) :: cgto_A, cgto_B
      real(kind=cfp), allocatable :: angular_integrals(:,:,:,:)
      real(kind=cfp), allocatable :: r_points(:)

      integer, parameter :: kode = 2
      real(kind=cfp), parameter :: half = 0.5_cfp

      integer :: i, j, err, l, m, lm, l_besi, nz, max_l_aux, n_threads, iam, n_cgto_A_m, n_cgto_B_m, n_Xlm, p, q, ij
      real(kind=cfp) :: r_square, arg, fnu, asym, term, R_A, R_B, start_t, end_t, R_AB(3), R_AB_square, fac
      real(kind=cfp), allocatable :: besi_vals(:,:), contraction_besi(:,:,:), contraction(:), prim_fac(:,:), besi_args(:,:), transl_cfs(:,:,:), c_lambda(:,:,:), prod_alp(:), K_pref(:)
      real(kind=cfp), allocatable :: prod_P(:,:), RP(:), exp_fac(:), prod_contr(:), Xlm_CGTO_A_center(:), Xlm_CGTO_B_center(:), Xlm_product_CGTO_center(:), transl_cfs_AB(:,:,:,:,:), c_pair_lambda(:,:,:,:,:)
      integer :: n_total_points, n_contr_pairs, Mg_A, Mg_B, Mg_A_ind, Mg_B_ind

         if (allocated(angular_integrals)) deallocate(angular_integrals)
         n_total_points = size(r_points)
         n_cgto_A_m = 2*cgto_A%l+1
         n_cgto_B_m = 2*cgto_B%l+1
         n_Xlm = (max_l+1)**2
         allocate(angular_integrals(n_total_points,n_cgto_B_m,n_cgto_A_m,n_Xlm),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_CGTO_pw_coefficients_analytic','Memory allocation 1 failed.',err,1)

         R_A = sqrt(dot_product(cgto_A%RA,cgto_A%RA))
         R_B = sqrt(dot_product(cgto_B%RA,cgto_B%RA))
         !todo if both CGTOs are at CMS then the present method will not work: I need a special case for this

!         !Calculate the angular integrals trivially if CGTO is at the CMS.
!         if (RA .eq. 0.0_cfp) then
!
!            allocate(contraction(n_total_points),stat=err)
!            if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_CGTO_pw_coefficients_analytic','Memory allocation 2 failed.',err,1)
!
!            do i=1,n_total_points
!               !Radial part of the CGTO
!               r_square = r_points(i)*r_points(i)
!               contraction(i) = 0.0_cfp
!               do j=1,cgto%n_prim
!                  contraction(i) = contraction(i) + cgto%ccf(j)*cgto%p_norm(j)*exp(-cgto%alp(j)*r_square)
!               enddo
!               contraction(i) = contraction(i)*cgto%norm*sqrt(fourpi/(2*cgto%l+1.0_cfp))*r_points(i)**cgto%l
!            enddo !i
!   
!            !The only projections that are non-zero are those where the real spherical harmonic has the same (l,m) values as the CGTO:
!            lm = cgto%l*cgto%l+cgto%l+1
!            do m=-cgto%l,cgto%l
!               term = (-1)**m
!               do i=1,n_total_points
!                  angular_integrals(i,cgto%l+m+1,lm+m) = term*contraction(i)
!               enddo !i
!            enddo
!            return
!
!         endif
!
!------- Use single centre expansion of the CGTO to evaluate the partial wave coefficients.
!        Angular integrals for s-type CGTOs can be always evaluated accurately using a single centre expansion; tests show that for higher-L CGTOs accuracy is preserved as well.
         !todo implement special case of both CGTOs sharing the same center -> SCE of a single CGTO.
         l_besi = max_l+cgto_A%l+cgto_B%l+1

         FNU = l_besi-1+half 
         !Compute the smallest value of the Bessel I argument for which the asymptotic expansion can be used
         !todo the asym value is not used anywhere at the moment
         if (cfp .eq. wp) then
            asym = MAX(17.00_cfp,0.550_cfp*FNU*FNU)
         else !quad precision case
            asym = MAX(76.00_cfp,0.550_cfp*FNU*FNU)
         endif 

         n_contr_pairs = cgto_A%n_prim*cgto_B%n_prim

         !Buffer for exponentially scaled modified Bessel values for each CGTO primitive exponent.
         allocate(besi_vals(1:l_besi,1:n_contr_pairs),contraction_besi(n_contr_pairs,0:l_besi,n_total_points),prod_P(3,n_contr_pairs),prod_alp(n_contr_pairs),&
                 &K_pref(n_contr_pairs),RP(n_contr_pairs),exp_fac(n_contr_pairs),prod_contr(n_contr_pairs),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','omp_calculate_CGTO_pw_coefficients_analytic','Memory allocation 3 failed.',err,1)

         R_AB(1:3) = cgto_A%RA(1:3)-cgto_B%RA(1:3)
         R_AB_square = dot_product(R_AB,R_AB)
         fac = cgto_A%norm*cgto_B%norm*twopi*sqrt(pi)
         ij = 0
         do p=1,cgto_A%n_prim
            do q=1,cgto_B%n_prim
               ij = ij + 1
               prod_alp(ij) = cgto_A%alp(p)+cgto_B%alp(q) !exponent of the product GTO
               prod_P(1:3,ij) = (cgto_A%alp(p)*cgto_A%RA(1:3)+cgto_B%alp(q)*cgto_B%RA(1:3))/prod_alp(ij) !center of the product GTO
               K_pref(ij) = exp(-cgto_A%alp(p)*cgto_B%alp(q)/prod_alp(ij)*R_AB_square) !exponential prefactor for the product GTO
               RP(ij) = sqrt(dot_product(prod_P(1:3,ij),prod_P(1:3,ij)))
               prod_contr(ij) = fac*cgto_A%ccf(p)*cgto_A%p_norm(p)*cgto_B%ccf(q)*cgto_B%p_norm(q)
            enddo !q
         enddo !p

         do i=1,n_total_points
            do ij=1,n_contr_pairs
               arg = 2.0_cfp*prod_alp(ij)*r_points(i)*RP(ij)
               call cfp_besi(arg,half,kode,l_besi,besi_vals(1:l_besi,ij),nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=l+1 is needed to get y_{alpha+l}.
            enddo !ij

            ij = 0
            do p=1,cgto_A%n_prim
               do q=1,cgto_B%n_prim
                  ij = ij + 1
                  r_square = (r_points(i)-RP(ij))**2
                  exp_fac(ij) = K_pref(ij)*prod_contr(ij)*exp(-prod_alp(ij)*r_square)/sqrt(prod_alp(ij)*r_points(i)*RP(ij))
               enddo !q
            enddo !p

            do l=0,l_besi-1
               !Radial parts of the product GTOs multiplied by the modified Bessel functions
               ij = 0
               do p=1,cgto_A%n_prim
                  do q=1,cgto_B%n_prim
                     ij = ij + 1
                     contraction_besi(ij,l,i) = exp_fac(ij)*besi_vals(l+1,ij)
                  enddo !q
               enddo !p
            enddo !l
         enddo !i

!         if (cgto%l .eq. 0) then
!            !Real spherical harmonics for the nuclei: result in the module array Xlm_CGTO_center
!            print *,'cent'
!            call precalculate_Xlm_for_CGTO_center(cgto%RA,max_l,Xlm_CGTO_center)
!
!            m_ind = 1
!            do l=0,max_l
!               lm = l*l+l+1
!               do m=-l,l
!                  do i=1,n_total_points
!                     angular_integrals(i,m_ind,lm+m) = Xlm_CGTO_center(lm+m)*contraction_besi(l,i)
!                  enddo !i
!               enddo !m
!            enddo !l
!         else

            !Real spherical harmonics for the CGTO A center: result in the module array Xlm_CGTO_A_center
            call precalculate_Xlm_for_CGTO_center(cgto_A%RA,max_l+cgto_A%l,Xlm_CGTO_A_center)

            !Real spherical harmonics for the CGTO B center: result in the module array Xlm_CGTO_B_center
            call precalculate_Xlm_for_CGTO_center(cgto_B%RA,max_l+cgto_B%l,Xlm_CGTO_B_center)

            max_l_aux = max_l+cgto_A%l+cgto_B%l
            call precalculate_Xlm_for_CGTO_product_center(n_contr_pairs,prod_P,max_l_aux,Xlm_product_CGTO_center)

            !Precalculate the coefficients in the translation formula for the pair of solid harmonics sitting on centers A,B: this requires Xlm_CGTO_center
            call precalculate_pair_solh_translation_coeffs(cgto_A%l,R_A,Xlm_CGTO_A_center,cgto_B%l,R_B,Xlm_CGTO_B_center,transl_cfs_AB)
!            print *,'prec done'

            start_t = omp_get_wtime()
            !Evaluate the SCE of the CGTO at the radial points
            !$OMP PARALLEL DEFAULT(NONE) PRIVATE(l,m,lm,Mg_A,Mg_A_ind,Mg_B,Mg_B_ind,i,n_threads,iam,c_pair_lambda) &
            !$OMP & SHARED(max_l,cgto_A,cgto_B,n_total_points,l_besi,angular_integrals,r_points,contraction_besi,n_contr_pairs,Xlm_product_CGTO_center,transl_cfs_AB)
            n_threads = omp_get_num_threads()
            iam = omp_get_thread_num()
            do l=0,max_l
               do m=-l,l
                  lm = l*l+l+m+1
                  if (mod(lm,n_threads) .ne. iam) cycle !work distribution
                  !Calculate the coupling coefficients needed to evaluate the projection using SCE of GTO.
                  call calculate_pair_lambda_couplings(l,m,cgto_A%l,cgto_B%l,n_contr_pairs,Xlm_product_CGTO_center,transl_cfs_AB,c_pair_lambda)
                  do Mg_A=-cgto_A%l,cgto_A%l
                     Mg_A_ind = Mg_A+cgto_A%l+1
                     do Mg_B=-cgto_B%l,cgto_B%l
                        Mg_B_ind = Mg_B+cgto_B%l+1
                        do i=1,n_total_points
                           angular_integrals(i,Mg_B_ind,Mg_A_ind,lm) = &
                           &CGTO_pair_pw_coefficient(r_points(i),l,cgto_A%l,Mg_A,cgto_B%l,Mg_B,c_pair_lambda,n_contr_pairs,contraction_besi(1:n_contr_pairs,0:l_besi,i))
                           !write(111,'(e25.15,3i,e25.15)') r_points(i),Mg_B,Mg_A,lm,angular_integrals(i,Mg_B_ind,Mg_A_ind,lm)
                        enddo !i
                     enddo !Mg_B
                  enddo !Mg_A
               enddo !m
            enddo !l
            !$OMP END PARALLEL
            end_t = omp_get_wtime()
!            print *,'pw took',end_t-start_t
!         endif

  end subroutine omp_calculate_CGTO_pair_pw_coefficients_analytic

  function CGTO_pair_pw_coefficient(r,l,Lg_A,Mg_A,Lg_B,Mg_B,c_pair_lambda,n_contr_pairs,contraction_besi)
      use special_functions, only: cfp_eval_poly_horner
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
  subroutine precalculate_pair_solh_translation_coeffs(CGTO_A_L,RA_A,Xlm_CGTO_A_center,CGTO_B_L,RA_B,Xlm_CGTO_B_center,transl_cfs_AB)
     use phys_const, only: fourpi
     implicit none
     real(kind=cfp), allocatable :: Xlm_CGTO_A_center(:), Xlm_CGTO_B_center(:), transl_cfs_AB(:,:,:,:,:)
     integer, intent(in) :: CGTO_A_L, CGTO_B_L
     real(kind=cfp), intent(in) :: RA_A, RA_B

     integer :: n_mp, err, CGTO_M, lp, mp, max_lp, lp_min, lp_max, CGTO_A_M, CGTO_A_M_ind, CGTO_B_M, CGTO_B_M_ind, la_p, lb_p, la_p_lb_p, ma_p, mb_p
     real(kind=cfp), allocatable :: transl_cfs_A(:,:,:), transl_cfs_B(:,:,:)
   
        !Translation coefficients for the individual CGTOs
        call precalculate_solh_translation_coeffs(CGTO_A_L,RA_A,Xlm_CGTO_A_center,transl_cfs_A)
        call precalculate_solh_translation_coeffs(CGTO_B_L,RA_B,Xlm_CGTO_B_center,transl_cfs_B)

        max_lp = CGTO_A_L+CGTO_B_L
        n_mp = 2*max_lp+1

        if (allocated(transl_cfs_AB)) deallocate(transl_cfs_AB)
        allocate(transl_cfs_AB(n_mp,0:max(max_lp,1),0:max(max_lp,1),2*CGTO_B_L+1,2*CGTO_A_L+1),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','precalculate_solh_translation_coeffs','Memory allocation failed.',err,1)

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
                                transl_cfs_AB(mp+lp+1,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) = transl_cfs_AB(mp+lp+1,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) + &
                                &transl_cfs_A(ma_p+la_p+1,la_p,CGTO_A_M_ind)*transl_cfs_B(mb_p+lb_p+1,lb_p,CGTO_B_M_ind)*cpl%rgaunt(lp,la_p,lb_p,mp,ma_p,mb_p)
                             enddo !mp
                          enddo !lp
                       enddo !mb_p
                    enddo !ma_p
                 enddo !lb_p
              enddo !la_p
           enddo !CGTO_B_M
        enddo !CGTO_A_M

  end subroutine precalculate_pair_solh_translation_coeffs

  !> \warning Requires precalculated values of the real spherical harmonics at the position of the CGTO nucleus. The coupling coefficients should also be precalculated for performance reasons.
  subroutine precalculate_solh_translation_coeffs(CGTO_L,RA,Xlm_CGTO_center,transl_cfs)
     use phys_const, only: fourpi
     implicit none
     real(kind=cfp), allocatable :: transl_cfs(:,:,:), Xlm_CGTO_center(:)
     integer, intent(in) :: CGTO_L
     real(kind=cfp), intent(in) :: RA

     integer :: n_mp, err, CGTO_M, lp, lpp, mp, mpp, ind
     real(kind=cfp) :: fac, sum_mpp
   
        n_mp = 2*CGTO_L+1

        if (allocated(transl_cfs)) deallocate(transl_cfs)
        allocate(transl_cfs(n_mp,0:max(CGTO_L,1),n_mp),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','precalculate_solh_translation_coeffs','Memory allocation failed.',err,1)

        call cpl%prec_G_cf(CGTO_L)

        fac = sqrt(fourpi/(2*CGTO_L+1.0_cfp))
        transl_cfs = 0.0_cfp
        do CGTO_M=-CGTO_L,CGTO_L
           do lp=0,CGTO_L
              lpp = CGTO_L-lp
              do mp=-lp,lp
                 sum_mpp = 0.0_cfp
                 do mpp=-lpp,lpp
                    ind = lpp*lpp+lpp+mpp+1
                    sum_mpp = sum_mpp + cpl%G_real_cf(CGTO_L,lp,CGTO_M,mp,mpp)*Xlm_CGTO_center(ind)*(-1)**mpp
                 enddo !mpp
                 transl_cfs(mp+lp+1,lp,CGTO_M+CGTO_L+1) = sum_mpp*fac*(-1)**(lpp+mp)*RA**lpp
              enddo !mp
           enddo !lp
        enddo !CGTO_M

  end subroutine precalculate_solh_translation_coeffs

  subroutine calculate_pair_lambda_couplings(l,m,CGTO_A_L,CGTO_B_L,n_contr_pairs,Xlm_product_CGTO_center,transl_cfs_AB,c_pair_lambda)
     implicit none
     real(kind=cfp), allocatable :: transl_cfs_AB(:,:,:,:,:), Xlm_product_CGTO_center(:)
     integer, intent(in) :: CGTO_A_L, CGTO_B_L, l, m, n_contr_pairs
     !OUTPUT:
     real(kind=cfp), allocatable :: c_pair_lambda(:,:,:,:,:)

     real(kind=cfp) :: coupling 
     integer :: CGTO_M, CGTO_M_ind, lm, lp, mp, lambda, mu, base, lambda_max, err, d1, d2, d3, d4, d5, mp_ind, CGTO_A_M, CGTO_B_M, la_p_lb_p, ij, CGTO_A_M_ind, CGTO_B_M_ind

        lambda_max = l+CGTO_A_L+CGTO_B_L
        d1 = n_contr_pairs
        d2 = max(lambda_max,1)
        d3 = max(CGTO_A_L+CGTO_B_L,1)
        d4 = 2*CGTO_B_L+1
        d5 = 2*CGTO_A_L+1
        if (ubound(c_pair_lambda,1) < d1 .or. ubound(c_pair_lambda,2) < d2 .or. ubound(c_pair_lambda,3) < d3 .or. ubound(c_pair_lambda,4) < d4 .or. ubound(c_pair_lambda,5) < d5) then
           if (allocated(c_pair_lambda)) deallocate(c_pair_lambda)
           allocate(c_pair_lambda(1:d1,0:d2,0:d3,1:d4,1:d5),stat=err)
           if (err .ne. 0) then
              print *,d1,d2,d3,d4,d5
              call xermsg('bto_gto_integrals_mod','calculate_pair_lambda_couplings','Memory allocation failed.',err,1)
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
                                c_pair_lambda(ij,lambda,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) = c_pair_lambda(ij,lambda,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind)&
                                & + coupling*transl_cfs_AB(mp_ind,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind)*Xlm_product_CGTO_center(ij+base)
                             enddo !ij
                          enddo !la_p_lb_p
                       enddo !CGTO_B_M
                    enddo !CGTO_A_M
                 enddo !mu
              enddo !mp
           enddo !lambda
        enddo !lp

  end subroutine calculate_pair_lambda_couplings

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
        if (ubound(c_lambda,1) < d1 .or. ubound(c_lambda,2) < d2 .or. ubound(c_lambda,3) < d3) then
           if (allocated(c_lambda)) deallocate(c_lambda)
           allocate(c_lambda(0:d1,0:d2,d3),stat=err)
           if (err .ne. 0) call xermsg('bto_gto_integrals_mod','calculate_lambda_couplings','Memory allocation failed.',err,1)
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
     use phys_const, only: twopi
     implicit none
     real(kind=cfp), intent(in) :: z, couplings(0:l)
     integer, intent(in) :: l, lmin,lmax
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
     use phys_const, only: twopi
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

  subroutine calculate_Tl_matrix(bspline_sol,Tl)
     use quadrature_module, only: cfp_bfqad
     implicit none
     type(bto_data), intent(inout) :: bspline_sol
     real(kind=cfp), allocatable :: Tl(:,:,:)

     type(radial_free_bspline) :: integrand
     integer :: err, l, i, j
     real(kind=cfp) :: r1, r2, quad

        err = bspline_sol%check()
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','calculate_Tl_matrix','Check of the B-spline grid has failed with an error.',err,1)

        integrand%bspline_radial_grid = bspline_sol

        if (allocated(Tl)) deallocate(Tl)
        allocate(Tl(bspline_sol%n,bspline_sol%n,0:max(1,bspline_sol%l)),stat=err)
        if (err .ne. 0) call xermsg('bto_gto_integrals_mod','calculate_Tl_matrix','Memory allocation failed.',err,1)
        Tl = 0.0_cfp

        bspline_sol%bcoef = 0.0_cfp

        do l=0,bspline_sol%l
           integrand%bspline_radial_grid%l = l
        
           do i=1,bspline_sol%n
              integrand%bspline_radial_grid%ind = i
              do j=1,bspline_sol%n
                 !integration limits: the lower limit is determined by the starting point of the B-spline that starts further to the right; the oposite holds for the upper integration limit. 
                 r1 = max(bspline_sol%knots(i),bspline_sol%knots(j))
                 r2 = min(bspline_sol%knots(i+bspline_sol%order),bspline_sol%knots(j+bspline_sol%order))
                 if (r1 .le. r2) then
                    !Calculate matrix element of the free radial schr. eq. in the B-spline basis
                    bspline_sol%bcoef(j) = 1.0_cfp
                    call cfp_bfqad(integrand, bspline_sol%knots, bspline_sol%bcoef, bspline_sol%n, bspline_sol%order, 0, r1, r2, bspline_sol%tol, quad, err, bspline_sol%work)
                    Tl(j,i,l) = quad
                    !write(stdout,'(3i5,e25.15)') j,i,l,quad
                    bspline_sol%bcoef(j) = 0.0_cfp
                 endif
              enddo !j
           enddo !i

        enddo !l

  end subroutine calculate_Tl_matrix

  real(wp) function wp_eval_radial_free_bspline(data,x) result(r)
      use bspline_base, only: bvalu
      implicit none
      class(radial_free_bspline) :: data
      real(wp), intent(in) :: x

      integer :: ind

         associate(grid => data%bspline_radial_grid)
         if (r_is_inside_A_B(grid%knots(grid%ind),grid%knots(grid%ind+grid%order),real(x,kind=cfp))) then
            grid%bcoef = 0.0_cfp
            grid%bcoef(grid%ind) = 1.0_cfp
            r = bvalu(grid%knots,grid%bcoef,grid%n,grid%order,2,real(x,kind=cfp),grid%inbv,grid%work)
            r = r - grid%l*(grid%l+1)/(x*x)*bvalu(grid%knots,grid%bcoef,grid%n,grid%order,0,real(x,kind=cfp),grid%inbv,grid%work)
            grid%bcoef(grid%ind) = 0.0_cfp
         endif
         end associate

  end function wp_eval_radial_free_bspline

  real(ep1) function ep_eval_radial_free_bspline(data,x) result(r)
      use bspline_base, only: bvalu
      implicit none
      class(radial_free_bspline) :: data
      real(ep1), intent(in) :: x

!      integer :: ind
!
!         associate(grid => data%bspline_radial_grid)
!         if (r_is_inside_A_B(grid%knots(grid%ind),grid%knots(grid%ind+grid%order),x)) then
!            grid%bcoef = 0.0_cfp
!            grid%bcoef(grid%ind) = 1.0_cfp
!            r = bvalu(grid%knots,grid%bcoef,grid%n,grid%order,2,x,grid%inbv,grid%work)
!            r = r - grid%l*(grid%l+1)/(x*x)*bvalu(grid%knots,grid%bcoef,grid%n,grid%order,0,x,grid%inbv,grid%work)
!            grid%bcoef(grid%ind) = 0.0_cfp
!         endif
!         end associate

  end function ep_eval_radial_free_bspline

  !> Solves the Poission equation for a given angular momentum value (l) and for the source term evaluated on a given quadrature grid.
  subroutine solve_poisson_equation(bspline_sol,Tl,source_term,r1,r2,r,w,n_total_points,max_l,l,A,B,bspline_cfs)
     use blas_lapack, only: trsm
     use bspline_base, only: bvalu
     implicit none
     type(bto_data), intent(inout) :: bspline_sol
     integer, intent(in) :: max_l, l, n_total_points
     real(kind=cfp), intent(in) :: source_term(n_total_points), r(n_total_points), w(n_total_points), Tl(bspline_sol%n,bspline_sol%n,0:max_l), r1, r2
     real(kind=cfp), intent(out) :: bspline_cfs(bspline_sol%n,1)
     real(kind=cfp), intent(out) :: A,B

     real(kind=cfp), parameter :: alpha = 1.0_cfp
     integer, parameter :: n = 1
     real(kind=cfp) :: BC1, BC2
     integer :: i, ind, ipiv(bspline_sol%n-2)

        !todo the l-value on input should rather be the maximum l value so here I should get solutions for all l-values at once!!!

        !Calculate the boundary conditions: integrate the source term
        BC1 = 0.0_cfp; BC2 = 0.0_cfp
        do i=1,n_total_points
           BC1 = BC1 + w(i)*source_term(i)*r(i)**(-l+1)
           BC2 = BC2 + w(i)*source_term(i)*r(i)**(l+2)
        enddo !i
        BC1 = BC1*r1**(l+1)
        BC2 = BC2/r2**l

        !Calculate the constants A,B defining the homogeneous solution:
        !todo careful if r1 ~ r2: rewrite the subtractions
        A = (-(BC1*r1**l) + BC2*r2**l)/(-r1**(1 + 2*l) + r2**(1 + 2*l))
        B = (r1**l*r2**l*(BC2*r1**(1 + l) - BC1*r2**(1 + l)))/(r1**(1 + 2*l) - r2**(1 + 2*l))
!        print *,'BC',BC1,BC2
!        print *,'A,B',A,B

        !Calculate projections of the source term on the B-spline basis omitting the first and the last B-splines in the basis:
        bspline_cfs = 0.0_cfp
        bspline_sol%bcoef = 0.0_cfp
        do i=1,n_total_points
           !todo use BSPLVD to calculate values of all B-splines at r(i)
           do ind=2,bspline_sol%n-1
              if (r_is_inside_A_B(bspline_sol%knots(ind),bspline_sol%knots(ind+bspline_sol%order),r(i))) then
                 bspline_sol%bcoef(ind) = 1.0_cfp
                 bspline_cfs(ind,1) = bspline_cfs(ind,1) + w(i)*source_term(i)*r(i)*bvalu(bspline_sol%knots,bspline_sol%bcoef,bspline_sol%n,bspline_sol%order,0,r(i),bspline_sol%inbv,bspline_sol%work)
                 bspline_sol%bcoef(ind) = 0.0_cfp
              endif
           enddo !ind
        enddo !i
        bspline_cfs = -(2*l+1)*bspline_cfs

        !Solve for the B-spline coefficients describing the solution with boundary conditions Y(r1) = Y(r2) = 0
        call dgesv(bspline_sol%n-2,n,Tl(2:bspline_sol%n-1,2:bspline_sol%n-1,l),bspline_sol%n-2,ipiv,bspline_cfs(2:bspline_sol%n-1,1:1),bspline_sol%n-2,i)

  end subroutine solve_poisson_equation
 
  !> Calculates partial wave expansion for the product of a pair of CGTOs.
  subroutine calculate_CGTO_pair_pw_coefficients_numerical(max_l,cgto_A,cgto_B,r,angular_integrals,epsrel,epsabs)
      use phys_const, only: fourpi, twopi, pi
      use special_functions, only: cfp_besi, cfp_eval_poly_horner
      use omp_lib
      implicit none
      integer, intent(in) :: max_l
      type(cgto_data), intent(in) :: cgto_A, cgto_B
      real(kind=cfp), allocatable :: angular_integrals(:,:,:,:)
      real(kind=cfp), intent(in) :: r(:), epsrel, epsabs

      integer, parameter :: kode = 2
      real(kind=cfp), parameter :: half = 0.5_cfp
      real(kind=cfp), parameter :: th_min = 0.0_cfp, th_max = pi, phi_min = 0.0_cfp, phi_max = twopi

      integer :: n, available, i, ind, order, largest_l, j, n_Xlm, last_Xlm, last_j, err, l, m, stride, n_cgto_A_m, n_cgto_B_m, m_A_ind, m_B_ind, lm, point
      real(kind=cfp) :: current_integral, r_A(3)
      real(kind=cfp), allocatable :: gto_vals(:,:,:), r1(:), r2(:), r3(:), w(:), aux(:,:)
      integer, allocatable :: to_converge(:,:)
      logical, allocatable :: converged(:,:,:)
      type(Xlm_x_pair_cgto_surface) :: xlm_pair_cgto_r
      integer :: my, rank, n_threads

         xlm_pair_cgto_r%cgto_A = cgto_A
         xlm_pair_cgto_r%cgto_B = cgto_B

         !Use quadrature on sphere to calculate the angular integral

         !check the largest L-value of the Xlm for which this quadrature can be used
         if (max_l > rule_max-(cgto_A%l+cgto_B%l+1)) call xermsg('bto_gto_integrals_mod','calculate_CGTO_pw_coefficients','max_l too large for the available Lebedev rules.',2,1)

         n_Xlm = (max_l+1)**2
         n_cgto_A_m = 2*cgto_A%l+1
         n_cgto_B_m = 2*cgto_B%l+1

         if (allocated(angular_integrals)) deallocate(angular_integrals)

         allocate(gto_vals(order_table(rule_max),n_cgto_B_m,n_cgto_A_m),to_converge(n_cgto_B_m,n_cgto_A_m),r1(mmax),r2(mmax),r3(mmax),w(mmax),aux(n_cgto_B_m,n_cgto_A_m),&
         &converged(n_Xlm,n_cgto_B_m,n_cgto_A_m),angular_integrals(size(r),n_cgto_B_m,n_cgto_A_m,n_Xlm),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','calculate_CGTO_pw_coefficients','Memory allocation failed.',err,1)

         !Precalculate spherical harmonics on the Lebedev grid
         call precalculate_Xlm_on_Lebedev_grid(max_l)

         !Over all radial points
         do point=1,size(r)

            to_converge(1:n_cgto_B_m,1:n_cgto_A_m) = n_Xlm !for each CGTO m we have to calculate n_Xlm angular integrals
            converged(1:n_Xlm,1:n_cgto_B_m,1:n_cgto_A_m) = .false.
   
            stride = (max_l_precalculated+1)**2
   
            last_j = 0
            !The minimum order of the Lebedev quadrature is given by the sum of orders of the spherical harmonics: their product is a polynomial of order cgto%l+cms_bto%l
            !TODO MAKE SURE WE BRANCH HERE INTO THE ADAPTIVE INTEGRATOR IF QUAD PRECISION IS BEING USED!!!
            do n=cgto_A%l+cgto_B%l+2,rule_max
   
               if (sum(to_converge) .eq. 0) exit !all angular integrals have been calculated and converged
   
               largest_l = n-(cgto_A%l+cgto_B%l+2) !largest L-value of the Xlm for which this quadrature can be used
               last_Xlm = (largest_l+1)**2
   
               !Loop over increasingly precise quadratures until convergence of all angular integrals is reached
               available = available_table(n)
            
               if (available == 1) then
                  !In case of CMS-only functions the order must be high enough to integrate the product of the spherical harmonics otherwise the rule can fall on the nodes of Xlm
                  order = order_table(n)
                  call ld_by_order(order,r1,r2,r3,w)
   
                  !Evaluate the CGTO for all m values at the quadrature points
                  do i=1,order
                     r_A(1:3) = (/r1(i),r2(i),r3(i)/)
                     aux(1:n_cgto_B_m,1:n_cgto_A_m) = gto_pair_eval_R_all_m(cgto_A,cgto_B,r(point),r_A)
                     gto_vals(i,1:n_cgto_B_m,1:n_cgto_A_m) = aux(1:n_cgto_B_m,1:n_cgto_A_m)
                  enddo !i
         
                  !For all CGTO m values
                  do m_A_ind=1,n_cgto_A_m
                     do m_B_ind=1,n_cgto_B_m
   
                        j=0
                        do l=0,min(max_l,largest_l)
                           do m=-l,l
         
                              j = j + 1
         
                              if (to_converge(m_B_ind,m_A_ind) .eq. 0) exit !all integrals have been calculated
                              if (converged(j,m_B_ind,m_A_ind)) cycle !skip the angular integrals that have been converged using a lower-order rule
           
                              !Apply the quadrature rule for the integral dOmega/(4*pi)*Xlm(Omega)*CGTO(Omega,r)
                              current_integral = 0.0_wp
                              do i=1,order
                                 !Use the precalculated real spherical harmonic stored in Xlm_Lebedev. todo transpose the Xlm_Lebedev array so I loop in array-element order
                                 ind = offset(n) + (i-1)*stride + j !j = l*l + l + m + 1
                                 current_integral = current_integral + w(i)*Xlm_Lebedev(ind)*gto_vals(i,m_B_ind,m_A_ind)
                              enddo !i
                              current_integral = current_integral*fourpi
                
                              !Is this the first estimate for the integral with Xlm (that has the index j)?
                              if (last_j .ge. j) then !NO
                                 !Check for relative precision 
                                 if (abs((current_integral-angular_integrals(point,m_B_ind,m_A_ind,j))/current_integral) .le. epsrel/10.0_cfp) then
                                    converged(j,m_B_ind,m_A_ind) = .true.
                                    to_converge(m_B_ind,m_A_ind) = to_converge(m_B_ind,m_A_ind)-1
                                 elseif (n .eq. rule_max) then
                                    !We need to resort to adaptive quadrature
                                    !write(100+iam,'("No convergence of the angular integral using Lebedev grid for r=",e,2i5,i)') r,cgto%l,cgto%m,j
                                    !todo modify so that the adaptive quadrature will converge at the same time integrals for the whole shell of CGTOs (i.e. for all GTO m).
                                    xlm_pair_cgto_r%neval = 0; xlm_pair_cgto_r%ndiv = 0
                                    xlm_pair_cgto_r%cgto_A%m = m_A_ind-cgto_A%l-1 !m value of the CGTO A
                                    xlm_pair_cgto_r%cgto_B%m = m_B_ind-cgto_B%l-1 !m value of the CGTO B
                                    xlm_pair_cgto_r%r = r(point)
                                    xlm_pair_cgto_r%l = l
                                    xlm_pair_cgto_r%m = m
                                    current_integral = quad2d(xlm_pair_cgto_r,th_min,th_max,phi_min,phi_max,epsrel/10.0_cfp)
                                    converged(j,m_B_ind,m_A_ind) = .true.
                                    to_converge(m_B_ind,m_A_ind) = to_converge(m_B_ind,m_A_ind)-1
                                    !write(100+iam,'("Adaptive quadrature finished. neval, ndiv=",2i,e)') xlm_pair_cgto_r%neval,xlm_pair_cgto_r%ndiv,current_integral
                                 endif
                              endif
       
                              angular_integrals(point,m_B_ind,m_A_ind,j) = current_integral
                              lm = l*l+l+m+1
                   
                           enddo !m
                        enddo !l
                        last_j = last_Xlm
   
                     enddo !m_B_ind
                  enddo !m_A_ind
   
               endif
          
            enddo !n

         enddo !point

         do lm=1,n_Xlm
            do m_A_ind=1,n_cgto_A_m
               do m_B_ind=1,n_cgto_B_m
                  do point=1,size(r)
                     write(112,'(e25.15,3i,e25.15)') r(point),m_B_ind-cgto_B%l-1,m_A_ind-cgto_A%l-1,lm,angular_integrals(point,m_B_ind,m_A_ind,lm)
                  enddo !i
               enddo !Mg_B
            enddo !Mg_A
         enddo !lm

  end subroutine calculate_CGTO_pair_pw_coefficients_numerical

  !> Evaluates the product of a pair of CGTOs for all m values at radial distance x and for angular direction given by vector R.
  function gto_pair_eval_R_all_m(cgto_A,cgto_B,x,R)
      implicit none
      type(cgto_data), intent(in) :: cgto_A, cgto_B
      real(kind=cfp), intent(in) :: R(3), x
      real(kind=cfp) :: gto_pair_eval_R_all_m(2*cgto_B%l+1,2*cgto_A%l+1)

      real(kind=cfp) :: sum_exp_A, sum_exp_B, r_A(3), r_B(3), r_A_square, r_B_square
      real(kind=cfp) :: SH_A(-cgto_A%l:cgto_A%l,0:cgto_A%l+1), SH_B(-cgto_B%l:cgto_B%l,0:cgto_B%l+1)
      integer :: j, m_A, m_B

         r_A(1:3) = (/x*R(1)-cgto_A%RA(1),x*R(2)-cgto_A%RA(2),x*R(3)-cgto_A%RA(3)/)
         r_B(1:3) = (/x*R(1)-cgto_B%RA(1),x*R(2)-cgto_B%RA(2),x*R(3)-cgto_B%RA(3)/)

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
         sum_exp_A = 0.0_wp
         do j=1,cgto_A%n_prim
            sum_exp_A = sum_exp_A + cgto_A%ccf(j)*cgto_A%p_norm(j)*exp(-cgto_A%alp(j)*r_A_square)
         enddo
         sum_exp_A = cgto_A%norm*sum_exp_A

         r_B_square = dot_product(r_B,r_B)
         sum_exp_B = 0.0_wp
         do j=1,cgto_B%n_prim
            sum_exp_B = sum_exp_B + cgto_B%ccf(j)*cgto_B%p_norm(j)*exp(-cgto_B%alp(j)*r_B_square)
         enddo
         sum_exp_B = cgto_B%norm*sum_exp_B

         do m_A=-cgto_A%l,cgto_A%l
            do m_B=-cgto_B%l,cgto_B%l
               gto_pair_eval_R_all_m(m_B+cgto_B%l+1,m_A+cgto_A%l+1) = SH_A(m_A,cgto_A%l)*sum_exp_A*SH_B(m_B,cgto_B%l)*sum_exp_B
            enddo !m_A
         enddo !m_b

  end function gto_pair_eval_R_all_m

  !> Evaluates the product of a pair of CGTOs for the given m values at radial distance x and for angular direction given by vector R.
  function gto_pair_eval_R(cgto_A,cgto_B,x,R)
      implicit none
      type(cgto_data), intent(in) :: cgto_A, cgto_B
      real(kind=cfp), intent(in) :: R(3), x
      real(kind=cfp) :: gto_pair_eval_R

      real(kind=cfp) :: sum_exp_A, sum_exp_B, r_A(3), r_B(3), r_A_square, r_B_square
      real(kind=cfp) :: SH_A(-cgto_A%l:cgto_A%l,0:cgto_A%l+1), SH_B(-cgto_B%l:cgto_B%l,0:cgto_B%l+1)
      integer :: j

         r_A(1:3) = (/x*R(1)-cgto_A%RA(1),x*R(2)-cgto_A%RA(2),x*R(3)-cgto_A%RA(3)/)
         r_B(1:3) = (/x*R(1)-cgto_B%RA(1),x*R(2)-cgto_B%RA(2),x*R(3)-cgto_B%RA(3)/)

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
         sum_exp_A = 0.0_wp
         do j=1,cgto_A%n_prim
            sum_exp_A = sum_exp_A + cgto_A%ccf(j)*cgto_A%p_norm(j)*exp(-cgto_A%alp(j)*r_A_square)
         enddo
         sum_exp_A = cgto_A%norm*sum_exp_A

         r_B_square = dot_product(r_B,r_B)
         sum_exp_B = 0.0_wp
         do j=1,cgto_B%n_prim
            sum_exp_B = sum_exp_B + cgto_B%ccf(j)*cgto_B%p_norm(j)*exp(-cgto_B%alp(j)*r_B_square)
         enddo
         sum_exp_B = cgto_B%norm*sum_exp_B

         gto_pair_eval_R = SH_A(cgto_A%m,cgto_A%l)*sum_exp_A*SH_B(cgto_B%m,cgto_B%l)*sum_exp_B

  end function gto_pair_eval_R

  !> Evaluates the product of CGTO and a real spherical harmonic at a given point in space. x,y are the spherical polar coordinates on the sphere.
  function eval_Xlm_x_pair_cgto_surface(this,x,y)
      use phys_const, only: fourpi
      implicit none
      class(Xlm_x_pair_cgto_surface) :: this
      real(kind=cfp), intent(in) :: x, y
      real(kind=cfp) :: eval_Xlm_x_pair_cgto_surface

      real(kind=wp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
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

         eval_Xlm_x_pair_cgto_surface = gto_pair_eval_R(this%cgto_A,this%cgto_B,this%r,r_A)*RH(this%m,this%l)*s

         this%neval = this%neval + 1

  end function eval_Xlm_x_pair_cgto_surface

  !> Assumes that the Xlm have been precalculated on the Lebedev grid and that the Xlm for the nuclei have been precalculated for all l up to l=max_l+cgto%l.
  !> This routine attempts to integrate using the Lebedev grid until the required relative precision. If the relative precision cannot be reached then we first check if the integral is smaller than epsabs.
  !> If it is then we regard it as converged otherwise we resort to the adaptive quadrature - now this is extremely slow and therefore it is used only as a last resort.
  subroutine calculate_CGTO_pw_coefficients_numerical(converged,max_l,cgto,RA,r,angular_integrals,epsrel,epsabs,nuclei,nuclei_RA)
      use phys_const, only: fourpi, twopi, pi
      use special_functions, only: cfp_besi, cfp_eval_poly_horner
      use omp_lib
      implicit none
      logical, allocatable :: converged(:,:)
      integer, intent(in) :: max_l
      type(cgto_data), intent(in) :: cgto
      real(kind=cfp), intent(out) :: angular_integrals(:,:)
      real(kind=cfp), intent(in) :: r, epsrel, epsabs, RA
      type(nucleus_type), optional, intent(in) :: nuclei(:)
      real(kind=cfp), optional, intent(in) :: nuclei_RA(:)

      integer, parameter :: kode = 2
      real(kind=cfp), parameter :: half = 0.5_cfp
      real(kind=cfp), parameter :: th_min = 0.0_cfp, th_max = pi, phi_min = 0.0_cfp, phi_max = twopi

      integer :: n, available, i, ind, order, largest_l, j, n_Xlm, last_Xlm, last_j, err, l, m, stride, n_cgto_m, m_ind, n_nuc, lm, nuc, lp, mp, mg, n_l, lp_mp, Lg_Mg
      real(kind=cfp) :: current_integral, r_A(3), r_square, fac, inv_r, arg, sum_mpp, lambda_cf(0:cgto%l+1), sum_mu, fnu, asym, val, contraction
      real(kind=cfp), allocatable :: gto_vals(:,:), r1(:), r2(:), r3(:), w(:), aux(:), fac_l(:,:)
      integer, allocatable :: to_converge(:)
      type(Xlm_x_cgto_surface) :: xlm_cgto_r
      logical :: do_nuclear_potential
      integer :: my, rank, n_threads

         xlm_cgto_r%cgto = cgto

         if (present(nuclei)) then
            do_nuclear_potential = .true.
            n_nuc = size(nuclei)
         else
            do_nuclear_potential = .false.
         endif

         if (ieor(present(nuclei),present(nuclei_RA))) call xermsg('bto_gto_integrals_mod','calculate_CGTO_pw_coefficients','On input both nuclei_RA and nuclei must be given if one of them is.',1,1)

         !todo the parts of the code within the if statements below should go into separate subroutines to improve clarity.
  
         !Calculate the angular integrals trivially if CGTO is at the CMS - this applies also to nuclear attraction integrals.
         if (RA .eq. 0.0_cfp) then
            angular_integrals = 0.0_cfp

            !Radial part of the CGTO
            r_square = r*r
            contraction = 0.0_cfp
            do j=1,cgto%n_prim
               contraction = contraction + cgto%ccf(j)*cgto%p_norm(j)*exp(-cgto%alp(j)*r_square)
            enddo
            contraction = contraction*cgto%norm*sqrt(fourpi/(2*cgto%l+1.0_cfp))*r**cgto%l

            if (do_nuclear_potential) then

               !Use the Legendre expansion at r to calculate the angular integral.
               n_l = cgto%l+max_l
               allocate(fac_l(0:n_l,n_nuc),stat=err)
               if (err .ne. 0) call xermsg('bto_gto_integrals_mod','calculate_CGTO_pw_coefficients','Memory allocation 2 failed.',err,1)

               !Radial part of the Legendre expansion
               inv_r = 1.0_cfp/r
               do nuc=1,n_nuc
                  do lp=0,cgto%l+max_l
                     fac_l(lp,nuc) = fourpi/(2*lp+1.0_cfp)*((nuclei_RA(nuc)*inv_r)**lp)*inv_r
                  enddo !lp
               enddo !nuc

               !Angular part for the whole shell of CGTOs and all (l,m) projections.
               do mg=-cgto%l,cgto%l
                  fac = (-1)**mg*contraction
                  do l=0,max_l
                     do m=-l,l
                        lm = l*l+l+m+1
                        do nuc=1,n_nuc
                           do lp=abs(cgto%l-l),cgto%l+l
                              do mp=-lp,lp
                                 lp_mp = lp*lp+lp+mp+1
                                 angular_integrals(lm,cgto%l+mg+1) = angular_integrals(lm,cgto%l+mg+1)&
                                                                    &-nuclei(nuc)%charge*fac*fac_l(lp,nuc)*Xlm_nuclei(lp_mp+n_Xlm_nuclei*(nuc-1))*cpl%rgaunt(cgto%l,lp,l,mg,mp,m)
                              enddo !mp
                           enddo !lp
                        enddo !nuc
                     enddo !m
                  enddo !l
               enddo !mg
            else
               !The only projections that are non-zero are those where the real spherical harmonic has the same (l,m) values as the CGTO:
               lm = cgto%l*cgto%l+cgto%l+1
               do m=-cgto%l,cgto%l
                  angular_integrals(lm+m,cgto%l+m+1) = (-1)**m*contraction
               enddo
            endif

            return
         endif

         !Use quadrature on sphere to calculate the angular integral

         !check the largest L-value of the Xlm for which this quadrature can be used
         if (max_l > rule_max-(cgto%l+1)) call xermsg('bto_gto_integrals_mod','calculate_CGTO_pw_coefficients','max_l too large for the available Lebedev rules.',2,1)

         n_Xlm = (max_l+1)**2
         n_cgto_m = 2*cgto%l+1

         converged(1:n_Xlm,1:n_cgto_m) = .false.

         allocate(gto_vals(order_table(rule_max),n_cgto_m),to_converge(n_cgto_m),r1(mmax),r2(mmax),r3(mmax),w(mmax),aux(2*cgto%l+1),stat=err)
         if (err .ne. 0) call xermsg('bto_gto_integrals_mod','calculate_CGTO_pw_coefficients','Memory allocation failed.',err,1)

         to_converge(:) = n_Xlm !for each CGTO m we have to calculate n_Xlm angular integrals

         stride = (max_l_precalculated+1)**2

         last_j = 0
         !The minimum order of the Lebedev quadrature is given by the sum of orders of the spherical harmonics: their product is a polynomial of order cgto%l+cms_bto%l
         !TODO MAKE SURE WE BRANCH HERE INTO THE ADAPTIVE INTEGRATOR IF QUAD PRECISION IS BEING USED!!!
         do n = cgto%l+2,rule_max

            if (sum(to_converge) .eq. 0) exit !all angular integrals have been calculated and converged

            largest_l = n-(cgto%l+2) !largest L-value of the Xlm for which this quadrature can be used
            last_Xlm = (largest_l+1)**2

            !Loop over increasingly precise quadratures until convergence of all angular integrals is reached
            available = available_table(n)
         
            if (available == 1) then
               !In case of CMS-only functions the order must be high enough to integrate the product of the spherical harmonics otherwise the rule can fall on the nodes of Xlm
               order = order_table(n)
               call ld_by_order(order,r1,r2,r3,w)

               !Evaluate the CGTO for all m values at the quadrature points
               if (do_nuclear_potential) then
                  do i=1,order
                     r_A(1:3) = (/r1(i),r2(i),r3(i)/)
                     aux(1:n_cgto_m) = gto_eval_R_all_m(cgto,r,r_A)
                     gto_vals(i,1:n_cgto_m) = aux(1:n_cgto_m)*nuclear_potential(nuclei,n_nuc,r,r_A)
                  enddo !i
               else
                  do i=1,order
                     r_A(1:3) = (/r1(i),r2(i),r3(i)/)
                     aux(1:n_cgto_m) = gto_eval_R_all_m(cgto,r,r_A)
                     gto_vals(i,1:n_cgto_m) = aux(1:n_cgto_m)
                  enddo !i
               endif
      
               !For all CGTO m values
               do m_ind=1,n_cgto_m 
                  j=0
                  do l=0,min(max_l,largest_l)
                     do m=-l,l
   
                        j = j + 1
   
                        if (to_converge(m_ind) .eq. 0) exit !all integrals have been calculated
                        if (converged(j,m_ind)) cycle !skip the angular integrals that have been converged using a lower-order rule
     
                        !Apply the quadrature rule for the integral dOmega/(4*pi)*Xlm(Omega)*CGTO(Omega,r)
                        current_integral = 0.0_wp
                        do i=1,order
                           !Use the precalculated real spherical harmonic stored in Xlm_Lebedev. todo transpose the Xlm_Lebedev array so I loop in array-element order
                           ind = offset(n) + (i-1)*stride + j !j = l*l + l + m + 1
                           current_integral = current_integral + w(i)*Xlm_Lebedev(ind)*gto_vals(i,m_ind)
                        enddo !i
                        current_integral = current_integral*fourpi
          
                        !Is this the first estimate for the integral with Xlm (that has the index j)?
                        if (last_j .ge. j) then !NO
                           !Check for relative precision 
                           if (abs((current_integral-angular_integrals(j,m_ind))/current_integral) .le. epsrel/10.0_cfp) then
                              converged(j,m_ind) = .true.
                              to_converge(m_ind) = to_converge(m_ind)-1
                           elseif (n .eq. rule_max) then
                              !We need to resort to adaptive quadrature
                              !write(100+iam,'("No convergence of the angular integral using Lebedev grid for r=",e,2i5,i)') r,cgto%l,cgto%m,j
                              !todo modify so that the adaptive quadrature will converge at the same time integrals for the whole shell of CGTOs (i.e. for all GTO m).
                              xlm_cgto_r%neval = 0; xlm_cgto_r%ndiv = 0
                              xlm_cgto_r%cgto%m = m_ind-cgto%l-1 !m value of the CGTO
                              xlm_cgto_r%r = r
                              xlm_cgto_r%l = l
                              xlm_cgto_r%m = m
                              current_integral = quad2d(xlm_cgto_r,th_min,th_max,phi_min,phi_max,epsrel/10.0_cfp)
                              converged(j,m_ind) = .true.
                              to_converge(m_ind) = to_converge(m_ind)-1
                              !write(100+iam,'("Adaptive quadrature finished. neval, ndiv=",2i,e)') xlm_cgto_r%neval,xlm_cgto_r%ndiv,current_integral
                           endif
                        endif
 
                        angular_integrals(j,m_ind) = current_integral
             
                     enddo !m
                  enddo !l
                  last_j = last_Xlm

               enddo !m_ind

            endif
       
         enddo !n

  end subroutine calculate_CGTO_pw_coefficients_numerical

  !> Evaluates CGTO for all m values at radial distance x and for angular direction given by vector R.
  function gto_eval_R_all_m(cgto_inp,x,R)
      use phys_const, only: fourpi
      implicit none

      type(cgto_data), intent(in) :: cgto_inp
      real(kind=cfp), intent(in) :: R(3), x
      real(kind=cfp) :: gto_eval_R_all_m(2*cgto_inp%l+1)

      real(kind=cfp) :: sum_exp, r_A(3), r_square
      real(kind=cfp) :: SH(-cgto_inp%l:cgto_inp%l,0:cgto_inp%l+1)
      integer :: j

         r_A(1:3) = (/x*R(1)-cgto_inp%RA(1),x*R(2)-cgto_inp%RA(2),x*R(3)-cgto_inp%RA(3)/)
         r_square = dot_product(r_A,r_A)
         if (cgto_inp%l > 0) then
            call solh(SH,r_A(1),r_A(2),r_A(3),cgto_inp%l)
         else
            SH(0,0) = 1.0_cfp
         endif
         sum_exp = 0.0_wp
         do j=1,cgto_inp%n_prim
            sum_exp = sum_exp + cgto_inp%ccf(j)*cgto_inp%p_norm(j)*exp(-cgto_inp%alp(j)*r_square)
         enddo
         sum_exp = cgto_inp%norm*sum_exp
         gto_eval_R_all_m(1:2*cgto_inp%l+1) = SH(-cgto_inp%l:cgto_inp%l,cgto_inp%l)*sum_exp

  end function gto_eval_R_all_m

  !> Evaluates the nuclear attraction potential for electron at point: x*R(1:3).
  function nuclear_potential(nuclei,n,x,R)
     use common_obj, only: nucleus_type
     implicit none
     integer, intent(in) :: n
     type(nucleus_type), intent(in) :: nuclei(n)
     real(kind=cfp), intent(in) :: R(3), x
     real(kind=cfp) :: nuclear_potential

     integer :: i
     real(kind=cfp) :: r_A(1:3)

        nuclear_potential = 0.0_cfp
        do i=1,n
           r_A(1:3) = (/x*R(1)-nuclei(i)%center(1),x*R(2)-nuclei(i)%center(2),x*R(3)-nuclei(i)%center(3)/)
           nuclear_potential = nuclear_potential - nuclei(i)%charge/sqrt(dot_product(r_A,r_A))
        enddo

  end function nuclear_potential

  !> Evaluates the product of CGTO and a real spherical harmonic at a given point in space. x,y are the spherical polar coordinates on the sphere.
  function eval_Xlm_x_cgto_surface(this,x,y)
      use phys_const, only: fourpi
      implicit none
      class(Xlm_x_cgto_surface) :: this
      real(kind=cfp), intent(in) :: x, y
      real(kind=cfp) :: eval_Xlm_x_cgto_surface

      real(kind=wp), parameter :: norm = 1.0_cfp/sqrt(fourpi)
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

         eval_Xlm_x_cgto_surface = gto_eval_R(this%cgto,this%r,r_A)*RH(this%m,this%l)*s

         this%neval = this%neval + 1

  end function eval_Xlm_x_cgto_surface

  !> Evaluates CGTO at radial distance x and for angular direction given by vector R.
  function gto_eval_R(cgto_inp,x,R)
      use phys_const, only: fourpi
      implicit none

      real(kind=cfp) :: gto_eval_R
      type(cgto_data), intent(in) :: cgto_inp
      real(kind=cfp), intent(in) :: R(3), x

      real(kind=cfp) :: sum_exp, r_A(3), r_square
      real(kind=cfp) :: SH(-cgto_inp%l:cgto_inp%l,0:cgto_inp%l+1)
      integer :: j

         r_A(1:3) = (/x*R(1)-cgto_inp%RA(1),x*R(2)-cgto_inp%RA(2),x*R(3)-cgto_inp%RA(3)/)
         r_square = dot_product(r_A,r_A)
         if (cgto_inp%l > 0) then
            call solh(SH,r_A(1),r_A(2),r_A(3),cgto_inp%l)
         else
            SH(0,0) = 1.0_cfp
         endif
         sum_exp = 0.0_wp
         do j=1,cgto_inp%n_prim
            sum_exp = sum_exp + cgto_inp%ccf(j)*cgto_inp%p_norm(j)*exp(-cgto_inp%alp(j)*r_square)
         enddo
         sum_exp = cgto_inp%norm*sum_exp
         gto_eval_R = SH(cgto_inp%m,cgto_inp%l)*sum_exp

  end function gto_eval_R

  subroutine calc_resh_coefficients(L)
    implicit none
    integer, intent(in) :: L

    integer :: l_it, m_it, ind

       if (L > max_l) then
          max_l = L
          if (allocated(a)) deallocate(a,b,c)
          allocate(a((L+1)**2),b((L+1)**2),c(L+1))
          ind = 0
          do l_it=1,L-1
             c(l_it) = -sqrt((2.0_cfp*l_it+3.0_cfp)/(2.0_cfp*l_it+2.0_cfp))
             do m_it=-l_it,l_it
                a(ind+m_it+l_it+1) = sqrt((4*l_it*l_it+8*l_it+3.0_cfp)/((l_it+1)**2-m_it*m_it))
                b(ind+m_it+l_it+1) = sqrt((l_it*l_it-m_it*m_it)/(4*l_it*l_it-1.0_cfp))
             enddo
             ind = ind + 2*l_it+1
          enddo
       endif

  end subroutine calc_resh_coefficients

  !> Calculates the real spherical harmonics assuming that the input values X,Y,Z lie on the unit sphere.
  !> \warning This routine is not threadsafe unless calc_resh_coefficients was called before the parallel region for high enough L to ensure the coefficients are never recalculated while using this routine.
  subroutine resh(SH,X,Y,Z,L)
    use precisn
    use phys_const, only: fourpi
    implicit none
    integer, intent(in) :: L
    real(kind=cfp), intent(out) :: SH(-L:L,0:L)
    real(kind=cfp), intent(in) :: x, y, z
  
    integer :: l_it, m_it, ind
    real(kind=cfp), parameter :: norm1 = sqrt(3.0_cfp/fourpi), norm0 = 1.0_cfp/sqrt(fourpi)
  
       SH = 0.0_cfp !vectorized
     
       !initialize the starting values
       SH(0,0) = norm0
       SH(-1,0) = 0.0_cfp
     
       !recursion
       !l_it = 0 case: L=1
       SH(1,1) = -norm1*x 
       SH(0,1) = norm1*z
       SH(-1,1)= -norm1*y
   
       !precalculate the coefficients a,b,c
       call calc_resh_coefficients(L)
   
       !diagonal recursions
       ind = 0
       do l_it = 1, L-1
          SH(l_it+1, l_it+1) = c(l_it)*(x*SH(l_it,l_it)-y*SH(-l_it,l_it))
          SH(-l_it-1,l_it+1) = c(l_it)*(y*SH(l_it,l_it)+x*SH(-l_it,l_it))
          !vertical recursions (vectorized)
          !DIR$ SIMD
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

  end subroutine calc_solh_coefficients

  !> Calculates the real solid harmonic.
  !> \warning This routine is not threadsafe unless calc_solh_coefficients was called before the parallel region for high enough L to ensure the coefficients are never recalculated while using this routine.
  subroutine solh(SH,x,y,z,L)
    use precisn
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
   
       call calc_solh_coefficients(L)
   
       !diagonal recursions
       ind = 0
       do l_it = 1, L-1
          SH(l_it+1, l_it+1) = cs(l_it)*(x*SH(l_it,l_it)-y*SH(-l_it,l_it))
          SH(-l_it-1,l_it+1) = cs(l_it)*(y*SH(l_it,l_it)+x*SH(-l_it,l_it))
          !vertical recursions (vectorized)
          !DIR$ SIMD
          do m_it = -l_it,l_it
             SH(m_it,l_it+1) = z*SH(m_it,l_it)*as(ind+m_it+l_it+1)-bs(ind+m_it+l_it+1)*rsq*SH(m_it,l_it-1)
          end do
          ind = ind + 2*l_it+1
       end do
  
  end subroutine solh

end module bto_gto_integrals_mod
