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
module bto_gto_test_routines_mod
  use general_quadrature
  use gto_data_obj
  use bto_data_obj
  use common_obj, only: nucleus_type
  use coupling_obj
  use bto_gto_integrals_mod, only: calc_resh_coefficients, calc_solh_coefficients, radial_grid_cgto_pair, omp_calculate_cgto_pair_pw_coefficients_analytic, omp_calculate_CGTO_pw_coefficients_analytic
  implicit none

  private

  !> Used to evaluate various coupling coefficients.
  type(couplings_type) :: cpl

  !> These are used to keep the largest L for which the coefficients a,b,c and as,bs,cs for the real spherical and real solid harmonic have been precalculated.
  integer :: max_l = -1, max_ls = -1
  !> Coefficients used for calculation of real spherical and real solid harmonics using the routines resh, solh.
  real(kind=cfp), allocatable :: a(:), b(:), c(:), as(:), bs(:), cs(:)

  type, extends(bound_user_function) :: BB_density_obj
     type(bto_data) :: cms_bto
     integer :: ind_A = 0, ind_B = 0
     integer :: l_A = 0, m_A = 0, l_B = 0, m_B = 0
     integer :: l = 0, m = 0
     real(kind=cfp), allocatable :: bcoef_A(:), bcoef_B(:)
     real(kind=cfp) :: norm_A, norm_B, A, B
  contains
     procedure :: init => init_BB_density
     procedure :: wp_eval => wp_eval_BB_density
     procedure :: ep_eval => ep_eval_BB_density
  end type BB_density_obj

  type, extends(bound_user_function) :: GG_density_obj
     type(cgto_data) :: cgto_A, cgto_B
     integer :: l = 0, m = 0
     real(kind=cfp), allocatable, private :: r_points(:), angular_integrals(:,:,:,:)
  contains
     procedure :: init => init_GG_density
     procedure :: wp_eval => wp_eval_GG_density
     procedure :: ep_eval => ep_eval_GG_density
  end type GG_density_obj

  type, extends(bound_user_function) :: BG_density_obj
     type(bto_data) :: cms_bto
     type(cgto_data) :: cgto
     integer :: l = 0, m = 0
     real(kind=cfp), allocatable, private :: r_points(:), angular_integrals(:,:,:)
  contains
     procedure :: init => init_BG_density
     procedure :: wp_eval => wp_eval_BG_density
     procedure :: ep_eval => ep_eval_BG_density
  end type BG_density_obj

  type, extends(function_2d) :: radial_leg_int_obj
     class(bound_user_function), pointer :: density_1 => null()
     class(bound_user_function), pointer :: density_2 => null()
     integer :: l = 0, m = 0
     real(kind=cfp), private :: prev_r1 = -1.0_cfp, prev_val_den_1 = 0.0_cfp
     real(kind=cfp), private :: prev_r2 = -1.0_cfp, prev_val_den_2 = 0.0_cfp
     integer, private :: l_old = -1, m_old = 0
  contains
     procedure :: eval => eval_radial_leg_int_obj
  end type radial_leg_int_obj

  public calc_2el_integral_BG_GG, calc_2el_integral_BG_BG, calc_2el_integral_BB_GG, GG_density_obj, BG_density_obj, BB_density_obj

contains

  subroutine calc_2el_integral_BB_GG(density_1,density_2,max_l_legendre,integral)
     use const, only: epsabs, epsrel
     use gto_function_obj
     use general_quadrature, only: n_10,w_10,x_10
     implicit none
     class(BB_density_obj), target :: density_1
     class(GG_density_obj), target :: density_2
     integer, intent(in) :: max_l_legendre
     real(kind=cfp), intent(out) :: integral

     type(radial_leg_int_obj) :: radial_leg_int
     real(kind=cfp) :: A_1,B_1,A_2,B_2, val, threshold
     real(kind=cfp), allocatable :: r1(:), w1(:)
     integer :: l, m, err, n_rng_knot = 1, n1_total_points

        call density_1%init(max_l_legendre)
        call density_2%init(max_l_legendre)

        radial_leg_int%density_1 => density_1
        radial_leg_int%density_2 => density_2

        A_1 = density_1%A
        B_1 = density_1%B
!
!------ Determine the radial grid appropriate for integrations involving this CGTO shell.
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give bet
        threshold = threshold*epsabs

        call radial_grid_CGTO_pair(density_2%cgto_A,density_2%cgto_B,threshold,density_1%B,x_10,w_10,n_10,n_rng_knot,r1,w1,n1_total_points,A_2,B_2)
        print *,'domain 1',A_1,B_1
        print *,'domain 2',A_2,B_2

        radial_leg_int%l_old = -1

        integral = 0.0_cfp
        do l=0,max_l_legendre
           do m=-l,l
              density_1%l = l; density_1%m = m
              density_2%l = l; density_2%m = m
              radial_leg_int%l = l; radial_leg_int%m = m
              val = cpl%rgaunt(l,density_1%l_A,density_1%l_B,m,density_1%m_A,density_1%m_B)
              if (val .ne. 0.0_cfp) then
                 !todo extrapolation
                 val = quad2d(radial_leg_int,A_1,B_1,A_2,B_2,10e-5_cfp) !10e-7_cfp)!epsrel)
                 integral = integral + val
                 print *,l,m,val !,integral
              endif
           enddo !m
        enddo !l

  end subroutine calc_2el_integral_BB_GG

  subroutine calc_2el_integral_BG_BG(density_1,density_2,max_l_legendre,integral)
     implicit none
     class(BG_density_obj), target :: density_1
     class(BG_density_obj), target :: density_2
     integer, intent(in) :: max_l_legendre
     real(kind=cfp), intent(out) :: integral

     type(radial_leg_int_obj) :: radial_leg_int
     real(kind=cfp) :: A_1,B_1,A_2,B_2, val, radius_1, radius_2
     real(kind=cfp), allocatable :: r1(:), w1(:)
     integer :: l, m, err

        call density_1%init(max_l_legendre)
        call density_2%init(max_l_legendre)

        radial_leg_int%density_1 => density_1
        radial_leg_int%density_2 => density_2

        A_1 = density_1%cms_bto%knots(density_1%cms_bto%ind); B_1 = density_1%cms_bto%knots(density_1%cms_bto%ind+density_1%cms_bto%order)
        A_2 = density_2%cms_bto%knots(density_2%cms_bto%ind); B_2 = density_2%cms_bto%knots(density_2%cms_bto%ind+density_2%cms_bto%order)

        print *,'domain 1',A_1,B_1
        print *,'domain 2',A_2,B_2

        radial_leg_int%l_old = -1

        integral = 0.0_cfp
        do l=0,max_l_legendre
           do m=-l,l
              density_1%l = l; density_1%m = m
              density_2%l = l; density_2%m = m
              radial_leg_int%l = l; radial_leg_int%m = m
              !todo extrapolation
              val = quad2d(radial_leg_int,A_1,B_1,A_2,B_2,10e-4_cfp) !10e-7_cfp)!epsrel)
              integral = integral + val
              print *,l,m,val !,integral
           enddo !m
        enddo !l

  end subroutine calc_2el_integral_BG_BG

  subroutine calc_2el_integral_BG_GG(density_1,density_2,max_l_legendre,integral)
     use const, only: epsabs, epsrel
     use gto_function_obj
     use general_quadrature, only: n_10,w_10,x_10
     implicit none
     class(GG_density_obj), target :: density_1
     class(BG_density_obj), target :: density_2
     integer, intent(in) :: max_l_legendre
     real(kind=cfp), intent(out) :: integral

     type(radial_leg_int_obj) :: radial_leg_int
     real(kind=cfp) :: A_1,B_1,A_2,B_2, val, radius_1, radius_2, threshold
     real(kind=cfp), allocatable :: r1(:), w1(:)
     integer :: l, m, err, n_rng_knot = 1, n1_total_points


        call density_1%init(max_l_legendre)
        call density_2%init(max_l_legendre)

        radial_leg_int%density_1 => density_1
        radial_leg_int%density_2 => density_2
!
!------ Determine the radial grid appropriate for integrations involving this CGTO shell.
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give bet
        threshold = threshold*epsabs

        call radial_grid_CGTO_pair(density_1%cgto_A,density_1%cgto_B,threshold,density_2%cms_bto%B,x_10,w_10,n_10,n_rng_knot,r1,w1,n1_total_points,A_1,B_1)
        print *,'pair domain=',A_1,B_1
        !todo branch into multipole expansion since the domains for GG and GB
        !don't overlap!!!
        if (B_1 .le. density_2%cms_bto%A) then
           print *, "multipole expansion should be used here"
           stop "error"
        endif

        A_2 = density_2%cms_bto%knots(density_2%cms_bto%ind); B_2 = density_2%cms_bto%knots(density_2%cms_bto%ind+density_2%cms_bto%order)

        print *,'domain 1',A_1,B_1
        print *,'domain 2',A_2,B_2

        radial_leg_int%l_old = -1

        integral = 0.0_cfp
        do l=0,max_l_legendre
           do m=-l,l
              density_1%l = l; density_1%m = m
              density_2%l = l; density_2%m = m
              radial_leg_int%l = l; radial_leg_int%m = m
              !todo extrapolation
              val = quad2d(radial_leg_int,A_1,B_1,A_2,B_2,10e-7_cfp)!epsrel)
              integral = integral + val
              print *,l,m,val !,integral
           enddo !m
        enddo !l

  end subroutine calc_2el_integral_BG_GG

!GG density
  subroutine init_GG_density(this,max_l_pw)
     implicit none
     class(GG_density_obj) :: this
     integer, intent(in) :: max_l_pw

     integer :: max_l_aux

        call this%cgto_A%normalize
        call this%cgto_B%normalize

        allocate(this%r_points(1))
        max_l_aux = max_l_pw+max(this%cgto_A%l,this%cgto_B%l)

        !The two routines below must be called before any OpenMP region to ensure thread safety of resh, solh.
        call calc_resh_coefficients(max_l_aux)
        call calc_solh_coefficients(max_l_aux)

        !Precalculate all couplings needed
        print *,'Gaunt'
        call cpl%prec_cgaunt(max_l_aux)

  end subroutine init_GG_density

  real(wp) function wp_eval_GG_density(data,x)
     implicit none
     class(GG_density_obj) :: data
     real(kind=wp), intent(in) :: x

     integer :: lm, Mg_B_ind, Mg_A_ind

        data%r_points = x

        call omp_calculate_CGTO_pair_pw_coefficients_analytic(data%l,data%cgto_A,data%cgto_B,data%r_points,data%angular_integrals)

        lm = data%l*data%l + data%l+data%m+1
        Mg_A_ind = data%cgto_A%l+data%cgto_A%m+1
        Mg_B_ind = data%cgto_B%l+data%cgto_B%m+1

        !This includes Jac:
        wp_eval_GG_density = x*x*data%angular_integrals(1,Mg_B_ind,Mg_A_ind,lm)

  end function wp_eval_GG_density

  real(ep1) function ep_eval_GG_density(data,x)
     implicit none
     class(GG_density_obj) :: data
     real(kind=ep1), intent(in) :: x

        stop "not implemented"

  end function ep_eval_GG_density
!
!BB density
  subroutine init_BB_density(this,max_l_pw)
     implicit none
     class(BB_density_obj) :: this
     integer, intent(in) :: max_l_pw

     integer :: max_l_aux

        this%cms_bto%ind = this%ind_A
        call this%cms_bto%init_grid
        allocate(this%bcoef_A,source=this%cms_bto%bcoef)

        call this%cms_bto%normalize
        this%norm_A = this%cms_bto%norm

        this%cms_bto%ind = this%ind_B
        call this%cms_bto%init_grid
        allocate(this%bcoef_B,source=this%cms_bto%bcoef)

        call this%cms_bto%normalize
        this%norm_B = this%cms_bto%norm

        this%A = max(this%cms_bto%knots(this%ind_A),this%cms_bto%knots(this%ind_B))
        this%B = min(this%cms_bto%knots(this%ind_A+this%cms_bto%order),this%cms_bto%knots(this%ind_B+this%cms_bto%order))

        max_l_aux = max(max_l_pw,this%l_A,this%l_B)

        !Precalculate all couplings needed
        print *,'Gaunt'
        call cpl%prec_cgaunt(max_l_aux)

  end subroutine init_BB_density

  real(wp) function wp_eval_BB_density(data,x)
     use bspline_base, only: bvalu
     implicit none
     class(BB_density_obj) :: data
     real(kind=wp), intent(in) :: x

     real(kind=cfp) :: x_cfp

        wp_eval_BB_density = 0.0_cfp
        if (x .ge. data%A .and. x .le. data%B) then

           wp_eval_BB_density = cpl%rgaunt(data%l,data%l_A,data%l_B,data%m,data%m_A,data%m_B)
           if (wp_eval_BB_density .ne. 0.0_cfp) then
              x_cfp = x
              !This includes Jac/r/r:
              wp_eval_BB_density = wp_eval_BB_density*data%norm_A*data%norm_B*&
              &bvalu(data%cms_bto%knots,data%bcoef_A,data%cms_bto%n,data%cms_bto%order,0,x_cfp,data%cms_bto%inbv,data%cms_bto%work)*&
              &bvalu(data%cms_bto%knots,data%bcoef_B,data%cms_bto%n,data%cms_bto%order,0,x_cfp,data%cms_bto%inbv,data%cms_bto%work)
           endif

        endif

  end function wp_eval_BB_density

  real(ep1) function ep_eval_BB_density(data,x)
     implicit none
     class(BB_density_obj) :: data
     real(kind=ep1), intent(in) :: x

        stop "not implemented"

  end function ep_eval_BB_density

!BG density
  subroutine init_BG_density(this,max_l_pw)
     implicit none
     class(BG_density_obj) :: this
     integer, intent(in) :: max_l_pw

     integer :: max_l_aux

        call this%cgto%normalize

        call this%cms_bto%init_grid
        call this%cms_bto%normalize

        this%cms_bto%bcoef = 0.0_cfp
        this%cms_bto%bcoef(this%cms_bto%ind) = 1.0_cfp

        allocate(this%r_points(1))
        max_l_aux = max_l_pw+this%cms_bto%l+this%cgto%l

        !The two routines below must be called before any OpenMP region to ensure thread safety of resh, solh.
        call calc_resh_coefficients(max_l_aux)
        call calc_solh_coefficients(max_l_aux)

        !Precalculate all couplings needed
        print *,'Gaunt'
        call cpl%prec_cgaunt(max_l_aux)

  end subroutine init_BG_density

  real(wp) function wp_eval_BG_density(data,x)
     use bspline_base, only: bvalu
     implicit none
     class(BG_density_obj) :: data
     real(kind=wp), intent(in) :: x

     integer :: lm, Mg_ind, lp, mp

        wp_eval_BG_density = 0.0_cfp
        if (x .ge. data%cms_bto%knots(data%cms_bto%ind) .and. x .le. data%cms_bto%knots(data%cms_bto%ind+data%cms_bto%order)) then

           data%r_points = x
   
           call omp_calculate_CGTO_pw_coefficients_analytic(data%l+data%cms_bto%l,data%cgto,data%r_points,data%angular_integrals)
           Mg_ind = data%cgto%l+data%cgto%m+1

           do lp=abs(data%l-data%cms_bto%l),data%l+data%cms_bto%l
              lm = lp*lp+lp+1
              do mp=-lp,lp
                 wp_eval_BG_density = wp_eval_BG_density + cpl%rgaunt(data%l,lp,data%cms_bto%l,data%m,mp,data%cms_bto%m)*data%angular_integrals(1,Mg_ind,lm+mp)
              enddo !mp
           enddo !lp
  
           !This includes Jac/r:
           wp_eval_BG_density = x*wp_eval_BG_density*data%cms_bto%norm*bvalu(data%cms_bto%knots,data%cms_bto%bcoef,data%cms_bto%n,data%cms_bto%order,0,data%r_points(1),data%cms_bto%inbv,data%cms_bto%work)
        endif

  end function wp_eval_BG_density

  real(ep1) function ep_eval_BG_density(data,x)
     implicit none
     class(BG_density_obj) :: data
     real(kind=ep1), intent(in) :: x

        stop "not implemented"

  end function ep_eval_BG_density
!
  real(cfp) function eval_radial_leg_int_obj(this,x,y)
     use phys_const, only: fourpi
     implicit none
     class(radial_leg_int_obj) :: this
     real(kind=cfp), intent(in) :: x, y

     real(kind=cfp) :: den_1, den_2

        !Make sure that when l,m values change that the prev_r1, prev_r2 are reset to negative values otherwise I may be retrieving old (but incompatible) previous values of the density.
        if (this%l_old .ne. this%l .or. this%m_old .ne. this%m) then
           this%prev_r1 = -1.0_cfp; this%prev_val_den_1 = 0.0_cfp
           this%prev_r2 = -1.0_cfp; this%prev_val_den_2 = 0.0_cfp
           this%l_old = this%l; this%m_old = this%m
        endif

        if (x .eq. this%prev_r1) then
           den_1 = this%prev_val_den_1
        else
           !the l,m values for this density must be set externally since density_1 is pointer to a generic bound_user_function
           den_1 = this%density_1%eval(x)
           !save the latest values of r1 and the corresponding density
           this%prev_r1 = x
           this%prev_val_den_1 = den_1
        endif

        if (y .eq. this%prev_r2) then
           den_2 = this%prev_val_den_2
        else
           !the l,m values for this density must be set externally since density_2 is pointer to a generic bound_user_function
           den_2 = this%density_2%eval(y)
           !save the latest values of r1 and the corresponding density
           this%prev_r2 = y
           this%prev_val_den_2 = den_2
        endif

        if (y .ge. x) then
           eval_radial_leg_int_obj = ((x/y)**this%l)/y
        else
           eval_radial_leg_int_obj = ((y/x)**this%l)/x
        endif
!        write(stdout,'(4e25.15)') x,y,den_1,den_2
        eval_radial_leg_int_obj = eval_radial_leg_int_obj*fourpi/(2*this%l+1.0_cfp)*den_1*den_2
!        write(stdout,'(3e25.15)') x,y,eval_radial_leg_int_obj
!        write(stdout,'(2i4,3e25.15)') this%l,this%m,x,y,eval_radial_leg_int_obj

  end function eval_radial_leg_int_obj

end module bto_gto_test_routines_mod
