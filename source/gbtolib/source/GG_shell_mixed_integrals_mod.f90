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
module GG_shell_mixed_integrals_gbl
  use basis_data_generic_gbl
  use coupling_obj_gbl
  use precisn_gbl
  use const_gbl, only: stdout, level2, level3
  use utils_gbl, only: xermsg
  use cgto_pw_expansions_gbl, only: CGTO_shell_pw_expansion_obj, pw_expansion_obj, calculate_lambda_couplings
  use special_functions_gbl, only: real_harmonics_obj
  use grid_gbl, only: grid_r1_r2, grid_r1_r2_obj
  implicit none
  
  private

  public GG_shell_integrals_obj, CMS_CGTO_on_grid

  !> Object that can be used to evaluate mixed BTO/CGTO 2-electron integrals involving at least two shells of CGTOs, i.e. BBGG, BGBG and BGGG classes of integrals.
  !> This object can only be used once the grid_r1_r2 from grid_mod has been initialized.
  type :: GG_shell_integrals_obj
     !> Data for the two CGTO shells which will be used in the (BB|GG), (BG|BG) and (BG|GG) integrals. For the BGGG class the two shells represent the pair of shells of the second electron.
     type(CGTO_shell_data_obj) :: cgto_shell_A, cgto_shell_B
     !> Starting indices for the functions in shells A and B.
     integer, private :: shell_A_starting_index = -1, shell_B_starting_index = -1
     !> Set to .true. by init if the respective shell is centered on the CMS.
     logical, private :: A_is_cms = .false., B_is_cms = .false.
     !> Final set of shell integrals for the BBGG or BGGG (BXGG) and BGBG classes.
     real(kind=cfp), allocatable, private :: BXGG(:,:), BGBG(:,:)
     !> Needed for the Lebedev method:
     real(kind=cfp), allocatable, private :: coulomb_integrals(:,:,:)
     !> Specifies in which order (ma first of mb first) are the integrals in coulomb_integrals stored.
     logical :: order_ab = .false.
     !> Needed for eval_BGGG class: it represents the CGTO shell for the first electron.
     type(CGTO_shell_data_obj), private :: cgto_shell
     !> Needed for eval_BGGG class: starting index for the functions in cgto_shell.
     integer, private :: shell_starting_index = -1
     !> Auxiliary needed for get_BBGG_shell_integrals.
     real(kind=cfp), allocatable, private :: couplings(:)
     !> Set to .true. following a call to init.
     logical, private :: initialized = .false.
  contains
     procedure :: init
     procedure :: final
     procedure :: eval_coulomb_integrals
     procedure :: eval_BGGG_integrals
     procedure :: eval_BBGG_integrals
     procedure :: eval_BGBG_integrals
     procedure :: eval_BGBG_integrals_direct
     procedure :: eval_GG_CCTT_prototype_integrals
     procedure :: eval_BG_CCTT_prototype_integrals
     procedure :: eval_BB_CCTT_prototype_integrals
     procedure :: get_BBGG_shell_integrals
     procedure :: get_BGGG_shell_integrals
     procedure :: get_BGBG_shell_integrals
     procedure, private :: eval_BGBG_integrals_cms_G
     procedure, private :: sph_BXGG_shell_integrals
     procedure, private :: Boys_projected_on_X_lm
     procedure, private :: Laguerre_GTO_projected_on_X_lm
     procedure, private :: Boys_projected_on_G_X_lm
     procedure, private :: Laguerre_GTO_projected_on_G_X_lm
     procedure, private :: precalculate_solh_translation_coeffs
     procedure, private :: precalculate_pair_solh_translation_coeffs
     !procedure, private :: calculate_lambda_couplings
  end type GG_shell_integrals_obj

  !> Used to evaluate the real spherical and solid harmonics. PRIVATE FOR ALL THREADS.
  type(real_harmonics_obj) :: real_harmonics

  !> Used to evaluate various coupling coefficients. PRIVATE FOR ALL THREADS.
  type(couplings_type) :: cpl

  !$OMP THREADPRIVATE(cpl,real_harmonics)

contains

  subroutine init(this,cgto_shell_A,shell_A_starting_index,cgto_shell_B,shell_B_starting_index)
     implicit none
     class(GG_shell_integrals_obj) :: this
     type(CGTO_shell_data_obj) :: cgto_shell_A, cgto_shell_B
     integer, intent(in) :: shell_A_starting_index, shell_B_starting_index

     real(kind=cfp) :: RA, RB, tol

        this%cgto_shell_A = cgto_shell_A
        this%shell_A_starting_index = shell_A_starting_index

        this%cgto_shell_B = cgto_shell_B
        this%shell_B_starting_index = shell_B_starting_index

        RA = sqrt(dot_product(cgto_shell_A%center,cgto_shell_A%center))
        RB = sqrt(dot_product(cgto_shell_B%center,cgto_shell_B%center))

        tol = F1MACH(4,cfp_dummy)

        this%A_is_cms = .false.
        this%B_is_cms = .false.

        if (RA .le. tol) this%A_is_cms = .true.
        if (RB .le. tol) this%B_is_cms = .true.

        !Make sure that the coupling coefficients will always get precalculated if requested for L larger than available from the precalculated buffer.
        cpl%always_precalculate = .true.

        this%initialized = .true.

  end subroutine init

  subroutine final(this)
     implicit none
     class(GG_shell_integrals_obj) :: this

        this%initialized = .false.

        this%shell_A_starting_index = -1
        this%shell_B_starting_index = -1
        this%shell_starting_index = -1
        this%A_is_cms = .false.
        this%B_is_cms = .false.

        if (allocated(this%BXGG)) deallocate(this%BXGG)
        if (allocated(this%BGBG)) deallocate(this%BGBG)
        if (allocated(this%coulomb_integrals)) deallocate(this%coulomb_integrals)
        if (allocated(this%couplings)) deallocate(this%couplings)

  end subroutine final

  subroutine get_BBGG_shell_integrals(this,bto_shell_1,bto_shell_2,starting_index_1,starting_index_2,&
                                      indexing_method,two_el_column,int_index,integrals,n_unique)
     use gto_routines_gbl, only: index_2el, reorder_and_index_2el
     implicit none
     class(GG_shell_integrals_obj) :: this
     type(BTO_shell_data_obj), intent(in) :: bto_shell_1,bto_shell_2
     integer, intent(in) :: starting_index_1, starting_index_2, two_el_column,indexing_method
     !We assume that these three arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)
     integer, intent(out) :: n_unique

     integer :: n_cgto_CD_m, lm, m, ma, mb, m_CD, l, err, ind, pair_index, l_min, l_max, lm_start, lm_end, &
                bto_ind_1, bto_ind_2, BB_index, n_bto_BA_m

        n_cgto_CD_m = (2*this%cgto_shell_A%l+1)*(2*this%cgto_shell_B%l+1)
        n_bto_BA_m = (2*bto_shell_2%l+1)*(2*bto_shell_1%l+1)
        n_unique = n_cgto_CD_m*n_bto_BA_m

        integrals(1:n_unique,two_el_column) = 0.0_cfp

        l_min = abs(bto_shell_1%l-bto_shell_2%l)
        l_max = bto_shell_1%l+bto_shell_2%l

        lm_start = (l_min-1 +1)**2+1
        lm_end =   (l_max   +1)**2

        if (allocated(this%couplings)) then
            if (size(this%couplings) < lm_end) deallocate(this%couplings)
        endif

        if (.not. allocated(this%couplings)) then
           allocate(this%couplings(lm_end),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','get_BBGG_shell_integrals','Memory allocation 1 failed.',err,1)
           endif
        endif

        bto_ind_1 = bto_shell_1%bspline_index
        bto_ind_2 = bto_shell_2%bspline_index

        pair_index = grid_r1_r2%pairs_to_unique_pairs(bto_ind_1,bto_ind_2)

        if (pair_index > 0) then !The pair index is > 0 only if the B-splines overlap.

           ind = 0
           do ma=-bto_shell_1%l,bto_shell_1%l
              do mb=-bto_shell_2%l,bto_shell_2%l
                 do l=l_min,l_max
                    do m=-l,l
                       lm = l*l+l+m+1
                       this%couplings(lm) = cpl%rgaunt(l,bto_shell_1%l,bto_shell_2%l,m,ma,mb)
                    enddo !m
                 enddo !l
                 do m_CD=1,n_cgto_CD_m
                    do lm=lm_start,lm_end
                       BB_index = pair_index + grid_r1_r2%n_unique_pairs*(lm-1)
                       integrals(ind+m_CD,two_el_column) = integrals(ind+m_CD,two_el_column)&
                                                         + this%couplings(lm)*this%BXGG(BB_index,m_CD)
                    enddo !lm
                 enddo !m_CD
                 ind = ind + n_cgto_CD_m
              enddo !mb
           enddo !ma

        endif

        !Compute indices
        if (indexing_method .eq. 2) then
           ! THERE IS A MISTAKE IN THE NEW VERSION OF REORDERING
           !call reorder_and_index_2el(this%cgto_shell_B%l,&
           !                           this%cgto_shell_A%l,&
           !                           bto_shell_2%l,&
           !                           bto_shell_1%l,&
           !                           this%shell_B_starting_index,&
           !                           this%shell_A_starting_index,&
           !                           starting_index_2,&
           !                           starting_index_1,&
           !                           two_el_column,&
           !                           0,&
           !                           integrals,&
           !                           n_unique,&
           !                           int_index)

           call reorder_and_index_2el(this%cgto_shell_B%l,&
                                      this%cgto_shell_A%l,&
                                      bto_shell_2%l,&
                                      bto_shell_1%l,&
                                      this%shell_B_starting_index,&
                                      this%shell_A_starting_index,&
                                      starting_index_2,&
                                      starting_index_1,&
                                      two_el_column,&
                                      int_index,&
                                      integrals)
        else
           call index_2el(this%cgto_shell_B%l,&
                          this%cgto_shell_A%l,&
                          bto_shell_2%l,&
                          bto_shell_1%l,&
                          this%shell_B_starting_index,&
                          this%shell_A_starting_index,&
                          starting_index_2,&
                          starting_index_1,&
                          int_index,.false.,.false.)
        endif

  end subroutine get_BBGG_shell_integrals

  subroutine get_BGGG_shell_integrals(this,bto_shell,starting_index,indexing_method,two_el_column,int_index,integrals,n_unique)
     use gto_routines_gbl, only: index_2el, reorder_and_index_2el
     implicit none
     class(GG_shell_integrals_obj) :: this
     type(BTO_shell_data_obj), intent(in) :: bto_shell
     integer, intent(in) :: starting_index, two_el_column, indexing_method
     !We assume that these two arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)
     integer, intent(out) :: n_unique

     integer :: max_l, i, k, n, ma,lp,mp, l,m,ind,err,n_shell_A,n_shell_B,n_shell_C,n_shell_D,n_shell_CD
     integer :: m_DC, lm, min_l, n1_points, cgto_m, base, BG_index, m_CD, l_min, l_max, lm_start, lm_end
     real(kind=cfp) :: RB, tol

        n_shell_A = 2*bto_shell%l+1
        n_shell_B = 2*this%cgto_shell%l+1
        n_shell_C = 2*this%cgto_shell_A%l+1
        n_shell_D = 2*this%cgto_shell_B%l+1
        n_shell_CD = n_shell_C*n_shell_D

        n_unique = n_shell_A*n_shell_B*n_shell_C*n_shell_D
        integrals(1:n_unique,two_el_column) = 0.0_cfp

        RB = sqrt(dot_product(this%cgto_shell%center,this%cgto_shell%center))
        tol = F1MACH(4,cfp_dummy)

        if (RB .le. tol) then !The CGTO shell B is sitting on the CMS so we compute the integral just like the BBGG class.

           l_min = abs(bto_shell%l-this%cgto_shell%l)
           l_max = bto_shell%l+this%cgto_shell%l
   
           lm_start = (l_min-1 +1)**2+1
           lm_end =   (l_max   +1)**2

           if (allocated(this%couplings)) then
              if (size(this%couplings) < lm_end) deallocate(this%couplings)
           endif
   
           if (.not. allocated(this%couplings)) then
              allocate(this%couplings(lm_end),stat=err)
              if (err /= 0) then
                 call xermsg('bto_gto_integrals_mod','get_BGGG_shell_integrals','Memory allocation 1 failed.',err,1)
              endif
           endif

           ind = 0
           do ma=-bto_shell%l,bto_shell%l
              do cgto_m=-this%cgto_shell%l,this%cgto_shell%l
                 do l=l_min,l_max
                    do m=-l,l
                       lm = l*l+l+m+1
                       this%couplings(lm) = cpl%rgaunt(l,bto_shell%l,this%cgto_shell%l,m,ma,cgto_m)
                    enddo !m
                 enddo !l
                 do m_CD=1,n_shell_CD
                    do lm=lm_start,lm_end
                       BG_index = bto_shell%bspline_index + grid_r1_r2%last_bspline_inner*(lm-1)
                       integrals(ind+m_CD,two_el_column) = integrals(ind+m_CD,two_el_column)&
                                                         + this%couplings(lm)*this%BXGG(BG_index,m_CD)
                    enddo !lm
                    !print *,ind+m_CD,integrals(ind+m_CD,two_el_column)
                 enddo !m_CD
                 ind = ind + n_shell_CD
              enddo !cgto_m
           enddo !ma

        else !We only need to retrieve the integrals that have been precalculated

           base = grid_r1_r2%last_bspline_inner*(grid_r1_r2%max_bspline_l+1)**2
           ind = 0
           do ma=-bto_shell%l,bto_shell%l
              lm = bto_shell%l*bto_shell%l+bto_shell%l+ma+1
              do cgto_m=1,n_shell_B
                 !The BG indices are effectively a 3D array: (bspline_index,bspline_lm,CGTO_B_m)
                 BG_index = bto_shell%bspline_index + grid_r1_r2%last_bspline_inner*(lm-1) + base*(cgto_m-1)
                 do m_CD=1,n_shell_CD
                    integrals(ind+m_CD,two_el_column) = this%BXGG(BG_index,m_CD)
                 enddo !m_CD
                 ind = ind + n_shell_CD
              enddo !cgto_m
           enddo !ma
        endif

        if (indexing_method .eq. 2) then
           ! THERE IS A MISTAKE IN THE NEW VERSION OF REORDERING
           !call reorder_and_index_2el(this%cgto_shell_B%l,&
           !                           this%cgto_shell_A%l,&
           !                           this%cgto_shell%l,&
           !                           bto_shell%l,&
           !                           this%shell_B_starting_index,&
           !                           this%shell_A_starting_index,&
           !                           this%shell_starting_index,&
           !                           starting_index,&
           !                           two_el_column,&
           !                           0,&
           !                           integrals,&
           !                           n_unique,&
           !                           int_index)

           call reorder_and_index_2el(this%cgto_shell_B%l,&
                                      this%cgto_shell_A%l,&
                                      this%cgto_shell%l,&
                                      bto_shell%l,&
                                      this%shell_B_starting_index,&
                                      this%shell_A_starting_index,&
                                      this%shell_starting_index,&
                                      starting_index,&
                                      two_el_column,&
                                      int_index,&
                                      integrals)
        else
           call index_2el(this%cgto_shell_B%l,&
                          this%cgto_shell_A%l,&
                          this%cgto_shell%l,&
                          bto_shell%l,&
                          this%shell_B_starting_index,&
                          this%shell_A_starting_index,&
                          this%shell_starting_index,&
                          starting_index,&
                          int_index,.false.,.false.)
        endif

  end subroutine get_BGGG_shell_integrals

  subroutine get_BGBG_shell_integrals(this,bto_shell_1,bto_shell_2,starting_index_1,starting_index_2,&
                                      indexing_method,two_el_column,int_index,integrals,n_unique)
     use gto_routines_gbl, only: index_2el, reorder_and_index_2el
     implicit none
     class(GG_shell_integrals_obj) :: this
     type(BTO_shell_data_obj), intent(in) :: bto_shell_1,bto_shell_2
     integer, intent(in) :: starting_index_1, starting_index_2, two_el_column,indexing_method
     !We assume that these three arrays have been allocated to the appropriate dimensions:
     integer, allocatable :: int_index(:,:)
     real(kind=cfp), allocatable :: integrals(:,:)
     integer, intent(out) :: n_unique

     integer :: n_shell_1,n_shell_B,n_shell_A,n_shell_2,j1,j2
     integer :: ind,ind_bspl,ind_cgto,n_lbmb,lb1mb1,lb2mb2,m1,m2,mb,ma

        n_shell_1 = 2*bto_shell_1%l+1
        n_shell_A = 2*this%cgto_shell_A%l+1
        n_shell_2 = 2*bto_shell_2%l+1
        n_shell_B = 2*this%cgto_shell_B%l+1
        j1 = bto_shell_1%bspline_index
        j2 = bto_shell_2%bspline_index

        n_unique = n_shell_1*n_shell_B*n_shell_A*n_shell_2
        integrals(1:n_unique,two_el_column) = 0.0_cfp

        n_lbmb = (grid_r1_r2%max_bspline_l+1)**2
        do m1=-bto_shell_1%l,bto_shell_1%l
           lb1mb1 = bto_shell_1%l*bto_shell_1%l+bto_shell_1%l+m1+1
           do m2=-bto_shell_2%l,bto_shell_2%l
              lb2mb2 = bto_shell_2%l*bto_shell_2%l+bto_shell_2%l+m2+1
              ind_bspl = j1 + grid_r1_r2%last_bspline_inner*(lb1mb1-1) + grid_r1_r2%last_bspline_inner*n_lbmb*(j2-1)&
                            + grid_r1_r2%last_bspline_inner**2*n_lbmb*(lb2mb2-1)
              do mb=1,n_shell_B
                 do ma=1,n_shell_A
                    ind_cgto = ma + n_shell_A*(mb-1)
                    ind = mb + n_shell_B*(m2+bto_shell_2%l) + n_shell_2*n_shell_B*(ma-1)&
                             + n_shell_2*n_shell_B*n_shell_A*(m1+bto_shell_1%l)
                    integrals(ind,two_el_column) = this%BGBG(ind_bspl,ind_cgto)
                 enddo !ma
              enddo !mb
           enddo !m2
        enddo !m1

        !Compute indices
        if (indexing_method .eq. 2) then
           ! THERE IS A MISTAKE IN THE NEW VERSION OF REORDERING
           !call reorder_and_index_2el(this%cgto_shell_B%l,&
           !                           bto_shell_2%l,&
           !                           this%cgto_shell_A%l,&
           !                           bto_shell_1%l,&
           !                           this%shell_B_starting_index,&
           !                           starting_index_2,&
           !                           this%shell_A_starting_index,&
           !                           starting_index_1,&
           !                           two_el_column,&
           !                           0,&
           !                           integrals,&
           !                           n_unique,&
           !                           int_index)

           call reorder_and_index_2el(this%cgto_shell_B%l,&
                                      bto_shell_2%l,&
                                      this%cgto_shell_A%l,&
                                      bto_shell_1%l,&
                                      this%shell_B_starting_index,&
                                      starting_index_2,&
                                      this%shell_A_starting_index,&
                                      starting_index_1,&
                                      two_el_column,&
                                      int_index,&
                                      integrals)
        else
           call index_2el(this%cgto_shell_B%l,&
                          bto_shell_2%l,&
                          this%cgto_shell_A%l,&
                          bto_shell_1%l,&
                          this%shell_B_starting_index,&
                          starting_index_2,&
                          this%shell_A_starting_index,&
                          starting_index_1,&
                          int_index,.false.,.false.)
        endif

  end subroutine get_BGBG_shell_integrals
 
  subroutine eval_coulomb_integrals(this)
     use const_gbl, only: epsrel, epsabs
     use phys_const_gbl, only: fourpi
     use cgto_hgp_gbl, only: sph_nari
     implicit none
     class(GG_shell_integrals_obj) :: this

     integer :: err, i, j, k, ind_a, ind_b, n_integrals
     real(kind=cfp), allocatable :: nari(:)
     real(kind=cfp) :: point(3)
     integer, allocatable :: int_index(:,:)

        if (.not. this%initialized) then
           call xermsg('GG_shell_mixed_integrals_mod','eval_coulomb_integrals','Object not initialized.',1,1)
        endif

        n_integrals = (2*this%cgto_shell_A%l+1)*(2*this%cgto_shell_B%l+1)

        if (allocated(this%coulomb_integrals)) deallocate(this%coulomb_integrals)
        allocate(nari(n_integrals),int_index(2,n_integrals),&
                 this%coulomb_integrals(grid_r1_r2%lebedev_order,grid_r1_r2%n1_total_points,n_integrals),stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_coulomb_integrals','Memory allocation failed.',err,1)

        if (this%cgto_shell_A%l .ge. this%cgto_shell_B%l) then
           this%order_AB = .true.
        else
           this%order_AB = .false.
        endif

        ind_a = 1
        ind_b = 1+2*this%cgto_shell_A%l+1
        do i=1,grid_r1_r2%n1_total_points
           do j=1,grid_r1_r2%lebedev_order
              point(1:3) = (/grid_r1_r2%r1(i)*grid_r1_r2%leb_r1(j),&
                             grid_r1_r2%r1(i)*grid_r1_r2%leb_r2(j),&
                             grid_r1_r2%r1(i)*grid_r1_r2%leb_r3(j)/)

              call sph_nari(this%cgto_shell_A%number_of_primitives,&
                            this%cgto_shell_A%center(1),&
                            this%cgto_shell_A%center(2),&
                            this%cgto_shell_A%center(3),&
                            this%cgto_shell_A%norm,&
                            this%cgto_shell_A%norms,&
                            this%cgto_shell_A%l,&
                            this%cgto_shell_A%exponents,&
                            this%cgto_shell_A%contractions,&
                            ind_a,&
                            this%cgto_shell_B%number_of_primitives,&
                            this%cgto_shell_B%center(1),&
                            this%cgto_shell_B%center(2),&
                            this%cgto_shell_B%center(3),&
                            this%cgto_shell_B%norm,&
                            this%cgto_shell_B%norms,&
                            this%cgto_shell_B%l,&
                            this%cgto_shell_B%exponents,&
                            this%cgto_shell_B%contractions,&
                            ind_b,&
                            point(1),point(2),point(3), nari,int_index)

              this%coulomb_integrals(j,i,1:n_integrals) = nari(1:n_integrals)*grid_r1_r2%leb_w(j)*fourpi
              
           enddo !j
        enddo !i

  end subroutine eval_coulomb_integrals

  subroutine eval_lambda_BG_cms_G(cgto_shell,l_max,lambda_klj)
     use phys_const_gbl, only: fourpi
     implicit none
     integer, intent(in) :: l_max
     type(CGTO_shell_data_obj) :: cgto_shell
     real(kind=cfp), allocatable :: lambda_klj(:,:,:) !(r1,l,j)

     integer :: err, l, j, k, k2, B_start, B_end
     real(kind=cfp) :: r12, fac
     real(kind=cfp), allocatable :: cgto_on_grid(:), f_l_w1_w2(:), r2_B(:,:), inv_r1(:), inv_r2(:)

        if (allocated(lambda_klj)) deallocate(lambda_klj)
        allocate(lambda_klj(grid_r1_r2%n1_total_points,0:max(l_max,1),&
                 grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner),&
                 f_l_w1_w2(grid_r1_r2%n2_total_points),&
                 r2_B(grid_r1_r2%n2_total_points,grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner),&
                 inv_r1(grid_r1_r2%n1_total_points),&
                 inv_r2(grid_r1_r2%n2_total_points),stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_lambda_BG_cms_G','Memory allocation failed.',err,1)

        !Generate the radial part of the CGTO
        call CMS_CGTO_on_grid(grid_r1_r2%r2,grid_r1_r2%n2_total_points,cgto_shell,cgto_on_grid)

        !Precalculate inverses of the r1 and r2 grid coordinates to replace the
        !division in the calculation of the radial part of the Legendre
        !resolution below.
        do k=1,grid_r1_r2%n1_total_points !r1
           inv_r1(k) = 1.0_cfp/grid_r1_r2%r1(k)
        enddo !k

        do k=1,grid_r1_r2%n2_total_points !r2
           inv_r2(k) = 1.0_cfp/grid_r1_r2%r2(k)
        enddo !k

        !Generate B-splines on the r2 grid multiplied by Jacobian and 1/r2 coming from the B-spline
        r2_B = 0.0_cfp
        do j=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
           B_start = grid_r1_r2%bspline_start_end_r2(1,j)
           B_end = grid_r1_r2%bspline_start_end_r2(2,j)
           r2_B(B_start:B_end,j) = grid_r1_r2%r2(B_start:B_end)*grid_r1_r2%B_vals_r2(B_start:B_end,j)
        enddo

        !$omp parallel do default(none) private(B_start, B_end, fac, j, k, l, k2, r12) &
        !$omp& firstprivate(f_l_w1_w2) shared(cgto_on_grid, grid_r1_r2, inv_r1, inv_r2, l_max, lambda_klj, r2_B)
        do k=1,grid_r1_r2%n1_total_points !r1

           !Calculate w1*w2*fac*r1*r<**l/r>**(l+1) from the Legendre expansion (l = 0)
           do k2=1,grid_r1_r2%n2_total_points !r2
              if (grid_r1_r2%r2(k2) .le. grid_r1_r2%r1(k)) then
                 f_l_w1_w2(k2) = grid_r1_r2%w1(k)*grid_r1_r2%w2(k2)
              else
                 r12 = grid_r1_r2%r1(k)*inv_r2(k2)
                 f_l_w1_w2(k2) = grid_r1_r2%w1(k)*grid_r1_r2%w2(k2)*r12
              endif
           enddo !k2

           do l=0,l_max

              !Generate the radial part of the Legendre expansion for the current l
              fac = fourpi/(2*l+1.0_cfp)

              !Integral over r2:
              do j=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                 B_start = grid_r1_r2%bspline_start_end_r2(1,j)
                 B_end = grid_r1_r2%bspline_start_end_r2(2,j)
                 lambda_klj(k,l,j) = fac*sum(f_l_w1_w2(B_start:B_end)*cgto_on_grid(B_start:B_end)*r2_B(B_start:B_end,j))
              enddo !j

              !Update w1*w2*fac*r1*r<**l/r>**(l+1) from the Legendre expansion (l > 0)
              do k2=1,grid_r1_r2%n2_total_points !r2
                 if (grid_r1_r2%r2(k2) .le. grid_r1_r2%r1(k)) then
                    r12 = grid_r1_r2%r2(k2)*inv_r1(k)
                 else
                    r12 = grid_r1_r2%r1(k)*inv_r2(k2)
                 endif
                 f_l_w1_w2(k2) = f_l_w1_w2(k2)*r12
              enddo !k2

           enddo !l
        enddo !k
        !$omp end parallel do

  end subroutine eval_lambda_BG_cms_G

  subroutine eval_BGBG_integrals_cms_G(this)
     use const_gbl, only: epsabs
     implicit none
     class(GG_shell_integrals_obj) :: this

     logical :: order_ab
     type(CGTO_shell_pw_expansion_obj) :: cgto_pw
     real(kind=cfp) :: threshold, cf
     integer :: lg,l_max,lpp_min,lpp_max,lbA,mbA,m,l,lbB,mbB,i,p,lbAmbA,lbBmbB,ind_cgto,ind_bto,jA,jB,B_start,B_end,&
                n_lbmb,k,nbto,ncgto,lpp,mpp,MgA_ind,MgB,MgB_ind,n_MgA,n_MgB
     real(kind=cfp), allocatable :: G_lm(:), r1_integral(:,:), lambda_klj(:,:,:)

        if (.not. this%initialized) then
           call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals_cms_G','Object not initialized.',1,1)
        endif

        if (.not.(this%A_is_cms .or. this%B_is_cms)) then
           call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals_cms_G','None of cgto_A, cgto_B sit on the CMS.',2,1)
        endif

        !order_ab = .true. <=> (BG_a|BG_b)
        !order_ab = .false. <=> (BG_b|BG_a)
        !The integrals must be always saved in the order (BG_a|BG_b) just like in eval_BGBG_integrals.

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs

        !Below the (BG_a|BG_b) integrals are evaluated so that G_b is the CMS CGTO
        !so we need to know how to map this order with the actual type of shells this%cgto_shell_A,this%cgto_shell_B.
        !The resulting integrals must be always ordered as (MgA,MgB) for the CGTO part, where MgA, MgB correspond to the shells this%cgto_shell_A,this%cgto_shell_B.
        if (this%B_is_cms) then
           order_ab = .true. !the required order of the CGTO shells and the order used in the evaluation are the same.
           lg = this%cgto_shell_B%l
           l_max = grid_r1_r2%max_bspline_l+lg !Maximum L in the Legendre expansion
           call eval_lambda_BG_cms_G(this%cgto_shell_B,l_max,lambda_klj)

           lpp_max = l_max+grid_r1_r2%max_bspline_l
           call cgto_pw%init_CGTO_shell_pw_expansion(this%cgto_shell_A,1)
           call cgto_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           call cgto_pw%eval_CGTO_single_projection_shell_pw_expansion(threshold,lpp_max)

           n_MgA = this%cgto_shell_A%number_of_functions
           n_MgB = this%cgto_shell_B%number_of_functions
        else
           order_ab = .false. !the required order of the CGTO shells and the order used in the evaluation are opposite.
           lg = this%cgto_shell_A%l
           l_max = grid_r1_r2%max_bspline_l+lg !Maximum L in the Legendre expansion
           call eval_lambda_BG_cms_G(this%cgto_shell_A,l_max,lambda_klj)

           lpp_max = l_max+grid_r1_r2%max_bspline_l
           call cgto_pw%init_CGTO_shell_pw_expansion(this%cgto_shell_B,1)
           call cgto_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           call cgto_pw%eval_CGTO_single_projection_shell_pw_expansion(threshold,lpp_max)

           n_MgA = this%cgto_shell_B%number_of_functions
           n_MgB = this%cgto_shell_A%number_of_functions
        endif

        !lg,MgB: angular momentum numbers of the CMS CGTO
        !n_MgA: number of M values for the non-CMS CGTO. MgA_ind: index of the M value of the non-CMS CGTO.
        !n_MgB: number of M values for the CMS CGTO. MgB, MgB_ind: M-values and indices of the CMS CGTO.

        nbto = (grid_r1_r2%max_bspline_l+1)**2*grid_r1_r2%last_bspline_inner
        nbto = nbto*nbto
        ncgto = this%cgto_shell_B%number_of_functions*this%cgto_shell_A%number_of_functions

        if (allocated(this%BGBG)) deallocate(this%BGBG)
        allocate(this%BGBG(nbto,ncgto),&
                 G_lm(grid_r1_r2%n1_total_points),&
                 r1_integral(grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner,&
                             grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner))

        this%BGBG = 0.0_cfp

        n_lbmb = (grid_r1_r2%max_bspline_l+1)**2

        !$omp parallel do schedule(dynamic, 1) default(none) &
        !$omp& private(MgA_ind, lbA, mbA, lbAmbA, lpp_min, lpp_max, l, m, G_lm, lpp, mpp, p, i, cf, k, jB, jA, B_start, B_end, &
        !$omp&         lbB, mbB, lbBmbB, MgB, MgB_ind, ind_cgto, ind_bto) &
        !$omp& firstprivate(r1_integral) &
        !$omp& shared(this, cgto_pw, l_max, grid_r1_r2, threshold, lambda_klj, order_ab, lg, n_MgA, n_MgB, n_lbmb)
        do MgA_ind=cgto_pw%cgto_shell%number_of_functions,1,-1

           do l=l_max,0,-1
              do lbA=0,grid_r1_r2%max_bspline_l
                 lpp_min = abs(l-lbA)
                 lpp_max = l+lbA
                 do mbA=-lbA,lbA
                    lbAmbA = lbA*lbA+lbA+mbA+1
                    do m=-l,l

                       !Construct the double angular projection of the CGTO
                       G_lm = 0.0_cfp
                       do lpp=lpp_min,lpp_max
                          do mpp=-lpp,lpp
                             p = lpp*lpp+mpp+lpp+1
                             i = cgto_pw%non_neg_indices_l(MgA_ind,p)
                             if (i > 0) then
                                cf = cpl%rgaunt(l,lbA,lpp,m,mbA,mpp)
                                if (cf /= 0.0_cfp) then
                                   do k=1,grid_r1_r2%n1_total_points
                                      G_lm(k) = G_lm(k) + cf*cgto_pw%angular_integrals(k,i)
                                   enddo !k
                                endif
                             endif
                          enddo !mpp
                       enddo !lpp

                       cf = maxval(abs(G_lm))
                       if (cf .le. threshold) cycle

                       !Compute the integral over r1: jB relates to the CMS CGTO B, jA relates to the CGTO A
                       do jB=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                          do jA=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                             B_start = grid_r1_r2%bspline_start_end_r1(1,jA)
                             B_end = grid_r1_r2%bspline_start_end_r1(2,jA)
                             r1_integral(jA,jB) = sum(G_lm(B_start:B_end)&
                                                      *grid_r1_r2%B_vals_r1(B_start:B_end,jA)&
                                                      *lambda_klj(B_start:B_end,l,jB))
                          enddo !jA
                       enddo !jB

                       do lbB=0,grid_r1_r2%max_bspline_l
                          do mbB=-lbB,lbB
                             lbBmbB = lbB*lbB+lbB+mbB+1
                             do MgB=-lg,lg
                                cf = cpl%rgaunt(l,lbB,lg,m,mbB,MgB)
                                if (cf /= 0.0_cfp) then
                                   MgB_ind = MgB+lg+1
                                   !compute the index and accumulate the final integral
                                   if (order_ab) then
                                      ind_cgto = MgA_ind + n_MgA*(MgB_ind-1) !save in order (MgA,MgB), where A,B are the actual cgto_shell_A, cgto_shell_B
                                      !save the BB part in order (jA,lAmA,jB,lBmB)
                                      do jB=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                                         do jA=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                                            ind_bto = jA + grid_r1_r2%last_bspline_inner*(lbAmbA-1)&
                                                         + grid_r1_r2%last_bspline_inner*n_lbmb*(jB-1)&
                                                         + grid_r1_r2%last_bspline_inner**2*n_lbmb*(lbBmbB-1)
                                            this%BGBG(ind_bto,ind_cgto) = this%BGBG(ind_bto,ind_cgto) + cf*r1_integral(jA,jB)
                                         enddo !jA
                                      enddo !jB
                                   else
                                      ind_cgto = MgB_ind + n_MgB*(MgA_ind-1) !save in order (MgA,MgB), where A,B are the actual cgto_shell_A, cgto_shell_B
                                      !save the BB part in order (jB,lBmB,jA,lAmA)
                                      do jB=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                                         do jA=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                                            ind_bto = jB + grid_r1_r2%last_bspline_inner*(lbBmbB-1)&
                                                         + grid_r1_r2%last_bspline_inner*n_lbmb*(jA-1)&
                                                         + grid_r1_r2%last_bspline_inner**2*n_lbmb*(lbAmbA-1)
                                            this%BGBG(ind_bto,ind_cgto) = this%BGBG(ind_bto,ind_cgto) + cf*r1_integral(jA,jB)
                                         enddo !jA
                                      enddo !jB
                                   endif
                                endif
                             enddo !MgB
                          enddo !mbB
                       enddo !lbB

                    enddo !m
                 enddo !mbA
              enddo !lbA
           enddo !l

        enddo !MgA_ind
        !$omp end parallel do

  end subroutine eval_BGBG_integrals_cms_G

  !> For efficiency the higher-L shell should be B.
  subroutine eval_BGBG_integrals(this,max_l_legendre)
     use const_gbl, only: epsabs
     use phys_const_gbl, only: fourpi
     use omp_lib
     implicit none
     class(GG_shell_integrals_obj) :: this
     integer, intent(in) :: max_l_legendre

     real(kind=cfp) :: threshold
     integer :: lpp, lpp_max, mpp, l,m,k,j1,j2,n_lbmb,n_m,p,p_min,p_max,k2,lb,mb,lbmb,p_offset,l_min,l_max,MgA,MgB,&
                i,ind,ind1,ind2,lb1mb1,lb2mb2,ind_AB,B_start,B_end,err,nbto,ncgto,n1,n2,lm
     integer :: p_start, p_end, iam, n_threads, i1,i2
     real(kind=cfp), allocatable :: couplings(:,:), G_lm(:,:,:), r1_integral(:), G_B_lm(:,:)
     integer(kind=1), allocatable :: G_lm_small(:,:)
     real(kind=cfp) :: fac, cf, last_contribution
     logical :: need_cf_B, need_cf_A
     type(CGTO_shell_pw_expansion_obj), save :: cgto_pw_B, cgto_pw_A
     real(kind=wp) :: t1,t2, tt1,tt2, startt,endt

        if (.not. this%initialized) call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals','Object not initialized.',1,1)

        write(level3,'("--------->GG_shell_integrals_obj:eval_BGBG_integrals")')

        startt = omp_get_wtime()

        !Compute the integrals easily if one of the CGTOs is sitting on CMS.
        if (this%A_is_cms .or. this%B_is_cms) then

           call this%eval_BGBG_integrals_cms_G
           endt = omp_get_wtime()

           write(level2,'("The Legendre expansion is an exact representation for this GG pair.")')
           write(level3,'("<---------GG_shell_integrals_obj:eval_BGBG_integrals and took [s]: ",f8.3)') endt-startt
           return
        endif

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs

        lpp_max = grid_r1_r2%max_bspline_l+grid_r1_r2%max_l_legendre

        if (allocated(cgto_pw_A%r_points)) then
           n2 = size(cgto_pw_A%r_points)
           if ((cgto_pw_A%cgto_shell_index /= this%shell_A_starting_index) .or. (n2 /= grid_r1_r2%n2_total_points)) then
              deallocate(cgto_pw_A%r_points)
           endif
        endif

        if (.not. allocated(cgto_pw_A%r_points)) then
           cgto_pw_A%cgto_shell_index = this%shell_A_starting_index
           write(level2,'(/,"Evaluating the Y_lm function for functions with starting index: ",i4)') this%shell_A_starting_index
           call cgto_pw_A%init_CGTO_shell_pw_expansion(this%cgto_shell_A,this%shell_A_starting_index)
           call cgto_pw_A%eval_BTO_CGTO_Y_lm(grid_r1_r2)
        endif

        if (allocated(cgto_pw_B%r_points)) then
           n1 = size(cgto_pw_B%r_points)
           if ((cgto_pw_B%cgto_shell_index /= this%shell_B_starting_index) .or. (n1 /= grid_r1_r2%n1_total_points)) then
              deallocate(cgto_pw_B%r_points)
           endif
        endif

        if (.not. allocated(cgto_pw_B%r_points)) then
           !The CGTO shell B needs projections only on the r1 grid.
           write(level2,'("Evaluating pw expansion for shell with starting index: ",i4)') this%shell_B_starting_index
           call cgto_pw_B%init_CGTO_shell_pw_expansion(this%cgto_shell_B,this%shell_B_starting_index)
           call cgto_pw_B%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           call cgto_pw_B%eval_CGTO_single_projection_shell_pw_expansion(threshold,lpp_max)
           cgto_pw_B%cgto_shell_index = this%shell_B_starting_index
        endif

        n_lbmb = (grid_r1_r2%max_bspline_l+1)**2
        n_m = 2*max_l_legendre+1

        i = 0
        do j1=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
           B_start = grid_r1_r2%bspline_start_end_r1(1,j1)
           B_end = grid_r1_r2%bspline_start_end_r1(2,j1)
           i = max(i,B_end-B_start+1)
        enddo

        nbto = (grid_r1_r2%max_bspline_l+1)**2*grid_r1_r2%last_bspline_inner
        nbto = nbto*nbto
        ncgto = cgto_pw_B%cgto_shell%number_of_functions*cgto_pw_A%cgto_shell%number_of_functions

        if (allocated(this%BGBG)) deallocate(this%BGBG)
        allocate(G_lm(grid_r1_r2%n1_total_points,n_lbmb*n_m,cgto_pw_B%cgto_shell%number_of_functions),this%BGBG(nbto,ncgto),&
                 G_lm_small(n_lbmb*n_m,cgto_pw_B%cgto_shell%number_of_functions),&
                 stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals','Memory allocation failed.',err,1)

        this%BGBG = 0.0_cfp

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(MgA,j1,B_start,B_end,i,ind,MgB,ind_AB,lb1mb1,j2,lb2mb2,ind2,m,r1_integral,&
        !$OMP & iam,n_threads,G_B_lm,fac,n_m,k,lb,mb,lpp,mpp,p,l_min,l_max,lbmb,cf,err,t1,t2,tt1,tt2,l,lm) &
        !$OMP & SHARED(cgto_pw_B,cgto_pw_A,G_lm,this,n_lbmb,grid_r1_r2,&
        !$OMP &        threshold,G_lm_small,last_contribution,i1,i2)
        i = 0
        i1 = 0
        i2 = 0
        do j2=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
           B_start = grid_r1_r2%bspline_start_end_r1(1,j2)
           B_end = grid_r1_r2%bspline_start_end_r1(2,j2)
           i = max(i,B_end-B_start+1)
        enddo
        n_m = 2*max_l_legendre+1
        allocate(G_B_lm(i,n_lbmb*n_m),r1_integral(n_m),stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals','Memory allocation failed.',err,1)

        !Precalculate the Gaunt coefficients
        l_max = max_l_legendre+grid_r1_r2%max_bspline_l
        call cpl%prec_cgaunt(l_max)

        last_contribution = 0.0_cfp

        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()

        do l=max_l_legendre,0,-1

           tt1 = omp_get_wtime()
           n_m = 2*l+1

           !Construct the double angular projections for the CGTO B pw projections
           t1 = omp_get_wtime()
           !$OMP SINGLE
           G_lm = 0.0_cfp
           !$OMP END SINGLE
           do MgB=1,cgto_pw_B%cgto_shell%number_of_functions
              do lb=0,grid_r1_r2%max_bspline_l
                 i = lb+1 + (grid_r1_r2%max_bspline_l+1)*(MgB-1)
                 if (mod(i,n_threads) /= iam) cycle !work redistribution
                 l_min = abs(l-lb)
                 l_max = l+lb
                 do mb=-lb,lb
                    lbmb = lb*lb+lb+mb+1
                    do m=-l,l
                       ind = m+l+1 + n_m*(lbmb-1)

                       do lpp=l_min,l_max !todo selection rules
                          do mpp=-lpp,lpp !todo selection rules
                             p = lpp*lpp+mpp+lpp+1
                             i = cgto_pw_B%non_neg_indices_l(MgB,p)
                             if (i > 0) then
                                cf = cpl%rgaunt(l,lb,lpp,m,mb,mpp)
                                if (cf /= 0.0_cfp) then
                                   do k=1,grid_r1_r2%n1_total_points
                                      G_lm(k,ind,MgB) = G_lm(k,ind,MgB) + cf*cgto_pw_B%angular_integrals(k,i)
                                   enddo !k
                                endif
                             endif
                          enddo !mpp
                       enddo !lpp

                       cf = maxval(abs(G_lm(1:grid_r1_r2%n1_total_points,ind,MgB)))
                       if (cf .le. threshold) then
                          G_lm_small(ind,MgB) = 1
                       else
                          G_lm_small(ind,MgB) = 0
                       endif

                    enddo !m
                 enddo !mb
              enddo !lb
           enddo !MgB
           !$OMP BARRIER
           !$OMP SINGLE
           t2 = omp_get_wtime()
           !print *,'part4',t2-t1
           !print *,'parts1-4',t2-tt1
           !$OMP END SINGLE

           !todo for Ga == Gb we need only the unique quartets of indices!!!
           t1 = omp_get_wtime()
           do MgB=1,cgto_pw_B%cgto_shell%number_of_functions
              do j2=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                 i = j2 + grid_r1_r2%last_bspline_inner*(MgB-1)
                 if (mod(i,n_threads) /= iam) cycle !work redistribution
                 B_start = grid_r1_r2%bspline_start_end_r1(1,j2)
                 B_end = grid_r1_r2%bspline_start_end_r1(2,j2)
                 i = B_end-B_start+1
                 do ind=1,n_lbmb*n_m
                    G_B_lm(1:i,ind) = G_lm(B_start:B_end,ind,MgB)*grid_r1_r2%B_vals_r1(B_start:B_end,j2)
                 enddo !ind
                 do MgA=1,cgto_pw_A%cgto_shell%number_of_functions
                    ind_AB = MgA + cgto_pw_A%cgto_shell%number_of_functions*(MgB-1) !save in order(MgA,MgB)
                    do j1=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                       do lb1mb1=1,n_lbmb
                          do lb2mb2=1,n_lbmb
                             !save the BB part in order (j1,l1m1,j2,l2m2)
                             ind = j1 + grid_r1_r2%last_bspline_inner*(lb1mb1-1)&
                                      + grid_r1_r2%last_bspline_inner*n_lbmb*(j2-1)&
                                      + grid_r1_r2%last_bspline_inner**2*n_lbmb*(lb2mb2-1)
                             ind2 = l+1 + n_m*(lb2mb2-1)
                             do m=-l,l
                                lm = l*l + l+m+1
                                p = CGTO_pw_A%Y_lm_non_neg_indices(MgA,lb1mb1,lm)
                                if (p == 0 .or. G_lm_small(ind2+m,MgB) == 1) then
                                   r1_integral(m+l+1) = 0.0_cfp
                                else
                                   r1_integral(m+l+1) = sum(G_B_lm(1:i,ind2+m)*CGTO_pw_A%Y_lm_mixed(B_start:B_end,p,j1))
                                endif
                             enddo !m
                             this%BGBG(ind,ind_AB) = this%BGBG(ind,ind_AB) + sum(r1_integral(1:n_m)) 
                          enddo !lb2mb2
                       enddo !j1
                    enddo !lb1mb1
                 enddo !MgA
              enddo !j2
           enddo !MgB
           !$OMP BARRIER
           !$OMP SINGLE
           t2 = omp_get_wtime()
           !print *,'part5',t2-t1
           !$OMP END SINGLE

           !$OMP SINGLE
           if (last_contribution .eq. 0.0_cfp) then
              !Find the largest contribution to any of the integrals. Since
              !we're looping over l in the inverse order this value corresponds
              !to the largest from the smallest contributions to the integrals.
              cf = 0.0_cfp
              do ind_AB=1,size(this%BGBG,2)
                 do ind=1,size(this%BGBG,1)
                    last_contribution = max(last_contribution,abs(this%BGBG(ind,ind_AB)))
                    if (last_contribution /= cf) then
                       i1 = ind
                       i2 = ind_AB
                       cf = last_contribution
                    endif
                 enddo
              enddo
           endif
           tt2 = omp_get_wtime()
           !print *,'total',tt2-tt1
           !$OMP END SINGLE
        enddo !l

        deallocate(G_B_lm)

        !$OMP END PARALLEL
        !
        !---- Analysis of the convergence
        !
        if (i1 > 0 .and. i2 > 0) then
           cf = abs(this%BGBG(i1,i2))
           fac = last_contribution/cf
           write(level2,'("Convergence (rel. prec.) of the Legendre expansion for the largest integral from the batch: "&
                         &,e8.2,1X,2e25.15)') fac, last_contribution, cf
        endif

        endt = omp_get_wtime()

        write(level3,'("<---------GG_shell_integrals_obj:eval_BGBG_integrals and took [s]: ",f8.3)') endt-startt

  end subroutine eval_BGBG_integrals

  !> For efficiency the higher-L shell should be B.
  subroutine eval_BGBG_integrals_direct(this,max_l_legendre)
     use const_gbl, only: epsabs
     use phys_const_gbl, only: fourpi
     use omp_lib
     implicit none
     class(GG_shell_integrals_obj) :: this
     integer, intent(in) :: max_l_legendre

     real(kind=cfp) :: threshold
     integer :: lpp, lpp_max, mpp, l,m,k,j1,j2,n_lbmb,n_m,p,p_min,p_max,k2,lb,mb,lbmb,p_offset,l_min,l_max,MgA,MgB,&
                i,ind,ind1,ind2,lb1mb1,lb2mb2,ind_AB,B_start,B_end,err,nbto,ncgto,n1,n2
     integer :: p_start, p_end, iam, n_threads, i1,i2
     real(kind=cfp), allocatable :: couplings(:,:), lambda_lpp_mpp(:,:,:,:), f_l_w1_w2(:,:), lambda_lm(:,:,:,:), tmp(:,:),&
                                    G_lm(:,:,:), r1_integral(:), inv_r1(:), inv_r2(:), bto_leg_part(:), G_B_lm(:,:)
     integer(kind=1), allocatable :: lambda_lm_small(:,:,:), G_lm_small(:,:)
     real(kind=cfp), allocatable :: r2_B(:,:)
     real(kind=cfp) :: r12, fac, cf, last_contribution
     logical :: need_cf_B, need_cf_A
     type(CGTO_shell_pw_expansion_obj), save :: cgto_pw_B, cgto_pw_A
     real(kind=wp) :: t1,t2, tt1,tt2, startt,endt

        if (.not. this%initialized) call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals_direct',&
           'Object not initialized.',1,1)

        write(stdout,'("--------->GG_shell_integrals_obj:eval_BGBG_integrals_direct")')

        startt = omp_get_wtime()

        !Compute the integrals easily if one of the CGTOs is sitting on CMS.
        if (this%A_is_cms .or. this%B_is_cms) then

           call this%eval_BGBG_integrals_cms_G
           endt = omp_get_wtime()

           write(stdout,'("The Legendre expansion is an exact representation for this GG pair.")')
           write(stdout,'("<---------GG_shell_integrals_obj:eval_BGBG_integrals_direct and took [s]: ",f8.3)') endt-startt
           return
        endif

        !The value below is the smallest GTO amplitude which is deemed significant.
        !We neglect integrals smaller than epsabs. Consider an integral with absolute value ~epsabs. 
        !If the required number of significant digits is N then we need to consider contributions of values to this integral not smaller than epsabs*10**(-N).
        threshold = 10**(-precision(cfp_dummy)+1.0_cfp) !we don't need the full relative precision since the quadrature rules don't give better relative precision than ~10e-10
        threshold = threshold*epsabs

        lpp_max = grid_r1_r2%max_bspline_l+max_l_legendre

        if (allocated(cgto_pw_A%r_points)) then
           n2 = size(cgto_pw_A%r_points)
           if ((cgto_pw_A%cgto_shell_index /= this%shell_A_starting_index) .or. (n2 /= grid_r1_r2%n2_total_points)) then
              deallocate(cgto_pw_A%r_points)
           endif
        endif

        if (.not. allocated(cgto_pw_A%r_points)) then
           !The CGTO shell B will be used to evaluate the lambda function so we need the projections for the r1 grid.
           write(stdout,'("Evaluating pw expansion for shell with starting index: ",i4)') this%shell_A_starting_index
           call cgto_pw_A%init_CGTO_shell_pw_expansion(this%cgto_shell_A,this%shell_A_starting_index)
           call cgto_pw_A%assign_grid(grid_r1_r2%r2,grid_r1_r2%w2)
           call cgto_pw_A%eval_CGTO_single_projection_shell_pw_expansion(threshold,lpp_max)
           cgto_pw_A%cgto_shell_index = this%shell_A_starting_index
        endif

        if (allocated(cgto_pw_B%r_points)) then
           n1 = size(cgto_pw_B%r_points)
           if ((cgto_pw_B%cgto_shell_index /= this%shell_B_starting_index) .or. (n1 /= grid_r1_r2%n1_total_points)) then
              deallocate(cgto_pw_B%r_points)
           endif
        endif

        if (.not. allocated(cgto_pw_B%r_points)) then
           !The CGTO shell A needs projections only on the r1 grid.
           write(stdout,'("Evaluating pw expansion for shell with starting index: ",i4)') this%shell_B_starting_index
           call cgto_pw_B%init_CGTO_shell_pw_expansion(this%cgto_shell_B,this%shell_B_starting_index)
           call cgto_pw_B%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
           call cgto_pw_B%eval_CGTO_single_projection_shell_pw_expansion(threshold,lpp_max)
           cgto_pw_B%cgto_shell_index = this%shell_B_starting_index
        endif

        n_lbmb = (grid_r1_r2%max_bspline_l+1)**2
        n_m = 2*max_l_legendre+1

        !maximum range of lppmpp indices needed from the angular projections
        p_min = min((max_l_legendre-grid_r1_r2%max_bspline_l)**2,max_l_legendre**2)
        p_max = (max_l_legendre+grid_r1_r2%max_bspline_l+1)**2
        p_offset = p_max-p_min

        i = 0
        do j1=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
           B_start = grid_r1_r2%bspline_start_end_r1(1,j1)
           B_end = grid_r1_r2%bspline_start_end_r1(2,j1)
           i = max(i,B_end-B_start+1)
        enddo

        nbto = (grid_r1_r2%max_bspline_l+1)**2*grid_r1_r2%last_bspline_inner
        nbto = nbto*nbto
        ncgto = cgto_pw_B%cgto_shell%number_of_functions*cgto_pw_A%cgto_shell%number_of_functions

        if (allocated(this%BGBG)) deallocate(this%BGBG)
        allocate(f_l_w1_w2(grid_r1_r2%n2_total_points,grid_r1_r2%n1_total_points),&
                 tmp(grid_r1_r2%n1_total_points,p_offset),&
                 inv_r1(grid_r1_r2%n1_total_points),&
                 inv_r2(grid_r1_r2%n2_total_points),r1_integral(n_m),couplings(p_offset,n_lbmb*n_m),&
                 lambda_lpp_mpp(p_offset,grid_r1_r2%n1_total_points,&
                                grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner,&
                                cgto_pw_A%cgto_shell%number_of_functions),&
                 lambda_lm(grid_r1_r2%n1_total_points,n_lbmb*n_m,&
                           grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner,&
                           cgto_pw_A%cgto_shell%number_of_functions),&
                 lambda_lm_small(n_lbmb*n_m,grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner,&
                                 cgto_pw_A%cgto_shell%number_of_functions),&
                 G_lm(grid_r1_r2%n1_total_points,n_lbmb*n_m,cgto_pw_B%cgto_shell%number_of_functions),this%BGBG(nbto,ncgto),&
                 G_lm_small(n_lbmb*n_m,cgto_pw_B%cgto_shell%number_of_functions),&
                 r2_B(grid_r1_r2%n2_total_points,grid_r1_r2%first_bspline_index:grid_r1_r2%last_bspline_inner),&
                 stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals_direct','Memory allocation failed.',err,1)

        !Precalculate inverses of the r1 and r2 grid coordinates to replace the
        !division in the calculation of the radial part of the Legendre
        !resolution below.
        do k=1,grid_r1_r2%n1_total_points !r1
           inv_r1(k) = 1.0_cfp/grid_r1_r2%r1(k)
        enddo !k

        do k=1,grid_r1_r2%n2_total_points !r2
           inv_r2(k) = 1.0_cfp/grid_r1_r2%r2(k)
        enddo !k

        r2_B = 0.0_cfp
        do j1=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
           B_start = grid_r1_r2%bspline_start_end_r2(1,j1)
           B_end = grid_r1_r2%bspline_start_end_r2(2,j1)
           r2_B(B_start:B_end,j1) = grid_r1_r2%r2(B_start:B_end)*grid_r1_r2%B_vals_r2(B_start:B_end,j1)
        enddo

        this%BGBG = 0.0_cfp

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(MgA,j1,B_start,B_end,i,ind,MgB,ind_AB,lb1mb1,j2,lb2mb2,ind1,ind2,m,r1_integral,&
        !$OMP & iam,n_threads,bto_leg_part,G_B_lm,fac,p_min,p_max,n_m,r12,&
        !$OMP & k,k2,lb,mb,lpp,mpp,p_offset,p,l_min,l_max,lbmb,cf,err,t1,t2,tt1,tt2,l) &
        !$OMP & SHARED(cgto_pw_B,cgto_pw_A,G_lm,this,n_lbmb,lambda_lm,f_l_w1_w2,lambda_lpp_mpp,inv_r1,inv_r2,r2_B,grid_r1_r2,&
        !$OMP &        lambda_lm_small,threshold,G_lm_small,max_l_legendre,last_contribution,i1,i2)
        i = 0
        i1 = 0
        i2 = 0
        do j2=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
           B_start = grid_r1_r2%bspline_start_end_r1(1,j2)
           B_end = grid_r1_r2%bspline_start_end_r1(2,j2)
           i = max(i,B_end-B_start+1)
        enddo
        n_m = 2*max_l_legendre+1
        allocate(bto_leg_part(grid_r1_r2%n2_total_points),G_B_lm(i,n_lbmb*n_m),stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_BGBG_integrals_direct','Memory allocation failed.',err,1)

        !Precalculate the Gaunt coefficients
        l_max = max_l_legendre+grid_r1_r2%max_bspline_l
        call cpl%prec_cgaunt(l_max)

        last_contribution = 0.0_cfp

        n_threads = omp_get_num_threads()
        iam = omp_get_thread_num()

        do l=max_l_legendre,0,-1

           tt1 = omp_get_wtime()
           n_m = 2*l+1

           !range of lppmpp indices needed from the angular projections: I need all lpp = abs(l-lb),...,l+lb where lb is the BTO
           !angular momentum: lb = 0,...,grid_r1_r2%max_bspline_l
           if (l .le. grid_r1_r2%max_bspline_l) then
              p_min = 0
           else
              p_min = (l-grid_r1_r2%max_bspline_l  )**2
           endif
           p_max = (l+grid_r1_r2%max_bspline_l+1)**2
           
           !Generate the radial part of the Legendre expansion for the current l
           t1 = omp_get_wtime()
           fac = fourpi/(2*l+1.0_cfp)
           !$OMP DO
           do k=1,grid_r1_r2%n1_total_points !r1
              do k2=1,grid_r1_r2%n2_total_points !r2
                 !Calculate w1*w2*fac*r1*r<**l/r>**(l+1) from the Legendre expansion
                 if (grid_r1_r2%r2(k2) .le. grid_r1_r2%r1(k)) then
                    r12 = grid_r1_r2%r2(k2)*inv_r1(k)
                    f_l_w1_w2(k2,k) = grid_r1_r2%w1(k)*grid_r1_r2%w2(k2)*fac*(r12)**l
                 else
                    r12 = grid_r1_r2%r1(k)*inv_r2(k2)
                    f_l_w1_w2(k2,k) = grid_r1_r2%w1(k)*grid_r1_r2%w2(k2)*fac*(r12)**(l+1)
                 endif
              enddo !k2
           enddo !k
           !$OMP END DO
           !$OMP SINGLE
           t2 = omp_get_wtime()
           !print *,'l',l
           !print *,'part1',t2-t1
           !$OMP END SINGLE

           !Generate lambda_lpp_mpp integrating over r2 all lppmpp projections with the radial part of the Legendre resolution and all radial B-splines.
           t1 = omp_get_wtime()
           !$OMP DO
           do j1=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
              !todo the sum over B_start:B_end takes the most time below so if the r2 grid is optimized the whole calculation can be much faster!
              B_start = grid_r1_r2%bspline_start_end_r2(1,j1)
              B_end = grid_r1_r2%bspline_start_end_r2(2,j1)
              do k=1,grid_r1_r2%n1_total_points
                 bto_leg_part(B_start:B_end) = f_l_w1_w2(B_start:B_end,k)*r2_B(B_start:B_end,j1)
                 do p=p_min+1,p_max
                    p_offset = p-p_min
                    do MgA=1,cgto_pw_A%cgto_shell%number_of_functions
                       i = cgto_pw_A%non_neg_indices_l(MgA,p)
                       if (i .eq. 0) then
                          lambda_lpp_mpp(p_offset,k,j1,MgA) = 0.0_cfp
                       else
                          lambda_lpp_mpp(p_offset,k,j1,MgA) = sum(bto_leg_part(B_start:B_end)&
                                                                  *cgto_pw_A%angular_integrals(B_start:B_end,i))
                       endif
                    enddo !cgto_m
                 enddo !p
              enddo !k
           enddo !j1
           !$OMP END DO
           !$OMP SINGLE
           t2 = omp_get_wtime()
           !print *,'part2',t2-t1
           !$OMP END SINGLE

           !Construct the double angular projections for the lambda function
           t1 = omp_get_wtime()
           !$OMP SINGLE
           lambda_lm = 0.0_cfp
           !$OMP END SINGLE
           do MgA=1,cgto_pw_A%cgto_shell%number_of_functions
              do j1=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                 i = j1 + grid_r1_r2%last_bspline_inner*(MgA-1)
                 if (mod(i,n_threads) /= iam) cycle !work redistribution
                 !todo transpose here the part lambda_lpp_mpp(:,:,j1,MgA) using abcd_to_cdab
                 do lb=0,grid_r1_r2%max_bspline_l
                    l_min = abs(l-lb)
                    l_max = l+lb
                    do mb=-lb,lb
                       lbmb = lb*lb+lb+mb+1
                       do m=-l,l
                          ind = m+l+1 + n_m*(lbmb-1)

                          do lpp=l_min,l_max,2
                             do mpp=-lpp,lpp !todo selection rules
                                cf = cpl%rgaunt(l,lb,lpp,m,mb,mpp)
                                if (cf /= 0.0_cfp) then
                                   p = lpp*lpp+mpp+lpp+1
                                   p_offset = p-p_min
                                   !todo change order of columns of lambda_lpp_mpp so k is first or transpose it after it has been determined above?
                                   !todo I can probably use cgto_pw_A%non_neg_indices_l(MgA,p) to figure out if lambda_lpp_mpp is significant!!!
                                   if (cgto_pw_A%non_neg_indices_l(MgA,p) > 0) then
                                      do k=1,grid_r1_r2%n1_total_points
                                          lambda_lm(k,ind,j1,MgA) = lambda_lm(k,ind,j1,MgA) + cf*lambda_lpp_mpp(p_offset,k,j1,MgA)
                                      enddo !k
                                   endif
                                endif
                             enddo !mpp
                          enddo !lpp

                          cf = maxval(abs(lambda_lm(1:grid_r1_r2%n1_total_points,ind,j1,MgA)))
                          if (cf .le. threshold) then
                             lambda_lm_small(ind,j1,MgA) = 1
                          else
                             lambda_lm_small(ind,j1,MgA) = 0
                          endif

                       enddo !m
                    enddo !mb
                 enddo !lb
              enddo !j1
           enddo !MgA
           !$OMP BARRIER
           !$OMP SINGLE
           t2 = omp_get_wtime()
           !print *,'part3',t2-t1
           !$OMP END SINGLE

           !Construct the double angular projections for the CGTO A pw projections
           t1 = omp_get_wtime()
           !$OMP SINGLE
           G_lm = 0.0_cfp
           !$OMP END SINGLE
           do MgB=1,cgto_pw_B%cgto_shell%number_of_functions
              do lb=0,grid_r1_r2%max_bspline_l
                 i = lb+1 + (grid_r1_r2%max_bspline_l+1)*(MgB-1)
                 if (mod(i,n_threads) /= iam) cycle !work redistribution
                 l_min = abs(l-lb)
                 l_max = l+lb
                 do mb=-lb,lb
                    lbmb = lb*lb+lb+mb+1
                    do m=-l,l
                       ind = m+l+1 + n_m*(lbmb-1)

                       do lpp=l_min,l_max !todo selection rules
                          do mpp=-lpp,lpp !todo selection rules
                             p = lpp*lpp+mpp+lpp+1
                             i = cgto_pw_B%non_neg_indices_l(MgB,p)
                             if (i > 0) then
                                cf = cpl%rgaunt(l,lb,lpp,m,mb,mpp)
                                if (cf /= 0.0_cfp) then
                                   do k=1,grid_r1_r2%n1_total_points
                                      G_lm(k,ind,MgB) = G_lm(k,ind,MgB) + cf*cgto_pw_B%angular_integrals(k,i)
                                   enddo !k
                                endif
                             endif
                          enddo !mpp
                       enddo !lpp

                       cf = maxval(abs(G_lm(1:grid_r1_r2%n1_total_points,ind,MgB)))
                       if (cf .le. threshold) then
                          G_lm_small(ind,MgB) = 1
                       else
                          G_lm_small(ind,MgB) = 0
                       endif

                    enddo !m
                 enddo !mb
              enddo !lb
           enddo !MgB
           !$OMP BARRIER
           !$OMP SINGLE
           t2 = omp_get_wtime()
           !print *,'part4',t2-t1
           !print *,'parts1-4',t2-tt1
           !$OMP END SINGLE

           !todo for Ga == Gb we need only the unique quartets of indices!!!
           t1 = omp_get_wtime()
           do MgB=1,cgto_pw_B%cgto_shell%number_of_functions
              do j2=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                 i = j2 + grid_r1_r2%last_bspline_inner*(MgB-1)
                 if (mod(i,n_threads) /= iam) cycle !work redistribution
                 B_start = grid_r1_r2%bspline_start_end_r1(1,j2)
                 B_end = grid_r1_r2%bspline_start_end_r1(2,j2)
                 i = B_end-B_start+1
                 do ind=1,n_lbmb*n_m
                    G_B_lm(1:i,ind) = G_lm(B_start:B_end,ind,MgB)*grid_r1_r2%B_vals_r1(B_start:B_end,j2)
                 enddo !ind
                 do MgA=1,cgto_pw_A%cgto_shell%number_of_functions
                    ind_AB = MgA + cgto_pw_A%cgto_shell%number_of_functions*(MgB-1) !save in order(MgA,MgB)
                    do j1=grid_r1_r2%first_bspline_index,grid_r1_r2%last_bspline_inner
                       do lb1mb1=1,n_lbmb
                          do lb2mb2=1,n_lbmb
                             !save the BB part in order (j1,l1m1,j2,l2m2)
                             ind = j1 + grid_r1_r2%last_bspline_inner*(lb1mb1-1)&
                                      + grid_r1_r2%last_bspline_inner*n_lbmb*(j2-1)&
                                      + grid_r1_r2%last_bspline_inner**2*n_lbmb*(lb2mb2-1)
                             ind2 = l+1 + n_m*(lb2mb2-1)
                             ind1 = l+1 + n_m*(lb1mb1-1)
                             do m=-l,l
                                if ((lambda_lm_small(ind1+m,j1,MgA) .eq. 1) .or. (G_lm_small(ind2+m,MgB) .eq. 1)) then
                                   r1_integral(m+l+1) = 0.0_cfp
                                else
                                   r1_integral(m+l+1) = sum(G_B_lm(1:i,ind2+m)*lambda_lm(B_start:B_end,ind1+m,j1,MgA))
                                endif
                             enddo !m
                             this%BGBG(ind,ind_AB) = this%BGBG(ind,ind_AB) + sum(r1_integral(1:n_m)) 
                          enddo !lbm2
                       enddo !j1
                    enddo !lb1mb1
                 enddo !MgA
              enddo !j2
           enddo !MgB
           !$OMP BARRIER
           !$OMP SINGLE
           t2 = omp_get_wtime()
           !print *,'part5',t2-t1
           !$OMP END SINGLE

           !$OMP SINGLE
           if (last_contribution .eq. 0.0_cfp) then
              !Find the largest contribution to any of the integrals. Since
              !we're looping over l in the inverse order this value corresponds
              !to the largest from the smallest contributions to the integrals.
              cf = 0.0_cfp
              do ind_AB=1,size(this%BGBG,2)
                 do ind=1,size(this%BGBG,1)
                    last_contribution = max(last_contribution,abs(this%BGBG(ind,ind_AB)))
                    if (last_contribution /= cf) then
                       i1 = ind
                       i2 = ind_AB
                       cf = last_contribution
                    endif
                 enddo
              enddo
           endif
           tt2 = omp_get_wtime()
           !print *,'total',tt2-tt1
           !$OMP END SINGLE
        enddo !l

        deallocate(bto_leg_part,G_B_lm)

        !$OMP END PARALLEL
        !
        !---- Analysis of the convergence
        !
        if (i1 > 0 .and. i2 > 0) then
           cf = abs(this%BGBG(i1,i2))
           fac = last_contribution/cf
           write(stdout,'("Convergence (rel. prec.) of the Legendre expansion for the largest integral from the batch: "&
                         &,e8.2,1X,2e25.15)') fac, last_contribution, cf
        endif

        endt = omp_get_wtime()

        write(stdout,'("<---------GG_shell_integrals_obj:eval_BGBG_integrals_direct and took [s]: ",f8.3)') endt-startt

  end subroutine eval_BGBG_integrals_direct

  subroutine eval_GG_CCTT_prototype_integrals(this,shell_1,shell_2)
     use cgto_pw_expansions_gbl, only: CGTO_shell_pair_pw_expansion_obj
     implicit none
     class(GG_shell_integrals_obj) :: this !=TT pair
     type(CGTO_shell_data_obj), intent(in) :: shell_1, shell_2 !=CC pair

     type(CGTO_shell_pair_pw_expansion_obj) :: cgto_pair_pw
     integer :: max_l

        if (.not. this%initialized) then
           call xermsg('GG_shell_mixed_integrals_mod','eval_GG_CCTT_prototype_integrals','Object not initialized.',1,1)
        endif

        call xermsg('GG_shell_mixed_integrals_mod','eval_GG_CCTT_prototype_integrals','Not implemented yet.',1,1)

        call cgto_pair_pw%init_CGTO_shell_pair_pw_expansion(this%cgto_shell_A,&
                                                            this%shell_A_starting_index,&
                                                            this%cgto_shell_B,&
                                                            this%shell_B_starting_index)

        !todo eval_regular_grid or use some other one?
        !call cgto_pair_pw%eval_regular_grid(0.0_cfp,rmat_radius,delta_r)

        max_l = shell_1%l + shell_2%l
        call cgto_pair_pw%eval_CGTO_shell_pair_pw_expansion(max_l)

        !todo implement:
        !call cgto_pair_pw%eval_radial_GG_CG_CG(grid_r1_r2) !result in cgto_pair_pw%radial_lm_BB_GG(lm,BA_ind,pair_index), BA_ind = Mg_B+cgto_B%l+1 + n_cgto_B_m*(Mg_A_ind-1)

  end subroutine eval_GG_CCTT_prototype_integrals

  subroutine eval_BG_CCTT_prototype_integrals(this,bspline_grid,CGTO_shell)
     use bspline_grid_gbl
     implicit none
     class(GG_shell_integrals_obj) :: this !=TT pair
     type(bspline_grid_obj), intent(in) :: bspline_grid !=C
     type(CGTO_shell_data_obj), intent(in) :: CGTO_shell !=C

        if (.not. this%initialized) then
           call xermsg('GG_shell_mixed_integrals_mod','eval_BG_CCTT_prototype_integrals','Object not initialized.',1,1)
        endif

        call xermsg('GG_shell_mixed_integrals_mod','eval_BG_CCTT_prototype_integrals','Not implemented yet.',1,1)

  end subroutine eval_BG_CCTT_prototype_integrals

  subroutine eval_BB_CCTT_prototype_integrals(this,grid_r1_r2)
     use bspline_grid_gbl
     use cgto_pw_expansions_gbl, only: CGTO_shell_pair_pw_expansion_obj
     use grid_gbl, only: grid_r1_r2_obj
     implicit none
     class(GG_shell_integrals_obj) :: this !=TT pair
     class(grid_r1_r2_obj), intent(in) :: grid_r1_r2 !grid_r1_r2%bspline_grid =CC

     type(CGTO_shell_pair_pw_expansion_obj) :: cgto_pair_pw
     integer :: max_l

        if (.not. this%initialized) then
           call xermsg('GG_shell_mixed_integrals_mod','eval_BB_CCTT_prototype_integrals','Object not initialized.',1,1)
        endif

        call cgto_pair_pw%init_CGTO_shell_pair_pw_expansion(this%cgto_shell_A,&
                                                            this%shell_A_starting_index,&
                                                            this%cgto_shell_B,&
                                                            this%shell_B_starting_index)

        call cgto_pair_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)

        max_l = 2*grid_r1_r2%max_bspline_l
        call cgto_pair_pw%eval_CGTO_shell_pair_pw_expansion(max_l)

        !todo the prototypes can be simplified even more for BB pairs from the outer region: uncoupled integrations over r1 and r2 with the corresponding uncoupled storage.
        call cgto_pair_pw%eval_radial_GG_BB(grid_r1_r2) !result in cgto_pair_pw%radial_lm_BB_GG(lm,BA_ind,pair_index), BA_ind = Mg_B+cgto_B%l+1 + n_cgto_B_m*(Mg_A_ind-1)

  end subroutine eval_BB_CCTT_prototype_integrals

  subroutine eval_BGGG_integrals(this,cgto_shell_C,shell_C_starting_index)
     use const_gbl, only: epsabs
     implicit none
     class(GG_shell_integrals_obj) :: this
     type(CGTO_shell_data_obj) :: cgto_shell_C
     integer, intent(in) :: shell_C_starting_index

     integer :: lena,lenb,la,lb,err,BG_size,GG_size,n_projections,max_l
     real(kind=cfp) :: xa,ya,za, xb,yb,zb, acnorm,bcnorm,threshold, RC, tol
     logical :: same_ab_shells

        if (.not. this%initialized) call xermsg('GG_shell_mixed_integrals_mod','eval_BGGG_integrals','Object not initialized.',1,1)

        this%cgto_shell = cgto_shell_C
        this%shell_starting_index = shell_C_starting_index

        tol = F1MACH(4,cfp_dummy)
        RC = sqrt(dot_product(cgto_shell_C%center,cgto_shell_C%center))

        xa = this%cgto_shell_A%center(1)
        ya = this%cgto_shell_A%center(2)
        za = this%cgto_shell_A%center(3)

        xb = this%cgto_shell_B%center(1)
        yb = this%cgto_shell_B%center(2)
        zb = this%cgto_shell_B%center(3)

        GG_size = this%cgto_shell_A%number_of_functions*this%cgto_shell_B%number_of_functions

        if (RC .le. tol) then
           !We need projections up to Lmax = Lmax_bspline+cgto_shell_C%l for the contraction with the Gaunt coefficients (see semi_analytic_BBGG_shell_integrals in bto_gto_integrals).
           n_projections = (grid_r1_r2%max_bspline_l+cgto_shell_C%l+1)**2
        else
           n_projections = (grid_r1_r2%max_bspline_l+1)**2
        endif

        BG_size = (grid_r1_r2%last_bspline_inner*n_projections)*(cgto_shell_C%number_of_functions)

        if (allocated(this%BXGG)) deallocate(this%BXGG)
        allocate(this%BXGG(BG_size,GG_size),stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_BGGG_integrals','Memory allocation failed.',err,1)

        same_ab_shells = this%shell_A_starting_index .eq. this%shell_B_starting_index

        !Evaluate integrals over all pairs of shells of BTOs
        call this%sph_BXGG_shell_integrals(this%cgto_shell_B%number_of_primitives,xb,yb,zb,this%cgto_shell_B%norm,&
                                           this%cgto_shell_B%norms,this%cgto_shell_B%l,this%cgto_shell_B%exponents,&
                                           this%cgto_shell_B%contractions,&
                                           this%cgto_shell_A%number_of_primitives,xa,ya,za,this%cgto_shell_A%norm,&
                                           this%cgto_shell_A%norms,this%cgto_shell_A%l,this%cgto_shell_A%exponents,&
                                           this%cgto_shell_A%contractions,&
                                           same_ab_shells, cgto_shell_C)

  end subroutine eval_BGGG_integrals

  subroutine eval_BBGG_integrals(this)
     implicit none
     class(GG_shell_integrals_obj) :: this

     integer :: lena,lenb,la,lb,err,BB_size,GG_size,n_projections
     real(kind=cfp) :: xa,ya,za, xb,yb,zb, acnorm,bcnorm
     logical :: same_ab_shells

        if (.not. this%initialized) call xermsg('GG_shell_mixed_integrals_mod','eval_BBGG_integrals','Object not initialized.',1,1)

        xa = this%cgto_shell_A%center(1)
        ya = this%cgto_shell_A%center(2)
        za = this%cgto_shell_A%center(3)

        xb = this%cgto_shell_B%center(1)
        yb = this%cgto_shell_B%center(2)
        zb = this%cgto_shell_B%center(3)

        GG_size = this%cgto_shell_A%number_of_functions*this%cgto_shell_B%number_of_functions

        n_projections = (2*grid_r1_r2%max_bspline_l+1)**2 !We need projections up to Lmax = 2*Lmax_bspline
        BB_size = grid_r1_r2%n_unique_pairs*n_projections

        if (allocated(this%BXGG)) deallocate(this%BXGG)
        allocate(this%BXGG(BB_size,GG_size),stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','eval_BBGG_integrals','Memory allocation failed.',err,1)

        same_ab_shells = this%shell_A_starting_index .eq. this%shell_B_starting_index

        !Evaluate integrals over all pairs of shells of BTOs
        call this%sph_BXGG_shell_integrals(this%cgto_shell_B%number_of_primitives,xb,yb,zb,this%cgto_shell_B%norm,&
                                           this%cgto_shell_B%norms,this%cgto_shell_B%l,this%cgto_shell_B%exponents,&
                                           this%cgto_shell_B%contractions,&
                                           this%cgto_shell_A%number_of_primitives,xa,ya,za,this%cgto_shell_A%norm,&
                                           this%cgto_shell_A%norms,this%cgto_shell_A%l,this%cgto_shell_A%exponents,&
                                           this%cgto_shell_A%contractions,&
                                           same_ab_shells)

  end subroutine eval_BBGG_integrals

  !> Capital X stands for one of: BTO, CGTO.
  subroutine sph_BXGG_shell_integrals(this, lena,xa,ya,za,acnorm,anorms,la,aexps,acoefs, &
                                      lenb,xb,yb,zb,bcnorm,bnorms,lb,bexps,bcoefs, same_ab_shells, cgto_shell_C)
     use const_gbl, only: mmax, imax_wp,boys_f_grid_step_wp,taylor_k_wp, imax_ep,boys_f_grid_step_ep,taylor_k_ep
     use phys_const_gbl, only: twopi, fourpi
     use eri_sph_coord_gbl, only: allocate_cf_space, Lag_n_hlf_k, cnla, cfp_solh_1d
     implicit none
     class(GG_shell_integrals_obj) :: this
     integer, intent(in) :: lena, lenb, la, lb
     real(kind=cfp), intent(in) :: xa,ya,za, xb,yb,zb
     real(kind=cfp), intent(in) :: acnorm,anorms(lena),aexps(lena),acoefs(lena)
     real(kind=cfp), intent(in) :: bcnorm,bnorms(lenb),bexps(lenb),bcoefs(lenb)
     logical, intent(in) :: same_ab_shells
     type(CGTO_shell_data_obj), optional :: cgto_shell_C

     real(kind=cfp), parameter :: inv_fourpi_sq = 1.0_cfp/(fourpi*fourpi)
     real(kind=cfp) :: Ra(3), Rb(3), R_prod(3), d, R_ba(3), d_ba, tol, alp, factor, norm, zeta, xsi, &
                       exp_T, exp_xsi, cnla_m1, T, Lag, rg, G, val
     integer :: n_l, i, j, sum_l, n_k, k, l, m, lap, map, ma, mb, lapp, mapp, lbpp, mbpp, Lboys, Mboys, ind, &
                n, lbp, mbp, base, err, n_lm, l_min,l_max, p,q, n_dummy, jj
     integer :: max_boys_l,max_laguerre_l,max_projection_l,d1,d2,d3,n_min,n_max,l_min_lag,l_max_lag,d4
     real(kind=cfp), allocatable :: solid_harmonics_ba(:), Boys_SLM(:,:), phi_nlm_xsi(:,:,:), pow_expa(:,:), pow_expb(:,:)
     real(kind=cfp), allocatable :: G_a_nz(:,:,:), G_b_nz(:,:,:), I_val(:,:), J_val(:,:), contracted(:)
     integer, allocatable :: G_a_ind(:,:,:), G_b_ind(:,:,:)
     real(kind=cfp) :: adj_norms_a(lena), adj_norms_b(lenb), prod_pow(lena*lenb), coupling, test, RC
     real(kind=cfp), allocatable :: Boys_Xlm(:,:,:), Boys_BX_Xlm(:,:,:), Laguerre_Xlm(:,:,:,:), Laguerre_BX_Xlm(:,:,:,:), &
                                    cgto_on_grid(:), Boys_G_Xlm(:,:,:,:), Laguerre_G_Xlm(:,:,:,:,:)
     type(pw_expansion_obj) :: u_grid
     logical :: I_is_zero, all_I_are_zero

        sum_l = la+lb
        n_l = sum_l+1
        n_max = (sum_l+2)/2

        !Account for possible equality of the A,B shells.
        if (same_ab_shells) then
           n_k = lena*(lena+1)/2
        else
           n_k = lena*lenb
        endif

        !Exponential grid in u for interval [0;1] with points accumulating towards 0.
        !val = 20.0_cfp !exponents larger than approx. 20 lead to a too small spacing between the first couple of integration intervals (for 50 intervals total).
        val = 20.0_cfp
        call u_grid%eval_exponential_grid(val,55,.false.)

        tol = F1MACH(4,cfp_dummy)

        if (present(cgto_shell_C)) then
           RC = sqrt(dot_product(cgto_shell_C%center,cgto_shell_C%center))
        endif

        max_boys_l = la + lb
        max_laguerre_l = max_boys_l - 2

        if (present(cgto_shell_C)) then
           !Is the CGTO sitting on the CMS?
           if (RC .le. tol) then
              !Evaluate the radial part of the CGTO on the radial grid.
              call CMS_CGTO_on_grid(grid_r1_r2%r1,grid_r1_r2%n1_total_points,cgto_shell_C,cgto_on_grid)

              !We need the X_lm projections for max_bspline_l+cgto_shell%l since we'll need to contract the final integrals with the coupling coefficients
              !to produce integrals for a pair of real spherical harms.
              max_projection_l = grid_r1_r2%max_bspline_l + cgto_shell_C%l
           else
              !Projections needed only up to the largest angular momentum of the B-spline.
              !The auxiliary projections needed to translate the CGTO_C are handled internally in the Boys and Laguerre routines.
              max_projection_l = grid_r1_r2%max_bspline_l
           endif
        else
           !We need the X_lm projections for 2*max_bspline_l since we'll need to contract the final integrals with the coupling coefficients to produce integrals for a pair of BTOs       
           max_projection_l = 2*grid_r1_r2%max_bspline_l
        endif

        !space for the Laguerre coefficients
        call allocate_cf_space(la,lb)

        n_lm = n_l**2
        i = max(la,1)
        j = max(lb,1)
        k = (2*i+1)*(2*j+1)
        l = max((la+lb)/2,1)
        m = max(n_max,1)
        allocate(solid_harmonics_ba(n_lm),Boys_SLM(n_lm,n_k),phi_nlm_xsi(n_lm,0:l,n_k),pow_expa(n_k,0:i),pow_expb(n_k,0:j),&
        &G_a_nz(2,-i:i,-i:i),G_b_nz(2,-j:j,-j:j),G_a_ind(2,-i:i,-i:i),G_b_ind(2,-j:j,-j:j),stat=err)
        if (err /= 0) stop "sph_BXGG_shell_integrals: memory allocation failed"

        factor = sqrt(fourpi/(2*la+1.0_cfp))
        do i=1,lena
           adj_norms_a(i) = acnorm*anorms(i)/((2*aexps(i))**(la+1.5_cfp))*factor*acoefs(i)
        enddo

        factor = sqrt(fourpi/(2*lb+1.0_cfp))
        do i=1,lenb
           adj_norms_b(i) = bcnorm*bnorms(i)/((2*bexps(i))**(lb+1.5_cfp))*factor*bcoefs(i)
        enddo

        factor = 8*(-1)**lb*(twopi)**3

        Ra(1:3) = (/xa,ya,za/)
        Rb(1:3) = (/xb,yb,zb/)
        R_ba = Rb - Ra
        d_ba = dot_product(R_ba,R_ba)

        !We calculate the solid harmonics including the leading factor sqrt(fourpi/(2*l+1)) but we remove it below since the formulae are for r**l*Xlm
        if (d_ba .le. tol) then
           solid_harmonics_ba = 0.0_cfp
           solid_harmonics_ba(1) = 1.0_cfp
        else
           call cfp_solh_1d(solid_harmonics_ba,R_ba(1),R_ba(2),R_ba(3),sum_l)
        endif

        do l=0,sum_l
           base = l*l
           m = 2*l+1
           solid_harmonics_ba(base+1:base+m) = sqrt((2*l+1.0_cfp)/fourpi)*solid_harmonics_ba(base+1:base+m)
        enddo

        n_min = -1
        n_max = -1
        l_min_lag = -1
        l_max_lag = -1
        do lap = 0,la
           do lbp = 0,lb
              l_min = abs(lap-lbp)
              l_max = lap+lbp-2
              if (l_max .ge. l_min) then
                 i = (lap+lbp-l_max)/2
                 if (n_min /= -1) then
                    n_min = min(n_min,i)
                 else
                    n_min = i
                 endif
                 n_max = max(n_max,(lap+lbp-l_min)/2)
                 if (l_min_lag /= -1) then
                    l_min_lag = min(l_min_lag,l_min)
                 else
                    l_min_lag = l_min
                 endif
                 l_max_lag = max(l_max_lag,l_max)
              endif
           enddo
        enddo

        !precalculate all G factors and Gaunt coefficients
        call cpl%prec_G_cf(sum_l+max_projection_l+max(n_max-1,0))

        !Calculate the phi_nlm function multipled by the corresponding solid harmonic and normalization coefficients
        !todo include in the input a flag telling me if the two shells are the same in which case I need to loop only over the unique combinations of contractions and multiply the off-diagonal by 2.
        k = 0
        do i=1,lena

           jj = lenb
           if (same_ab_shells) jj = i

           do j=1,jj
              k = k + 1

              norm = adj_norms_a(i)*adj_norms_b(j)*factor
              if (same_ab_shells .and. i /= j) norm = 2.0_cfp*norm !Account for equivalence of the (i,j) and (j,i) contractios.

              zeta = aexps(i) + bexps(j)
              xsi = aexps(i)*bexps(j)/zeta
              R_prod = (aexps(i)*Ra+bexps(j)*Rb)/zeta

              if (present(cgto_shell_C)) then !BGGG class

                 if (RC .le. tol) then
                    !For CGTO sitting on the CMS the BG projections are the same as for the BB case. Only the radial integration must be done differently.

                    !Obtain the values of the Boys function projected on the real spherical harmonics X_lm for all r1 grid points.
                    call this%Boys_projected_on_X_lm(u_grid,norm,zeta,R_prod,grid_r1_r2%r1,max_projection_l,max_boys_l,&
                                                     Boys_Xlm)

                    !Integrate the projected Boys function over all B-splines and the radial part of the CGTO on the r1 grid.
                    call r1_integrate_Boys_B_cms_G(grid_r1_r2,Boys_Xlm,k,n_k,cgto_on_grid,Boys_BX_Xlm)

                    if (n_min /= -1) then
                       call this%Laguerre_GTO_projected_on_X_lm(norm,zeta,R_prod,grid_r1_r2%r1,&
                                                                max_projection_l,l_min_lag,l_max_lag,n_min,n_max,Laguerre_Xlm)
                       call r1_integrate_Laguerre_GTO_B_cms_G(grid_r1_r2,Laguerre_Xlm,k,n_k,n_min,n_max,&
                                                              cgto_on_grid,Laguerre_BX_Xlm)
                    endif
                 else
                    !We need to project angularly the CGTO on the Boys function and on the Laguerre GTO.

                    !Obtain the values of the Boys function projected on the real spherical harmonics X_lm for all r1 grid points.
                    call this%Boys_projected_on_G_X_lm(u_grid,norm,zeta,R_prod,cgto_shell_C,grid_r1_r2%r1,&
                                                       max_projection_l,max_boys_l,Boys_G_Xlm)
   
                    !Integrate the projected Boys function over all B-splines and the angular projections of the CGTO on the r1 grid.
                    call r1_integrate_Boys_BG(grid_r1_r2,Boys_G_Xlm,k,n_k,Boys_BX_Xlm)
                    
                    if (n_min /= -1) then
                       call this%Laguerre_GTO_projected_on_G_X_lm(norm,zeta,R_prod,cgto_shell_C,grid_r1_r2%r1,&
                                                                  max_projection_l,l_min_lag,l_max_lag,n_min,n_max,Laguerre_G_Xlm)
                       call r1_integrate_Laguerre_GTO_BG(grid_r1_r2,Laguerre_G_Xlm,k,n_k,n_min,n_max,Laguerre_BX_Xlm)
                    endif
                 endif !RC .le. tol

              else !BBGG class

                 !Obtain the values of the Boys function projected on the real spherical harmonics X_lm for all r1 grid points.
                 call this%Boys_projected_on_X_lm(u_grid,norm,zeta,R_prod,grid_r1_r2%r1,max_projection_l,max_boys_l,Boys_Xlm)

                 !Integrate the projected Boys function over all pairs of B-splines on the r1 grid.
                 call r1_integrate_Boys_BB(grid_r1_r2,Boys_Xlm,k,n_k,Boys_BX_Xlm)
                 
                 if (n_min /= -1) then
                    call this%Laguerre_GTO_projected_on_X_lm(norm,zeta,R_prod,grid_r1_r2%r1,&
                                                             max_projection_l,l_min_lag,l_max_lag,n_min,n_max,Laguerre_Xlm)
                    call r1_integrate_Laguerre_GTO_BB(grid_r1_r2,Laguerre_Xlm,k,n_k,n_min,n_max,Laguerre_BX_Xlm)
                 endif
              endif

              do lap=0,la
                 pow_expa(k,lap) = (-aexps(i)/zeta)**lap
              enddo

              do lbp=0,lb
                 pow_expb(k,lbp) = (bexps(j)/zeta)**lbp
              enddo

              alp = 1.0_cfp/(4.0_cfp*xsi)
              xsi = xsi*d_ba
              exp_xsi = exp(-xsi) !*norm
              do l=0,sum_l
                 do n=0,(sum_l/2)
                    Lag = cnla(n,l,alp)*exp_xsi*Lag_n_hlf_k(n,l,xsi)*(-1)**n
                    do m=-l,l
                       ind = l*l+l+m+1
                       phi_nlm_xsi(ind,n,k) = Lag*solid_harmonics_ba(ind)
                    enddo !m
                 enddo !n
              enddo

           enddo !j
        enddo !i

        this%BXGG = 0.0_cfp

        !Below the loops over i or j up to n_dummy correspond to loops over the indices of the BX pair of functions.
        !The algorithm is otherwise identical to the spherical algorithm for NAI for a pair GG shells.
        n_dummy = size(Boys_BX_Xlm,1)
        i = max(la,1)
        j = max(lb,1)
        k = (2*i+1)*(2*j+1)
        allocate(I_val(n_k,k),J_val(n_dummy,n_k),contracted(n_dummy),stat=err)
        if (err /= 0) stop "sph_BXGG_shell_integrals: memory allocation 2 failed"

        do lbp=0,lb
           do lap=0,la

              lapp = la-lap
              lbpp = lb-lbp

              do k=1,n_k
                 prod_pow(k) = pow_expa(k,lap)*pow_expb(k,lbp)*inv_fourpi_sq
              enddo

              !Preload the coupling coefficients
              G_a_nz = 0.0_cfp
              G_a_ind = 0
              do mapp=-lapp,lapp
                 do map=-lap,lap
                    i = 0
                    do ma = -la,la
                       G = cpl%G_real_cf(la,lap,ma,map,mapp)
                       if (G /= 0.0_cfp) then
                          i=i+1
                          G_a_nz(i,map,mapp) = G
                          G_a_ind(i,map,mapp) = ma
                       endif
                    enddo !ma
                 enddo !map
              enddo !mapp

              G_b_nz = 0.0_cfp
              G_b_ind = 0
              do mbpp=-lbpp,lbpp
                 do mbp=-lbp,lbp
                    i = 0
                    do mb = -lb,lb
                       G = cpl%G_real_cf(lb,lbp,mb,mbp,mbpp)
                       if (G /= 0.0_cfp) then
                          i=i+1
                          G_b_nz(i,mbp,mbpp) = G
                          G_b_ind(i,mbp,mbpp) = mb
                       endif
                    enddo !mb
                 enddo !mbp
              enddo !mbpp

              !Construct the I-terms
              all_I_are_zero = .true.
              i = 0
              do mapp=-lapp,lapp
                 do mbpp=-lbpp,lbpp
                    i = i + 1

                    I_is_zero = .true.

                    call cpl%bounds_rg(lapp,lbpp,mapp,mbpp,l_min,l_max)
                    I_val(1:n_k,i) = 0.0_cfp
                    !todo loop only over those l,m,n which have non-zero phi_nlm_xsi 
                    do l=l_max,l_min,-2
                       n = (lapp+lbpp-l)/2
                       do m=-l,l
                          ind = l*l+l+m+1
                          rg = cpl%rgaunt(lapp,lbpp,l,mapp,mbpp,m)
                          test = rg*solid_harmonics_ba(ind)
                          if (test /= 0.0_cfp) then
                             do k=1,n_k
                                I_val(k,i) = I_val(k,i) + rg*phi_nlm_xsi(ind,n,k)
                             enddo
                          endif
                       enddo !m
                    enddo !l

                    do k=1,n_k
                       if (I_val(k,i) /= 0.0_cfp) then
                          I_is_zero = .false.
                          exit
                       endif
                    enddo

                    if (.not.(I_is_zero)) then
                       I_val(1:n_k,i) = I_val(1:n_k,i)*prod_pow(1:n_k)
                       all_I_are_zero = .false.
                    endif

                 enddo !mbpp
              enddo !mapp

              if (all_I_are_zero) cycle

              !Contract over the primed indices
              do mbp=-lbp,lbp
                 do map=-lap,lap

                    Lboys = lap+lbp
           
                    J_val(1:n_dummy,1:n_k) = 0.0_cfp
                    do Mboys = -Lboys,Lboys
                       rg = cpl%rgaunt(lbp,lap,Lboys,mbp,map,Mboys)
                       if (rg /= 0.0_cfp) then
                          !todo save in Boys_BX_Xlm only those L,M which are non-zero similar to angular_integrals.
                          ind = Lboys*Lboys+Lboys+Mboys+1
                          do k=1,n_k
                             do i=1,n_dummy
                                J_val(i,k) = J_val(i,k) + rg*Boys_BX_Xlm(i,k,ind)
                             enddo !i
                          enddo
                       endif
                    enddo

                    call cpl%bounds_rg(lbp,lap,mbp,map,l_min,l_max)
                    do l=min(Lboys-2,l_max),l_min,-2 !Lboys-2,abs(lap-lbp),-1
                       n = (lap+lbp-l)/2                  
                       if (n < n_min .or. n > n_max) stop "error: n out of range"
                       if (l < l_min_lag .or. l > l_max_lag) stop "error: l out of range"
                       do m=-l,l
                          rg = cpl%rgaunt(lbp,lap,l,mbp,map,m)
                          if (rg /= 0.0_cfp) then
                             ind = l*l+l+m+1
                             !todo save in Laguerre_BX_Xlm only those l,m which are non-zero similar to angular_integrals.
                             do k=1,n_k
                                do i=1,n_dummy
                                   J_val(i,k) = J_val(i,k) + rg*Laguerre_BX_Xlm(i,k,ind,n) !n must correspond to n-1 being used to evaluate Laguerre_X_X_nlm!!
                                enddo !i
                             enddo
                          endif
                       enddo !m
                    enddo !l

                    !Perform the contraction over primitives and contraction over the pp indices
                    i = 0
                    do mapp=-lapp,lapp
                       do mbpp=-lbpp,lbpp
                          i = i + 1
                          do j=1,n_dummy
                             contracted(j) = sum(I_val(1:n_k,i)*J_val(j,1:n_k)) !todo this is a matrix*matrix multiplication: I should do gemm(J_val,I_val,contr); contr(1:n_dummy,1:n_pp) !!!
                          enddo !j
                          !if (contracted .eq. 0.0_cfp) cycle

                          !At most two values of ma and mb are allowed
                          !todo this should be rewritten so that for each ma,mb combination I'll have a matrix Gab(1:n_pp,1:n_ma_mb) and then do gemm(contr,Gab,BXGG) adding to BXGG
                          do p=1,2
                             mb = G_b_ind(p,mbp,mbpp)
                             ind = (2*la+1)*(mb+lb) + la+1
                             do q=1,2
                                ma = G_a_ind(q,map,mapp)
                                coupling = G_a_nz(q,map,mapp)*G_b_nz(p,mbp,mbpp)
                                if (coupling .eq. 0.0_cfp) cycle
                                do j=1,n_dummy
                                   this%BXGG(j,ind+ma) = this%BXGG(j,ind+ma) + contracted(j)*coupling
                                enddo !j
                             enddo !q
                          enddo !p

                       enddo !mbpp
                    enddo !mapp

                 enddo !map
              enddo !mbp

           enddo !lap
        enddo !lbp

  end subroutine sph_BXGG_shell_integrals

  !> Calculates the angular integrals \int d\Omega(r_{1}) X_{l,m}(r_{1}) CGTO_{Lg,Mg}(r_{1}-A)*C_{n-1,l}^{a}(1/(4*zeta))*\phi_{n-1,lp,mp}^{b}(zeta,(zeta_center-r_{1})),
  !> for the grid of r1 points, l,m projections up to max_projection_l, lp,mp between l_min_lag and l_max_lag and n values between n_min,n_max. 
  !> \phi_{n-1,lp,mp}^{b}(zeta,(zeta_center-r_{1})) = exp(-zeta*(r_{1}-zeta_center)**2)*L_{n-1}^{lp+1/2}(zeta*(r_{1}-zeta_center)**2)*S_{lp,mp}(zeta_center-r_{1})*sqrt((2*lp+1)/(4*pi)).
  subroutine Laguerre_GTO_projected_on_G_X_lm(this,norm,zeta,zeta_center,cgto_shell,r,&
                                              max_projection_l,l_min_lag,l_max_lag,n_min,n_max,Laguerre_G_X_nlm)
     use eri_sph_coord_gbl, only: cnla
     use phys_const_gbl, only: fourpi, pi
     use special_functions_gbl, only: cfp_besi, cfp_eval_poly_horner_many
     use lag_cfs_gbl, only: Laguerre_poly_factorization_cfs, l_lag_lim, n_lag_lim
     implicit none
     class(GG_shell_integrals_obj) :: this
     real(kind=cfp), intent(in) :: zeta, zeta_center(3), norm
     real(kind=cfp), allocatable :: r(:)
     integer, intent(in) :: max_projection_l, n_min, n_max, l_min_lag, l_max_lag
     type(CGTO_shell_data_obj) :: cgto_shell
     !OUT:
     real(kind=cfp), allocatable :: Laguerre_G_X_nlm(:,:,:,:,:) !radial_point,Mg,lm,lpmp,n

     integer, parameter :: kode = 2
     real(kind=cfp), parameter :: half = 0.5_cfp

     integer :: n, i, j, err, n_lm, n_lpmp, max_l_aux, l, m, lp, besi_dim, n_lag, lm, nz, l_lag, m_lag, m_ind, &
                lpmp, l_min, l_max, p, n_mg, n_lm_aux
     integer :: ll, ml, lpp, mpp, lag_lm, l_aux, m_aux, lambda, mu, lm_aux, n_m1, k, mp, ind, lppmpp, Mg, sum_l, n_lm_sum, &
                base, LbMb, lpmup, lambdap,mup, lmu, la_p_lb
     real(kind=cfp) :: R_dist, R_zeta, arg, exp_arg, tol, val, exp_fac, bessel_fac, exp_part, alp, cf, prod_center(3), prod_exp, &
                       preexp_factor, Rab_sq, Rab(3), Rp, RC
     real(kind=cfp), allocatable :: Laguerre_aux(:,:,:), Xlm_CGTO_center(:), cnla_fac(:,:), y(:), d_cfs(:,:), c_lambda(:,:,:), &
                                    transl_cfs(:,:,:), c(:,:,:,:,:), d(:,:,:), P_lag(:,:), tmp(:)
     real(kind=cfp), allocatable :: r1_x_r2(:), x_xx_p(:), outer_poly(:,:), Xlm_product_center(:,:), R_aux(:,:), &
                                    Xlm_CGTO_C_center(:), transl_cfs_all_AB(:,:,:,:,:), transl_cfs_AB(:,:,:,:,:)
     real(kind=cfp), allocatable :: Laguerre_aux_contr(:,:), dpp(:,:,:), transl_cf(:,:)
     logical :: is_zero

        n = size(r)
        n_lm = (max_projection_l+1)**2
        n_lpmp = (l_max_lag+1)**2
        n_mg = 2*cgto_shell%l+1

        max_l_aux = max_projection_l + l_max_lag + max(n_max-1,0) + cgto_shell%l
        alp = 1.0_cfp/(4.0_cfp*zeta)
        R_dist = sqrt(dot_product(zeta_center,zeta_center))
        R_zeta = zeta*R_dist

        l = max(1,max_projection_l)
        lp = max(1,l_max_lag)
        besi_dim = max_l_aux + 1
        n_lm_aux = besi_dim**2
        allocate(Laguerre_aux(0:besi_dim,n,cgto_shell%number_of_primitives),cnla_fac(n_min:n_max,0:lp),y(besi_dim),&
                 Xlm_product_center(n_lm_aux,cgto_shell%number_of_primitives),stat=err)
        if (err /= 0) then
           call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_G_X_lm','Memory allocation 1 failed.',err,1)
        endif

        tol = F1MACH(4,cfp_dummy)
 
        if (allocated(Laguerre_G_X_nlm)) then
           if (size(Laguerre_G_X_nlm,1) /= n .or. size(Laguerre_G_X_nlm,2) /= n_mg .or. &
               size(Laguerre_G_X_nlm,3) /= n_lm .or. size(Laguerre_G_X_nlm,4) /= n_lpmp .or. &
               size(Laguerre_G_X_nlm,5) /= n_max) deallocate(Laguerre_G_X_nlm)
        endif
 
        if (.not. allocated(Laguerre_G_X_nlm)) then
           allocate(Laguerre_G_X_nlm(n,n_mg,n_lm,n_lpmp,n_min:n_max),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_G_X_lm','Memory allocation 2 failed.',err,1)
           endif
        endif

        Laguerre_G_X_nlm = 0.0_cfp

        cnla_fac = 0.0_cfp
        do l_lag=l_min_lag,l_max_lag
           val = sqrt((2*l_lag+1.0_cfp)/fourpi) !get rid of the leading factor sqrt(4pi/(2*l+1)) for the solid harmonics since the NAI formulae are only for r**l*Xlm
           do n_lag=n_min,n_max
              cnla_fac(n_lag,l_lag) = norm*val*fourpi*(-1)**n_lag*cnla(n_lag-1,l_lag,alp)*(-1)**l_lag
           enddo
        enddo

        Rab = zeta_center-cgto_shell%center
        Rab_sq = dot_product(Rab,Rab)

        Laguerre_aux = 0.0_cfp
        do j=1,cgto_shell%number_of_primitives

           prod_exp = zeta + cgto_shell%exponents(j)
           prod_center = (zeta*zeta_center + cgto_shell%exponents(j)*cgto_shell%center)/prod_exp
           Rp = sqrt(dot_product(prod_center,prod_center))
           preexp_factor = cgto_shell%norm*cgto_shell%norms(j)*cgto_shell%contractions(j)&
                           *exp(-zeta*cgto_shell%exponents(j)/prod_exp*Rab_sq)

           !Real spherical harmonics for the product GTO center: result in the array Xlm_CGTO_center
           call real_harmonics%precalculate_Xlm_for_CGTO_center(prod_center,max_l_aux,Xlm_CGTO_center)

           Xlm_product_center(1:n_lm_aux,j) = Xlm_CGTO_center(1:n_lm_aux)

           !Evaluate the partial wave projections of the GTO for the partial waves up to max_l_aux as needed later to perform the translations of the polynomials.
           do i=1,n
              arg = 2.0_cfp*r(i)*prod_exp*Rp
              exp_arg = prod_exp*(r(i)-Rp)**2                 
              exp_fac = preexp_factor*exp(-exp_arg)
   
              if (arg .le. tol) then !Evaluate limit for arg -> 0
                 l = 0
                 bessel_fac = 1.0_cfp !only l=0 Bessel function is non-zero
                 exp_part = bessel_fac*exp_fac
                 Laguerre_aux(l,i,j) = exp_part
              else
                 call cfp_besi(arg, half, kode, besi_dim, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
                 val = sqrt(pi/(2.0_cfp*arg))
                 do l=0,max_l_aux
                    bessel_fac = y(l+1)*val
                    exp_part = bessel_fac*exp_fac
                    Laguerre_aux(l,i,j) = exp_part
                 enddo !l
              endif
           enddo !i
        enddo !j

        if (n_max .eq. 1) then !L_{0}^{lp+1/2}(x) = 1: the polynomial is trivial

           if (l_max_lag .eq. 0) then 
              !s-type Laguerre GTO: projection of exp(-zeta*(r_{1}-zeta_center)**2) so no translation of S_lpmp needed.

              if (cgto_shell%l .eq. 0) then

                 do j=1,cgto_shell%number_of_primitives
                    do l=0,max_projection_l
                       do m=-l,l
                          lm = l*l+l+m+1
                          Laguerre_G_X_nlm(1:n,1,lm,1,1) = Laguerre_G_X_nlm(1:n,1,lm,1,1)&
                                                         + Laguerre_aux(l,1:n,j)*Xlm_product_center(lm,j)*cnla_fac(1,0)
                       enddo !m
                    enddo !l
                 enddo !j

              else
                 !The CGTO solid harmonic must be translated

                 !todo allocate only once and use the save attribute
                 allocate(d_cfs(n,1:cgto_shell%l+1),R_aux(n,(cgto_shell%l+1)**2),tmp(n),stat=err)
                 if (err /= 0) then
                    call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_G_X_lm',&
                                'Memory allocation 3 failed.',err,1)
                 endif
    
                 R_dist = sqrt(dot_product(cgto_shell%center,cgto_shell%center))
    
                 !Real spherical harmonics for the CGTO center: result in the array Xlm_CGTO_center
                 call real_harmonics%precalculate_Xlm_for_CGTO_center(cgto_shell%center,cgto_shell%l,Xlm_CGTO_center)

                 !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
                 call this%precalculate_solh_translation_coeffs(cgto_shell%l,R_dist,Xlm_CGTO_center,transl_cfs)

                 do l=0,max_projection_l
                    do m=-l,l
                       lm = l*l+l+m+1
       
                       R_aux = 0.0_cfp
                       do lp=0,cgto_shell%l
                          do mp=-lp,lp
                             lpmp = lp*lp+lp+mp+1
       
                             !Couple to the spherical harmonic X_{l,m}
                             do lpp=abs(lp-l),lp+l
                                lppmpp = lpp*lpp+lpp+1
                                do mpp=-lpp,lpp
                                   cf = cpl%rgaunt(l,lpp,lp,m,mpp,mp)*cnla_fac(1,0)
                                   if (cf /= 0.0_cfp) then
                                      do j=1,cgto_shell%number_of_primitives
                                         if (Xlm_product_center(lppmpp+mpp,j) /= 0.0_cfp) then
                                            R_aux(1:n,lpmp) = R_aux(1:n,lpmp)&
                                                            + cf*Laguerre_aux(lpp,1:n,j)*Xlm_product_center(lppmpp+mpp,j)
                                         endif
                                      enddo !j
                                   endif
                                enddo !mpp
                             enddo !lpp
       
                          enddo !mp
                       enddo !lp
       
                       do Mg=1,cgto_shell%number_of_functions
       
                          !Translate the solid harmonic S_{cgto_shell%l,M} of the CGTO
                          d_cfs = 0.0_cfp
                          do lp=0,cgto_shell%l
                             lpmp = lp*lp+lp+1
       
                             do mp=-lp,lp
                                if (transl_cfs(mp+lp+1,lp,Mg) /= 0.0_cfp) then
                                   d_cfs(1:n,lp+1) = d_cfs(1:n,lp+1) + transl_cfs(mp+lp+1,lp,Mg)*R_aux(1:n,lpmp+mp)
                                endif
                             enddo !mp
                          enddo !lp
       
                          !Evaluate the polynomial for each r(i)
                          call cfp_eval_poly_horner_many(cgto_shell%l,r,n,d_cfs,tmp(1:n))

                          !Contraction over primitives
                          Laguerre_G_X_nlm(1:n,Mg,lm,1,1) = Laguerre_G_X_nlm(1:n,Mg,lm,1,1) + tmp(1:n)
       
                       enddo !Mg
        
                    enddo !m
                 enddo !l

              endif

           else  !l_max_lag > 0
              !Laguerre GTO with L > 0: we need to translate the S_lpmp solid harmonic part. We also include translation of the CGTO solid harmonic part.
              sum_l = l_max_lag+cgto_shell%l
              n_lm_sum = (sum_l+1)**2
              i = 2*sum_l+1
              ind = (2*l_max_lag+1)*(2*cgto_shell%l+1)
              !todo allocate only once and use the save attribute
              allocate(d_cfs(n,1:sum_l+1),R_aux(n,cgto_shell%number_of_primitives*n_lm_sum*(l_max_lag+1)),&
                       transl_cfs_all_AB(i,0:sum_l,0:sum_l,ind,l_max_lag+1),tmp(n),stat=err)
              if (err /= 0) then
                 call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_G_X_lm','Memory allocation 3 failed.',err,1)
              endif
 
              R_dist = sqrt(dot_product(zeta_center,zeta_center))
              RC = sqrt(dot_product(cgto_shell%center,cgto_shell%center)) 
 
              !Real spherical harmonics for the CGTO center: result in the module array Xlm_CGTO_C_center
              call real_harmonics%precalculate_Xlm_for_CGTO_center(cgto_shell%center,cgto_shell%l,Xlm_CGTO_C_center)
 
              !Real spherical harmonics for the Laguerre GTO center: result in the module array Xlm_CGTO_center
              call real_harmonics%precalculate_Xlm_for_CGTO_center(zeta_center,l_max_lag,Xlm_CGTO_center)
 
              !For every l_lag precalculate the coefficients in the translation formula for the pair of solid harmonics
              transl_cfs_all_AB = 0.0_cfp
              do l_lag=l_min_lag,l_max_lag
                 call this%precalculate_pair_solh_translation_coeffs(l_lag,R_dist,Xlm_CGTO_center,cgto_shell%l,&
                                                                     RC,Xlm_CGTO_C_center,transl_cfs_AB)
                 ind = 0
                 do m_lag=1,2*l_lag+1
                    do Mg=1,cgto_shell%number_of_functions
                       ind = ind + 1
                       do la_p_lb=0,l_lag+cgto_shell%l
                          do lp=0,l_lag+cgto_shell%l
                             transl_cfs_all_AB(1:2*lp+1,lp,la_p_lb,ind,l_lag+1) = transl_cfs_AB(1:2*lp+1,lp,la_p_lb,Mg,m_lag)
                          enddo !lp
                       enddo !la_p_lb
                    enddo !Mg
                 enddo !m_lag
              enddo !l_lag
 
              do l=0,max_projection_l
                 do m=-l,l
                    lm = l*l+l+m+1
  
                    !R_aux is effectively a 4D array (i,lpmp,l_lag,j) with the last
                    !three dimensions contracted into one.
                    R_aux = 0.0_cfp
                    do lp=0,sum_l
                       do mp=-lp,lp
                          lpmp = lp*lp+lp+mp+1
    
                          !Couple to the spherical harmonic X_{l,m}
                          do lpp=abs(lp-l),lp+l
                             lppmpp = lpp*lpp+lpp+1
                             do mpp=-lpp,lpp
                                cf = cpl%rgaunt(l,lpp,lp,m,mpp,mp)
                                if (cf /= 0.0_cfp) then
                                   do l_lag=l_min_lag,l_max_lag
                                      do j=1,cgto_shell%number_of_primitives
                                         if (Xlm_product_center(lppmpp+mpp,j) /= 0.0_cfp) then
                                            base = n_lm_sum*l_lag + n_lm_sum*(l_max_lag+1)*(j-1)
                                            R_aux(1:n,base+lpmp) = R_aux(1:n,base+lpmp)&
                                                                 +cf&
                                                                 *Laguerre_aux(lpp,1:n,j)&
                                                                 *cnla_fac(1,l_lag)&
                                                                 *Xlm_product_center(lppmpp+mpp,j)
                                         endif
                                      enddo !j
                                   enddo
                                endif
                             enddo !mpp
                          enddo !lpp
   
                       enddo !mp
                    enddo !lp
 
                    do l_lag=l_min_lag,l_max_lag
                       base = n_lm_sum*l_lag
                       ind = 0
                       do m_lag=1,2*l_lag+1
                          LbMb = l_lag*l_lag+m_lag
                          do Mg=1,cgto_shell%number_of_functions
                             ind = ind + 1
    
                             do j=1,cgto_shell%number_of_primitives
                                base = n_lm_sum*l_lag + n_lm_sum*(l_max_lag+1)*(j-1)

                                !Translate the solid harmonics S_{l_lag,M}*S_{Lg,Mg}
                                d_cfs = 0.0_cfp
                                do la_p_lb=0,sum_l
                                   do lp=0,sum_l 
                                      lpmp = base+lp*lp+lp+1
                                      do mp=-lp,lp
                                         if (transl_cfs_all_AB(mp+lp+1,lp,la_p_lb,ind,l_lag+1) /= 0.0_cfp) then
                                            !d_cfs(1:n,la_p_lb+1) = d_cfs(1:n,la_p_lb+1) + transl_cfs_all_AB(mp+lp+1,lp,la_p_lb,ind,l_lag+1)*R_aux(1:n,lpmp+mp)
                                            do i=1,n !in quad precision ifort 16.0.1 fails here but changing the line above into the explicit loop resolves the problem
                                               d_cfs(i,la_p_lb+1) = d_cfs(i,la_p_lb+1)&
                                                                  + transl_cfs_all_AB(mp+lp+1,lp,la_p_lb,ind,l_lag+1)&
                                                                  * R_aux(i,lpmp+mp)
                                            enddo
                                         endif
                                      enddo !mp
                                   enddo !lp
                                enddo !la_p_lb
       
                                !Evaluate the polynomial for each r(i)
                                call cfp_eval_poly_horner_many(sum_l,r,n,d_cfs,tmp(1:n))
                                
                                !Contraction over primitives
                                Laguerre_G_X_nlm(1:n,Mg,lm,LbMb,1) = Laguerre_G_X_nlm(1:n,Mg,lm,LbMb,1) + tmp(1:n)
                             enddo !j
 
                          enddo !Mg
                       enddo !m_lag
                    enddo !l_lag
 
                 enddo !m
              enddo !l
           endif

        else !n_max /= 1: L_{0}^{lp+1/2}(x) is non-trivial

           if (l_max_lag > l_lag_lim .or. n_max-1 > n_lag_lim) then
              call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_G_X_lm',&
                          'Need Laguerre coefficients which have not been generated. See module lag_cfs.',3,1)
           endif

           !Get the coefficients needed to factorize L(x+y) (x = scalar, y = dot_product(r1,r2)) in terms of the real spherical harmonics.
           call Laguerre_poly_factorization_cfs(c)

           if (l_max_lag .eq. 0) then 
              !s-type Laguerre GTO so we need to translate only the Laguerre polynomial. Does this ever happen???
              !The implementation of this case is easy: see the calculation of P_lag in the branch below.
              call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_G_X_lm',&
                          'Not implemented Laguerre translation for the case of l_max_lag .eq. 0.',2,1)

            else

              !Laguerre GTO with L > 0: we need to translate both the Laguerre polynomial and the S_lpmp part.
              n_m1 = n_max-1
              i = (n_m1 + l_max_lag +1)**2
              j = max(l_max_lag+1,n_max)
              lpmp = (n_m1+l_max_lag+cgto_shell%l+max_projection_l+1)**2
              allocate(P_lag(1:n,0:n_max-1),tmp(n),r1_X_r2(n),x_xx_p(n),transl_cf(n,(l_max_lag+cgto_shell%l+1)**2),&
                       Laguerre_aux_contr(n,lpmp),dpp(n,0:n_m1,(max_projection_l+l_max_lag+cgto_shell%l+1)**2),&
                       d_cfs(n,1:cgto_shell%l+l_max_lag+1),outer_poly(n,1:n_max),stat=err)
              if (err /= 0) then
                 call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_G_X_lm','Memory allocation 4 failed.',err,1)
              endif

              R_dist = sqrt(dot_product(zeta_center,zeta_center))
              RC = sqrt(dot_product(cgto_shell%center,cgto_shell%center)) 
 
              !Real spherical harmonics for the CGTO center: result in the module array Xlm_CGTO_C_center
              call real_harmonics%precalculate_Xlm_for_CGTO_center(cgto_shell%center,cgto_shell%l,Xlm_CGTO_C_center)
 
              !Real spherical harmonics for the Laguerre GTO center: result in the module array Xlm_CGTO_center
              call real_harmonics%precalculate_Xlm_for_CGTO_center(zeta_center,j,Xlm_CGTO_center)
 
              do i=1,n
                 x_xx_p(i) = zeta*(r(i)**2+R_dist**2)
                 r1_X_r2(i) = 2*zeta*r(i)*R_dist
              enddo !i

              !Contraction over primitives
              Laguerre_aux_contr = 0.0_cfp
              do lambdap=0,n_m1+l_max_lag+cgto_shell%l+max_projection_l
                 do mup=-lambdap,lambdap
                    lpmup = lambdap*lambdap+lambdap+mup+1
                    do j=1,cgto_shell%number_of_primitives
                       if (Xlm_product_center(lpmup,j) /= 0.0_cfp) then
                          Laguerre_aux_contr(1:n,lpmup) = Laguerre_aux_contr(1:n,lpmup)&
                                                        + Laguerre_aux(lambdap,1:n,j)&
                                                        * Xlm_product_center(lpmup,j)
                       endif
                    enddo !j
                 enddo !mup
              enddo !lambdap

              !Couple the Laguerre translation harmonics to the spherical harmonics of the solid harmonic translation and projection
              do lpp=0,max_projection_l+l_max_lag+cgto_shell%l

                 do mpp=-lpp,lpp
                    lppmpp = lpp*lpp+lpp+mpp+1
               
                    !Contraction over mu and lambdap,mup     
                    do lambda=0,n_m1
                       dpp(1:n,lambda,lppmpp) = 0.0_cfp
                       do mu=-lambda,lambda
                          lmu = lambda*lambda+lambda+mu+1
                          tmp(1:n) = Xlm_CGTO_center(lmu)
                          do lambdap=abs(lambda-lpp),lambda+lpp
                             do mup=-lambdap,lambdap
                                cf = cpl%rgaunt(lpp,lambda,lambdap,mpp,mu,mup)
                                if (cf /= 0.0_cfp) then
                                   lpmup = lambdap*lambdap+lambdap+mup+1
                                   dpp(1:n,lambda,lppmpp) = dpp(1:n,lambda,lppmpp) + cf*tmp(1:n)*Laguerre_aux_contr(1:n,lpmup)
                                endif
                             enddo !mup
                          enddo !lambdap
                       enddo !mu
                    enddo !lambda

                 enddo !mpp
              enddo !lpp

              do l_lag=l_min_lag,l_max_lag

                 call this%precalculate_pair_solh_translation_coeffs(l_lag,R_dist,Xlm_CGTO_center,cgto_shell%l,&
                                                                     RC,Xlm_CGTO_C_center,transl_cfs_AB)

                 do n_lag=n_min,n_max

                    !Ensure that phi_n-1,lm^{b} is taken:
                    n_m1 = n_lag-1

                    if (n_m1 .eq. 0) then 
                       !Laguerre polynomial for n-1 == 0 requested; this is 1; 
                       !the factor of 4pi compensates for the product X_00(r)*X_00(zeta_center)=1/(4*pi) coming from the expansion of L_{0}^{l+1/2}(x+y)
                       P_lag(1:n,0) = 1.0_cfp*cnla_fac(n_lag,l_lag)*fourpi
                    else
                       do lp=0,n_m1
                          !Generate P_{lp}^{n_lag-1,l_lag}*(-1)**lp*cnla_fac(n_lag,l_lag)
                          !Use the Horner scheme twice to evaluate the 2D polynomial:
                          !P_{lp}^{n_lag-1,l_lag} = sum_{p} (zeta*(r(i)**2+R_dist**2))**p * outer_poly_{p,r(i)}; outer_poly_{p,r(i)} = sum_{k} c(k,p,lp,l_lag,n_m1)* (2*zeta*r(i)*R_dist)
                          do p = 0,n_m1
                             do k=0,n_m1
                                d_cfs(1:n,k+1) = c(k,p,lp,l_lag,n_m1)
                             enddo !k
                             !Evaluate the inner polynomial in (2*zeta*r(i)*R_dist) for each value of r(i): the result are coefficients for the outer polynomial
                             call cfp_eval_poly_horner_many(n_m1,r1_X_r2,n,d_cfs,outer_poly(1:n,p+1))
                          enddo !p
                          !Evaluate the outer polynomial in (zeta*(r(i)**2+R_dist**2)) for each value of r(i):
                          call cfp_eval_poly_horner_many(n_m1,x_xx_p,n,outer_poly,P_lag(1:n,lp))
                          P_lag(1:n,lp) = P_lag(1:n,lp)*cnla_fac(n_lag,l_lag)*(-1)**lp
                       enddo !lp
                    endif

                    !Couple the Laguerre translation harmonics to the spherical harmonics of the solid harmonic translation and projection
                    Laguerre_aux_contr = 0.0_cfp !re-use this array
                    do lpp=0,max_projection_l+l_lag+cgto_shell%l

                       do mpp=-lpp,lpp
                          lppmpp = lpp*lpp+lpp+mpp+1
                     
                          !Contraction over lambda
                          do lambda=0,n_m1
                             Laguerre_aux_contr(1:n,lppmpp) = Laguerre_aux_contr(1:n,lppmpp)&
                                                            + dpp(1:n,lambda,lppmpp)&
                                                            * P_lag(1:n,lambda)
                          enddo !lambda

                       enddo !mpp
                    enddo !lpp
                          
                    do l=0,max_projection_l
                       do m=-l,l
                          lm = l*l+l+m+1

                          do lp=0,l_lag+cgto_shell%l
                             do mp=-lp,lp
                                lpmp = lp*lp+lp+mp+1

                                !Contraction over lpp,mpp
                                transl_cf(1:n,lpmp) = 0.0_cfp
                                do lpp=abs(lp-l),lp+l !,2
                                   do mpp=-lpp,lpp
                                      lppmpp = lpp*lpp+lpp+mpp+1
                                      cf = cpl%rgaunt(l,lp,lpp,m,mp,mpp)
                                      if (cf /= 0.0_cfp) then
                                         transl_cf(1:n,lpmp) = transl_cf(1:n,lpmp) + cf*Laguerre_aux_contr(1:n,lppmpp)
                                      endif
                                   enddo !mpp
                                enddo !lpp

                             enddo !mp
                          enddo !lp

                          do ml=-l_lag,l_lag
                             m_lag = ml+l_lag+1
                             lag_lm = l_lag*l_lag+m_lag
                             do Mg=1,cgto_shell%number_of_functions

                                !Translate the product of the two solid harmonics with the pair of M-values ml,Mg:
                                ll = l_lag+cgto_shell%l
                                do la_p_lb=0,ll

                                   d_cfs(1:n,la_p_lb+1) = 0.0_cfp
                                   do lp=0,l_lag+cgto_shell%l
                                      do mp=-lp,lp
                                         lpmp = lp*lp+lp+mp+1
                                         cf = transl_cfs_AB(lp+mp+1,lp,la_p_lb,Mg,m_lag)
                                         if (cf /= 0.0_cfp) then
                                            !d_cfs(1:n,la_p_lb+1) = d_cfs(1:n,la_p_lb+1) + cf*transl_cf(1:n,lpmp)
                                            do i=1,n !in quad precision ifort 16.0.1 fails here but changing the line above into the explicit loop resolves the problem
                                               d_cfs(i,la_p_lb+1) = d_cfs(i,la_p_lb+1) + cf*transl_cf(i,lpmp)
                                            enddo
                                         endif
                                      enddo !mp
                                   enddo !lp

                                enddo !la_p_lb

                                !Evaluate the polynomial for each r(i)
                                call cfp_eval_poly_horner_many(ll,r,n,d_cfs,Laguerre_G_X_nlm(1:n,Mg,lm,lag_lm,n_lag))

                             enddo !Mg
                          enddo !ml
                          
                       enddo !m
                    enddo !l

                 enddo !n_lag

              enddo !l_lag

           endif !l_max_lag .eq. 0

        endif !n_max .eq. 1

  end subroutine Laguerre_GTO_projected_on_G_X_lm

  !> Calculates the angular integrals \int d\Omega(r_{1}) X_{l,m}(r_{1}) C_{n-1,l}^{a}(1/(4*zeta))*\phi_{n-1,lp,mp}^{b}(zeta,(zeta_center-r_{1})) for the grid of r1 points, l,m projections up to max_projection_l, 
  !> lp,mp between l_min_lag and l_max_lag and n values between n_min,n_max. 
  !> \phi_{n-1,lp,mp}^{b}(zeta,(zeta_center-r_{1})) = exp(-zeta*(r_{1}-zeta_center)**2)*L_{n-1}^{lp+1/2}(zeta*(r_{1}-zeta_center)**2)*S_{lp,mp}(zeta_center-r_{1})*sqrt((2*lp+1)/(4*pi)).
  subroutine Laguerre_GTO_projected_on_X_lm(this,norm,zeta,zeta_center,r,&
                                            max_projection_l,l_min_lag,l_max_lag,n_min,n_max,Laguerre_X_nlm)
     use eri_sph_coord_gbl, only: cnla
     use phys_const_gbl, only: fourpi, pi
     use special_functions_gbl, only: cfp_besi, cfp_eval_poly_horner_many
     use lag_cfs_gbl, only: Laguerre_poly_factorization_cfs, l_lag_lim, n_lag_lim
     implicit none
     class(GG_shell_integrals_obj) :: this
     real(kind=cfp), intent(in) :: zeta, zeta_center(3), norm
     real(kind=cfp), allocatable :: r(:)
     integer, intent(in) :: max_projection_l, n_min, n_max, l_min_lag, l_max_lag
     !OUT:
     real(kind=cfp), allocatable :: Laguerre_X_nlm(:,:,:,:) !radial_point,lm,lpmp,n

     integer, parameter :: kode = 2
     real(kind=cfp), parameter :: half = 0.5_cfp

     integer :: n, i, j, err, n_lm, n_lpmp, max_l_aux, l, m, lp, besi_dim, n_lag, lm, nz, l_lag, m_lag, m_ind,&
                lpmp, l_min, l_max, p
     integer :: ll, ml, lpp, mpp, lag_lm, l_aux, m_aux, lambda, mu, lm_aux, n_m1, k, mp, ind
     real(kind=cfp) :: R_dist, R_zeta, arg, exp_arg, tol, val, exp_fac, bessel_fac, exp_part, alp, cf
     real(kind=cfp), allocatable :: Laguerre_aux(:,:), Xlm_CGTO_center(:), cnla_fac(:,:), y(:), d_cfs(:,:), c_lambda(:,:,:),&
                                    transl_cfs(:,:,:), c(:,:,:,:,:), d(:,:,:), P_lag(:,:), tmp(:)
     real(kind=cfp), allocatable :: r1_x_r2(:), x_xx_p(:), outer_poly(:,:)
     logical :: is_zero

        n = size(r)
        n_lm = (max_projection_l+1)**2
        n_lpmp = (l_max_lag+1)**2

        max_l_aux = max_projection_l + l_max_lag + max(n_max-1,0)
        alp = 1.0_cfp/(4.0_cfp*zeta)
        R_dist = sqrt(dot_product(zeta_center,zeta_center))
        R_zeta = zeta*R_dist

        l = max(1,max_projection_l)
        lp = max(1,l_max_lag)
        besi_dim = max_l_aux + 1
        allocate(Laguerre_aux(0:besi_dim,n),cnla_fac(n_min:n_max,0:lp),y(besi_dim),stat=err)
        if (err /= 0) then
           call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_X_lm','Memory allocation 1 failed.',err,1)
        endif

        !Real spherical harmonics for the GTO center: result in the array Xlm_CGTO_center
        call real_harmonics%precalculate_Xlm_for_CGTO_center(zeta_center,max_l_aux,Xlm_CGTO_center)

        tol = F1MACH(4,cfp_dummy)
 
        if (allocated(Laguerre_X_nlm)) then
           if (size(Laguerre_X_nlm,1) /= n .or. size(Laguerre_X_nlm,2) /= n_lm .or.&
               size(Laguerre_X_nlm,3) /= n_lpmp .or. size(Laguerre_X_nlm,4) /= n_max) deallocate(Laguerre_X_nlm)
        endif
 
        if (.not. allocated(Laguerre_X_nlm)) then
           allocate(Laguerre_X_nlm(n,n_lm,n_lpmp,n_min:n_max),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_X_lm','Memory allocation 2 failed.',err,1)
           endif
        endif

        Laguerre_X_nlm = 0.0_cfp

        cnla_fac = 0.0_cfp
        do l_lag=l_min_lag,l_max_lag
           val = sqrt((2*l_lag+1.0_cfp)/fourpi) !get rid of the leading factor sqrt(4pi/(2*l+1)) for the solid harmonics since the NAI formulae are only for r**l*Xlm
           do n_lag=n_min,n_max
              cnla_fac(n_lag,l_lag) = norm*val*fourpi*(-1)**n_lag*cnla(n_lag-1,l_lag,alp)*(-1)**l_lag
           enddo
        enddo

        !Evaluate the partial wave projections of the GTO for the partial waves up to max_l_aux as needed later to perform the translations of the polynomials.
        Laguerre_aux = 0.0_cfp
        do i=1,n
           arg = 2.0_cfp*r(i)*R_zeta
           exp_arg = zeta*(r(i)-R_dist)**2                 
           exp_fac = exp(-exp_arg)

           if (arg .le. tol) then !Evaluate limit for arg -> 0
              l = 0
              bessel_fac = 1.0_cfp !only l=0 Bessel function is non-zero
              exp_part = bessel_fac*exp_fac
              Laguerre_aux(l,i) = exp_part
           else
              call cfp_besi(arg, half, kode, besi_dim, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
              val = sqrt(pi/(2.0_cfp*arg))
              do l=0,max_l_aux
                 bessel_fac = y(l+1)*val
                 exp_part = bessel_fac*exp_fac
                 Laguerre_aux(l,i) = exp_part
              enddo !l
           endif
        enddo !i

        if (n_max .eq. 1) then !L_{0}^{lp+1/2}(x) = 1: the polynomial is trivial

           if (l_max_lag .eq. 0) then 
              !s-type Laguerre GTO: projection of exp(-zeta*(r_{1}-zeta_center)**2) so no translation of S_lpmp needed.

              do l=0,max_projection_l
                 do m=-l,l
                    lm = l*l+l+m+1
                    if (Xlm_CGTO_center(lm) /= 0.0_cfp) then
                       Laguerre_X_nlm(1:n,lm,1,1) = Laguerre_aux(l,1:n)*Xlm_CGTO_center(lm)*cnla_fac(1,0)
                    endif
                 enddo !m
              enddo !l

           else 
              !Laguerre GTO with L > 0: we need to translate the S_lpmp solid harmonic part.

              !Projections for a standard GTO with arbitrary lp,mp so we need to translate the solid harmonic S_{lp,mp}.
              allocate(d_cfs(n,1:l_max_lag+1),stat=err)
              if (err /= 0) then
                 call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_X_lm','Memory allocation 3 failed.',err,1)
              endif

              do l_lag=l_min_lag,l_max_lag
    
                 !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
                 call this%precalculate_solh_translation_coeffs(l_lag,R_dist,Xlm_CGTO_center,transl_cfs)
    
                 do l=0,max_projection_l
                    do m=-l,l
                       lm = l*l+l+m+1
    
                       !Calculate the coupling coefficients needed to evaluate the projection using SCE of GTO.
                       !call this%calculate_lambda_couplings(l_lag,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
                       call calculate_lambda_couplings(l_lag,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
         
                       do m_lag=-l_lag,l_lag
                          m_ind = m_lag+l_lag+1
                          lpmp = l_lag*l_lag+l_lag+m_lag+1

                          !Skip the symmetry-forbidden projections
                          is_zero = .true.
                          do lambda=0,l+l_lag
                             do lp=0,l_lag
                                if (c_lambda(lambda,lp,m_ind) /= 0.0_cfp) then
                                   is_zero = .false.
                                   exit
                                endif
                             enddo !lp
                          enddo !lambda
                          if (is_zero) cycle

                           !Evaluate the polynomial coefficients
                           do lp=0,l_lag
                              l_min = abs(l-lp)
                              l_max = l+lp
                              j = lp+1
                              do i=1,n
                                 d_cfs(i,j) = sum(c_lambda(l_min:l_max,lp,m_ind)*Laguerre_aux(l_min:l_max,i))
                              enddo
                           enddo !lp

                           !Evaluate the polynomial for each r(i)
                           d_cfs = d_cfs*cnla_fac(1,l_lag)
                           call cfp_eval_poly_horner_many(l_lag,r,n,d_cfs,Laguerre_X_nlm(1:n,lm,lpmp,1))

                       enddo !m_lag
         
                    enddo !m
                 enddo !l
    
              enddo !l_lag

           endif

        else !n_max /= 1: L_{0}^{lp+1/2}(x) is non-trivial

           if (l_max_lag > l_lag_lim .or. n_max-1 > n_lag_lim) then
              call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_X_lm',&
                          'Need Laguerre coefficients which have not been generated. See module lag_cfs.',3,1)
           endif

           !Get the coefficients needed to factorize L(x+y) (x = scalar, y = dot_product(r1,r2)) in terms of the real spherical harmonics.
           call Laguerre_poly_factorization_cfs(c)

           if (l_max_lag .eq. 0) then 
              !s-type Laguerre GTO so we need to translate only the Laguerre polynomial. Does this ever happen???
              !The implementation of this case is easy: see the calculation of P_lag in the branch below.
              call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_X_lm',&
                          'Not implemented Laguerre translation for the case of l_max_lag .eq. 0.',2,1)

            else

              !Laguerre GTO with L > 0: we need to translate both the Laguerre polynomial and the S_lpmp part.
              n_m1 = n_max-1
              i = (n_m1 + l_max_lag +1)**2
              j = max(l_max_lag+1,n_max)
              allocate(P_lag(1:n,0:n_max-1),tmp(n),r1_X_r2(n),x_xx_p(n),d(i,0:n_m1,0:l_max_lag),d_cfs(n,1:j),&
                       outer_poly(n,1:n_max),stat=err)
              if (err /= 0) then
                 call xermsg('GG_shell_mixed_integrals_mod','Laguerre_GTO_projected_on_X_lm','Memory allocation 4 failed.',err,1)
              endif
 
              do i=1,n
                 x_xx_p(i) = zeta*(r(i)**2+R_dist**2)
                 r1_X_r2(i) = 2*zeta*r(i)*R_dist
              enddo !i

              do l_lag=l_min_lag,l_max_lag

                 !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
                 call this%precalculate_solh_translation_coeffs(l_lag,R_dist,Xlm_CGTO_center,transl_cfs)

                 do n_lag=n_min,n_max

                    !Ensure that phi_n-1,lm^{b} is taken:
                    n_m1 = n_lag-1

                    if (n_m1 .eq. 0) then 
                       !Laguerre polynomial for n-1 == 0 requested; this is 1;
                       !the factor of 4pi compensates for the product X_00(r)*X_00(zeta_center)=1/(4*pi) coming from the expansion of L_{0}^{l+1/2}(x+y)
                       P_lag(1:n,0) = 1.0_cfp*cnla_fac(n_lag,l_lag)*fourpi
                    else
                       do lp=0,n_m1
                          !Generate P_{lp}^{n_lag-1,l_lag}*(-1)**lp*cnla_fac(n_lag,l_lag)
                          !Use the Horner scheme twice to evaluate the 2D polynomial:
                          !P_{lp}^{n_lag-1,l_lag} = sum_{p} (zeta*(r(i)**2+R_dist**2))**p * outer_poly_{p,r(i)}; outer_poly_{p,r(i)} = sum_{k} c(k,p,lp,l_lag,n_m1)* (2*zeta*r(i)*R_dist)
                          do p = 0,n_m1
                             do k=0,n_m1
                                d_cfs(1:n,k+1) = c(k,p,lp,l_lag,n_m1)
                             enddo !k
                             !Evaluate the inner polynomial in (2*zeta*r(i)*R_dist) for each value of r(i): the result are coefficients for the outer polynomial
                             call cfp_eval_poly_horner_many(n_m1,r1_X_r2,n,d_cfs,outer_poly(1:n,p+1))
                          enddo !p
                          !Evaluate the outer polynomial in (zeta*(r(i)**2+R_dist**2)) for each value of r(i):
                          call cfp_eval_poly_horner_many(n_m1,x_xx_p,n,outer_poly,P_lag(1:n,lp))
                          P_lag(1:n,lp) = P_lag(1:n,lp)*cnla_fac(n_lag,l_lag)*(-1)**lp
                       enddo !lp
                    endif

                    do ml=-l_lag,l_lag
                       m_ind = ml+l_lag+1
                       lag_lm = l_lag*l_lag+m_ind

                       !Couple the spherical harmonics from the Laguerre translation with the ones from the solid harmonic translation.
                       d = 0.0_cfp
                       do lpp=0,l_lag
                          do lp=0,n_m1
                             do l_aux=abs(lp-lpp),lp+lpp
                                do m_aux=-l_aux,l_aux
                                   lm_aux = l_aux*l_aux+m_aux+l_aux+1
   
                                   do mp=-lp,lp
                                      ind = lp*lp+lp+mp+1
                                      do mpp=-lpp,lpp
                                         d(lm_aux,lp,lpp) = d(lm_aux,lp,lpp)&
                                                          + cpl%rgaunt(l_aux,lp,lpp,m_aux,mp,mpp)&
                                                          * Xlm_CGTO_center(ind)&
                                                          * transl_cfs(mpp+lpp+1,lpp,m_ind)
                                      enddo !mpp
                                   enddo !mp
   
                                enddo !M
                             enddo !L
                          enddo !lp
                       enddo !lpp
                    
                       do l=0,max_projection_l
                          do m=-l,l
                             lm = l*l+l+m+1
   
                             do lpp=0,l_lag

                                !Sum over lp and multiply in the factor P coming from the translation of the Laguerre polynomial.
                                d_cfs(1:n,lpp+1) = 0.0_cfp
                                do lp=0,n_m1
                                   
                                   !Couple to the spherical harmonics of the projection.
                                   tmp = 0.0_cfp
                                   do l_aux=abs(lp-lpp),lp+lpp
                                      do m_aux=-l_aux,l_aux
                                         lm_aux = l_aux*l_aux+m_aux+l_aux+1
                                         do lambda=abs(l-l_aux),l+l_aux
                                            do mu=-lambda,lambda
                                               ind = lambda*lambda+lambda+mu+1
                                               cf = cpl%rgaunt(lambda,l_aux,l,mu,m_aux,m)
                                               if (cf /= 0.0_cfp) then
                                                  tmp(1:n) = tmp(1:n)&
                                                           + cf&
                                                           * Laguerre_aux(lambda,1:n)&
                                                           * Xlm_CGTO_center(ind)&
                                                           * d(lm_aux,lp,lpp)
                                               endif
                                            enddo !mu
                                         enddo !lambda
                                      enddo !m_aux
                                   enddo !l_aux

                                   d_cfs(1:n,lpp+1) = d_cfs(1:n,lpp+1) + tmp(1:n)*P_lag(1:n,lp)

                                enddo !lp
                             enddo !lpp

                             !Evaluate the polynomial for each r(i)
                             call cfp_eval_poly_horner_many(l_lag,r,n,d_cfs,Laguerre_X_nlm(1:n,lm,lag_lm,n_lag))
   
                          enddo !m
                       enddo !l

                    enddo !ml

                 enddo !n_lag

              enddo !l_lag

           endif !l_max_lag .eq. 0

        endif !n_max .eq. 1

  end subroutine Laguerre_GTO_projected_on_X_lm

  !> Assumes Laguerre_BG_Xlm already allocated
  subroutine r1_integrate_Laguerre_GTO_BG(grid,Laguerre_G_X_nlm,k,n_k,n_min,n_max,Laguerre_BG_Xlm)
     implicit none
     class(grid_r1_r2_obj) :: grid
     real(kind=cfp), allocatable :: Laguerre_G_X_nlm(:,:,:,:,:), Laguerre_BG_Xlm(:,:,:,:)
     integer, intent(in) :: k,n_min,n_max,n_k

     integer :: d2,d3,lm,lpmp,n,i,ind,s_point,e_point,err,d1,base,n_mg,mg

        n_mg = size(Laguerre_G_X_nlm,2)
        d2 = size(Laguerre_G_X_nlm,3)
        d3 = size(Laguerre_G_X_nlm,4)

        if (.not.(allocated(Laguerre_BG_Xlm))) then
           d1 = grid%last_bspline_inner*d2*n_mg
           allocate(Laguerre_BG_Xlm(d1,n_k,d3,n_min:n_max),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','r1_integrate_Laguerre_GTO_BB','Memory allocation failed.',err,1)
           endif
        endif

        base = grid%last_bspline_inner*d2

        !Integrate Laguerre_G_X_nlm over all B-splines.
        !todo swap the loops so the one over i is the outer most one
        do n=n_min,n_max
           do lpmp=1,d3
              do lm=1,d2
                 do mg=1,n_mg
                    do i=1,grid%last_bspline_inner
                       s_point = grid%bspline_start_end_r1(1,i)
                       e_point = grid%bspline_start_end_r1(2,i)
                       !We put the 2D array of pairs and lm indices into a 1D
                       !array with the view that the BBGG and BGGG algorithms
                       !are generic and do not need to know what is the nature
                       !of the dummy indices that are being looped over.
                       ind = i + grid%last_bspline_inner*(lm-1) + base*(mg-1)
                       Laguerre_BG_Xlm(ind,k,lpmp,n) = sum(grid%B_vals_r1(s_point:e_point,i)&
                                                          *Laguerre_G_X_nlm(s_point:e_point,mg,lm,lpmp,n)&
                                                          *grid%r1(s_point:e_point)&
                                                          *grid%w1(s_point:e_point))
                    enddo !i
                 enddo !mg
              enddo !lm
           enddo !lpmp
        enddo !d4
 
  end subroutine r1_integrate_Laguerre_GTO_BG

  !> Assumes Laguerre_BG_Xlm already allocated
  !> Routine assuming the CGTO lies on the CMS.
  subroutine r1_integrate_Laguerre_GTO_B_cms_G(grid,Laguerre_X_nlm,k,n_k,n_min,n_max,cgto_on_grid,Laguerre_BG_Xlm)
     implicit none
     class(grid_r1_r2_obj) :: grid
     real(kind=cfp), allocatable :: Laguerre_X_nlm(:,:,:,:), Laguerre_BG_Xlm(:,:,:,:), cgto_on_grid(:)
     integer, intent(in) :: k,n_min,n_max,n_k

     integer :: d2,d3,lm,lpmp,n,pair,ind,s_point,e_point,err,d1,i

        d2 = size(Laguerre_X_nlm,2)
        d3 = size(Laguerre_X_nlm,3)

        if (.not.(allocated(Laguerre_BG_Xlm))) then
           d1 = grid%last_bspline_inner*d2
           allocate(Laguerre_BG_Xlm(d1,n_k,d3,n_min:n_max),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','r1_integrate_Laguerre_GTO_B_cms_G','Memory allocation failed.',err,1)
           endif
        endif

        !Integrate Laguerre_X_nlm over unique pairs of B-splines.
        do n=n_min,n_max
           do lpmp=1,d3
              do lm=1,d2
                 do i=1,grid%last_bspline_inner
                    s_point = grid%bspline_start_end_r1(1,i)
                    e_point = grid%bspline_start_end_r1(2,i)
                    !We put the 2D array of pairs and lm indices into a 1D
                    !array with the view that the BBGG and BGGG algorithms
                    !are generic and do not need to know what is the nature
                    !of the dummy indices that are being looped over.
                    ind = i + grid%last_bspline_inner*(lm-1)
                    Laguerre_BG_Xlm(ind,k,lpmp,n) = sum(grid%B_vals_r1(s_point:e_point,i)&
                                                    *cgto_on_grid(s_point:e_point)&
                                                    *Laguerre_X_nlm(s_point:e_point,lm,lpmp,n)&
                                                    *grid%r1(s_point:e_point)&
                                                    *grid%w1(s_point:e_point))
                 enddo !i
              enddo !lm
           enddo !lpmp
        enddo !d4
 
  end subroutine r1_integrate_Laguerre_GTO_B_cms_G

  !> Assumes Laguerre_BB_Xlm already allocated
  subroutine r1_integrate_Laguerre_GTO_BB(grid,Laguerre_X_nlm,k,n_k,n_min,n_max,Laguerre_BB_Xlm)
     implicit none
     class(grid_r1_r2_obj) :: grid
     real(kind=cfp), allocatable :: Laguerre_X_nlm(:,:,:,:), Laguerre_BB_Xlm(:,:,:,:)
     integer, intent(in) :: k,n_min,n_max,n_k

     integer :: d2,d3,lm,lpmp,n,pair,ind,s_point,e_point,err,d1

        d2 = size(Laguerre_X_nlm,2)
        d3 = size(Laguerre_X_nlm,3)

        if (.not.(allocated(Laguerre_BB_Xlm))) then
           d1 = grid%n_unique_pairs*d2
           allocate(Laguerre_BB_Xlm(d1,n_k,d3,n_min:n_max),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','r1_integrate_Laguerre_GTO_BB','Memory allocation failed.',err,1)
           endif
        endif

        !Integrate Laguerre_X_nlm over unique pairs of B-splines.
        do n=n_min,n_max
           do lpmp=1,d3
              do lm=1,d2
                 do pair=1,grid%n_unique_pairs
                    s_point = grid%BB_start_end_r1(1,pair)
                    e_point = grid%BB_start_end_r1(2,pair)
                    !We put the 2D array of pairs and lm indices into a 1D
                    !array with the view that the BBGG and BGGG algorithms
                    !are generic and do not need to know what is the nature
                    !of the dummy indices that are being looped over.
                    ind = pair + grid%n_unique_pairs*(lm-1)
                    Laguerre_BB_Xlm(ind,k,lpmp,n) = sum(grid%BB_vals_r1(s_point:e_point,pair)&
                                                       *Laguerre_X_nlm(s_point:e_point,lm,lpmp,n)&
                                                       *grid%w1(s_point:e_point))
                 enddo !pair
              enddo !lm
           enddo !lpmp
        enddo !d4
 
  end subroutine r1_integrate_Laguerre_GTO_BB

  !> Assumes Boys_BX_Xlm already allocated
  !> Routine assuming the CGTO lies on the CMS.
  subroutine r1_integrate_Boys_B_cms_G(grid,Boys_Xlm,k,n_k,cgto_on_grid,Boys_BG_Xlm)
     implicit none
     class(grid_r1_r2_obj) :: grid
     real(kind=cfp), allocatable :: Boys_Xlm(:,:,:), Boys_BG_Xlm(:,:,:), cgto_on_grid(:)
     integer, intent(in) :: k,n_k

     integer :: d2,d3,lm,LbMb,ind,s_point,e_point,d1,err,i

        d2 = size(Boys_Xlm,2)
        d3 = size(Boys_Xlm,3)

        if (.not.(allocated(Boys_BG_Xlm))) then
           d1 = grid%last_bspline_inner*d2
           allocate(Boys_BG_Xlm(d1,n_k,d3),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','r1_integrate_Boys_B_cms_G','Memory allocation 1 failed.',err,1)
           endif
        endif

        !Integrate Boys_Xlm over unique the B-splines and the radial CGTO.
        do LbMb=1,d3
           do lm=1,d2
              do i=1,grid%last_bspline_inner
                 s_point = grid%bspline_start_end_r1(1,i)
                 e_point = grid%bspline_start_end_r1(2,i)
                 !We put the 2D array of pairs and lm indices into a 1D
                 !array with the view that the BBGG and BGGG algorithms
                 !are generic and do not need to know what is the nature
                 !of the dummy indices that are being looped over.
                 ind = i + grid%last_bspline_inner*(lm-1)
                 Boys_BG_Xlm(ind,k,LbMb) = sum(grid%B_vals_r1(s_point:e_point,i)&
                                              *cgto_on_grid(s_point:e_point)&
                                              *Boys_Xlm(s_point:e_point,lm,LbMb)&
                                              *grid%r1(s_point:e_point)&
                                              *grid%w1(s_point:e_point))
              enddo !pair
           enddo !lm
        enddo !LbMb
 
  end subroutine r1_integrate_Boys_B_cms_G

  !> Assumes Boys_BX_Xlm already allocated
  subroutine r1_integrate_Boys_BG(grid,Boys_G_Xlm,k,n_k,Boys_BG_Xlm)
     implicit none
     class(grid_r1_r2_obj) :: grid
     real(kind=cfp), allocatable :: Boys_G_Xlm(:,:,:,:), Boys_BG_Xlm(:,:,:)
     integer, intent(in) :: k,n_k

     integer :: d2,d3,lm,LbMb,ind,s_point,e_point,d1,err,n_mg,base,mg,i

        n_mg = size(Boys_G_Xlm,2)
        d2 = size(Boys_G_Xlm,3)
        d3 = size(Boys_G_Xlm,4)

        if (.not.(allocated(Boys_BG_Xlm))) then
           d1 = grid%last_bspline_inner*d2*n_mg
           allocate(Boys_BG_Xlm(d1,n_k,d3),stat=err)
           if (err /= 0) then
              call xermsg('GG_shell_mixed_integrals_mod','r1_integrate_Boys_BG','Memory allocation failed.',err,1)
           endif
        endif

        base = grid%last_bspline_inner*d2

        !Integrate Boys_G_Xlm over B-splines and radial parts of the CGTO.
        !todo swap the loops so the one over i is the outer most one
        do LbMb=1,d3
           do lm=1,d2
              do mg=1,n_mg
                 do i=1,grid%last_bspline_inner
                    s_point = grid%bspline_start_end_r1(1,i)
                    e_point = grid%bspline_start_end_r1(2,i)
                    !We put the 3D array of B-spline indices, lm indices and CGTO m-values into a 1D
                    !array in the order (bspline_index,bspline_lm,CGTO_B_m) with the view that the BBGG and BGGG algorithms
                    !are generic and do not need to know what is the nature of the dummy indices that are being looped over.
                    ind = i + grid%last_bspline_inner*(lm-1) + base*(mg-1)
                    Boys_BG_Xlm(ind,k,LbMb) = sum(grid%B_vals_r1(s_point:e_point,i)&
                                                 *Boys_G_Xlm(s_point:e_point,mg,lm,LbMb)&
                                                 *grid%r1(s_point:e_point)&
                                                 *grid%w1(s_point:e_point))
                 enddo !i
              enddo !mg
           enddo !lm
        enddo !LbMb

  end subroutine r1_integrate_Boys_BG

  !> Assumes Boys_BX_Xlm already allocated
  subroutine r1_integrate_Boys_BB(grid,Boys_Xlm,k,n_k,Boys_BB_Xlm)
     implicit none
     class(grid_r1_r2_obj) :: grid
     real(kind=cfp), allocatable :: Boys_Xlm(:,:,:), Boys_BB_Xlm(:,:,:)
     integer, intent(in) :: k,n_k

     integer :: d2,d3,lm,LbMb,pair,ind,s_point,e_point,d1,err

        d2 = size(Boys_Xlm,2)
        d3 = size(Boys_Xlm,3)

        if (.not.(allocated(Boys_BB_Xlm))) then
           d1 = grid%n_unique_pairs*d2
           allocate(Boys_BB_Xlm(d1,n_k,d3),stat=err)
           if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','r1_integrate_Boys_BB','Memory allocation failed.',err,1)
        endif

        !Integrate Boys_Xlm over unique pairs of B-splines.
        do LbMb=1,d3
           do lm=1,d2
              do pair=1,grid%n_unique_pairs
                 s_point = grid%BB_start_end_r1(1,pair)
                 e_point = grid%BB_start_end_r1(2,pair)
                 !We put the 2D array of pairs and lm indices into a 1D
                 !array with the view that the BBGG and BGGG algorithms
                 !are generic and do not need to know what is the nature
                 !of the dummy indices that are being looped over.
                 ind = pair + grid%n_unique_pairs*(lm-1)
                 Boys_BB_Xlm(ind,k,LbMb) = sum(grid%BB_vals_r1(s_point:e_point,pair)&
                                              *Boys_Xlm(s_point:e_point,lm,LbMb)&
                                              *grid%w1(s_point:e_point))
              enddo !pair
           enddo !lm
        enddo !LbMb
 
  end subroutine r1_integrate_Boys_BB

  !> Evaluates the radial part of a CGTO which is assumed to be sitting on the CMS on the grid of radial points.
  subroutine CMS_CGTO_on_grid(r,n,cgto_shell,cgto_on_grid)
     use phys_const_gbl, only: fourpi
     implicit none
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: r(n)
     type(CGTO_shell_data_obj), intent(in) :: cgto_shell
     real(kind=cfp), allocatable :: cgto_on_grid(:)

     integer :: i,j,err
     real(kind=cfp) :: fac
     real(kind=cfp), allocatable :: cgto_solh(:)

        if (allocated(cgto_on_grid)) deallocate(cgto_on_grid)

        allocate(cgto_solh(n),cgto_on_grid(n),stat=err)
        if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','CMS_CGTO_on_grid','Memory allocation failed.',err,1)

        !Calculate the radial part of the CGTO on the r1 grid.
        fac = sqrt(fourpi/(2*cgto_shell%l+1.0_cfp))
        do i=1,n
           cgto_solh(i) = fac*r(i)**cgto_shell%l
        enddo

        cgto_on_grid = 0.0_cfp
        do j=1,cgto_shell%number_of_primitives
           fac = cgto_shell%norm*cgto_shell%norms(j)*cgto_shell%contractions(j)
           do i=1,n
              cgto_on_grid(i) = cgto_on_grid(i) + fac*exp(-cgto_shell%exponents(j)*r(i)**2)*cgto_solh(i)
           enddo
        enddo

  end subroutine CMS_CGTO_on_grid

  !> Calculates the angular integrals \int d\Omega(r_{1}) X_{l,m}(r_{1}) CGTO_{Lg,Mg}(r_{1}-A)*norm*C_{-1,L}^{a}F_{L}(zeta*(r_{1}-zeta_center)**2)*S_{L,M}(r_{1}-zeta_center)*sqrt((2*L+1)/(4*pi))
  !> for the grid of r1 points, l,m projections up to max_projection_l, L,M up to max_boys_l and all Mg values corresponding to Lg.
  !> F_{L}(x) = \int_{0}^{1} du u**2*exp(-u**2*x)
  subroutine Boys_projected_on_G_X_lm(this,u_grid,norm,zeta,zeta_center,cgto_shell,r,max_projection_l,max_boys_l,Boys_G_Xlm)
     use phys_const_gbl, only: fourpi, pi
     use special_functions_gbl, only: cfp_besi, cfp_eval_poly_horner_many
     use const_gbl, only: epsabs
     implicit none
     class(GG_shell_integrals_obj) :: this
     type(pw_expansion_obj) :: u_grid
     real(kind=cfp), intent(in) :: norm, zeta, zeta_center(3)
     real(kind=cfp), allocatable :: r(:)
     integer, intent(in) :: max_projection_l,max_boys_l
     type(CGTO_shell_data_obj) :: cgto_shell
     !OUT:
     real(kind=cfp), allocatable :: Boys_G_Xlm(:,:,:,:) !Boys_G_Xlm(radial point,CGTO m-value,lm,LbMb)

     integer, parameter :: kode = 2, n_u_points = 18
     real(kind=cfp), parameter :: half = 0.5_cfp

     integer :: err, n, n_lm, LbMb, lm, i, u, Lboys, Mboys, max_l_aux, besi_dim, l, lp, mp, lpp, mpp, j, nz, m, n_lm_boys, m_ind,&
                lap, map, lambda, CGTO_L, l_min,l_max, n_mg, base, n_lm_aux, ind, Mg
     integer :: lppmpp, lpmp, sum_l, n_lm_sum, la_p_lb, u_lim, u_lim_saved
     real(kind=cfp) :: alp, cnla(0:max_boys_l+1), R_dist, bessel_fac, u_sq, val, u_sq_zeta, fac, prec, max_u_fac
     real(kind=cfp) :: preexp_factor, cf, Rab(3), Rab_sq, prod_center(3), contr, RC
     real(kind=cfp), allocatable :: y(:), u_fac(:,:), d_cfs(:,:), Boys_aux(:,:,:,:), Xlm_CGTO_center(:), transl_cfs(:,:,:), &
                                    c_lambda(:,:,:), exp_part(:), prod_exponent(:), R_distprod(:), contr_fac(:)
     real(kind=cfp), allocatable :: Xlm_product_center(:), R_aux(:,:), transl_cfs_all(:,:,:,:), Xlm_CGTO_C_center(:), &
                                    transl_cfs_all_AB(:,:,:,:,:), transl_cfs_AB(:,:,:,:,:), arg_u(:), exp_fac(:)
     logical :: non_zero, is_zero
     real(kind=cfp), save :: tol = 0.0_cfp

     !debug
     !real(kind=cfp) :: term1, term2, term3, interterm, cnorm, second_der, der_tol, arg_0, correction, arg_a, arg_b, u_min
     real(kind=cfp), parameter :: exponent_del_thrs = 30.0_cfp !todo move to the parameters
     integer, parameter :: bisection_iterations = 15 !todo move to the parameters
     !real(kind=cfp) :: exponent_0, exponent_1, exponent_u, u_grid_max, u_1, u_2, u_sq_loc, u_sq_zeta_loc, max_u_fac_loc, contr_loc
     real(kind=cfp) :: exponent_0, exponent_1, exponent_u, u_grid_max, u_sq_loc, u_sq_zeta_loc, max_u_fac_loc, contr_loc
     real(kind=cfp) :: u_grid_cut
     integer :: u_int_grid, k, u_lim_loc, u_int_cut
     type(pw_expansion_obj) :: u_grid_loc, u_grid_temp
     real(kind=cfp), allocatable :: u_fac_loc(:,:), Xlm_CGTO_center_loc(:), &
                                    prod_exponent_loc(:), R_distprod_loc(:), contr_fac_loc(:)
     real(kind=cfp), allocatable :: Xlm_product_center_loc(:), arg_u_loc(:), exp_fac_loc(:)
     real(kind=cfp) :: prod_center_loc(3), preexp_factor_loc

       n = size(r)
       n_mg = cgto_shell%number_of_functions
       n_lm = (max_projection_l+1)**2
       n_lm_boys = (max_boys_l+1)**2

       if (tol .eq. 0.0_cfp) then
          tol = F1MACH(4,cfp_dummy)
       endif

       if (allocated(Boys_G_Xlm)) then
          if (size(Boys_G_Xlm,1) /= n .or. size(Boys_G_Xlm,2) /= n_mg .or. size(Boys_G_Xlm,4) /= n_lm_boys .or. &
               size(Boys_G_Xlm,3) /= n_lm) deallocate(Boys_G_Xlm)
       endif

       if (.not. allocated(Boys_G_Xlm)) then
          allocate(Boys_G_Xlm(n,n_mg,n_lm,n_lm_boys),stat=err)
          if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_G_X_lm','Memory allocation failed.',1,1)
       endif

       Boys_G_Xlm(:,:,:,:) = 0.0_cfp

       alp = 1.0_cfp/(4.0_cfp*zeta)
       do Lboys=0,max_boys_l
          val = sqrt((2*Lboys+1.0_cfp)/fourpi) !get rid of the leading factor sqrt(4pi/(2*l+1)) for the solid harmonics since the NAI formulae are only for r**l*Xlm
          cnla(Lboys) = val*norm*fourpi/(2.0_cfp*alp)**(Lboys+0.5_cfp)*(-1)**Lboys
       enddo

       max_l_aux = max_projection_l + max_boys_l + cgto_shell%l

       besi_dim = max_l_aux+1
       i = max(max_l_aux,1)
       j = max(max_boys_l,1)
       n_lm_aux = (max_l_aux+1)**2
       !todo allocate only once and use the save attribute
       allocate(y(besi_dim),u_fac(0:j,u_grid%n_total_points),Boys_aux(n_lm_aux,n,0:j,cgto_shell%number_of_primitives),&
                Xlm_product_center(u_grid%n_total_points*n_lm_aux),exp_part(0:i),prod_exponent(u_grid%n_total_points),&
                R_distprod(u_grid%n_total_points),contr_fac(u_grid%n_total_points),&
                arg_u(u_grid%n_total_points),exp_fac(u_grid%n_total_points),stat=err)
       if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_G_X_lm','Memory allocation 2 failed.',err,1)

       !todo allocate only once and use the save attribute
       allocate(u_fac_loc(0:j,u_grid%n_total_points),&
                Xlm_product_center_loc(u_grid%n_total_points*n_lm_aux),prod_exponent_loc(u_grid%n_total_points),&
                R_distprod_loc(u_grid%n_total_points),contr_fac_loc(u_grid%n_total_points),&
                arg_u_loc(u_grid%n_total_points),exp_fac_loc(u_grid%n_total_points),stat=err)
       if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_G_X_lm','Memory allocation 3 failed.',err,1)
       val = 20.0_cfp
       call u_grid_temp%eval_exponential_grid(val,28,.false.)

       do u=1,u_grid%n_total_points
          do L=0,max_boys_l
             u_fac(L,u) = u_grid%r_points(u)**(2*L)*u_grid%weights(u)
          enddo
       enddo

       max_u_fac = maxval(u_fac)*maxval(abs(cnla))

       Rab(1:3) = zeta_center(1:3)-cgto_shell%center(1:3)
       Rab_sq = dot_product(Rab,Rab)

       prec = tol*epsabs

       !Note on the contraction: in principle the contraction could be done
       !directly in the loop over j,i.e. summing into the array Boys_aux.
       !However, this could potentially lead to loss of precision due to
       !accumulation of small contributions to possibly large values stored in
       !Boys_aux. Therefore the contraction is performed at the end of the
       !j-loop.
       Boys_aux = 0.0_cfp
       do j=1,cgto_shell%number_of_primitives

          contr = cgto_shell%norm*cgto_shell%norms(j)*cgto_shell%contractions(j)

          !Calculate the product GTOs for each u: exp(-u^2*zeta(r-zeta_center)**2)*exp(-alp_{j}*(r-A)**2)
          i = 0
          u_lim = u_grid%n_total_points
          do u=1,u_grid%n_total_points
             u_sq = u_grid%r_points(u)**2
             u_sq_zeta = u_sq*zeta
   
             prod_exponent(u) = u_sq_zeta+cgto_shell%exponents(j)
             prod_center(1:3) = (u_sq_zeta*zeta_center(1:3) + cgto_shell%exponents(j)*cgto_shell%center(1:3))/prod_exponent(u)
             preexp_factor = exp(-u_sq_zeta*cgto_shell%exponents(j)/prod_exponent(u)*Rab_sq)
   
             R_distprod(u) = sqrt(dot_product(prod_center(1:3),prod_center(1:3)))
   
             !The preexponential factor coming from the product of the Boys GTO and the CGTO and the normalization factors
             contr_fac(u) = contr*preexp_factor

             !Real spherical harmonics for the product center: result in the array Xlm_CGTO_center
             call real_harmonics%precalculate_Xlm_for_CGTO_center(prod_center(1:3),max_l_aux,Xlm_CGTO_center)
             Xlm_product_center(i+1:i+n_lm_aux) = Xlm_CGTO_center(1:n_lm_aux)
             i = i + n_lm_aux

             !Find the u for which the prefactor is too small.
             if (u_sq**(max_boys_l)*preexp_factor < prec) then
                if (u_lim .eq. u_grid%n_total_points) u_lim = u
             endif
          enddo

          !r1 points:
          u_lim_saved = u_lim
          do i=1,n

             !Integral over u:
             !todo for very large zeta exponent the potential of the GG pair is
             !essentially a point charge: the integration using a fixed u-grid may be inaccurate!!!
             !The same happens for very large r(i). This should be dealt with
             !zooming onto the u-interval where the integrand is non-zero,
             !rescaling it onto [0;1] and integrating then. This would require
             !recalculating all u-dependent values.

             !here we decide, if we use the u_grid given on input or some other one
             !first: no regard to the possible point of unsmoothness, if the exponent of exp_fac in u=0 minus exponent of exp_fac in
             !u=1 is smaller than 23 (approx log(1E+10)) then use the u_grid on input, else find value of u_min where the difference of
             !exponents is approximately 23, then construct a new quadrature grid on the interval from 0 to u_min otherwise the
             !quadrature method will be the same

             ! first determine whether to integrate on whole (0,1) or on a subinterval
             exponent_0 = exp_fac_exponent(0.0_cfp,zeta,cgto_shell%exponents(j),zeta_center,cgto_shell%center,r(i))
             exponent_1 = exp_fac_exponent(1.0_cfp,zeta,cgto_shell%exponents(j),zeta_center,cgto_shell%center,r(i))
             if (exponent_0 - exponent_1 < exponent_del_thrs) then
                u_int_grid = 1 ! use the u_grid given on input
                
                u_grid_max = 1.0_cfp
             else ! generate new u_grid
                u_int_grid = 2

                u_grid_max = 0.5_cfp
                exponent_u = exp_fac_exponent(u_grid_max,zeta,cgto_shell%exponents(j),zeta_center,cgto_shell%center,r(i))
                do k=1,bisection_iterations
                   if (exponent_0 - exponent_u < exponent_del_thrs) then
                      u_grid_max = u_grid_max + 0.5_cfp**(k+1)
                   else
                      u_grid_max = u_grid_max - 0.5_cfp**(k+1)
                   endif
                   exponent_u = exp_fac_exponent(u_grid_max,zeta,cgto_shell%exponents(j),zeta_center,cgto_shell%center,r(i))
                enddo
             endif

             ! second, determine whether to cut the quadrature in two
             ! i.e. determine, whether the argument has a point of unsmoothness, and whether it is within the quadrature interval
             u_grid_cut = arg_min(zeta,cgto_shell%exponents(j),zeta_center,cgto_shell%center)
             u_int_cut = 1
             if (0 < u_grid_cut .and. u_grid_cut < u_grid_max) then
                u_int_grid = 2
                u_int_cut = 2
             endif

             !finally, generate the quadrature, if requested
             if (u_int_grid == 2) then
                u_grid_loc = u_grid
                if (u_int_cut == 1) then
                   u_grid_loc%r_points(:) = u_grid_loc%r_points(:) * u_grid_max
                   u_grid_loc%weights(:) = u_grid_loc%weights(:) * u_grid_max
                elseif(u_int_cut == 2) then
                   u_grid_loc%r_points(1:u_grid_temp%n_total_points) = u_grid_temp%r_points(:) * u_grid_cut
                   u_grid_loc%weights(1:u_grid_temp%n_total_points) = u_grid_temp%weights(:) * u_grid_cut
                   u_grid_loc%r_points(u_grid_temp%n_total_points+1:u_grid_loc%n_total_points) = & 
                      u_grid_temp%r_points(:) * (u_grid_max-u_grid_cut) + u_grid_cut
                   u_grid_loc%weights(u_grid_temp%n_total_points+1:u_grid_loc%n_total_points) = &
                      u_grid_temp%weights(:) * (u_grid_max-u_grid_cut)
                endif
             endif
             
             !Precalculate the arguments needed for the integration
             do u=1,u_grid%n_total_points
                arg_u(u) = 2.0_cfp*r(i)*R_distprod(u)*prod_exponent(u)
                exp_fac(u) = contr_fac(u)*exp(-prod_exponent(u)*(r(i)-R_distprod(u))**2)
             enddo !u

             !Update the value of u_lim based on the value of exp_fac:
             u_lim = u_lim_saved
             do u=u_lim_saved,u_grid%n_total_points
                if (abs(exp_fac(u)) > prec) u_lim=u
             enddo !u

             !skip those r-points for which the u-integral would definitely be too small.
             if (abs(max_u_fac)*maxval(abs(exp_fac)) < prec) then
                !print *,'skipping',r(i)
                cycle
             endif
 
             !if (u_lim /= u_grid%n_total_points) print *,'u_adj',u_grid%n_total_points,u_lim
             !todo rescale the u-interval

             !Integral over u:
             do u=1,u_grid%n_total_points
                if (arg_u(u) .le. tol) then !Evaluate limit for S*r -> 0
                   bessel_fac = 1.0_cfp !only lpp=0 Bessel function is non-zero
                   exp_part(0) = bessel_fac*exp_fac(u)
                   do Lboys=0,max_boys_l
                      Boys_aux(1,i,Lboys,j) = Boys_aux(1,i,Lboys,j) + exp_part(0)*u_fac(Lboys,u)*cnla(Lboys)
                   enddo
                else
                   call cfp_besi(arg_u(u), half, kode, besi_dim, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
                   val = sqrt(pi/(2.0_cfp*arg_u(u)))
                   base = (u-1)*n_lm_aux
                   do lpp=0,max_l_aux
                      bessel_fac = y(lpp+1)*val
                      exp_part(lpp) = bessel_fac*exp_fac(u)
                   enddo
                   do Lboys=0,max_boys_l
                      val = u_fac(Lboys,u)*cnla(Lboys)
                      do lpp=0,max_l_aux
                         lppmpp = lpp*lpp+lpp+1
                         ind = base + lppmpp
                         fac = exp_part(lpp)*val
                         do mpp=-lpp,lpp
                            Boys_aux(lppmpp+mpp,i,Lboys,j) = Boys_aux(lppmpp+mpp,i,Lboys,j) + fac*Xlm_product_center(ind+mpp)
                         enddo
                      enddo
                   enddo !Lboys
                endif
             enddo !u

             if (u_int_grid == 2) then ! redo the whole u-integral, requires recalculating all u-dependent values

                do u=1,u_grid_loc%n_total_points
                   do L=0,max_boys_l
                      u_fac_loc(L,u) = u_grid_loc%r_points(u)**(2*L)*u_grid_loc%weights(u)
                   enddo
                enddo

                max_u_fac_loc = maxval(u_fac_loc)*maxval(abs(cnla))

                contr_loc = cgto_shell%norm*cgto_shell%norms(j)*cgto_shell%contractions(j)

                !Calculate the product GTOs for each u: exp(-u^2*zeta(r-zeta_center)**2)*exp(-alp_{j}*(r-A)**2)
                k = 0
                u_lim_loc = u_grid_loc%n_total_points
                do u=1,u_grid_loc%n_total_points
                   u_sq_loc = u_grid_loc%r_points(u)**2
                   u_sq_zeta_loc = u_sq_loc*zeta
   
                   prod_exponent_loc(u) = u_sq_zeta_loc+cgto_shell%exponents(j)
                   prod_center_loc(1:3) = (u_sq_zeta_loc*zeta_center(1:3) + &
                                          cgto_shell%exponents(j)*cgto_shell%center(1:3))/prod_exponent_loc(u)
                   preexp_factor_loc = exp(-u_sq_zeta_loc*cgto_shell%exponents(j)/prod_exponent_loc(u)*Rab_sq)
   
                   R_distprod_loc(u) = sqrt(dot_product(prod_center_loc(1:3),prod_center_loc(1:3)))
   
                   !The preexponential factor coming from the product of the Boys GTO and the CGTO and the normalization factors
                   contr_fac_loc(u) = contr_loc*preexp_factor_loc

                   !Real spherical harmonics for the product center: result in the array Xlm_CGTO_center
                   call real_harmonics%precalculate_Xlm_for_CGTO_center(prod_center_loc(1:3),max_l_aux,Xlm_CGTO_center_loc)
                   Xlm_product_center_loc(k+1:k+n_lm_aux) = Xlm_CGTO_center_loc(1:n_lm_aux)
                   k = k + n_lm_aux

                   !Find the u for which the prefactor is too small.
                   if (u_sq_loc**(max_boys_l)*preexp_factor_loc < prec) then
                      if (u_lim_loc .eq. u_grid_loc%n_total_points) u_lim_loc = u
                   endif
                enddo

                !Precalculate the arguments needed for the integration
                do u=1,u_grid_loc%n_total_points
                   arg_u_loc(u) = 2.0_cfp*r(i)*R_distprod_loc(u)*prod_exponent_loc(u)
                   exp_fac_loc(u) = contr_fac_loc(u)*exp(-prod_exponent_loc(u)*(r(i)-R_distprod_loc(u))**2)
                enddo !u

                !Update the value of u_lim based on the value of exp_fac:
                !todo is this OK?
                u_lim_loc = u_lim_saved
                do u=u_lim_saved,u_grid_loc%n_total_points
                   if (abs(exp_fac_loc(u)) > prec) u_lim_loc=u
                enddo !u

                !skip those r-points for which the u-integral would definitely be too small.
                if (abs(max_u_fac_loc)*maxval(abs(exp_fac_loc)) < prec) then
                   !print *,'skipping',r(i)
                   cycle
                endif

 
                Boys_aux(:,i,:,j) = 0.0_cfp
                !Integral over u:
                do u=1,u_grid_loc%n_total_points
                   if (arg_u_loc(u) .le. tol) then !Evaluate limit for S*r -> 0
                      bessel_fac = 1.0_cfp !only lpp=0 Bessel function is non-zero
                      exp_part(0) = bessel_fac*exp_fac_loc(u)
                      do Lboys=0,max_boys_l
                         Boys_aux(1,i,Lboys,j) = Boys_aux(1,i,Lboys,j) + exp_part(0)*u_fac_loc(Lboys,u)*cnla(Lboys)
                      enddo
                   else
                      call cfp_besi(arg_u_loc(u), half, kode, besi_dim, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
                      val = sqrt(pi/(2.0_cfp*arg_u_loc(u)))
                      base = (u-1)*n_lm_aux
                      do lpp=0,max_l_aux
                         bessel_fac = y(lpp+1)*val
                         exp_part(lpp) = bessel_fac*exp_fac_loc(u)
                      enddo
                      do Lboys=0,max_boys_l
                         val = u_fac_loc(Lboys,u)*cnla(Lboys)
                         do lpp=0,max_l_aux
                            lppmpp = lpp*lpp+lpp+1
                            ind = base + lppmpp
                            fac = exp_part(lpp)*val
                            do mpp=-lpp,lpp
                               Boys_aux(lppmpp+mpp,i,Lboys,j) = Boys_aux(lppmpp+mpp,i,Lboys,j) + fac*Xlm_product_center_loc(ind+mpp)
                            enddo
                         enddo
                      enddo !Lboys
                   endif
                enddo !u
             endif !u_int_grid == 2
          enddo !i
       enddo !j

       !Perform the CGTO contraction: the result is stored in Boys_aux(:,:,:,1)
       do j=2,cgto_shell%number_of_primitives
          do Lboys=0,max_boys_l
             do i=1,n
                Boys_aux(1:n_lm_aux,i,Lboys,1) = Boys_aux(1:n_lm_aux,i,Lboys,1) + Boys_aux(1:n_lm_aux,i,Lboys,j)
             enddo
          enddo !Lboys
       enddo !j

       !Finally, perform translations of the solid harmonics where needed:

       if (max_boys_l > 0) then !We need to translate the solid harmonic S_{L,M} coupled to the Boys function

          if (cgto_shell%l .eq. 0) then
             !No need to translate the CGTO

             i = 2*max_boys_l+1
             !todo allocate only once and use the save attribute
             allocate(d_cfs(n,1:max_boys_l+1),R_aux(n,n_lm_boys*(max_boys_l+1)),&
                      transl_cfs_all(i,0:max_boys_l,i,max_boys_l+1),stat=err)
             if (err /= 0) then
                call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_G_X_lm','Memory allocation 3 failed.',err,1)
             endif

             R_dist = sqrt(dot_product(zeta_center,zeta_center))

             !Real spherical harmonics for the CGTO center: result in the array Xlm_CGTO_center
             call real_harmonics%precalculate_Xlm_for_CGTO_center(zeta_center,max_boys_l,Xlm_CGTO_center)
             
             !For every Lboys precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
             transl_cfs_all = 0.0_cfp
             do Lboys=0,max_boys_l

                !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
                call this%precalculate_solh_translation_coeffs(Lboys,R_dist,Xlm_CGTO_center,transl_cfs)

                transl_cfs_all(1:2*Lboys+1,0:Lboys,1:2*Lboys+1,Lboys+1) = transl_cfs(1:2*Lboys+1,0:Lboys,1:2*Lboys+1)
             enddo !Lboys

             do l=0,max_projection_l
                do m=-l,l
                   lm = l*l+l+m+1
 
                   !R_aux is effectively a 3D array (i,lpmp,Lboys) with the last
                   !two dimensions contracted into one.
                   R_aux = 0.0_cfp
                   do lp=0,max_boys_l
                      do mp=-lp,lp
                         lpmp = lp*lp+lp+mp+1
   
                         !Couple to the spherical harmonic X_{l,m}
                         do lpp=abs(lp-l),lp+l
                            lppmpp = lpp*lpp+lpp+1
                            do mpp=-lpp,lpp
                               cf = cpl%rgaunt(l,lpp,lp,m,mpp,mp)
                               if (cf /= 0.0_cfp) then
                                  do Lboys=0,max_boys_l
                                     base = n_lm_boys*Lboys
                                     R_aux(1:n,base+lpmp) = R_aux(1:n,base+lpmp) + cf*Boys_aux(lppmpp+mpp,1:n,Lboys,1)
                                  enddo
                               endif
                            enddo !mpp
                         enddo !lpp
  
                      enddo !mp
                   enddo !lp

                   do Lboys=0,max_boys_l
                      base = n_lm_boys*Lboys
                      do Mboys=1,2*Lboys+1
                         LbMb = Lboys*Lboys+Mboys
   
                         !Translate the solid harmonic S_{Lboys,M} 
                         d_cfs = 0.0_cfp
                         do lp=0,Lboys
                            lpmp = base+lp*lp+lp+1
   
                            do mp=-lp,lp
                               if (transl_cfs_all(mp+lp+1,lp,Mboys,Lboys+1) /= 0.0_cfp) then
                                  d_cfs(1:n,lp+1) = d_cfs(1:n,lp+1) + transl_cfs_all(mp+lp+1,lp,Mboys,Lboys+1)*R_aux(1:n,lpmp+mp)
                               endif
                            enddo !mp
                         enddo !lp
   
                         !Evaluate the polynomial for each r(i)
                         call cfp_eval_poly_horner_many(Lboys,r,n,d_cfs,Boys_G_Xlm(1:n,1,lm,LbMb))
   
                      enddo !Mboys
                   enddo !Lboys

                enddo !m
             enddo !l

          else
             !Translation of the CGTO needed too

             sum_l = max_boys_l+cgto_shell%l
             n_lm_sum = (sum_l+1)**2
             i = 2*sum_l+1
             ind = (2*max_boys_l+1)*(2*cgto_shell%l+1)
             !todo allocate only once and use the save attribute
             allocate(d_cfs(n,1:sum_l+1),R_aux(n,n_lm_sum*(max_boys_l+1)),&
                      transl_cfs_all_AB(i,0:sum_l,0:sum_l,ind,max_boys_l+1),stat=err)
             if (err /= 0) then
                call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_G_X_lm','Memory allocation 3 failed.',err,1)
             endif

             R_dist = sqrt(dot_product(zeta_center,zeta_center))
             RC = sqrt(dot_product(cgto_shell%center,cgto_shell%center)) 

             !Real spherical harmonics for the CGTO center: result in the module array Xlm_CGTO_C_center
             call real_harmonics%precalculate_Xlm_for_CGTO_center(cgto_shell%center,cgto_shell%l,Xlm_CGTO_C_center)

             !Real spherical harmonics for the Boys center: result in the module array Xlm_CGTO_center
             call real_harmonics%precalculate_Xlm_for_CGTO_center(zeta_center,max_boys_l,Xlm_CGTO_center)

             !For every Lboys precalculate the coefficients in the translation formula for the pair of solid harmonics
             transl_cfs_all_AB = 0.0_cfp
             do Lboys=0,max_boys_l
                call this%precalculate_pair_solh_translation_coeffs(Lboys,R_dist,Xlm_CGTO_center,cgto_shell%l,&
                                                                    RC,Xlm_CGTO_C_center,transl_cfs_AB)
                ind = 0
                do Mboys=1,2*Lboys+1
                   do Mg=1,cgto_shell%number_of_functions
                      ind = ind + 1
                      do la_p_lb=0,Lboys+cgto_shell%l
                         do lp=0,Lboys+cgto_shell%l
                            transl_cfs_all_AB(1:2*lp+1,lp,la_p_lb,ind,Lboys+1) = transl_cfs_AB(1:2*lp+1,lp,la_p_lb,Mg,Mboys)
                         enddo !lp
                      enddo !la_p_lb
                   enddo !Mg
                enddo !Mboys
             enddo !Lboys

             do l=0,max_projection_l
                do m=-l,l
                   lm = l*l+l+m+1
 
                   !R_aux is effectively a 3D array (i,lpmp,Lboys) with the last
                   !two dimensions contracted into one.
                   R_aux = 0.0_cfp
                   do lp=0,sum_l
                      do mp=-lp,lp
                         lpmp = lp*lp+lp+mp+1
   
                         !Couple to the spherical harmonic X_{l,m}
                         do lpp=abs(lp-l),lp+l
                            lppmpp = lpp*lpp+lpp+1
                            do mpp=-lpp,lpp
                               cf = cpl%rgaunt(l,lpp,lp,m,mpp,mp)
                               if (cf /= 0.0_cfp) then
                                  do Lboys=0,max_boys_l
                                     base = n_lm_sum*Lboys
                                     R_aux(1:n,base+lpmp) = R_aux(1:n,base+lpmp) + cf*Boys_aux(lppmpp+mpp,1:n,Lboys,1)
                                  enddo
                               endif
                            enddo !mpp
                         enddo !lpp
  
                      enddo !mp
                   enddo !lp

                   do Lboys=0,max_boys_l
                      base = n_lm_sum*Lboys
                      ind = 0
                      do Mboys=1,2*Lboys+1
                         LbMb = Lboys*Lboys+Mboys
                         do Mg=1,cgto_shell%number_of_functions
                            ind = ind + 1
   
                            !Translate the solid harmonics S_{Lboys,M}*S_{Lg,Mg}
                            d_cfs = 0.0_cfp
                            do la_p_lb=0,sum_l
                               do lp=0,sum_l 
                                  lpmp = base+lp*lp+lp+1
                                  do mp=-lp,lp
                                     if (transl_cfs_all_AB(mp+lp+1,lp,la_p_lb,ind,Lboys+1) /= 0.0_cfp) then
                                        d_cfs(1:n,la_p_lb+1) = d_cfs(1:n,la_p_lb+1)&
                                                             + transl_cfs_all_AB(mp+lp+1,lp,la_p_lb,ind,Lboys+1)&
                                                             * R_aux(1:n,lpmp+mp)
                                     endif
                                  enddo !mp
                               enddo !lp
                            enddo !la_p_lb
      
                            !Evaluate the polynomial for each r(i)
                            call cfp_eval_poly_horner_many(sum_l,r,n,d_cfs,Boys_G_Xlm(1:n,Mg,lm,LbMb))

                         enddo !Mg
                      enddo !Mboys
                   enddo !Lboys

                enddo !m
             enddo !l

          endif

       else !max_l_boys == 0

          if (cgto_shell%l .eq. 0) then
             !No need to translate the CGTO

             do lm=1,n_lm
                Boys_G_Xlm(1:n,1,lm,1) = Boys_aux(lm,1:n,0,1)
             enddo

          else
             !Translation of the CGTO needed.

             !todo allocate only once and use the save attribute
             allocate(d_cfs(n,1:cgto_shell%l+1),R_aux(n,(cgto_shell%l+1)**2),stat=err)
             if (err /= 0) then
                call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_G_X_lm','Memory allocation 3 failed.',err,1)
             endif

             R_dist = sqrt(dot_product(cgto_shell%center,cgto_shell%center))

             !Real spherical harmonics for the CGTO center: result in the array Xlm_CGTO_center
             call real_harmonics%precalculate_Xlm_for_CGTO_center(cgto_shell%center,cgto_shell%l,Xlm_CGTO_center)

             !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
             call this%precalculate_solh_translation_coeffs(cgto_shell%l,R_dist,Xlm_CGTO_center,transl_cfs)

             do l=0,max_projection_l
                do m=-l,l
                   lm = l*l+l+m+1

                   R_aux = 0.0_cfp
                   do lp=0,cgto_shell%l
                      do mp=-lp,lp
                         lpmp = lp*lp+lp+mp+1

                         !Couple to the spherical harmonic X_{l,m}
                         do lpp=abs(lp-l),lp+l
                            lppmpp = lpp*lpp+lpp+1
                            do mpp=-lpp,lpp
                               cf = cpl%rgaunt(l,lpp,lp,m,mpp,mp)
                               if (cf /= 0.0_cfp) then
                                  R_aux(1:n,lpmp) = R_aux(1:n,lpmp) + cf*Boys_aux(lppmpp+mpp,1:n,0,1)
                               endif
                            enddo !mpp
                         enddo !lpp

                      enddo !mp
                   enddo !lp

                   do Mg=1,cgto_shell%number_of_functions

                      !Translate the solid harmonic S_{cgto_shell%l,M} of the CGTO
                      d_cfs = 0.0_cfp
                      do lp=0,cgto_shell%l
                         lpmp = lp*lp+lp+1

                         do mp=-lp,lp
                            if (transl_cfs(mp+lp+1,lp,Mg) /= 0.0_cfp) then
                               d_cfs(1:n,lp+1) = d_cfs(1:n,lp+1) + transl_cfs(mp+lp+1,lp,Mg)*R_aux(1:n,lpmp+mp)
                            endif
                         enddo !mp
                      enddo !lp

                      !Evaluate the polynomial for each r(i)
                      call cfp_eval_poly_horner_many(cgto_shell%l,r,n,d_cfs,Boys_G_Xlm(1:n,Mg,lm,1))

                   enddo !Mg
 
                enddo !m
             enddo !l

          endif

       endif

  end subroutine Boys_projected_on_G_X_lm

  !> Calculates the angular integrals \int d\Omega(r_{1}) X_{l,m}(r_{1}) norm*C_{-1,L}^{a}F_{L}(zeta*(r_{1}-zeta_center)**2)*S_{L,M}(r_{1}-zeta_center)*sqrt((2*L+1)/(4*pi))
  !> for the grid of r1 points and l,m projections up to max_projection_l and L,M up to max_boys_l.
  !> F_{L}(x) = \int_{0}^{1} du u**2*exp(-u**2*x)
  subroutine Boys_projected_on_X_lm(this,u_grid,norm,zeta,zeta_center,r,max_projection_l,max_boys_l,Boys_Xlm)
     use phys_const_gbl, only: fourpi, pi
     use const_gbl, only: exp_arg_asymptotic_wp, exp_arg_asymptotic_ep
     use special_functions_gbl, only: cfp_besi, cfp_eval_poly_horner_many
     use general_quadrature_gbl, only: quadrature_u_integral
     implicit none
     class(GG_shell_integrals_obj) :: this
     type(pw_expansion_obj) :: u_grid
     real(kind=cfp), intent(in) :: norm, zeta, zeta_center(3)
     real(kind=cfp), allocatable :: r(:)
     integer, intent(in) :: max_projection_l,max_boys_l
     !OUT:
     real(kind=cfp), allocatable :: Boys_Xlm(:,:,:) !Boys_Xlm(radial point,lm,LbMb)

     integer, parameter :: kode = 2, n_u_points = 18
     real(kind=cfp), parameter :: half = 0.5_cfp

     integer :: err, n, n_lm, LbMb, lm, i, u, Lboys, Mboys, max_l_aux, besi_dim, l, lp, mp, lpp, mpp, j, nz, m, n_lm_boys,&
                m_ind, lap, map, lambda, CGTO_L, l_min,l_max
     real(kind=cfp) :: alp, cnla(0:max_boys_l+1), R_zeta, R_dist, bessel_fac, arg, exp_arg, arg_u, exp_arg_u, exp_fac, exp_part,&
                       u_sq, tol, val, cf
     real(kind=cfp) :: exp_arg_asymptotic, x(n_u_points), w(n_u_points), scale_factor
     real(kind=cfp), allocatable :: y(:), u_fac(:,:), d_cfs(:,:), Boys_aux(:,:,:), Xlm_CGTO_center(:), transl_cfs(:,:,:),&
                                    u_fac_asym(:,:), u_points(:,:), c_lambda(:,:,:)
     logical :: non_zero, is_zero

       n = size(r)
       n_lm = (max_projection_l+1)**2
       n_lm_boys = (max_boys_l+1)**2

       tol = F1MACH(4,cfp_dummy)

       if (allocated(Boys_Xlm)) then
          if (size(Boys_Xlm,1) /= n .or. size(Boys_Xlm,3) /= n_lm_boys .or. size(Boys_Xlm,2) /= n_lm) deallocate(Boys_Xlm)
       endif

       if (.not. allocated(Boys_Xlm)) then
          allocate(Boys_Xlm(n,n_lm,n_lm_boys),stat=err)
          if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_X_lm','Memory allocation failed.',1,1)
       endif

       Boys_Xlm(:,:,:) = 0.0_cfp

       R_dist = sqrt(dot_product(zeta_center,zeta_center))
       R_zeta = zeta*R_dist

       alp = 1.0_cfp/(4.0_cfp*zeta)
       do Lboys=0,max_boys_l
          val = sqrt((2*Lboys+1.0_cfp)/fourpi) !get rid of the leading factor sqrt(4pi/(2*l+1)) for the solid harmonics since the NAI formulae are only for r**l*Xlm
          cnla(Lboys) = val*norm*fourpi/(2.0_cfp*alp)**(Lboys+0.5_cfp)*(-1)**Lboys
       enddo

       max_l_aux = max_projection_l + max_boys_l

       besi_dim = max_l_aux+1
       i = max(max_l_aux,1)
       j = max(max_boys_l,1)
       allocate(y(besi_dim),u_fac(0:j,u_grid%n_total_points),Boys_aux(0:i,n,0:j),u_fac_asym(n_u_points,0:j),&
                u_points(n_u_points,0:j),stat=err)
       if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_X_lm','Memory allocation 2 failed.',err,1)

       do u=1,u_grid%n_total_points
          do L=0,max_boys_l
             u_fac(L,u) = u_grid%r_points(u)**(2*L)*u_grid%weights(u)
          enddo
       enddo

       do L=0,max_boys_l
          call quadrature_u_integral(x,w,n_u_points,L)
          do u=1,n_u_points
             u_fac_asym(u,L) = w(u)
             u_points(u,L) = x(u)
          enddo
       enddo

       !Real spherical harmonics for the GTO center: result in the array Xlm_CGTO_center
       call real_harmonics%precalculate_Xlm_for_CGTO_center(zeta_center,max_l_aux,Xlm_CGTO_center)

       if (cfp .eq. wp) then
          exp_arg_asymptotic = exp_arg_asymptotic_wp
       elseif (cfp .eq. ep) then
          exp_arg_asymptotic = exp_arg_asymptotic_ep
       else
          stop "error: unsupported value of cfp"
       endif

       Boys_aux = 0.0_cfp
       do i=1,n

          arg = 2.0_cfp*r(i)*R_zeta
          exp_arg = zeta*(r(i)-R_dist)**2

          !todo if R_dist == 0 then the u-integral is just the standard Boys function and lpp_max = 0.
          if (exp_arg .ge. exp_arg_asymptotic) then
             !The u-integral can be calculated asymptotically using a qudrature on the interval [0;infinity] for the weight function u**(2*L)*exp(-u**2).
             exp_fac = 1.0_cfp/exp_arg
             do Lboys=0,max_boys_l
                scale_factor = exp_arg**(-(Lboys+0.5_cfp))
                do u=n_u_points,1,-1 !1,n_u_points
                   arg_u = arg*u_points(u,Lboys)**2*exp_fac
                   if (arg_u .le. tol) then !Evaluate limit for S*r -> 0
                      bessel_fac = 1.0_cfp !only lpp=0 Bessel function is non-zero
                      Boys_aux(0,i,Lboys) = Boys_aux(0,i,Lboys) + scale_factor*bessel_fac*u_fac_asym(u,Lboys)*cnla(Lboys)
                   else
                      call cfp_besi(arg_u, half, kode, besi_dim, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
                      val = sqrt(pi/(2.0_cfp*arg_u))
                      do lpp=0,max_l_aux
                         bessel_fac = y(lpp+1)*val
                         Boys_aux(lpp,i,Lboys) = Boys_aux(lpp,i,Lboys) + scale_factor*bessel_fac*u_fac_asym(u,Lboys)*cnla(Lboys)
                      enddo !lpp
                   endif
                enddo !u
             enddo !Lboys
          else
             !The u-integral must be calculated using Gauss-Legendre quadratures on the interval [0;1]
             do u=1,u_grid%n_total_points
                u_sq = u_grid%r_points(u)**2
   
                arg_u = arg*u_sq
                exp_arg_u = exp_arg*u_sq
   
                exp_fac = exp(-exp_arg_u)
   
                if (arg_u .le. tol) then !Evaluate limit for S*r -> 0
                   bessel_fac = 1.0_cfp !only lpp=0 Bessel function is non-zero
                   exp_part = bessel_fac*exp_fac
                   do Lboys=0,max_boys_l
                      Boys_aux(0,i,Lboys) = Boys_aux(0,i,Lboys) + exp_part*u_fac(Lboys,u)*cnla(Lboys)
                   enddo
                else
                   call cfp_besi(arg_u, half, kode, besi_dim, y, nz) !cfp_besi gives: y_{alpha+k-1}, k=1,...,N. Hence N=data%l+1 is needed to get y_{data%l}.
                   val = sqrt(pi/(2.0_cfp*arg_u))
                   do lpp=0,max_l_aux
                      bessel_fac = y(lpp+1)*val
                      exp_part = bessel_fac*exp_fac
                      do Lboys=0,max_boys_l
                         Boys_aux(lpp,i,Lboys) = Boys_aux(lpp,i,Lboys) + exp_part*u_fac(Lboys,u)*cnla(Lboys)
                      enddo
                   enddo !l
                endif
             enddo !u
          endif

       enddo !i

       if (max_boys_l > 0) then !Adapted from omp_calculate_CGTO_pw_coefficients_analytic

          allocate(d_cfs(n,1:max_boys_l+1),stat=err)
          if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','Boys_projected_on_X_lm','Memory allocation 3 failed.',err,1)

          do Lboys=0,max_boys_l

             !Precalculate the coefficients in the translation formula for the solid harmonics: this requires Xlm_CGTO_center
             call this%precalculate_solh_translation_coeffs(Lboys,R_dist,Xlm_CGTO_center,transl_cfs)

             do l=0,max_projection_l
                do m=-l,l
                   lm = l*l+l+m+1

                   !Calculate the coupling coefficients needed to evaluate the projection using SCE of GTO.
                   !call this%calculate_lambda_couplings(Lboys,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
                   call calculate_lambda_couplings(Lboys,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
     
                   do Mboys=-Lboys,Lboys
                      m_ind = Mboys+Lboys+1
                      LbMb = Lboys*Lboys+Lboys+Mboys+1


                      !Skip the symmetry-forbidden projections
                      is_zero = .true.
                      do lambda=0,l+Lboys
                         do lp=0,Lboys
                            if (c_lambda(lambda,lp,m_ind) /= 0.0_cfp) then
                               is_zero = .false.
                               exit
                            endif
                         enddo !lp
                      enddo !lambda
                      if (is_zero) cycle

                       !Evaluate the polynomial coefficients
                       do lp=0,Lboys
                          l_min = abs(l-lp)
                          l_max = l+lp
                          j = lp+1
                          do i=1,n
                             d_cfs(i,j) = sum(c_lambda(l_min:l_max,lp,m_ind)*Boys_aux(l_min:l_max,i,Lboys))
                          enddo
                       enddo !lp
 
                       !Evaluate the polynomial for each r(i)
                       call cfp_eval_poly_horner_many(Lboys,r,n,d_cfs,Boys_Xlm(1:n,lm,LbMb))
     
                   enddo !Mboys
     
                enddo !m
             enddo !l

          enddo !Lboys

       else

          do l=0,max_projection_l
             do m=-l,l
                lm = l*l+l+m+1
                if (Xlm_CGTO_center(lm) /= 0.0_cfp) then
                   Boys_Xlm(1:n,lm,1) = Boys_aux(l,1:n,0)*Xlm_CGTO_center(lm)
                endif
             enddo !m
          enddo !l

       endif

  end subroutine Boys_projected_on_X_lm

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
         if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','radial_grid_CGTO_pair','Memory allocation 1 failed.',err,1)

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
         if (err /= 0) call xermsg('GG_shell_mixed_integrals_mod','radial_grid_CGTO_pair','Memory allocation 2 failed.',err,1)

         !Construct quadratures for the individual ranges:
         cnt = 0
         do i=1,n_ranges
           !Prepare the quadrature grid for the interval [range_start,range_end] expanding the canonical G-L quadrature.
           call gl_expand_A_B(x,w,n_points,tmp_r,tmp_w,range_start_end(1,i),range_start_end(2,i))
           r_points(cnt+1:cnt+n_points) = tmp_r(1:n_points)
           weights(cnt+1:cnt+n_points) = tmp_w(1:n_points)
           cnt = cnt + n_points
           !print *,'gto',range_start,range_end
         enddo
!
         if (cnt /= n_total_points) then
            call xermsg('GG_shell_mixed_integrals_mod','radial_grid_CGTO_pair','Error constructing the radial grid.',3,1)
         endif

  end subroutine radial_grid_CGTO_pair

  !> \warning Requires precalculated values of the real spherical harmonics at the position of the CGTO nucleus. The coupling coefficients should also be precalculated for performance reasons.
  subroutine precalculate_pair_solh_translation_coeffs(this,CGTO_A_L,RA_A,Xlm_CGTO_A_center,&
                                                            CGTO_B_L,RA_B,Xlm_CGTO_B_center,transl_cfs_AB)
     use phys_const_gbl, only: fourpi
     implicit none
     class(GG_shell_integrals_obj) :: this
     real(kind=cfp), allocatable :: Xlm_CGTO_A_center(:), Xlm_CGTO_B_center(:), transl_cfs_AB(:,:,:,:,:)
     integer, intent(in) :: CGTO_A_L, CGTO_B_L
     real(kind=cfp), intent(in) :: RA_A, RA_B

     integer :: n_mp, err, CGTO_M, lp, mp, max_lp, lp_min, lp_max, CGTO_A_M, CGTO_A_M_ind, CGTO_B_M, CGTO_B_M_ind, &
                la_p, lb_p, la_p_lb_p, ma_p, mb_p
     real(kind=cfp), allocatable :: transl_cfs_A(:,:,:), transl_cfs_B(:,:,:)
   
        !Translation coefficients for the individual CGTOs
        call this%precalculate_solh_translation_coeffs(CGTO_A_L,RA_A,Xlm_CGTO_A_center,transl_cfs_A)
        call this%precalculate_solh_translation_coeffs(CGTO_B_L,RA_B,Xlm_CGTO_B_center,transl_cfs_B)

        max_lp = CGTO_A_L+CGTO_B_L
        n_mp = 2*max_lp+1

        if (allocated(transl_cfs_AB)) deallocate(transl_cfs_AB)
        allocate(transl_cfs_AB(n_mp,0:max(max_lp,1),0:max(max_lp,1),2*CGTO_B_L+1,2*CGTO_A_L+1),stat=err)
        if (err /= 0) then
           call xermsg('cgto_pw_expansions_mod','precalculate_pair_solh_translation_coeffs','Memory allocation failed.',err,1)
        endif

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
                             if (mod(lp+la_p+lb_p,2) /= 0) cycle !selection rule for Gaunt coefficients
                             do mp=-lp,lp
                                transl_cfs_AB(mp+lp+1,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) = &
                                transl_cfs_AB(mp+lp+1,lp,la_p_lb_p,CGTO_B_M_ind,CGTO_A_M_ind) &
                                + transl_cfs_A(ma_p+la_p+1,la_p,CGTO_A_M_ind)&
                                * transl_cfs_B(mb_p+lb_p+1,lb_p,CGTO_B_M_ind)&
                                * cpl%rgaunt(lp,la_p,lb_p,mp,ma_p,mb_p)
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
  subroutine precalculate_solh_translation_coeffs(this,CGTO_L,RA,Xlm_CGTO_center,transl_cfs)
     use phys_const_gbl, only: fourpi
     implicit none
     class(GG_shell_integrals_obj) :: this
     real(kind=cfp), allocatable :: transl_cfs(:,:,:), Xlm_CGTO_center(:)
     integer, intent(in) :: CGTO_L
     real(kind=cfp), intent(in) :: RA

     integer :: n_mp, err, CGTO_M, lp, lpp, mp, mpp, ind
     real(kind=cfp) :: fac, sum_mpp
   
        n_mp = 2*CGTO_L+1

        if (allocated(transl_cfs)) deallocate(transl_cfs)
        allocate(transl_cfs(n_mp,0:max(CGTO_L,1),n_mp),stat=err)
        if (err /= 0) then
           call xermsg('GG_shell_mixed_integrals_mod','precalculate_solh_translation_coeffs','Memory allocation failed.',err,1)
        endif

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
                    sum_mpp = sum_mpp + cpl%G_real_cf(CGTO_L,lp,CGTO_M,mp,mpp)*Xlm_CGTO_center(ind)
                 enddo !mpp
                 transl_cfs(mp+lp+1,lp,CGTO_M+CGTO_L+1) = sum_mpp*fac*(-1)**(lpp)*RA**lpp
              enddo !mp
           enddo !lp
        enddo !CGTO_M

  end subroutine precalculate_solh_translation_coeffs

  !> \warning Requires precalculated values of the real spherical harmonics at the position of the CGTO nucleus. The coupling coefficients should also be precalculated for performance reasons.
  !todo This routine is exactly the same as a routine with same name in cgto_pw_expansions_gbl, except here is one additional input
  !parameter - this - GG_shell_integrals_obj, but it is not used
!  subroutine calculate_lambda_couplings(this,CGTO_L,l,m,Xlm_CGTO_center,transl_cfs,c_lambda)
!     implicit none
!     class(GG_shell_integrals_obj) :: this
!     real(kind=cfp), allocatable :: transl_cfs(:,:,:), Xlm_CGTO_center(:)
!     integer, intent(in) :: CGTO_L, l, m
!     !OUTPUT:
!     real(kind=cfp), allocatable :: c_lambda(:,:,:)
!
!     real(kind=cfp) :: cf, sum_mu
!     integer :: CGTO_M, CGTO_M_ind, lm, lp, lpp, mp, mu, lambda, ind, lambda_max, err, d1, d2, d3
!
!        lambda_max = l+CGTO_L
!
!        d1 = max(lambda_max,1)
!        d2 = max(CGTO_L,1)
!        d3 = 2*CGTO_L+1
!        if (allocated(c_lambda)) then
!           if (ubound(c_lambda,1) < d1 .or. ubound(c_lambda,2) < d2 .or. ubound(c_lambda,3) < d3) then
!              deallocate (c_lambda)
!           end if
!        end if
!        if (.not. allocated(c_lambda)) then
!           allocate(c_lambda(0:d1,0:d2,d3),stat=err)
!           if (err .ne. 0) call xermsg('GG_shell_mixed_integrals_mod','calculate_lambda_couplings','Memory allocation failed.',&
!                                       err,1)
!        endif
!
!        do CGTO_M=-CGTO_L,CGTO_L
!           CGTO_M_ind = CGTO_L+CGTO_M+1
!           lm = l*l+l+m+1
!           do lp=0,CGTO_L
!              lpp = CGTO_L-lp
!
!              c_lambda(0:lambda_max,lp,CGTO_M_ind) = 0.0_cfp
!              do lambda=abs(l-lp),l+lp
!                 if (mod(l+lp+lambda,2) /= 0) cycle !selection rule for the Gaunt coefficients
!                 do mp=-lp,lp
!                    cf = transl_cfs(mp+lp+1,lp,CGTO_M_ind)
!                    if (cf .eq. 0.0_cfp) cycle
!                    sum_mu = 0.0_cfp
!                    do mu=-lambda,lambda
!                       ind = lambda*lambda+lambda+mu+1
!                       sum_mu = sum_mu + cf*cpl%rgaunt(l,lp,lambda,m,mp,mu)*Xlm_CGTO_center(ind)
!                    enddo !mu
!                    c_lambda(lambda,lp,CGTO_M_ind) = c_lambda(lambda,lp,CGTO_M_ind) + sum_mu
!                 enddo !mp
!              enddo !lambda
!
!           enddo !lp
!        enddo !CGTO_M
!     
!  end subroutine calculate_lambda_couplings

  pure function exp_fac_exponent(u,zeta,cgto_exp,zeta_center,cgto_center,r1)
     implicit none
     real(kind=cfp), intent(in) :: u, zeta, cgto_exp, zeta_center(3), cgto_center(3), r1
     real(kind=cfp) :: exp_fac_exponent
     real(kind=cfp) :: u_sq_zeta, u_sq_zeta_offset, norm_1_sq, norm_2, vec(3)

        u_sq_zeta = u*u*zeta
        u_sq_zeta_offset = u_sq_zeta + cgto_exp
        vec = zeta_center - cgto_center
        norm_1_sq = dot_product(vec,vec)
        vec = u_sq_zeta*zeta_center + cgto_exp*cgto_center
        norm_2 = sqrt(dot_product(vec,vec))

        exp_fac_exponent = - (u_sq_zeta*cgto_exp*norm_1_sq + (r1*u_sq_zeta_offset - norm_2)**2)/u_sq_zeta_offset

  end function exp_fac_exponent

  pure function arg_min(zeta,cgto_exp,zeta_center,cgto_center)
     implicit none
     real(kind=cfp), intent(in) :: zeta, cgto_exp, zeta_center(3), cgto_center(3)
     real(kind=cfp) :: arg_min
     real(kind=cfp) :: prod, zeta_norm_sq

        arg_min = -1.0_cfp

        prod = dot_product(zeta_center, cgto_center)
        zeta_norm_sq = zeta*dot_product(zeta_center, zeta_center)
        if (prod < 0) then
           arg_min = sqrt(- prod * cgto_exp / zeta_norm_sq)
        endif

  end function arg_min

end module GG_shell_mixed_integrals_gbl
