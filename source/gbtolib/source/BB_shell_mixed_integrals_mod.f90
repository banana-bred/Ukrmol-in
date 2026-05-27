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
module BB_shell_mixed_integrals_gbl
  use basis_data_generic_gbl
  use coupling_obj_gbl
  use precisn_gbl
  use const_gbl, only: stdout
  use utils_gbl, only: xermsg
  use cgto_pw_expansions_gbl, only: CGTO_shell_pw_expansion_obj
  use special_functions_gbl, only: real_harmonics_obj 
  use grid_gbl, only: grid_r1_r2
  implicit none
  
  private

  public BB_shell_integrals_obj

  !we replaced all the grid_r1_r2 grid with grid_r1_r2, so that it does not conflict with the grid_r1_r2 of type
  !(legendre_grid_r1_r2) from cgto_pw_expansions_gbl

  !> Object that can be used to evaluate the (BB|BB) and the mixed BTO/CGTO 2-electron integrals (BG|GG).
  !> This object can only be used once the grid_r1_r2 from grid_mod has been initialized.
  type :: BB_shell_integrals_obj
     !> Data for the CGTO shell which will be used in the (BG|BB) integrals
     type(CGTO_shell_data_obj) :: cgto_shell
     !> Starting indices for the functions in cgto_shell.
     integer, private :: cgto_shell_starting_index = -1
     !> Set to .true. by init if the CGTO shell is centered on the CMS.
     logical, private :: cgto_shell_is_cms = .false.
     !> Final set of integrals of the BGBB type for all unique shells of BTOs.
     real(kind=cfp), allocatable, private :: BGBB(:,:)
     !> Needed by eval_radial_BBBB_integrals.
     real(kind=cfp), allocatable, private :: radial_BBBB(:,:,:)
     !> Maximum L for which the radial_BBBB integrals have been calculated and the maximum L in the BTO basis.
     integer, private :: max_l = -1, max_bspline_l = -1
     !> Needed for BBBB and BGBB integral reconstruction.
     integer, allocatable, private :: pairs_to_unique_pairs(:,:)
     !> Needed for BBBB and BGBB integral reconstruction.
     integer, private :: n_unique_pairs = 0, n_bsplines = 0, n_cgto_shells = 0, n_lbmb = 0, n_lm = 0, offset_BBBB = 0
     !> Needed for BGBB and BBBB class reconstruction.
     integer, allocatable, private :: BGBB_aux_offset(:), bto_starting_indices(:,:), bspline_range(:,:), &
                                      cgto_starting_indices(:), cgto_shell_n_fns(:)
     !> Set to .true. following a call to eval_radial_BBBB_integrals.
     logical, private :: radial_BBBB_calculated = .false.
     !> Set to .true. following a call to prepare_for_reconstruction.
     logical, private :: reconstruction_prepared = .false.
  contains
     procedure :: eval_radial_BBBB_integrals
     procedure :: eval_BGBB_integrals
     procedure :: eval_BGBB_auxiliary_integrals
     procedure :: get_BGBB_shell_integrals
     procedure :: get_BBBB_shell_integrals
     procedure :: get_BBBB_integrals_all_pq_fixed_BB
     procedure :: get_BGBB_integrals_all_pq_fixed_BG
     procedure :: get_BGBB_integrals_all_pq_fixed_BB
     procedure :: get_lambda_function
     procedure :: get_radial_BBBB
     procedure :: prepare_for_reconstruction
     procedure :: reconstruct_BBBB_shell_integrals
     procedure :: final
  end type BB_shell_integrals_obj

  !> Used to evaluate various coupling coefficients. PRIVATE FOR ALL THREADS.
  type(couplings_type) :: cpl

  !$OMP THREADPRIVATE(cpl)

contains

  subroutine get_lambda_function(this)
     implicit none
     class(BB_shell_integrals_obj) :: this

     integer :: max_l

        !Precalculate the lambda function. Do not forget that it must be multiplied by 1/r1 when integrating over r1!
        max_l = 2*grid_r1_r2%max_bspline_l
        if (.not.(allocated(grid_r1_r2%lambda_BB_r1_l_ij))) call grid_r1_r2%eval_lambda_l_BB(max_l)

  end subroutine get_lambda_function

  subroutine final(this)
     implicit none
     class(BB_shell_integrals_obj) :: this

        this%cgto_shell_starting_index = -1
        this%cgto_shell_is_cms = .false.
        if (allocated(this%radial_BBBB)) deallocate(this%radial_BBBB)
        if (allocated(this%BGBB_aux_offset)) deallocate(this%BGBB_aux_offset)
        if (allocated(this%pairs_to_unique_pairs)) deallocate(this%pairs_to_unique_pairs)
        if (allocated(this%bto_starting_indices)) deallocate(this%bto_starting_indices)
        if (allocated(this%bspline_range)) deallocate(this%bspline_range)
        this%radial_BBBB_calculated = .false.
        this%max_l = -1
        this%max_bspline_l = -1
        this%n_unique_pairs = 0
        this%n_bsplines = 0
        this%reconstruction_prepared = .false.

  end subroutine final

  !> Assumes that the lambda function has been precalculated before using a call to this%get_lambda function.
  subroutine eval_BGBB_auxiliary_integrals(this,cgto_shell,cgto_shell_starting_index,aux_radial_BGBB)
     use special_functions_gbl, only: ipair
     implicit none
     class(BB_shell_integrals_obj) :: this
     type(CGTO_shell_data_obj) :: cgto_shell
     integer, intent(in) :: cgto_shell_starting_index
     real(kind=cfp), allocatable :: aux_radial_BGBB(:)

     real(kind=cfp) :: RA, tol
     integer :: max_l, l, m, l1,m1,Mg,i,j,lm,l1m1,si,ei,n_lbmb,err,p,n_lm,base,n,n1,n2,n3,n4
     type(CGTO_shell_pw_expansion_obj) :: cgto_pw

        this%cgto_shell = cgto_shell
        this%cgto_shell_starting_index = cgto_shell_starting_index

        RA = sqrt(dot_product(cgto_shell%center,cgto_shell%center))

        tol = F1MACH(4,cfp_dummy)

        this%cgto_shell_is_cms = .false.

        if (RA .le. tol) this%cgto_shell_is_cms = .true.

        n_lbmb = (grid_r1_r2%max_bspline_l+1)**2
        max_l = 2*grid_r1_r2%max_bspline_l

        !Check that lambda function has been precalculated using call grid_r1_r2%eval_lambda_l_BB(max_l). Do not forget that it must be multiplied by 1/r1 when integrating over r1!
        if (.not. allocated(grid_r1_r2%lambda_BB_r1_l_ij)) then
           call xermsg('BB_shell_mixed_integrals_mod','eval_BGBB_auxiliary_integrals',&
                       'The lambda function has not been precalculated!',1,1)
        endif

        !Evaluate the double angular integrals for the CGTO on the r1 grid
        call cgto_pw%init_CGTO_shell_pw_expansion(this%cgto_shell,this%cgto_shell_starting_index)
        call cgto_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
        call cgto_pw%eval_CGTO_shell_pw_expansion(grid_r1_r2%bspline_grid%knots,grid_r1_r2%max_bspline_l,0,max_l)

        n_lm = (max_l+1)**2
        n = grid_r1_r2%n_unique_pairs*n_lm*n_lbmb*grid_r1_r2%last_bspline_inner*this%cgto_shell%number_of_functions
        if (size(aux_radial_BGBB) .ne. n) then
           if (allocated(aux_radial_BGBB)) deallocate(aux_radial_BGBB)
           allocate(aux_radial_BGBB(n),stat=err)
           if (err /= 0) then
              call xermsg('BB_shell_mixed_integrals_mod','eval_BGBB_auxiliary_integrals','Memory allocation failed.',err,1)
           endif
        endif

        !Evaluate the BGBB auxiliary integrals for all shells of BTOs in the first coordinate and all unique combinations of BTO indices in the second coordinate:
        n1 = grid_r1_r2%n_unique_pairs
        n2 = n1*n_lm
        n3 = n2*grid_r1_r2%last_bspline_inner
        n4 = n3*n_lbmb
        aux_radial_BGBB = 0.0_cfp
        do l1=0,grid_r1_r2%max_bspline_l
           do m1=-l1,l1
              l1m1 = l1*l1+l1+m1+1
              do Mg=1,this%cgto_shell%number_of_functions

                 !Evaluate the radial integrals:
                 !todo if the CGTO is on CMS then the radial integrals depend only on l and should calculated outside of the l1,m1,Mg loops.
                 do l=0,max_l
                    do m=-l,l
                       lm = l*l+l+m+1

                       p = cgto_pw%non_neg_indices_l_lp(Mg,l1m1,lm)
                       if (p > 0) then
                          do i=1,grid_r1_r2%last_bspline_inner
                             si = grid_r1_r2%bspline_start_end_r1(1,i)
                             ei = grid_r1_r2%bspline_start_end_r1(2,i)
                             base = n1*(lm-1) + n2*(i-1) + n3*(l1m1-1) + n4*(Mg-1)
                             do j=1,grid_r1_r2%n_unique_pairs
                                !int_{r1} r1**2 * 1/r1*lambda_BB_r1_l_ij * B(r1)/r1 * G_{lm,l1m1}(r1)
                                !store the integrals as a 5D array (j,lm,i,l1m1,Mg)
                                aux_radial_BGBB(base+j) = sum(grid_r1_r2%B_vals_r1(si:ei,i)&
                                                              *cgto_pw%gaunt_angular_integrals(si:ei,p)&
                                                              *grid_r1_r2%lambda_BB_r1_l_ij(si:ei,l,j))
                             enddo !j
                          enddo !i
                       endif

                    enddo !m
                 enddo !l

              enddo !Mg
           enddo !m1
        enddo !l1

  end subroutine eval_BGBB_auxiliary_integrals

  !> Assumes that the lambda function has been precalculated before using a call to this%get_lambda function.
  subroutine eval_BGBB_integrals(this,cgto_shell,cgto_shell_starting_index)
     use special_functions_gbl, only: ipair
     implicit none
     class(BB_shell_integrals_obj) :: this
     type(CGTO_shell_data_obj) :: cgto_shell
     integer, intent(in) :: cgto_shell_starting_index

     real(kind=cfp) :: RA, tol, cf
     integer :: max_l, l, m, l1,m1,l3,m3,l4,m4,Mg,i,j,lm,ind_21,ind_43,l1m1,l3m3,l4m4,si,ei,n_lbmb,ind_lm,err,p,n_lm
     integer(kind=1), allocatable :: radial_lm_is_significant(:)
     real(kind=cfp), allocatable :: radial_BGBB(:,:,:)
     type(CGTO_shell_pw_expansion_obj) :: cgto_pw

        this%cgto_shell = cgto_shell
        this%cgto_shell_starting_index = cgto_shell_starting_index

        RA = sqrt(dot_product(cgto_shell%center,cgto_shell%center))

        tol = F1MACH(4,cfp_dummy)

        this%cgto_shell_is_cms = .false.

        if (RA .le. tol) this%cgto_shell_is_cms = .true.

        n_lbmb = (grid_r1_r2%max_bspline_l+1)**2
        max_l = 2*grid_r1_r2%max_bspline_l

        !Check that lambda function has been precalculated using call grid_r1_r2%eval_lambda_l_BB(max_l). Do not forget that it must be multiplied by 1/r1 when integrating over r1!
        if (.not. allocated(grid_r1_r2%lambda_BB_r1_l_ij)) then
           call xermsg('BB_shell_mixed_integrals_mod','eval_BGBB_integrals','The lambda function has not been precalculated!',1,1)
        endif

        !Evaluate the double angular integrals for the CGTO on the r1 grid
        call cgto_pw%init_CGTO_shell_pw_expansion(this%cgto_shell,this%cgto_shell_starting_index)
        call cgto_pw%assign_grid(grid_r1_r2%r1,grid_r1_r2%w1)
        call cgto_pw%eval_CGTO_shell_pw_expansion(grid_r1_r2%bspline_grid%knots,grid_r1_r2%max_bspline_l,0,max_l)

        if (allocated(this%BGBB)) deallocate(this%BGBB)
        ind_21 = grid_r1_r2%last_bspline_inner*n_lbmb*this%cgto_shell%number_of_functions      
        ind_43 = (ipair(n_lbmb) + n_lbmb)*grid_r1_r2%n_unique_pairs
        n_lm = (max_l+1)**2
        allocate(this%BGBB(ind_43,ind_21),&
                 radial_BGBB(grid_r1_r2%n_unique_pairs,grid_r1_r2%last_bspline_inner,n_lm),&
                 radial_lm_is_significant(n_lm),stat=err)
        if (err /= 0) call xermsg('BB_shell_mixed_integrals_mod','eval_BGBB_integrals','Memory allocation failed.',err,1)

        !Evaluate the BGBB integrals for all shells of BTOs in the first coordinate and all unique combinations of BTO shells in the second coordinate:
        this%BGBB = 0.0_cfp
        do l1=0,grid_r1_r2%max_bspline_l
           do m1=-l1,l1
              l1m1 = l1*l1+l1+m1+1
              do Mg=1,this%cgto_shell%number_of_functions

                 !Evaluate the radial integrals:
                 !todo if the CGTO is on CMS then the radial integrals depend only on l and should calculated outside of the l1,m1,Mg loops.
                 radial_lm_is_significant = 0
                 do l=0,max_l
                    do m=-l,l
                       lm = l*l+l+m+1

                       p = cgto_pw%non_neg_indices_l_lp(Mg,l1m1,lm)
                       if (p > 0) then
                          radial_lm_is_significant(lm) = 1
                          do i=1,grid_r1_r2%last_bspline_inner
                             si = grid_r1_r2%bspline_start_end_r1(1,i)
                             ei = grid_r1_r2%bspline_start_end_r1(2,i)
                             do j=1,grid_r1_r2%n_unique_pairs
                                !int_{r1} r1**2 * 1/r1*lambda_BB_r1_l_ij * B(r1)/r1 * G_{lm,l1m1}(r1)
                                radial_BGBB(j,i,lm) = sum(grid_r1_r2%B_vals_r1(si:ei,i)&
                                                          *cgto_pw%gaunt_angular_integrals(si:ei,p)&
                                                          *grid_r1_r2%lambda_BB_r1_l_ij(si:ei,l,j))
                             enddo !j
                          enddo !i
                       endif

                    enddo !m
                 enddo !l

                 !Loop over the unique combinations of angular parts of the BB pair:
                 do l3=0,grid_r1_r2%max_bspline_l
                    do l4=0,l3

                       do l=abs(l3-l4),l3+l4
                          do m=-l,l
                             lm = l*l+l+m+1
                             if (radial_lm_is_significant(lm) .eq. 1) then
     
                                do m3=-l3,l3
                                   l3m3 = l3*l3+l3+m3+1
                                   do m4=-l4,l4
                                      l4m4 = l4*l4+l4+m4+1
                                      if (l4m4 > l3m3) cycle
                                      cf = cpl%rgaunt(l3,l4,l,m3,m4,m)
                                      if (cf .ne. 0.0_cfp) then
                                         l4m4 = l4*l4+l4+m4+1
                                         ind_lm = ipair(l3m3)+l4m4
                                         do i=1,grid_r1_r2%last_bspline_inner
                                            ind_21 = i + grid_r1_r2%last_bspline_inner*(Mg-1)&
                                                       + grid_r1_r2%last_bspline_inner*this%cgto_shell%number_of_functions*&
                                                         (l1m1-1) !3D array (i,Mg,l1m1)
                                            do j=1,grid_r1_r2%n_unique_pairs
                                               ind_43 = j + grid_r1_r2%n_unique_pairs*(ind_lm-1) !2D array (j,ind_lm)
                                               this%BGBB(ind_43,ind_21) = this%BGBB(ind_43,ind_21) + cf*radial_BGBB(j,i,lm)
                                            enddo !j
                                         enddo !i
                                      endif
                                   enddo !m4
                                enddo !m3
                             endif

                          enddo !m
                       enddo !l

                    enddo !l4
                 enddo !l3

              enddo !Mg
           enddo !m1
        enddo !l1

  end subroutine eval_BGBB_integrals

  !> Obtain the 2-electron integral (BG|BB). Assumes that eval_BGBB_integrals was called first.
  subroutine get_BGBB_shell_integrals(this,&
                                      bto_shell_1,&
                                      bto_shell_3,&
                                      bto_shell_4,&
                                      starting_index_1,&
                                      starting_index_3,&
                                      starting_index_4,&
                                      indexing_method,&
                                      two_el_column,&
                                      int_index,&
                                      integrals,&
                                      n_unique)
    use gto_routines_gbl, only: reorder_and_index_2el, index_2el
    use special_functions_gbl, only: ipair
    implicit none
    class(BB_shell_integrals_obj) :: this
    type(BTO_shell_data_obj), intent(in) :: bto_shell_1,bto_shell_3,bto_shell_4
    integer, intent(in) :: starting_index_1, starting_index_3, starting_index_4, two_el_column,indexing_method
    !We assume that these three arrays have been allocated to the appropriate dimensions:
    integer, allocatable :: int_index(:,:)
    real(kind=cfp), allocatable :: integrals(:,:)
    integer, intent(out) :: n_unique

    integer :: i, j, l1m1_min,l1m1_max,l3m3_min,l3m3_max,l4m4_min,l4m4_max,n_l4m4,n_l3m3,l1m1,Mg,ind_21,base,&
               l3m3,l4m4,ind_lm,ind,ind_43
    
       i = bto_shell_1%bspline_index
       j = grid_r1_r2%pairs_to_unique_pairs(bto_shell_3%bspline_index,bto_shell_4%bspline_index)

       n_l4m4 = 2*bto_shell_4%l+1
       n_l3m3 = 2*bto_shell_3%l+1

       n_unique = n_l4m4*n_l3m3*this%cgto_shell%number_of_functions*(2*bto_shell_1%l+1)

       if (j .ne. 0) then

          l1m1_min = bto_shell_1%l**2 + 1
          l1m1_max = (bto_shell_1%l+1)**2
   
          l3m3_min = bto_shell_3%l**2 + 1
          l3m3_max = (bto_shell_3%l+1)**2
   
          l4m4_min = bto_shell_4%l**2 + 1
          l4m4_max = (bto_shell_4%l+1)**2
   
          do l1m1=l1m1_min,l1m1_max
             do Mg=1,this%cgto_shell%number_of_functions
                ind_21 = i + grid_r1_r2%last_bspline_inner*(Mg-1)&
                           + grid_r1_r2%last_bspline_inner*this%cgto_shell%number_of_functions*(l1m1-1) !3D array (i,Mg,l1m1)
                base = n_l3m3*n_l4m4*(Mg-1) + n_l3m3*n_l4m4*this%cgto_shell%number_of_functions*(l1m1-l1m1_min)
   
                do l3m3=l3m3_min,l3m3_max
                   do l4m4=l4m4_min,l4m4_max
                      ind_lm = ipair(max(l3m3,l4m4)) + min(l3m3,l4m4)
                      ind_43 = j + grid_r1_r2%n_unique_pairs*(ind_lm-1) !2D array (j,ind_lm)
                      ind = l4m4-l4m4_min + 1 + n_l4m4*(l3m3-l3m3_min) + base !save in order (m4,m3,Mg,m1)
                      integrals(ind,two_el_column) = this%BGBB(ind_43,ind_21)
                   enddo !l4m4
                enddo !l3m3
             enddo !mg
          enddo !l1m1
       else
          !The pair of shells 3 and 4 does not overlap
          integrals(1:n_unique,two_el_column) = 0.0_cfp
       endif

       !Compute indices
       if (indexing_method .eq. 2) then
          ! THERE IS A MISTAKE IN THE NEW VERSION OF REORDERING
          !call reorder_and_index_2el(bto_shell_4%l,&
          !                           bto_shell_3%l,&
          !                           this%cgto_shell%l,&
          !                           bto_shell_1%l,&
          !                           starting_index_4,&
          !                           starting_index_3,&
          !                           this%cgto_shell_starting_index,&
          !                           starting_index_1,&
          !                           two_el_column,&
          !                           0,&
          !                           integrals,&
          !                           n_unique,&
          !                           int_index)
          call reorder_and_index_2el(bto_shell_4%l,&
                                     bto_shell_3%l,&
                                     this%cgto_shell%l,&
                                     bto_shell_1%l,&
                                     starting_index_4,&
                                     starting_index_3,&
                                     this%cgto_shell_starting_index,&
                                     starting_index_1,&
                                     two_el_column,&
                                     int_index,&
                                     integrals)
       else
          call index_2el(bto_shell_4%l,&
                         bto_shell_3%l,&
                         this%cgto_shell%l,&
                         bto_shell_1%l,&
                         starting_index_4,&
                         starting_index_3,&
                         this%cgto_shell_starting_index,&
                         starting_index_1,&
                         int_index,&
                         .false.,.false.)
       endif

  end subroutine get_BGBB_shell_integrals

  !> Assumes that the lambda function has been precalculated before using a call to this%get_lambda function.
  subroutine eval_radial_BBBB_integrals(this)
     use omp_lib
     implicit none
     class(BB_shell_integrals_obj) :: this

     integer :: err, i, j, l, si, ei
     real(kind=cfp), allocatable :: inv_r1(:)
     real(kind=wp) :: t1, t2
     logical :: in_parallel

        in_parallel = omp_in_parallel()

        if (.not.(in_parallel)) write(stdout,'(/,"--------->","BB_shell_mixed_integrals_mod:eval_radial_BBBB_integrals")')

        t1 = omp_get_wtime()
 
        this%max_l = 2*grid_r1_r2%max_bspline_l
        this%max_bspline_l = grid_r1_r2%max_bspline_l
        this%n_bsplines = grid_r1_r2%last_bspline_inner

        call cpl%prec_cgaunt(this%max_l)

        !Check that lambda function has been precalculated using call grid_r1_r2%eval_lambda_l_BB(this%max_l). Do not forget that below it must be multiplied by 1/r1 when integrating over r1!
        if (.not. allocated(grid_r1_r2%lambda_BB_r1_l_ij)) then
           call xermsg('BB_shell_mixed_integrals_mod','eval_radial_BBBB_integrals',&
                       'The lambda function has not been precalculated!',1,1)
        endif

        if (allocated(this%radial_BBBB)) deallocate(this%radial_BBBB)
        allocate(inv_r1(grid_r1_r2%n1_total_points),this%radial_BBBB(0:max(this%max_l,1),&
                 grid_r1_r2%n_unique_pairs,grid_r1_r2%n_unique_pairs),stat=err)
        if (err /= 0) call xermsg('BB_shell_mixed_integrals_mod','eval_radial_BBBB_integrals','Memory allocation failed.',err,1)

        do i=1,grid_r1_r2%n1_total_points
           inv_r1(i) = 1.0_cfp/grid_r1_r2%r1(i)
        enddo !i

        !Loop over all unique pairs of B-splines in the r1 coordinate
        this%radial_BBBB = 0.0_cfp
        !$OMP PARALLEL IF (.not.(in_parallel)) DEFAULT(NONE) PRIVATE(i,si,ei,l,j) SHARED(this,grid_r1_r2,inv_r1)
        !$OMP DO 
        do i=1,grid_r1_r2%n_unique_pairs
           si = grid_r1_r2%BB_start_end_r1(1,i)
           ei = grid_r1_r2%BB_start_end_r1(2,i)
           do l=0,this%max_l
             do j=1,grid_r1_r2%n_unique_pairs
                !int_{r1} r1**2 * 1/r1*lambda_BB_r1_l_ij * B(r1)/r1 * B(r1)/r1
                this%radial_BBBB(l,i,j) = sum(grid_r1_r2%BB_vals_r1(si:ei,i)&
                                              *grid_r1_r2%lambda_BB_r1_l_ij(si:ei,l,j)&
                                              *inv_r1(si:ei))
             enddo !j
           enddo !l
        enddo !i
        !$OMP END DO
        !$OMP END PARALLEL

        t2 = omp_get_wtime()

        this%radial_BBBB_calculated = .true.

        if (.not.(in_parallel)) then
           write(stdout,'("<---------","BB_shell_mixed_integrals_mod:eval_radial_BBBB_integrals took [s]: ",f15.5)') t2-t1
        endif
     
  end subroutine eval_radial_BBBB_integrals

  subroutine get_radial_BBBB(this,radial_BBBB)
    implicit none
    class(BB_shell_integrals_obj) :: this
    real(kind=cfp), allocatable :: radial_BBBB(:)

    integer :: err, i, j, l, ind

      if (.not. this%radial_BBBB_calculated) then
         call xermsg('BB_shell_mixed_integrals_mod','get_radial_BBBB',&
                     'The radial integrals have not been calculated. Run eval_radial_BBBB_integrals first.',1,1)
      endif

      if (allocated(radial_BBBB)) deallocate(radial_BBBB)
      i = size(this%radial_BBBB,1)*size(this%radial_BBBB,2)*size(this%radial_BBBB,3)

      allocate(radial_BBBB(i),stat=err)
      if (err /= 0) call xermsg('BB_shell_mixed_integrals_mod','get_radial_BBBB','Memory allocation failed.',err,1)

      !Transform the 3D array into a 1D array
      ind = 0
      do j=1,size(this%radial_BBBB,3)
         do i=1,size(this%radial_BBBB,2)
            do l=0,size(this%radial_BBBB,1)-1
               ind = ind + 1
               radial_BBBB(ind) = this%radial_BBBB(l,i,j)
            enddo
         enddo
      enddo

  end subroutine get_radial_BBBB

  !> Obtain the 2-electron integral (BB|BB). Assumes that radial_BBBB_integrals was called first.
  subroutine get_BBBB_shell_integrals(this,&
                                      bto_shell_1,&
                                      bto_shell_2,&
                                      bto_shell_3,&
                                      bto_shell_4,&
                                      starting_index_1,&
                                      starting_index_2,&
                                      starting_index_3,&
                                      starting_index_4,&
                                      indexing_method,&
                                      two_el_column,&
                                      int_index,&
                                      integrals,&
                                      n_unique)
    use gto_routines_gbl, only: reorder_and_index_2el, index_2el
    use const_gbl, only: level3 !debug
    implicit none
    class(BB_shell_integrals_obj) :: this
    type(BTO_shell_data_obj), intent(in) :: bto_shell_1,bto_shell_2,bto_shell_3,bto_shell_4
    integer, intent(in) :: starting_index_1, starting_index_2, starting_index_3, starting_index_4, two_el_column,indexing_method
    !We assume that these three arrays have been allocated to the appropriate dimensions:
    integer, allocatable :: int_index(:,:)
    real(kind=cfp), allocatable :: integrals(:,:)
    integer, intent(out) :: n_unique

    integer :: min_l, max_l, i, j, l, m, lm, lm_max, m1m2, m3m4, ind, n_shell_12, n_shell_34, err, m1, m2, m3, m4
    real(kind=cfp), allocatable :: couplings12(:,:), couplings34(:,:)

      if (.not. this%radial_BBBB_calculated) then
         call xermsg('BB_shell_mixed_integrals_mod','get_BBBB_shell_integrals',&
                     'The radial integrals have not been calculated. Run eval_radial_BBBB_integrals first.',1,1)
      endif

      min_l = max(abs(bto_shell_1%l-bto_shell_2%l),abs(bto_shell_3%l-bto_shell_4%l))
      max_l = min(bto_shell_1%l+bto_shell_2%l,bto_shell_3%l+bto_shell_4%l)

      i = grid_r1_r2%pairs_to_unique_pairs(bto_shell_1%bspline_index,bto_shell_2%bspline_index)
      j = grid_r1_r2%pairs_to_unique_pairs(bto_shell_3%bspline_index,bto_shell_4%bspline_index)

      n_unique =  bto_shell_1%number_of_functions&
                 *bto_shell_2%number_of_functions&
                 *bto_shell_3%number_of_functions&
                 *bto_shell_4%number_of_functions

      if (((i .ne. 0) .and. (j .ne. 0)) .and. (min_l .le. max_l)) then

         n_shell_12 = bto_shell_1%number_of_functions*bto_shell_2%number_of_functions
         n_shell_34 = bto_shell_3%number_of_functions*bto_shell_4%number_of_functions
         lm_max = (max_l+1)**2 - min_l**2
         allocate(couplings12(lm_max,n_shell_12),couplings34(lm_max,n_shell_34),stat=err)
         if (err /= 0) call xermsg('BB_shell_mixed_integrals_mod','get_BBBB_shell_integrals','Memory allocation failed.',err,1)

         m1m2 = 0
         do m1=-bto_shell_1%l,bto_shell_1%l
            do m2=-bto_shell_2%l,bto_shell_2%l
               m1m2 = m1m2 + 1
  
               lm = 0
               do l=min_l,max_l
                  do m=-l,l
                     lm = lm + 1
                     couplings12(lm,m1m2) = cpl%rgaunt(bto_shell_1%l,bto_shell_2%l,l,m1,m2,m)*this%radial_BBBB(l,i,j)
                  enddo !m
               enddo !l
  
            enddo !m2
         enddo !m1
  
         m3m4 = 0
         do m3=-bto_shell_3%l,bto_shell_3%l
            do m4=-bto_shell_4%l,bto_shell_4%l
               m3m4 = m3m4 + 1
  
               lm = 0
               do l=min_l,max_l
                  do m=-l,l
                     lm = lm + 1
                     couplings34(lm,m3m4) = cpl%rgaunt(bto_shell_3%l,bto_shell_4%l,l,m3,m4,m)
                  enddo !m
               enddo !l
  
            enddo !m4
         enddo !m3

         !todo matrix * matrix multiplication
         ind = 0
         m1m2 = 0
         do m1=-bto_shell_1%l,bto_shell_1%l
            do m2=-bto_shell_2%l,bto_shell_2%l
               m1m2 = m1m2 + 1
               m3m4 = 0
               do m3=-bto_shell_3%l,bto_shell_3%l
                  do m4=-bto_shell_4%l,bto_shell_4%l
                     m3m4 = m3m4 + 1

                     ind = ind + 1

                     integrals(ind,two_el_column) = sum(couplings12(1:lm_max,m1m2)*couplings34(1:lm_max,m3m4))
  
                  enddo !m4
               enddo !m3
            enddo !m2
         enddo !m1

      else !at least one of the two pairs of B-splines don't overlap or angular momentum selection rules forbid the integral
         integrals = 0.0_cfp
      endif

      !debug
      !write(level3,'("Inside get_BBBB_shell_integrals")')
      !write(level3,'("B-spline indices",4i5)') bto_shell_1%bspline_index, bto_shell_2%bspline_index, bto_shell_3%bspline_index, &
      !                                         bto_shell_4%bspline_index
      !write(level3,'("The integrals: ",100e25.15)') integrals(:,two_el_column)                                            

      !Compute indices
      if (indexing_method .eq. 2) then
         ! THERE IS A MISTAKE IN THE NEW VERSION OF REORDERING
         !call reorder_and_index_2el(bto_shell_4%l,bto_shell_3%l,bto_shell_2%l,bto_shell_1%l,&
         !&starting_index_4,starting_index_3,starting_index_2,starting_index_1,two_el_column,0,integrals,n_unique,int_index)
         call reorder_and_index_2el(bto_shell_4%l,&
                                    bto_shell_3%l,&
                                    bto_shell_2%l,&
                                    bto_shell_1%l,&
                                    starting_index_4,&
                                    starting_index_3,&
                                    starting_index_2,&
                                    starting_index_1,&
                                    two_el_column,&
                                    int_index,&
                                    integrals)
      else
         call index_2el(bto_shell_4%l,bto_shell_3%l,bto_shell_2%l,bto_shell_1%l,&
         &starting_index_4,starting_index_3,starting_index_2,starting_index_1,int_index,.false.,.false.)
      endif

  end subroutine get_BBBB_shell_integrals

  subroutine prepare_for_reconstruction(this,&
                                        bspline_grid,&
                                        max_bspline_l,&
                                        offset_BBBB,&
                                        n_cgto_shells,&
                                        BGBB_aux_offset,&
                                        bto_starting_indices,&
                                        bspline_range,&
                                        cgto_starting_indices,&
                                        cgto_shell_n_fns)
    use bspline_grid_gbl
    implicit none
    class(BB_shell_integrals_obj) :: this
    class(bspline_grid_obj) :: bspline_grid
    integer, intent(in) :: max_bspline_l, n_cgto_shells, BGBB_aux_offset(n_cgto_shells), bto_starting_indices(:,0:), &
                           bspline_range(:,0:), cgto_starting_indices(:), cgto_shell_n_fns(:), offset_BBBB

    integer :: i, j, l, err, ind

       !Get the mapping between each pair of B-splines and the unique set of overlapping B-splines
       call bspline_grid%get_unique_pairs(this%pairs_to_unique_pairs,this%n_unique_pairs)

       this%max_l = 2*max_bspline_l
       this%max_bspline_l = max_bspline_l
       this%n_bsplines = bspline_grid%get_last_inner_bspline(bspline_grid%C) !%n
       call cpl%prec_cgaunt(this%max_l)
       cpl%always_precalculate = .true.

       this%offset_BBBB = offset_BBBB
       this%n_cgto_shells = n_cgto_shells

       this%n_lbmb = (this%max_bspline_l+1)**2
       this%n_lm = (this%max_l+1)**2

       !Offsets for each CGTO shell of the start of auxiliary integrals in the integrals array needed to assemble the BGBB class of integrals.
       if (n_cgto_shells > 0) then
          if (allocated(this%BGBB_aux_offset)) deallocate(this%BGBB_aux_offset)
          if (allocated(this%cgto_starting_indices)) deallocate(this%cgto_starting_indices)
          if (allocated(this%cgto_shell_n_fns)) deallocate(this%cgto_shell_n_fns)
          allocate(this%BGBB_aux_offset,source=BGBB_aux_offset,stat=err)
          if (err /= 0) then
             call xermsg('BB_shell_mixed_integrals_mod','prepare_for_reconstruction','Memory allocation 2a failed.',err,1)
          endif
          allocate(this%cgto_starting_indices,source=cgto_starting_indices,stat=err)
          if (err /= 0) then
             call xermsg('BB_shell_mixed_integrals_mod','prepare_for_reconstruction','Memory allocation 2b failed.',err,1)
          endif
          allocate(this%cgto_shell_n_fns,source=cgto_shell_n_fns,stat=err)
          if (err /= 0) then
             call xermsg('BB_shell_mixed_integrals_mod','prepare_for_reconstruction','Memory allocation 2c failed.',err,1)
          endif
       endif

       if (allocated(this%bto_starting_indices)) deallocate(this%bto_starting_indices)
       allocate(this%bto_starting_indices,source=bto_starting_indices,stat=err)
       if (err /= 0) call xermsg('BB_shell_mixed_integrals_mod','prepare_for_reconstruction','Memory allocation 3 failed.',err,1)

       if (allocated(this%bspline_range)) deallocate(this%bspline_range)
       allocate(this%bspline_range,source=bspline_range,stat=err)
       if (err /= 0) call xermsg('BB_shell_mixed_integrals_mod','prepare_for_reconstruction','Memory allocation 4 failed.',err,1)

       this%reconstruction_prepared = .true.

  end subroutine prepare_for_reconstruction

  !> WARNING: Note that this routine does not zero-out the arrays reconstructed and integrals since it is assumed this might be used also in a routine which reconstructs
  !> the BBBB class evaluating other elements of the array integrals, see routine omp_two_p_transform_pqrs_block_to_ijrs_AO.
  subroutine get_BGBB_integrals_all_pq_fixed_BG(this,m_r,bto_shell_r,m_s,cgto_shell_s,s_shell,ao_integrals,two_el_column,&
                                                reconstructed,all_zero,integrals)
    use special_functions_gbl, only: ipair
    implicit none
    class(BB_shell_integrals_obj) :: this
    type(BTO_shell_data_obj), intent(in) :: bto_shell_r
    type(CGTO_shell_data_obj), intent(in) :: cgto_shell_s
    integer, intent(in) :: m_r,m_s,s_shell,two_el_column
    logical, intent(out) :: all_zero
    real(kind=cfp), pointer :: ao_integrals(:,:)
    real(kind=cfp), intent(out) :: integrals(:)
    logical, intent(out) :: reconstructed(:)

    integer :: l1m1,offset,i,n1,n2,n3,n4,base,lm,l3,l4,l,m,m3,m4,l3m3,l4m4,i3,i4,i3i4,p,q,pq
    real(kind=cfp) :: cf
    logical :: first_lm

       if (.not. this%reconstruction_prepared) then
          call xermsg('BB_shell_mixed_integrals_mod','get_BBBB_integrals_all_pq_fixed_BG',&
                      'Reconstruction not prepared: call prepare_for_reconstruction first.',1,1)
       endif

       !starting position in the ao_integrals array of the auxiliary integrals needed to assemble the BGBB class for the given CGTO m.
       offset = this%BGBB_aux_offset(s_shell)

       !Evaluate the BGBB auxiliary integrals for a fixed combination of BTO/CGTO indices for one electron and all unique combinations of a pair of BTO indices for the second electron.
       l1m1 = bto_shell_r%l*bto_shell_r%l+bto_shell_r%l+m_r+1
       i = bto_shell_r%bspline_index
       n1 = this%n_unique_pairs
       n2 = n1*this%n_lm
       n3 = n2*this%n_bsplines
       n4 = n3*this%n_lbmb

       !Loop over the unique combinations of angular parts of the BB pair:
       all_zero = .false.
       do l3=0,this%max_bspline_l
          do l4=0,l3

             first_lm = .true.
             do l=abs(l3-l4),l3+l4
                do m=-l,l
                   lm = l*l+l+m+1
                   base = n1*(lm-1) + n2*(i-1) + n3*(l1m1-1) + n4*(m_s+cgto_shell_s%l) + offset

                   do m3=-l3,l3
                      l3m3 = l3*l3+l3+m3+1
                      do m4=-l4,l4
                         l4m4 = l4*l4+l4+m4+1
                         if (l4m4 > l3m3) cycle
                         cf = cpl%rgaunt(l3,l4,l,m3,m4,m)
                         if (cf .ne. 0.0_cfp) then
                            !todo the list of i3,i4 can be determined at the top of the l3,l4 loops!
                            do i3=this%bspline_range(1,l3),this%bspline_range(2,l3)
                               p = this%bto_starting_indices(i3,l3) + m3+l3
                               do i4=this%bspline_range(1,l4),this%bspline_range(2,l4)
                                  q = this%bto_starting_indices(i4,l4) + m4+l4
                                  if (q > p) cycle !this is important since we need to sum only the unique contributions below
                                  pq = ipair(p)+q
                                  i3i4 = this%pairs_to_unique_pairs(i3,i4)
                                  if (first_lm) integrals(pq) = 0.0_cfp !Make sure the integrals are zeroed-out before accumulation
                                  if (i3i4 > 0) then !the B-splines overlap
                                     !Note that here we rely on the pq elements being zeroed out outside of this routine.
                                     integrals(pq) = integrals(pq) + cf*ao_integrals(base+i3i4,two_el_column)
                                  endif
                                  reconstructed(pq) = .true.
                               enddo !j
                            enddo !i
                         endif
                      enddo !m4
                   enddo !m3

                   if (first_lm) first_lm = .false.
                enddo !m
             enddo !l

          enddo !l4
       enddo !l3

  end subroutine get_BGBB_integrals_all_pq_fixed_BG

  !> WARNING: Note that this routine does not zero-out the arrays reconstructed and integrals since it is assumed this might be used also in a routine which reconstructs
  !> the BBBB class evaluating other elements of the array integrals, see routine omp_two_p_transform_pqrs_block_to_ijrs_AO.
  subroutine get_BGBB_integrals_all_pq_fixed_BB(this,m_r,bto_shell_r,m_s,bto_shell_s,ao_integrals,two_el_column,&
                                                reconstructed,all_zero,integrals)
    use special_functions_gbl, only: ipair
    implicit none
    class(BB_shell_integrals_obj) :: this
    type(BTO_shell_data_obj), intent(in) :: bto_shell_r,bto_shell_s
    integer, intent(in) :: m_r,m_s,two_el_column
    logical, intent(out) :: all_zero
    real(kind=cfp), pointer :: ao_integrals(:,:)
    real(kind=cfp), intent(out) :: integrals(:)
    logical, intent(out) :: reconstructed(:)

    integer :: l1m1,offset,i,n1,n2,n3,n4,base,lm,l,m,p,q,pq,l_min,l_max,l1,m1,Mg,cgto_shell,rs_pair,err,nz,j
    real(kind=cfp) :: cf
    real(kind=cfp), allocatable :: couplings(:)
    integer, allocatable :: lm_list(:)

       if (.not. this%reconstruction_prepared) then
          call xermsg('BB_shell_mixed_integrals_mod','get_BGBB_integrals_all_pq_fixed_BB',&
                      'Reconstruction not prepared: call prepare_for_reconstruction first.',1,1)
       endif

       !Evaluate the BGBB auxiliary integrals for a fixed pair of BTO indices for one electron and all combinations of BTO/CGTO indices for the second electron.
       n1 = this%n_unique_pairs
       n2 = n1*this%n_lm
       n3 = n2*this%n_bsplines
       n4 = n3*this%n_lbmb

       rs_pair = this%pairs_to_unique_pairs(bto_shell_r%bspline_index,bto_shell_s%bspline_index)
       l_min = abs(bto_shell_r%l-bto_shell_s%l)
       l_max = bto_shell_r%l+bto_shell_s%l
       i = (l_max+1)**2 - l_min**2
       allocate(couplings(i),lm_list(i),stat=err)
       if (err /= 0) then
          call xermsg('BB_shell_mixed_integrals_mod','get_BGBB_integrals_all_pq_fixed_BB','Memory allocation failed.',err,1)
       endif

       nz = 0
       lm = 0
       do l=l_min,l_max
          do m=-l,l
             lm = lm + 1
             cf = cpl%rgaunt(bto_shell_r%l,bto_shell_s%l,l,m_r,m_s,m)
             if (cf .ne. 0.0_cfp) then
                nz = nz + 1
                couplings(nz) = cf
                lm_list(nz) = lm
             endif
          enddo !m
       enddo !l

       all_zero = .false.
       do cgto_shell=1,this%n_cgto_shells
          offset = this%BGBB_aux_offset(cgto_shell)
          do l1=0,this%max_bspline_l
             do m1=-l1,l1
                l1m1 = l1*l1+l1+m1+1
                do Mg=1,this%cgto_shell_n_fns(cgto_shell)
                   q = this%cgto_starting_indices(cgto_shell) + Mg-1

                   do j=1,nz
                      cf = couplings(j)
                      lm = lm_list(j)
                      base = rs_pair + n1*(lm-1) + n3*(l1m1-1) + n4*(Mg-1) + offset
                      if (j .eq. 1) then
                         do i=this%bspline_range(1,l1),this%bspline_range(2,l1)
                            !Note that the case q > p will never arise since CGTOs and BTOs are properly ordered in the basis so CGTOs are always before BTOs 
                            p = this%bto_starting_indices(i,l1) + m1+l1
                            pq = ipair(p)+q
                            !The first term so don't add: make sure the integrals array are accumulated with zero starting value.
                            integrals(pq) = cf*ao_integrals(base+n2*(i-1),two_el_column)
                            reconstructed(pq) = .true.
                         enddo !i
                      else
                         do i=this%bspline_range(1,l1),this%bspline_range(2,l1)
                            !Note that the case q > p will never arise since CGTOs and BTOs are properly ordered in the basis so CGTOs are always before BTOs 
                            p = this%bto_starting_indices(i,l1) + m1+l1
                            pq = ipair(p)+q
                            !Accumulate the contributions:
                            integrals(pq) = integrals(pq) + cf*ao_integrals(base+n2*(i-1),two_el_column)
                            reconstructed(pq) = .true.
                         enddo !i
                      endif
                   enddo !j

                enddo !Mg
             enddo !m1
          enddo !l1
       enddo !cgto_shell

  end subroutine get_BGBB_integrals_all_pq_fixed_BB

  !> WARNING: Note that this routine does not zero-out the arrays reconstructed and integrals since it is assumed this might be used also in a routine which reconstructs
  !> the BGBB class evaluating other elements of the array integrals, see routine omp_two_p_transform_pqrs_block_to_ijrs_AO.
  subroutine get_BBBB_integrals_all_pq_fixed_BB(this,m_r,m_s,bto_shell_r,bto_shell_s,ao_integrals,two_el_column,&
                                                reconstructed,all_zero,integrals)
    use special_functions_gbl, only: ipair
    implicit none
    class(BB_shell_integrals_obj) :: this
    type(BTO_shell_data_obj), intent(in) :: bto_shell_r,bto_shell_s
    integer, intent(in) :: m_r,m_s,two_el_column
    real(kind=cfp), pointer :: ao_integrals(:,:)
    logical, intent(out) :: all_zero
    real(kind=cfp), intent(out) :: integrals(:)
    logical, intent(out) :: reconstructed(:)

    integer :: rs_pair, min_l, max_l, l, m, lm, lp, lq, mp, mq, mpmq, i, j, ij, p, q, pq, min_l_new,max_l_new,lm_min_rs,&
               err,n1,n2,base,nz
    real(kind=cfp), allocatable :: couplings_rs(:,:), couplings_pq(:,:)
    integer, allocatable :: lm_list(:,:)
    real(kind=cfp) :: cf
    logical :: zero

       if (.not. this%reconstruction_prepared) then
          call xermsg('BB_shell_mixed_integrals_mod','get_BBBB_integrals_all_pq_fixed_BB',&
                      'Reconstruction not prepared: call prepare_for_reconstruction first.',1,1)
       endif
    
       rs_pair = this%pairs_to_unique_pairs(bto_shell_r%bspline_index,bto_shell_s%bspline_index)

       all_zero = .false.
       if (rs_pair .eq. 0) then !The rs-pair of radial B-splines does not overlap: all integrals are zero
          all_zero = .true.

          do lp=0,this%max_bspline_l
             do lq=0,lp
                do i=this%bspline_range(1,lp),this%bspline_range(2,lp)
                   do j=this%bspline_range(1,lq),this%bspline_range(2,lq)
                      do mp=0,2*lp
                         p = this%bto_starting_indices(i,lp) + mp
                         do mq=0,2*lq
                            q = this%bto_starting_indices(j,lq) + mq
                            pq = ipair(max(p,q))+min(p,q)
                            integrals(pq) = 0.0_cfp
                            reconstructed(pq) = .true.
                         enddo !mq
                      enddo !mp
                   enddo !j
                enddo !i
             enddo
          enddo

          return
       endif

       !Below the element this%bto_starting_indices(i,l) contains the starting index in the shell of BTOs with radial B-spline i and angular momentum l.

       min_l = abs(bto_shell_r%l-bto_shell_s%l)
       max_l = bto_shell_r%l+bto_shell_s%l

       mpmq = (2*this%max_bspline_l+1)**2
       nz = (max_l-min_l+1)*4
       allocate(couplings_rs(nz,this%n_unique_pairs),lm_list(2,nz),couplings_pq(nz,mpmq),stat=err)
       if (err /= 0) then
          call xermsg('BB_shell_mixed_integrals_mod','get_BBBB_integrals_all_pq_fixed_BB','Memory allocation failed.',err,1)
       endif

       n1 = this%max_l+1
       n2 = n1*this%n_unique_pairs
       nz = 0
       do l=min_l,max_l
          base = l+1 + n2*(rs_pair-1) + this%offset_BBBB
          do m=-l,l
             !For a fixed set of bto_shell_r%l,bto_shell_s%l,l,m_r,m_s there can be at most 4 allowed values of m hence the factor 4 in the allocation statement above.
             cf = cpl%rgaunt(bto_shell_r%l,bto_shell_s%l,l,m_r,m_s,m)
             if (cf .ne. 0.0_cfp) then
                nz = nz + 1
                lm_list(1:2,nz) = (/l,m/)
                do ij=1,this%n_unique_pairs
                   couplings_rs(nz,ij) = cf*ao_integrals(base+n1*(ij-1),two_el_column) !3D array indexing: offset + (l,ij,rs_pair) = l+1 + n1*(ij-1) + n2*(rs_pair-1) + offset
                enddo !ij
             endif
          enddo !m
       enddo !l

       do lp=0,this%max_bspline_l
          do lq=0,lp

             min_l_new = max(min_l,abs(lp-lq))
             max_l_new = min(max_l,lp+lq)

             if (min_l_new > max_l_new) then !the pq integrals are zero due to angular momentum selection rules
                do i=this%bspline_range(1,lp),this%bspline_range(2,lp)
                   do j=this%bspline_range(1,lq),this%bspline_range(2,lq)
                      do mp=0,2*lp
                         p = this%bto_starting_indices(i,lp) + mp
                         do mq=0,2*lq
                            q = this%bto_starting_indices(j,lq) + mq
                            pq = ipair(max(p,q))+min(p,q)
                            integrals(pq) = 0.0_cfp
                            reconstructed(pq) = .true.
                         enddo !mq
                      enddo !mp
                   enddo !j
                enddo !i

                cycle
             endif

             mpmq = 0
             do mp=-lp,lp
                do mq=-lq,lq
                   mpmq = mpmq + 1
        
                   do lm=1,nz
                      l = lm_list(1,lm)
                      m = lm_list(2,lm)
                      if (l .ge. min_l_new .and. l .le. max_l_new) then
                         couplings_pq(lm,mpmq) = cpl%rgaunt(lp,lq,l,mp,mq,m)
                      else
                         couplings_pq(lm,mpmq) = 0.0_cfp
                      endif
                   enddo

                enddo !m2
             enddo !m1

             ij = 0
             do i=this%bspline_range(1,lp),this%bspline_range(2,lp)
                do j=this%bspline_range(1,lq),this%bspline_range(2,lq)
                   ij = this%pairs_to_unique_pairs(i,j)
                   if (ij .le. 0) then
                      zero = .true. !the B-splines don't overlap
                   else
                      zero = .false.
                   endif
                   mpmq = 0
                   do mp=0,2*lp
                      p = this%bto_starting_indices(i,lp) + mp
                      do mq=0,2*lq
                         mpmq = mpmq + 1
                         q = this%bto_starting_indices(j,lq) + mq
                         pq = ipair(max(p,q))+min(p,q)
                         if (zero) then
                            integrals(pq) = 0.0_cfp
                         else
                            !todo matrix*matrix multiplication
                            integrals(pq) = sum(couplings_pq(1:nz,mpmq)*couplings_rs(1:nz,ij))
                         endif
                         reconstructed(pq) = .true.
                      enddo !mq
                   enddo !mp
                enddo !j
             enddo !i

          enddo !lq
       enddo !lp

  end subroutine get_BBBB_integrals_all_pq_fixed_BB

  !> Obtain the shell of 2-electron integrals (BB|BB). Assumes that prepare_for_reconstruction was called first.
  subroutine reconstruct_BBBB_shell_integrals(this,bto_shell_1,bto_shell_2,bto_shell_3,bto_shell_4,&
                                              starting_index_1,starting_index_2,starting_index_3,starting_index_4,&
                                              offset,integrals,n_unique)
    use gto_routines_gbl, only: reorder_and_index_2el
    implicit none
    class(BB_shell_integrals_obj) :: this
    type(BTO_shell_data_obj), intent(in) :: bto_shell_1,bto_shell_2,bto_shell_3,bto_shell_4
    integer, intent(in) :: starting_index_1, starting_index_2, starting_index_3, starting_index_4,offset
    real(kind=cfp), allocatable :: integrals(:,:)
    integer, intent(out) :: n_unique

    integer :: min_l, max_l, i, j, l, m, lm, lm_max, m1m2, m3m4, ind, n_shell_12, n_shell_34, err, m1, m2, m3, m4
    real(kind=cfp), allocatable :: couplings12(:,:), couplings34(:,:)
    integer, parameter :: two_el_column = 1

      if (.not. this%reconstruction_prepared) then
         call xermsg('BB_shell_mixed_integrals_mod','reconstruct_BBBB_shell_integrals',&
                     'Reconstruction not prepared: call prepare_for_reconstruction first.',1,1)
      endif

      stop "wouldnt work since this%radial_BBBB is no longer being constructed in prepare_for_reconstruction."

      min_l = max(abs(bto_shell_1%l-bto_shell_2%l),abs(bto_shell_3%l-bto_shell_4%l))
      max_l = min(bto_shell_1%l+bto_shell_2%l,bto_shell_3%l+bto_shell_4%l)

      if (max_l > this%max_l) call xermsg('BB_shell_mixed_integrals_mod','reconstruct_BBBB_shell_integrals',&
      &'Attempt to get integrals for BTO angular momentum for which the reconstruction has not been prepared .',2,1)

      i = this%pairs_to_unique_pairs(bto_shell_1%bspline_index,bto_shell_2%bspline_index)
      j = this%pairs_to_unique_pairs(bto_shell_3%bspline_index,bto_shell_4%bspline_index)

      n_unique = bto_shell_1%number_of_functions&
                 *bto_shell_2%number_of_functions&
                 *bto_shell_3%number_of_functions&
                 *bto_shell_4%number_of_functions

      if (((i .ne. 0) .and. (j .ne. 0)) .and. (min_l .le. max_l)) then

         n_shell_12 = bto_shell_1%number_of_functions*bto_shell_2%number_of_functions
         n_shell_34 = bto_shell_3%number_of_functions*bto_shell_4%number_of_functions
         lm_max = (max_l+1)**2 - min_l**2
         allocate(couplings12(lm_max,n_shell_12),couplings34(lm_max,n_shell_34),stat=err)
         if (err /= 0) then
            call xermsg('BB_shell_mixed_integrals_mod','reconstruct_BBBB_shell_integrals','Memory allocation failed.',err,1)
         endif

         m1m2 = 0
         do m1=-bto_shell_1%l,bto_shell_1%l
            do m2=-bto_shell_2%l,bto_shell_2%l
               m1m2 = m1m2 + 1
  
               lm = 0
               do l=min_l,max_l
                  do m=-l,l
                     lm = lm + 1
                     couplings12(lm,m1m2) = cpl%rgaunt(bto_shell_1%l,bto_shell_2%l,l,m1,m2,m)*this%radial_BBBB(l,i,j)
                  enddo !m
               enddo !l
  
            enddo !m2
         enddo !m1
  
         m3m4 = 0
         do m3=-bto_shell_3%l,bto_shell_3%l
            do m4=-bto_shell_4%l,bto_shell_4%l
               m3m4 = m3m4 + 1
  
               lm = 0
               do l=min_l,max_l
                  do m=-l,l
                     lm = lm + 1
                     couplings34(lm,m3m4) = cpl%rgaunt(bto_shell_3%l,bto_shell_4%l,l,m3,m4,m)
                  enddo !m
               enddo !l
  
            enddo !m4
         enddo !m3

         !todo matrix * matrix multiplication
         ind = 0
         m1m2 = 0
         do m1=-bto_shell_1%l,bto_shell_1%l
            do m2=-bto_shell_2%l,bto_shell_2%l
               m1m2 = m1m2 + 1
               m3m4 = 0
               do m3=-bto_shell_3%l,bto_shell_3%l
                  do m4=-bto_shell_4%l,bto_shell_4%l
                     m3m4 = m3m4 + 1

                     ind = ind + 1

                     integrals(offset+ind,two_el_column) = sum(couplings12(1:lm_max,m1m2)*couplings34(1:lm_max,m3m4))
  
                  enddo !m4
               enddo !m3
            enddo !m2
         enddo !m1

      else !at least one of the two pairs of B-splines don't overlap or angular momentum selection rules forbid the integral
         integrals(offset+1:offset+n_unique,two_el_column) = 0.0_cfp
         n_unique = 0
      endif

      !Reorder the integrals
      ! THERE IS A MISTAKE IN THE NEW VERSION OF REORDERING
      call reorder_and_index_2el(bto_shell_4%l,&
                                 bto_shell_3%l,&
                                 bto_shell_2%l,&
                                 bto_shell_1%l,&
                                 starting_index_4,&
                                 starting_index_3,&
                                 starting_index_2,&
                                 starting_index_1,&
                                 two_el_column,&
                                 offset,&
                                 integrals,&
                                 n_unique)

  end subroutine reconstruct_BBBB_shell_integrals

end module BB_shell_mixed_integrals_gbl
