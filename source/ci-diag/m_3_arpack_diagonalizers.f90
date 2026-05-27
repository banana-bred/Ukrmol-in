! Copyright 2019
!
! For a comprehensive list of the developers that contributed to these codes
! see the UK-AMOR website.
!
! This file is part of UKRmol-in (UKRmol+ suite).
!
!     UKRmol-in is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     UKRmol-in is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  UKRmol-in (in source/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!====================================================================!
! Author: Pavlos G. Galiatsatos.                                     !
! Corrections, modifications and improvements: Michael Lysaght.      !
! Date: 2011/06/30.                                                  !
!====================================================================!

module m_3_arpack_diagonalizers

  use precisn
  use m_2_matrix_sparse

  implicit none

  contains

!====================================================================!
!                                                                    !
! subroutine: "arpack_coo_half_auto".                                !
!                                                                    !
! a. input: row_index_array,                                         !
!           column_index_array,                                      !
!           ham_elem_array,                                          !
!           which,                                                   !
!           ncv_coeff,                                               !
!           nev,                                                     !
!           maxitr,                                                  !
!           tol                                                      !
!                                                                    !
! b. output: vectors_array,                                          !
!            values_array                                            !
!                                                                    !
! c. purpose: diagonalization of the hamiltonian matrix.             !                    
!                                                                    !
!====================================================================!

  subroutine arpack_coo_half_auto( row_index_array,      &  !   1. input & deallocated.
                                   column_index_array,   &  !   2. input & deallocated.
                                   ham_elem_array,       &  !   3. input & deallocated.
                                   matrix_dimension,     &  !   4. input.
                                   non_zero_ham_elem,    &  !   5. input.
                                   which,                &  !   6. input.
                                   ncv_coeff,            &  !   7. input.
                                   num_eigvals,          &  !   8. input & output.
                                   maxitr,               &  !   9. input & output.
                                   tol,                  &  !  10. input & output.
                                   vectors_array,        &  !  11. output.
                                   values_array )           !  12. output.

  use blas_lapack_gbl, only: blasint
  use omp_lib

  implicit none

!====================================================================!
!                                                                    !
! 1. standard symmetric eigenvalue problem.                          !
!                                                                    !
!    a. we want to solve a*x = lambda*x in regular mode.             !
!                                                                    !
!    b. dsaupd = arpack reverse communication interface routine.     !
!                                                                    !
!    c. dseupd = arpack routine that returns ritz values             !
!                and ritz vectors.                                   !
!                                                                    !
!    e. leading dimensions for all arrays are :                      !
!                                                                    !
!       maxn   : maximum dimension of the matrix allowed.            !
!       maxnev : maximum "nev" allowed.                              !         
!       maxncv : maximum "ncv" allowed.                              !
!                                                                    !
!====================================================================!

! 2. input variables.

  integer(kind=shortint), intent(in), allocatable, dimension(:) :: row_index_array        !  1.
  integer(kind=shortint), intent(in), allocatable, dimension(:) :: column_index_array     !  2.
  real(kind=wp)   , intent(in), allocatable, dimension(:) :: ham_elem_array               !  3.
  character(len=2), intent(in) :: which                                                   !  4.
  integer         , intent(in) :: ncv_coeff                                               !  5.

! 3. input - output variables.
  
  integer               , intent(in) :: matrix_dimension                                  !  6.
  integer(kind=longint) , intent(in) :: non_zero_ham_elem                                 !  7.
  integer               , intent(inout) :: num_eigvals                                    !  8.
  integer               , intent(inout) :: maxitr                                         !  9.
  real(kind=wp), intent(inout) :: tol                                                     !  10.
  real(kind=wp), intent(out), allocatable, dimension(:,:) :: vectors_array                !  11.
  real(kind=wp), intent(out), allocatable, dimension(:,:) :: values_array                 !  12.

! 4. local variables. part 1/2.
  
  real(kind=wp), allocatable, dimension(:) :: x
  real(kind=wp), allocatable, dimension(:) :: y
  integer :: nt

!==============================================================================
! openmp_auto help variables.

  real(kind=wp), allocatable, dimension(:,:) ::  x1
  real(kind=wp), allocatable, dimension(:,:) ::  y1
  integer(kind=longint)  :: i2

!==============================================================================

  real(kind=wp), allocatable, dimension(:) :: workl 
  real(kind=wp), allocatable, dimension(:) :: workd 
  real(kind=wp), allocatable, dimension(:) :: resid
  real(kind=wp), allocatable, dimension(:) :: ax
  real(kind=wp) :: sigma

  logical(kind=blasint), allocatable, dimension(:) :: select_arpack
  logical(kind=blasint) :: rvec
  
  character(len=1) :: bmat
 
  integer(kind=blasint), dimension(11) :: iparam
  integer(kind=blasint), dimension(11) :: ipntr
  integer(kind=blasint) :: ido
  integer(kind=blasint) :: ncv
  integer(kind=blasint) :: nev
  integer(kind=blasint) :: lworkl
  integer(kind=blasint) :: info
  integer(kind=blasint) :: ierr
  integer(kind=blasint) :: mode
  integer(kind=blasint) :: ishfts
  integer(kind=blasint) :: maxnev
  integer(kind=blasint) :: maxncv
  integer(kind=blasint) :: matrix_dimen
  
! 5. local variables. part 2/2.
 
  real(kind=wp) :: t1                 ! timing variable.
  real(kind=wp) :: t2                 ! timing variable.
  real(kind=wp) :: total_time_mv      ! timing variable.
  real(kind=wp) :: total_time_dsaupd  ! timing variable.
  real(kind=wp) :: total_time_dseupd  ! timing variable.
  integer(kind=shortint) :: i1
  integer(kind=shortint), allocatable ,dimension(:) :: irow, icolumn

!$OMP PARALLEL

  nt = omp_get_num_threads()

!$OMP END PARALLEL


! 6. executable statements.
!    a standard eigenvalue problem is solved, bmat = 'i'.
!    "nev" is the number of eigenvalues to be approximated.
!    the following conditions must be satisfied:
!
!    nev <= maxnev,
!    nev + 1 <= ncv <= maxncv 

!========================================================================================!
!
! pgg note:
!  
! The following is a very important setting.
!
! 1. If the "ncv" is small (a typical small value is ncv = 2 * nev + 1) 
! then we get back less matrix-vectors operations
! and more restarted arnoldi operations.  
!
! 2. If the "ncv" is large (a tyical large value is ncv = 10 * nev +1)
! then the situation is the oposite tahn the above one.
! 
! 3. As a conclusions adjusting this value maybe we mabnaged to speed up thearpcak.
! 

  nev = num_eigvals
  ncv = ncv_coeff * nev + 1    ! the num of lanczos' vectors to be generated and used
                               ! at each restarted arnoldi algorithm.
!
!========================================================================================!

  matrix_dimen = matrix_dimension

  if ( ncv >= matrix_dimen ) then  ! adjusting the ncv.
  
     ncv = matrix_dimen      
            
  end if

  maxncv = ncv        ! this is my setting.
  maxnev = nev + 10   ! this is my setting.  

  allocate( x(1:matrix_dimen) )
  allocate( y(1:matrix_dimen) )  

!========================================================================================!
! Allocating ram space for the par_auto help variables.

  allocate( x1( 1:matrix_dimen, 1:nt ) )      ! x1
  allocate( y1( 1:matrix_dimen, 1:nt ) )      ! y1
  allocate( irow(1:nt) )                      ! irow
  allocate( icolumn(1:nt) )                   ! icolumn

!========================================================================================!

  allocate( vectors_array(matrix_dimen,maxncv) )
  allocate( workl(maxncv*(maxncv+8)) )
  allocate( workd(3*matrix_dimen) )
  allocate( values_array(maxncv,2) )
  allocate( resid(matrix_dimen) )
  allocate( ax(matrix_dimen) )  
  allocate( select_arpack(maxncv) )     

  ! 7. output to console 1.
    
  if ( nev .gt. maxnev ) then
    print *, ' error with _sdrv1: nev is greater than maxnev '
    go to 9000
  else if ( ncv .gt. maxncv ) then
    print *, ' error with _sdrv1: ncv is greater than maxncv '
    go to 9000
  end if
  
  bmat  = 'I'  ! normal problem.

! 8. setting the working space arrays.
!
!    the work array workl is used in dsaupd as
!    workspace.  its dimension lworkl is set as
!    illustrated below.  the parameter tol determines
!    the stopping criterion.  if tol<=0, machine
!    precision is used.  the variable ido is used for
!    reverse communication and is initially set to 0.
!    setting info=0 indicates that a random vector is
!    generated in dsaupd to start the arnoldi
!    iteration.
   
  lworkl = ncv*(ncv+8) ! workspace for dsaupd. a typical value here.
  
  info   = 0           ! random vector in dsaupd for a start of arnoldi.
  
  ido    = 0           ! must value for a start.

! 9. setting the remaining parameters.
!
!    this program uses exact shifts with respect to 
!    the current hessenberg matrix, iparam(1) = 1.
!    iparam(3) specifies the maximum number of arnoldi
!    iterations allowed. mode 1 of dsaupd is used,
!    iparam(7) = 1. all these options may be
!    changed by the user. for details, see the
!    documentation in dsaupd.

  ishfts    = 1      ! exact shifts.
   
  mode      = 1      ! typical problem.

  iparam    = 0

  iparam(1) = ishfts ! exact shifts to the matrix.
   
  iparam(3) = maxitr ! maximum number of arnoldi iterations allowed.
    
  iparam(7) = mode   ! mode 1 of dsaupd.

! 10. the main loop. reverse communication.
 
  ! setting clock to zero.
 
  total_time_mv     = 0.0_wp
  total_time_dsaupd = 0.0_wp
  total_time_dseupd = 0.0_wp

 10  continue

  ! 10.a. repeatedly call the routine dsaupd and take  
  !       actions indicated by parameter ido until    
  !       either convergence is indicated or maxitr   
  !       has been exceeded.         


  call cpu_time(t1)  ! timing dsaupd.
  
  call dsaupd ( ido, bmat, matrix_dimen, which,  &
                nev, tol, resid,                 & 
                ncv, vectors_array,              &
                matrix_dimen, iparam, ipntr,     &
                workd, workl,                    &
                lworkl, info )

  call cpu_time(t2)  ! timing dsaupd.

  total_time_dsaupd = total_time_dsaupd + (t2-t1)  ! timing dsaupd.

  if (ido .eq. -1 .or. ido .eq. 1) then ! if --> alpha.
         
  ! 10.b. perform matrix vector multiplication y <--- op*x.
  !       x=workd(ipntr(1)) as the input,
  !       and return the result to y=workd(ipntr(2)). 
            
    call cpu_time(t1)  ! timing matrix-vector.
 

!=======================================================================================!
! Initializing openmp help variables.
  
    do i1 = 1, nt
      x1(:,i1) = workd(ipntr(1) : ipntr(1) + matrix_dimen - 1)            !  1.
    end do

    y1 = 0.0_wp

!=======================================================================================!

! omp_auto do-loop.

!$OMP PARALLEL DO                     &
!$OMP DEFAULT( NONE )                 &
!$OMP PRIVATE( i1 )                   & 
!$OMP PRIVATE( i2 )                   &
!$OMP SHARED( non_zero_ham_elem )     &
!$OMP SHARED( nt )                    &
!$OMP SHARED( irow )                  &
!$OMP SHARED( icolumn )               &
!$OMP SHARED( y1 )                    &
!$OMP SHARED( x1 )                    &
!$OMP SHARED( ham_elem_array )        &
!$OMP SHARED( row_index_array )       &
!$OMP SHARED( column_index_array )    &
!$OMP NUM_THREADS( nt )

    do i1 = 1, nt

      do i2 = (i1-1) * non_zero_ham_elem / nt + 1 , (i1) * non_zero_ham_elem / nt     
 
        irow(i1) = row_index_array(i2)
 
        icolumn(i1) = column_index_array(i2)
 
        y1(irow(i1),i1) = y1(irow(i1),i1) + ham_elem_array(i2) * x1(icolumn(i1),i1)
     
        if (irow(i1).ne.icolumn(i1)) y1(icolumn(i1),i1) = y1(icolumn(i1),i1)+ham_elem_array(i2)*x1(irow(i1),i1)
        
      end do
 
    end do 

!$OMP END PARALLEL DO

    !====================================================================================!
    !# WARNING : DANGER OF ROUND-OFF ERRORS AT THE FOLLOWING STATEMENT.

    y = 0.0_wp

    do i1 = 1, nt

      y = y + y1(:,i1) 

    end do      
  
    !====================================================================================!
  
    workd(ipntr(2):2*matrix_dimen) = y

    call cpu_time(t2)  ! timing matrix-vector.

    total_time_mv = total_time_mv + (t2-t1) ! timing matrix-vector.

    go to 10 ! the most beautiful goto statement.

  end if ! end if --> alpha.


  ! 10.c. either we have convergence or there is an error.

  if ( info .lt. 0 ) then ! if --> beta  


  ! 10.d. error message. check the documentation in dsaupd.
        
    print *, ' '
    print *, ' error with _saupd, info = ', info
    print *, ' check documentation in _saupd '
    print *, ' '

  else ! if --> beta.


  ! 10.e. no fatal errors occurred. 
  !       post-process using dseupd.
  !       computed eigenvalues may be extracted.
  !       eigenvectors may also be computed now if desired,
  !       indicated by rvec = .true.   

    rvec = .true.

    call cpu_time(t1)  ! timing dseupd.

    call dseupd ( rvec, 'All', select_arpack,      &
                  values_array, vectors_array,     &
                  matrix_dimen, sigma, bmat,       &
                  matrix_dimen, which, nev,        &
                  tol, resid, ncv,                 &
                  vectors_array, matrix_dimen,     &
                  iparam, ipntr, workd, workl,     &
                  lworkl, ierr )

    call cpu_time(t2)  ! timing dseupd.
  
    total_time_dseupd = total_time_dseupd + (t2-t1) ! timing dseupd.

    write(*,*) "==========================================================================="
    write(*,*) " Inside subroutine: 'arpack_coo_openmp_auto_fortran'. "
    write(*,*) 
    write(*,*) " total_time_matrix_vector_products (seconds) --> ", total_time_mv
    write(*,*) " total_time_dsaupd_subroutine (seconds) --> ", total_time_dsaupd
    write(*,*) " total_time_dseupd --> ", total_time_dseupd
    write(*,*) 
    write(*,*) "==========================================================================="

  ! 10.f. eigenvalues are returned in the first column 
  !       of the two dimensional array "values_array" and the       
  !       corresponding eigenvectors are returned in   
  !       the first nev columns of the two dimensional 
  !       array "vectors_array" if requested.  otherwise, an         
  !       orthogonal basis for the invariant subspace  
  !       corresponding to the eigenvalues in "values_array" is     
  !       returned in "vectors_array".
                               
    if (ierr .ne. 0) then ! if --> gamma.

  ! 10.g. error condition:  
  !       check the documentation of dseupd. 

      print *, ' '
      print *, ' error with _seupd, info = ', ierr
      print *, ' check the documentation of _seupd. '
      print *, ' '
  
    end if ! end if --> gamma.

  ! 10.h. print additional convergence information.

    if (info .eq. 1) then ! if --> delta.
      print *, ' '
      print *, ' maximum number of iterations reached.'
      print *, ' '
    else if ( info .eq. 3) then ! if --> delta.
      print *, ' ' 
      print *, ' no shifts could be applied during implicit', &
               ' arnoldi update, try increasing ncv.'
      print *, ' '
    end if ! end if --> delta. 

    write(*,*) "==========================================================================="
    write(*,*) " Inside subroutine: 'arpack_coo_openmp_auto_fotran'. "
    write(*,*)
    write(*,*) " _sdrv1 "
    write(*,*) " ====== "
    write(*,*) 
    write(*,*) " size of the matrix is --> ", matrix_dimen
    write(*,*) " the number of ritz values requested is --> ", nev 
    write(*,*) " the number of arnoldi vectors generated (ncv) is --> ", ncv 
    write(*,*) " the portion of the spectrum is --> ", which
    write(*,*) " the number of converged ritz values is --> ", iparam(5) 
    write(*,*) " the number of implicit arnoldi update iterations taken is -->", iparam(3) 
    write(*,*) " the number of op*x is --> ", iparam(9)
    write(*,*) " the convergence criterion is --> ", tol 
    write(*,*)
    write(*,*) "==========================================================================="

  end if ! end if --> beta.
   
!  done with the subroutine.

 9000 continue

    num_eigvals = nev

  end subroutine arpack_coo_half_auto

! FINI.

 end module m_3_arpack_diagonalizers
