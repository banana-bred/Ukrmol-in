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
      module cdenprop_defs
      use precisn_gbl
      use blas_lapack_gbl, only: blasint
      implicit none
      !*********************************************************!
      !
      !Module containing parameter and derived type definitions
      !
      !*********************************************************!

      integer, parameter :: idp = wp
      integer, parameter :: ir_max=8 !Maximum no of IR per point group
      integer, parameter :: tsym_max=1000 !Maximum number of target symmetries.

      integer, parameter :: maxnuc=20
      integer, parameter :: maxprop=2  !Highest l to be read from target properties

      !Parameters used in formatting of prop.out for scattering calculations
      ! (I.E. when prop.out need to be read by swinterf)
      integer, parameter :: small_swintf_no_states = 100 !Number of states used in "small, swinterf-compatible"
      integer, parameter :: large_swintf_prop_states_nuclei = 1000000 !Sum of non-zero properties, number of states and number of nuclei. Used in "large swinterf-compatible"
      integer, parameter :: large_swintf_no_nuclei = 1000 !Number of nuclei used in "large swinterf-compatible"
      integer, parameter :: large_swintf_no_states = 10000 !Number of states used in "large swinterf-compatible"
      integer, parameter :: large_swintf_non_zero_properties = 1000000 !Number of non-zero properties used in "large swinterf-compatible"
!     ---------------------------

!     D2h group multiplication table
!     Note the change in dimensions from 8,8 to maxsym, maxsym
!     to achieve consistency via the global parameters.
      integer, dimension(ir_max,ir_max), parameter :: IPD2H=RESHAPE( (/ &
     &        1,2,3,4,5,6,7,8, &
     &        2,1,4,3,6,5,8,7, &
     &        3,4,1,2,7,8,5,6, &
     &        4,3,2,1,8,7,6,5, &
     &        5,6,7,8,1,2,3,4, &
     &        6,5,8,7,2,1,4,3, &
     &        7,8,5,6,3,4,1,2, &
     &        8,7,6,5,4,3,2,1/) , (/ ir_max, ir_max /) )


      character(len=5), dimension(8), parameter :: CD2H_symb=(/ &
     & ' AG  ', ' B3U ', ' B2U ', ' B1G ', ' B1U ', ' B2G ', ' B3G ', ' AU  ' /)
      character(len=5), dimension(4), parameter :: CC2V_symb=(/ &
     & ' A1  ', ' B1  ', ' B2  ', ' A2  '/)


      character(len=10), dimension(8), parameter ::spin_symb=(/&
     & 'Singlet   ', 'Doublet   ', 'Triplet   ', 'Quadruplet',      &
     & 'Quintuplet', 'Sextuplet ', 'Septuplet ', 'Octuplet  ' /)

!
!     Define the symbols for the elements in the periodic table
!
      character(len=2), dimension(103) :: ASYMB=(/&
     &     'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ',                   &
     &     'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', &
     &     'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', &
     &     'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', &
     &     'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', &
     &     'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', &
     &     'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', &
     &     'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', &
     &     'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', &
     &     'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', &
     &     'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'/)
!
!     Define the masses of the elements in the periodic table
!
      real(kind=idp), dimension(103) ::AMASS=(/&
     &      1.0078246_idp, 4.002601_idp, 7.01600_idp,                            &
     &      9.01218_idp, 11.009307_idp, 12.000000_idp, 14.0030738_idp,           &
     &     15.9949141_idp, 18.9984022_idp, 19.992441_idp, 22.9898_idp,           &
     &     23.98504_idp, 26.98153_idp, 27.976929_idp, 30.973764_idp,             &
     &     31.9720727_idp, 34.9688531_idp, 39.962386_idp, 38.96371_idp,          &
     &     39.96259_idp, 44.95592_idp, 48._idp, 50.9440_idp, 51.9405_idp,        &
     &     54.9380_idp, 55.9349_idp, 58.9332_idp, 57.9353_idp, 62.9296_idp,      &
     &     63.9291_idp, 68.9257_idp, 73.9219_idp, 74.9216_idp, 79.9165_idp,      &
     &     78.91839_idp, 83.91151_idp, 84.9117_idp, 87.9056_idp, 88.9059_idp,    &
     &     89.9043_idp, 92.9060_idp, 97.9055_idp, 98._idp, 101.9037_idp,         &
     &    102.9048_idp, 107.90389_idp, 106.90509_idp, 113.9036_idp, 114.9041_idp,&
     &    120._idp, 120.9038_idp, 129.9067_idp, 126.90466_idp, 131.90416_idp,    &
     &    132.9051_idp, 137.9050_idp, 138.9061_idp, 139.9053_idp, 140.9074_idp,  &
     &    141.9075_idp, 145._idp, 151.9195_idp, 152.9209_idp, 157.9241_idp,      &
     &    159.9250_idp, 163.9288_idp, 164.9303_idp, 165.9304_idp, 168.9344_idp,  &
     &    173.9390_idp, 174.9409_idp, 179.9468_idp, 180.9480_idp, 183.9510_idp,  &
     &    186.9560_idp, 192._idp, 192.9633_idp, 194.9648_idp, 196.9666_idp,      &
     &    201.970625_idp, 204.9745_idp, 207.9766_idp, 208.9804_idp, 209._idp,    &
     &    210._idp, 222._idp, 223._idp, 226._idp, 227._idp, 232._idp, 231._idp,  &
     &    238._idp,237._idp, 244._idp, 243._idp, 247._idp, 247._idp, 251._idp,   &
     &    252._idp, 257._idp, 258._idp, 259._idp, 260._idp/)

      type CSFheader
         CHARACTER(LEN=120) :: NAME
         integer :: MGVN
         double precision :: S, SZ, R, PIN
         integer :: NORB, NSRB,NOCSF,NELT, lcdof, IDIAG,NSYM, SYMTYP, lndof
         integer, dimension (:) :: npflg(6)
         double precision ::  thres
         integer :: NCTARG,NTGSYM
         integer, allocatable, dimension (:) :: iphz, nctgt, notgt, mcont, gucont, nob, ndtrf, nodo, numtgt
         integer, allocatable, dimension (:) :: itarget_overall_phase, idtarg
         integer, allocatable, dimension (:,:) ::  itarget_symmetry_order
         integer :: iposit
         integer, allocatable, dimension (:) :: nob0, nobl, nob0l,no_l2_virtuals
         integer :: l2nocsf,last_continuum_csf
      contains
         procedure :: dealloc=>dealloc_CSFheader
      end type CSFheader

      type CSFbody
         ! MPI shared memory windows
         integer :: icdo_window = -1
         integer :: indo_window = -1
         integer :: ndo_window  = -1
         integer :: cdo_window  = -1
         ! pointers to MPI shared memory
         integer,   pointer :: icdo(:) => null()
         integer,   pointer :: indo(:) => null()
         integer,   pointer :: ndo(:)  => null()
         real(idp), pointer :: cdo(:)  => null()
      contains
         procedure :: dealloc => dealloc_CSFbody
         final :: finalize_CSFbody
      end type CSFbody

      type CIvect
         integer :: nset, nrec
         CHARACTER(len=120) :: name
         integer :: nnuc,nocsf,nstat,mgvn
         real(kind=idp) :: S, SZ
         integer :: nelt
         real(kind=idp) :: e0
         character(len=8), dimension(:) :: cname(maxnuc)
         real(kind=idp), dimension(:) :: charge(maxnuc),xnuc(maxnuc), ynuc(maxnuc), znuc(maxnuc)
         real(kind=idp), allocatable, dimension(:) ::  ei
         integer , allocatable, dimension (:) :: iphz
         real(kind=idp), allocatable, dimension(:,:) :: CV
      !!!DESCRIPTORS FOR THE DISTRIBUTED SCALAPACK MATRIX, CV (if used)
         integer(blasint) :: mat_dimen,mat_dimen_r,mat_dimen_c
         integer(blasint) :: blacs_context
         integer(blasint) :: myrow,mycol
         integer(blasint) :: nprow,npcol
         integer(blasint) :: scal_block_size
         integer(blasint) :: local_row_dimen,local_col_dimen
         integer(blasint) :: descr_CV_mat(50)
         integer(blasint) :: lda
         logical :: CV_is_scalapack = .false.
      !!!DESCRIPTORS FOR THE DISTRIBUTED SCALAPACK MATRIX, CV (if used)
      contains
         procedure :: dealloc=>dealloc_CIvect
!         final :: destroy_CIvect !Use when intel finally fully implement 2003 standard
!         procedure :: destroy => destroy_CIvect
         procedure :: init_CV
         procedure :: final_CV
         procedure :: A_B_matmul
         procedure :: redistribute
         procedure :: gather_vectors
         procedure :: set_CV_element
         procedure :: add_to_CV_element
         procedure :: local_to_global
         procedure :: global_to_local

         procedure :: formatted_write
         procedure :: unformatted_write
         procedure :: stream_write

         ! generic :: write(formatted) => formatted_write
         ! generic :: write(unformatted) => unformatted_write
      end type CIvect

      type property_integrals
         integer :: no_of_integrals,no_of_properties
         integer, allocatable, dimension(:) :: lp, mp, nilmq, qp
         character(len=8), allocatable, dimension(:) :: property_name
         integer, allocatable, dimension(:) :: nob,mob, mpob
         integer, allocatable, dimension(:,:) :: indexv,inverted_indexv,istart
         real(kind=idp), allocatable, dimension(:,:) :: xintegrals
         logical :: ukrmolp_ints = .false.
      end type property_integrals

#if defined(scalapack) && defined(usempi)
      ! explicit interfaces for the ScaLAPACK routines used in this module
      interface
         subroutine blacs_get (context, what, val)
            import blasint
            integer(blasint), intent(in)  :: context, what
            integer(blasint), intent(out) :: val
         end subroutine blacs_get

         subroutine blacs_gridinit (context, layout, nprows, npcols)
            import blasint
            integer(blasint), intent(out) :: context
            integer(blasint), intent(in)  :: nprows, npcols
            character(len=1), intent(in)  :: layout
         end subroutine blacs_gridinit

         subroutine blacs_gridexit (context)
            import blasint
            integer(blasint), intent(in) :: context
         end subroutine blacs_gridexit

         subroutine blacs_gridinfo (context, nprows, npcols, myrow, mycol)
            import blasint
            integer(blasint), intent(in)  :: context
            integer(blasint), intent(out) :: nprows, npcols, myrow, mycol
         end subroutine blacs_gridinfo

         subroutine descinit (descr, mat_dimen_r, mat_dimen_c, block_size_r, block_size_c, first_r, first_c, context, lda, info)
            import blasint
            integer(blasint), intent(out) :: descr(*), info
            integer(blasint), intent(in)  :: mat_dimen_r, mat_dimen_c, block_size_r, block_size_c, first_r, first_c, context, lda
         end subroutine descinit

         integer(blasint) function numroc (mat_dimen, block_size, myproc, first_proc, nprocs)
            import blasint
            integer(blasint), intent(in) :: mat_dimen, block_size, myproc, first_proc, nprocs
         end function numroc

         integer(blasint) function indxl2g (indxloc, nb, iproc, isrcproc, nprocs)
            import blasint
            integer(blasint), intent(in) :: indxloc, nb, iproc, isrcproc, nprocs
         end function indxl2g

         subroutine infog1l (gindx, nb, nprocs, myroc, isrcproc, lindx, rocsrc)
            import blasint
            integer(blasint), intent(in)  :: gindx, nb, nprocs, myroc, isrcproc
            integer(blasint), intent(out) :: lindx, rocsrc
         end subroutine infog1l

         subroutine infog2l (grindx, gcindx, desc, nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc)
            import blasint
            integer(blasint), intent(in)  :: grindx, gcindx, desc(*), nprow, npcol, myrow, mycol
            integer(blasint), intent(out) :: lrindx, lcindx, rsrc, csrc
         end subroutine infog2l

         subroutine dgsum2d (icontxt, scope, top, m, n, a, lda, rdest, cdest)
            import blasint, wp
            integer(blasint), intent(in)    :: icontxt, m, n, lda, rdest, cdest
            character(len=1), intent(in)    :: scope, top
            real(wp),         intent(inout) :: a(*)
         end subroutine dgsum2d

         subroutine pdgemm (transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
            import blasint, wp
            character(len=1), intent(in)    :: transa, transb
            integer(blasint), intent(in)    :: m, n, k, ia, ja, ib, jb, ic, jc, desca(*), descb(*), descc(*)
            real(wp),         intent(in)    :: alpha, beta, a(*), b(*)
            real(wp),         intent(inout) :: c(*)
         end subroutine pdgemm

         subroutine pdgemr2d (m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)
            import blasint, wp
            integer(blasint), intent(in)    :: m, n, ia, ja, ib, jb, ictxt, desca(*), descb(*)
            real(wp),         intent(in)    :: a(*)
            real(wp),         intent(inout) :: b(*)
         end subroutine pdgemr2d
      end interface
#endif

      contains

!     Destructors for allocatable arrays in type definitions
!     ------------------------------------------------------------------
      subroutine dealloc_CIvect(this)
      implicit none

      class(CIvect), intent(inout) :: this

      if (allocated(this%ei)) deallocate(this%ei)
      if (allocated(this%iphz)) deallocate(this%iphz)
      if (allocated(this%CV)) deallocate(this%CV)

      end subroutine dealloc_CIvect

      subroutine dealloc_CSFheader(this)
      implicit none

      class(CSFheader), intent(inout) :: this

      if (allocated(this%iphz)) deallocate(this%iphz)
      if (allocated(this%nctgt)) deallocate(this%nctgt)
      if (allocated(this%notgt)) deallocate(this%notgt)
      if (allocated(this%numtgt)) deallocate(this%numtgt)
      if (allocated(this%itarget_symmetry_order)) deallocate(this%itarget_symmetry_order)
      if (allocated(this%idtarg)) deallocate(this%idtarg)
      if (allocated(this%mcont)) deallocate(this%mcont)
      if (allocated(this%gucont)) deallocate(this%gucont)
      if (allocated(this%nob)) deallocate(this%nob)
      if (allocated(this%ndtrf)) deallocate(this%ndtrf)
      if (allocated(this%nodo)) deallocate(this%nodo)
      if (allocated(this%nob0)) deallocate(this%nob0)
      if (allocated(this%nob0l)) deallocate(this%nob0l)
      if (allocated(this%nobl)) deallocate(this%nobl)
      if (allocated(this%no_l2_virtuals)) deallocate(this%no_l2_virtuals)

      end subroutine dealloc_CSFheader


    !> \brief  Memory cleanup for CSFbody
    !> \author J Benda
    !> \date   2010
    !>
    !> Performs automatic deallocation of members (if not already done).
    !>
    subroutine dealloc_CSFbody (this)

        use mpi_gbl, only: shared_communicator
        use mpi_memory_gbl, only: mpi_memory_deallocate_integer, mpi_memory_deallocate_real

        class(CSFbody), intent(inout) :: this

        if (associated(this % icdo)) then
            call mpi_memory_deallocate_integer(this % icdo, size(this % icdo), this % icdo_window, shared_communicator)
            nullify (this % icdo)
        end if

        if (associated(this % indo)) then
            call mpi_memory_deallocate_integer(this % indo, size(this % indo), this % indo_window, shared_communicator)
            nullify (this % indo)
        end if

        if (associated(this % ndo)) then
            call mpi_memory_deallocate_integer(this % ndo, size(this % ndo), this % ndo_window, shared_communicator)
            nullify (this % ndo)
        end if

        if (associated(this % cdo)) then
            call mpi_memory_deallocate_real(this % cdo, size(this % cdo), this % cdo_window, shared_communicator)
            nullify (this % cdo)
        end if

    end subroutine dealloc_CSFbody


    !> \brief  Class destructor for CSFbody
    !> \author J Benda
    !> \date   2020
    !>
    !> Performs automatic deallocation of members (if not already done).
    !>
    subroutine finalize_CSFbody (this)

        type(CSFbody), intent(inout) :: this

        call this % dealloc

    end subroutine finalize_CSFbody


    !> \brief   Get preferred maximal ScaLAPACK block size
    !> \author  J Benda
    !> \date    2022
    !>
    !> The default block size is 64. However, it can be overriden by the used by means of the environment
    !> variable CDENPROP_SCALAPACK_BLOCK_SIZE.
    !>
    integer function get_scalapack_max_block_size () result (block_size)

        use const_gbl, only: stdout

        character(len=1024) :: text
        integer             :: length, status, number

        block_size = 64

        call get_environment_variable('CDENPROP_SCALAPACK_BLOCK_SIZE', text, length, status)

        if (status == 0) then
            read (text, *) number
            if (number < 1) then
                write (stdout, '(3A)') 'Warning: Ignoring invalid CDENPROP_SCALAPACK_BLOCK_SIZE with value "', text, '"'
            else
                block_size = number
            end if
        end if

    end function get_scalapack_max_block_size


      subroutine init_CV(this,mat_dimen_r,mat_dimen_c)
      use mpi_gbl, only: nprocs, mpi_xermsg
      use const_gbl, only: stdout
      implicit none
      class(CIvect), intent(inout) :: this
      integer, intent(in) :: mat_dimen_r, mat_dimen_c

      integer :: ido, ifail
      integer(blasint) :: info, scalapack_min_block_size, scalapack_max_block_size

         write(stdout,'(/," CIvect%init_CV start")')
#if defined(scalapack) && defined(usempi)
         if (this%CV_is_scalapack) then !SCALAPACK ARRAY

            this%mat_dimen_r = mat_dimen_r
            this%mat_dimen_c = mat_dimen_c

            this%mat_dimen = 0
            if (this%mat_dimen_r .eq. mat_dimen_c) this%mat_dimen = this%mat_dimen_r

!We assume that BLACS has been initialized like this somewhere else.
!            do ido=1,int( sqrt( real(nprocs) ) + 1 )
!              if(mod(nprocs,ido) .eq. 0) this%nprow = ido
!            end do
!
!            this%npcol = nprocs/this%nprow
!
!            call blacs_get( -1, 0, this%blacs_context )
!            call blacs_gridinit(  this%blacs_context , 'r', this%nprow, this%npcol )
            call blacs_gridinfo(  this%blacs_context , this%nprow, this%npcol, this%myrow, this%mycol )

            write(stdout,"('context = ',i4,' nprocs = ',i4,' matdimen = ',2i8,&
                &' nrow = ',i8,' ncol = ',i8,' myrow = ',i8,' mycol = ',i8)") &
                this % blacs_context, nprocs, this % mat_dimen_r, this % mat_dimen_c, this % nprow, &
                this % npcol, this % myrow, this % mycol

            scalapack_min_block_size = 1
            scalapack_max_block_size = get_scalapack_max_block_size()

            this%scal_block_size = min ( this%mat_dimen_r/this%nprow, this%mat_dimen_c/this%npcol )
            this%scal_block_size = min(this%scal_block_size, scalapack_max_block_size)
            this%scal_block_size = max(this%scal_block_size, scalapack_min_block_size)

            this%local_row_dimen = numroc(this%mat_dimen_r,this%scal_block_size,this%myrow,0_blasint,this%nprow)
            this%local_col_dimen = numroc(this%mat_dimen_c,this%scal_block_size,this%mycol,0_blasint,this%npcol)

            this%lda = max (1_blasint,this%local_row_dimen)

            write(stdout,"('block_size = ',i4,' local_row_size = ',i8,' local_col_size = ',i8,' lda = ',i8)") &
                this % scal_block_size, this % local_row_dimen, this % local_col_dimen, this % lda

            if (this % myrow >= 0 .and. this % mycol >= 0) then
               call descinit (this % descr_CV_mat, this % mat_dimen_r, this % mat_dimen_c, this % scal_block_size, &
                              this % scal_block_size, 0_blasint, 0_blasint, this % blacs_context, this % lda, info)

               if (info /= 0) then
                  call mpi_xermsg('CIvect', 'init_CV', 'Error in getting description for A', int(info), 1)
               end if
            else
               this % descr_CV_mat(:) = -1
            end if

            if(allocated(this%CV)) deallocate(this%CV)
            allocate(this%CV(this%lda,this%local_col_dimen),stat=ifail)
            if (ifail /= 0) then
               call mpi_xermsg('CIvect', 'init_CV', 'Error in LOCAL this%CV allocation', ifail, 1)
            end if

            this%CV = 0.0_idp

         else !STANDARD ARRAY
#endif
            this%CV_is_scalapack = .false.

            this%mat_dimen_r = mat_dimen_r
            this%mat_dimen_c = mat_dimen_c

            this % local_row_dimen = mat_dimen_r
            this % local_col_dimen = mat_dimen_c

            this%mat_dimen = 0
            if (this%mat_dimen_r .eq. mat_dimen_c) this%mat_dimen = this%mat_dimen_r

            if(allocated(this%CV)) deallocate(this%CV)
            allocate(this%CV(mat_dimen_r,mat_dimen_c),stat=ifail)
            if (ifail /= 0) then
               call mpi_xermsg('CIvect', 'init_CV', 'Error in this%CV allocation', ifail, 1)
            end if

            this%CV = 0.0_idp
#if defined(scalapack) && defined(usempi)
         endif
#endif
         write(stdout,'(/," CIvect%init_CV finished")')

      end subroutine init_CV


      subroutine final_CV(this)
      use const_gbl, only: stdout
      implicit none
      class(CIvect), intent(inout) :: this

         if (allocated(this%CV)) deallocate(this%CV)

         this%mat_dimen = 0
         this%mat_dimen_r = 0
         this%mat_dimen_c = 0
         this%local_row_dimen = 0
         this%local_col_dimen = 0
         this%lda = 0

      end subroutine final_CV

      !>  this%CV = matmul(A%CV,B%CV)
      subroutine A_B_matmul(this,A,B,TRANSA,TRANSB)
      use maths, only: maths_dmatrix_multiply_blas95
      class(CIvect) :: this
      class(CIvect), intent(in) :: A, B
      character(len=1), intent(in) :: TRANSA,TRANSB
      logical :: ta, tb
      real(kind=idp), parameter :: alpha = 1.0_idp, beta = 0.0_idp
      integer(blasint), parameter :: one = 1
      integer(blasint) :: m, n, k
#if defined(scalapack) && defined(usempi)
         if (this%CV_is_scalapack) then !SCALAPACK ARRAY
            if (A%CV_is_scalapack .and. B%CV_is_scalapack) then
               m = this%mat_dimen_r
               n = this%mat_dimen_c
               if ((TRANSA .eq. 'T') .or. (TRANSA .eq. 't')) then
                  k = A%mat_dimen_r
               else
                  k = A%mat_dimen_c
               endif
               call pdgemm (transa, transb, m, n, k, alpha, A % CV, one, one, A % descr_CV_mat(1:9), B % CV, one, one, &
                            B % descr_CV_mat(1:9), beta, this % CV, one, one, this % descr_CV_mat(1:9))
            else
               print *,'A_B_matmul: at least one of A,B are not in SCALAPACK format!!!'
               stop
            endif
         else !STANDARD ARRAY
#endif
            ta = (TRANSA == 'T' .or. TRANSA == 't')
            tb = (TRANSB == 'T' .or. TRANSB == 't')
            call maths_dmatrix_multiply_blas95(A % CV, B % CV, this % CV, ta, tb)
#if defined(scalapack) && defined(usempi)
         endif
#endif

      end subroutine A_B_matmul

      subroutine gather_vectors(this,vec,nreq,rdest,cdest)
      implicit none
      class(CIvect), intent(in) :: this
      real(kind=idp), allocatable :: vec(:,:)
      integer, intent(in) :: nreq,rdest,cdest

      integer(blasint) :: err, ig, jg, il, jl, iprow, ipcol, irow, icol, ncols
#if defined(scalapack) && defined(usempi)
         if (this%CV_is_scalapack) then !SCALAPACK ARRAY

            call blacs_barrier(this%blacs_context, 'a')

            if (nreq > this%mat_dimen_c) then
               print *,'nreq outside of allowed range',nreq,this%mat_dimen_c
               stop
            endif

            if (allocated(vec)) deallocate(vec)
            allocate(vec(this%mat_dimen_r,nreq),stat=err)
            if (err .ne. 0) then
               print *,'allocation error in gather_vectors'
               stop
            endif

            vec = 0.0_idp
            do jg=1,nreq
               do ig=1, this%mat_dimen_r
                  call infog2l(ig, jg, this%descr_CV_mat(1:9), this%nprow, this%npcol, this%myrow, this%mycol, il, jl, iprow, ipcol)
                  if (this%myrow==iprow .and. this%mycol==ipcol) then
                     vec(ig,jg) = this%CV(il,jl)
                  endif
               enddo
            enddo

            !Gather the vectors on process with coordinates (rdest, cdest) if (rdest=-1, cdest=-1) then the vector is gathered on all processes
            irow = rdest
            icol = cdest
            ncols = nreq
            call dgsum2d(this%blacs_context, 'all', ' ', this%mat_dimen_r, ncols, vec, this%mat_dimen_r, irow, icol)

         else !STANDARD ARRAY
#endif
            if (allocated(vec)) deallocate(vec)
            allocate(vec(size(this%CV,1),nreq),stat=err)
            if (err .ne. 0) then
               print *,'allocation error in gather_vectors standard',size(this%CV,1),nreq
               stop
            endif

            vec(:,1:nreq) = this%CV(:,1:nreq)
#if defined(scalapack) && defined(usempi)
         endif
#endif
      end subroutine gather_vectors


        !> \brief   Redistribute matrix between two BLACS contexts
        !> \authors J Benda
        !> \date    2019 - 2022
        !>
        !> Redistribute matrix, present in one BLACS context, to another BLACS context. The source and target
        !> matrix sizes do not have to be equal. In that case, only the largest common submatrix will be redistributed.
        !>
        !> \warning Typical legacy ScaLAPACK implementations use 4-byte integers internally. When redistributing
        !>          a large matrix using PDGEMR2D, integer overflow can occur, which often manifests as a memory allocation
        !>          failure.
        !>
        !> \param[inout] this         Distributed matrix to write.
        !> \param[in]    that         Distributed matrix to read.
        !> \param[in]    context_opt  BLACS context that contains both all source and all destination MPI processes.
        !>                            If not provided, the subroutine will use this % blacs_context instead.
        !> \param[in]    nrows_opt    Number of rows to copy. If not provided, all rows will be copied.
        !> \param[in]    ncols_opt    Number of columns to copy. If not provided, all columns will be copied.
        !> \param[in]    arow_opt     Initial source row index. Default: 1.
        !> \param[in]    acol_opt     Initial source column index. Default: 1.
        !> \param[in]    brow_opt     Initial destination row index. Default: 1.
        !> \param[in]    bcol_opt     Initial destination column index. Default: 1.
        !>
        subroutine redistribute (this, that, context_opt, nrows_opt, ncols_opt, arow_opt, acol_opt, brow_opt, bcol_opt)

            class(CIvect),    intent(inout)        :: this
            type(CIvect),     intent(in)           :: that
            integer(blasint), intent(in), optional :: context_opt
            integer,          intent(in), optional :: nrows_opt, ncols_opt, arow_opt, acol_opt, brow_opt, bcol_opt

            integer(blasint) :: context, m, n, ai, aj, bi, bj

            m = min(this % mat_dimen_r, that % mat_dimen_r)
            n = min(this % mat_dimen_c, that % mat_dimen_c)

            if (present(nrows_opt)) m = nrows_opt
            if (present(ncols_opt)) n = ncols_opt

            ai = 1;  if (present(arow_opt)) ai = arow_opt
            aj = 1;  if (present(acol_opt)) aj = acol_opt
            bi = 1;  if (present(brow_opt)) bi = brow_opt
            bj = 1;  if (present(bcol_opt)) bj = bcol_opt

#if defined(scalapack) && defined(usempi)
            if (this % CV_is_scalapack) then
                if (present(context_opt)) then
                    context = context_opt
                else
                    context = this % blacs_context
                end if

                call pdgemr2d(m, n, that % CV, ai, aj, that % descr_CV_mat(1:9), &
                                    this % CV, bi, bj, this % descr_CV_mat(1:9), context)
            else
#endif
                this % CV(ai : ai + m - 1, aj : aj + n - 1) = that % CV(bi : bi + m - 1, bj : bj + n - 1)
#if defined(scalapack) && defined(usempi)
            end if
#endif
        end subroutine redistribute


        !> \brief   Convert local matrix indices to global ones
        !> \authors J Benda
        !> \date    2019
        !>
        !> In case of distributed matrix, get the real coordinates (in the global matrix).
        !> Trivially returns the same indices in non-distributed case.
        !>
        !> \param[in]  this   Distributed matrix to read.
        !> \param[in]  i_loc  Row index in the local portion of the matrix.
        !> \param[in]  j_loc  Column index in the local portion of the matrix.
        !> \param[out] i      Global row index.
        !> \param[out] j      Global column index.
        !>
        subroutine local_to_global (this, i_loc, j_loc, i, j)

            class(CIvect), intent(in)  :: this
            integer,       intent(in)  :: i_loc, j_loc
            integer,       intent(out) :: i, j

            if (i_loc < 1 .or. this % local_row_dimen < i_loc .or. &
                j_loc < 1 .or. this % local_col_dimen < j_loc) then
                i = -1
                j = -1
                return
            end if

#if defined(usempi) && defined(scalapack)
            if (this % CV_is_scalapack) then
                i = indxl2g(int(i_loc, blasint), this % scal_block_size, this % myrow, 0_blasint, this % nprow)
                j = indxl2g(int(j_loc, blasint), this % scal_block_size, this % mycol, 0_blasint, this % npcol)
            else
#endif
                i = i_loc
                j = j_loc
#if defined(usempi) && defined(scalapack)
            end if
#endif
        end subroutine local_to_global


        !> \brief   Convert global matrix indices to local ones
        !> \authors J Benda
        !> \date    2019
        !>
        !> In case of distributed matrix, get the local coordinates corresponding to
        !> the global ones. Trivially returns the same indices in non-distributed case.
        !> If the specified global matrix row/column is not located at the current processor,
        !> -1 is returned in the corresponding index. Only when both `i_loc` and
        !> `j_loc` are positive, the element specified by the `i`,`j` pair is located
        !> at the current processor.
        !>
        !> \param[in]  this   Distributed matrix to read.
        !> \param[in]  i      Global row index.
        !> \param[in]  j      Global column index.
        !> \param[out] i_loc  Row index in the local portion of the matrix.
        !> \param[out] j_loc  Column index in the local portion of the matrix.
        !>
        subroutine global_to_local (this, i, j, i_loc, j_loc)

            class(CIvect), intent(in)  :: this
            integer,       intent(in)  :: i, j
            integer,       intent(out) :: i_loc, j_loc
#if defined(usempi) && defined(scalapack)
            integer(blasint) :: proc_row, proc_col, zero = 0, gi, gj, li, lj
#endif
            if (i < 1 .or. this % mat_dimen_r < i .or. &
                j < 1 .or. this % mat_dimen_c < j) then
                i_loc = -1
                j_loc = -1
                return
            end if

#if defined(usempi) && defined(scalapack)
            if (this % CV_is_scalapack) then
                gi = i; call infog1l(gi, this % scal_block_size, this % nprow, this % myrow, zero, li, proc_row)
                gj = j; call infog1l(gj, this % scal_block_size, this % npcol, this % mycol, zero, lj, proc_col)
                i_loc = li; if (proc_row /= this % myrow) i_loc = -1
                j_loc = lj; if (proc_col /= this % mycol) j_loc = -1
            else
#endif
                i_loc = i
                j_loc = j
#if defined(usempi) && defined(scalapack)
            end if
#endif
        end subroutine global_to_local


        !> \brief  Write distributed matrix to a stream file
        !> \author J Benda
        !> \date   2023
        !>
        !> Write the (potentially distributed) CI matrix to a stream file identified by the MPI-IO handle `fh`. This subroutine
        !> use the standard MPI-IO routines. If the matrix is distributed (i.e. `CV_is_scalapack` is true), then all processes
        !> in the BLACS context of this matrix need to call this.
        !>
        !> The default high-performance MPI-IO module "Vulcan" of Open MPI begins with allocation of a descriptor data structure
        !> where the matrix to be written is first indexed by extracting the beginning and the size of every contiguous segment
        !> of data in memory of any process. In case of a distributed matrix with a given ScaLAPACK block size n, these contiguous
        !> segments are vectors of size n, into which each column of the matrix is decomposed. For an N-by-N matrix, there will be
        !> N*(N/n) such segments, each of size 2×8 bytes. For matrix of rank N = 150,000 and ScaLAPACK block of size n = 64,
        !> this would mean allocation of 5 GiB *by every process*! To avoid such an excessive (and non-scalable!) memory use,
        !> this subroutine first redistributes the matrix to be written in such a way that it is distributed among processes
        !> by whole contiguous columns. This reduces the memory need by the factor N/n. Other MPI distributions *may* not have this
        !> issue, but this reorganization also should not harm anything.
        !>
        !> \param[in] this  CIvect matrix to write.
        !> \param[in] fh    A valid MPI-IO file handle of a file opened for writing.
        !>
        subroutine stream_write (this, fh)

            use mpi_gbl, only: mpiint, mpi_mod_barrier, mpi_mod_file_write

            class(CIvect),   intent(in) :: this
            integer(mpiint), intent(in) :: fh

            logical, parameter :: recombine_columns_first = .true.  ! redistribute matrix by columns before starting MPI-IO

            type(CIvect)     :: m
            integer(blasint) :: info, zero = 0, one = 1
            integer          :: ierr

#if defined(scalapack) && defined(usempi)
            if (this % CV_is_scalapack) then
                if (recombine_columns_first) then
                    ! create BLACS context with one process row and many columns
                    call blacs_get(-one, zero, m % blacs_context)
                    call blacs_gridinit(m % blacs_context, 'r', one, this % nprow * this % npcol)
                    call blacs_gridinfo(m % blacs_context, m % nprow, m % npcol, m % myrow, m % mycol)

                    ! initialize the auxiliary distributed matrix m
                    m % mat_dimen_r     = this % mat_dimen_r
                    m % mat_dimen_c     = this % mat_dimen_c
                    m % lda             = this % mat_dimen_r
                    m % CV_is_scalapack = this % CV_is_scalapack
                    m % local_row_dimen = m % mat_dimen_r
                    call descinit(m % descr_CV_mat, m % mat_dimen_r, m % mat_dimen_c, m % mat_dimen_r, one, &
                                  zero, zero, m % blacs_context, m % mat_dimen_r, info)
                    m % local_col_dimen = numroc(m % mat_dimen_c, one, m % mycol, zero, m % npcol)
                    allocate (m % CV(m % lda, m % local_col_dimen))

                    ! copy the original distributed matrix elements to the auxiliary distributed matrix using ScaLAPACK subroutine
                    call m % redistribute(this)

                    ! release the BLACS context (no longer needed)
                    call blacs_gridexit(m % blacs_context)

                    ! write the auxiliary matrix to disk
                    call mpi_mod_file_write(fh, int(m % mat_dimen_r), int(m % mat_dimen_c), &
                                                int(m % nprow), int(m % npcol),&
                                                int(m % mat_dimen_r), 1, m % CV, &
                                                int(m % local_row_dimen), int(m % local_col_dimen))
                else
                    ! or just write the original matrix as it is (may blow up with Open MPI)
                    call mpi_mod_file_write(fh, int(this % mat_dimen_r), int(this % mat_dimen_c), &
                                                int(this % nprow), int(this % npcol),&
                                                int(this % scal_block_size), int(this % scal_block_size), this % CV, &
                                                int(this % local_row_dimen), int(this % local_col_dimen))
                end if
                call mpi_mod_barrier(ierr)
            else
#endif
                call mpi_mod_file_write(fh, this % CV, int(this % mat_dimen_r), int(this % mat_dimen_c))
#if defined(scalapack) && defined(usempi)
            end if
#endif

        end subroutine stream_write


        !> \brief   Write distributed matrix to file (formatted)
        !> \author  J Benda
        !> \date    2023
        !>
        !> Derived-type I/O overload for formatted files. Will write all columns of the (potentially) distributed matrix
        !> into a single record on the file specified by unit `lu`. Only the matrix elements are written to file, but none
        !> of the other information about its distribution.
        !>
        !> The format of the output is '(xDy.z)', where the integer number 'x', 'y' and 'z' are given in the derived type
        !> format specifier as
        !> ```
        !>    write (u, '(dt(x,y,z))') cdenprop_matrix
        !> ```
        !>
        !> \param[in]    this        Distributed matrix to write to file.
        !> \param[in]    lu          Unit number of unformatted output.
        !> \param[in]    iotype      I/O type (not used).
        !> \param[in]    v_list      I/O format specification.
        !> \param[out]   io_stat     Exit success indicator.
        !> \param[inout] io_message  Status message to be filled in on return.
        !>
        subroutine formatted_write (this, lu, io_type, v_list, io_stat, io_message)

            use iso_c_binding,   only: c_f_pointer, c_loc
            use iso_fortran_env, only: int32
            use mpi_gbl,         only: master, myrank, mpi_mod_bcast

            class(CIvect), intent(in)    :: this
            integer,       intent(in)    :: lu
            character(*),  intent(in)    :: io_type
            integer,       intent(in)    :: v_list(:)
            integer,       intent(out)   :: io_stat
            character(*),  intent(inout) :: io_message

            target :: v_list

            real(wp), allocatable :: column(:)
            integer(int32), pointer :: v_list32(:)

            character(len=20) :: number_format
            integer(blasint)  :: iprow, ipcol, i, j, k, l, remaining, nelem, numbers_per_line

            io_stat = 0
            remaining = 0
            nelem = this%mat_dimen_r
            numbers_per_line = 4
            number_format = '(4d20.13,/)'

            ! build the requested number format
            if (size(v_list) == 3) then
#ifdef __GFORTRAN__
                ! GCC bug #108680: gfortran actually gives us array of *short* integers regardless of -fdefault-integer-8
                ! (affected versions: at least 11.3, 12.2 and 13.0)
                call c_f_pointer(c_loc(v_list), v_list32, [ size(v_list) ])  ! reinterpret array as 32-bit integers
                numbers_per_line = v_list32(1)
                write (number_format, '(a,i0,a,i0,a,i0,a)') '(', v_list32(1), 'd', v_list32(2), '.', v_list32(3), ',/)'
#else
                numbers_per_line = v_list(1)
                write (number_format, '(a,i0,a,i0,a,i0,a)') '(', v_list(1), 'd', v_list(2), '.', v_list(3), ',/)'
#endif
            end if

            allocate (column(this%mat_dimen_r + numbers_per_line - 1))

            do l = 1, this%mat_dimen_c
                ! obtain the next column
#if defined(scalapack) && defined(usempi)
                if (this % CV_is_scalapack) then
                    ! copy out all local elements of this matrix column
                    do k = 1, this%mat_dimen_r
                        call infog2l(k, l, this%descr_CV_mat, this%nprow, this%npcol, &
                                     this%myrow, this%mycol, i, j, iprow, ipcol)
                        if (this%myrow == iprow .and. this%mycol == ipcol) then
                            column(remaining + k) = this%CV(i, j)
                        else
                            column(remaining + k) = 0
                        end if
                     end do
                     ! sum the partial columns from all tasks to master
                     call dgsum2d(this%blacs_context, 'all', ' ', this%mat_dimen_r, 1_blasint, &
                                  column(remaining + 1:remaining + this%mat_dimen_r), this%mat_dimen_r, 0_blasint, 0_blasint)
                else
#endif
                    column(remaining + 1:remaining + this%mat_dimen_r) = this%CV(1:this%mat_dimen_r, l)
#if defined(scalapack) && defined(usempi)
                end if
#endif

                do k = 1, nelem / numbers_per_line
                    ! master will write the vector to file, placing `numbers_per_line` on each line
                    if (myrank == master) then
                        write (lu, number_format, iostat=io_stat, iomsg=io_message) &
                            column((k - 1)*numbers_per_line + 1 : k*numbers_per_line)
                    end if

                    ! if there was any problem with writing, notify other tasks and exit
                    call mpi_mod_bcast(io_stat, master)
                    if (io_stat /= 0) return
                end do

                ! move non-written data to the beginning of the vector
                remaining = mod(nelem, numbers_per_line)
                if (remaining > 0) then
                    column(1 : remaining) = column(nelem - remaining + 1 : nelem)
                end if
                nelem = remaining + this%mat_dimen_r
            end do

            ! write the remaining numbers
            if (remaining > 0) then
                if (myrank == master) then
                    write (lu, number_format, iostat=io_stat, iomsg=io_message) column(1 : remaining)
                end if
                call mpi_mod_bcast(io_stat, master)
                if (io_stat /= 0) return
            end if

        end subroutine formatted_write


        !> \brief   Write distributed matrix to file (unformatted)
        !> \author  J Benda
        !> \date    2023
        !>
        !> Derived-type I/O overload for unformatted files. Will write all columns of the (potentially) distributed matrix
        !> into a single record on the file specified by unit `lu`. Only the matrix elements are written to file, but none
        !> of the other information about its distribution.
        !>
        !> \param[in]    this        Distributed matrix to write to file.
        !> \param[in]    lu          Unit number of unformatted output.
        !> \param[out]   io_stat     Exit success indicator.
        !> \param[inout] io_message  Status message to be filled in on return.
        !>
        subroutine unformatted_write (this, lu, io_stat, io_message)

            use mpi_gbl, only: master, myrank, mpi_mod_bcast

            class(CIvect), intent(in)    :: this
            integer,       intent(in)    :: lu
            integer,       intent(out)   :: io_stat
            character(*),  intent(inout) :: io_message

            real(wp), allocatable :: column(:)
            integer(blasint) :: i, j, mat_row, mat_col, iprow, ipcol

            io_stat = 0

            allocate (column(this%mat_dimen_r))

            do mat_col = 1, this%mat_dimen_c
                ! obtain the next column
#if defined(scalapack) && defined(usempi)
                if (this % CV_is_scalapack) then
                    ! copy out all local elements of this matrix column
                    do mat_row = 1, this%mat_dimen_r
                        call infog2l(mat_row, mat_col, this%descr_CV_mat, &
                                     this%nprow, this%npcol, this%myrow, this%mycol, i, j, iprow, ipcol)
                        if (this%myrow == iprow .and. this%mycol == ipcol) then
                           column(mat_row) = this%CV(i, j)
                        else
                           column(mat_row) = 0
                        end if
                     end do
                     ! sum the partial columns from all tasks to master
                     call dgsum2d(this%blacs_context, 'all', ' ', this%mat_dimen_r, 1_blasint, &
                                  column, this%mat_dimen_r, 0_blasint, 0_blasint)
                else
#endif
                    ! just pick the column of the local matrix
                    column = this%CV(:, mat_col)
#if defined(scalapack) && defined(usempi)
                end if
#endif

                ! master will write the vector to file
                if (myrank == master) then
                   write (lu, iostat=io_stat, iomsg=io_message) column
                end if

                ! if there was any problem with writing, notify other tasks and exit
                call mpi_mod_bcast(io_stat, master)
                if (io_stat /= 0) return
            end do

        end subroutine unformatted_write


      subroutine set_CV_element(this,val,i,j)
      implicit none
      class(CIvect), intent(inout) :: this
      real(kind=idp), intent(in) :: val
      integer, intent(in) :: i,j

      integer(blasint) :: i_loc, j_loc, proc_row, proc_col
#if defined(scalapack) && defined(usempi)
         if (this%CV_is_scalapack) then !SCALAPACK ARRAY
            !Figure out which proc it belongs to and the local matrix index
            call infog2l(int(i, blasint), int(j, blasint), this % descr_CV_mat(1:9), &
                         this % nprow, this % npcol, this % myrow, this % mycol, i_loc, j_loc, proc_row, proc_col)
            if ((this%myrow .eq. proc_row) .and. (this%mycol .eq. proc_col)) then
               if ((i_loc > this%local_row_dimen) .or. (j_loc > this%local_col_dimen)) then
                  print *,'error inserting element',i,j,i_loc,j_loc,this%local_row_dimen,this%local_col_dimen
               endif
               this%CV(i_loc, j_loc) = val
            endif
         else !STANDARD ARRAY
#endif
            this%CV(i,j) = val
#if defined(scalapack) && defined(usempi)
         endif
#endif
      end subroutine set_CV_element

      subroutine add_to_CV_element(this,val,i,j)
      implicit none
      class(CIvect), intent(inout) :: this
      real(kind=idp), intent(in) :: val
      integer, intent(in) :: i,j

      integer(blasint) :: i_loc, j_loc, proc_row, proc_col
#if defined(scalapack) && defined(usempi)
         if (this%CV_is_scalapack) then !SCALAPACK ARRAY
            !Figure out which proc it belongs to and the local matrix index
            call infog2l(int(i, blasint), int(j, blasint), this % descr_CV_mat(1:9), &
                         this % nprow, this % npcol, this % myrow, this % mycol, i_loc, j_loc, proc_row, proc_col)
            if ((this%myrow .eq. proc_row) .and. (this%mycol .eq. proc_col)) then
               if ((i_loc > this%local_row_dimen) .or. (j_loc > this%local_col_dimen)) then
                  print *,'error adding into element',i,j,i_loc,j_loc,this%local_row_dimen,this%local_col_dimen
               endif
               this%CV(i_loc, j_loc) = this%CV(i_loc, j_loc) + val
            endif
         else !STANDARD ARRAY
#endif
            this%CV(i,j) = this%CV(i,j) + val
#if defined(scalapack) && defined(usempi)
         endif
#endif
      end subroutine add_to_CV_element


      end module
