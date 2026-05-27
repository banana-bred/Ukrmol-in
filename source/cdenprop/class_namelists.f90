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
      module class_namelists
      use cdenprop_defs, only: tsym_max
      implicit none

      integer, parameter :: idp = kind(1d0)

      private

      type, public :: denprop_namelist
!        Namelist variables for &INPUT
         character(len=80) :: NAME
         integer :: NTGT,  NFTG, IWRITE, LNGTH, LENGTH, NFTD, NFTA, NFDA, LCOF, NFTMT, NMSET, NFTDL,&
     &              LDAR, NFTINT, iplotfg, iqdfg, ISW, IPOL, NFTMOM, NMSET2, iprop, idiag, grounden,&
     &              geocfn, ksym, gflag
         integer :: lutdip, verbosity  ! used only in CDENPROP_ALL
         integer, allocatable :: NFTSOR(:), NTGTF(:), NTGTS(:), NTGTL(:), NPFLG(:), numtgt(:)
         integer, allocatable :: itarget_spin(:)  ! used only in cdenprop_all
         logical :: zlast, ukrmolp_ints, QMOLN
!        arrays sizes
         integer,public :: maxtgt, maxpfg
      contains
         procedure, public :: dealloc=>dealloc_target_namelist
         procedure, public :: init=>init_target_namelist
         procedure, public :: read_namelist
         procedure, public :: write_namelist
      end type denprop_namelist

      type, public :: cdenprop_namelists
!~          namelist /DENINP/ nciset, lucsf, lucivec, luneutralci, numtgt, luprop, lupropw,          &
!~      &                  lutargci, nstat, lutdip, ntset, threshold, lunp1targ, no_l2_virtuals,      &
!~      &                  ludyson, ludyson_gcsf, itarget_spin, lutarg_phase, max_multipole,          &
!~      &                  iuse_bound, ibra_or_ket_bound, ukrmolp_ints,isw
!~
!~          namelist /L2CONTRACT/ nciset_l2contract, lucivec_l2contract,lutransmat_l2contract

         integer :: lucivec, lutdip, ntset, luneutralci, ludyson, ludyson_gcsf, lucivec_l2contract,&
                    lutarg_phase, max_multipole, iuse_bound, ibra_or_ket_bound, isw, luprop, lupropw,&
                    nbset, gflag
         character(len=11) :: bform
         real(kind=idp) :: threshold
         integer, allocatable, dimension(:) :: lucsf, nciset, nstat, numtgt, notgt, lunp1targ,     &
     &                                         no_l2_virtuals, nciset_l2contract,                  &
     &                                         lutransmat_l2contract, itarget_spin
         logical :: ukrmolp_ints

!        arrays sizes
         integer :: maxtgt
      contains
         procedure, public :: dealloc=>dealloc_cdenprop_namelists
         procedure, public :: init=>init_cdenprop_namelists
         procedure, public :: construct_namelists_from_denprop_input
         procedure, public :: return_namelist_variables

      end type cdenprop_namelists

      contains

!     *********************************
!     Subroutines for cdenrpo_namelists
!     *********************************

      subroutine init_cdenprop_namelists(this, maxtgt)
      implicit none

!     Arguments
      class(cdenprop_namelists), intent(inout) :: this
      integer :: maxtgt

      if (.not. allocated(this%lucsf ) ) allocate( this%lucsf(maxtgt)  )
      if (.not. allocated(this%nciset ) ) allocate( this%nciset(maxtgt)  )
      if (.not. allocated(this%nstat ) ) allocate( this%nstat(maxtgt)  )
      if (.not. allocated(this%numtgt ) ) allocate( this%numtgt(tsym_max)  )
      if (.not. allocated(this%notgt ) ) allocate( this%notgt(maxtgt)  )
      if (.not. allocated(this%lunp1targ ) ) allocate( this%lunp1targ(maxtgt)  )
      if (.not. allocated(this%no_l2_virtuals ) ) allocate( this%no_l2_virtuals(maxtgt)  )
      if (.not. allocated(this%nciset_l2contract ) ) allocate( this%nciset_l2contract(maxtgt)  )
      if (.not. allocated(this%lutransmat_l2contract ) ) allocate( this%lutransmat_l2contract(maxtgt)  )
      if (.not. allocated(this%itarget_spin ) ) allocate( this%itarget_spin(tsym_max)  )

!~       this%isw=0
!~       this%luprop=17                ! File unit for input of property integrals.
!~       this%ludyson=123              ! File unit for dyson orbital output.
!~       this%ludyson_gcsf=9978        ! Output information needed to contract the L2 configurations.
!~       this%lupropw=24               ! File unit for dipole output.
!~       this%lutarg_phase=333         ! File unit for scattering run target vectors.
!~       this%nstat=0                  ! Number of states to calculate dipoles for defaults to all (0,0).
!~       this%lutdip=24                ! File unit for bound state transition dipoles.
!~       this%ntset=1                  ! Set for bound state transtion dipoles.
!~       this%threshold =1.E-10_idp    ! Printing threshold for dipoles.
!~       this%no_l2_virtuals=0         ! Don't change.
!~       this%itarget_spin=2           ! Array of dim no_target_symmetries: defaults to doublet.
!~       this%max_multipole=1          ! 1: Dipole transitions to the continuum.
!~                                     ! 2: Quad transitions.
!~       this%iuse_bound=0             ! Use bound states from BOUND.
!~       this%ibra_or_ket_bound=0      ! Set which are the bound states - only needed for Dyson orbital
!~                                     ! output when iuse_bound=0 and both bra and ket contain continuum
!~                                     ! orbitals.
!~       this%ukrmolp_ints = .false.   ! Are we using property integrals generated using the new integral code?

!     Defaults for L2CONTRACT
      this%nciset_l2contract=0
      this%lucivec_l2contract=0
      this%lutransmat_l2contract(1)=3001
      this%lutransmat_l2contract(2)=3002

 !    Set parameter values
      this%maxtgt=maxtgt

      this%nbset = 1
      this%gflag = 1
      this%bform = 'UNFORMATTED'

      end subroutine init_cdenprop_namelists

      subroutine dealloc_cdenprop_namelists(this)
      implicit none

      class(cdenprop_namelists), intent(inout) :: this

      if (allocated(this%lucsf ) ) deallocate( this%lucsf  )
      if (allocated(this%nciset ) ) deallocate( this%nciset  )
      if (allocated(this%nstat ) ) deallocate( this%nstat  )
      if (allocated(this%numtgt ) ) deallocate( this%numtgt  )
      if (allocated(this%notgt ) ) deallocate( this%notgt  )
      if (allocated(this%lunp1targ ) ) deallocate( this%lunp1targ  )
      if (allocated(this%no_l2_virtuals ) ) deallocate( this%no_l2_virtuals  )
      if (allocated(this%nciset_l2contract ) ) deallocate( this%nciset_l2contract  )
      if (allocated(this%lutransmat_l2contract ) ) deallocate( this%lutransmat_l2contract  )
      if (allocated(this%itarget_spin ) ) deallocate( this%itarget_spin  )

      end subroutine dealloc_cdenprop_namelists

      subroutine construct_namelists_from_denprop_input(this,denprop_input, isym_one, isym_two)
      implicit none
!     Arguments
      class(cdenprop_namelists), intent(inout) :: this
      class(denprop_namelist), intent(inout) :: denprop_input
      integer :: isym_one, isym_two

      this%lucsf(1)=denprop_input%NFTSOR(isym_one)
      this%lucsf(2)=denprop_input%NFTSOR(isym_two)

      this%lucivec=denprop_input%NFTG

      this%nciset(1)= denprop_input%NTGTF(isym_one)
      this%nciset(2)= denprop_input%NTGTF(isym_two)

      this%nstat(1)= denprop_input%NTGTL(isym_one)
      this%nstat(2)= denprop_input%NTGTL(isym_two)

      this%luprop=denprop_input%NFTINT
      this%lupropw=denprop_input%NFTMT
      this%lutdip=denprop_input%lutdip

      this%numtgt=denprop_input%numtgt
      this%itarget_spin=denprop_input%itarget_spin
      this%max_multipole=2
      this%isw=denprop_input%isw

      this%ukrmolp_ints=denprop_input%ukrmolp_ints

      end subroutine construct_namelists_from_denprop_input

      subroutine return_namelist_variables (this, lucsf, lucivec, nciset, nstat, numtgt, &
                                            lutdip, luprop, lupropw, max_multipole, ukrmolp_ints, isw, gflag, itarget_spin)
      implicit none
!     Arguments
      class(cdenprop_namelists), intent(inout) :: this
      integer :: lucsf(:), nciset(:), nstat(:), numtgt(tsym_max), lucivec, lutdip, luprop, lupropw, max_multipole, isw, gflag
      integer :: itarget_spin(:)
      logical :: ukrmolp_ints

      lucsf(1)=this%lucsf(1)
      lucsf(2)=this%lucsf(2)

      nciset(1)=this%nciset(1)
      nciset(2)=this%nciset(2)

      nstat(1)=this%nstat(1)
      nstat(2)=this%nstat(2)

      numtgt=this%numtgt
      itarget_spin=this%itarget_spin
      lucivec=this%lucivec
      lutdip=this%lutdip
      luprop=this%luprop
      lupropw=this%lupropw
      ukrmolp_ints=this%ukrmolp_ints
      max_multipole=this%max_multipole
      isw = this%isw
      gflag = this%gflag

      end subroutine return_namelist_variables

!     *******************************
!     Subroutines for denprop_namelist
!     *******************************
      subroutine init_target_namelist(this, maxtgt, maxpfg)
      implicit none

!     Arguments
      class(denprop_namelist), intent(inout) :: this
      integer :: maxtgt, maxpfg

      if (.not. allocated(this%NFTSOR) ) allocate( this%NFTSOR(maxtgt) )
      if (.not. allocated(this%NTGTF ) ) allocate( this%NTGTF(maxtgt)  )
      if (.not. allocated(this%NTGTS ) ) allocate( this%NTGTS(maxtgt)  )
      if (.not. allocated(this%NTGTL ) ) allocate( this%NTGTL(maxtgt)  )
      if (.not. allocated(this%NPFLG ) ) allocate( this%NPFLG(maxpfg)  )
      if (.not. allocated(this%numtgt) ) allocate( this%numtgt(maxtgt) )
      if (.not. allocated(this%itarget_spin)) allocate(this%itarget_spin(maxtgt))

 !    Set parameter values
      this%maxtgt=maxtgt
      this%maxpfg=maxpfg

      end subroutine init_target_namelist

      subroutine dealloc_target_namelist(this)
      implicit none

      class(denprop_namelist), intent(inout) :: this

      if (allocated(this%NFTSOR)) deallocate(this%NFTSOR)
      if (allocated(this%NTGTF) ) deallocate(this%NTGTF)
      if (allocated(this%NTGTS) ) deallocate(this%NTGTS)
      if (allocated(this%NTGTL) ) deallocate(this%NTGTL)
      if (allocated(this%NPFLG) ) deallocate(this%NPFLG)
      if (allocated(this%numtgt)) deallocate(this%numtgt)
      if (allocated(this%itarget_spin)) deallocate(this%itarget_spin)

      end subroutine dealloc_target_namelist

      subroutine read_namelist(this, lu_namelist)
      implicit none
!     Arguments
      class(denprop_namelist), intent(inout) :: this
      integer :: lu_namelist

!     Local
      integer :: i

!     denprop namelist /INPUT/
      character(len=80) :: NAME
      integer :: NTGT, NFTG, IWRITE, LNGTH, LENGTH, NFTD, NFTA, NFDA, LCOF, NFTMT, NMSET, NFTDL,   &
     &           LDAR, NFTINT, iplotfg, iqdfg, ISW, IPOL, NFTMOM, NMSET2, iprop, idiag, grounden,  &
     &           geocfn, ksym, gflag, lutdip, verbosity
      logical :: zlast, ukrmolp_ints, QMOLN
      integer, allocatable :: NFTSOR(:),NTGTF(:), NTGTS(:), NTGTL(:), NPFLG(:), numtgt(:), itarget_spin(:)

      NAMELIST /INPUT / NAME, NTGT, NFTSOR, NFTG, NTGTF, NTGTS, NTGTL,   &
     &   NPFLG, IWRITE, LNGTH, LENGTH, NFTD, NFTA, NFDA, LCOF, NFTMT,    &
     &   NMSET, NFTDL, LDAR, NFTINT, iplotfg, iqdfg, ISW, IPOL, NFTMOM,  &
     &   NMSET2, iprop, idiag, grounden, ZLAST, geocfn, ksym, QMOLN,     &
     &   ukrmolp_ints, numtgt, gflag, lutdip, itarget_spin

      allocate(NFTSOR(this%maxtgt),NTGTF(this%maxtgt), NTGTS(this%maxtgt), NTGTL(this%maxtgt), NPFLG(this%maxpfg), &
       numtgt(this%maxtgt), itarget_spin(this%maxtgt))

!     Set default values

      DO I=1, this%MAXTGT
         NFTSOR(I)=0
         NTGTF(I)=i
         NTGTS(I)=1
         NTGTL(I)=1
      END DO

      NTGT=0
      numtgt=0
      itarget_spin=2
      iprop=1
      idiag=0
      iplotfg=0
      iqdfg=0
      grounden=0.0_idp
      ksym=1
      qmoln=.false.
      ukrmolp_ints = .true.
      geocfn=1
      verbosity=1
!
!...... Logical unit numbers
!
      IWRITE=6
      NFTINT=17
      NFDA=20
      NFTD=21
      lutdip=24
      NFTG=25
      NFTA=30
      NFTMT=24
      NFTDL=60
      NFTMOM=61
!
!..... Record sizes for some logical units
!
      LDAR=620
      LNGTH=4095
      LENGTH=819
!
!..... Miscellaneous switches
!
      DO I=1, this%MAXPFG
         NPFLG(I)=0
      END DO
!
!...... No of CI coefficients per box
!
      LCOF=10000
      NMSET=1
      NMSET2=1
      ISW=0
      gflag=1
      IPOL=0
!
!..... End of a stack of &INPUT cards
!
      ZLAST=.FALSE.

      read(lu_namelist,INPUT)

      this%geocfn       = geocfn
      this%name         = name
      this%NFTSOR       = NFTSOR
      this%NTGTF        = NTGTF
      this%NTGTS        = NTGTS
      this%NTGTL        = NTGTL
      this%NTGT         = NTGT
      this%numtgt       = numtgt
      this%itarget_spin = itarget_spin
      this%iprop        = iprop
      this%idiag        = idiag
      this%iplotfg      = iplotfg
      this%iqdfg        = iqdfg
      this%grounden     = grounden
      this%ksym         = ksym
      this%qmoln        = qmoln
      this%ukrmolp_ints = ukrmolp_ints
      this%verbosity    = verbosity
      this%IWRITE       = IWRITE
      this%NFTINT       = NFTINT
      this%NFDA         = NFDA
      this%NFTD         = NFTD
      this%NFTG         = NFTG
      this%NFTA         = NFTA
      this%NFTMT        = NFTMT
      this%NFTDL        = NFTDL
      this%NFTMOM       = NFTMOM
      this%lutdip       = lutdip
      this%LDAR         = LDAR
      this%LNGTH        = LNGTH
      this%LENGTH       = LENGTH
      this%NPFLG        = NPFLG
      this%LCOF         = LCOF
      this%NMSET        = NMSET
      this%NMSET2       = NMSET2
      this%ISW          = ISW
      this%gflag        = gflag
      this%IPOL         = IPOL
      this%ZLAST        = ZLAST

      end subroutine read_namelist

      subroutine write_namelist(this, lu_namelist)
      implicit none
!     Arguments
      class(denprop_namelist), intent(inout) :: this
      integer :: lu_namelist

!     denprop namelist /INPUT/
      character(len=80) :: NAME
      integer :: NTGT, NFTG, IWRITE, LNGTH, LENGTH, NFTD, NFTA, NFDA, LCOF, NFTMT, NMSET, NFTDL,   &
     &           LDAR, NFTINT, iplotfg, iqdfg, ISW, IPOL, NFTMOM, NMSET2, iprop, idiag, grounden,  &
     &           geocfn, ksym, gflag
      logical :: zlast, ukrmolp_ints, QMOLN
      integer, allocatable :: NFTSOR(:),NTGTF(:), NTGTS(:), NTGTL(:), NPFLG(:), numtgt(:)

      NAMELIST /INPUT / NAME, NTGT, NFTSOR, NFTG, NTGTF, NTGTS, NTGTL,   &
     &   NPFLG, IWRITE, LNGTH, LENGTH, NFTD, NFTA, NFDA, LCOF, NFTMT,    &
     &   NMSET, NFTDL, LDAR, NFTINT, iplotfg, iqdfg, ISW, IPOL, NFTMOM,  &
     &   NMSET2, iprop, idiag, grounden, ZLAST, geocfn, ksym, QMOLN,     &
     &   ukrmolp_ints, numtgt, gflag

      allocate(NFTSOR(this%maxtgt), NTGTF(this%maxtgt), NTGTS(this%maxtgt), NTGTL(this%maxtgt), NPFLG(this%maxpfg), &
        numtgt(this%maxtgt))

      geocfn       = this%geocfn
      name         = this%name
      NFTSOR       = this%NFTSOR
      NTGTF        = this%NTGTF
      NTGTS        = this%NTGTS
      NTGTL        = this%NTGTL
      NTGT         = this%NTGT
      numtgt       = this%numtgt
      iprop        = this%iprop
      idiag        = this%idiag
      iplotfg      = this%iplotfg
      iqdfg        = this%iqdfg
      grounden     = this%grounden
      ksym         = this%ksym
      qmoln        = this%qmoln
      ukrmolp_ints = this%ukrmolp_ints
      IWRITE       = this%IWRITE
      NFTINT       = this%NFTINT
      NFDA         = this%NFDA
      NFTD         = this%NFTD
      NFTG         = this%NFTG
      NFTA         = this%NFTA
      NFTMT        = this%NFTMT
      NFTDL        = this%NFTDL
      NFTMOM       = this%NFTMOM
      LDAR         = this%LDAR
      LNGTH        = this%LNGTH
      LENGTH       = this%LENGTH
      NPFLG        = this%NPFLG
      LCOF         = this%LCOF
      NMSET        = this%NMSET
      NMSET2       = this%NMSET2
      ISW          = this%ISW
      gflag        = this%gflag
      IPOL         = this%IPOL
      ZLAST        = this%ZLAST

      write(lu_namelist,INPUT)

      end subroutine write_namelist

      end module class_namelists
