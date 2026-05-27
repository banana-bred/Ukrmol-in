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
module common_obj_gbl
use precisn_gbl
use utils_gbl, only: xermsg
use const_gbl, only: level2, level3, nuc_nam_len, nlmax_ecp, lmax_ecp

   !> \class <darray_1d>
   !> 1D array of reals
   type darray_1d
      !> The array values.
      real(kind=cfp), allocatable :: a(:)
      !> Dimension of the array.
      integer :: d1 = 0
   end type darray_1d

   !> \class <array_2d>
   !> 2D array of reals
   type darray_2d
      !> The array values.
      real(kind=cfp), allocatable :: a(:,:)
      !> First dimension of the array.
      integer :: d1 = 0
      !> Second dimension of the array.
      integer :: d2 = 0
   contains
      !> Allocates the space for the array a. The values d1, d2 must be set before. Returns a non-zero number if allocation was unsucessfull.
      procedure :: init => init_darray_2d
      !> Deallocates space. Returns a non-zero number if the array has not been allocated before.
      procedure :: final => final_darray_2d
   end type darray_2d

   !> \class <iarray_1d>
   !> 1D array of integers
   type iarray_1d
      integer, allocatable :: a(:)
      !> Dimension of the array.
      integer :: d1 = 0
   end type iarray_1d

   !> \class <ecp>
   !>
   !> Parameters for an Effective Core Potential for one nucleus.
   !>
   !> https://www.tc.uni-koeln.de/cgi-bin/pp.pl?language=en,job=getreadme
   !> https://www.molpro.net/manual/doku.php?id=effective_core_potentials
   !>
   !> For each pseudopotential, there are 4 basic parameters:
   !> the number of core electrons, the number of l-projectors (lmax) in the one-
   !> component (non-relativistic or scalar-relativistic) ECP, the number of l-
   !> projectors (lmax') of the SO potential (if given; lmax'=0 otherwise), and
   !> the total number of parameters listed below the commentary line. The latter
   !> parameters provide information on V(lmax) first, and then for the semi-local
   !> one-component and SO potentials, V(l) and V'(l') respectively, in the order
   !> l=0, 1, 2, ..., lmax-1; l'=1, 2, ..., lmax'. For each V(l) or V'(l'), the
   !> number of terms of the form A(i)*r**(n(i)-2)*exp(-a(i)*r**2) is given first,
   !> and then the parameters specifying the individual terms in the sequence
   !> n(1),a(1),A(1);n(2),a(2),A(2);..... Note that the V'(l') are defined as
   !> radial prefactors of l*s terms, i.e., the difference of l+1/2 and l-1/2
   !> potentials, for a given l, is multiplied by 2/(2l+1). In other words the
   !> SO potentials assumed here are just the radial prefactors of the l*s term
   !> without further prefactors. This is consistent with the basis sets from:
   !> https://www.tc.uni-koeln.de/PP/clickpse.en.html. Also see the notes from
   !> NWChem website: https://nwchemgit.github.io/ECP.html#spin-orbit-ecps
   !> 
   !> For each valence basis set of a specified symmetry (s, p, d, ...),
   !> the number of exponents is specified first, then the number of recom-
   !> mended contractions and the contraction patterns (n.m defines the range
   !> of primitives to be contracted). On the following lines, the exponents
   !> of the primitives are given first, and afterwards the sets of contraction
   !> coefficients.
   !>
   !> Currently Core Polarization Potentials (CPPs) and SO ECPs are not implemented.
   !>
   !> All potentials have the generic radial form: V_{l}(r) = \sum_{j}^{n_{l}} c_{j}^{l}*r**(m_{j}^{l}-2)*exp(-g_{j}^{l}*r**2) 
   type ecp_obj
      !> ID of the nucleus, i.e. its number. This must be the same as the 'nuc' value in the parent object nucleus_type.
      integer :: nuc = -1
      !> Number of core electrons replaced by ECPs on this nucleus. Setting n_core = 0 means no ECPs are used.
      integer :: n_core = 0
      !> lmax for the semi-local terms (with angular-momentum projectors P_{l}): \sum_{l=0}^{lmax-1} V_{loc,l}(r) P_{l}
      integer :: lmax = 0
      !> lmax' for the SO terms: \sum_{l=1}^{lmaxp} \Delta V_{SO,l}(r) P_{l} l.s P_{l}
      integer :: lmaxp = 0
      !> n_{lmax} for the local term.
      integer :: nlmax = 0
      !> COEFFICIENTS OF THE LOCAL TERMS. Intended dimension (1:nlmax).
      real(kind=cfp), allocatable :: cloc(:), gloc(:)
      integer, allocatable :: mloc(:)
      !> n_{l} for the semi-local terms. Intended dimension (0:lmax-1).
      integer, allocatable :: nl(:)
      !> COEFFICIENTS OF THE SEMI-LOCAL TERMS. Dimension (1:max(nl),0:lmax-1)
      real(kind=cfp), allocatable :: c(:,:), g(:,:)
      integer, allocatable :: m(:,:)
      !> n_{l}^{'} for the SO terms. Intended dimension (1:lmaxp).
      integer, allocatable :: nlp(:)
      !> COEFFICIENTS OF THE SO TERMS. Dimension (1:max(nlp),1:lmaxp)
      real(kind=cfp), allocatable :: cSO(:,:), gSO(:,:)
      integer, allocatable :: mSO(:,:)
      !>
      logical :: have_ecp = .false., have_semi_local_ecp = .false., have_so_ecp = .false.
   contains
      !> Reads the ECP parameters from the namelist &ECP specified in this routine.
      procedure read_ecp_nml
      !> Evaluates the local part of the ECP on a given radial grid.
      !> Note that the Molpro website [1] includes the screened nuclear potential -(Z-this%n_core)/r as part of the ECP.
      !> However, this term is actually the standard nuclear repulsion integral with the effective charge Z-this%n_core, see cgto_integrals.f90. For this
      !> reason we don't include this potential in our definition of the ECP.
      !> [1] https://www.molpro.net/manual/doku.php?id=effective_core_potentials
      procedure eval_local_ecp
      !> Evaluates the radial part of the semi-local part of the ECP for a given angular momentum and a given radial grid.
      procedure eval_semi_local_ecp
      !> Evaluates the radial part of the SO part of the ECP for a given angular momentum and a given radial grid.
      procedure eval_SO_ecp
      !> Writes the ECP data to a given position in an open file. Only the master writes.
      procedure write_ecp
      !> Reads the ECP data from a given position in an open file. Only the master reads.
      procedure read_ecp
      !> Broadcasts the ECP data from the master rank to all other ranks.
      procedure :: mpi_ecp_bcast
   end type ecp_obj

   !> \class <nucleus_type>
   type nucleus_type
      !> Cartesian coordinates of the center of the nucleus.
      real(kind=cfp) :: center(1:3)
      !> Charge on the nucleus. In theory this is integer, but it may be useful to allow non-integer values for some experiments.
      !> If this "nucleus" represents the center for the basis of the continuum, then Z must be set to 0.
      real(kind=cfp) :: charge
      !> ID of the nucleus, i.e. its number.
      integer :: nuc
      !> Nucleus name, i.e. atomic symbol.
      character(len=nuc_nam_len) :: name
      !> Effective Core Potential (ECP) on this nucleus.
      type(ecp_obj) :: ecp
   contains
      !> This function checks that if nuc = 0 then the corresponding distance from the origin is 0, i.e. nuc = 0 is reserved for CMS functions. 
      procedure :: check => check_nucleus_type
      !> Print the nuclear data.
      procedure :: print => print_nucleus_type
      !> Returns .true. if the nucleus has all attributes of the scattering centre: nuc = 0; Z = 0.0_cfp; nam = 'sc'; center = 0.0_cfp
      procedure :: is_continuum => is_scattering_center
      !> Returns the distance from the CMS of the nucleus
      procedure :: r_CMS => distance_CMS
   end type nucleus_type

   interface resize_array
      module procedure resize_2d_array_cfp
      module procedure resize_3d_array_cfp
      module procedure resize_2d_array_int
      module procedure resize_3d_array_int
   end interface resize_array

   !private routines
   private print_nucleus_type, check_nucleus_type, is_scattering_center, distance_CMS

contains

   subroutine print_orbital_table(energies,num_sym,n_orbs,n_sym,convert_to_eV,print_level)
      use sort_gbl, only: cfp_sort_float_int_1d
      use phys_const_gbl, only: to_eV
      implicit none
      integer, intent(in) :: n_orbs, n_sym, num_sym(2,n_orbs), print_level
      real(kind=cfp) :: energies(n_orbs)
      logical, intent(in) :: convert_to_eV

      integer, allocatable :: permutation(:), nob(:)
      integer :: err, i, num, sym

         allocate(permutation(n_orbs),nob(n_sym),stat=err)
         if (err .ne. 0) call xermsg ('common_obj', 'print_orbital_table', 'Memory allocation failed.',err,1)

         nob = 0

         do i=1,n_orbs
            permutation(i) = i
         enddo !i

         !sort the orbitals according to their energies (if required)
         write(print_level,'(/,10X,"Sorting orbital energies...")')
         call cfp_sort_float_int_1d(n_orbs,energies,permutation)
         write(print_level,'("...done")')
   
         !convert the energies from H to eV
         if (convert_to_eV) energies = energies*to_eV
   
         write(print_level,'(/,5X,"List of Energy-sorted molecular orbitals (number.symmetry) follows:",/)')
         write(print_level,'(5X,"|No |NUM.SYMM|ENERGY[eV]|",8(" ",i1," |"))') (i,i=1,n_sym)
   
         do i=1,n_orbs
            num = num_sym(1,permutation(i))
            sym = num_sym(2,permutation(i))
            nob(sym) = nob(sym) + 1
            if (i > 1) then
               if (energies(i)*energies(i-1) < 0.0_cfp) then
                  write(print_level,'(5X,25("-"))')
               endif
            endif
            write(print_level,'(5X,"|",i3,"|",i6,".",i1,"|",F10.3,"|",8(i3,"|"))') i,num,sym, energies(i), nob(1:n_sym)
         enddo

   end subroutine print_orbital_table

   function init_darray_2d(this)
      implicit none
      class(darray_2d) :: this
      integer :: init_darray_2d

      integer :: err
      
         init_darray_2d = 0

         if (this%d1 .le. 0) then 
            init_darray_2d = 1
            return
         endif

         if (this%d2 .le. 0) then
            init_darray_2d = 2
            return
         endif

         if (allocated(this%a)) deallocate(this%a)

         allocate(this%a(this%d1,this%d2),stat=err)
         if (err .ne. 0) init_darray_2d = 3
      
   end function init_darray_2d

   function final_darray_2d(this)
      implicit none
      class(darray_2d) :: this
      integer :: final_darray_2d

         final_darray_2d = 0

         if (allocated(this%a)) then
            deallocate(this%a)
            this%d1 = 0
            this%d2 = 0
         else
            final_darray_2d = 1
            return
         endif

   end function final_darray_2d

   subroutine read_ecp_nml(this,lunit)
      use const_gbl, only: level1, level3
      use iso_fortran_env,      only: iostat_end
      implicit none
      class(ecp_obj) :: this
      integer, intent(in) :: lunit

      logical :: have_ecp = .false.
      integer :: nuc = 0, n_core = 0, lmax = 0, lmaxp = 0, nlmax = 0
      integer :: nl(0:lmax_ecp) = 0, nlp(1:lmax_ecp) = 0
      real(kind=cfp) :: local_term(3*nlmax_ecp) = 0.0_cfp
      real(kind=cfp) :: semilocal_term(3*nlmax_ecp,0:lmax_ecp) = 0.0_cfp
      real(kind=cfp) :: so_term(3*nlmax_ecp,0:lmax_ecp) = 0.0_cfp

      integer :: i, j, l, eof

      namelist/ecp/nuc,n_core,lmax,lmaxp,nlmax,nl,nlp,local_term,semilocal_term,so_term,have_ecp

      write(level3,'("--------->","read_ecp_nml")')

      this%have_ecp = .false.

      read(unit=lunit,nml=ecp,iostat=eof)

      if (eof == 0) this%have_ecp = have_ecp !namelist read-in correctly
      this%have_semi_local_ecp = .false.
      this%have_so_ecp = .false.

      if (this%have_ecp) then
      if (eof == 0) then !namelist read-in: an ECP is present

         write(level1,'(/,10X,i0,"   CORE POTENTIAL",/)') nuc
   
         if (any([nuc,n_core,lmax,lmaxp,nlmax] < 0)) then
            print *,nuc,n_core,lmax,lmaxp,nlmax
            call xermsg ('common_obj', 'read_ecp_nml', &
                         'At least one of nuc,n_core,lmax,lmaxp,nlmax is negative. Invalid input.', 1, 1)
         endif
   
         this%nuc = nuc
         this%n_core = n_core
         this%lmax = lmax
         this%lmaxp = lmaxp
         this%nlmax = nlmax

         !Extract parameters of the local term
         write(level1,'(/,10X,"g POTENTIAL",/)')
   
         if (this%nlmax == 0) then
            call xermsg ('common_obj', 'read_ecp_nml', &
                         'nlmax is zero. If you wish to exclude the local term set cloc(1) = 0.', 2, 1)
         endif
  
         if (allocated(this%cloc)) deallocate(this%cloc)
         if (allocated(this%mloc)) deallocate(this%mloc)
         if (allocated(this%gloc)) deallocate(this%gloc)
         allocate(this%cloc(this%nlmax), this%mloc(this%nlmax), this%gloc(this%nlmax))

         j = 0
         do i = 1, 3*this%nlmax, 3
            j = j + 1
            this%mloc(j) = nint(local_term(i))
            this%gloc(j) = local_term(i+1)
            this%cloc(j) = local_term(i+2)
         enddo
   
         write(level1,'(1X,"POWERS      ",3X,20(i8,7X))') int(this%mloc(1:this%nlmax))
         write(level1,'(1X,"EXPONENTIALS",3X,20(f8.4,7X))') this%gloc(1:this%nlmax)
         write(level1,'(1X,"COEFFICIENTS",3X,20(f8.4,7X))') this%cloc(1:this%nlmax)
   
         !Extract parameters of the semi-local term
         if (this%lmax > 0) then

            this%have_semi_local_ecp = .true.

            if (any(nl(0:this%lmax-1) <= 0)) then
               print *,nl
               call xermsg ('common_obj', 'read_ecp_nml', 'At least one of nl(l) is zero.', 4, 1)
            endif
      
            if (allocated(this%nl)) deallocate(this%nl)
            if (allocated(this%c)) deallocate(this%c)
            if (allocated(this%m)) deallocate(this%m)
            if (allocated(this%g)) deallocate(this%g)
            i = 3*maxval(nl(0:this%lmax-1))
            allocate(this%m(i,0:this%lmax-1),this%g(i,0:this%lmax-1),this%c(i,0:this%lmax-1),this%nl(0:this%lmax-1))

            this%nl(0:this%lmax-1) = nl(0:this%lmax-1)
      
            do l = 0, this%lmax-1
   
               write(level1,'(/,10X,"L=",i0," - g POTENTIAL",/)') l
      
               j = 0
               do i = 1, 3*this%nl(l), 3
                  j = j + 1
                  this%m(j,l) = nint(semilocal_term(i,l))
                  this%g(j,l) = semilocal_term(i+1,l)
                  this%c(j,l) = semilocal_term(i+2,l)
               enddo
   
               write(level1,'(1X,"POWERS      ",3X,20(i8,7X))') int(this%m(1:this%nl(l),l))
               write(level1,'(1X,"EXPONENTIALS",3X,20(f8.4,7X))') this%g(1:this%nl(l),l)
               write(level1,'(1X,"COEFFICIENTS",3X,20(f8.4,7X))') this%c(1:this%nl(l),l)            
            enddo
         endif !lmax
   
         !Extract parameters of the SO term
         if (this%lmaxp > 0) then
   
            this%have_so_ecp = .true.

            write(level1,'(/,10X,i0,3X,"SPIN-ORBIT POTENTIAL",/)') nuc
   
            if (any(nlp(1:this%lmaxp) <= 0)) then
               print *,nlp
               call xermsg ('common_obj', 'read_ecp_nml', 'At least one of nlp(l) is zero.', 4, 1)
            endif
      
            if (allocated(this%nlp)) deallocate(this%nlp)
            if (allocated(this%cSO)) deallocate(this%cSO)
            if (allocated(this%mSO)) deallocate(this%mSO)
            if (allocated(this%gSO)) deallocate(this%gSO)
            i = 3*maxval(nlp(1:this%lmaxp))
            allocate(this%mSO(i,1:this%lmaxp),this%gSO(i,1:this%lmaxp),this%cSO(i,1:this%lmaxp),this%nlp(1:this%lmaxp))

            this%nlp(1:this%lmaxp) = nlp(1:this%lmaxp)
      
            do l = 1, this%lmaxp
      
               write(level1,'(/,10X,"L=",i0," POTENTIAL",/)') l
   
               j = 0
               do i = 1, 3*this%nlp(l), 3
                  j = j + 1
                  this%mSO(j,l) = nint(so_term(i,l))
                  this%gSO(j,l) = so_term(i+1,l)
                  this%cSO(j,l) = so_term(i+2,l)
               enddo
               
               write(level1,'(1X,"POWERS      ",3X,20(i8,7X))') int(this%mSO(1:this%nlp(l),l))
               write(level1,'(1X,"EXPONENTIALS",3X,20(f8.4,7X))') this%gSO(1:this%nlp(l),l)
               write(level1,'(1X,"COEFFICIENTS",3X,20(f8.4,7X))') this%cSO(1:this%nlp(l),l)            
            enddo
         endif !lmaxp

      else if (eof == iostat_end) then

         write(level3,'("Searching for another ECP...")')

      else !If there is an issue not related to end of file then check with a non iostat call to read the nml

         rewind(lunit)
         do 
            read(unit=lunit,nml=ecp)
         enddo
      endif
      else
         write(level1,'(/,10X,"   NO ECP ON THIS NUCLEUS",/)')
      endif !have_ecp

      write(level3,'("<---------","read_ecp_nml")')

   end subroutine read_ecp_nml

   subroutine write_ecp(this,lunit,posit)
      implicit none
      class(ecp_obj) :: this
      integer, intent(in) :: lunit, posit

         write(level3,'("--------->","write_ecp")')

         write(lunit,pos=posit,err=10) merge(1, 0, this%have_ecp), merge(1, 0, this%have_semi_local_ecp), &
                                       merge(1, 0, this%have_so_ecp)

         write(lunit,err=10) this%nuc, this%n_core, this%lmax, this%lmaxp, this%nlmax

         if (this%have_ecp) then
            write(lunit,err=10) this%cloc(1:this%nlmax), this%mloc(1:this%nlmax), this%gloc(1:this%nlmax)
         endif

         if (this%have_semi_local_ecp) then
            write(lunit,err=10) this%nl(0:this%lmax-1)
            write(lunit,err=10) this%m(:,0:this%lmax-1),this%g(:,0:this%lmax-1),this%c(:,0:this%lmax-1)
         endif

         if (this%have_so_ecp) then
            write(lunit,err=10) this%nlp(1:this%lmaxp)
            write(lunit,err=10) this%mSO(:,1:this%lmaxp),this%gSO(:,1:this%lmaxp),this%cSO(:,1:this%lmaxp)
         endif

         write(level3,'("<---------","write_ecp")')

         return

 10      call xermsg('common_obj','write_ecp','Error writing the ecp_obj data into the file and position given.',2,1)

   end subroutine write_ecp

   subroutine read_ecp(this,lunit,posit)
      implicit none
      class(ecp_obj) :: this
      integer, intent(in) :: lunit, posit

      integer :: il(3), i, err

         write(level3,'("--------->","read_ecp")')

         read(lunit,pos=posit,err=10) il(1:3)

         this%have_ecp = (il(1) == 1)
         this%have_semi_local_ecp = (il(2) == 1)
         this%have_so_ecp = (il(3) == 1)

         read(lunit,err=10) this%nuc, this%n_core, this%lmax, this%lmaxp, this%nlmax

         if (this%have_ecp) then
            if (allocated(this%cloc)) deallocate(this%cloc)
            if (allocated(this%mloc)) deallocate(this%mloc)
            if (allocated(this%gloc)) deallocate(this%gloc)
            allocate(this%cloc(1:this%nlmax), this%mloc(1:this%nlmax), this%gloc(1:this%nlmax),stat=err)
            if (err .ne. 0) call xermsg('common_obj','read_ecp', 'Memory allocation 1 failed.',err,1)
            read(lunit,err=10) this%cloc(1:this%nlmax), this%mloc(1:this%nlmax), this%gloc(1:this%nlmax)
         endif

         if (this%have_semi_local_ecp) then
            if (allocated(this%m)) deallocate(this%m)
            if (allocated(this%g)) deallocate(this%g)
            if (allocated(this%c)) deallocate(this%c)
            if (allocated(this%nl)) deallocate(this%nl)

            allocate(this%nl(0:this%lmax-1),stat=err)
            if (err .ne. 0) call xermsg('common_obj','read_ecp', 'Memory allocation 2 failed.',err,1)
            read(lunit,err=10) this%nl(0:this%lmax-1)

            i = 3*maxval(this%nl(0:this%lmax-1))
            allocate(this%m(i,0:this%lmax-1),this%g(i,0:this%lmax-1),this%c(i,0:this%lmax-1),stat=err)
            if (err .ne. 0) call xermsg('common_obj','read_ecp', 'Memory allocation 3 failed.',err,1)

            read(lunit,err=10) this%m(1:i,0:this%lmax-1),this%g(1:i,0:this%lmax-1),this%c(1:i,0:this%lmax-1)
         endif

         if (this%have_so_ecp) then
            if (allocated(this%mSO)) deallocate(this%mSO)
            if (allocated(this%gSO)) deallocate(this%gSO)
            if (allocated(this%cSO)) deallocate(this%cSO)
            if (allocated(this%nlp)) deallocate(this%nlp)

            allocate(this%nlp(1:this%lmaxp),stat=err)
            if (err .ne. 0) call xermsg('common_obj','read_ecp', 'Memory allocation 4 failed.',err,1)
            read(lunit,err=10) this%nlp(1:this%lmaxp)

            i = 3*maxval(this%nlp(1:this%lmaxp))
            allocate(this%mSO(i,1:this%lmaxp),this%gSO(i,1:this%lmaxp),this%cSO(i,1:this%lmaxp),stat=err)
            if (err .ne. 0) call xermsg('common_obj','read_ecp', 'Memory allocation 5 failed.',err,1)

            read(lunit,err=10) this%mSO(1:i,1:this%lmaxp),this%gSO(1:i,1:this%lmaxp),this%cSO(1:i,1:this%lmaxp)
         endif

         write(level3,'("<---------","read_ecp")')

         return

 10      call xermsg('common_obj','read_ecp','Error reading the ecp_obj data from the file and position given.',2,1)

   end subroutine read_ecp

   subroutine mpi_ecp_bcast(this)
      use mpi_gbl
      implicit none
      class(ecp_obj) :: this

      integer :: err, i

         !master broadcasts all its data to the other processes
         call mpi_mod_bcast(this%have_ecp,master)
         call mpi_mod_bcast(this%have_semi_local_ecp,master)
         call mpi_mod_bcast(this%have_so_ecp,master)
         call mpi_mod_bcast(this%nuc,master)
         call mpi_mod_bcast(this%n_core,master)
         call mpi_mod_bcast(this%lmax,master)
         call mpi_mod_bcast(this%lmaxp,master)
         call mpi_mod_bcast(this%nlmax,master)

         if (myrank .ne. master) then

            if (this%have_ecp) then
               if (allocated(this%cloc)) deallocate(this%cloc)
               if (allocated(this%mloc)) deallocate(this%mloc)
               if (allocated(this%gloc)) deallocate(this%gloc)
               allocate(this%cloc(1:this%nlmax), this%mloc(1:this%nlmax), this%gloc(1:this%nlmax),stat=err)
               if (err .ne. 0) call xermsg('common_obj','mpi_ecp_bcast', 'Memory allocation 1 failed.',err,1)
            endif

            if (this%have_semi_local_ecp) then
               if (allocated(this%m)) deallocate(this%m)
               if (allocated(this%g)) deallocate(this%g)
               if (allocated(this%c)) deallocate(this%c)
               if (allocated(this%nl)) deallocate(this%nl)
   
               allocate(this%nl(0:this%lmax-1),stat=err)
               if (err .ne. 0) call xermsg('common_obj','mpi_ecp_bcast', 'Memory allocation 2 failed.',err,1)
            endif

            if (this%have_so_ecp) then
               if (allocated(this%mSO)) deallocate(this%mSO)
               if (allocated(this%gSO)) deallocate(this%gSO)
               if (allocated(this%cSO)) deallocate(this%cSO)
               if (allocated(this%nlp)) deallocate(this%nlp)
   
               allocate(this%nlp(1:this%lmaxp),stat=err)
               if (err .ne. 0) call xermsg('common_obj','mpi_ecp_bcast', 'Memory allocation 3 failed.',err,1)
            endif

         endif

         if (this%have_ecp) then
            call mpi_mod_bcast(this%cloc,master)
            call mpi_mod_bcast(this%mloc,master)
            call mpi_mod_bcast(this%gloc,master)
         endif

         if (this%have_semi_local_ecp) then
            call mpi_mod_bcast(this%nl,master)
            i = 3*maxval(this%nl(0:this%lmax-1))
            if (myrank .ne. master) then
               allocate(this%m(i,0:this%lmax-1),this%g(i,0:this%lmax-1),this%c(i,0:this%lmax-1),stat=err)
               if (err .ne. 0) call xermsg('common_obj','mpi_ecp_bcast', 'Memory allocation 4 failed.',err,1)
            endif
            call mpi_mod_bcast(this%c,master)
            call mpi_mod_bcast(this%m,master)
            call mpi_mod_bcast(this%g,master)
         endif

         if (this%have_so_ecp) then
            call mpi_mod_bcast(this%nlp,master)
            i = 3*maxval(this%nlp(1:this%lmaxp))
            if (myrank .ne. master) then
               allocate(this%mSO(i,1:this%lmaxp),this%gSO(i,1:this%lmaxp),this%cSO(i,1:this%lmaxp),stat=err)
               if (err .ne. 0) call xermsg('common_obj','mpi_ecp_bcast', 'Memory allocation 5 failed.',err,1)
            endif
            call mpi_mod_bcast(this%cSO,master)
            call mpi_mod_bcast(this%mSO,master)
            call mpi_mod_bcast(this%gSO,master)
         endif

   end subroutine mpi_ecp_bcast

   subroutine eval_local_ecp(this,r,local_ecp)
      use const_gbl, only: int_rel_prec
      implicit none
      class(ecp_obj) :: this
      real(kind=cfp), intent(in) :: r(:)
      real(kind=cfp), allocatable :: local_ecp(:)

      integer :: n, i, j

         n = size(r)

         if (this%n_core <= 0) then
            call xermsg ('common_obj', 'eval_local_ecp', 'ECP is not initialized.', 2, 1)
         endif

         if (allocated(local_ecp)) deallocate(local_ecp)
         allocate(local_ecp(n))

         local_ecp = 0.0_cfp
         do i = 1, n
            do j = 1, this%nlmax
               local_ecp(i) = local_ecp(i) + this%cloc(j)*r(i)**(this%mloc(j)-2)*exp(-this%gloc(j)*r(i)*r(i))
            enddo !j
         enddo !i

         do i = n, 1, -1
            if (abs(local_ecp(i)) > int_rel_prec) then
               write(level3,'(10X,"Local ECP range approx.: ",e25.15," Bohr")') r(i)
               exit
            endif
         enddo

   end subroutine eval_local_ecp

   subroutine eval_semi_local_ecp(this,r,semi_local_ecp)
      use const_gbl, only: int_rel_prec
      implicit none
      class(ecp_obj) :: this
      real(kind=cfp), intent(in) :: r(:)
      real(kind=cfp), allocatable :: semi_local_ecp(:,:)

      integer :: n, i, j, l

         n = size(r)

         if (this%lmax <= 0) then
            call xermsg ('common_obj', 'eval_semi_local_ecp', 'Semi-local ECP is not initialized.', 1, 1)
         endif

         if (allocated(semi_local_ecp)) deallocate(semi_local_ecp)
         allocate(semi_local_ecp(n,0:this%lmax-1))

         semi_local_ecp = 0.0_cfp

         do l = 0, this%lmax-1
            do i = 1, n
               do j = 1, this%nl(l)
                  semi_local_ecp(i,l) = semi_local_ecp(i,l) + this%c(j,l)*r(i)**(this%m(j,l)-2)*exp(-this%g(j,l)*r(i)*r(i))
               enddo !j
            enddo !i
         enddo !l

         do l = 0, this%lmax-1
            do i = n, 1, -1
               if (abs(semi_local_ecp(i,l)) > int_rel_prec) then
                  write(level3,'(10X,"Semi-local ECP range approx.: ",e25.15," Bohr"," L = ",i2)') r(i), l
                  exit
               endif
            enddo
         enddo !l

   end subroutine eval_semi_local_ecp

   subroutine eval_SO_ecp(this,r,SO_ecp)
      use const_gbl, only: int_rel_prec
      implicit none
      class(ecp_obj) :: this
      real(kind=cfp), intent(in) :: r(:)
      real(kind=cfp), allocatable :: SO_ecp(:,:)

      integer :: n, i, j, l

         n = size(r)

         if (this%lmaxp <= 0) then
            call xermsg ('common_obj', 'eval_SO_ecp', 'SO ECP is not initialized.', 1, 1)
         endif

         if (allocated(SO_ecp)) deallocate(SO_ecp)
         allocate(SO_ecp(n,1:this%lmaxp))

         SO_ecp = 0.0_cfp

         do l = 1, this%lmaxp
            do i = 1, n
               do j = 1, this%nlp(l)
                  SO_ecp(i,l) = SO_ecp(i,l) + this%cSO(j,l)*r(i)**(this%mSO(j,l)-2)*exp(-this%gSO(j,l)*r(i)*r(i))
               enddo !j
            enddo !i
         enddo !l

         do l = 0, this%lmaxp
            do i = n, 1, -1
               if (abs(SO_ecp(i,l)) > int_rel_prec) then
                  write(level3,'(10X,"Spin-orbit ECP range approx.: ",e25.15," Bohr"," L = ",i2)') r(i), l
                  exit
               endif
            enddo
         enddo !l

   end subroutine eval_SO_ecp

   !> \todo Sort out the nuc == 0 question.
   function check_nucleus_type(this)
      implicit none
      class(nucleus_type) :: this
      integer :: check_nucleus_type

         check_nucleus_type = 0

         !nuc = 0 is reserved only for nuclei centered on the CMS. todo This is not true and it should be probably removed.
!         if (dot_product(this%center,this%center) == 0.0_cfp) then
!            if (this%nuc .ne. 0) then
!               check_nucleus_type = 1
!               call this%print
!               call xermsg('common_obj','check_nucleus_type','nuc = 0 must be used for nuclei centered on the CMS.',1,1)
!            endif
!         endif

         if (this%nuc == 0 .and. dot_product(this%center,this%center) /= 0.0_cfp) then
            check_nucleus_type = 2
            call this%print
            call xermsg ('common_obj', 'check_nucleus_type', &
                         'nuc = 0 (i.e. nucleus at the CMS), but the distance from the origin is not 0.', 2, 1)
         endif

         !if we're using ECPs we must ensure that the intended potential input matches the actual nucleus
         if (this%ecp%n_core > 0 .and. this%nuc /= this%ecp%nuc) then
            check_nucleus_type = 3
            call this%print
            call xermsg ('common_obj', 'check_nucleus_type', &
                         'nuc attributes of nucleus_type and of the corresponding ecp dont match.', 3, 1)
         endif

   end function check_nucleus_type

   subroutine print_nucleus_type(this)
      implicit none
      class(nucleus_type) :: this

         write(level2,'("Nucleus name: ",a)') adjustl(this%name)
         write(level2,'("Nucleus number: ",i5)') this%nuc
         write(level2,'("Nucleus Z: ",f20.15)') this%charge
         write(level2,'("Nucleus center: ",3f25.15)') this%center(:)

   end subroutine print_nucleus_type

   function is_scattering_center(this)
      implicit none
      class(nucleus_type) :: this
      logical :: is_scattering_center

      real(kind=cfp) :: r

         is_scattering_center = .false.

         r = dot_product(this%center,this%center)
         if (this%name .eq. 'sc' .and. this%nuc .eq. 0 .and. r .eq. 0.0_cfp) is_scattering_center = .true.

   end function

   function distance_CMS(this)
      implicit none
      class(nucleus_type) :: this
      real(kind=cfp) :: distance_CMS

         distance_CMS = sqrt(dot_product(this%center,this%center))

   end function

   function resize_2d_array_cfp(array,d1,d2)
      implicit none
      integer, intent(in) :: d1, d2
      real(kind=cfp), allocatable :: array(:,:)
      integer :: resize_2d_array_cfp

      logical :: do_allocation
      integer :: err

         resize_2d_array_cfp = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2),stat=err)
            resize_2d_array_cfp = err
         endif

   end function resize_2d_array_cfp

   function resize_3d_array_cfp(array,d1,d2,d3)
      implicit none
      integer, intent(in) :: d1, d2, d3
      real(kind=cfp), allocatable :: array(:,:,:)
      integer :: resize_3d_array_cfp

      logical :: do_allocation
      integer :: err

         resize_3d_array_cfp = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2 .or. size(array,3) < d3) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2,d3),stat=err)
            resize_3d_array_cfp = err
         endif

   end function resize_3d_array_cfp

   function resize_2d_array_int(array,d1,d2)
      implicit none
      integer, intent(in) :: d1, d2
      integer, allocatable :: array(:,:)
      integer :: resize_2d_array_int

      logical :: do_allocation
      integer :: err

         resize_2d_array_int = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2),stat=err)
            resize_2d_array_int = err
         endif

   end function resize_2d_array_int

   function resize_3d_array_int(array,d1,d2,d3)
      implicit none
      integer, intent(in) :: d1, d2, d3
      integer, allocatable :: array(:,:,:)
      integer :: resize_3d_array_int

      logical :: do_allocation
      integer :: err

         resize_3d_array_int = 0

         if (.not.allocated(array)) then
            do_allocation = .true.
         elseif (size(array,1) < d1 .or. size(array,2) < d2 .or. size(array,3) < d3) then
            do_allocation = .true.
            deallocate(array)
         else
            do_allocation = .false.
         endif
   
         if (do_allocation) then
            allocate(array(d1,d2,d3),stat=err)
            resize_3d_array_int = err
         endif

   end function resize_3d_array_int

   subroutine resize_copy_2d_array(array,d1,d2)
       implicit none
       real(kind=cfp), allocatable :: array(:,:)
       integer, intent(in) :: d1,d2

       real(kind=cfp), allocatable :: tmp(:,:)
       integer :: err, d1_al, d2_al, d1_old, d2_old

          d1_al = d1
          d2_al = d2

          d1_old = size(array,1)
          d2_old = size(array,2)

          if (.not.(allocated(array))) then
             allocate(array(d1_al,d2_al),stat=err)
             if (err .ne. 0) call xermsg('common_obj','resize_copy_2d_array','Memory allocation 1 failed.',err,1)
             array = 0.0_cfp
          else
             d1_al = max(d1_al,d1_old)
             d2_al = max(d2_al,d2_old)

             call move_alloc(array,tmp)

             allocate(array(d1_al,d2_al),stat=err)
             if (err .ne. 0) call xermsg('common_obj','resize_copy_2d_array','Memory allocation 2 failed.',err,1)
             array = 0.0_cfp

             array(1:d1_old,1:d2_old) = tmp(1:d1_old,1:d2_old)
             deallocate(tmp)
          endif

   end subroutine resize_copy_2d_array

end module common_obj_gbl
