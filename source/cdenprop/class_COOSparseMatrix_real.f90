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
module class_COOSparseMatrix_real

  use iso_fortran_env, only: output_unit
  use mpi_gbl,         only: mpi_xermsg
  use precisn

  implicit none

  integer, parameter :: idp = wp
  integer, parameter :: ip = shortint !longint
  integer(ip), private :: test
  integer, parameter :: maxint = huge(test)-2

  private


  type, public :: COOMatrix_real
     integer, private  :: M , N
     integer, private :: estimated_number_of_nnz
     integer :: nnz
     integer, private :: allocated_size
     logical, private :: is_matrix_initialized = .false.
     real(kind=idp), allocatable  :: value_list(:)
     integer(ip), allocatable :: index_i_list(:),index_j_list(:)
     integer :: iwrite = output_unit
  contains
     procedure, public :: init_matrix
     procedure, public :: end_matrix_construction
     procedure, public :: insert_single
     procedure, public :: print_structure
     procedure, public :: get_indexed_element
     procedure, public :: search_by_i_j
     procedure, public :: set_iwrite
  end type COOMatrix_real

  contains

    subroutine set_iwrite(this,iwrite_inp)
    class(COOMatrix_real), intent(inout) :: this
    integer, intent(in) :: iwrite_inp

       this%iwrite = iwrite_inp

    end subroutine set_iwrite

    subroutine init_matrix(this,M,N,est_nnz)
    class(COOMatrix_real), intent(inout) :: this
    integer, intent(in) :: M,N,est_nnz
    integer :: err

    if (this % is_matrix_initialized) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'init_matrix', &
                       'Trying to initialize a matrix already initialized', 1, 1)
    endif

    if (est_nnz > (M*N)) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'init_matrix', &
                       'Trying to allocate more space than elements in the matrix', 2, 1)
    endif

    write(this%iwrite,*) "SparseMatrix: Creating a matrix of dimension", M," ", N, " ", est_nnz
    if (ip .eq. shortint) write(this%iwrite,*) "SparseMatrix: integers will be stored as SHORTINT"
    if (ip .eq. longint) write(this%iwrite,*) "SparseMatrix: integers will be stored as LONGINT"

    if (est_nnz < 10) then
       this%estimated_number_of_nnz = M*N
    else
       this%estimated_number_of_nnz = est_nnz
    endif

!     this%estimated_number_of_nnz = est_nnz
    this%M = M
    this%N = N

    allocate(this%index_i_list(this%estimated_number_of_nnz),stat=err)
    if (err .ne. 0) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'init_matrix', &
                       'allocate this%index_i_list failed', this%estimated_number_of_nnz, 1)
    endif
    allocate(this%index_j_list(this%estimated_number_of_nnz),stat=err)
    if (err .ne. 0) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'init_matrix', &
                       'allocate this%index_j_list failed', this%estimated_number_of_nnz, 1)
    endif
    allocate(this%value_list(this%estimated_number_of_nnz),stat=err)
    if (err .ne. 0) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'init_matrix', &
                       'allocate this%index_list failed', this%estimated_number_of_nnz, 1)
    endif

    this%index_i_list = 0
    this%index_j_list = 0
    this%value_list = 0.d0

    this%allocated_size = this%estimated_number_of_nnz
    this%nnz = 0
    this%is_matrix_initialized = .true.

    end subroutine init_matrix
!
! When we are not inserting more elements it may be useful to resize all structures just to the number of non-zeros
    subroutine end_matrix_construction(this)
    class(COOMatrix_real), intent(inout) :: this
    integer :: new_size, old_size
    real(kind=idp), allocatable :: tmp_values(:)
    integer(ip), allocatable :: tmp_i(:),tmp_j(:)

!     Reallocate all the structures only to the nnz
    new_size = this%nnz
    old_size = this%allocated_size
    write(this%iwrite,*) "SparseMatrix: Memory saving of ",(old_size - new_size)*(3*8), " bytes"
    allocate(tmp_values(new_size),tmp_i(new_size),tmp_j(new_size))
    tmp_values = 0.d0
    tmp_i = 0
    tmp_j = 0


! COPY all arrays to the temporal ones
    tmp_values(1:new_size) = this%value_list(1:new_size)
    tmp_i(1:new_size) = this%index_i_list(1:new_size)
    tmp_j(1:new_size) = this%index_j_list(1:new_size)

! Deallocate structure arrays, and allocate again with new_size
    deallocate(this%index_i_list)
    deallocate(this%index_j_list)
    deallocate(this%value_list)

    allocate(this%index_i_list(new_size))
    allocate(this%index_j_list(new_size))
    allocate(this%value_list(new_size))

    this%index_i_list = 0
    this%index_j_list = 0
    this%value_list = 0.d0

! COPY THE VALUES BACK INTO THE MATRIX STRUCTURE
    this%allocated_size = new_size
    this%value_list(1:new_size) = tmp_values(1:new_size)
    this%index_i_list(1:new_size) = tmp_i(1:new_size)
    this%index_j_list(1:new_size) = tmp_j(1:new_size)
    deallocate(tmp_i,tmp_j,tmp_values)

    end subroutine end_matrix_construction

    subroutine insert_single(this, pos_i, pos_j, value)
    class(COOMatrix_real), intent(inout) :: this
    integer, intent(in) :: pos_i,pos_j
    real(kind=idp), intent(in) :: value
    integer :: position, new_size, old_size
    real(kind=idp), allocatable :: tmp_values(:)
    integer(ip), allocatable :: tmp_i(:),tmp_j(:)

    if (.not. this % is_matrix_initialized) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'insert_single', 'Trying to insert before initializing the matrix', 1, 1)
    endif

    if ((pos_i > this%M).or.(pos_j > this%N)) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'insert_single', 'Out of place index', 2, 1)
    endif

    position = this%nnz+1

! if number of nnz is smaller than the allocated size
! we need to allocate more memory :: first guess add another estimated_number_of_nnz
    if ((position > this%allocated_size)) then
         write(this%iwrite,*) &
            "COOMatrix_real PENALTY IN MATRIX INSERTION, ALLOCATING DOUBLE THE MEMORY, please increase nnz :: ", &
            this % nnz, &
            this % allocated_size
! Calculate new size, allocate temporal arrays and initialize to 0
       old_size = this%allocated_size
       new_size = this%allocated_size * 2
       allocate(tmp_values(new_size),tmp_i(new_size),tmp_j(new_size))
       tmp_values = 0.d0
       tmp_i = 0
       tmp_j = 0

! COPY all arrays to the temporal ones
       tmp_values(1:old_size) = this%value_list(1:old_size)
       tmp_i(1:old_size) = this%index_i_list(1:old_size)
       tmp_j(1:old_size) = this%index_j_list(1:old_size)

! Deallocate structure arrays, and allocate again with new_size
       deallocate(this%index_i_list)
       deallocate(this%index_j_list)
       deallocate(this%value_list)

       allocate(this%index_i_list(new_size))
       allocate(this%index_j_list(new_size))
       allocate(this%value_list(new_size))

       this%index_i_list = 0
       this%index_j_list = 0
       this%value_list = 0.d0

! COPY THE VALUES BACK INTO THE MATRIX STRUCTURE
       this%allocated_size = new_size
       this%value_list(1:old_size) = tmp_values(1:old_size)
       this%index_i_list(1:old_size) = tmp_i(1:old_size)
       this%index_j_list(1:old_size) = tmp_j(1:old_size)
! DEALLOCATE TEMPORAL ARRAYS
       deallocate(tmp_i,tmp_j,tmp_values)

!UPDATE THE ALLOCATED SIZE
       this%allocated_size = new_size
    endif

!  Insert values and update nnz
    if (abs(pos_i) > maxint) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'insert_single', 'pos_i out of range', pos_i, 1)
    else if (abs(pos_j) > maxint) then
       call mpi_xermsg('class_COOSparseMatrix_real', 'insert_single', 'pos_j out of range', pos_j, 1)
    endif

    this%index_i_list(position) = pos_i
    this%index_j_list(position) = pos_j
    this%value_list(position) = value
    this%nnz = this%nnz + 1

    if (value.eq.0.d0) then
       write(this%iwrite,*) "INSERTING A 0 into the sparse matrix"
    endif


    end subroutine insert_single

    subroutine get_indexed_element(this,index_in,i,j,value)
    class(COOMatrix_real), intent(inout) :: this
    integer, intent(in) :: index_in
    integer, intent(out):: i,j
    real(kind=idp), intent(out):: value

    if (.not. this % is_matrix_initialized) then
       write(this%iwrite,*) "COOMatrix_real ERROR ACCESING:: Trying to access before initializing the matrix"
       i=0
       j=0
       value = 0.d0
       return

    endif
    if (this%nnz < index_in) then
       write(this%iwrite,*) "COOMatrix_real ERROR ACCESING:: Trying to access an element with index bigger than nnz"
       i=0
       j=0
       value = 0.d0
       return

    endif

    i = this%index_i_list(index_in)
    j = this%index_j_list(index_in)
    value = this%value_list(index_in)

    end subroutine get_indexed_element

    subroutine search_by_i_j(this,i,j,index_out,value)
    class(COOMatrix_real), intent(inout) :: this
    integer, intent(out) :: index_out
    integer, intent(in):: i,j
    integer :: i_iterate
    real(kind=idp), intent(out):: value

    index_out = 0
    value = 0.d0

    if (.not. this % is_matrix_initialized) then
       write(this%iwrite,*) "COOMatrix_real ERROR SEARCHING:: Trying to access before initializing the matrix"
       return

    endif

    if ((this%M<i).or.(this%N<j)) then
       write(this%iwrite,*) "COOMatrix_real ERROR SEARCHING:: Trying to access an element with ij bigger MN"
       return

    endif


    do i_iterate=1,this%nnz
       if ((this%index_i_list(i_iterate) == i).and.(this%index_j_list(i_iterate) == j)) then
          value = this%value_list(i_iterate)
          index_out = i_iterate
          exit

       endif

    enddo

    end subroutine search_by_i_j

    subroutine print_structure(this)
    class(COOMatrix_real), intent(inout) :: this
    integer :: i

    if (.not. this % is_matrix_initialized) then
        write(this%iwrite,*) "COOMatrix_real ERROR PRINTING:: Trying to print before initializing the matrix"
        return

    endif

    write(this%iwrite,*) "COOMatrix_real STRUCTURE------------------------------------"
    write(this%iwrite,*) " DIM M,N, SIZE :: " , this%M, this%N, this%M*this%N
    write(this%iwrite,*) " NNZ :: " , this%nnz
    write(this%iwrite,*) " Matrix occupancy :: ", 100.d0*this%nnz/float(this%M*this%N) , " %"
    write(this%iwrite,*) " index_i / index_j / value :: "

    write(this%iwrite, '((4x,I9,I9,E25.15))') (this%index_i_list(i), this%index_j_list(i), this%value_list(i), i=1,this%nnz)

    write(this%iwrite,*) " Sizes of arrays :: ",size(this%index_i_list),size(this%index_j_list),size(this%value_list)
    write(this%iwrite,*) " Allocated size :: ",this%allocated_size
    if (ip .eq. shortint) write(this%iwrite,*) "Integers are stored as SHORTINT"
    if (ip .eq. longint) write(this%iwrite,*) "Integers are stored as LONGINT"
    write(this%iwrite,*) "COOMatrix_real STRUCTURE------------------------------------"

    end subroutine print_structure

end module class_COOSparseMatrix_real
