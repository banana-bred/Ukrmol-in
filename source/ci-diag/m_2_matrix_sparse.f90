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
!=======================================================================================!
! Author: Pavlos G. Galiatsatos.                                                        !
! Corrections, modifications and improvements: Michael Lysaght .                        !
! Date: 2011/06/30.                                                                     !
!=======================================================================================!

module m_2_matrix_sparse

  use precisn 

  contains

!=======================================================================================!
!                                                                                       !
! Subroutine: "matrix_coo_ram_half".                                                    !
!                                                                                       !
! a. Input:                                                                             !
!                                                                                       !
! integer, intent(in) :: file_number_from_scatci                                        !         
!                                                                                       !              
! b. Output:                                                                            ! 
!                                                                                       !
! integer(kind=shortint), intent(out) :: dimen_matrix                                   !
! integer(kind=longint) , intent(out) :: non_zero_ham_elem_half                         !
! integer(kind=shortint), intent(out), allocatable, dimension(:) :: row_index_array     !
! integer(kind=shortint), intent(out), allocatable, dimension(:) :: column_index_array  ! 
! real(kind=wp), intent(out), allocatable, dimension(:) :: ham_elem_array               !
!                                                                                       !
! c. Purpose:                                                                           !
!                                                                                       !
! Writing in RAM the indices of rows, columns                                           !
! and the correspoding hamiltonian matrix elements.                                     !
! Only the symmetric unique elements triangle are strored in RAM.                       ! 
!                                                                                       !
!=======================================================================================!

  subroutine matrix_coo_ram_half( file_number_from_scatci,   &
                                  dimen_matrix,              &
                                  non_zero_ham_elem_half,    &
                                  row_index_array,           & 
                                  column_index_array,        & 
                                  ham_elem_array )

  implicit none

!=======================================================================================!
! A. Declaration part of the code.                                                      !
!=======================================================================================!

!=======================================================================================!
!                                                                                       ! 
! 1. Declaration of the constant variables of the code. Part 1/2.                       !
!                                                                                       !
! a. The phrase "constant variables of the code" meens that                             !
!    since they get a value, this value must not be changed at                          !
!    any step of the rest of the program. In this sense they are constants,             !
!    but not "parameters". This is the only proper way, since                           !
!    all these variables are getting values during the execution                        !
!    of the program, while the "fort.26" is read.                                       !
!                                                                                       ! 
! b. "dimen_matrix" is the dimension of the hamiltonian matrix.                         !
!    The default type in scatci is single integer precision.                            !
!                                                                                       !  
! c. "non_zero_ham_elem_half" is the total numbers of non-zero                          !
!    hamiltonian elements of the lower(upper) triangle only.                            !
!    The kind type is double precision  integer, since this variable                    !
!    is getting large values. This fact does not confilct with the                      !
!    kind type definition of the same variable inside scatci,                           !
!    which is single precision integer                                                  ! 
!    (see note on 'nelm' in scatci_diag). MAL: 15/06/2011                               !
!                                                                                       !
! d. These variables have intent(out) attribute, and they are passed                    !
!    to any subroutine which diagonalizes the hamiltonian matrix.                       !
!                                                                                       !
!=======================================================================================!

  integer, intent(out) :: dimen_matrix
  integer(kind=longint), intent(out) :: non_zero_ham_elem_half

!=======================================================================================!
!                                                                                       !
! 2. Declaration of the arrays for storing in ram the hamiltonian matrix,               !
!    as a lower(upper) triangular sparse matrix.                                        !
!                                                                                       !
! a. "row_index_array" is the single integer precision array                            !
!    for storing the numbers of rows which hold the non-zero                            !
!    hamiltonian matrix elements.                                                       !
!                                                                                       !
! b. "column_index_array" is the single integer precision array                         !
!    for storing the numbers of columns which hold the non-zero                         !
!    hamiltonian matrix elements.                                                       !
!                                                                                       ! 
! c. "ham_elem_array" is the double real precision array                                !
!    for storing the hamiltonian matrix elements.                                       !
!                                                                                       ! 
! d. The integer precision of a. and b. is single. With the current kind type           !
!    we are able to handle a matrix of a dimension 2*(10**9), which is huge.            !
!    So, there is no reason to alter the kind type to any higher kind type.             !
!                                                                                       !
! e. These variables have intent(out) attribute, and they are exported                  !
!    to any subroutine which diagonalize the hamiltonian matrix.                        !
!                                                                                       !
!=======================================================================================!

  integer(kind=shortint), intent(out), allocatable, dimension(:) :: row_index_array
  integer(kind=shortint), intent(out), allocatable, dimension(:) :: column_index_array
  real(kind=wp), intent(out), allocatable, dimension(:) :: ham_elem_array

!=======================================================================================!
!                                                                                       !
! 3. Declaration of the file (holding the Hamiltonian) unit number (usually 26).        !
!                                                                                       !  
! a. This variable has intent(in) attribute and is the only input to                    !
!    the "matrix_coo_ram_half" subroutine.                                              !
!                                                                                       !
!=======================================================================================!

  integer, intent(in) :: file_number_from_scatci

!=======================================================================================!
!                                                                                       !
! 4. The Input Unit.                                                                    !
!                                                                                       !
! a. Definition of the input unit number.                                               !
!                                                                                       ! 
! b. The unit number is local to this subroutine, is arbitrarily choosen,               !
!    and does not conflict with any other unit number definition inside                 !
!    the main scatci program.                                                           !
!                                                                                       !
! c. This variable is local to the subroutine.                                          !
!                                                                                       ! 
!=======================================================================================!
 
   character(len=9) :: filename
   
!=======================================================================================!
!                                                                                       !
! 5. Declaration of some temporary help variables.                                      !
!                                                                                       !
! a. "i1" is a loop variable.                                                           !
!                                                                                       !
! b. "itmp1" is to store temporarly the buffers. There is no any                        !
!    further process of those values. It is a trash variable.                           !
!                                                                                       !
! c. "error_check" is to determine when the reading has                                 !
!    reached the end of the "fort.26" file and to exit the                              !
!    infinite reading do-loop. It is a must variable.                                   !
!                                                                                       !
! d. These variables are local to the subroutine.                                       !
!                                                                                       !
!=======================================================================================!

   integer :: i1 
   integer :: itmp1, error_check

!=======================================================================================!
!                                                                                       !
! 6. Declaration of some temporary help arrays.                                         !
!                                                                                       !
! a. The "ij_tmp" tmp help array is a single precision integer array                    !
!    and must remain as single precision integer array, since the                       !
!    related output of scatci is an integer single precision output.                    !
!                                                                                       ! 
! b. If the related output of scatci, written in the original                           !
!    scatci file "fort.26", is modified and the single type becomes                     !
!    double integer precision, then I have to modify the "ij_tmp"                       !
!    integer type to double precision integer type.                                     !
!                                                                                       ! 
! c. The "emx_tmp" tmp help array is a double precision real array,                     !
!    and must remain as such since the related output of the scatci                     !
!    is a double real precision output.                                                 !
!    There in no reason this to be changed in the future.                               !
!                                                                                       ! 
! d. The purpose of the "ij_tmp" array is to store the indices of                       !
!    rows and columns of hamiltonian matrix elememnts.                                  !
!    The length of this array is the length of the buffer that the                      !
!    scatci uses to write the "fort.26" file. This buffer length                        !
!    is read and used by this program automatically.                                    !
!                                                                                       ! 
! e. The purpose of the "emx_tmp" array is to store the corresponding                   !
!    hamiltonian elements. The length of this array is the length of                    !
!    the buffer that the scatci uses to write the "fort.26" file. This                  !
!    buffer is read and used by this program automatically.                             !
!                                                                                       !
! f. These variables are local to the subroutine.                                       !
!                                                                                       !
!=======================================================================================!
 
  integer, allocatable, dimension(:,:) :: ij_tmp
  real(kind=wp), allocatable, dimension(:) :: emx_tmp 
  
!=======================================================================================!
!                                                                                       ! 
! 7. Declaration of the constant variables of the code. Part 2/2.                       !
!                                                                                       ! 
! a. The phrase "constant variables of the code" meens that                             !
!    since they get a value this value must not be changed at                           !
!    any step of the rest of the program. In this sense they are constant               !
!    but not "parameter". This is the only proper way since                             !
!    all these variables are getting values during the execution                        !
!    of the program, while the "fort.26" is being read.                                 !
!                                                                                       !
! b. "num_buffer" and "num_buffer_tmp" are variables for the                            !
!    first value buffer and for all the rest, correspondingly.                          !
!                                                                                       !
! c. "num_buffer_last" is the dimension of the last buffer.                             !
!                                                                                       !
! d. "num_total_buffers" is the total number of buffers inside the "fort.26" file.      !
!                                                                                       !
! e. These variables are local to the subroutine.                                       !
!                                                                                       !
!=======================================================================================!
 
  integer :: num_buffer
  integer :: num_buffer_tmp
  integer :: num_buffer_last
  integer :: num_total_buffers 

!=======================================================================================!
! B. Execution part of the code.                                                        !
!=======================================================================================!

!=======================================================================================!
!                                                                                       !
! 1. Opening the unit.                                                                  ! 
!                                                                                       !
! a. "matrix_from_scatci" is the unit opened for reading                                !
!    the original "fort.26" scatci produced unformatted file.                           !
!                                                                                       !
!=======================================================================================!

  write(filename,'(i0)') file_number_from_scatci

  open( unit=file_number_from_scatci, &
        file='fort.'//trim(filename), &
        form="unformatted",           &
        status="old")                

!=======================================================================================!
!                                                                                       !
! 2. Reading the dimension of the matrix and the buffer length.                         !
!                                                                                       ! 
! a. First we rewind the "matrix_from_scatci" unit.                                     !
!                                                                                       ! 
! b. "dimen_matrix" is the dimension of the hamiltonian matrix.                         !
!                                                                                       ! 
! c. "num_buffer" is the length of the main buffer.                                     !
!    The length of the last buffer will be computed later.                              !
!                                                                                       !
!=======================================================================================!
 
  rewind(unit=file_number_from_scatci)
  read(unit=file_number_from_scatci) dimen_matrix, num_buffer 

!=======================================================================================!
!                                                                                       !
! 3. Evaluating the length of the last buffer and the number of the total buffers.      !
!                                                                                       !
! a. "num_buffer_last" is the length of the last buffer.                                !
!                                                                                       !
! b. "num_total_buffers" is the variable which its final value                          !
!    is the total number of buffers inside the "fort.26" file.                          !
!                                                                                       !
!=======================================================================================!

  num_total_buffers=0 ! Setting to zero.

  do ! Do-Loop: ALPHA.

    num_total_buffers=num_total_buffers+1

    read(unit=file_number_from_scatci, iostat=error_check) num_buffer_tmp

    if ((num_buffer_tmp /= num_buffer) .and. (error_check == 0)) then ! If: ALPHA.

      num_buffer_last=num_buffer_tmp
 
      exit 

    else if (error_check /= 0) then ! Else-If: ALPHA.

      exit

    end if ! End If: ALPHA.

  end do  ! End of Do-Loop: ALPHA.

!=======================================================================================!
!                                                                                       !
! 4. Evaluating the number of non-zero hamiltonian matrix elements.                     !
!                                                                                       !
! a. "non_zero_ham_elem_half" is the number of total upper(lower)                       !
!    triangle non-zero matrix hamiltonian elements.                                     !
!                                                                                       !
!=======================================================================================!
 
  non_zero_ham_elem_half = (num_total_buffers-1)*num_buffer + num_buffer_last
 
!=======================================================================================!
!                                                                                       ! 
! 5. Allocating the temporary arrays.                                                   !
!                                                                                       !
! a. An important note is that the last "ij_tmp" must                                   !
!    have the common buffer length size. Otherwise the                                  !
!    reading "fort.26" process fails. This is because                                   !
!    in such a way the "fort.26" is built inside scatci.                                !
!                                                                                       !
! b. The following arrays are tmp.                                                      !
!    They are deallocated as soon as they become useless.                               !
!                                                                                       !
!=======================================================================================!

  allocate(ij_tmp(1:2,1:num_buffer)) 
  allocate(emx_tmp(1:num_buffer))    

!=======================================================================================!
!                                                                                       !
! 6. Allocating the arrays for rows, columns and matrix elements.                       !
!                                                                                       !
! a. The following arrays are non-local and have to be passed to                        !
!    any subroutine which diagonalize the matrix. Do not deallocate them here.          !
!    Deallocation must be done inside the calling subroutine.                           !
!                                                                                       !
!========================================================================================

  allocate(row_index_array(1:non_zero_ham_elem_half))
  allocate(column_index_array(1:non_zero_ham_elem_half))
  allocate(ham_elem_array(1:non_zero_ham_elem_half))

!=======================================================================================!
!                                                                                       !
! 7. Rewinding the "matrix_from_scatci" unit.                                           !
!                                                                                       ! 
! a. The purpose is to rewind the unit and to skip the first line.                      !
!                                                                                       !
!=======================================================================================! 

  rewind(unit=file_number_from_scatci) 
  read(unit=file_number_from_scatci) 

!=======================================================================================!
!                                                                                       !
! 8. Reading all the buffers except the last and writing the infos.                     !
!                                                                                       !
! a. The purpose is to read all the elements up to the previous                         !
!    to the last buffer of "fort.26" file, and on the fly, for                          !
!    each buffer, to store its infos to the non-local arrays                            !
!    "row_index_array", "column_index_array", "ham_elem_array".                         !
!                                                                                       !
! b. The reading and writing of the last buffer is the next step.                       !
!                                                                                       !
!=======================================================================================!
 
 do i1=1, num_total_buffers-1 ! Do-Loop: BETA.

    read(unit=file_number_from_scatci) itmp1, ij_tmp(:,:), emx_tmp(:)
   
    row_index_array(1+(i1-1)*num_buffer : i1*num_buffer) =     &
                                   
      ij_tmp(1,1:num_buffer)
    
    column_index_array(1+(i1-1)*num_buffer : i1*num_buffer) =  &

      ij_tmp(2,1:num_buffer)
    
    ham_elem_array(1+(i1-1)*num_buffer : i1*num_buffer) =      &
         
      emx_tmp(1:num_buffer)

 end do ! End of Do-Loop: BETA.

!=======================================================================================!
!                                                                                       !
! 9. Reading the last buffer and writting the infos.                                    !
!                                                                                       !
! a. The purpose is to read the last buffer and to put the data                         !
!    to the arrays "row_index_array", "column_index_array" and                          !
!    "ham_elem_array".                                                                  !
!                                                                                       !
!=======================================================================================!

  read(unit=file_number_from_scatci) itmp1, ij_tmp(:,:), emx_tmp(:)

  row_index_array(num_buffer*(num_total_buffers-1)+1 : non_zero_ham_elem_half) =     &
                             
    ij_tmp(1,1:num_buffer_last)

  column_index_array(num_buffer*(num_total_buffers-1)+1 : non_zero_ham_elem_half) =  &

    ij_tmp(2,1:num_buffer_last)

  ham_elem_array(num_buffer*(num_total_buffers-1)+1 : non_zero_ham_elem_half) =      & 
  
    emx_tmp(1:num_buffer_last)

!=======================================================================================!
!                                                                                       !
! 10. Deallocating the tmp help arrays.                                                 !
!                                                                                       !
! a. The "ij_tmp(:,:)" and "emx_tmp(:)" arrays are no more needed.                      !
!                                                                                       !
!=======================================================================================!

  deallocate(ij_tmp)  
  deallocate(emx_tmp) 

!=======================================================================================!
!                                                                                       !
! 11. Closing the unit and keeping the file.                                            !
!                                                                                       !
!=======================================================================================!

  close(unit=file_number_from_scatci, status="keep")

  end subroutine matrix_coo_ram_half

! FINI.

end module m_2_matrix_sparse
