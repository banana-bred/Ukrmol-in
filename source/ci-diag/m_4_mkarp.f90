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

module m_4_mkarp

implicit none

integer :: default_arpack_max_iter = 1000000

contains

!=====================================================================!
!=====================================================================!  
!                                                                     !
! 1a. Author: Pavlos G. Galiatsatos.                                  !
! 1b. Modifications, Corrections, Improvements: Michael Lysaght       !
!                                                                     !
! 2. Date: 2011/06/30.                                                !
!                                                                     !
! 3. Subroutine:  mkarp.                                              !
!                                                                     !
! 4. Details:                                                         !
!                                                                     !
!    a. I deleted the previous version.                               !
!    b. I putted a new version.                                       !
!    c. The interface has been simplified.                            !
!                                                                     !
! 5. Purpose:                                                         !
!                                                                     !
!    a. Reading the "fort.26" file.                                   !
!    b. Building the hamiltonian matrix in coordinate format.         !
!    c. Storing in ram ONLY the half elements of the                  !
!       hamiltonian matrix in coordinate format.                      !
!    d. Use ARPACK to diagonalize.                                    !
!    e. Not "save-state" implemented yet, in ARPACK.                  !
!                                                                     !
!=====================================================================!
!=====================================================================!
 
  subroutine mkarp(nhdim,nstat,maxiter,crite,file_name_from_scatci) 
                                                                  
    use precisn 
    use m_2_matrix_sparse
    use m_3_arpack_diagonalizers

    implicit none 
                                                                        
!====================================================================!             
! 1. Declaration of the interface arguments.                         !             
!====================================================================!                                                                                    
    integer, intent(out)                     :: nhdim
    integer, intent(inout)                   :: nstat
    integer, intent(inout)                   :: maxiter
    real(wp),intent(inout)                   :: crite
    integer, intent(in)                      :: file_name_from_scatci 

    integer(kind=longint) :: nelm  ! mal 10/06/2011 : as defined in pgg version

!MAL 10/06/2011, Note on 'nelm': 'nelm' is no longer included in the 
!argument list of subroutine 'mkarp', due to the fact that it is calculated
!upon reading in the fort.26 file. 'nelm' is the value of for the number of
!non-zero elements in the sparse N+1 Hamiltonian stored in the fort.26 file.
!Because, nelm can be a "large" integer value, it is declared here as
!integer(kind=longint). This means that, irrespective of whether one compiles
!with the -i4 or -i8 flag, 'nelm' will be able to store 64-bit integer values.
!The advantage of declaring 'nelm' in this way is that for certain problems
!one may need to only compile with the -i4 flag (i.e., the construction of
!the N+1 Hamiltonian does not require the -i8 flag), but can still deal
!with a large number of non-zero N+1 Hamiltonian elements. 

    integer                               :: ncv_coeff
    character(len=2)                      :: which

!=====================================================================!
! 2. Input/Output for "arpack_coo_half".                              !
!=====================================================================!
                                                                        

    integer(kind=shortint), allocatable, dimension(:) :: row_index_array       
    integer(kind=shortint), allocatable, dimension(:) :: column_index_array
    real(kind=wp), allocatable, dimension(:)    :: ham_elem_array 
    real(kind=wp), allocatable, dimension(:,:)  :: vectors_array
    real(kind=wp), allocatable, dimension(:,:)  :: values_array  

!======================================================================!
!  3. Local variables.                                                 !
!======================================================================!

    integer :: i1, i2

!======================================================================!             
! 4. Executable part.                                                  !             
!======================================================================!
     
   write(*,*) "======================================================"
   write(*,*) " inside subroutine --> 'mkarp' --> begins. "
   write(*,*) " "

    if (maxiter == -1) then                        !if clause added by mal 13/06/2011
      maxiter = default_arpack_max_iter            !Unlike pgg version, maxiter is now
    else if (maxiter < 300) then                   !passed to mkarp from namelist
      maxiter = 300                                !maxiter should be chosen to be at 
      print*, 'maxiter has been set to 300'        !at least 300 for arpack.
    end if

    call matrix_coo_ram_half( file_name_from_scatci, & 
                              nhdim,                 &
                              nelm,                  &
                              row_index_array,       &
                              column_index_array,    &
                              ham_elem_array )
  
    write(*,*) " 1. nstat --> ", nstat 
    write(*,*) " 2. nhdim --> ", nhdim 
    write(*,*) " 3. nelm  --> ", nelm 
    write(*,*) " 4. row_index_array(10)    --> ", row_index_array(10)
    write(*,*) " 5. column_index_array(10) --> ", column_index_array(10)
    write(*,*) " 6. ham_elem_array(10)     --> ", ham_elem_array(10)
    
    which = 'SA'
    ncv_coeff = 2     
     
    call arpack_coo_half_auto( row_index_array,                  &
                               column_index_array,               &
                               ham_elem_array,                   &
                               nhdim,                            &
                               nelm,                             &
                               which,                            &
                               ncv_coeff,                        &
                               nstat,                            &
                               maxiter,                          &
                               crite,                            &
                               vectors_array,                    &
                               values_array )

! opening the unit to write the eigensystem.

    open( unit=11,                                    &
          file="___tmp_eigensys",                     & 
          status="replace",                           &
          action="write",                             &
          form="unformatted" )                    

! writing to hard disk the eigenvalues.

    do i1 = 1, nstat
      write(unit=11) values_array(i1,1)
    end do

! writing to hard disk the eigenvwectors. 

    do i1 = 1, nstat
      do i2 = 1, nhdim
        write(unit=11) vectors_array(i2,i1)
      end do 
    end do
    
! closing the unit.

    close( unit=11, status="keep" )

! deallocating the ram space.

    deallocate( row_index_array ) 
    deallocate( column_index_array ) 
    deallocate( ham_elem_array ) 
    deallocate( values_array ) 
    deallocate( vectors_array )  
    
    write(*,*) " "  
    write(*,*) " inside subroutine --> 'mkarp' --> ends. "
    write(*,*) "====================================================="
    write(*,*) " "

  end subroutine mkarp

! FINI.

end module m_4_mkarp 

