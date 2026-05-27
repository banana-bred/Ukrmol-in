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

!> \brief   Matrix element module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module MatrixElement_module

    use precisn,            only: wp
    use const_gbl,          only: stdout
    use MatrixCache_module, only: MatrixCache
    use BaseMatrix_module,  only: BaseMatrix

    implicit none

    public MatrixElementVector

    private

    !> \brief This handles the matrix elements and also expands the vector size if we have reached max capacity
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    type, extends(BaseMatrix) :: MatrixElementVector
        type(MatrixCache) :: matrix_cache
    contains
        procedure, public :: sort                   => sort_matelem
       !procedure, public :: prune_threshold
        procedure, public :: initialize_struct_self => initialize_struct_standard
        procedure, public :: construct_self         => construct_mat_standard
        procedure, public :: insert_matelem_self    => insert_matelem_standard
        procedure, public :: get_matelem_self       => get_matelem_standard
        procedure, public :: clear_self             => clear_standard
        procedure, public :: expand_capacity        => expand_standard
        procedure, public :: destroy_self           => destroy_standard
       !procedure, public :: count_occurance
       !procedure, public :: construct
    end type MatrixElementVector

contains

    !> \brief   Check that type is constructed
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Checks if the type is constructed and hard aborts when not.
    !>
    subroutine check_constructed (this)
        class(MatrixElementVector) :: this

        if (.not. this % constructed) then
            write (stdout, "('Vector::constructed - Vector is not constructed')")
            stop "Vector::constructed - Vector is not constructed"
        end if

    end subroutine check_constructed


    !> \brief   Construct the type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine construct_mat_standard (this)
        class(MatrixElementVector) :: this

        print *, 'MatrixElement::construct_mat_standard'

        call this % matrix_cache % construct
        call this % matrix_cache % clear

    end subroutine construct_mat_standard


    !> \brief   Initialize the type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine initialize_struct_standard (this, matrix_size, matrix_type, block_size)
        class(MatrixElementVector) :: this
        integer, intent(in)        :: matrix_size, matrix_type, block_size
        integer                    :: ido

        this % initialized = .true.
       !this % matrix_dimension = matrix_size

    end subroutine initialize_struct_standard


    !> \brief   Insert a new element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine insert_matelem_standard (this, i, j, coefficient, class, thresh)
        class(MatrixElementVector) :: this
        integer,  intent(in)       :: i, j, class
        real(wp), intent(in)       :: coefficient, thresh

        !this%matrix_size = max(this%matrix_size,i,j)
        if (abs(coefficient) < thresh) return

        call this % matrix_cache % insert_into_cache(i, j, coefficient)

        this % n = this % matrix_cache % n

    end subroutine insert_matelem_standard


    !> \brief   Retrieve element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine get_matelem_standard (this, idx, i, j, coeff)
        class(MatrixElementVector) :: this
        integer,  intent(in)       :: idx
        integer,  intent(out)      :: i, j
        real(wp), intent(out)      :: coeff

        call this % matrix_cache % get_from_cache(idx, i, j, coeff)

    end subroutine get_matelem_standard


    !> \brief   Retrieve element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine expand_standard (this, capacity)
        class(MatrixElementVector) :: this
        integer, intent(in)        :: capacity

        call this % matrix_cache % expand_capacity(capacity)

        write (stdout, "('Current Matrix array size is ',i4)") this % matrix_cache % max_capacity

    end subroutine expand_standard


    !> \brief   Clear elements
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine clear_standard (this)
        class(MatrixElementVector) :: this
        integer                    :: ido

        !do ido=1,this%max_capacity
        !    this%matrix_elements(ido)%i = -1
        !    this%matrix_elements(ido)%j = -1
        !    this%matrix_elements(ido)%coefficient = 0.0_wp
        !enddo
        !this%matrix_size = 0

        call this % matrix_cache % clear

    end subroutine clear_standard


    !> \brief   Sort elements
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine sort_matelem (this)
        class(MatrixElementVector) :: this

        call this % matrix_cache % sort

    end subroutine sort_matelem


    !> \brief   Destroy elements
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine destroy_standard (this)
        class(MatrixElementVector) :: this

        call this % matrix_cache % destroy

    end subroutine destroy_standard

end module MatrixElement_module
