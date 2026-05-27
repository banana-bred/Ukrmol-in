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

!> \brief   Base matrix module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module BaseMatrix_module

    use const_gbl,              only: stdout
    use consts_mpi_ci,          only: MATRIX_INDEX_FORTRAN, DEFAULT_MATELEM_THRESHOLD
    use Options_module,         only: Options
    use precisn,                only: wp

    implicit none

    public BaseMatrix

    private

    !> \brief   Base matrix type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This handles the matrix elements and also expands the vector size if we have reached max capacity.
    !>
    type, abstract :: BaseMatrix
        real(wp)              :: threshold
        integer               :: matrix_type
        integer               :: matrix_indexing = MATRIX_INDEX_FORTRAN
        integer               :: matrix_mode
        integer               :: excluded_rowcol = 2000000000
        !>Number of elements in the array of both the integral and coefficients
        logical               :: constructed = .false.
        integer               :: matrix_dimension
        integer               :: n = 0
        integer               :: block_size
        logical               :: remove_row_column = .false.
        real(wp), allocatable :: diagonal(:)
        logical               :: initialized = .false.
    contains
        procedure, public :: initialize_matrix_structure
        procedure, public :: construct
        procedure, public :: insert_matrix_element
        procedure, public :: get_matrix_element
        procedure, public :: exclude_row_column
        procedure, public :: is_empty
        procedure, public :: get_size
        procedure, public :: get_matrix_size

        procedure, public :: clear
        procedure, public :: destroy

        procedure, public :: update_continuum
        procedure, public :: set_options
        procedure, public :: update_pure_L2
        procedure, public :: finalize_matrix
        procedure, public :: store_diagonal
        procedure, public :: initialize_struct_self

        procedure(generic_construct), deferred :: construct_self
        procedure(generic_insert),    deferred :: insert_matelem_self
        procedure(generic_get),       deferred :: get_matelem_self
        procedure(generic_clear),     deferred :: clear_self
        procedure(generic_destroy),   deferred :: destroy_self

        procedure, public  :: expand_capacity
        procedure, public  :: print => print_mat

        procedure, private :: check_bounds
        procedure, private :: check_constructed
    end type

    abstract interface
        subroutine generic_construct (this)
            import :: BaseMatrix
            class(BaseMatrix) :: this
        end subroutine generic_construct
    end interface

    abstract interface
        subroutine generic_insert(this, i, j, coefficient, class, thresh)
           use precisn, only: wp
            import :: BaseMatrix
            class(BaseMatrix)    :: this
            integer,  intent(in) :: i, j
            real(wp), intent(in) :: coefficient
            integer,  intent(in) :: class
            real(wp), intent(in) :: thresh
        end subroutine generic_insert
    end interface

    abstract interface
        subroutine generic_get (this, idx, i, j, coeff)
            use precisn, only: wp
            import :: BaseMatrix
            class(BaseMatrix)     :: this
            integer,  intent(in)  :: idx
            integer,  intent(out) :: i, j
            real(wp), intent(out) :: coeff
        end subroutine generic_get
    end interface

    abstract interface
        subroutine generic_clear (this)
            import :: BaseMatrix
            class(BaseMatrix) :: this
        end subroutine generic_clear
    end interface

    abstract interface
       subroutine generic_destroy (this)
        import :: BaseMatrix
        class(BaseMatrix) :: this
        end subroutine generic_destroy
    end interface

contains

    !> \brief   Check that matrix is constructed
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Check that matrix is constructed; hard stop if not.
    !>
    subroutine check_constructed (this)
        class(BaseMatrix) :: this

        if (.not. this % constructed) then
            write (stdout, "('Vector::constructed - Vector is not constructed')")
            stop "Vector::constructed - Vector is not constructed"
        end if

    end subroutine check_constructed


    !> \brief   Construct the matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine construct (this, mat_mode, threshold)
        class(BaseMatrix)              :: this
        integer,            intent(in) :: mat_mode
        real(wp), optional, intent(in) :: threshold

        integer :: err

        if (present(threshold)) then
            this % threshold = threshold
        else
            this % threshold = DEFAULT_MATELEM_THRESHOLD
        end if

        this % matrix_mode = mat_mode

        call this % construct_self
        this % remove_row_column = .false.
        this % constructed = .true.

    end subroutine construct


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine exclude_row_column (this, row_column)
        class(BaseMatrix)   :: this
        integer, intent(in) :: row_column

        this % excluded_rowcol = row_column
        if (this % excluded_rowcol > 0) then
            this % remove_row_column = .true.
        end if

    end subroutine exclude_row_column


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine initialize_matrix_structure (this, matrix_size, matrix_type, block_size)
        class(BaseMatrix)   :: this
        integer, intent(in) :: matrix_size, matrix_type, block_size
        integer             :: ido, ierr

        this % initialized = .true.
        this % matrix_dimension = matrix_size
        this % matrix_type = matrix_type
        this % block_size  = block_size

        if (this % remove_row_column .and. this % excluded_rowcol <= this % matrix_dimension) then
            this % matrix_dimension = this % matrix_dimension - 1
        end if

        if (allocated(this % diagonal)) deallocate(this % diagonal)
        write (stdout, "('Allocating dimension ',i12)") this % matrix_dimension

        allocate(this % diagonal(this % matrix_dimension), stat = ierr)
        this % diagonal = 0
        if (ierr /= 0) then
            stop "ERROR"
        end if

        call this % initialize_struct_self(this % matrix_dimension, this % matrix_type, this % block_size)

    end subroutine initialize_matrix_structure


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Usually does nothing but some matrix formats may need more parameters in order to function.
    !>
    subroutine set_options (this, option_val)
        class(BaseMatrix)          :: this
        class(Options), intent(in) :: option_val
    end subroutine set_options


    !> \brief   Initialize the type
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Do nothing
    !>
    subroutine initialize_struct_self (this, matrix_size, matrix_type, block_size)
        class(BaseMatrix)   :: this
        integer, intent(in) :: matrix_size, matrix_type, block_size
    end subroutine initialize_struct_self


    !> \brief   Set matrix element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine insert_matrix_element (this, i, j, coefficient, class_, thresh_)
        class(BaseMatrix)              :: this
        integer,  intent(in)           :: i, j
        real(wp), intent(in)           :: coefficient
        integer,  intent(in), optional :: class_
        real(wp), intent(in), optional :: thresh_
        real(wp) :: thresh
        integer  :: class, i_, j_

        !this%matrix_size = max(this%matrix_size,i,j)

        if (.not. this % initialized) then
            stop "Matrix structure not initialized"
        end if

        if (present(thresh_)) then
            thresh = thresh_
        else
            thresh = this % threshold
        endif

        if (present(class_)) then
            class = class_
        else
            class = -1
        end if

        !Make sure we are in the lower triangular
        if (i >= j) then
            i_ = i
            j_ = j
        else
            i_ = j
            j_ = i
        end if

        !Exclude if needed
        if (this % remove_row_column) then
            if (i_ == this % excluded_rowcol) return
            if (j_ == this % excluded_rowcol) return
            if (i_ > this % excluded_rowcol) i_ = i_ - 1
            if (j_ > this % excluded_rowcol) j_ = j_ - 1
        end if
        if (i == j) call this % store_diagonal(i, coefficient)
        if (abs(coefficient) < thresh) return

        call this % insert_matelem_self(i_, j_, coefficient, class, thresh)

    end subroutine insert_matrix_element


    !> \brief   Set diagonal element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine store_diagonal (this, i, coeff)
        class(BaseMatrix)    :: this
        integer,  intent(in) :: i
        real(wp), intent(in) :: coeff

        !if(myrank /= master) return
        this % diagonal(i) = coeff

    end subroutine store_diagonal


    !> \brief   Retrieve matrix element
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine get_matrix_element (this, idx, i, j, coeff)
        class(BaseMatrix)     :: this
        integer,  intent(in)  :: idx
        integer,  intent(out) :: i, j
        real(wp), intent(out) :: coeff

        if (this % check_bounds(idx)) then
            call this % get_matelem_self(idx, i, j, coeff)
        end if

    end subroutine get_matrix_element


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine update_continuum (this, force_update)
        class(BaseMatrix)   :: this
        logical, intent(in) :: force_update
    end subroutine update_continuum


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine update_pure_L2 (this, force_update, count_)
        class(BaseMatrix)   :: this
        logical, intent(in) :: force_update
        integer, optional   :: count_
    end subroutine update_pure_L2


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine finalize_matrix (this)
        class(BaseMatrix) :: this
    end subroutine finalize_matrix


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine expand_capacity (this, capacity)
        class(BaseMatrix)   :: this
        integer, intent(in) :: capacity
    end subroutine expand_capacity


    !> \brief   ?
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    logical function check_bounds (this, i)
        class(BaseMatrix)   :: this
        integer, intent(in) :: i

        if (i <= 0 .or. i > this % n) then
            write (stdout, "('Vector::check_bounds - Out of Bounds access')")
            stop "Vector::check_bounds - Out of Bounds access"
            check_bounds = .false.
        else
            check_bounds = .true.
        end if

    end function check_bounds


    !> \brief   Determine if matrix is empty
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Check if there are no elements in the mastrix.
    !>
    logical function is_empty(this)
        class(BaseMatrix) :: this

        is_empty = (this % n == 0)

    end function is_empty


    !> \brief   Clear matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine clear (this)
        class(BaseMatrix) :: this
        integer           :: ido

        !do ido=1,this%max_capacity
        !    this%matrix_elements(ido)%i = -1
        !    this%matrix_elements(ido)%j = -1
        !    this%matrix_elements(ido)%coefficient = 0.0_wp
        !enddo
        !this%matrix_size = 0

        this % n = 0
        call this % clear_self
        this % initialized = .false.

    end subroutine clear


    !> \brief   Get matrix size (rank)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer function get_size (this)
        class(BaseMatrix) :: this

        get_size = this % n

    end function get_size


    !> \brief   Get matrix size (number of elements)
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    integer function get_matrix_size(this)
        class(BaseMatrix) :: this

        if (.not. this % initialized) then
            stop "Matrix structure not initialized"
        end if

        get_matrix_size = this % matrix_dimension

    end function get_matrix_size


    !> \brief   Print matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine print_mat (this)
        class(BaseMatrix) :: this
        integer           :: i, j, ido
        real(wp)          :: coeff

        do ido = 1, this % n
            call this % get_matrix_element(ido, i, j, coeff)
            write (stdout, "(i8,i8,D16.8)") i, j, coeff
        end do

    end subroutine print_mat


    !> \brief   Destroy matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine destroy (this)
        class(BaseMatrix) :: this
        integer           :: num_arrays,ido

        if (allocated(this % diagonal)) deallocate(this % diagonal)

        call this % destroy_self

        this % constructed = .false.

    end subroutine destroy


    !> \brief   Update matrix
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine update_matrix (this)
        class(BaseMatrix) :: this
    end subroutine update_matrix

end module BaseMatrix_module
