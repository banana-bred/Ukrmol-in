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

!> \brief   Determinant binary search tree
!> \authors J Benda
!> \date    2018 - 2019
!>
!> This module contains a special data structure, \ref det_tree, used (e.g.) by CONGEN to carry out quick look-ups in the
!> storage of determinants. It is based on the type "bstree" from libAMOR.
!>
module congen_bstree

    use containers, only: bstree

    implicit none

    !> \brief   Determinant binary search tree
    !> \authors J Benda
    !> \date    2019
    !>
    !> This is an extension of the plain integer binary search tree to one that operates on the determinant
    !> storage. The relation used in the ordering is the lexicographical order of the determinants as provided
    !> by the overriden subroutine "compare".
    !>
    type, extends(bstree) :: det_tree
        integer, pointer :: det(:) => null()
        integer, pointer :: ndo(:) => null()
        integer          :: nelt = 0
    contains
        procedure, public :: init       => init_det_tree
        procedure, public :: compare    => compare_determinats
        procedure, public :: locate_det => locate_in_det_tree
        procedure, public :: dtwrite    => write_determinant
        final :: final_det_tree
    end type det_tree

contains

    !> \brief   Initialize determinant tree
    !> \authors J Benda
    !> \date    2019
    !>
    !> Stores pointer to the determinant storage, so that it can be used when comparing determinant indices
    !> in the binary search tree subroutines.
    !>
    !> \param this    Tree object to initialize
    !> \param ndo     Determinant storage
    !> \param nelt    Size of a determinant
    !>
    subroutine init_det_tree (this, ndo, nelt)

        class(det_tree),  intent(inout) :: this
        integer, pointer, intent(in)    :: ndo(:)
        integer,          intent(in)    :: nelt

        this % ndo  => ndo
        this % nelt = nelt

        allocate (this % det(nelt))

    end subroutine init_det_tree


    !> \brief   Finalize determinant tree
    !> \authors J Benda
    !> \date    2019
    !>
    !> Releases all allocated memory.
    !>
    !> \param this Tree object to finalize
    !>
    subroutine final_det_tree (this)

        type(det_tree), intent(inout) :: this

        if (associated(this % det)) then
            deallocate (this % det)
        end if

    end subroutine final_det_tree


    !> \brief   Find a specific determinant in the storage
    !> \authors J Benda
    !> \date    2019
    !>
    !> Return index in the determinant storage, where the given determinant is located. Return -1 when there is no such
    !> determinant.
    !>
    !> \param this  Tree object to search
    !> \param det   Determinant to find (length nelt)
    !>
    integer function locate_in_det_tree (this, det) result (res)

        class(det_tree),                 intent(inout) :: this
        integer, dimension(this % nelt), intent(in)    :: det

        this % det = det

        res = this % locate(-1)  ! negative id signalizes unary comparison, see "compare_determinats"

    end function locate_in_det_tree


    !> \brief   Write determinat
    !> \authors J Benda
    !> \date    2019
    !>
    !> Used in debuggind output of the whole tree.
    !>
    !> \param this  Binary search tree
    !> \param lu    Unit for writing
    !> \param id    Tree node id
    !>
    subroutine write_determinant (this, lu, id)

        class(det_tree), intent(in) :: this
        integer,         intent(in) :: lu, id

        integer :: i

        write (lu, '("[")', advance = 'no')

        do i = 1, this % nelt
            write (lu, '(1x,I0)', advance = 'no') this % ndo((id - 1) * this % nelt + i)
        end do

        write (lu, '(1x,"]")')

    end subroutine write_determinant


    !> \brief   Lexicographically compare two determinats
    !> \authors J Benda
    !> \date    2018 - 2019
    !>
    !> Binary tree needs a notion of order of its elements to be able to work with the data
    !> so efficiently. In case of arrays of numbers, the typical order is the "lexicographical"
    !> order. The array A is less than B if for the first position I for which A(I) /= B(I)
    !> holds that A(I) < B(I).
    !>
    !> If any of the given indices is negative, the reference determinant stored within the
    !> predicate type is used instead.
    !>
    !> \param this  Determinant tree object to use
    !> \param i     Index of the first determinant
    !> \param j     Index of the second determinant
    !> \param data  Additional payload required by base class interface. Not used in congen.
    !>
    !> \return -1 if the first determinant is lexicographically less,
    !>          0 if they are equal,
    !>         +1 if the first is greater.
    !>
    integer function compare_determinats (this, i, j, data) result (verdict)

        use iso_c_binding, only: c_ptr

        class(det_tree),       intent(in) :: this
        integer,               intent(in) :: i, j
        type(c_ptr), optional, intent(in) :: data

        integer, pointer :: d1(:), d2(:)
        integer :: k

        ! get first determinant view
        if (i >= 0) then
            d1 => this % ndo((i - 1) * this % nelt + 1 : i * this % nelt)
        else
            d1 => this % det(1:this % nelt)
        end if

        ! get second determinant view
        if (j >= 0) then
            d2 => this % ndo((j - 1) * this % nelt + 1 : j * this % nelt)
        else
            d2 => this % det(1:this % nelt)
        end if

        ! compare determinants element by element
        do k = 1, this % nelt

            ! is first less?
            if (d1(k) < d2(k)) then
                verdict = -1
                return
            end if

            ! is first greater?
            if (d1(k) > d2(k)) then
                verdict =  1
                return
            end if

        end do

        ! they are equal
        verdict = 0

    end function compare_determinats

end module congen_bstree
