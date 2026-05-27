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

!> \brief   Elementary algorithms for use in the UKRmol+ suite
!> \authors J Benda
!> \date    2019
!>
!> This module provides a simple algorithm library for reuse in UKRmol+ programs.
!>
module algorithms

    implicit none

    private

    interface findloc
        module procedure findloc_i_1d
    end interface findloc

    interface insertion_sort
        module procedure insertion_sort_i
    end interface insertion_sort

    public findloc, indexx, indexx_stable, insertion_sort, rank

contains

    !> \brief   Linear search in integer array
    !> \authors J Benda
    !> \date    2019
    !>
    !> A convenience subroutine defined in Fortran 2008. Can be safely erased once Fortran 2008 is fully embraced.
    !>
    !> \param[in] array  One-dimensional array to search.
    !> \param[in] value  Value to find.
    !> \param[in] dimen  Dummy integer.
    !>
    !> \return Zero if not found or index in array otherwise.
    !>
    integer function findloc_i_1d (array, value, dimen) result (position)

        integer, intent(in) :: array(:), value, dimen
        integer             :: i

        position = 0
        do i = 1, size(array)
            if (array(i) == value) then
                position = i
                return
            end if
        end do

    end function findloc_i_1d


    !> \brief   In-place insertion sort of integer array
    !> \authors J Benda
    !> \date    2019
    !>
    !> Insertion sort iterates over elements from left to right and places them on the right place within the left,
    !> already sorted part of the array. Produces a non-descending array.
    !>
    subroutine insertion_sort_i (array)

        integer, intent(inout) :: array(:)
        integer                :: i, j, m, n

        n = size(array)
        do i = 2, n
            do j = i - 1, 1, -1
                if (array(j) > array(j + 1)) then
                    m = array(j)
                    array(j) = array(j + 1)
                    array(j + 1) = m
                end if
            end do
        end do

    end subroutine insertion_sort_i


    !> \brief   Create a sort index for a real array
    !> \authors C J Gillan
    !>
    !> Takes an array and produces a set of indices such that ARRIN(INDX(J)) is in ascending order J=1,2,..,N
    !>
    !> \param[in]  n      Number of elements in the array to be ordered.
    !> \param[in]  arrin  Real array which is to placed in ascending order.
    !> \param[out] indx   A set of indices for ascending indices.
    !>
    !> This routine is taken from the book Numerical Receipes by Press, Flannery, Teukolsky and Vetterling Chapter 8 p. 233.
    !> ISBN 0-521-30811-9 pub. Cambridge University Press (1986) QA297.N866
    !>
    !> This routines has been adapted by Charles J Gillan for use in the R-matrix codes.
    !>
    subroutine indexx (n, arrin, indx)

        use precisn, only: wp

        real(wp), parameter :: VSMALL = 1.0e-20_wp

        integer,  intent(in)  :: n
        integer,  intent(out) :: indx(n)
        real(wp), intent(in)  :: arrin(n)

        integer  :: i, indxt, ir, j, l
        real(wp) :: q

        ! Initialize the index array with consecutive integers
        do j = 1, n
            indx(j) = j
        end do
        if (n <= 1) return

        l = n/2 + 1
        ir = n

        !> From here on the algorithm is HEAPSORT wit indirect addressing
        !> through INDX in all references to ARRIN

    10  continue

        if (l > 1) then
            l = l - 1
            indxt = indx(l)
            q = arrin(indxt)
        else
            indxt = indx(ir)
            q = arrin(indxt)
            indx(ir) = indx(1)
            ir = ir - 1
            if (ir == 1) then
                indx(1) = indxt
                go to 800
            end if
        end if

        i = l
        j = l + l

    20  continue

        if (j <= ir) then
            if (j < ir) then
                if (arrin(indx(j)) < arrin(indx(j + 1)) + VSMALL) j = j + 1
            end if
            if (q < arrin(indx(j)) + VSMALL) then
                indx(i) = indx(j)
                i = j
                j = j + j
            else
                j = ir + 1
            end if
            go to 20
        end if
        indx(i) = indxt

        ! loop back to process another element
        go to 10

    800 continue

    end subroutine indexx


    !> \brief   Create a sort index for a real array using merge sort
    !> \authors J Benda
    !> \date    2025
    !>
    !> Takes an array and produces a set of indices such that ARRIN(INDX(J)) is in non-descending order for J=1,2,..,N.
    !> Unlike \ref indexx, this subroutine makes sure that the relative order of equal elements is not changed.
    !> This is guaranteed by using a stable sort algorithm (merge sort).
    !>
    !> \param[in]  n      Number of elements in the array to be ordered.
    !> \param[in]  arrin  Real array which is to placed in non-descending order.
    !> \param[out] indx   A set of indices for non-descending indices.
    !>
    subroutine indexx_stable (n, arrin, indx)

        use precisn, only: wp

        integer,  intent(in)  :: n
        integer,  intent(out) :: indx(n)
        real(wp), intent(in)  :: arrin(n)

        integer,  allocatable :: indx0(:)

        integer :: i, j, jmax, k, kmax, m, p

        allocate (indx0(n))

        ! Initialize the index array with consecutive integers
        do i = 1, n
            indx(i) = i
        end do

        ! Perform merge sort
        m = 1
        do while (m < n)
            indx0 = indx
            p = 1
            do while (p + m <= n)
                i = p
                j = p
                k = p + m
                jmax = j + m - 1
                kmax = min(n, k + m - 1)
                do while (j <= jmax .or. k <= kmax)
                    if (j <= jmax .and. k <= kmax) then
                        if (arrin(indx0(j)) <= arrin(indx0(k))) then
                            indx(i) = indx0(j)
                            j = j + 1
                        else
                            indx(i) = indx0(k)
                            k = k + 1
                        end if
                    else if (j <= jmax) then
                        indx(i) = indx0(j)
                        j = j + 1
                    else
                        indx(i) = indx0(k)
                        k = k + 1
                    end if
                    i = i + 1
                end do
                p = p + 2*m
            end do
            m = 2*m
        end do

    end subroutine indexx_stable


    !> \brief   Construct inverse permutation
    !> \authors C J Gillan
    !>
    !> Takes an array of indices as output by the routine \ref indexx and returns a table of ranks.
    !>
    !> \param[in]  n      Number of elements in the array to be ordered.
    !> \param[in]  indx   Array of indices created by routine \ref indexx.
    !> \param[out] irank  A set ranks corresponding to the indices.
    !>
    !> This routine is taken from the book Numerical Receipes by Press, Flannery, Teukolsky and Vetterling Chapter 8 p. 233.
    !> ISBN 0-521-30811-9 pub. Cambridge University Press (1986) QA297.N866
    !>
    !> This routines has been adapted by Charles J Gillan for use in the R-matrix codes.
    !>
    subroutine rank (n, indx, irank)

        integer, intent(in)  :: n
        integer, intent(in)  :: indx(n)
        integer, intent(out) :: irank(n)

        integer :: j

        ! Convert the indices into ranks
        do j = 1, n
            irank(indx(j)) = j
        end do

    end subroutine rank

end module algorithms
