! Copyright 2025
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
!> \brief   Unit test for `sort_gbl::multiway_merge_float_int`
!> \authors J Benda
!> \date    2025
!>
!> Prepares an array with presorted chunks, applies the m-way merge routine, and tests the resulting order.
!> Checks both in-memory and scratch-backed versions. Uses the sequences from the sieve of Eratosthenes, which
!> have very non-uniform lengths.
!>
program test_multiway_merge_float_int

    use precisn_gbl, only: cfp
    use sort_gbl,    only: multiway_merge_float_int

    implicit none

    integer,   parameter   :: n = 1000
    integer,   allocatable :: ints0(:), ints(:), sieve(:), sizes(:)
    real(cfp), allocatable :: reals(:)

    integer :: i, j, m, k

    allocate (ints0(n), sieve(n), sizes(n))

    do i = 1, n
        sieve(i) = i
    end do

    ! Eratosthenes sieve sequences:
    !   1) 1
    !   2) multiples of 2
    !   3) multiples of 3 (but not of 2)
    !   4) multiples of 5 (but not of 2 and 3)
    !   5) ... etc ...

    i = 1           ! position in the output array `ints0`
    ints0(1) = 1    ! manually added the first one-element sequence [1]
    sizes(1) = 1    ! length of the first sequence
    m = 1           ! number of sequences
    do j = 2, n
        do k = 2, n
            if (sieve(k) /= 0 .and. mod(sieve(k), j) == 0) then
                i = i + 1
                ints0(i) = sieve(k)
                sieve(k) = 0
                sizes(m + 1) = sizes(m + 1) + 1
            end if
        end do
        if (sizes(m + 1) > 0) then
            m = m + 1
        end if
    end do

    ! test without scratch
    ints = ints0
    reals = ints0
    call multiway_merge_float_int(n, reals, ints, sizes(1:m))
    call verify_order(ints, reals)

    ! test with scratch (in the working directory)
    ints = ints0
    reals = ints0
    call multiway_merge_float_int(n, reals, ints, sizes(1:m), './')
    call verify_order(ints, reals)

contains

    subroutine verify_order(ints, reals)

        integer,   intent(in) :: ints(:)
        real(cfp), intent(in) :: reals(:)

        integer :: i

        do i = 1, size(ints)
            if (ints(i) /= i .or. reals(i) /= i) then
                print '(*(a,i0))', 'Sort failed at position ', i, ': ', ints(i)
                stop 1
            end if
        end do

    end subroutine verify_order

end program
