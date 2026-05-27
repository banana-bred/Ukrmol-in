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

!> \brief   Utility module
!> \authors A Al-Refaie
!> \date    2017
!>
!> \note 16/01/2019 - Jakub Benda: Unifom coding style and expanded documentation.
!>
module Utility_module

    use consts_mpi_ci, only: NIDX
    use precisn,       only: longint, wp
    use mpi_gbl,       only: mpi_mod_wtime, mpi_xermsg

    implicit none

    public string_hash, get_real_time, get_cpu_time, compute_total_triangular, triangular_index_to_ij
    public compute_total_box, box_index_to_ij

contains

    !> \brief   Calculate triangular area
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Calculate n*(n+1)/2.
    !>
    integer(longint) function compute_total_triangular (n)
        integer, intent(in) :: n

        compute_total_triangular = n * (n + 1_longint) / 2

    end function compute_total_triangular


    !> \brief   Calculate rectangular product
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Calculate width*height.
    !>
    integer function compute_total_box (width, height)
        integer, intent(in) :: width, height

        compute_total_box = width * height

    end function compute_total_box


    !> \brief   Extract indices from rectangular multi-index
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine box_index_to_ij (idx, height, i, j)
        integer, intent(in)  :: idx, height
        integer, intent(out) :: i, j

        i = mod(idx - 1, height) + 1
        j = (idx - 1) / height + 1

    end subroutine box_index_to_ij


    !> \brief   Extract indices from triangular multi-index
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    subroutine triangular_index_to_ij (idx_f, N, row, column)
        integer, intent(in)  :: idx_f, N
        integer, intent(out) :: row, column
        integer              :: idx, ii, K, jj

        idx = idx_f - 1
        ii = N * (N + 1) / 2 - 1 - idx
        K = int((sqrt(8.0_wp*real(ii,wp) + 1.0_wp) - 1.0_wp) / 2.0_wp, longint)
        jj = ii - K * (K + 1) / 2

        row = N - 1 - K
        column = N - 1 - jj

    end subroutine triangular_index_to_ij


    !> \brief   Calculate a string hash
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> This hash assumes at least 29-bit integers. It is supposedly
    !> documented in Aho, Sethi, and Ullman, pp. 434-438
    !>
    integer function string_hash (str, table_size) result(h)
        character(len=*), intent(in) :: str
        integer,          intent(in) :: table_size
        integer :: i, chr, g, mask = int(Z"1FFFFFF")

        h = 0
        do i = 1, len_trim(str)
            chr = ichar(str(i:i))
            h   = ishft(h,  4) + chr
            g   = ishft(h,-24)
            h   = iand(ieor(h,g), mask)
        end do
        h = 1 + modulo(h, table_size)

    end function string_hash


    !> \brief   Get current (real) time
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    !> Uses a function from GBTOlib.
    !>
    function get_real_time () result(t)
        real(wp) :: t

        t = mpi_mod_wtime()

    end function get_real_time


    !> \brief   Get current (CPU) time
    !> \authors A Al-Refaie
    !> \date    2017
    !>
    function get_cpu_time () result(t)
        real(wp) :: t

        call cpu_time(t)

    end function get_cpu_time


    !> \brief   Store 8 integers in given number of integers
    !> \authors J Benda
    !> \date    2023
    !>
    !> Copy or pack 8 integers to the given packed storage.
    !>
    subroutine pack_ints (i1, i2, i3, i4, i5, i6, i7, i8, z)

        integer,          intent(in)  :: i1, i2, i3, i4, i5, i6, i7, i8
        integer(longint), intent(out) :: z(:)

        integer(longint), parameter :: max16 = 2_longint**16 - 1
        integer(longint), parameter :: max32 = 2_longint**32 - 1

        select case (NIDX)

            case (2)

                if (i1 > max16 .or. i2 > max16 .or. i3 > max16 .or. i4 > max16 .or. &
                    i5 > max16 .or. i6 > max16 .or. i7 > max16 .or. i8 > max16) then
                    call mpi_xermsg('Utility_module', 'pack_ints', &
                                    'Cannot pack integers, values too large. Try increasing NIDX.', 2, 1)
                end if

                z(1) = i1
                z(1) = ior(ishft(z(1), 16_longint), int(i2, longint))
                z(1) = ior(ishft(z(1), 16_longint), int(i3, longint))
                z(1) = ior(ishft(z(1), 16_longint), int(i4, longint))
                z(2) = i5
                z(2) = ior(ishft(z(2), 16_longint), int(i6, longint))
                z(2) = ior(ishft(z(2), 16_longint), int(i7, longint))
                z(2) = ior(ishft(z(2), 16_longint), int(i8, longint))

            case (4)

                if (i1 > max32 .or. i2 > max32 .or. i3 > max32 .or. i4 > max32 .or. &
                    i5 > max32 .or. i6 > max32 .or. i7 > max32 .or. i8 > max32) then
                    call mpi_xermsg('Utility_module', 'pack_ints', &
                                    'Cannot pack integers, values too large. Try increasing NIDX.', 4, 1)
                end if

                z(1) = i1
                z(1) = ior(ishft(z(1), 32_longint), int(i2, longint))
                z(2) = i3
                z(2) = ior(ishft(z(2), 32_longint), int(i4, longint))
                z(3) = i5
                z(3) = ior(ishft(z(3), 32_longint), int(i6, longint))
                z(4) = i7
                z(4) = ior(ishft(z(4), 32_longint), int(i8, longint))

            case (8)

                z(1) = i1
                z(2) = i2
                z(3) = i3
                z(4) = i4
                z(5) = i5
                z(6) = i6
                z(7) = i7
                z(8) = i8

            case default

                call mpi_xermsg('Utility_module',  'pack_ints', 'Packing implemented only for NIDX = 2, 4 or 8.', 1, 1)

        end select

    end subroutine pack_ints


    !> \brief   Retrieve 8 integers from given number of integers
    !> \authors J Benda
    !> \date    2023
    !>
    !> Copy or unpack 8 default integers from the provided packed storage.
    !>
    subroutine unpack_ints (z, u)

        integer(longint), intent(in)  :: z(:)
        integer,          intent(out) :: u(8)

        integer(longint), parameter :: mask16 = 2_longint**16 - 1
        integer(longint), parameter :: mask32 = 2_longint**32 - 1

        integer(longint) :: n

        select case (NIDX)

            case (2)

                n = z(1);                   u(4) = int(iand(n, mask16))
                n = ishft(n, -16_longint);  u(3) = int(iand(n, mask16))
                n = ishft(n, -16_longint);  u(2) = int(iand(n, mask16))
                n = ishft(n, -16_longint);  u(1) = int(iand(n, mask16))
                n = z(2);                   u(8) = int(iand(n, mask16))
                n = ishft(n, -16_longint);  u(7) = int(iand(n, mask16))
                n = ishft(n, -16_longint);  u(6) = int(iand(n, mask16))
                n = ishft(n, -16_longint);  u(5) = int(iand(n, mask16))

            case (4)

                n = z(1);                   u(2) = int(iand(n, mask32))
                n = ishft(n, -32_longint);  u(1) = int(iand(n, mask32))
                n = z(2);                   u(4) = int(iand(n, mask32))
                n = ishft(n, -32_longint);  u(3) = int(iand(n, mask32))
                n = z(3);                   u(6) = int(iand(n, mask32))
                n = ishft(n, -32_longint);  u(5) = int(iand(n, mask32))
                n = z(4);                   u(8) = int(iand(n, mask32))
                n = ishft(n, -32_longint);  u(7) = int(iand(n, mask32))

            case (8)

                u(1) = int(z(1))
                u(2) = int(z(2))
                u(3) = int(z(3))
                u(4) = int(z(4))
                u(5) = int(z(5))
                u(6) = int(z(6))
                u(7) = int(z(7))
                u(8) = int(z(8))

            case default

                call mpi_xermsg('Utility_module',  'unpack_ints', 'Unpacking implemented only for NIDX = 2, 4 or 8.', 1, 1)

        end select

    end subroutine unpack_ints


    !> \brief   Compare two arrays
    !> \authors J Benda
    !> \date    2023
    !>
    !> Lexicographically compare two integer arrays of equal length NIDX. Return 0 if the arrays are indentical, -1 if the first
    !> unequal element is smaller in the first array that in the second array, and +1 otherwise.
    !>
    integer function lexicographical_compare (a, b) result (verdict)

        integer(longint), intent(in) :: a(NIDX), b(NIDX)

        integer :: i

        verdict = 0

        do i = 1, NIDX
            if (a(i) < b(i)) then
                verdict = -1
                exit
            end if
            if (a(i) > b(i)) then
                verdict = +1
                exit
            end if
        end do

    end function lexicographical_compare

end module Utility_module
