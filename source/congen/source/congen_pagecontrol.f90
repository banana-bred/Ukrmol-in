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
!> \brief Page control.
!>
!> Set of routines for pagination of the text output from the program.
!> Every couple of lines a page header is produced.
!>
module congen_pagecontrol

    implicit none

    ! subroutines used from other modules
    public addl
    public ctlpg1
    public newpg
    public space
    public taddl
    public taddl1

contains

    !> \brief Add blank lines to output
    !>
    !> Increments the output line counter by the given number of lines. If the result exceeds the built-in
    !> limit, a new page header is printed. Note that the lines are virtual, none is actually printed if 
    !> the pagination threshold is not reached.
    !>
    !> This function is mostly used before writing a known number of lines, so that all lines to be written
    !> appear on the same page.
    !>
    subroutine addl (lines)

        use global_utils, only : getin
        use congen_data,  only : head, lppr, lppmax, npage, nftw

        integer, intent(in) :: lines

        lppr = lppr + lines

        if (lppr > lppmax) then
            lppr = lines
            npage = mod(npage + 1, 1000)
            
            ! JMC this is putting the max. 3-digit number npage into the last 3 elements of HEAD
            call getin (npage, 3, head(size(head) - 2), 1)
            write(nftw, '("1",132A1,/)') head
        end if

    end subroutine addl


    !> \brief Controls page layout and counting
    !>
    subroutine ctlpg1 (lpp, h1, nh1, h2, nh2)

        use congen_data, only : head, lppr, lppmax, npage

        integer :: lpp, nh1, nh2
        character(len=1), dimension(nh1) :: h1
        character(len=1), dimension(nh2) :: h2
        intent (in) h1, h2, lpp, nh1, nh2

        character(len=1) :: blank = ' '
        integer :: i, il, im, ip
        character(len=1), dimension(4) :: page = (/ 'P', 'A', 'G', 'E' /)

        lppmax = lpp
        if (lpp <= 0 .or. lpp > 60) lppmax = 60 ! JMC this line sets the default value for LPPMAX (see item LPP in the documentation)
        lppr = lppmax
        npage = 0
        do i = 1, size(head)
            head(i) = blank
        end do

        im = min(nh1, 120)
        head(1:im) = h1(1:im)

        if (im /= 120) then
            il = im + 1
            im = min(im + nh2, 120)
            ip = 0
            do i = il, im
                ip = ip + 1
                head(i) = h2(ip)
            end do
        end if

        im = size(head) - 8 ! JMC changing from 124, as the last 8 elements of HEAD will contain 'PAGE' then ' ' then a 3-digit number 

        do i = 1, 4
            im = im + 1
            head(im) = page(i)
        end do

    end subroutine ctlpg1


    !> \brief Start new page of output
    !>
    !> Resets the line counter, increments the page counter and prints the page header on standard output.
    !>
    subroutine newpg

        use global_utils, only : getin
        use congen_data,  only : head, lppr, npage, nftw

        lppr  = 0
        npage = mod(npage + 1, 1000)

        ! JMC this is putting the max. 3-digit number npage into the last 3 elements of HEAD
        call getin (npage, 3, head(size(head) - 2), 1)

        write(nftw, '("1",132A1,/)') head

    end subroutine newpg


    !> \brief Add blank lines to output
    !>
    !> Increments the global line counter by the given number of lines. If the pagination limit is not reached,
    !> it will actually write those empty lines (containing one space each) to output.
    !>
    !> \param lines Number of blank lines
    !>
    subroutine space (lines)

        use congen_data, only : lppr, lppmax, nftw

        integer, intent(in) :: lines

        integer :: j

        lppr = lppr + lines

        if (lppr < lppmax) then
            do j = 1, lines
                write(nftw, '(" ")')
            end do
        end if

    end subroutine space


    !> \brief Add blank lines to output
    !>
    !> Advance the internal output line counter and check if pagination limit has been hit.
    !>
    !> \param lines Number of lines to add
    !> \param lt    Output status: 0 for no lines added (would overflow page), 1 for all lines added.
    !>
    subroutine taddl (lines, lt)

        use congen_data, only : lppr, lppmax

        integer, intent(in)  :: lines
        integer, intent(out) :: lt

        lt = 0
        if (lppr + lines <= lppmax) then
            lppr = lppr + lines
            lt = 1
        end if

    end subroutine taddl


    !> \brief Check available lines.
    !>
    !> Checks if the given number of lines fits to page.
    !>
    !> \param lines  Number of lines to test.
    !> \param lt     On return 0 if lines do not fit on page, 1 if they do.
    !>
    subroutine taddl1 (lines, lt)

        use congen_data, only : lppr, lppmax

        integer :: lines, lt
        intent (in) lines
        intent (out) lt

        lt = 0
        if (lppr + lines > lppmax) return
        lt = 1

    end subroutine taddl1


end module
