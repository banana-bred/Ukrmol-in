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
module integer_packing

! Module containing data and subprograms for the integer packing / unpacking.
! Joanne Carr June 2010

implicit none
private
save

! Data definitions - local to this module so the default setting is PRIVATE

integer, parameter :: itwo     = 2
integer, parameter :: ithree   = 3
integer, parameter :: ifour    = 4
integer, parameter :: isix     = 6
integer, parameter :: iten     = 10
integer, parameter :: ififteen = 15
integer, parameter :: isixteen = 16
integer, parameter :: i32      = 32

integer, parameter :: two_two     = 2**itwo ! 4
integer, parameter :: two_three   = 2**ithree ! 8
integer, parameter :: two_four    = 2**ifour ! 16
integer, parameter :: two_six     = 2**isix ! 64
integer, parameter :: two_ten     = 2**iten ! 1024
integer, parameter :: two_fifteen = 2**ififteen ! 32768
integer, parameter :: two_sixteen = 2**isixteen ! 65536

! expose the subprogram names
public :: swmol3_pack2
public :: swmol3_pack4
public :: sword_unpack
public :: pack8ints
public :: unpack8ints
public :: ipack9
public :: unpak9
public :: imask
public :: jmask

contains

      integer function swmol3_pack2(i, j, lut)
      ! First packing routine from swmol3 that puts two integers into one.
      ! The first integer is bit-shifted "up" 10 places then the second is added to it.
      implicit none
      
      integer, intent(in) :: i, j ! the two integers to be packed
      integer, intent(in) :: lut  ! the unit number for printing

! Some checks on the input data
      if (i .ge. two_four) then ! this limit arises partly from the second packing routine 
                                ! from swmol3 (see below) and partly from the unpacking as 
                                ! coded in sword originally.
         write(lut,*) ' ERROR in swmol3_pack2: the first integer exceeds the allowed range: '
         write(lut,*) ' i = ',i,' >= ',two_four
         stop
      else if (j .ge. two_ten) then
         write(lut,*) ' ERROR in swmol3_pack2: the second integer exceeds the allowed range: '
         write(lut,*) ' j = ',j,' >= ',two_ten
         stop
      end if

      swmol3_pack2 = i*two_ten + j

      end function swmol3_pack2

      integer function swmol3_pack4(i, j, lut)
      ! Second packing routine from swmol3 that puts two integers into one.
      ! The first integer is bit-shifted "up" 15 places then the second is added to it.
      ! The two integers to be packed are both packed integers already from swmol3_pack2.
      ! The corresponding unpacking is done in sword: see sword_unpack below.
      implicit none

      integer, intent(in) :: i, j ! the two integers to be packed
      integer, intent(in) :: lut  ! the unit number for printing

! No checks on the input data here - to save CPU time we're relying on the checks in swmol3_pack2

      swmol3_pack4 = i*two_fifteen + j

      end function swmol3_pack4

      subroutine sword_unpack(ipacked, i, j, k, l)
      ! Routine that unpacks the integers packed in swmol3_pack4 
      ! We are pulling out bits 25 to 28 into i
      !                    bits 15 to 24 into j
      !                    bits 10 to 13 into k
      !                    bits 0  to 9  into l
      implicit none
  
      integer, intent(in) :: ipacked
      integer, intent(out) :: i, j, k, l

      ! The bit shifting operation ISHFT with a negative 2nd argument returns the 1st argument logically
      ! shifted right, so that the bits we're interested in are in the rightmost (lowest) positions.
      ! Note that this will only work correctly with a little-endian architecture!!!
      ! IAND (logical AND intrinsic function) masks the 1st argument with the 2nd in a bitwise fashion.
      i=IAND(ISHFT(ipacked,-(iten+ififteen)),two_four-1) ! two_four - 1 is 15, i.e. the lowest four bits are 1's; rest 0
      j=IAND(ISHFT(ipacked,-ififteen),two_ten-1) ! two_ten - 1 is 1023 i.e. the lowest 10 bits are 1's; rest 0
      k=IAND(ISHFT(ipacked,-iten),two_four-1)
      l=IAND(ipacked,two_ten-1)

      end subroutine sword_unpack

!*==unpack8ints.spg  processed by SPAG 6.56Rc at 17:23 on  5 Nov 2010
! Author: Dermot Madden
! scatci - unpacking routine in scatci
! Packs eight 16-bit integers from an array of two 64-bit integers
      SUBROUTINE unpack8ints(Z,U)
      USE precisn, ONLY : longint ! for 64-bit integers
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER, DIMENSION(8) :: U
      INTEGER(KIND=longint), DIMENSION(2) :: Z
      INTENT (IN) Z
      INTENT (OUT) U
!
! Local variables
!
      INTEGER(KIND=longint) :: MASK=INT(two_sixteen - 1,KIND=longint) ! 65535, i.e. the number whose lowest 16 bits are all 1; 
      ! the rest are 0.  This is a longint so that the arguments to IAND are of the same type, as required by the F95 standard.
      INTEGER(KIND=longint) :: NX
!
!*** End of declarations rewritten by SPAG
!
      nx=Z(1)
      U(4)=iand(nx,mask)
      nx=ishft(nx,-isixteen)
      U(3)=iand(nx,mask)
      nx=ishft(nx,-isixteen)
      U(2)=iand(nx,mask)
      nx=ishft(nx,-isixteen)
      U(1)=iand(nx,mask)
 
 
      nx=Z(2)
      U(8)=iand(nx,mask)
      nx=ishft(nx,-isixteen)
      U(7)=iand(nx,mask)
      nx=ishft(nx,-isixteen)
      U(6)=iand(nx,mask)
      nx=ishft(nx,-isixteen)
      U(5)=iand(nx,mask)
 
      RETURN
      END SUBROUTINE UNPACK8INTS

!*==pack8ints.spg  processed by SPAG 6.56Rc at 17:23 on  5 Nov 2010
! Author: Dermot Madden
! scatci - packing routine in scatci
! Packs eight 16-bit integers into an array of two 64-bit integers
! Note that the integers passes to the funcation are 32bit, but the
! routine assumes that they have values below 65535 (ie assumed to be
! at most 16bits)
! JMC addendum to the above - care must be taken using the sign 
! bit of integers for packing when default 32-bit integers are converted 
! to 64-bit integers to match kinds.  If the 32-bit integer is negative, 
! then the corresponding 64-bit number will have all the extra bits set to 
! 1, i.e bits 32-63.  This will be fine when ix and jx are IOR'd and put 
! into Z(1) because Z(1) is then left-shifted by 32 bits.  However when 
! kx and lx are IOR'd and the result is IOR'd with the current Z(1) then 
! the 32 left-most bits of Z(1) will all be set to 1 when kx is negative.
! The same is obviously true for Z(2).
! So the conclusion is that I3 and I7 must be less than 2^15, while the 
! other arguments must be less than 2^16.  I've added a test on the integers 
! but it might cause a performance hit...
      SUBROUTINE pack8ints(I1,I2,I3,I4,I5,I6,I7,I8,Z)
      USE precisn, ONLY : longint ! for 64-bit integers
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: I1, I2, I3, I4, I5, I6, I7, I8
      INTEGER(KIND=longint), DIMENSION(2) :: Z
      INTENT (IN) I1, I2, I3, I4, I5, I6, I7, I8
      INTENT (INOUT) Z
!
! Local variables
!
      INTEGER :: IX, JX, KX, LX
!
!*** End of declarations rewritten by SPAG
!
! Some checks on the input data
      if (I1 .ge. two_sixteen .or. I3 .ge. two_fifteen .or. &
          I5 .ge. two_sixteen .or. I7 .ge. two_fifteen .or. &
          I2 .ge. two_sixteen .or. I4 .ge. two_sixteen .or. &
          I6 .ge. two_sixteen .or. I8 .ge. two_sixteen) then
         write(*,*) ' ERROR in pack8ints: an integer exceeds the allowed range'
         write(*,*) ' I3,I7 = ',I3,I7,' but the limit is ',two_fifteen
         write(*,*) ' I1,I2,I4,I5,I6,I8 = ',I1,I2,I4,I5,I6,I8,' but the limit is ',two_sixteen
         stop
      end if

      ix=ishft(I1,isixteen)
      jx=I2
      kx=ishft(I3,isixteen)
      lx=I4
      Z(1)=IOR(ix,jx)
      Z(1)=ishft(Z(1),i32)
      Z(1)=IOR(Z(1),INT(IOR(kx,lx),KIND=longint))
 
      ix=ishft(I5,isixteen)
      jx=I6
      kx=ishft(I7,isixteen)
      lx=I8
      Z(2)=IOR(ix,jx)
      Z(2)=ishft(Z(2),i32)
      Z(2)=IOR(Z(2),INT(IOR(kx,lx),KIND=longint))
 
      RETURN
      END SUBROUTINE PACK8INTS

!*==ipack9.spg  processed by SPAG 6.56Rc at 14:09 on 18 Nov 2010
! Author: Dermot Madden
! denprop - packing routine in denprop
! packing is specific to the data being packed in this routine
! so was not edited
      FUNCTION IPACK9(IOPCDE,I,MGVNI,J,MGVNJ,NWORD1)
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: I, J, MGVNI, MGVNJ, NWORD1
      INTEGER, DIMENSION(7) :: IOPCDE
      INTEGER :: IPACK9
      INTENT (IN) I, IOPCDE, J, MGVNI, MGVNJ
      INTENT (INOUT) NWORD1
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     IPACK9 - Integer PACKing of 9 variables
!
!     Stores into the integer word, IPACK9, the operator and transition
!     information; state id information is packed into NWORD1 too.
!
!     INPUT DATA:
!         IOPCDE THE FIRST SEVEN NUMERALS DEFINING THE PROPERTY OPERATOR
!              I THE NUMBER OF THE VECTOR WITHIN MANIFOLD I
!          MGVNI LAMBDA QUANTUM NUMBER DEFINING MANIFOLD I
!              J THE NUMBER OF THE VECTOR WITHIN MANIFOLD J
!          MGVNJ LAMBDA QUANTUM NUMBER DEFINING MANIFOLD J
!
!     OUTPUT DATA:
!          IPACK9 COMPACTED FORM OF ALL OF THE ABOVE DATA, EXCLUDING STA
!           NWORD COMPACTED FORM OF STATE IDS.
!
!     Notes:
!       It is assumed that the machine word length is at least 32 bits
!     for an integer. This means that this routine will work on IBM,
!     Vax and Cray computers.
!
!     This packing is used by TMT to store moments and is thus required
!     by codes reading the results.
!
!     Author: Charles J Gillan, February 1989 and March 1994
!
!     Copyright 1994 (c)  Charles J Gillan
!     All rights reserved
!
!***********************************************************************
!
!---- Initialize the two output words now
!
      IPACK9=0
      NWORD1=0
!
!=======================================================================
!
!     IBM System 370 and AIX/6000, Vax/VMS, Microsoft Fortran Version
!
!=======================================================================
!
!---- STORE THE OPERATOR CODE AND THE LAMBDA VALUES INTO IPACK9
!
      NWORD1=ISHFT(IOPCDE(1),29)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(IOPCDE(2),27)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(IOPCDE(3),25)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(IOPCDE(4),23)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(IOPCDE(5),17)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(IOPCDE(6),15)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(IOPCDE(7),13)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(MGVNI,10)
      IPACK9=IOR(IPACK9,NWORD1)
!
      NWORD1=ISHFT(MGVNJ,7)
      IPACK9=IOR(IPACK9,NWORD1)
!
!---- STORE THE STATE VECTOR IDS INTO NWORD1
!
      NWORD1=ISHFT(I,16)
      NWORD1=IOR(NWORD1,J)
!
      RETURN
!
!=======================================================================
!
!     Cray COS and UNICOS version
!
!=======================================================================
!
!---- STORE THE OPERATOR CODE AND THE LAMBDA VALUES INTO IPACK9
!
!      NWORD1=SHIFTL(IOPCDE(1),29)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(IOPCDE(2),27)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(IOPCDE(3),25)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(IOPCDE(4),23)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(IOPCDE(5),17)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(IOPCDE(6),15)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(IOPCDE(7),13)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(MGVNI,10)
!      IPACK9=OR(IPACK9,NWORD1)
!
!      NWORD1=SHIFTL(MGVNJ,7)
!      IPACK9=OR(IPACK9,NWORD1)
!
!---- STORE THE STATE VECTOR IDS INTO NWORD1
!
!      NWORD1=SHIFTL(I,16)
!      NWORD1=OR(NWORD1,J)
!
!      RETURN
!
      END FUNCTION IPACK9

!*==unpak9.spg  processed by SPAG 6.56Rc at 14:09 on 18 Nov 2010
! denprop - unpacking routine in denprop
! unpacking is specific to the data being packed in this routine
! so was not edited
      SUBROUTINE UNPAK9(IOPCDE,MGVNI,MGVNJ,I,J,NWORD,NSTATE)
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: I, J, MGVNI, MGVNJ, NSTATE, NWORD
      INTEGER, DIMENSION(7) :: IOPCDE
      INTENT (IN) NSTATE, NWORD
      INTENT (OUT) I, IOPCDE, J, MGVNI, MGVNJ
!
! Local variables
!
      INTEGER :: MASK=two_two-1, MASK1=two_three-1, &
                 MASK2=two_sixteen-1, MASK3=two_six-1, NX
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************
!
!     UNPACKS IPACK9 THE OPERATOR AND TRANSITION INFORMATION. IT IS
!     ASSUMED THAT THE MACHINE WORD LENGTH IS AT LEAST 32 BITS FOR AN
!     INTEGER. THE ROUTINE WILL THEREFORE WORK ON IBM, CRAY AND VAX COMP
!
!     INPUT DATA:
!         NSTATE PACKED STATE VECTOR IDS
!          NWORD PACKED OPERATOR CODE AND STATE LAMBDA VALUES
!
!     OUTPUT DATA:
!         IOPCDE THE FIRST SEVEN NUMERALS DEFINING THE PROPERTY OPERATOR
!              I THE NUMBER OF THE VECTOR WITHIN MANIFOLD I
!          MGVNI LAMBDA QUANTUM NUMBER DEFINING MANIFOLD I
 
!              J THE NUMBER OF THE VECTOR WITHIN MANIFOLD J
!          MGVNJ LAMBDA QUANTUM NUMBER DEFINING MANIFOLD J
!
!     This packing is used by TMT to store moments and is thus required
!     by codes reading the results.
!
!***********************************************************************
!
!---- MASKS
!
!      DATA MASK/3/MASK1/7/MASK2/65535/, mask3/63/
!
!**** IBM AND VAX VERSION *************
!
!     UNPACK THE OPERATOR CODE AND LAMBDA VALUES FROM NWORD.
!     NOTE THAT IOPCDE(1) IS SPECIAL GIVEN THAT IT CAN BE 4 !
!     and THAT IOPCDE(5) (the label of the scattering centre) can be eve
!     larger, so we allocate 6 bits, allowing up to 63 centres.
!
      NX=ISHFT(NWORD,-29)
      IOPCDE(1)=IAND(NX,MASK1)
!
      NX=ISHFT(NWORD,-27)
      IOPCDE(2)=IAND(NX,MASK)
!
      NX=ISHFT(NWORD,-25)
      IOPCDE(3)=IAND(NX,MASK)
!
      NX=ISHFT(NWORD,-23)
      IOPCDE(4)=IAND(NX,MASK)
!
      NX=ISHFT(NWORD,-17)
      IOPCDE(5)=IAND(NX,MASK3)
!
      NX=ISHFT(NWORD,-15)
      IOPCDE(6)=IAND(NX,MASK)
!
      NX=ISHFT(NWORD,-13)
      IOPCDE(7)=IAND(NX,MASK)
!
! Michal Bug fix- incorrect unpacking of mgvni and mgvnj (unpacked 2-bits in
! 3-bits packed variable)
      NX=ISHFT(NWORD,-10)
      MGVNI=IAND(NX,MASK1)
!
      NX=ISHFT(NWORD,-7)
      MGVNJ=IAND(NX,MASK1)
!
!---- UNPACK THE STATE VECTOR IDS FROM NWORD1
!
      I=ISHFT(NSTATE,-16)
      J=IAND(NSTATE,MASK2)
!
!**** CRAY VERSION *************
!
!     UNPACK THE OPERATOR CODE AND LAMBDA VALUES FROM NWORD.
!     NOTE THAT IOPCDE(1) IS SPECIAL GIVEN THAT IT CAN BE 4 !
!
!      NX=SHIFTR(NWORD,29)
!      IOPCDE(1)=AND(NX,MASK1)
!
!      NX=SHIFTR(NWORD,27)
!      IOPCDE(2)=AND(NX,MASK)
!
!      NX=SHIFTR(NWORD,25)
!      IOPCDE(3)=AND(NX,MASK)
!
!      NX=SHIFTR(NWORD,23)
!      IOPCDE(4)=AND(NX,MASK)
!
!      NX=SHIFTR(NWORD,17)
!      IOPCDE(5)=AND(NX,MASK3)
!
!      NX=SHIFTR(NWORD,15)
!      IOPCDE(6)=AND(NX,MASK)
!
!      NX=SHIFTR(NWORD,13)
!      IOPCDE(7)=AND(NX,MASK)
!
!      NX=SHIFTR(NWORD,10)
!      MGVNI=AND(NX,MASK)
!
!      NX=SHIFTR(NWORD,7)
!      MGVNJ=AND(NX,MASK)
!
!---- UNPACK THE STATE VECTOR IDS FROM NWORD1
!
!      I=SHIFTR(NSTATE,16)
!      J=AND(NSTATE,MASK2)
!
      RETURN
!
      END SUBROUTINE UNPAK9

!*==imask.spg  processed by SPAG 6.56Rc at 10:27 on 24 Nov 2010
! Author: Dermot Madden
!     pword can be considered as 4 blocks of 16bits. This function
!     isolates the 1st block (the lowest order 16 bit) converts
!     those 16bits to an integer (i.e. same bit pattern, different
!     type of number) and returns that value
! JMC Since this routine is hardwired to select only the lowest 16
! JMC bits, it will be safe as-is for both 32-bit and 64-bit integers (etc).
! JMC There was an identical routine to this one in swtrmo.f
 
      FUNCTION IMASK(pword)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      REAL(KIND=wp) :: PWORD
      INTEGER :: IMASK
      INTENT (IN) PWORD
!
! Local variables
!
      INTEGER :: IWORD, MASK=two_sixteen-1 ! 65535, i.e. the number whose lowest 16 bits are all 1; the rest are 0
      REAL(KIND=wp) :: TWORD
!
!*** End of declarations rewritten by SPAG
!
      EQUIVALENCE(tword,iword)
!
!---- Set the value of MASK using Mil-Std-53 function
!
!      MASK=IZERO
!
!      DO I=15,0,-1
!         MASK=IBSET(MASK,I)
!      END DO
!
      tword=pword
      IMASK=iand(iword,MASK)
 
      RETURN
      END FUNCTION IMASK

!*==jmask.spg  processed by SPAG 6.56Rc at 10:27 on 24 Nov 2010
! Author: Dermot Madden
!      This function appears to be broken.
!      It returns an integer, whose first 16bits are the
!      first 16bits of kl, and whose second 16bits are an
!      Inclusive OR operation on the second 16bits of pword
!      and the second 16bits of kl. This sequence of 32bits is
!      then converted to an integer and returned.
!      The function is called in an IF statement, so it is possible
!      that the condition is never met.
! JMC There was an identical routine to this one in swtrmo.f
! JMC depending on how the unpacking is done, this may not work correctly 
! for 64-bit integers.
 
      FUNCTION JMASK(pword,kl)
      USE precisn, ONLY : wp ! for specifying the kind of reals                        
      IMPLICIT NONE
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
      INTEGER :: KL
      REAL(KIND=wp) :: PWORD
      INTEGER :: JMASK
      INTENT (IN) KL, PWORD
!
! Local variables
!
      INTEGER :: IWORD, LL, MASK=two_sixteen-1 ! 65535, i.e. the number whose lowest 16 bits are all 1; the rest are 0
      REAL(KIND=wp) :: TWORD
!
!*** End of declarations rewritten by SPAG
!
      EQUIVALENCE(tword,iword)
!
      tword=pword
      ll=kl
      jmask=ior(iand(iword,not(MASK)),ll)
 
      RETURN
      END FUNCTION JMASK

end module integer_packing
