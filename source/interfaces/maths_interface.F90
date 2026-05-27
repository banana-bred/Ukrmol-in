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

!> \brief   Interfacing module for library mathematical routines
!> \author  Joanne Carr, J Benda
!> \date    2010 - 2019
!>
!> Some subroutines can optionally interface blas95/lapack95 routines from Intel MKL (when compiled with -DBLA_F95).
!>
module maths

    use precisn, only: wp ! for specifying the kind of reals

    implicit none

    private

    public :: maths_erfc ! calculates the complement of the error function erf(x)
    public :: diag       ! matrix diagonalizer for psn
    public :: maths_tcoef ! simple matrix transposition
    public :: maths_cveca_plus_vecb ! constant*veca + vecb
    public :: maths_mxm ! matrix multplication in swtrmo; should be merged with maths_matrix_multiply at some point??

    public :: maths_matrix_multiply_blas95
    public :: maths_dmatrix_multiply_blas95
    public :: maths_zmatrix_multiply_blas95

    interface  maths_matrix_multiply_blas95
        module procedure  maths_dmatrix_multiply_blas95
        module procedure  maths_zmatrix_multiply_blas95
    end interface

contains

    !> Calculates the complement of the error function erf(x)
    !> Note that ifail is always set to 0 in s15erfc, i.e.
    !> no failure reported on exit
    real(wp) function maths_erfc (x) result (y)

        use global_utils, only: s15erfc ! the location of the actual function

        real(kind=wp), intent(in) :: x
        integer :: ifail

        y = s15erfc(x, ifail)

    end function maths_erfc


    !> DIAG DIAGONALIZES THE DENSITY MATRIX BLOCKS. THIS USED A NAG
    !> LIBRARY ROUTINE F02FAF TO PERFORM THE DIAGONALIZATION.
    !> This routine has been replaced by the lapack double
    !> precision routine dsyev (nv 16/5/02)
    !>
    !> INPUT DATA :
    !>          H  LOWER TRIANGLE OF THE DENSITY MATRIX (stored as a
    !>             linear array with MK*(MK+1)/2 elements)
    !>         MK  SIZE OF THE DENSITY MATRIX IS MK BY MK
    !>      IPRNT  SWITCH FOR DEBUG PRINTING
    !>
    !> OUTPUT DATA :
    !>         EIG  EIGEN-ENERGIES OF THE DENSITY MATRIX
    !>         VEC  EIGEN-VECTORS OF THE DENSITY MATRIX
    !>
    !> TO SAVE SPACE AND TO COMPLY WITH THE REQUIREMENTS OF THE NAG
    !> ROUTINE, THE LOWER TRIANGLE OF THE DENSITY MATRIX IS COPIED
    !> INTO THE ARRAY VEC BEFORE THE CALL. IT IS OVERWRITTEN ON OUTPUT
    !> BY THE EIGEN-VECTORS.
    SUBROUTINE DIAG (MK, H, EIG, VEC, IPRNT)

        ! Dummy arguments
        INTEGER :: IPRNT, MK
        REAL(KIND=wp), DIMENSION(MK) :: EIG
        REAL(KIND=wp), DIMENSION(*) :: H
        REAL(KIND=wp), DIMENSION(MK,MK) :: VEC
        INTENT (IN) H, IPRNT
        INTENT (INOUT) VEC

        ! Local variables
        INTEGER :: I, IFAIL, J, K, LX, NDIM
        REAL(KIND=wp), DIMENSION(3*mk) :: X

 9000   FORMAT(/,10X,' DEBUG INFORMATION IN SUBROUTINE DIAG',//)
 9010   FORMAT(/,10X,' HALF TRIANGLE DIMENSION = ',I5,/)
 9020   FORMAT(/,10X,' H-MATRIX ',/)
 9030   FORMAT(2X,6F12.8)
 9040   FORMAT(/,10X,'VEC MATRIX ',/)
 9050   FORMAT(/,10X,'EIGENVALUES ',/)
 9060   FORMAT((/6(I3,F10.5)))
 9070   FORMAT(/,10X,'EIGENVECTORS ',/)
 9080   FORMAT(2X,I5,6F12.5,/,(7X,6F12.5))
 9999   FORMAT(/,10X,'DIAG : LIBRARY ROUTINE RETURNS IFAIL = ',I2,/)
 9888   FORMAT(/,10X,'ERROR IN DIAG : INPUT MATRIX DIMENSION MK = ',I5,/)

        IFAIL = 0

        ! COPY LOWER TRIANGLE INTO LOWER HALF OF VEC.
        K=0
        DO I=1, MK
            DO J=1, I
                K=K+1
                VEC(I,J)=H(K)
                VEC(J,I)=H(K)
            END DO
        END DO


        ! DEBUG CODE. CHECKS IF EVERYTHING IS STORED OK.
        IF (IPRNT > 1) THEN
            NDIM=(MK*(MK+1))/2
            WRITE(*,9000)
            WRITE(*,9010)NDIM
            WRITE(*,9020)
            WRITE(*,9030)(H(K),K=1,NDIM)
            WRITE(*,9040)
            DO I=1, MK
                WRITE(*,9030)(VEC(J,I),J=1,MK)
            END DO
        END IF

        lx=3*mk

        IF (mk > 0) THEN
            CALL dsyev('V','L',mk,vec,mk,eig,x,lx,ifail)
            IF (IFAIL /= 0) THEN
                WRITE(*,9999)IFAIL
                STOP
            END IF
        ELSE
            WRITE(*,9888)MK
            STOP
        END IF

        ! DEBUG CODE PRINTS OUT EIGEN VALUES AND VECTORS.
        IF (IPRNT > 1) THEN
            WRITE(*,9050)
            WRITE(*,9060)(I,EIG(I),I=1,MK)
            WRITE(*,9070)
            DO I=1, MK
                WRITE(*,9080)I, (VEC(J,I),J=1,MK)
            END DO
        END IF

        RETURN

    END SUBROUTINE DIAG


    ! (m by n) matrix C = (m by k) matrix A * (k by n) matrix B
    subroutine maths_dmatrix_multiply_blas95 (a, b, c, transpose_A, transpose_B, coeff_A, coeff_C)
        use blas_lapack_gbl, only: blasint
#ifdef BLA_F95
        use blas95, only: gemm
#endif
        integer(blasint) :: m, n, k, OP_a_rows, OP_a_cols, OP_b_rows, OP_b_cols, OP_c_rows, OP_c_cols
        integer(blasint) :: lda, ldb, ldc ! leading dimensions of A, B and C, respectively
                                          ! (should be >= m, k, m otherwise it's an error).
                                          ! In dgemm, the A array is declared with dimensions
                                          ! (lda,*) and equivalently for B and C.
        real(wp), intent(in)    :: a(:,:)
        real(wp), intent(in)    :: b(:,:)
        real(wp), intent(inout) :: c(:,:)
        real(wp), intent(in), optional :: coeff_A, coeff_C
        logical,  intent(in), optional :: transpose_A, transpose_B

        real(wp)         :: coeff_A_local, coeff_C_local
        character(len=1) :: transpose_A_local, transpose_B_local

        ! Default values
        transpose_A_local = 'N'
        transpose_B_local = 'N'
        coeff_A_local = 1.0_wp
        coeff_C_local = 0.0_wp

        ! Adjust for optional parameters
        if (present(transpose_A)) then
            if (transpose_A) then
                transpose_A_local = 'T'
            end if
        end if
        if (present(transpose_B)) then
            if (transpose_B) then
                transpose_B_local = 'T'
            end if
        end if
        if (present(coeff_A)) coeff_A_local = coeff_A
        if (present(coeff_C)) coeff_C_local = coeff_C

#ifdef BLA_F95
        call gemm(a, b, c, transpose_A_local, transpose_B_local, coeff_A_local, coeff_C_local)
#else
        ! Calculate shapes
        OP_a_rows = merge(size(a,1), size(a,2), transpose_A_local == 'N')
        OP_a_cols = merge(size(a,2), size(a,1), transpose_A_local == 'N')
        OP_b_rows = merge(size(b,1), size(b,2), transpose_B_local == 'N')
        OP_b_cols = merge(size(b,2), size(b,1), transpose_B_local == 'N')
        OP_c_rows = size(c,1)
        OP_c_cols = size(c,2)

        m = min(OP_a_rows,OP_c_rows)
        k = min(OP_a_cols,OP_b_rows)
        n = min(OP_b_cols,OP_c_cols)

        lda = size(a,1)
        ldb = size(b,1)
        ldc = size(c,1)

        ! call the Level 3 blas routine, assuming that we don't need to
        ! transpose either input matrix (hence the two 'n's), and that we
        ! simply want to calculate c = a*b  i.e. c = 1.0_wp*a*b + 0.0_wp*c
        !
        ! dgemm will also do some error checking.  For example, it will stop
        ! if any of m,n,k are < 0

       call dgemm(transpose_A_local, transpose_B_local, m, n, k, coeff_A_local, a, lda, b, ldb, coeff_C_local, c, ldc)
#endif
    end subroutine maths_dmatrix_multiply_blas95


    ! (m by n) matrix C = (m by k) matrix A * (k by n) matrix B
    subroutine maths_zmatrix_multiply_blas95 (a, b, c, transpose_A, transpose_B, coeff_A, coeff_C)
        use blas_lapack_gbl, only: blasint
#ifdef BLA_F95
        use blas95, only: gemm
#endif
        integer(blasint) :: m, n, k, OP_a_rows, OP_a_cols, OP_b_rows, OP_b_cols, OP_c_rows, OP_c_cols
        integer(blasint) :: lda, ldb, ldc ! leading dimensions of A, B and C, respectively
                                          ! (should be >= m, k, m otherwise it's an error).
                                          ! In dgemm, the A array is declared with dimensions
                                          ! (lda,*) and equivalently for B and C.
        complex(wp), intent(in)    :: a(:,:)
        complex(wp), intent(in)    :: b(:,:)
        complex(wp), intent(inout) :: c(:,:)
        complex(wp), intent(in), optional :: coeff_A, coeff_C
        logical,     intent(in), optional :: transpose_A, transpose_B

        complex(wp)      :: coeff_A_local, coeff_C_local
        character(len=1) :: transpose_A_local, transpose_B_local

        ! Default values
        transpose_A_local = 'N'
        transpose_B_local = 'N'
        coeff_A_local = 1.0_wp
        coeff_C_local = 0.0_wp

        ! Adjust for optional parameters
        if (present(transpose_A)) then
            if (transpose_A) then
                transpose_A_local = 'T'
            end if
        end if
        if (present(transpose_B)) then
            if (transpose_B) then
                transpose_B_local = 'T'
            end if
        end if
        if (present(coeff_A)) coeff_A_local = coeff_A
        if (present(coeff_C)) coeff_C_local = coeff_C

#ifdef BLA_F95
        call gemm(a, b, c, transpose_A_local, transpose_B_local, coeff_A_local, coeff_C_local)
#else
        ! Calculate shapes
        OP_a_rows = merge(size(a,1), size(a,2), transpose_A_local == 'N')
        OP_a_cols = merge(size(a,2), size(a,1), transpose_A_local == 'N')
        OP_b_rows = merge(size(b,1), size(b,2), transpose_B_local == 'N')
        OP_b_cols = merge(size(b,2), size(b,1), transpose_B_local == 'N')
        OP_c_rows = size(c,1)
        OP_c_cols = size(c,2)

        m = min(OP_a_rows,OP_c_rows)
        k = min(OP_a_cols,OP_b_rows)
        n = min(OP_b_cols,OP_c_cols)

        lda = size(a,1)
        ldb = size(b,1)
        ldc = size(c,1)

        ! call the Level 3 blas routine, assuming that we don't need to
        ! transpose either input matrix (hence the two 'n's), and that we
        ! simply want to calculate c = a*b  i.e. c = 1.0_wp*a*b + 0.0_wp*c
        !
        ! dgemm will also do some error checking.  For example, it will stop
        ! if any of m,n,k are < 0

        call zgemm(transpose_A_local, transpose_B_local, m, n, k, coeff_A_local, a, lda, b, ldb, coeff_C_local, c, ldc)
#endif
    end subroutine maths_zmatrix_multiply_blas95


    ! Utility to transpose a matrix. A is transposed and written to B leaving A unchanged.
    subroutine maths_tcoef(a,b,nao,nmo)

        ! dummy arguments
        integer :: nao, nmo
        real(kind=wp), dimension(nao,nmo) :: a
        real(kind=wp), dimension(nmo,nao) :: b
        intent (in) a, nao, nmo
        intent (out) b

        ! local variables
        integer :: i, j

        do i=1, nmo
            do j=1, nao
                b(i,j)=a(j,i)
            end do
        end do

        return

    end subroutine maths_tcoef


    ! Utility to interface to the BLAS routine DAXPY, which evaluates a
    ! constant (X) times a vector (VECA) plus a vector (VECB).
    subroutine maths_cveca_plus_vecb(n,x,veca,inca,vecb,incb)

        integer :: n, inca, incb
        real(kind=wp) :: x
        real(kind=wp), dimension(*) :: veca, vecb

        call daxpy(n,x,veca,inca,vecb,incb)

    end subroutine maths_cveca_plus_vecb


    ! Routine to emulate CRAY SCILIB routine, MXM
    subroutine maths_mxm(a,nar,b,nac,c,nbc)

        use consts, only : zero=>xzero

        ! dummy arguments
        integer :: nac, nar, nbc
        real(kind=wp), dimension(nar,*) :: a, c
        real(kind=wp), dimension(nac,*) :: b
        intent (in) a, b, nac, nar, nbc
        intent (inout) c

        ! local variables
        integer :: i, j, k

        do j=1, nbc
            do i=1, nar
                c(i,j)=zero
                do k=1, nac
                c(i,j)=c(i,j)+a(i,k)*b(k,j)
                end do
            end do
        end do

        return

    end subroutine maths_mxm

end module maths
