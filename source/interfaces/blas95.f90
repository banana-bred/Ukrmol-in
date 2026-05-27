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

!> \brief   Intel-MKL-like BLAS95 wrappers
!> \authors J Benda
!> \date    2018 - 2020
!>
!> Provides type-generic interfaces for several BLAS subroutines used in UKRmol+.
!> This is only needed when not using Intel MKL (e.g. when using purely free software
!> toolchain). The wrappers are written to behave exactly as Intel MKL BLAS95.
!> An underlying BLAS library (e.g. Netlib BLAS, or OpenBLAS) is needed to
!> build this library.
!>
module blas95

    use blas_lapack_gbl, only: blasint
    use iso_fortran_env, only: real64

    implicit none

    !> \brief Matrix-martrix multiplication
    !>
    !> Multiply two dense matrices:
    !> \f[
    !>      \mathbf{c} \leftarrow \alpha \mathbf{a} \mathbf{b} + \beta \mathbf{c}
    !> \f]
    !>
    !> Only the two input and one output matrices are required, with optional flags that
    !> specify transposition, and linear combination coefficients.
    !>
    !> When the matrices A and B are of incompatible size (i.e. A has more columns than B
    !> has rows, or similarly for transposed matrices), only subset of the larger matrix is
    !> used in the multiplication (keeping the "leading dimension" parameter passed to xGEMM
    !> intact).
    !>
    !> \param[in]  a        First input matrix.
    !> \param[in]  b        Second input matrix.
    !> \param[in]  transa   Transposition state of \c a, default 'N'.
    !> \param[in]  transb   Transposition state of \c b, default 'N'.
    !> \param[in]  alpha    Multiplication factor for a*b, default 1.
    !> \param[in]  beta     Multiplication factor for c, default 0.
    !> \param[out] c        Output matrix.
    !>
    interface gemm
        module procedure dgemm_wrapper
        module procedure zgemm_wrapper
    end interface

    !> \brief Matrix-vector.
    !>
    !> Multiply vector by a matrix.
    !>
    interface gemv
        module procedure dgemv_wrapper
        module procedure zgemv_wrapper
    end interface

contains

    subroutine dgemm_wrapper (a, b, c, transa, transb, alpha, beta)

        real(real64), dimension(:,:), intent(in)  :: a, b
        real(real64), dimension(:,:), intent(out) :: c
        character(len=1), intent(in), optional    :: transa, transb
        real(real64),     intent(in), optional    :: alpha, beta

        real(real64), parameter :: one  = 1
        real(real64), parameter :: zero = 0

        integer(blasint) :: m, n, k, lda, ldb, ldc, OP_a_rows, OP_a_cols, OP_b_rows, OP_b_cols, OP_c_rows, OP_c_cols
        real(real64)     :: alpha1, beta1
        character(len=1) :: transa1, transb1

        transa1 = 'N'
        transb1 = 'N'
        alpha1  = one
        beta1   = zero

        if (present(transa)) transa1 = transa
        if (present(transb)) transb1 = transb
        if (present(alpha))  alpha1  = alpha
        if (present(beta))   beta1   = beta

        OP_a_rows = merge(size(a,1), size(a,2), transa1 == 'N')
        OP_a_cols = merge(size(a,2), size(a,1), transa1 == 'N')
        OP_b_rows = merge(size(b,1), size(b,2), transb1 == 'N')
        OP_b_cols = merge(size(b,2), size(b,1), transb1 == 'N')
        OP_c_rows = size(c,1)
        OP_c_cols = size(c,2)

        m = min(OP_a_rows,OP_c_rows)
        k = min(OP_a_cols,OP_b_rows)
        n = min(OP_b_cols,OP_c_cols)

        lda = size(a,1)
        ldb = size(b,1)
        ldc = size(c,1)

        call dgemm (transa1, transb1, m, n, k, alpha1, a, lda, b, ldb, beta1, c, ldc)

    end subroutine dgemm_wrapper


    subroutine zgemm_wrapper (a, b, c, transa, transb, alpha, beta)

        complex(real64), dimension(:,:), intent(in)  :: a, b
        complex(real64), dimension(:,:), intent(out) :: c
        character(len=1), intent(in), optional       :: transa, transb
        complex(real64), intent(in), optional        :: alpha, beta

        complex(real64), parameter :: one  = 1
        complex(real64), parameter :: zero = 0

        integer(blasint) :: m, n, k, lda, ldb, ldc, OP_a_rows, OP_a_cols, OP_b_rows, OP_b_cols, OP_c_rows, OP_c_cols
        complex(real64)  :: alpha1, beta1
        character(len=1) :: transa1, transb1

        transa1 = 'N'
        transb1 = 'N'
        alpha1  = one
        beta1   = zero

        if (present(transa)) transa1 = transa
        if (present(transb)) transb1 = transb
        if (present(alpha))  alpha1  = alpha
        if (present(beta))   beta1   = beta

        OP_a_rows = merge(size(a,1), size(a,2), transa1 == 'N')
        OP_a_cols = merge(size(a,2), size(a,1), transa1 == 'N')
        OP_b_rows = merge(size(b,1), size(b,2), transb1 == 'N')
        OP_b_cols = merge(size(b,2), size(b,1), transb1 == 'N')
        OP_c_rows = size(c,1)
        OP_c_cols = size(c,2)

        m = min(OP_a_rows,OP_c_rows)
        k = min(OP_a_cols,OP_b_rows)
        n = min(OP_b_cols,OP_c_cols)

        lda = size(a,1)
        ldb = size(b,1)
        ldc = size(c,1)

        call zgemm (transa1, transb1, m, n, k, alpha1, a, lda, b, ldb, beta1, c, ldc)

    end subroutine zgemm_wrapper


    subroutine dgemv_wrapper (a, x, y, alpha, beta, trans)

        real(real64), dimension(:,:), intent(in) :: a
        real(real64), dimension(:), intent(in)   :: x
        real(real64), dimension(:), intent(out)  :: y
        real(real64), intent(in), optional       :: alpha, beta
        character(len=1), intent(in), optional   :: trans

        real(real64), parameter :: one  = 1
        real(real64), parameter :: zero = 0

        integer(blasint) :: m, n, lda, incx, incy
        real(real64)     :: alpha1, beta1
        character(len=1) :: trans1

        trans1 = 'N'
        alpha1 = one
        beta1  = zero

        if (present(trans)) trans1 = trans
        if (present(alpha)) alpha1 = alpha
        if (present(beta))  beta1  = beta

        m = size(a, 1)
        n = size(a, 2)

        lda = merge(m, n, trans1 == 'N')

        incx = 1
        incy = 1

        call dgemv (trans1, m, n, alpha1, a, lda, x, incx, beta1, y, incy)

    end subroutine dgemv_wrapper


    subroutine zgemv_wrapper (a, x, y, alpha, beta, trans)

        complex(real64), dimension(:,:), intent(in) :: a
        complex(real64), dimension(:), intent(in)   :: x
        complex(real64), dimension(:), intent(out)  :: y
        complex(real64), intent(in), optional       :: alpha, beta
        character(len=1), intent(in), optional      :: trans

        complex(real64), parameter :: one  = 1
        complex(real64), parameter :: zero = 0

        integer(blasint) :: m, n, lda, incx, incy
        complex(real64)  :: alpha1, beta1
        character(len=1) :: trans1

        trans1 = 'N'
        alpha1 = one
        beta1  = zero

        if (present(trans)) trans1 = trans
        if (present(alpha)) alpha1 = alpha
        if (present(beta))  beta1  = beta

        m = size(a, 1)
        n = size(a, 2)

        lda = merge(m, n, trans1 == 'N')

        incx = 1
        incy = 1

        call zgemv (trans1, m, n, alpha1, a, lda, x, incx, beta1, y, incy)

    end subroutine zgemv_wrapper

end module blas95
