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

!> \brief   Intel-MKL-like LAPACK95 wrappers
!> \authors J Benda
!> \date    2018 - 2020
!>
!> Provides type-generic interfaces for several LAPACK subroutines used in UKRmol+.
!> This is only needed when not using Intel MKL (e.g. when using purely free software
!> toolchain). The wrappers are written to behave exactly as Intel MKL LAPACK95.
!> An underlying LAPACK library (e.g. Netlib LAPACK, or OpenBLAS) is needed to
!> build this library.
!>
module lapack95

    use blas_lapack_gbl, only: blasint
    use iso_fortran_env, only: real64

    implicit none

    !> \brief LU decomposition
    !>
    !> xGETRF computes an LU factorization of a general M-by-N matrix A
    !> using partial pivoting with row interchanges.
    !>
    !> The factorization has the form
    !>     A = P * L * U
    !> where P is a permutation matrix, L is lower triangular with unit
    !> diagonal elements (lower trapezoidal if m > n), and U is upper
    !> triangular (upper trapezoidal if m < n).
    !>
    !> This is the right-looking Level 3 BLAS version of the algorithm.
    !>
    interface getrf
        module procedure dgetrf_wrapper
        module procedure zgetrf_wrapper
    end interface getrf


    !> \brief LU back-substitution
    !>
    !> xGETRS solves a system of linear equations
    !>     A * X = B  or  A**T * X = B
    !> with a general N-by-N matrix A using the LU factorization computed
    !> by xGETRF.
    !>
    interface getrs
        module procedure dgetrs_wrapper
        module procedure zgetrs_wrapper
    end interface getrs


    !> \brief Iterative refinement of a solution
    !>
    !> xGERFS improves the computed solution to a system of linear
    !> equations and provides error bounds and backward error estimates for
    !> the solution.
    !>
    interface gerfs
        module procedure dgerfs_wrapper
    end interface gerfs


    !> \brief Get condition number
    !>
    !> xGECON estimates the reciprocal of the condition number of a general
    !> real matrix A, in either the 1-norm or the infinity-norm, using
    !> the LU factorization computed by xGETRF.
    !>
    !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
    !> condition number is computed as
    !>     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
    !>
    interface gecon
        module procedure dgecon_wrapper
    end interface gecon


    !> \brief Reduction to bi-diagonal matrix
    !>
    !> xGEBRD reduces a general real M-by-N matrix A to upper or lower
    !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
    !>
    !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    !>
    interface gebrd
        module procedure dgebrd_wrapper
    end interface gebrd


    !> \brief Singular value decomposition of bi-diagonal matrix
    !>
    !> xBDSQR computes the singular values and, optionally, the right and/or
    !> left singular vectors from the singular value decomposition (SVD) of
    !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
    !> zero-shift QR algorithm.  The SVD of B has the form
    !>
    !>     B = Q * S * P**T
    !>
    !> where S is the diagonal matrix of singular values, Q is an orthogonal
    !> matrix of left singular vectors, and P is an orthogonal matrix of
    !> right singular vectors.  If left singular vectors are requested, this
    !> subroutine actually returns U*Q instead of Q, and, if right singular
    !> vectors are requested, this subroutine returns P**T*VT instead of
    !> P**T, for given real input matrices U and VT.  When U and VT are the
    !> orthogonal matrices that reduce a general matrix A to bidiagonal
    !> form:  A = U*B*VT, as computed by xGEBRD, then
    !>
    !>     A = (U*Q) * S * (P**T*VT)
    !>
    !> is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
    !> for a given real input matrix C.
    !>
    !> See "Computing  Small Singular Values of Bidiagonal Matrices With
    !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
    !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
    !> no. 5, pp. 873-912, Sept 1990) and
    !> "Accurate singular values and differential qd algorithms," by
    !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
    !> Department, University of California at Berkeley, July 1992
    !> for a detailed description of the algorithm.
    !>
    interface bdsqr
        module procedure dbdsqr_wrapper
    end interface bdsqr


    !> \brief Singular value decomposition
    !>
    !> xGESVD computes the singular value decomposition (SVD) of a real
    !> M-by-N matrix A, optionally computing the left and/or right singular
    !> vectors. The SVD is written
    !>
    !>     A = U * SIGMA * transpose(V)
    !>
    !> where SIGMA is an M-by-N matrix which is zero except for its
    !> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
    !> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
    !> are the singular values of A; they are real and non-negative, and
    !> are returned in descending order.  The first min(m,n) columns of
    !> U and V are the left and right singular vectors of A.
    !>
    !> Note that the routine returns V**T, not V.
    !>
    interface gesvd
        module procedure dgesvd_wrapper
    end interface gesvd

contains

    subroutine dgetrf_wrapper (a, ipiv, info)

        real(real64), dimension(:,:), intent(inout)  :: a
        integer(blasint), intent(out), optional      :: ipiv(:)
        integer(blasint), intent(out), optional      :: info

        integer(blasint)                       :: m, n, lda, info1
        integer(blasint), dimension(size(a,1)) :: ipiv1

        m = size(a, 1)
        n = size(a, 2)

        lda = m

        call dgetrf (m, n, a, lda, ipiv1, info1)

        if (present(info)) info = info1
        if (present(ipiv)) ipiv = ipiv1

    end subroutine dgetrf_wrapper


    subroutine zgetrf_wrapper (a, ipiv, info)

        complex(real64), dimension(:,:), intent(inout) :: a
        integer(blasint), intent(out), optional        :: ipiv(:)
        integer(blasint), intent(out), optional        :: info

        integer(blasint)                       :: m, n, lda, info1
        integer(blasint), dimension(size(a,1)) :: ipiv1

        m = size(a, 1)
        n = size(a, 2)

        lda = m

        call zgetrf (m, n, a, lda, ipiv1, info1)

        if (present(info)) info = info1
        if (present(ipiv)) ipiv = ipiv1

    end subroutine zgetrf_wrapper


    subroutine dgetrs_wrapper (a, ipiv, b, trans, info)

        real(real64), dimension(:,:), intent(in)    :: a
        real(real64), dimension(:,:), intent(inout) :: b
        integer(blasint), intent(in)                :: ipiv(:)
        character(len=1), intent(in), optional      :: trans
        integer(blasint), intent(out), optional     :: info

        integer(blasint), allocatable :: ipiv1(:)

        integer(blasint) :: n, nrhs, lda, ldb, info1
        character(len=1) :: trans1

        trans1 = 'N'
        ipiv1 = ipiv

        if (present(trans)) trans1 = trans

        n    = size(a, 1)
        nrhs = size(b, 2)

        lda = n
        ldb = n

        call dgetrs(trans1, n, nrhs, a, lda, ipiv1, b, ldb, info1)

        if (present(info)) info = info1

    end subroutine dgetrs_wrapper


    subroutine zgetrs_wrapper (a, ipiv, b, trans, info)

        complex(real64), dimension(:,:), intent(in)    :: a
        complex(real64), dimension(:,:), intent(inout) :: b
        integer(blasint), dimension(:), intent(in)     :: ipiv
        character(len=1), intent(in), optional         :: trans
        integer(blasint), intent(out), optional        :: info

        integer(blasint), allocatable :: ipiv1(:)

        integer(blasint) :: n, nrhs, lda, ldb, info1
        character(len=1) :: trans1

        trans1 = 'N'
        ipiv1 = ipiv

        if (present(trans)) trans1 = trans

        n    = size(a, 1)
        nrhs = size(b, 2)

        lda = n
        ldb = n

        call zgetrs(trans1, n, nrhs, a, lda, ipiv1, b, ldb, info1)

        if (present(info)) info = info1

    end subroutine zgetrs_wrapper


    subroutine dgerfs_wrapper (a, af, ipiv, b, x, trans, ferr, berr, info)

        real(real64),     intent(in)    :: a(:,:), af(:,:), b(:,:)
        real(real64),     intent(inout) :: x(:,:)
        integer(blasint), intent(in)    :: ipiv(:)

        real(real64),     intent(inout), optional :: ferr(:), berr(:)
        character(len=1), intent(in),    optional :: trans
        integer(blasint), intent(out),   optional :: info

        real(real64), allocatable :: work(:), ferr1(:), berr1(:)
        integer(blasint), allocatable :: iwork(:), ipiv1(:)

        character(len=1) :: trans1
        integer(blasint) :: n, nrhs, info1, lda, ldaf, ldb, ldx

        trans1 = 'N'
        ipiv1 = ipiv

        if (present(trans)) trans1 = trans

        n    = size(ipiv)
        nrhs = size(b, 2)

        lda  = n
        ldaf = n
        ldb  = n
        ldx  = n

        allocate (work(n), iwork(n), ferr1(n), berr1(n))

        call dgerfs(trans1, n, nrhs, a, lda, af, ldaf, ipiv1, b, ldb, x, ldx, ferr1, berr1, work, iwork, info1)

        if (present(info)) info = info1
        if (present(ferr)) ferr = ferr1
        if (present(berr)) berr = berr1

    end subroutine dgerfs_wrapper


    subroutine dgecon_wrapper (a, anorm, rcond, norm, info)

        real(real64),     intent(in)    :: a(:,:)
        real(real64),     intent(out)   :: anorm, rcond

        character(len=1), intent(in),  optional :: norm
        integer(blasint), intent(out), optional :: info

        real(real64), allocatable :: work(:)
        integer(blasint), allocatable :: iwork(:)

        integer(blasint) :: n, lda, info1, one = 1
        character(len=1) :: norm1

        norm1 = '1'
        n     = size(a, 1)
        lda   = max(one, n)

        if (present(norm)) norm1 = norm

        allocate (work(4*n), iwork(n))

        call dgecon (norm1, n, a, lda, anorm, rcond, work, iwork, info1)

        if (present(info)) info = info1

    end subroutine dgecon_wrapper


    subroutine dgebrd_wrapper (a, d, e, tauq, taup, info)

        real(real64), intent(inout) :: a(:,:)
        real(real64), intent(inout), optional :: d(:), e(:), tauq(:), taup(:)
        integer(blasint), intent(out), optional :: info

        real(real64), allocatable :: d1(:), e1(:), tauq1(:), taup1(:), work(:)

        integer(blasint) :: info1, m, n, sz, lwork, lda, one = 1

        m = size(a, 1)
        n = size(a, 2)

        lda   = max(one, m)
        sz    = max(one, min(m, n))
        lwork = max(one, max(m, n))

        allocate(d1(sz), e1(sz), tauq1(sz), taup1(sz), work(lwork))

        call dgebrd (m, n, a, lda, d1, e1, tauq1, taup1, work, lwork, info1)

        if (present(info)) info = info1

        if (min(m, n) > 0) then
            if (present(d))    d = d1
            if (present(e))    e = e1
            if (present(tauq)) tauq = tauq1
            if (present(taup)) taup = taup1
        end if

    end subroutine dgebrd_wrapper


    subroutine dbdsqr_wrapper (d, e, vt, u, c, uplo, info)

        real(real64),     intent(inout)           :: d(:), e(:)
        integer(blasint), intent(out),   optional :: info
        real(real64),     intent(inout), optional :: vt(:,:), u(:,:), c(:,:)
        character(len=1), intent(in),    optional :: uplo

        real(real64), allocatable :: work(:)
        integer(blasint) :: info1, n, ldvt, nru, ncc, ldu, ldc, ncvt, one = 1
        character(len=1) :: uplo1

        n = size(d)

        uplo1 = 'N'
        ldvt  = 1
        ncvt  = 0
        nru   = 0
        ncc   = 0

        if (present(uplo)) uplo1 = uplo
        if (present(vt)) ldvt = size(vt, 1)
        if (present(vt)) ncvt = size(vt, 2)
        if (present(u))  nru  = size(u,  1)
        if (present(c))  ncc  = size(c,  2)

        ldu = max(one, nru)
        ldc = max(one, n)

        allocate (work(max(1, 4*n)))

        call dbdsqr (uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info1)

        if (present(info)) info = info1

    end subroutine dbdsqr_wrapper


    subroutine dgesvd_wrapper (a, s, u, vt, ww, job, info)

        real(real64), intent(inout) :: a(:,:), s(:)

        real(real64),     intent(inout), optional :: u(:,:), vt(:,:), ww(:,:)
        character(len=1), intent(in),    optional :: job
        integer(blasint), intent(out),   optional :: info

        real(real64), allocatable :: work(:), ww1(:,:)
        integer(blasint) :: info1, m, n, lda, ldu, ldvt, lwork, one = 1
        character(len=1) :: jobu, jobvt

        m = size(a, 1)
        n = size(a, 2)

        lda  = max(one, m)
        ldu  = max(one, m)
        ldvt = max(one, n)

        lwork = max(3 * min(m, n) + max(m, n), 5 * min(m, n))

        allocate (work(max(one, lwork)))

        if (present(u))  then ; jobu  = 'A' ; else if (job == 'U') then ; jobu  = 'O' ; else ; jobu  = 'N' ; end if
        if (present(vt)) then ; jobvt = 'A' ; else if (job == 'V') then ; jobvt = 'O' ; else ; jobvt = 'N' ; end if

        call dgesvd (jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info1)

        if (present(info)) info = info1

        ! TODO: ww

    end subroutine dgesvd_wrapper

end module lapack95
