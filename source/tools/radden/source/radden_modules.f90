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
MODULE constants
use precisn
IMPLICIT NONE
INTEGER, PARAMETER :: line_len=132 !length of one line

REAL(kind=wp), PARAMETER :: pi=3.14159265358979323_wp
REAL(kind=wp), PARAMETER :: to_au=1.88972599_wp !Conversion constant from angstroms to a.u. This is the value which Molpro uses when writing the Molden file.
REAL(kind=wp), PARAMETER :: r_inf=20.0_wp       !Radius of the sphere in which the basis functions are integrated numerically when checking their norms.

CHARACTER(len=5),PARAMETER :: BFTYPE(1:2)=(/'SA-BF','CGTO '/) !two possible headers of the density matrix input file saying which basis set has to be used.

END MODULE constants

MODULE prefactors
!This module declares angular factors used in definitions of primitive gaussians. These definitions are compatible with the latest (4.8 version) of Molden.
IMPLICIT NONE
CHARACTER(len=4) :: ang_fact(1:5,1:15) !definitions of the angular parts of gaussians. The order of the angular factors corresponds to the Molden format.
CHARACTER(len=1), PARAMETER :: bas_typ(1:5)=(/'s','p','d','f','g'/)  !names of angular types of basis functions
INTEGER, PARAMETER :: bf_num(1:5)=(/1,3,6,10,15/)                    !number of cartesian gaussians corresponding to each angular type

DATA ang_fact(1,1)      /'s'/
DATA ang_fact(2:2,1:3)  /'x','y','z'/
DATA ang_fact(3:3,1:6)  /'xx','yy','zz','xy','xz','yz'/
DATA ang_fact(4:4,1:10) /'xxx','yyy','zzz','xyy','xxy','xxz','xzz','yzz','yyz','xyz'/
DATA ang_fact(5:5,1:15) /'xxxx','yyyy','zzzz','xxxy','xxxz','yyyx','yyyz','zzzx','zzzy','xxyy','xxzz','yyzz','xxyz','yyxz','zzxy'/

END MODULE prefactors

MODULE saved_vars 
!the variables declared here are used by routines integrating the basis functions, i.e. routines in the module 'basis_fns_num'
use constants
use precisn
IMPLICIT NONE

SAVE

REAL(kind=wp) :: r, theta, phi  !cooridnates used during integration of the basis functions
REAL(kind=wp) :: prev_s_th, prev_c_th !the most recent values of sin(theta) and cos(theta): saves time in evaluating GET_X, Y, Z functions
INTEGER :: prim_i               !index of the primitive basis function to be evaluated
INTEGER :: contr_i,contr_j      !index of the contracted basis function to be evaluated
INTEGER :: orb_i                !index of the orbital whose density we want to calculate
INTEGER :: sabf_i,sabf_j        !index of the symmetry adapted basis function to be evaluated

LOGICAL :: iscgto               !which basis we have - cgto or symmetry adapted basis. This is used in integration of densities/orbitals.

END MODULE saved_vars

MODULE basis_atoms
!This module declares the two most important variables 'basis' and 'orbital' which hold the complete basis set information 
!and information about the orbitals and (or) density matrices processed. These two variables are 'saved' therefore modifiable
!by any surboutine which uses this module. Therefore any subroutine with 'use basis_atoms' either reads some information from
!these two variables or it modifies them, e.g. when reading the basis set or the orbitals from a file.
use constants
use precisn

IMPLICIT NONE

INTEGER, PARAMETER :: max_at  = 100   !maximum number of atoms allowed
INTEGER, PARAMETER :: maxsym  = 8     !maximum number of symmetries (used only when the R-matrix symmetry-adapted basis is used)
INTEGER, PARAMETER :: max_orb = 1000  !maximum number of orbitals allowed

TYPE basis_type

    INTEGER :: ncgto  !number of contracted gaussians
    INTEGER :: nprim  !number of primary gaussians
    INTEGER :: nadapt !number of symmetry contractions
    INTEGER :: maxall !maximum number of contractions used in one contracted function (needed for optimal allocations of some arrays)
    INTEGER :: no_at  !number of atoms
    INTEGER :: bind   !number of lines in the molden file corresponding to basis set information
    
    REAL(kind=wp), ALLOCATABLE :: alpha(:)          !(nprim) exponents of primitive gaussians
    REAL(kind=wp), ALLOCATABLE :: contractions(:,:) !(ncgto,maxall) contraction coefficients (function:coefficients)
    REAL(kind=wp), ALLOCATABLE :: at_coord(:,:)     !(3,no_at) coordinates of atoms (coordinate 1-3,atom)
    REAL(kind=wp), ALLOCATABLE :: overlaps_cgto(:,:)!(ncgto,ncgto) overlap integrals of contracted gaussians as calculated analytically
    
    INTEGER, ALLOCATABLE :: at_num(:)        !(no_at) indices of atoms
    INTEGER, ALLOCATABLE :: at_ind(:)        !(nprim) number of an atom at which a primitive gaussian function is centred
    INTEGER, ALLOCATABLE :: contr_alpha(:,:) !(ncgto,maxall) index to an exponent of a prim. gaussian in the array 'alpha' corresponding to a given 
                                             !contraction coefficient in a contr. basis function
    INTEGER, ALLOCATABLE :: ang_exp(:,:)     !(nprim,3) holds exponents i,j,k of prefactors x^i*y^j*z^k (angular behaviours of primitive GTOs)
    INTEGER, ALLOCATABLE :: no_contr(:)      !(ncgto) number of primitives building up a contraction
    
    CHARACTER(len=4), ALLOCATABLE :: at_sym(:)  !(no_at) atomic symbols
    CHARACTER(len=4), ALLOCATABLE :: prefact(:) !(nprim) contains angular mom. factor symbol, e.g. 'xx', 'xy', 'xxx', etc. for all primitive gaussians

!!!!additional variables used only when the R-matrix basis (symmetry-adapted) is used:
    INTEGER :: nsabf                               !total number of symmetry adapted basis functions
    INTEGER :: nsymt                               !number of symmetries
    INTEGER :: nbft(1:maxsym)                      !number of basis functions per symmetry
    INTEGER, ALLOCATABLE :: no_sabf(:)             !(nsabf) number of contractions building up a given symmetry adapted function
    REAL(kind=wp), ALLOCATABLE :: sa_contr(:,:)    !(nsabf,ncgto) coefficients of contracted functions building up a given symmetry adapted function
    INTEGER, ALLOCATABLE :: contr_ind(:,:)         !(nsabf,ncgto) index of a contracted gaussian corresponding to a given coefficient 
                                                   !in a symmetry adapted basis function
!!!!

!!!! Spherical gaussian basis information, i.e. GTO(r-A) = Norm*r(r-A)^{l}*Y(l,m)*exp(-alpha*(r-A)^2); m = 2*l+1
    INTEGER :: sph_ncgto  !number of contracted gaussians
    INTEGER :: sph_nprim  !number of primary gaussians
    INTEGER :: sph_maxall !maximum number of contractions used in one contracted function (needed for optimal allocations of some arrays)

    REAL(kind=wp), ALLOCATABLE :: sph_a(:)          !(sph_nprim) exponents of primitive gaussians
    REAL(kind=wp), ALLOCATABLE :: sph_contr(:,:)    !(sph_ncgto,maxall) contraction coefficients (function:coefficients)
    REAL(kind=wp), ALLOCATABLE :: sph_olap_cgto(:,:)!(sph_ncgto,sph_ncgto) overlap integrals of contracted gaussians as calculated analytically

    INTEGER, ALLOCATABLE :: sph_at_ind(:)    !(sph_nprim) number of an atom at which a primitive gaussian function is centred
    INTEGER, ALLOCATABLE :: sph_contr_a(:,:) !(sph_ncgto,maxall) index to an exponent of a prim. gaussian in the array 'alpha' corresponding to a given 
                                             !contraction coefficient in a contr. basis function
    INTEGER, ALLOCATABLE :: l_ang(:), m_ang(:) !(sph_nprim) 'l' and 'm' angular parts of the primitive gaussians
    INTEGER, ALLOCATABLE :: sph_no_contr(:)  !(sph_ncgto) number of primitives building up a contraction

END TYPE basis_type

TYPE orbital_type

     REAL(kind=wp), ALLOCATABLE :: coeff(:) !coefficients of primitive baussians for this orbital
     REAL(kind=wp), ALLOCATABLE :: dm(:,:)  !density matrix of this orbital
     REAL(kind=wp) :: sym                   !symmetry specification: 'number.symmetry' NOT IMPLEMENTED YET (provided by Molpro from v. 2010.something)
     REAL(kind=wp) :: ene                   !energy of the orbital [H]                 NOT IMPLEMENTED YET
     CHARACTER(len=5) :: spin               !string defining spin ('Alpha' or 'Beta')  NOT IMPLEMENTED YET
     REAL(kind=wp) :: occup                 !occupation number                         NOT IMPLEMENTED YET
     REAL(kind=wp) :: tot_ch                !total charge in the orbital (as calculated analytically during construction of the density matrix)

END TYPE orbital_type

SAVE 
!these variables need to be accessible by subroutines which enter integration of the orbitals/density matrices

TYPE(basis_type) :: basis 
TYPE(orbital_type), ALLOCATABLE :: orbital(:)

END MODULE basis_atoms

MODULE el_fns
!This module contains some elementray functions like factorials, mapping from spherical to cartesian coordinates, etc.

CONTAINS

  FUNCTION f_cf(k,l,m,a,b)
  use precisn
  !computes the coefficient 'f_k' from the formula (2.5) from the Taketa paper; f_cf(k,l,m,a,b) = f_k(l,m,a,b) from the Taketa paper
  IMPLICIT NONE
  REAL(kind=wp) :: f_cf, a, b
  INTEGER :: k, l, m, j, p, q
  
  f_cf=0.0_wp
  DO p=0,l
     DO q=0,m
        IF (p+q .eq. k) THEN  !this is the straightforward way
           f_cf=f_cf+Binom(l,p)*Binom(m,q)*a**(l-p)*b**(m-q) 
        END IF
     END DO
  END DO

  END FUNCTION f_cf

  FUNCTION Binom(m,n)
  use precisn
  !binomial coefficient m over n
  IMPLICIT NONE
  REAL(kind=wp) :: Binom
  INTEGER :: n,m

  Binom=FACTORIAL(m)/(FACTORIAL(n)*FACTORIAL(m-n))

  END FUNCTION Binom

  FUNCTION FACTORIAL(f)
  use precisn
  INTEGER :: f, FACTORIAL, i

  FACTORIAL=1

  IF (f < 0) THEN
     RETURN
  END IF

  DO i=2,f
     FACTORIAL=FACTORIAL*i
  END DO

  END FUNCTION FACTORIAL

  FUNCTION D_FACTORIAL(f)
  use precisn
  !double factorial
  IMPLICIT NONE
  INTEGER :: f, D_FACTORIAL, i

  D_FACTORIAL=1

  IF (f < 0) THEN
     RETURN
  END IF

  DO i=1,f,2
     D_FACTORIAL=D_FACTORIAL*i
  END DO

  END FUNCTION D_FACTORIAL

  FUNCTION GET_X(rd,ph)
  !BEWARE: this function is optimized for calculation of the orbital densities:
  !it uses the saved variable prev_s_th. It is therefore NOT a general spherical
  !-> cartesian coordinates mapping function!!
  use constants
  use precisn
  use saved_vars
  IMPLICIT NONE
  REAL(kind=wp) :: GET_X, rd, ph

  GET_X=rd*prev_s_th*cos(ph) !rd*sin(th)*cos(ph)
  
  END FUNCTION GET_X

  FUNCTION GET_Y(rd,ph)
  !BEWARE: this function is optimized for calculation of the orbital densities:
  !it uses the saved variable prev_s_th. It is therefore NOT a general spherical
  !-> cartesian coordinates mapping function!!
  use constants
  use precisn
  use saved_vars
  IMPLICIT NONE
  REAL(kind=wp) :: GET_Y, rd, ph

  GET_Y=rd*prev_s_th*sin(ph) !sin(th)*sin(ph)

  END FUNCTION GET_Y

  FUNCTION GET_Z(rd,ph)
  !BEWARE: this function is optimized for calculation of the orbital densities:
  !it uses the saved variable prev_c_th. It is therefore NOT a general spherical
  !-> cartesian coordinates mapping function!!
  use constants
  use precisn
  use saved_vars
  IMPLICIT NONE
  REAL(kind=wp) :: GET_Z, rd, ph

  GET_Z=rd*prev_c_th !cos(th)

  END FUNCTION GET_Z

  FUNCTION R_SQ(x,y,z)
  use precisn
  IMPLICIT NONE
  REAL(kind=wp) :: x, y, z, R_SQ

  R_SQ = x*x + y*y + z*z  

  END FUNCTION R_SQ

  FUNCTION U_mu_l_m(mu,l,m)
  !calculates elements of the unitary matrix transforming the complex spherical harmonics to the real spherical harmonics. See below:
  !The real spherical harmonics are defined as: X_{l}^{mu}=sum_{m}U^{mu}_{l,m}Y_{l}^{m}, where Y_{l}^{m} are the complex spherical harmonics
  !the matrix U^{mu}_{l,m} is defined here according to the expression (12) in: 
  !Homeier, H.H.H. & Steinborn, E.O. Some properties of the coupling coefficients of real spherical harmonics and their relation to Gaunt coefficients. 
  !Journal of Molecular Structure: THEOCHEM 368, 31-37 (1996).
  !The matrix U^{mu}_{l,m} corresponds to the usual transformation from complex to real spherical harmonics
  use precisn
  IMPLICIT NONE
  COMPLEX(kind=wp) :: U_mu_l_m
  INTEGER, INTENT(IN) :: mu, l, m
  INTEGER :: d_m_mu, d_m_mmu, H_mu, H_mmu !delta_{m,mu}, delta_{m,-mu}, Heaviside(mu), Heaviside(-mu)
  COMPLEX(kind=wp), PARAMETER :: imu = CMPLX(0d0,1d0) !0 + i*1

  IF (abs(m) > l .or. abs(mu) > l) STOP "Invalid value of mu or m passed to the function U_mu_l_m! Program terminated."

  IF (abs(mu) .ne. abs(m)) THEN !complex spherical harmonics with different 'm' don't contribute
     U_mu_l_m = 0d0
     RETURN
  ENDIF

  IF (m == 0 .and. mu == 0) THEN !X_{l}^{0} = Y_{l}^{0}
     U_mu_l_m = 1d0
     RETURN
  ENDIF

  d_m_mu = 0
  d_m_mmu= 0
  IF (m == mu)  d_m_mu=1
  IF (m == -mu) d_m_mmu=1
  
  H_mu = 0
  H_mmu = 0
  IF (mu > 0) THEN
     H_mu = 1
  ELSE !mu .le. 0
     H_mmu = 1
  ENDIF

  U_mu_l_m = 1d0/sqrt(2d0)*(H_mu*d_m_mu + H_mmu*imu*d_m_mu*(-1)**m + H_mmu*(-imu)*d_m_mmu + H_mu*d_m_mmu*(-1)**m)

  END FUNCTION U_mu_l_m

  FUNCTION C_G_coeff(l1,m1,l2,m2,L,M)
  !calculates the usual Clebsch-Gordan coefficient <l1,m1,l2,m2|L,M>. All arguments = WHOLE NUMBERS. Allowance for half-integer arguments can be achieved
  !simply by adjusting the 'factorial' function.
  !the formulas: Edmonds - Angular momentum in QM
  use precisn
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l1,m1,l2,m2,L,M
  REAL(kind=wp) :: C_G_coeff
  REAL(kind=wp) :: ro, sigma, tau
  INTEGER :: r, r_min, r_max

  C_G_coeff = 0d0

  IF (m1 + m2 .ne. M) RETURN
  IF (L < abs(l1-l2) .or. L > l1+l2) RETURN

  ro = sqrt((2d0*L+1d0)*factorial(l1+l2-L)*factorial(L+l1-l2)*factorial(l2+L-l1)/factorial(l1+l2+L+1))
  sigma=sqrt(1d0*factorial(L+M)*factorial(L-M)*factorial(l1+m1)*factorial(l1-m1)*factorial(l2+m2)*factorial(l2-m2))

  r_min=MAX(0,l1+l2-L+1-(l1+m1+1),l1+l2-L+1-(l2-m2+1))
  r_max=MIN(l1+l2-L,l1-m1,l2+m2)

  tau=0d0 
  DO r=r_min,r_max
     tau = tau + ((-1d0)**r) / (factorial(l1-m1-r) * factorial(l2+m2-r) * factorial(L-l2+m1+r) * factorial(L-l1-m2+r) &
                                * factorial(l1+l2-L-r) * factorial(r))
  ENDDO

  C_G_coeff = ro*sigma*tau

  END FUNCTION C_G_coeff

  FUNCTION W_3j(j1,j2,j3,m1,m2,m3)
  !calculates the Wigner 3j symbol using the usual C-G coefficients. All arguments = WHOLE NUMBERS.
  !Checked with Mathematica, that it gives correct values.
  !For real spherical harmonics the C_G_coeff needs to be replaced by R_G_coeff.
  use precisn
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j1,j2,j3,m1,m2,m3
  REAL(kind=wp) :: W_3j

  W_3j = (-1)**(j1-j2-m3)/sqrt(2d0*j3+1)*C_G_coeff(j1,m1,j2,m2,j3,-m3)

  END FUNCTION

  FUNCTION W_6j(j1,j2,j3,j4,j5,j6)
  !calculates the Wigner's 6-j symbol through calculation of the Racah W-coefficients.
  !DOES NOT ALWAYS RETURN THE CORRECT VALUE! The problem lies in the calculation of the sign in the sum over z. The absolute values of the terms are fine.
  use precisn
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j1,j2,j3,j4,j5,j6
  REAL(kind=wp) :: W_6j
  INTEGER :: z, z_min, z_max, a1, a2, a3, a4, b1, b2, b3, b4
  REAL(kind=wp) :: w

  a1=j1+j2+j3
  a2=j5+j4+j3
  a3=j1+j5+j6
  a4=j2+j4+j6
  b1=j1+j2+j5+j4
  b2=j1+j4+j3+j6
  b3=j2+j5+j3+j6

  z_min=MAX(0,a1,a2,a3,a4)
  z_max=MIN(b1,b2,b3)

  w=0d0
  DO z=z_min,z_max
     w = w + (-1d0)**(z+b1) * factorial(z+1) / (factorial(z-a1) * factorial(z-a2) * factorial(z-a3) * factorial(z-a4) &
                                                 * factorial(b1-z) * factorial(b2-z) * factorial(b3-z))
  ENDDO
  PRINT *,z_min,z_max

  W_6j=delta(j1,j2,j3)*delta(j5,j4,j3)*delta(j1,j5,j6)*delta(j2,j4,j6)*w
  PRINT *,W_6j,delta(j1,j2,j3),delta(j5,j4,j3),delta(j1,j5,j6),delta(j2,j4,j6)

  END FUNCTION W_6j

  FUNCTION delta(a,b,c)
  !auxiliary function used in the calculation of the Wigner's 6_j symbol. CAN CAUSE PROBLEMS FOR HIGH VALUES OF THE ARGUMENTS!!
  use precisn
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a,b,c
  REAL(kind=wp) :: delta

  delta = sqrt(1d0*factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1))

  END FUNCTION delta

  FUNCTION G_coeff(l1,m1,l2,m2,l3,m3)
  !calculates the Gaunt coefficient, i.e. the coupling coefficient for complex spherical harmonics. All arguments = WHOLE NUMBERS.
  !the formula: Edmonds - Angular momentum in QM
  use precisn
  use constants, only: pi
  INTEGER, INTENT(IN) :: l1,m1,l2,m2,l3,m3
  REAL(kind=wp) :: G_coeff

  G_coeff = sqrt((2d0*l1+1)*(2d0*l2+1)*(2d0*l3+1)/(4d0*pi))*C_G_coeff(l1,0,l2,0,l3,0)*C_G_coeff(l1,m1,l2,m2,l3,m3)

  END FUNCTION G_coeff

  FUNCTION R_G_coeff(l1,m1,l2,m2,l3,m3)
  !Formula (27) of the Homeier et al. paper. See subroutine U_mu_l_m for full citation.
  !calculates the coupling coefficients of real spherical harmonics (or R-Gaunt coefficients). All arguments = WHOLE NUMBERS.
  !these are defined as <l1,m1,l2,m2,l3,m3>=int_{theta,phi} X_{l3,m3}*X_{l2,m2}*X_{l1,m1}, where X_{l,m} are REAL SPHERICAL HARMONICS as defined in
  !the subroutine U_mu_l_m. The R-Gaunt coefficients are expressed in terms of the Gaunt coefficients as defined for the complex spherical harmonics.
  use precisn
  INTEGER, INTENT(IN) :: l1,m1,l2,m2,l3,m3
  REAL(kind=wp) :: R_G_coeff
  INTEGER :: mt2, mt3, b2, b3

  !formula (27) of the Homeier et al. paper. See subroutine U_mu_l_m for full citation.

  !before the loop itself there might be the IF statements similar to the ones used in the C_G_coeff routine. Additional selection rules exist.

  R_G_coeff = 0d0

  b3=l3
  b2=max(l1,l2)

  DO mt2=-b2,b2
     DO mt3=-b3,b3
        R_G_coeff = R_G_coeff + real(conjg(U_mu_l_m(m1,l1,mt2+mt3)) * U_mu_l_m(m2,l2,mt2) * U_mu_l_m(m3,l3,mt3)) &
                                 * G_coeff(l2,mt2,l3,mt3,l1,mt2+mt3)
     ENDDO
  ENDDO

  END FUNCTION R_G_coeff

END MODULE el_fns

MODULE basis_fns
!This module contains routines for analytic normalization of the basis functions and for evaluation of basis functions at an arbitrary point (x,y,z)

CONTAINS

  FUNCTION OVERLAP_CGTO(c_i,c_j)
  !calculates analytically the overlap integral between contracted gaussians with indices 'c_i' and 'c_j'
  use precisn
  use basis_atoms
  use el_fns, only: f_cf, D_FACTORIAL, R_SQ
  !INPUT: the common variable 'basis', indices 'c_i' and 'c_j' of contracted gaussians
  IMPLICIT NONE
  INTEGER :: i, j, exp_x_i,exp_x_j,exp_y_i,exp_y_j,exp_z_i,exp_z_j,c_i,c_j, t, u, v
  REAL(kind=wp) :: tmp, OVERLAP_CGTO, gam, px, py, pz, ABsq, ol, f_x,f_y,f_z, ni, nj, n

  tmp=0.0_wp                    !this holds the final overlap integral of the contracted gaussians
  DO i=1,basis%no_contr(c_i)    !over all primitives of the contraction 'c_i'
     DO j=1,basis%no_contr(c_j) !over all primitives of the contraction 'c_j'
           !now calculate the overlap integral for this pair of primitives
          
           gam=basis%alpha(basis%contr_alpha(c_i,i))+basis%alpha(basis%contr_alpha(c_j,j))  !exponent of the product gaussian

           !exponents of the angular factors of the current primitive gaussians
           exp_x_i=basis%ang_exp(basis%contr_alpha(c_i,i),1)
           exp_y_i=basis%ang_exp(basis%contr_alpha(c_i,i),2)
           exp_z_i=basis%ang_exp(basis%contr_alpha(c_i,i),3)
           exp_x_j=basis%ang_exp(basis%contr_alpha(c_j,j),1)
           exp_y_j=basis%ang_exp(basis%contr_alpha(c_j,j),2)
           exp_z_j=basis%ang_exp(basis%contr_alpha(c_j,j),3)

           !centre of the product gaussian
           px=(basis%alpha(basis%contr_alpha(c_i,i))*basis%at_coord(1,basis%at_ind(basis%contr_alpha(c_i,i)))+&
             &basis%alpha(basis%contr_alpha(c_j,j))*basis%at_coord(1,basis%at_ind(basis%contr_alpha(c_j,j))))/gam
           py=(basis%alpha(basis%contr_alpha(c_i,i))*basis%at_coord(2,basis%at_ind(basis%contr_alpha(c_i,i)))+&
             &basis%alpha(basis%contr_alpha(c_j,j))*basis%at_coord(2,basis%at_ind(basis%contr_alpha(c_j,j))))/gam
           pz=(basis%alpha(basis%contr_alpha(c_i,i))*basis%at_coord(3,basis%at_ind(basis%contr_alpha(c_i,i)))+&
             &basis%alpha(basis%contr_alpha(c_j,j))*basis%at_coord(3,basis%at_ind(basis%contr_alpha(c_j,j))))/gam

           !square of distance between the centers of the contracted gaussians
           ABsq = R_SQ(basis % at_coord(1, basis % at_ind(basis % contr_alpha(c_i,i))) - &
                       basis % at_coord(1, basis % at_ind(basis % contr_alpha(c_j,j))),  &
                       basis % at_coord(2, basis % at_ind(basis % contr_alpha(c_i,i))) - &
                       basis % at_coord(2, basis % at_ind(basis % contr_alpha(c_j,j))),  &
                       basis % at_coord(3, basis % at_ind(basis % contr_alpha(c_i,i))) - &
                       basis % at_coord(3, basis % at_ind(basis % contr_alpha(c_j,j))))

           !integral of the product gaussian; from formula (2.12) from Taketa, Huzinaga and O-Hata, J. of the Phys. Soc. of Japan, vol. 21, No. 11, 1966 
           ol=0.0_wp
           f_x=0.0_wp; f_y=0.0_wp; f_z=0.0_wp;
           DO t=0,floor((exp_x_i+exp_x_j)/2.0_wp) !i
              f_x=f_x+f_cf(2*t,exp_x_i,exp_x_j,px-basis%at_coord(1,basis%at_ind(basis%contr_alpha(c_i,i))),&
                          &px-basis%at_coord(1,basis%at_ind(basis%contr_alpha(c_j,j))))*D_FACTORIAL(2*t-1)/((2.0_wp*gam)**t)
           END DO
           DO u=0,floor((exp_y_i+exp_y_j)/2.0_wp) !j
              f_y=f_y+f_cf(2*u,exp_y_i,exp_y_j,py-basis%at_coord(2,basis%at_ind(basis%contr_alpha(c_i,i))),&
                          &py-basis%at_coord(2,basis%at_ind(basis%contr_alpha(c_j,j))))*D_FACTORIAL(2*u-1)/((2.0_wp*gam)**u)
           END DO
           DO v=0,floor((exp_z_i+exp_z_j)/2.0_wp) !k
              f_z=f_z+f_cf(2*v,exp_z_i,exp_z_j,pz-basis%at_coord(3,basis%at_ind(basis%contr_alpha(c_i,i))),&
                          &pz-basis%at_coord(3,basis%at_ind(basis%contr_alpha(c_j,j))))*D_FACTORIAL(2*v-1)/((2.0_wp*gam)**v)
           END DO

           !we don't have to calculate the norm of the primitive gaussians here since it has been already included in the 
           !contraction coefficients during the normalization of primitives!! If you want to use 'ni', 'nj' then you need to remove the normalization of the 
           !primitive gaussians in the subroutine 'normalize'

           ni=1.0_wp!NORM_PRIM_GTO(basis%alpha(basis%contr_alpha(c_i,i)),exp_x_i,exp_y_i,exp_z_i)
           nj=1.0_wp!NORM_PRIM_GTO(basis%alpha(basis%contr_alpha(c_j,j)),exp_x_j,exp_y_j,exp_z_j)

           !assemble all the calculated factors above to the formula (2.12)
           n=(pi/gam)**(3.0_wp/2.0_wp)*exp(-basis%alpha(basis%contr_alpha(c_i,i))*basis%alpha(basis%contr_alpha(c_j,j))/gam*ABsq)
           ol=n*f_x*f_y*f_z*ni*nj

           ol=ol*basis%contractions(c_i,i)*basis%contractions(c_j,j) !multiply the overlap by the corresponding contraction coefficients
           tmp=tmp+ol                                                !increment the final integral

     END DO
  END DO

  OVERLAP_CGTO=tmp

  END FUNCTION OVERLAP_CGTO

  SUBROUTINE NORMALIZE_CGTO(iprint)
  !This normalizes all primitive and contracted basis functions by multiplying the contraction coefficients by calculated normalization factors
  use basis_atoms
  use precisn
  use saved_vars
  IMPLICIT NONE
  !INPUT: the common variable 'basis'
  LOGICAL :: iprint  !do we want a detailed printout of the contraction coefficients?
  !auxiliary
  INTEGER :: i, j, counter, c, exp_x_i,exp_x_j,exp_y_i,exp_y_j,exp_z_i,exp_z_j
  REAL(kind=wp) :: test, tmp, gn
  
  !Normalize all primitives (see the comments in the function overlap_cgto)!
  WRITE (*,'(/,5X,"Primitive cartesian GTOs normalized first.")')
  DO c=1,basis%ncgto          !over all contracted gaussians
     DO i=1,basis%no_contr(c) !over all primitives of the contraction 'c'
           exp_x_i=basis%ang_exp(basis%contr_alpha(c,i),1) !angular exponents of the primitive gaussian
           exp_y_i=basis%ang_exp(basis%contr_alpha(c,i),2)
           exp_z_i=basis%ang_exp(basis%contr_alpha(c,i),3)
           basis % contractions(c,i) = basis % contractions(c,i) &
                                     * NORM_PRIM_GTO(basis % alpha(basis % contr_alpha(c,i)), exp_x_i, exp_y_i, exp_z_i)
     END DO
  END DO
  
  IF (iprint) THEN !detailed printout requested
     WRITE (*, '(/,5X,"Contractions of normalized basis functions:")')
     WRITE (*,'(5X,"contr|Norm of the contr.func.|Normalized contractions")')
  END IF

  DO c=1,basis%ncgto
     tmp=OVERLAP_CGTO(c,c) !calculate the overlap integral
     basis%contractions(c,1:basis%no_contr(c))=basis%contractions(c,1:basis%no_contr(c))/sqrt(tmp) !multiply the contraction coefficients through
     IF (iprint) then
        WRITE (*,'(5X,i5,"|",f20.10,3X,"|",10f15.10)') c,1_wp/sqrt(tmp),(basis%contractions(c,counter),counter=1,basis%no_contr(c))
     end if
  END DO

  WRITE (*,'(/,5X,"Contracted cartesian GTOs normalized.")')

  END SUBROUTINE NORMALIZE_CGTO

  SUBROUTINE NORMALIZE_SABF(iprint)
  !This normalizes all symmetry adapted functions
  use basis_atoms
  use precisn
  use saved_vars
  IMPLICIT NONE
  !INPUT: the common variable 'basis' 
  LOGICAL :: iprint  !do we want a detailed printout?
  !auxiliary
  INTEGER :: i, j, k, counter
  REAL(kind=wp) :: tmp

  IF (iprint) THEN !detailed printout requested
     WRITE (*, '(/,5X,"Coefficients for normalized symmetry adapted basis functions:")')
     WRITE (*,'(5X,"sa_in|Norm of the sym.ad.fun.|Normalized coefficients")')
  END IF

  DO i=1,basis%nsabf              !over all symmetry-adapted functions
     tmp=0.0_wp
     DO j=1,basis%no_sabf(i)      !over all contractions in the symmetry adapted function 'i'
        DO k=1,basis%no_sabf(i)   !over all contractions in the symmetry adapted function 'i'
           tmp=tmp+basis%sa_contr(i,j)*basis%sa_contr(i,k)*OVERLAP_CGTO(basis%contr_ind(i,j),basis%contr_ind(i,k))
        END DO
     END DO
     basis%sa_contr(i,1:basis%no_sabf(i))=basis%sa_contr(i,1:basis%no_sabf(i))/sqrt(tmp)
     IF (iprint) then
        WRITE (*,'(5X,i5,"|",f20.10,3X,"|",10f15.10)') i,1_wp/sqrt(tmp),(basis%sa_contr(i,counter),counter=1,basis%no_sabf(i))
     end if
  END DO

  WRITE (*,'(/,5X,"Symmetry adapted functions normalized.")')

  END SUBROUTINE NORMALIZE_SABF

  FUNCTION OVERLAP_SABF(c_i,c_j)
  !calculates analytically the overlap integral between the symmetry adapted functions c_i and c_j
  use basis_atoms
  use precisn
  !INPUT: the common variable 'basis'
  IMPLICIT NONE
  INTEGER :: i, j, k, c_i, c_j
  REAL(kind=wp) :: OVERLAP_SABF, tmp

  !no need to optimize the sum below - it is quick enough as it is
  tmp=0.0_wp
  DO j=1,basis%no_sabf(c_i)    !over all contractions in the symmetry adapted function 'c_i'
     DO k=1,basis%no_sabf(c_j) !over all contractions in the symmetry adapted function 'c_j'
        tmp=tmp+basis%sa_contr(c_i,j)*basis%sa_contr(c_j,k)*OVERLAP_CGTO(basis%contr_ind(c_i,j),basis%contr_ind(c_j,k))
     END DO
  END DO
  OVERLAP_SABF=tmp

  END FUNCTION OVERLAP_SABF

  FUNCTION NORM_PRIM_GTO(alp,i,j,k)
  use el_fns, only: FACTORIAL
  use precisn
  use constants, only: pi
  !Returns the normalization factor for a primitive GTO: x^i*y^j*z^k*exp[-alp*(x^2+y^2+z^2)]
  IMPLICIT NONE
  REAL(kind=wp) :: alp, NORM_PRIM_GTO
  INTEGER :: i,j,k

  NORM_PRIM_GTO = (2.0_wp*alp/pi)**(0.75_wp) * sqrt((8_wp*alp)**(i+j+k) * FACTORIAL(i) * FACTORIAL(j) * FACTORIAL(k) &
                / (FACTORIAL(2*i) * FACTORIAL(2*j) * FACTORIAL(2*k)))
  
  END FUNCTION NORM_PRIM_GTO

  FUNCTION PRIM_GAUSS(ind,x,y,z)
  use el_fns, only: R_SQ
  use basis_atoms
  use prefactors
  use precisn
  !returns the value of the primitive gaussian 'ind' at x,y,z; The coordinates x,y,z are assumed to be measured from the molecular origin.
  IMPLICIT NONE
  REAL(kind=wp) :: PRIM_GAUSS
  !INPUT: the common variable 'basis';'ind','x','y','z' explained above
  INTEGER :: ind             !index of the primitive gaussian
  REAL(kind=wp) :: x, y, z   !coordinates (in a.u.) of the point where I want to evaluate the primitive gaussian
  REAL(kind=wp) :: x_tmp, y_tmp, z_tmp

  x_tmp=x-basis%at_coord(1,basis%at_ind(ind))
  y_tmp=y-basis%at_coord(2,basis%at_ind(ind))
  z_tmp=z-basis%at_coord(3,basis%at_ind(ind))

  IF (basis % ang_exp(ind,1) == 0 .and. &
      basis % ang_exp(ind,1) == basis % ang_exp(ind,2) .and. &
      basis % ang_exp(ind,1) == basis % ang_exp(ind,3)) THEN !avoid slow evaluation of the power functions for s-type gaussians

    PRIM_GAUSS = exp(-basis%alpha(ind)*R_SQ(x_tmp,y_tmp,z_tmp))

  ELSE

    PRIM_GAUSS = (x_tmp)**basis%ang_exp(ind,1)*(y_tmp)**basis%ang_exp(ind,2)*(z_tmp)**basis%ang_exp(ind,3)*&
                &exp(-basis%alpha(ind)*R_SQ(x_tmp,y_tmp,z_tmp))
  ENDIF

  END FUNCTION PRIM_GAUSS

  FUNCTION CONTR_GAUSS(ic,x,y,z)
  use basis_atoms
  use precisn
  !returns the value of a CONTRACTED gaussian 'ic' at x,y,z; The coordinates x,y,z are assumed to be measured from the molecular origin.
  !INPUT: the common variable 'basis';'ic','x','y','z' explained above
  IMPLICIT NONE
  INTEGER :: ic, i, j
  REAL(kind=wp) :: x, y, z, CONTR_GAUSS

  CONTR_GAUSS=0.0_wp
  DO i=1,basis%no_contr(ic) !loop over primitives
        CONTR_GAUSS=CONTR_GAUSS + basis%contractions(ic,i)*PRIM_GAUSS(basis%contr_alpha(ic,i),x,y,z)
  END DO
  END FUNCTION CONTR_GAUSS

  FUNCTION SA_BF(isa,x,y,z)
  use basis_atoms
  use precisn
  !returns the value of a symmetry adapted function 'is' at x,y,z; The coordinates x,y,z are assumed to be measured from the molecular origin.
  !INPUT: the common variable 'basis';'is','x','y','z' explained above
  IMPLICIT NONE
  INTEGER :: isa, i, j
  REAL(kind=wp) :: x, y, z, SA_BF

  SA_BF=0.0_wp
  DO i=1,basis%no_sabf(isa) !loop over contractions
        SA_BF=SA_BF + basis%sa_contr(isa,i)*CONTR_GAUSS(basis%contr_ind(isa,i),x,y,z)
  END DO
  END FUNCTION SA_BF

END MODULE basis_fns

MODULE basis_fns_num
!This module contains routines neccessary for numerical integration of basis functions and orbitals

CONTAINS

  FUNCTION INT_FUNC_SIMPLE(F, A, B)
  use precisn
  !integrates an arbitrary real function 'F' over an interval [A,B]
  !'n_int' is an index 1->3 saying which integration is being performed (over 'r' or 'theta' or 'phi'). The rest are auxiliary variables required by DQAGS
  IMPLICIT NONE
  INTEGER, PARAMETER :: LIMIT=5000
  INTEGER, PARAMETER :: LENW=LIMIT*4
  REAL(kind=wp), external :: F
  REAL(kind=wp), PARAMETER :: EPSABS=1d-09 !controls precision of each integration
  REAL(kind=wp), PARAMETER :: EPSREL=1d-09
  INTEGER :: NEVAL, IER, LAST, IWORK(1:LIMIT)
  REAL(kind=wp) :: A, B, INT_FUNC_SIMPLE, ABSERR, WORK(1:LENW)

  !adaptive integrator from the Slatec library  
  CALL DQAGS (F, A, B, EPSABS, EPSREL, INT_FUNC_SIMPLE, ABSERR, NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)

  END FUNCTION INT_FUNC_SIMPLE

  FUNCTION INT_FUNC(F, A, B, n_int)
  use precisn
  use saved_vars
  !integrates an arbitrary function 'F' over the interval (A,B); if an imprecision of integration is detected the interval (A,B) will be subdivided up to
  !10 times until convergence is reached. This is here in order to mitigate possible imprecisions arising from the extrapolation used by the DQAGS routine.
  CHARACTER(5), PARAMETER :: str(1:3) = (/'r    ','theta','phi  '/) !string used in printout of the error message 4
  INTEGER, PARAMETER :: num_int=3 !number of integrations performed at once - each requires a stand alone (i)work arrays
  INTEGER, PARAMETER :: LIMIT=5000
  INTEGER, PARAMETER :: SUBDIV=100 !maximum number of subdivisions of the starting interval (A,B)
  INTEGER, PARAMETER :: LENW=LIMIT*4
  REAL(kind=wp), external :: F
  REAL(kind=wp), PARAMETER :: EPSABS(1:num_int)=(/1d-09,1d-09,1d-09/)
  REAL(kind=wp), PARAMETER :: EPSREL(1:num_int)=(/1d-09,1d-09,1d-09/)
  INTEGER :: NEVAL, LAST, IWORK(1:LIMIT,1:num_int), n_int, test_a(1:SUBDIV,1:num_int), test_b(1:SUBDIV,1:num_int)
  INTEGER :: now_chck(1:num_int), i, no_int(1:num_int), IER(1:num_int), gone_th(1:SUBDIV,1:num_int)
  REAL(kind=wp) :: A, B, INT_FUNC, ABSERR, WORK(1:LENW,1:num_int), INT_TMP(1:num_int), lim_A, lim_B, middle, INT_ALL(1:num_int)
  LOGICAL :: err_det(1:num_int)
  
  IER(n_int) = 0
  err_det(n_int)=.false. !did we detect an insufficient precision during this run?
  INT_FUNC=0d0
  no_int(n_int)=1 !number of intervals which were not integrated to the desired precision. no_int(n_int)=1 is the necessary initial condition.
  lim_A = A
  lim_B = B
  gone_th(:,n_int)=0 !did I already go through this node?
  test_a(:,n_int)=0d0; test_b(:,n_int)=0d0 !array of intervals to test
  now_chck(n_int)=1   !index of the node (interval) which will be tested for converged integration
  INT_TMP(n_int)=0d0
  INT_ALL(n_int)=0d0
  
  !CALL DQAGS (F, lim_A, lim_B, EPSABS(n_int),EPSREL(n_int), INT_FUNC, ABSERR, NEVAL, IER(n_int), LIMIT, LENW, LAST, IWORK(:,n_int), WORK(:,n_int)) 
  !above: the original one-go integration only
  
  DO !the loop of the subdivision algorithm

    !try to integrate the function over the test interval lim_A, lim_B
    CALL DQAGS (F, lim_A, lim_B, EPSABS(n_int), EPSREL(n_int), INT_TMP(n_int), ABSERR, NEVAL, IER(n_int), &
                LIMIT, LENW, LAST, IWORK(:,n_int), WORK(:,n_int))
    IF (IER(n_int) .eq. 4) THEN  !the integration did not converge to the desired limit? -> divide the integration into two pieces and continue

      err_det(n_int)=.true. !an imprecision in integration was detected. Used in the ELSE part of this IF.

      IF (now_chck(n_int) .eq. 1) INT_ALL(n_int)=INT_TMP(n_int) !store integral over the original interval. Used if the subdivision does not converge.
!      WRITE (*,'(a,a5,a,f10.5,";",f10.5,"; INT = ",e22.15)'),&
!                &'#ERROR NUMBER = 4: Insufficient precision detected in integration over ',str(n_int),' along the interval:',lim_A, lim_B,INT_TMP(n_int)
!      WRITE(*,'("#Trying to subdivide the interval and achieve the required precision...")')
      gone_th(now_chck(n_int),n_int)=1  !I just got through the node with the index now_chck(n_int) and limits lim_A, lim_B
      no_int(n_int)=no_int(n_int)+2     !total number of intervals to check: +2 = the original interval was checked, but two new subintervals emerge.

      IF (no_int(n_int) > SUBDIV) THEN  !the maximum number of allowed subdivisions was reached
         WRITE(*, '("#The number of subdivisions exceeded SUBDIV! The desired precision of integration over ",a5 &
              &," cannot be reached by subdivision.")') str(n_int)
         INT_FUNC=INT_ALL(n_int)  !return as the result the integral over the original interval.
         EXIT
      ENDIF

      !calculate limits for the two subintervals and add these intervals on the list of the intervals which have to be checked
      middle=min(lim_A,lim_B)+abs(lim_B-lim_A)/2d0
      !first subinterval
      test_a(no_int(n_int)-1,n_int)=min(lim_A,lim_B)
      test_b(no_int(n_int)-1,n_int)=middle
      !second subinterval
      test_a(no_int(n_int),n_int)=middle
      test_b(no_int(n_int),n_int)=max(lim_A,lim_B)
      !I've not gone through these two nodes yet
      gone_th(no_int(n_int),n_int)=0
      gone_th(no_int(n_int)-1,n_int)=0
      !Choose the next interval to be checked as the first (the choice is arbitrary) of the two just calculated
      lim_A=test_a(no_int(n_int)-1,n_int)
      lim_B=test_b(no_int(n_int)-1,n_int)
      now_chck(n_int)=no_int(n_int)-1     !this is the index of the node which will be tested now
  !    WRITE (*,'(a,2f20.15)'),'#Integrating along the interval:',lim_A, lim_B

    ELSE !integration over the current interval was OK

       gone_th(now_chck(n_int),n_int)=1  !I just went through the now_chck(n_int) node
       INT_FUNC=INT_FUNC+INT_TMP(n_int)  !add the integral over this subinterval to the toal integral. I am adding only those integrals which were OK.

       !search for a node which has not been tried yet:
       now_chck(n_int)=0 !this marks that no intervals are left for integration
       IF (err_det(n_int)) THEN !only search for another intervals if we've detected an imprecision in the first integration.
                                            !We don't want to waste time if the first integral is OK.
          DO i=1,no_int(n_int)
             IF (gone_th(i,n_int) .eq. 0) THEN !if such a node (i-th) has been found then check it
                now_chck(n_int)=i     !this interval has not been tested yet
                lim_A=test_a(i,n_int)
                lim_B=test_b(i,n_int)
      !          WRITE (*,'(a,2f20.15)'),'#Integrating along the interval:',lim_A, lim_B
                EXIT
             ENDIF
          ENDDO

       ENDIF

       !were all the nodes integrated successfuly? -> EXIT
       IF (now_chck(n_int) .eq. 0) THEN
          !only write this message if we tried to subdivide the integration limits
          IF (err_det(n_int)) &
             &WRITE(*,'("#Subdivision of the parent interval over ",a5," was needed. r, theta, phi: ",f10.5,";",f10.5,";",f10.5,&
             &"; Number of integrations: ",i5)') str(n_int), r,theta, phi, no_int(n_int)
          EXIT
       ENDIF

    ENDIF

  ENDDO
  
  END FUNCTION INT_FUNC

  FUNCTION CONTR_G_r(rad)
  !Returns the value of the overlap function between the contracted GTOs 'contr_i' and 'contr_j' integrated over theta, phi; evaluated at r=rad
  use precisn
  use basis_atoms
  use saved_vars
  !INPUT: the common variables 'contr_i', 'contr_j', radius 'rad'
  IMPLICIT NONE
  REAL(kind=wp) :: rad !radius at which I want to evaluate the primary gaussian
  REAL(kind=wp) :: CONTR_G_r
  
  !evaluate the gaussian*Jac(r) 'contr_i' at saved values of 'theta' and 'phi'
  r=rad
  CONTR_G_r = (r*r)*INT_FUNC(CONTR_G_theta,0.0_wp,pi,2)

  END FUNCTION

  FUNCTION CONTR_G_theta(th)
  !Returns the value of the overlap function between the contracted GTOs 'contr_i' and 'contr_j' integrated over phi; evaluated at r=rad, theta=th
  use precisn
  use basis_atoms
  use saved_vars
  !INPUT: the common variables 'contr_i', 'contr_j', angle theta 'th'
  IMPLICIT NONE
  REAL(kind=wp) :: th !theta at which I want to evaluate the primary gaussian integrated over phi
  REAL(kind=wp) :: CONTR_G_theta
   
  !evaluate the gaussian*Jac(theta) 'contr_i' at saved values of 'r' and 'phi'
  theta=th
  prev_s_th = sin(theta) !saves time in evaluating the corresponding X, Y, Z coordinates for the same value of theta in CONTR_G_phi. Also see GET_X, Y, Z
  prev_c_th = cos(theta) 
  CONTR_G_theta = prev_s_th*INT_FUNC(CONTR_G_phi,0.0_wp,2.0_wp*pi,3)

  END FUNCTION

  FUNCTION CONTR_G_phi(ph)
  !Returns the value of the overlap function between the contracted GTOs 'contr_i' and 'contr_j' evaluated at r=rad, theta=th, phi=ph
  use precisn
  use basis_atoms
  use saved_vars
  use el_fns, only: GET_X,GET_Y,GET_Z
  use basis_fns, only: CONTR_GAUSS
  !INPUT: the common variables 'contr_i', 'contr_j', angle phi 'ph'
  IMPLICIT NONE
  REAL(kind=wp) :: ph !phi at which I want to evaluate the primary gaussian evaluated at 'rad', 'th'
  REAL(kind=wp) :: CONTR_G_phi,tmp1,tmp2
  REAL(kind=wp) :: x_tmp, y_tmp, z_tmp

  !evaluate the overlap function between contractions 'contr_i' and 'contr_j' at saved values of 'r' and 'theta'
  phi=ph
  x_tmp=GET_X(r,ph)
  y_tmp=GET_Y(r,ph)
  z_tmp=GET_Z(r,ph)

  tmp1=CONTR_GAUSS(contr_i,x_tmp,y_tmp,z_tmp)
  tmp2=CONTR_GAUSS(contr_j,x_tmp,y_tmp,z_tmp)
  CONTR_G_phi = tmp1*tmp2
  
  END FUNCTION

  FUNCTION SABF_r(rad)
  !Returns the value of the overlap function between the symmetry adapted functions 'sabf_i' and 'sabf_j' integrated over theta, phi; evaluated at r=rad
  use precisn
  use basis_atoms
  use saved_vars
  IMPLICIT NONE
  REAL(kind=wp) :: rad !radius at which I want to evaluate the symmetry-adapted function
  REAL(kind=wp) :: SABF_r

  r=rad
  SABF_r = (r*r)*INT_FUNC(SABF_theta,0.0_wp,pi,2)

  END FUNCTION

  FUNCTION SABF_theta(th)
  !Returns the value of the overlap function between the functions 'sabf_i' and 'sabf_j' integrated over phi; evaluated at r=rad, theta=th
  use precisn
  use basis_atoms
  use saved_vars
  IMPLICIT NONE
  REAL(kind=wp) :: th !theta at which I want to evaluate the symmetry-adapted function
  REAL(kind=wp) :: SABF_theta

  theta=th
  prev_s_th = sin(theta) !saves time in evaluating GET_X, Y, Z for the same value of theta in CONTR_G_phi
  prev_c_th = cos(theta)
  SABF_theta = prev_s_th*INT_FUNC(SABF_phi,0.0_wp,2.0_wp*pi,3)

  END FUNCTION

  FUNCTION SABF_phi(ph)
  !Returns the value of the overlap function between the functions 'sabf_i' and 'sabf_j' evaluated at r=rad, theta=th, phi=ph
  use precisn
  use basis_atoms
  use saved_vars
  use el_fns, only: GET_X,GET_Y,GET_Z
  use basis_fns, only: SA_BF
  IMPLICIT NONE
  REAL(kind=wp) :: ph !phi at which I want to evaluate the symmetry-adapted function
  REAL(kind=wp) :: SABF_phi,tmp1,tmp2
  REAL(kind=wp) :: x_tmp, y_tmp, z_tmp

  !evaluate the overlap function between the symmetry ad. fns. 'sabf_i' and 'sabf_j' at saved values of 'r' and 'theta' and 'phi'
  phi=ph
  x_tmp=GET_X(r,ph)
  y_tmp=GET_Y(r,ph)
  z_tmp=GET_Z(r,ph)

  tmp1=SA_BF(sabf_i,x_tmp,y_tmp,z_tmp)
  tmp2=SA_BF(sabf_j,x_tmp,y_tmp,z_tmp)
  SABF_phi = tmp1*tmp2

  END FUNCTION

  FUNCTION ORB_DEN_r(rad)
  !Returns charge density of the orbital/density 'orb_i' at r=rad
  use precisn
  use basis_atoms
  use saved_vars  !this is where value of the 'orb_i' comes from
  IMPLICIT NONE
  REAL(kind=wp) :: ORB_DEN_r
  !INPUT:
  REAL(kind=wp) :: rad       !radius at which we want the density
  !value of the 'orb_i' variable which is set externally

  r=rad
  ORB_DEN_r = (r*r)*INT_FUNC(ORB_DEN_theta,0.0_wp,pi,2)
  !radius in Bohrs and Angstroms,charge density,charge density normalized to 1
  WRITE (*,'(f20.10, f20.10, f20.10, f20.10)') r, r/to_au, ORB_DEN_r,ORB_DEN_r/orbital(orb_i)%tot_ch

  END FUNCTION ORB_DEN_r

  FUNCTION ORB_DEN_theta(th)
  !Returns charge density of the orbital 'orb_i' at r=rad, theta=th
  use precisn
  use basis_atoms
  use saved_vars !this is where value of the 'orb_i' comes from
  IMPLICIT NONE
  REAL(kind=wp) :: ORB_DEN_theta
  !INPUT:
  REAL(kind=wp) :: th !theta at which we want the density
  !value of the 'orb_i' variable which is set externally

  theta=th
  prev_s_th = sin(theta) !saves time in evaluating GET_X, Y, Z for the same value of theta in CONTR_G_phi
  prev_c_th = cos(theta)
  ORB_DEN_theta = prev_s_th*INT_FUNC(ORB_DEN_phi,0.0_wp,2.0_wp*pi,3)

  END FUNCTION ORB_DEN_theta

  FUNCTION ORB_DEN_phi(ph)
  !Returns charge density of the orbital 'orb_i' at r=rad, theta=th, phi=ph
  use precisn
  use basis_atoms
  use saved_vars
  use el_fns, only: GET_X,GET_Y,GET_Z
  use basis_fns, only: CONTR_GAUSS, SA_BF
  IMPLICIT NONE
  REAL(kind=wp) :: ORB_DEN_phi
  !INPUT:
  REAL(kind=wp) :: ph          !phi at which we want to evaluate the primary gaussian
  !value of the 'orb_i' variable which is set externally and the common variable 'orbital'
  !auxiliary:
  INTEGER :: i, j
  REAL(kind=wp) :: x_tmp, y_tmp, z_tmp, bf_i_tmp

  x_tmp=GET_X(r,ph)
  y_tmp=GET_Y(r,ph)
  z_tmp=GET_Z(r,ph)

  !evaluate the density of orbital 'orb_i' at 'r' and 'theta' and 'phi'
  phi=ph
  ORB_DEN_phi = 0.0_wp
  IF (iscgto) THEN !the basis functions are contracted gaussians
     DO i=1,basis%ncgto
        bf_i_tmp=CONTR_GAUSS(i,x_tmp,y_tmp,z_tmp)
        DO j=1,i
           !it actually happens quite often that dm(i,j)=0.0_wp; inclusion of this conditioning pays off very much (despite the break up of efficient compilation of the double do loop): speed up of about 2x was achieved for the 
           !1st orbital of water. For molecules with high symmetry inclusion of the conditioning should result in a better scaling when comparing the speed of calculations done in a small and then in a bigger basis set.
           IF (orbital(orb_i)%dm(i,j) .ne. 0.0_wp) ORB_DEN_phi = & 
             &ORB_DEN_phi + orbital(orb_i)%dm(i,j)*bf_i_tmp*CONTR_GAUSS(j,x_tmp,y_tmp,z_tmp)
        END DO
     END DO
  ELSE                         !the basis functions are symmetry adapted functions
     DO i=1,basis%nsabf
        bf_i_tmp=SA_BF(i,x_tmp,y_tmp,z_tmp)
        DO j=1,i
           IF (orbital(orb_i)%dm(i,j) .ne. 0.0_wp) ORB_DEN_phi = & !it actually happens quite often that dm(i,j)=0.0_wp; see the comment above.
             &ORB_DEN_phi + orbital(orb_i)%dm(i,j)*bf_i_tmp*SA_BF(j,x_tmp,y_tmp,z_tmp) 
        END DO
     END DO
  END IF
  
  END FUNCTION ORB_DEN_phi

  SUBROUTINE CHECK_GTO_NORMS
  !Numerically integrates the overlaps of contracted basis functions and compares them with analytical expressions
  use basis_atoms
  use saved_vars
  use constants
  use basis_fns, only: OVERLAP_CGTO
  !INPUT: the common variable 'basis'
  IMPLICIT NONE
  INTEGER :: i, j
  REAL(kind=wp) :: test,ol_anl
 
  WRITE (*,'(/,5X,"Checking norms of contracted GTOs")')
  WRITE (*,'(5X,"c_i  |c_j |Numerical integral  |Analytic           |Difference ")')
  DO i=1,basis%ncgto
     DO j=1,basis%ncgto
        contr_i=i     !index of the first CGTO
        contr_j=j     !index of the second CGTO
        test=INT_FUNC(CONTR_G_r,0.0_wp,r_inf,1)   !overlap <i|j>
        !basis%overlaps_cgto(i,j)=test
        ol_anl=OVERLAP_CGTO(i,j)
        WRITE (*,'(5X,2i5,"|",3e20.11)') i,j,test,ol_anl,test-ol_anl
     END DO
  END DO

  END SUBROUTINE CHECK_GTO_NORMS

  SUBROUTINE CHECK_SABF_NORMS
  !Numerically integrates the overlaps of contracted basis functions and compares them with analytical expressions
  use basis_atoms
  use saved_vars
  use constants
  use basis_fns, only: OVERLAP_CGTO
  !INPUT: the common variable 'basis'
  IMPLICIT NONE
  INTEGER :: i, j, p, q
  REAL(kind=wp) :: test,ol_anl
 
  WRITE (*,'(/,5X,"Checking norms of symmetry adapted functions")')
  WRITE (*,'(5X,"s_i  |s_j |Numerical integral  |Analytic           |Difference ")')
  DO i=1,basis%nsabf
     DO j=1,basis%nsabf
        sabf_i=i         !index of the first symmetry adapted function
        sabf_j=j         !index of the second symmetry adapted function
        test=INT_FUNC(SABF_r,0.0_wp,r_inf,1) !numerical overlap <i|j>
        ol_anl=0.0_wp
        DO p=1,basis%no_sabf(i)
           DO q=1,basis%no_sabf(j)
              ol_anl=ol_anl+basis%sa_contr(i,p)*basis%sa_contr(j,q)*OVERLAP_CGTO(basis%contr_ind(i,p),basis%contr_ind(j,q)) !analytic overlap <i|j>
           END DO
        END DO
        WRITE (*,'(5X,2i5,"|",3e20.11)') i,j,test,ol_anl,test-ol_anl
     END DO
  END DO
  
  ENDSUBROUTINE CHECK_SABF_NORMS

  SUBROUTINE CALC_DENSITIES(req_dm,rstart,rfinish)
  !This is the main routine which calculates the radial densities for all required orbitals/density matrices
  use saved_vars
  use constants
  use precisn
  IMPLICIT NONE
  !INPUT:
  INTEGER :: req_dm                   !number of required orbitals/density matrices
  REAL(kind=wp) :: rstart, rfinish    !interval for radial integration
  !auxiliary
  INTEGER :: i
  REAL(kind=wp) :: integral

  IF (rstart .eq. rfinish) THEN
     WRITE(*,'(/,5X,a,f12.5)') "Radial charge density will be calculated only for one point [Bohr] : r = ",rstart
  ELSE
     WRITE(*,'(/,5X,a,"[",f12.5,",",f12.5,"]")') "Radial densities will be calculated on interval [Bohr] : ",rstart,rfinish
  ENDIF

  DO orb_i=1,req_dm   !for all required orbitals/density matrices
     WRITE (*,'(/,5X,a,i5)') 'Calculating densities for density matrix: ',orb_i
     WRITE (*,'(5X,a)') 'Radius [Bohr]      |Radius [Angstrom]  |Charge density     |Charge density/total charge'
     if (rstart .eq. rfinish) then !density only at one point is required
        integral = ORB_DEN_r(rstart)
     else !we require integral over a certain radial range
        integral=INT_FUNC(ORB_DEN_r,rstart,rfinish,1)
        WRITE (*,'(5X,"Integrated charge density for density matrix ",i5," = ",f20.15)') orb_i, integral
     endif
  END DO

  WRITE (*,'(/,5X,a)') 'Calculation of charge densities finished.'
 
  END SUBROUTINE CALC_DENSITIES

END MODULE basis_fns_num

MODULE molden

CONTAINS

  SUBROUTINE READ_MOLDEN(mldu,mld_file,dmfromf)
  !reads all basis set information from the Molden file (mldu) to the common variable 'basis'
  use basis_atoms
  IMPLICIT NONE
  !INPUT:
  INTEGER :: mldu                      !unit number for the Molden file
  LOGICAL :: dmfromf                   !do we want to have density matrices read from a file?
  CHARACTER(len=line_len) :: mld_file  !path to the Molden file
  !auxiliary:
  INTEGER :: io

  !make sure that the input file contains information on basis set and molecular orbitals
  CALL CHECK_FILE(mldu)
  
  !how many atoms do we have?
  basis%no_at=NOATOMS(mldu)
  IF (basis%no_at > max_at) THEN
     STOP "Number of atoms exceeds the current limit set in the code. Change the parameter 'max_at' and recompile."
  END IF
  ALLOCATE(basis%at_coord(1:3,1:basis%no_at),basis%at_num(1:basis%no_at),basis%at_sym(1:basis%no_at))
  
  !get information on coordinates of atoms
  CALL GETATOMS(mldu)
  
  !get bas set dimensions and other parameters needed for allocation of arrays
  CALL GETBASISINFO(mldu)
  
  !now allocate arrays for the cartesian basis set
  ALLOCATE(basis % alpha(1 : basis % nprim), &
           basis % contractions(1 : basis % ncgto, 1 : basis % maxall), &
           basis % at_ind(1 : basis % nprim), &
           basis % prefact(1 : basis % nprim), &
           basis % contr_alpha(1 : basis % ncgto, 1 : basis % maxall), &
           basis % ang_exp(1 : basis % nprim, 1 : 3), &
           basis % overlaps_cgto(1 : basis % ncgto, 1 : basis % ncgto), &
           basis % no_contr(1 : basis % ncgto))
  
  basis % alpha(:) = 0.0_wp
  basis % contractions(:,:) = 0.0_wp
  basis % at_ind(:) = 0
  basis % prefact(:) = ''
  basis % contr_alpha(:,:) = 0
  basis % ang_exp(:,:) = 0
  basis % overlaps_cgto(:,:) = 0.0_wp
  basis % no_contr(:) = 0
  
  !alloacte arrays for the spherical gaussian basis set
  ALLOCATE(basis % sph_a(1 : basis % sph_nprim), &
           basis % sph_contr(1 : basis % sph_ncgto, 1 : basis % sph_maxall), &
           basis % sph_at_ind(1 : basis % sph_nprim), &
           basis % sph_contr_a(1 : basis % sph_ncgto, 1 : basis % sph_maxall), &
           basis % sph_olap_cgto(1 : basis % sph_ncgto, 1 : basis % sph_ncgto), &
           basis % sph_no_contr(1 : basis % sph_ncgto), &
           basis % l_ang(1 : basis % sph_nprim), &
           basis % m_ang(1 : basis % sph_nprim))

  !now we can read all basis set information to the variable 'basis'
  CALL GETBASIS(mldu)

  END SUBROUTINE READ_MOLDEN

  SUBROUTINE CHECK_FILE(IUNIT)
  !Checks that the file contains the information in the Molden format
  IMPLICIT NONE
  !INPUT:
  INTEGER :: IUNIT !unit number of the Molden file
  !auxiliary:
  INTEGER :: ifail
  
  WRITE (*,'(5X,a)') 'Attempt to find orbital coefficients...'
  
  CALL SEARCH(IUNIT,'[MO]',len('[MO]'),ifail,.true.)
  IF (ifail .ne. 0) THEN
     STOP "The Molden file does not contain orbital coefficients. Program terminated."
  END IF
  
  WRITE (*,'(a)') 'OK'
  
  WRITE (*,'(5X,a)') 'Attempt to find basis set information...'
  
  CALL SEARCH(IUNIT,'[GTO]',len('[GTO]'),ifail,.true.)
  IF (ifail .ne. 0) THEN
     STOP "The Molden file does not contain basis set information. Program terminated."
  END IF
  
  WRITE (*,'(a)') 'OK'

  END SUBROUTINE CHECK_FILE

  SUBROUTINE CALC_DM(molunit,req_dm,ilist,dm_list,norb,dmunit,iprint,header_base)
  !This routine reads the orbitals from the Molden file, constructs density matrices and calculates total charge in all space, i.e. fills in the common
  !variable 'orbital(:)'.
  use constants
  use precisn
  use basis_atoms
  use basis_fns, only: OVERLAP_CGTO
  IMPLICIT NONE
  !INPUT:
  INTEGER :: molunit            !unit number of the Molden file
  INTEGER :: dmunit             !unit number for the file where the calculated density matrices will be written
  INTEGER :: req_dm             !number of orbitals we want to read in
  INTEGER :: norb               !total number of orbitals on the Molden file
  INTEGER :: dm_list(1:max_orb) !list of selected orbitals we want to process
  LOGICAL :: ilist,iprint       !do we want to use the 'dm_list' or read all orbitals instead?;detail printou requested?
  CHARACTER(len=line_len) :: header_base !header for the saved density matrices
  !auxiliary:
  INTEGER :: ifail, i, j, k, m, n, io
  REAL(kind=wp) :: tmp, tot_ch
  CHARACTER(len=line_len) :: line
  
  CALL SEARCH(molunit,'[MO]',len('[MO]'),ifail,.true.)  !header for orbitals
  
  WRITE (*,'(/,5X,a)') 'Searching for available orbitals...'
  
  norb=0
  DO
    CALL SEARCH(molunit,' Ene=',len(' Ene='),ifail,.false.) !scan the file for orbitals
    IF (ifail .ne. 0) THEN
       EXIT
    ELSE
       norb=norb+1
    END IF
  END DO
  
  WRITE (*,'(5X,a,i5)')'Number of orbitals found: ',norb
  
  IF (req_dm < 0) THEN !all orbitals requested
     req_dm=norb
  END IF

  IF (req_dm > norb) THEN
     WRITE(*,'(5X,a)') 'Required number of orbitals higher than available! Reseting to the number of available orbitals.'
     req_dm=norb
  END IF
  
  WRITE (*,'(5X,a,i5)')'Required number of orbitals: ',req_dm

  IF (ilist) THEN !we want to read in only some orbitals
     WRITE (*,'(5X,a,100i5)') "Requested sequence numbers of density matrices: ",(dm_list(i),i=1,req_dm)
     DO i=1,size(dm_list)
        IF (dm_list(i) > norb) then
           STOP "A requested sequence number in 'dm_list' higher than available number of density matrices. Program terminated."
        end if
     END DO
  END IF

  ALLOCATE(orbital(1:req_dm))
  DO i=1,req_dm
     ALLOCATE(orbital(i)%coeff(1:basis%ncgto))
     ALLOCATE(orbital(i)%dm(1:basis%ncgto,1:basis%ncgto))
  END DO
  
  CALL SEARCH(molunit,'[MO]',len('[MO]'),ifail,.true.) !get on the start of molecular orbitals information

  !formats of orbital information as written by Molpro
  !IMPLEMENT!!
  !write (ifile,'(a5,3x,i4,".",i1)') 'Sym=',i,isk
  !write (ifile,'(a5,3x,F10.4)') 'Ene=',eig(i+nts(isk))
  !write(ifile,'(a12)')'Spin= Alpha'
  !write(ifile,'(1x,a6,3x,f8.6)')'Occup=',anocc(i+nts(isk))
  
  WRITE (*,'(5X,a)') 'Reading the required orbitals...'
  
  DO i=1,req_dm

     IF (ilist) THEN !if only some orbitals are required then read the next one from the array 'dm_list'
        REWIND(molunit)
        DO j=1,dm_list(i) !skip the 'dm_list(i)-1' orbitals
           CALL SEARCH(molunit,' Occup=',len(' Occup='),ifail,.false.) !This is an assumption that the orbital coefficients follow this line!
        END DO     
     ELSE                 !we want the orbitals one by one -> read the next orbital
        CALL SEARCH(molunit,' Occup=',len(' Occup='),ifail,.false.) !This is an assumption that the orbital coefficients follow this line!
     END IF

     DO j=1,basis%ncgto   !read the orbital coefficients
        READ(molunit,'(a)') line  !number, coefficient
        READ(line,*) k, orbital(i)%coeff(k) 
     END DO

  END DO

  WRITE (*,'(a)') 'OK'
  
  OPEN(unit=dmunit,form='formatted',iostat=io) !output of the density matrices
  IF (io .ne. 0) THEN
     STOP "Problem opening the DMUINT for density matrices. Program terminated."
  END IF
  
  WRITE (*,'(/,5X,a)') 'Calculating density matrices for the required orbitals...'
  WRITE (*,'(5X,a)') 'Assuming double occupancy of all orbitals.'
  IF (iprint) WRITE (*,'(/,5X,a)') 'xyz c_i|xyz c_j|Overlap integral'  !detailed printout requested
  DO i=1,req_dm
     orbital(i)%tot_ch=0.0_wp  !total charge in the orbital
     DO m=1,basis%ncgto
        DO n=1,m
           IF (m .eq. n) THEN !diagonal elements remain the same
              orbital(i)%dm(m,n)=2.0_wp*orbital(i)%coeff(m)*orbital(i)%coeff(n)
           ELSE               !extra factor for the off-diagonal elements
              orbital(i)%dm(m,n)=2.0_wp*2.0_wp*orbital(i)%coeff(m)*orbital(i)%coeff(n)
           END IF
           !overlap <m|n> needed for calculation of the total charge
           tmp=OVERLAP_CGTO(m,n)
           basis%overlaps_cgto(m,n)=tmp
           IF (i .eq. 1 .and. iprint) THEN !print out the normalized basis set if required
              WRITE (*,'(5X,3i1,i4,"|",3i1,i4,"|",2f20.15)') &
                basis % ang_exp(basis % contr_alpha(m,1), 1), &
                basis % ang_exp(basis % contr_alpha(m,1), 2), &
                basis % ang_exp(basis % contr_alpha(m,1), 3), m, &
                basis % ang_exp(basis % contr_alpha(n,1), 1), &
                basis % ang_exp(basis % contr_alpha(n,1), 2), &
                basis % ang_exp(basis % contr_alpha(n,1), 3), n, tmp
           END IF
           orbital(i)%tot_ch=orbital(i)%tot_ch+orbital(i)%dm(m,n)*tmp !total charge in the orbital is calculated here
        END DO
     END DO
     IF (i.eq.1) WRITE(*,'(/,5X,"Total charge in the required orbitals:")')
     WRITE(*,'(5X,i5,"|",f20.15)') i,orbital(i)%tot_ch
  END DO
  
  WRITE (*,'(/,5X,a)')'Saving density matrices...'
  
  !save the calculated density matrices to the unit 'dmunit'
  WRITE(dmunit, '(a)') BFTYPE(2) !'Contracted GTO basis'     !type of the basis set corresponding to the saved density matrix elements
  WRITE(dmunit, '(i5)') norb                                 !number of density matrices on this file
  WRITE(dmunit, '(i5)') basis%ncgto                          !dimension of the density matrix
     DO i=1,req_dm
        WRITE(dmunit,'(a,1X,i5)') trim(header_base),i        !title of the current density matrix
        DO m=1,basis%ncgto
           DO n=1,m
              WRITE(dmunit,*) orbital(i)%dm(m,n)             !lower triangle of the symmetric density matrix
           END DO
        END DO
     END DO
  WRITE(dmunit,*)
  CLOSE(dmunit)
  WRITE (*,'(5X,a)')'Density matrices saved successfully.'
  
  END SUBROUTINE CALC_DM

  FUNCTION NOATOMS(IUNIT)
  !finds from the opened molden file (INUIT) the number of atoms the molecule has
  use constants
  IMPLICIT NONE
  INTEGER :: IUNIT, NOATOMS, ifail, i
  CHARACTER(len=line_len) :: line

  CALL SEARCH(1,'[Atoms]',len('[Atoms]'),ifail,.true.)
  IF (ifail .ne. 0) THEN
     STOP "Molden file does not contain information on atoms. Program terminated."
  END IF 
  
  noatoms=0 !get the number of lines which build up the basis set information
  DO
    READ (1,'(a)') line
    IF (line(1:5) .eq. '[GTO]') THEN  !the assumption is that the information on atoms is immediately followed by information on the basis set!
       EXIT
    ELSE
       noatoms=noatoms+1
    END IF
  END DO

  END FUNCTION NOATOMS

  SUBROUTINE GETATOMS(IUNIT)
  !Reads coordinates of atoms, their numbers (indices) and symbols of them from the Molden file (IUNIT)
  use basis_atoms
  !INPUT: the common variable 'basis'
  IMPLICIT NONE
  INTEGER :: IUNIT, i, ifail, iat, k
  CHARACTER(len=line_len) :: line

  CALL SEARCH(1,'[Atoms]',len('[Atoms]'),ifail,.true.)
  DO i=1,size(basis%at_num)
     READ (IUNIT,'(a)') line
     READ(line,*) basis%at_sym(i), basis%at_num(i), iat, (basis%at_coord(k,i),k=1,3)
  END DO

  WRITE (*,'(/,5X,a)')'Coordinates of atoms [Angstroms]: '
  DO i=1,basis%no_at
     WRITE (*,'(5X,i2,1X,a,1X,3f20.10)') i,basis%at_sym(i),(basis%at_coord(k,i),k=1,3)
  END DO

  basis%at_coord(:,:)=basis%at_coord(:,:)*to_au !convert atomic coordinates to a.u.
  
  END SUBROUTINE GETATOMS

  SUBROUTINE GETBASISINFO(IUNIT)
  !returns the number of contracted (basis%ncgto) gaussians and the number of primitive (basis%nprim) gaussians of the basis on the Molden file
  !INPUT: the common variable 'basis'
  use prefactors
  use basis_atoms
  IMPLICIT NONE
  INTEGER :: IUNIT, tmp, ios, i, j, ang_typ, ifail
  INTEGER :: np(1:5), gf(1:5), sph_gtos
  REAL :: one
  CHARACTER(len=1) :: letter
  CHARACTER(len=line_len) :: line

  CALL SEARCH(1,'[GTO]',len('[GTO]'),ifail,.true.)

  basis%bind=0 !get the number of lines which build up the basis set information
  DO
    READ (1,'(a)') line
    IF (line(1:4) .eq. '[MO]') THEN
       EXIT
    ELSE
       basis%bind=basis%bind+1
    END IF
  END DO

  CALL SEARCH(1,'[GTO]',len('[GTO]'),ifail,.true.) !get back on the position of the start of the basis set information
  
  gf(:)=0 !number of s-type contracted gaussians, etc.

  basis%maxall=-1 !maximum number of contractions across all types of functions

  basis%nprim=0  !number of primitives
  
  np(:)=0 !number of primitive s-type functions, etc.
  basis%sph_nprim=0 !total number of primitive spherical gaussians
  basis%sph_ncgto=0 !total number of contracted spherical gaussians
  
  DO j=1,(basis%bind-1)
     READ (1,'(a)') line

     DO ang_typ=1,size(bas_typ) !angular types of basis (s, p, d, f)

        IF (line(2:2) .eq. trim(bas_typ(ang_typ))) THEN
           READ(line,*) letter, tmp, one !basis type, no. of contractions, 1.00
           np(ang_typ)=np(ang_typ)+tmp*bf_num(ang_typ)
           gf(ang_typ)=gf(ang_typ)+1*bf_num(ang_typ)
           basis%sph_nprim=basis%sph_nprim+tmp*(2*(ang_typ-1)+1)  !=... + number of contractions*(2*l+1)
           basis%sph_ncgto=basis%sph_ncgto+1*(2*(ang_typ-1)+1)    !=... + (2*l+1)
           IF (tmp > basis%maxall) basis%maxall=tmp
        END IF

    END DO

  END DO   

  basis%ncgto=sum(gf)
  basis%nprim=sum(np)
  basis%sph_maxall=basis%maxall

  WRITE (*,'(/,5X,a,i5)')'Number of primitive gaussians:  ',basis%nprim
  WRITE (*,'(5X,a,i5)')  'Number of contracted gaussians: ',basis%ncgto
  WRITE (*,'(/,5X,a,i5)')'Number of primitive spherical gaussians:  ',basis%sph_nprim
  WRITE (*,'(5X,a,i5)')  'Number of contracted spherical gaussians: ',basis%sph_ncgto

  END SUBROUTINE GETBASISINFO

  SUBROUTINE GETBASIS(IUNIT)
  !Reads all basis set information and constructs all neccessary index arrays as well
  !INPUT: the common variable 'basis'
  use prefactors
  use basis_atoms
  IMPLICIT NONE

  REAL(kind=wp) :: alp, contr
  INTEGER :: IUNIT, i, j, at_i, tmp, k, c_ind, ind, p, ang_typ, q, r, ifail, l, m, sph_k, sph_c_ind
  REAL :: one
  CHARACTER(len=line_len) :: line, at_char(1:max_at)
  CHARACTER(len=1) :: letter
  CHARACTER(len=1) :: crds(1:3)
  DATA crds/'x','y','z'/

  CALL SEARCH(IUNIT,'[GTO]',len('[GTO]'),ifail,.true.) !get back on the position of the start of the basis set information
  
  k=0
  c_ind=0
  sph_c_ind=0
  sph_k = 0

  DO i=1,basis%no_at
     at_char(i)=''
     WRITE(at_char(i),'(1x,i2,1x,i1)') i, 0 !prepare headers which correspond to the beginnings of basis sets on individual atoms
  END DO
  
  WRITE(*,'(/,5X,"Basis set details follow:")')
  WRITE(*,'(5X,"CARTESIAN GAUSSIANS====================================================||&
       &SPHERICAL GAUSSIANS======================================================")')
  WRITE(*,'(5X,a6,a6,a5,a19,a19,a,a6,a5,a2,a6,a6,a7,a19,a19,a,a6,a5)') &
       'Unct.|','Cont.|','Ang.|',' Gauss. exponent |',' Contraction cf. |','Atom|',&
       'Ct.In|','No.P.','||','Unct.|','Cont.|',' L| M |',' Gauss. exponent |',' Contraction cf. |','Atom|','Ct.In|','No.P.'

  j=0
  DO !loop over all lines containing the basis set information (1 -> basis%bind) 
     READ (IUNIT,'(a)') line
     j=j+1 !we read one line; 'j' counts the number of read lines
     IF (j .eq. basis%bind) EXIT !we have reached the end of the basis set information

     !check if basis set information for another atom starts
     DO i=1,basis%no_at
        IF (trim(line) .eq. trim(at_char(i))) THEN !yes we are on the beginning of a new basis set for an atom 'i'
           at_i=i                                  !current number of atom
        END IF
     END DO
  
     DO ang_typ=1,size(bas_typ) !angular types of basis (1=s, 2=p, 3=d, 4=f, 5=g)

        IF (line(2:2) .eq. trim(bas_typ(ang_typ))) THEN !which angular type are we reading in? s, or p, or ...

        READ(line,*) letter, tmp, one !letter=basis set type, tmp=number of contractions, one=1 (redundant)
        c_ind=c_ind+1  !current index of the contracted cartesian gaussian
        sph_c_ind=sph_c_ind+1 !current index of the contracted spherical gaussian
        DO i=1,tmp     !tmp=number of contracted gaussians following
           READ (IUNIT,'(a)') line
           j=j+1       !we just read one line
           READ(line,*) alp, contr !exponent, contraction coeff
           !construction of the spherical gaussians
           l=ang_typ-1
           DO m=1,2*l+1                            !construct all primitive and contracted GTOs of the current angular type from the read exponent
              sph_k=sph_k+1                        !this counts primitive spherical gaussian functions
              basis%sph_a(sph_k)=alp               !exponent of the gaussian
              basis%l_ang(sph_k)=l                 !l value of the corresponding spherical harmonic
              basis%m_ang(sph_k)=m-l-1             !m value of the corresponding spherical harmonic
              basis%sph_contr_a(sph_c_ind+m-1,i)=sph_k !index of the primitive spherical gaussian funtion building up the current contracted function
              basis%sph_contr(sph_c_ind+m-1,i)=contr   !contraction index of the current contracted spherical gaussian function
              basis%sph_at_ind(sph_k)=at_i             !at which atom is this basis function centred
              basis%sph_no_contr(sph_c_ind+m-1)=tmp    !number of primitives building up the current contraction
           ENDDO
           !construction of the cartesian gaussians
           DO p=0,(bf_num(ang_typ)-1)              !construct all primitive and contracted GTOs of the current angular type from the read exponent
              k=k+1                                !this counts primitive functions 
              basis%prefact(k)=ang_fact(ang_typ,p+1) !prefactor string defining angular behaviour of the current GTO
              basis%ang_exp(k,1)=0                   !extract the integer exponents from the prefactor string assigned above:
              basis%ang_exp(k,2)=0
              basis%ang_exp(k,3)=0
              DO q=1,3                                           !for all coordinates
                 DO r=1,len(basis%prefact(k))                    !scan all characters of the prefactor string
                    IF (basis%prefact(k)(r:r) .eq. crds(q)) THEN !see if the current coordinate is at the current position in the prefactor string
                       basis%ang_exp(k,q)=basis%ang_exp(k,q)+1   !increase the exponent for this coordinate if it is present in the string
                    END IF
                 END DO
              END DO
              basis%alpha(k)=alp                   !exponent of the current primitive basis function
              basis%at_ind(k)=at_i                 !at which atom is this basis function centred
              basis%contr_alpha(c_ind+p,i)=k       !index of the primitive funtion building up the current contracted function
              basis%contractions(c_ind+p,i)=contr  !contraction index of the current contracted function
              basis%no_contr(c_ind+p)=tmp          !number of primitives building up the current contraction
              IF (p+1 .le. 2*l+1) THEN
                 WRITE(*,'(5X,i5,"|",i5,"|",a,"|",d18.10,"|",d18.10,"|",i2,2X,"|",i5,"|",i5,"||",i5,"|",i5,"|",i2," ",i3,"|",&
                 &d18.10,"|",d18.10,"|",i2,2X,"|",i5,"|",i5)') k, c_ind + p, basis % prefact(k), alp, contr, at_i, &
                    basis % contr_alpha(c_ind + p, i), basis % no_contr(c_ind + p), basis % sph_contr_a(sph_c_ind + p, i), &
                    sph_c_ind + p, basis % l_ang(basis % sph_contr_a(sph_c_ind + p, i)), &
                    basis % m_ang(basis % sph_contr_a(sph_c_ind + p, i)), &
                    basis % sph_a(basis % sph_contr_a(sph_c_ind + p, i)), &
                    basis % sph_contr(sph_c_ind + p, i), at_i, &
                    basis % contr_alpha(c_ind + p, i), &
                    basis % no_contr(c_ind + p)
              ELSE
                 WRITE(*,'(5X,i5,"|",i5,"|",a,"|",d18.10,"|",d18.10,"|",i2,2X,"|",i5,"|",i5,"||",a)') &
                    k, c_ind + p, basis % prefact(k), alp, contr, at_i, basis % contr_alpha(c_ind + p, i), &
                    basis % no_contr(c_ind + p), "-------------------------------------------------------------------------"
              ENDIF
           END DO
        END DO
        c_ind=c_ind+(bf_num(ang_typ)-1) !increment contraction index by the number of angular behaviours which were now processed
        sph_c_ind=sph_c_ind+(2*l)   !increment contraction index by the number of angular behaviours of the spherical gaussians which were now processed

        END IF

     END DO

  END DO

  END SUBROUTINE GETBASIS

      SUBROUTINE SEARCH(IUNIT,A,ln,ifail,rew)
!***********************************************************************
!
!     Utility to search a dataset IUNIT for a header A, given the length
!     'ln' of the string A. If rew is set to true, the unit will be rewinded
!***********************************************************************
!
      CHARACTER(len=ln)  A,B
!
      INTEGER, PARAMETER :: IWRITE=6
      INTEGER :: IUNIT, ln, ifail
      LOGICAL :: rew
!
      IF (rew) THEN
         REWIND IUNIT
      END IF

 1    READ(IUNIT,'(a)',END=990) B
!      WRITE(IWRITE,1500) B
!
      IF (B .eq. A) then
        CONTINUE
      ELSE
        GOTO 1
      END IF
!
!..... We reach this point if the search has been successful
!
      ifail = 0
      RETURN
!
!---- Process error condition namely, header not found by end of file.
!
  990 CONTINUE
!
      !WRITE(IWRITE,9900) A,IUNIT
      ifail = 1
      return
!
!---- Format Statements
!
 1500 FORMAT(10X,' Data = ',4A8)
 9900 FORMAT(/' **** Error in SEARCH: ',//, ' Attempt to find header (A) = ',A,' on unit = ',I3, ' has failed.'/)
!
      END SUBROUTINE SEARCH

END MODULE molden

MODULE rmat_rd
!This module contains routines for reading in basis set information from the R-matrix fort.2 file to the structure 'basis'.
!Only the subroutine READ_RMAT is new, all the other routines come from the R-matrix codes.

CONTAINS

  SUBROUTINE READ_RMAT(RMU)
  use basis_atoms
  !reads all basis set information from the R-matrix fort.2 (RMU) file to the common variable 'basis'.
  IMPLICIT NONE
  !INPUT:
  INTEGER :: RMU  !unit number for the R-matrix file containing basis set information (2 is the default)
  !auxiliray:
  integer, parameter :: IWRITE=6,MAXCFN=1000
  integer, parameter :: LTRAN=500
  
  integer :: NPRIMS,nsabf,nsymt,nsymabf,nbft(1:maxsym)
  integer :: i,j,k,ind,counter,LSABFTAB,isymfn,ipntr,itab,icount
  
  character(len=4) :: CNUCNAME(max_at),CSYMBOL(MAXCFN),CANGULAR(MAXCFN)
  integer :: NUCNUMB(max_at),CHARGNUC(max_at)
  integer :: itrans(LTRAN),ics(LTRAN)
  real(8) :: ctrans(LTRAN)

  integer,allocatable :: icontno(:)
  real(8),allocatable :: coef(:)
  integer, allocatable :: lmn(:,:)

  integer, allocatable :: itran(:)
  real(8), allocatable :: ctran(:)
  integer, allocatable :: ISYMABF(:)

  !Print the header from the fort.2 file and get basic information about the basis set
  CALL PRTINF(RMU,IWRITE,basis%nadapt,ICs,ITRANs,CTRANs,basis%nsymt,basis%nbft)
  !
  !     NSABF counts the total number of symmetry adapted basis
  !     functions here.
  !
  basis%nsabf = 0
  DO I=1,basis%nsymt
     basis%nsabf = basis%nsabf+basis%nbft(I)
  enddo
  ALLOCATE(itran(basis%nadapt),ctran(basis%nadapt),isymabf(basis%nadapt))

  !Read the nuclear information
  REWIND RMU
  READ(RMU)
  READ(RMU)
  READ(RMU) basis%no_at  !number of atoms

  IF(basis%no_at.GT.max_at)THEN
     WRITE(IWRITE,9900)
     WRITE(IWRITE,9910) basis%no_at,max_at
     STOP 999
  END IF

9900 FORMAT(//' **** Error in READORBS: ',//)
9910 FORMAT(10X,'Parameter value of MAXNUC is too small in this case ',&
          10X,'NNUC - number of nuclei = ',I5,/,&
          10X,'MAXNUC limit in code    = ',I5,/,&
          10X,'Recompilation is needed',//)

  ALLOCATE(basis%at_sym(1:basis%no_at),basis%at_coord(1:3,1:basis%no_at),basis%at_num(1:basis%no_at))
  !read cooridinates of atoms and their symbols
  DO I=1,basis%no_at
     READ(RMU) basis%at_sym(i),basis%at_num(i),basis%at_coord(1,i),basis%at_coord(2,i),basis%at_coord(3,i),CHARGNUC(I)
  END DO
  !
  !---- Read the Gaussian basis function information from the integrals
  !     file now.
  !
  !       NPRIMS - Number of primitive Gaussian functions
  !       NCGTO - Number of contracted Gaussian functions
  !
  READ(RMU) basis%nprim,basis%ncgto !this can possibly contain scattering functions as well - add reordering accoring to raddensity?

  !now allocate arrays for the basis set
  ALLOCATE(basis % alpha(1 : basis % nprim), &
           basis % contractions(1 : basis % ncgto, 1 : basis % nprim), &
           basis % at_ind(1 : basis % nprim), &
           basis % prefact(1 : basis % nprim), &
           basis % contr_alpha(1 : basis % ncgto, 1 : basis % nprim), &
           basis % ang_exp(1 : basis % nprim, 1 : 3), &
           basis % overlaps_cgto(1 : basis % ncgto, 1 : basis % ncgto), &
           basis % no_contr(1 : basis % ncgto))

  basis % alpha(:) = 0.0_wp
  basis % contractions(:,:) = 0.0_wp
  basis % at_ind(:) = 0
  basis % prefact(:) = ''
  basis % contr_alpha(:,:) = 0
  basis % ang_exp(:,:) = 0
  basis % overlaps_cgto(:,:) = 0.0_wp
  basis % no_contr(:) = 0
  
  READ(RMU) (basis%no_contr(i),i=1,basis%ncgto)  !number of primitives per contraction

  !
  !...... Allocate storage space for the primitive basis function table 
  !
  allocate (icontno(basis%nprim),coef(basis%nprim),lmn(3,basis%nprim))

  !
  !...... Invoke the routine which reads the table of primitive Gaussians
  !       and their identifications into core
  !
  CALL GETAB(IWRITE,RMU,basis%ncgto,basis%no_contr,icontno,csymbol,basis%at_ind,&
       basis%prefact,basis%no_at,basis%at_sym,basis%alpha,coef,&
       basis%nsabf,isymabf,itran,ctran,basis%nadapt,LSABFTAB)

  !---- Each primitive Gaussian has associated with it an angular term in
  !     the form x^i * y^j * z^k and this is represented in MOLECULE by
  !     a character string. 
  !
  DO I=1,basis%nprim
     CALL INDEX(basis%prefact(i),basis%ang_exp(i,1),basis%ang_exp(i,2),basis%ang_exp(i,3))
  enddo

  !finally transform the arrays read from fort.2 to the convention used here

  !mapping primitive -> contracted
  ind=0
  DO i=1,basis%nprim
     IF (icontno(i) .eq. ind) THEN !are we still processing the same contracted function 'icontno(i)'?
        counter=counter+1
        !index of a primitive function which is a part of the contracted function 'icontno(i)' as contraction number 'counter'
        basis%contr_alpha(icontno(i),counter)=i     
        !contraction coefficient number 'counter' in contracted function 'icontno(i)'
        basis%contractions(icontno(i),counter)=coef(i)
     ELSE !new contracted function information starts
        counter=1
        ind=icontno(i)
        basis%contr_alpha(icontno(i),counter)=i
        basis%contractions(icontno(i),counter)=coef(i)
     END IF
  END DO

  !mapping contracted -> symmetry adapted
  ALLOCATE(basis%no_sabf(1:basis%nsabf),basis%sa_contr(1:basis%nsabf,1:basis%ncgto),basis%contr_ind(1:basis%nsabf,1:basis%ncgto))
  ind=0
  DO i=1,basis%nadapt
     IF (isymabf(i) .eq. ind) THEN !are we still processing the same symmetry adapted function 'icontno(i)'?
        counter=counter+1
        !the number ('counter') of contracted gaussians forming the symmetry adapted function 'isymabf(i)'
        basis%no_sabf(isymabf(i))=counter
        !index number 'counter' of the contracted gaussian in linear combination of contracted gaussians 
        !corresponding to the symmetry adapted function 'isymabf(i)'
        basis%contr_ind(isymabf(i),counter)=itran(i)
        !coefficient in linear combination of contracted gaussians
        basis%sa_contr(isymabf(i),counter)=ctran(i)
     ELSE !new symmetry adapted function information starts
        counter=1
        ind=isymabf(i)
        basis%no_sabf(isymabf(i))=counter
        basis%contr_ind(isymabf(i),counter)=itran(i)
        basis%sa_contr(isymabf(i),counter)=ctran(i)
     END IF
  END DO
 
  END SUBROUTINE READ_RMAT

      SUBROUTINE PRTINF(ITAPE,IWRITE,MAXTRAN,IC,ITRAN,CTRAN,&
           NSYMT,NBFT)
!***********************************************************************
!
!     Prints the header from the dataset of Gaussian integrals that
!     are computed by the code MOLECULE - this means that all data
!     prior to the records of integrals is printed. Note that some data
!     is returned to the caller here so that this is more than just a
!     utility routine.
!
!     Input data:
!          ITAPE  Logical unit for the MOLECULE output
!         IWRITE  Logical unit for the printer
!             IC  Workspace array for no. of primitives per contracted
!                 Gaussian function
!          ITRAN  Workspace array for integers
!          CTRAN  Workspace array for real*8 variables
!
!     Output data:
!           NSYMT Number of symmetries (IRRs) in this molecular point
!                 group.
!            NBFT Number of basis functions per symmetry
!
!***********************************************************************
      IMPLICIT double precision (a-h,o-z)
!
      character(len=4), allocatable :: cnucnam(:)
      CHARACTER(LEN=8)    LABEL(4),CL,CNAKO,CNAME
      CHARACTER(LEN=192)  CHEADER
      integer :: l,nmax,k,itran,maxtran,kb,ii,ic,iabas,jend,ibbas,j,nnuc
      integer :: idosph,i,nbft,itape,nsymt,iwrite
      real(8) :: ctran,cx,xx,ch,x,potnuc
      
!
      DIMENSION  IC(*),ITRAN(*),CTRAN(*),NBFT(*),X(3)
      double precision, allocatable :: ctran1(:)
!
      DATA CL/'        '/,CNAKO/'        '/
!
      WRITE(IWRITE,1000)
!
!---- Prepare the output file of tail integrals.
!
      open(unit=ITAPE,form='unformatted')
!
!---- Data on this record is as follows:
!
!      CHEADER - Character header from input to MOLECULE code
!        NSYMT - No. of Irreducible Representations (Symmetries) in the
!                set of basis functions
!         NBFT - No. of atomic basis functions per symmetry
!       POTNUC - Nuclear potential energy
!       IDOSPH - Flag stating whether or not we have spherical harmonic
!                basis
!
      READ(ITAPE) CHEADER,NSYMT,(NBFT(I),I=1,NSYMT),POTNUC,IDOSPH
!
      WRITE(IWRITE,500) CHEADER(1:65)
      WRITE(IWRITE,510) NSYMT,(I,NBFT(I),I=1,NSYMT)
      WRITE(IWRITE,520) POTNUC
      if(IDOSPH.EQ.1) then
         ! spherical harmonics basis
         WRITE(IWRITE,530) IDOSPH
      else
         ! cartesian basis
         WRITE(IWRITE,540) IDOSPH
      endif
!
!=======================================================================
!
!     1st set of data is on the nuclear configuration
!
!=======================================================================
!
      READ(ITAPE) (LABEL(I),I=1,4)
      READ(ITAPE)  NNUC
!
      WRITE(IWRITE,2000) (LABEL(I),I=1,4)
      WRITE(IWRITE,2010) NNUC
!
!...... Loop over nuclei - obtain positions and charges
!
      WRITE(IWRITE,2015)
!
      DO I=1,LEN(CNAME)
         CNAME(I:I)=' '
      END DO
!
      DO I=1,NNUC
         
         READ(ITAPE)    CNAME(1:4),II,(X(J),J=1,3),CH
         WRITE(IWRITE,2020) CNAME,II,(X(J),J=1,3),CH
         
      END DO
!
!=======================================================================
!
!     2nd set of data is on Primitive and Contracted Gaussian functions
!
!=======================================================================
!
!---- IBBAS is the number of primitive functions; IABAS the contracted
!
      READ(ITAPE)    IBBAS,IABAS
!
      WRITE(IWRITE,2500) IBBAS,IABAS
!
!---- Array IC has one entry for each contracted function. It gives the
!     number of primitives in that contraction.
!
      READ(ITAPE)    (IC(I),I=1,IABAS)
!
!---- For each contracted function we loop over the number of primitives
!     within it and obtain data on each.
!
      DO I=1,IABAS
         JEND=IC(I)

         DO  J=1,JEND
            READ(ITAPE)    CL(1:4),KB,CNAKO(1:4),XX,CX
         enddo
      enddo
!
!=======================================================================
!
!     3rd set of data is on symmetry information
!
!=======================================================================
!
      READ(ITAPE) (LABEL(I),I=1,4)
      READ(ITAPE) IABAS
!
      maxtran = 0
      DO I=1,IABAS
         READ(ITAPE)     J,(ITRAN(K),CTRAN(K),K=1,J)
         maxtran = maxtran+j
      enddo
!
!=======================================================================
!
!     4th set of data is on MULLIKEN Population information
!
!     This contains information on the nuclear charges and on the
!     angular behaviour of the symmetry adaped basis functions.
!
!=======================================================================
!
      READ(ITAPE) (LABEL(I),I=1,4)
      READ(ITAPE) NMAX,IABAS
!
      allocate (cnucnam(nmax),ctran1(nmax))

      READ(ITAPE) NMAX,(CNUCNAM(I),CTRAN1(I),I=1,NMAX)

      DO I=1,IABAS
         READ(ITAPE)    J,K,L
      enddo
      
      deallocate (cnucnam,ctran1)
!
      RETURN
!
!---- Format Statements
!
  500 FORMAT('',5X,'Header Card (1:65) : ',A,/,'')
  510 FORMAT('',5X,'No. of symmetries in basis set = ',I3,/,'',/,&
           '',5X,'No. of basis functions per symmetry: ',/,&
           ('',5X,I2,1X,I3))
  520 FORMAT('',/,'',5X,'Nuclear Potential Energy = ',F15.7,&
           ' (Hartrees) ')
  530 FORMAT ('',/,'',5X,'Spherical Basis, IDOSPH=',I12,/,'')
  540 FORMAT ('',/,'',5X,'Cartesian Basis, IDOSPH=',I12,/,'')
!
 1000 FORMAT('',/,'',/,'',/,'',10X,&
           'Atomic Integrals File: Header Records ',/,'')
!
 2000 FORMAT('',5X,'Section 1 Header = ',4A,/,'')
 2010 FORMAT('',5X,'Number of nuclear centers = ',I3,/,'')
 2015 FORMAT('',5X,' Symbol ',1X,'No.',' (X,Y,Z) Co-ordinates ',T50,&
           ' Charge ',/,'')
 2020 FORMAT('',5X,A,1X,I2,1X,3(F10.6,1X),F6.3)
 2500 FORMAT('',/,'',/,'',5X,&
           'Number of primitive  Gaussian functions = ',I4,/,'',&
           5X,'Number "" contracted  " " "    " " " "  = ',I4,/,'',/,'')
 2550 FORMAT(5X,'Primitive functions per contracted function:',/)
 2560 FORMAT((5X,8(I3,'.',1X,I3,1X)))
 2600 FORMAT(/,5X,'A N A L Y S I S  of  each  C O N T R A C T I O N',/)
 2610 FORMAT(/,5X,'Contracted function no. = ',I3,/,&
           5X,'No. of primitives       = ',I3,/)
 2620 FORMAT(1X,'Primitive = ',I2,1X,'CL = ',A,1X,'KB = ',I3,1X,&
           'CNAKO = ',A,1X,'XX = ',D12.5,1X,'CX = ',D12.5,1X)
!
 3000 FORMAT(/,5X,'Section 2 Header = ',4A,/)
 3010 FORMAT(5X,'Number of symmetry adapted functions = ',I3,/)
 3020 FORMAT('',/,'',5X,'Symmetry function no. = ',I3,/,&
           5X,'No. of components     = ',I3,/,&
           5X,' (ITRAN,CTRAN) pairs: ',//,&
           (5X,4(I4,1X,D12.5,1X)))
!
 4000 FORMAT(/,5X,'Section 3 Header = ',4A,/)
 4010 FORMAT(5X,'Number of atomic nuclei (NMAX) = ',I3,/,'',/,&
           5X,'Number ""   ""   basis functions (IABAS) = ',I5,/,'',/)
 4020 FORMAT(5X,'NMAX = ',I4,' Nuclei/Charge  pair definitions follow: ',/)
 4030 FORMAT((5X,3(I3,'.',1X,A,1X,D12.5,1X)))
 4040 FORMAT(/,5X,'No. of items (IABAS) = ',I5,/)
 4050 FORMAT(5X,'Item = ',I3,' J = ',I3,' K = ',I3,' L = ',I3,1X,&
           'Angular function = ',A)
!  
      END SUBROUTINE PRTINF


SUBROUTINE GETAB(IWRITE,ITAPE,NCONTRAS,NPRCONT,&
     ICONTRNO,CSYMBOL,ngnuc,CANGULAR,&
     nnuc,cnucname,XPONENT,XCONTCOF,&
     NSYMABF,ISYMABF,ITRAN,CTRAN,LTRAN,LSABFTAB)
!***********************************************************************
!
!     GETAB - GeT primitive and symmetry adapted Gaussian TABles
!
!     This routine reads data from the atomic integrals file into a
!     tabular format that is used later by the code. The first table
!     relates primitive Gaussian functions to contracted ones while the
!     second relates contracted Gaussian functions to symmetry adapted
!     basis functions.
!
!     Input data:
!         IWRITE Logical unit for the printer
!          ITAPE Logical unit holding the data on the functions
!                i.e. Output by the MOLECULE code
!       NCONTRAS Number of contracted Gaussian functions
!        NPRCONT Table of the no. of primitives per contracted function
!          LTRAN Size of the arrays ISYMABF,ITRAN and CTRAN
!
!     Output data:
!        ICONTRNO Sequence number of the contracted Gaussian function
!                 to which this primitive belongs
!         CSYMBOL Unique symbol which identifies the symmetry independen
!                 atom to which this basis function belongs
!        CANGULAR Designation of the angular term for this primitive
!                 Gaussian
!         XPONENT The exponent for this primitive Gaussian function
!        XCONTCOF The contraction coefficient for this primitive Gaussia
!         NSYMABF Number of symmetry adapted basis functions
!         ISYMABF Pointer array defining symmetry adapted basis
!                 function numbers
!           ITRAN Pointer to contracted Gaussian function which belongs
!                 to the symmetry adapted basis function in the
!                 corresponfing entry in ISYMABF
!           CTRAN Coefficient partnering the corresponding entries in
!                 ISYMABF and ITRAN.
!        LSABFTAB Length of the symmetry adapted basis function table
!
!     Notes:
!
!        The contracted Gaussian expansions are normalized after
!     being read in. This is because the output from MOLECULE may be
!     un-normalized.
!
!     Linkage:
!
!        INDEX, NORM
!
!***********************************************************************
  IMPLICIT double precision (a-h,o-z)
!
  CHARACTER(LEN=4)  CANGULAR,csymbol,cnucname
  integer :: nnuc,ltran,k,jnlast,ilast,i,ncontras,j,nprcont,itran
  integer :: icontrno,itape,isymbnum,jnuc,inuc,mult,jnc,ngnuc
  integer :: lsabftab,iwrite,isymabf,kbias,nsymabf,l,jcode
  integer :: ixpower,iypower,izpower
  real(8) :: ctran,xponent,xcontcof
  
!
  DIMENSION  NPRCONT(*),cnucname(nnuc),ngnuc(*)
  DIMENSION  ICONTRNO(*),CSYMBOL(*),CANGULAR(*)
  DIMENSION  XPONENT(*),XCONTCOF(*)
  DIMENSION  ISYMABF(*),ITRAN(ltran),CTRAN(ltran)
!
!---- We assume that ITAPE is positioned so that we may begin reading.
!
      K = 0
      jnlast = 0
      ilast = 0
      DO 10 I=1,NCONTRAS
      DO 20 J=1,NPRCONT(I)
      K=K+1
      ICONTRNO(K)=I
      READ(ITAPE,ERR=900) CSYMBOL(K),ISYMBNUM,&
           CANGULAR(K),XPONENT(K),XCONTCOF(K)
!
      if(j.eq.1) then
! I THINK THIS IS JUST COSMETIC
        if(i.gt.1) then
         if(csymbol(k).ne.csymbol(k-1)) then
          ilast = i-1
          jnlast = jnuc
         endif
        end if
        do 21 inuc=1,nnuc
        if(csymbol(k).eq.cnucname(inuc)) jnuc=inuc
 21     continue
!
!---- Must check for 'equivalent atoms' with non-unique labels.
        if(jnuc.ne.jnlast.and.jnuc.ne.jnlast+1) then
          mult = jnuc-jnlast
          jnc = mult-1-mod(i-ilast-1,mult)
          jnuc = jnuc-jnc
        endif
      endif
      ngnuc(k) = jnuc
 20   CONTINUE
!
!.... Build the JCODE value by looking at the angular terms and
!     then invoke the normalization routine for this contraction.
!
!        At this stage K is a high water marker and will point to the
!        last primitive read. SInce all primitives within a contraction
!        are of the same angular behaviour then we may use this to work
!        out the JCODE value. However a new pointer must be set to
!        give access to the exponents and coefficients.
!
      CALL INDEX(CANGULAR(K),IXPOWER,IYPOWER,IZPOWER)
      JCODE=IXPOWER+IYPOWER+IZPOWER+1
!
      L=K-NPRCONT(I)+1
 10   CONTINUE
!
!---- Having read the contracted Gaussian information we must now
!     skip the header that follows it and read the symmetry adapted
!     basis function information. It defines linear combinations of
!     contracted Gaussians which make symmetry adapted functions
!
      READ(ITAPE,ERR=910)
      READ(ITAPE,ERR=920) NSYMABF
!
      KBIAS=0
!
      DO 50 I=1,NSYMABF
      READ(ITAPE,ERR=930) J,(ITRAN(KBIAS+K),CTRAN(KBIAS+K),K=1,J)
      DO 60 K=1,J
      ISYMABF(KBIAS+K)=I
 60   CONTINUE
!
      KBIAS=KBIAS+J
!
      IF(KBIAS.GT.LTRAN)THEN
        WRITE(IWRITE,9900)
        WRITE(IWRITE,9950) KBIAS,LTRAN
        STOP 999
      END IF
 50   CONTINUE
!
!---- Print out the symmetry adapted basis function table
!
      WRITE(IWRITE,5000)
      
      DO K=1,KBIAS
         WRITE(IWRITE,5010) K,ISYMABF(K),ITRAN(K),CTRAN(K)
      END DO
!
!---- Store the length of the table for later use in the transform
!     code.
!
      LSABFTAB=KBIAS
!
      RETURN
!
!---- Error condition handlers
!
!..... Error during of a record of primitive basis function data
!
  900 CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9910) K
      STOP 999
!
!..... Error reading header of SYMTRANS data
!
  910 CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9920)
      STOP 999
!
!..... Error reading number of symmetry adapted functions
!
  920 CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9930)
      STOP 999
!
!..... Error reading a record for one symmetry adapted function
!
 930  CONTINUE
!
      WRITE(IWRITE,9900)
      WRITE(IWRITE,9940) I
      STOP 999
!
!---- Format Statements
!
 1000 FORMAT(//,15X,'====> GETAB - READ PRIMITIVE FN. TABLE <====',/)
 1010 FORMAT(/,15X,'Logical unit for atomic ints     (ITAPE) = ',I5,/,&
           15X,'Total number of contracted fns. (NCONTRAS) = ',I5)
!
 5000 FORMAT('',/,'',10X,'Symmetry Adapted Function Definition Table',/,&
           '',10X,'==========================================',/,'',/,&
           '',10X,'This relates contracted Gaussian functions',/,&
           '',10X,'to symmetry adapted basis functions for the',/,&
           '',10X,'molecular point group. ',/,'',/,'',/,&
           '',10X,'Row  Basis Function  Cont. Gaussian    Coefficient',/,&
           '',10X,'---  --------------  --------------    -----------',&
           /,'',/,'')
 5010 FORMAT('',5X,I7,3X,I10,8X,I10,7X,F10.6,1X)
 9900 FORMAT(//,'',10X,'**** Error in GETAB: ',//)
 9910 FORMAT('',10X,'Reading of record number = ',I5,' has failed ',/)
 9930 FORMAT('',10X,'Reading NSYMABF has failed ',/)
 9920 FORMAT('',10X,'Reading header for SYMTRANS has failed',/)
 9940 FORMAT('',10X,'Reading of SYMTRANS record no. = ',I5,' failed',/)
 9950 FORMAT('',10X,'Not enough space in the symmetry adapted function',&
           ' table.',//,&
           '',10X,'Need at present = ',I10,/,&
           '',10X,'Given value     = ',I10,/)
!
      END SUBROUTINE GETAB

      SUBROUTINE INDEX(CANGULAR,IXPOWER,IYPOWER,IZPOWER)
!***********************************************************************
!
!     INDEX - convert angular string to INDEX values
!
!     Takes the string representing the powers of x,y and z in the
!     Cartesian Gaussian generated by MOLECULE and builds integer
!     representations of these.
!
!     Input data:
!       CANGULAR the string representation of the angular part of the
!                primitive Gaussian. eg 'S       ' means that there is
!                no angular part.
!
!     Output data:
!         IXPOWER The power associated with the X co-ordinate
!         IYPOWER The power associated with the Y co-ordinate
!         IZPOWER The power associated with the Z co-ordinate
!
!     Notes:
!
!     This routine is essentially a lexical analyser for the output from
!     MOLECULE. Read the data cards in subroutine READIN in MOLECULE to
!     understand what is being generated there. It is READIN that
!     generates the data that MOLECULE writes and that we are analysing
!     here.
!
!     Author: Charles J Gillan
!
!     Copyright (c) 1995 Charles J Gillan
!     All rights reserved
!
!***********************************************************************
!
      CHARACTER(LEN=4)  CANGULAR,BUFF
      integer :: iypower,ixpower,izpower
!
!---- Initialization to case for 'S'
!
      IXPOWER=0
      IYPOWER=0
      IZPOWER=0
!
!---- Case 1. 'S' type primitive Gaussian
!
      IF(CANGULAR(1:1).EQ.'S'.AND.CANGULAR(2:4).EQ.'   ') return
!
!---- Case 2. 'X' type primitive Gaussian
!
      IF(CANGULAR(1:1).EQ.'X'.AND.CANGULAR(2:4).EQ.'   ')THEN
        IXPOWER=1
        return
      END IF
!
!---- Case 3. 'Y' type primitive Gaussian
!
      IF(CANGULAR(1:1).EQ.'Y'.AND.CANGULAR(2:4).EQ.'   ')THEN
        IYPOWER=1
        return
      END IF
!
!---- Case 4. 'Z' type primitive Gaussian
!
      IF(CANGULAR(1:1).EQ.'Z'.AND.CANGULAR(2:4).EQ.'   ')THEN
        IZPOWER=1
        return
      END IF
!
!---- Case 5. 'XX' type primitive Gaussian
!
      IF(CANGULAR(1:2).EQ.'XX'.AND.CANGULAR(3:4).EQ.'  ')THEN
        IXPOWER=2
        return
      END IF
!
!---- Case 6. 'XY' type primitive Gaussian
!
      IF(CANGULAR(1:2).EQ.'XY'.AND.CANGULAR(3:4).EQ.'  ')THEN
        IXPOWER=1
        IYPOWER=1
        return
      END IF
!
!---- Case 7. 'XZ' type primitive Gaussian
!
      IF(CANGULAR(1:2).EQ.'XZ'.AND.CANGULAR(3:4).EQ.'  ')THEN
        IXPOWER=1
        IZPOWER=1
        return
      END IF
!
!---- Case 8. 'YY' type primitive Gaussian
!
      IF(CANGULAR(1:2).EQ.'YY'.AND.CANGULAR(3:4).EQ.'  ')THEN
        IYPOWER=2
        return
      END IF
!
!---- Case 9. 'YZ' type primitive Gaussian
!
      IF(CANGULAR(1:2).EQ.'YZ'.AND.CANGULAR(3:4).EQ.'  ')THEN
        IYPOWER=1
        IZPOWER=1
        return
      END IF
!
!---- Case 10. 'ZZ' type primitive Gaussian
!
      IF(CANGULAR(1:2).EQ.'ZZ'.AND.CANGULAR(3:4).EQ.'  ')THEN
        IZPOWER=2
        return
      END IF
!
!---- Case 11. l=3 and higher (use internal i/o to convert)
!
      write(buff,2000) cangular
 2000 format(a4)
      read (buff,2001) ixpower,iypower,izpower
 2001 format(1x,3i1)
!
      RETURN
!
      END SUBROUTINE INDEX

!      SUBROUTINE READ_RMAT_DM(bas)
!      use constants
!      use precisn
!      use basis_atoms
!      IMPLICIT NONE
!      INTEGER :: i
!      TYPE(basis_type) :: bas
!      REAL(kind=wp), ALLOCATABLE :: denmat(:)
!
!!---- Arguments of the subroutine DENRD
!      character(len=80)  CHEAD
!      INTEGER idenget, NFTDEN, LDENMAT, ICODE, maxden, NNUC, NSYM, NCI, NCJ, NOB(1:maxsym), NREC
!      INTEGER ISPINI,ISZI,LAMDAI,ISPINJ,ISZJ,LAMDAJ,IGUI,IGUJ,NOCSFI,NOCSFJ,IREFLI,IREFLJ,ISYMTYPI,ISYMTYPJ
!      double precision ECI,ECJ
!      double precision RGEOM(3,max_at)
!
!!
!!---- Compute absolute maximum value for the size of the density
!!     matrix
!!
!      MAXDEN=0
!      DO i=1,bas.nsymt
!         MAXDEN=MAXDEN + bas.nbft(i)*bas.nbft(i)
!      END DO   
!!
!      allocate(denmat(maxden))
!
!      NFTDEN=60
!      IDENGET=1
!
!      CALL DENRD(NFTDEN,IDENGET,DENMAT,LDENMAT,&
!     &                 ICODE,CHEAD,&
!     &                 ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI,&
!     &                 ISPINJ,ISZJ,LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ,&
!     &                 NSYM,NOB,NOCSFI,NOCSFJ,NCI,NCJ,ECI,ECJ,&
!     &                 NNUC,RGEOM,maxden)
!
!      deallocate(denmat)
!
!      END SUBROUTINE READ_RMAT_DM
!
!
!!!Taken from denprop.f, release-0.1
!      SUBROUTINE DENRD(NFTDEN,IDENGET,DENMAT,LDENMAT,&
!     &                 ICODE,CHEAD,&
!     &                 ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI,&
!     &                 ISPINJ,ISZJ,LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ,&
!     &                 NSYM,NOB,NOCSFI,NOCSFJ,NCI,NCJ,ECI,ECJ,&
!     &                 NNUC,RGEOM,maxden)
!!***********************************************************************
!!
!!     DENRD - Reads a Density Matrix from the library file on
!!             unit NFTDEN
!!
!!     Input data:
!!         NFTDEN Fortran logical unit number for the density matrix
!!                library dataset
!!        IDENSET Set number, on NFTDEN, at which this data will be
!!                written
!!         IWRITE Logical unit for the printer
!!
!!     Output data:
!!         DENMAT The complete density matrix to be written out.
!!        LDENMAT Total number of elements in the density matrix
!!          ICODE Code defining the symmetry make up of the density
!!                matrix:
!!                = 1 means a state with itself i.e. ground state
!!                    hence delta lambda=0 and matrix is triangular
!!                = 2 means that states are of the same symmetry but
!!                    are not the same state - matrix triangular.
!!                = 3 means that states have different symmetries and
!!                    so delta lambda is not zero.
!!          CHEAD Character header describing the data
!!         ISPINI 2*S+1 quantum number for first state (J for second)
!!           ISZI Z-projection of S for first state    (J for second)
!!         LAMDAI Z-projection of angular momentum (I=first,J=second)
!!           IGUI G/U quantum number (if any)      (I=first,J=second)
!!         IREFLI Reflection symmetry (if any)     (I=first,J=second)
!!       ISYMTYPI Linear/Abelian flag for the molecular point group
!!           NSYM Number of symmetries in the orbital set (C-inf-v)
!!            NOB No. of orbitals per C-inf-v symmetry
!!         NOCSFI No. of CSFs defining wavefunction of first state
!!         NOCSFJ No. of CSFs  " " "    " " " " " " "  second state
!!            NCI Which CI vector from the Hamiltonian was used for
!!                the first state
!!            NCJ Which CI vector from the Hamiltonian was used for
!!                the second state.
!!            ECI Absolute energy, in Hartrees, of the wavefunction I
!!            ECJ   " "     " "    "    " " "   "   "    "  "  "    J
!!           NNUC Number of nuclei in the system
!!          RGEOM Nuclear configuration at which this density matrix was
!!                generated: (X,Y,Z) co-ordinates for each nucleus.
!!         IWRITE Logical unit for the printer
!!
!!     Linkage:
!!
!!         DENGET
!!
!!     Note:
!!
!!        The format of each member of the density matrix library is
!!
!!     Record 1:   Header defining the library member
!!
!!     Record 2:   Symmetry data for first wavefunction
!!
!!                 ISPIN,ISZ,LAMDA,IGU,IREFL,ISYMTYPI
!!
!!     Record 2a:  Symmetry data for the second wavefunction
!!
!!     Record 3:   Orbital set, CSF and eigenvector information
!!
!!                 NSYM,NOB,NOCSFI,NOCSFJ
!!
!!     Record 4:   Geometry information
!!
!!                 NNUC
!!
!!                 followed by NNUC records of the form
!!
!!                 RGEOM(1,I),RGEOM(2,I),RGEOM(3,I)
!!
!!     Record 5:   Density matrix elements
!!
!!                 NCI,NCJ,ECI,ECJ,LDENMAT,(DENMAT(I),I=1,LDENMAT)
!!
!!***********************************************************************
!      IMPLICIT double precision (a-h,o-z)
!!
!      character(len=80)  CHEAD
!!
!      PARAMETER IDENKEY=60
!      PARAMETER IWRITE=6
!      DIMENSION  DENMAT(maxden),RGEOM(3,*)
!      INTEGER j, idenset, ireset
!!
!!---- Arguments of the subroutine
!      INTEGER idenget, NFTDEN, LDENMAT, ICODE, maxden, NNUC, NSYM, NCI, NCJ, NOB(*), NREC
!      INTEGER ISPINI,ISZI,LAMDAI,ISPINJ,ISZJ,LAMDAJ,IGUI,IGUJ,NOCSFI,NOCSFJ,IREFLI,IREFLJ,ISYMTYPI,ISYMTYPJ
!      double precision ECI,ECJ
!!---- 
!      save idenset
!!
!!---- Position the unit in order to read the set number IDENGET.
!!
!      ireset = 0
!!
! 42   if(idenget.le.1.or.ireset.ne.0) then
!        if(idenget.le.1) idenset=1
!!
!        CALL DENGET(NFTDEN,IDENSET,IDENKEY,'UNFORMATTED')
!!
!!---- Record 1:
!!
!!......... The header record
!!
!        READ(NFTDEN) I,IDENSET,NREC
!        READ(NFTDEN) icode,CHEAD
!!
!!---- Record 2 and possibly 2a:
!!
!!......... Data defining the symmetries of the state(s) involved
!!
!        READ(NFTDEN) ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI
!!
!        IF(ICODE.EQ.3)THEN
!          READ(NFTDEN) ISPINJ,ISZJ,LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ
!        ELSE
!          ISPINJ=ISPINI
!          ISZJ=ISZI
!          LAMDAJ=LAMDAI
!          IGUJ=IGUI
!          IREFLJ=IREFLI
!          ISYMTYPJ=ISYMTYPI
!        END IF
!!
!!---- Record 3:
!!
!!.......... Data defining the orbital set and the CSF expansions
!!           for each of the states
!!
!        READ(NFTDEN) NSYM,(NOB(I),I=1,NSYM),NOCSFI,NOCSFJ
!!
!!---- Record 4:
!!
!!.......... The nuclear geometry information
!!
!        READ(NFTDEN) NNUC
!!
!        DO 40 I=1,NNUC
!        READ(NFTDEN) (RGEOM(J,I),J=1,3)
!  40    CONTINUE
!!
!      ireset = 0
!      endif
!!
!!---- Record 5:
!!
!!.......... The density matrix itself
!!
!      if(idenget.eq.0) return
!!
!      ldenmat = 0
!      READ(NFTDEN,err=41) icode,NCI,NCJ,ECI,ECJ,LDENMAT,&
!     &                    (DENMAT(I),I=1,LDENMAT)
!!
!!---- Print the data that has just been read.
!!
!        WRITE(IWRITE,3000)
!        WRITE(IWRITE,3010) NFTDEN,IDENSET,NREC,ICODE
!        WRITE(IWRITE,3020) CHEAD(1:40)
!        WRITE(IWRITE,3030) CI,ISPINI,ISZI,LAMDAI,IGUI,IREFLI,ISYMTYPI
!        WRITE(IWRITE,3030) CJ,ISPINJ,ISZJ,LAMDAJ,IGUJ,IREFLJ,ISYMTYPJ
!        WRITE(IWRITE,3040) NSYM,(NOB(I),I=1,NSYM)
!        WRITE(IWRITE,3050) NOCSFI,NOCSFJ,NCI,NCJ,ECI,ECJ
!        WRITE(IWRITE,3060) NNUC
!        DO 1 J=1,NNUC
!           WRITE(IWRITE,3070) J,(RGEOM(I,J),I=1,3)
!   1    CONTINUE
!!
!      RETURN
!!
! 41   ireset = idenset
!      idenset = ireset+1
!      backspace nftden
!      go to 42
!!
!!---- Format Statements
!!
! 3000 FORMAT(///,10X,'Data read from Density Matrix Library',//)
! 3010 FORMAT(10X,'Logical unit number for the library   = ',I5,/,&
!     &       10X,'Density matrix set number to be read  = ',I5,/,&
!     &       10X,'Number of records in the set          = ',I5,/,&
!     &       10X,'Format code for type of matrix        = ',I5,//)
! 3020 FORMAT(10X,'Set character header (1:40)           = ',A)
! 3030 FORMAT(/,10X,'Symmetry details for wavefunction ',A,/,&
!     &       10X,'-----------------------------------',//,&
!     &       10X,'Total spin quantum number  = ',I5,/,&
!     &       10X,'Z-projection of spin       = ',I5,/,&
!     &       10X,'Lambda or Irreducible Rep. = ',I5,/,&
!     &       10X,'G/U flag for D-inf-h only  = ',I5,/,&
!     &       10X,'Sigma reflection (C-inf-v) = ',I5,/,&
!     &       10X,'C-inf-v or Abelian flag    = ',I5,/)
! 3040 FORMAT(/,10X,'Orbital Set details: ',/,&
!     &       10X,'-------------------- ',//,&
!     &       10X,'No. of symmetries in the set = ',I5,/,&
!     &       10X,'Orbital per symmetry = ',20(I3,1X),/)
! 3050 FORMAT(/,10X,'CSF details: ',/,&
!     &       10X,'-----------  ',//,&
!     &       10X,'No. of CSFs in wavefunction I = ',I5,/,&
!     &       10X,'No. of CSFs in wavefunction J = ',I5,/,&
!     &       10X,'Root number used for wfn I    = ',I5,/,&
!     &       10X,' ""   """    ""   "  wfn J    = ',I5,/,&
!     &       10X,'Energy (Hartrees)    wfn I    = ',F20.12,/,&
!     &       10X,' " "   (Hartrees)    wfn J    = ',F20.12,/)
!3060  FORMAT(/,10X,'Nuclear Configuration: ',/,&
!     &       10X,'---------------------- ',//,&
!     &       10X,'Number of nuclei = ',I5,//,&
!     &       10X,' No. ',1X,5X,'X',11X,'Y',11X,'Z',/)
! 3070 FORMAT(10X,  I5   ,1X,3(F10.7,1X))
!!
!      END SUBROUTINE DENRD
!
!      SUBROUTINE DENGET(LUNIT,INSET,INKEY,FORM)
!!***********************************************************************
!!
!!     DENGET locates set number INSET on unit LUNIT with KEY = INKEY
!!
!!     If NSET = 0 file is positioned at end-of-information
!!             = 1 the file is opened
!!             = n file is positioned at the beginning of set number n
!!
!!     On return INSET = sequence number of current set
!!
!!***********************************************************************
!      IMPLICIT NONE
!!
!!..... Integer variables passed in the argument list
!!
!      INTEGER  LUNIT,INSET,INKEY
!!
!!..... Local integer variables
!!
!      INTEGER  KEY,NSET,NREC,I,MSET
!!
!!..... Character variables passed in the argument list
!!
!      character(len=11) FORM
!!
!!..... Local logical variables including the debug
!!
!      LOGICAL  OP
!!
!!---- Enquire about the status of the dataset on unit LUNIT,
!!
!      INQUIRE(UNIT=LUNIT,OPENED=OP)
!!
!!---- If the data set is not OPEN then we must open it and rewind.
!!
!      IF(.NOT.OP) OPEN(UNIT=LUNIT,ERR=99,FORM=FORM,STATUS='UNKNOWN')
!      REWIND(UNIT=LUNIT,ERR=100)
!!
!!---- If we were simply asked to open the dataset, INSET=1, then we
!!     may exit now.
!!
!      IF (INSET.EQ.1) GOTO 800
!!
!!---- Locate set number INSET
!!
! 5    CONTINUE
!!
!      IF (FORM.EQ.'FORMATTED') THEN
!        READ(LUNIT,*,END=9) KEY,NSET,NREC
!        IF (NSET.EQ.INSET .AND. KEY.EQ.INKEY) THEN
!          BACKSPACE LUNIT
!          GOTO 800
!        ELSE
!          DO 2 I=1,NREC
!          READ(LUNIT,*,END=199)
! 2        continue
!        END IF
!      ELSE
!        READ(LUNIT,END=9) KEY,NSET,NREC
!        IF (NSET.EQ.INSET .AND. KEY.EQ.INKEY) THEN
!          BACKSPACE LUNIT
!          GOTO 800
!        ELSE
!          DO 3 I=1,NREC
!          READ(LUNIT,END=199)
! 3        continue
!        END IF
!      END IF
!!
!      IF (NSET+1.EQ.INSET) GOTO 800
!!
!      MSET=NSET
!!
!      GOTO 5
!!
!!---- At END of file on read we branch here.
!!
! 9    CONTINUE
!!
!      IF (INSET.EQ.0) THEN
!        BACKSPACE LUNIT
!        INSET=MSET+1
!        GOTO 800
!      ELSE
!        GOTO 99
!      END IF
!!
!  800 CONTINUE
!      RETURN
!!
!!---- Error handler - set number NSET not found on unit LUNIT
!!
!   99 CONTINUE
!!
!      WRITE(6,9900)
!      WRITE(6,9920) LUNIT,INSET
!!
!      STOP
!!
!!---- Error on rewind of the file
!!
!  100 CONTINUE
!!
!      WRITE(6,9900)
!      WRITE(6,9930) LUNIT
!!
!      STOP
!!
!!---- End of file while reading contents of a dataset
!!
!  199 CONTINUE
!!
!      WRITE(6,9900)
!      WRITE(6,9940) LUNIT
!!
!!---- Format Statements
!!
! 9900 FORMAT(//,5X,'**** Error in DENGET ',//)
! 9920 FORMAT(5X,'File, unit =',I3,', does not contain dataset NSET =',&
!     &       I5,' or data header error')
! 9930 FORMAT(5X,'Rewind failed on unit = ',I3,/)
! 9940 FORMAT(5X,'Error on read in middle of a set of data ',i5,/)
!!
!      END SUBROUTINE DENGET

END MODULE rmat_rd
