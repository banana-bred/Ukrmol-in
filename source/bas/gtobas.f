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
!*==gtobas_data.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
      MODULE GTOBAS_DATA
      USE PRECISN, ONLY : WP
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER, PARAMETER :: NX=1000, MX=20, NMAX=5000

      INTEGER, PARAMETER :: IOUT=6 ! Logical unit for output printing
 
      REAL(KIND=wp), PARAMETER :: SMALL=1.E-4_wp

C Contents of common blocks:
      INTEGER, SAVE :: IQQQ, MA, NH, NL, NPT, NCOM
      REAL(KIND=wp), SAVE :: RADMAT
      REAL(KIND=wp), DIMENSION(nx), SAVE :: X
      REAL(KIND=wp), DIMENSION(nx,mx), SAVE :: YPO
      REAL(KIND=wp), DIMENSION(mx,mx), SAVE :: CH
      REAL(KIND=wp), DIMENSION(nmax), SAVE :: PCOM, XICOM

C      COMMON /coefff/ CH, NH
C      COMMON /DATA  / NPT, MA, MP
C      COMMON /LL    / NL
C      COMMON /POINTS/ X, YPO, RADMAT
C      COMMON /QUES2 / IQQQ
C      COMMON /F1COM / PCOM, XICOM, NCOM

C
C*** End of declarations rewritten by SPAG
C
      END MODULE GTOBAS_DATA

!*==gtobas.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
 
      PROGRAM GTOBAS
 
c#########################################################################
c
c     GTOBAS is a program for fitting Gaussian-type orbitals (GTOs) to
c     continuum functions. The Exponents of GTOs are optimized using the
c     method of Nestmann and Peyerimhoff [1]. The numerical routines are
c     from Numerical Recipes [2]. The numerical continuum functions, which
c     provide the input for the program, must be generated from an
c     external module (e.g. NUMBAS).
c
c     [1] Nestmann B. M. and Peyerimhoff S. D., J. Phys. B. (1990) L773
c     [2] Press W. H., Teukolski S. A., Vetterling W. T. and Flannery
c     B. P., Numerical Recipes in Fortran 77, Cambridge University Press
c     (1986-1992).
c
c#########################################################################
c
      USE GTOBAS_DATA, only:nx, mx, IOUT, SMALL, NH, NPT, MA, NL, 
     &                      X, YPO, RADMAT, IQQQ
      USE PRECISN, ONLY : WP
      USE GLOBAL_UTILS, ONLY : PRINT_UKRMOL_HEADER
      USE CONSTS, ONLY : ZERO=>XZERO, ONE=>XONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(mx) :: A, EXPO, FIX, W
      REAL(KIND=wp) :: BETA=0.016_wp, CHISQ, ERROR=99999.0_wp, EXCH, 
     &                FINALFIT, FRET, FTOL=1.E-9_wp, GAMMA=1.39_wp, 
     &                RAND, RDLOW=0.01_wp, RDUP=0.49_wp, XP
      INTEGER :: I, IGUESS=2, II, IPRINT=0, ISWMOL3=18, ITER, J, 
     &           JJ, K, LUINP=5, LUNUMB=13, LUPLOT=17, LVAL=0, NNPLOT, 
     &           NOEXP=0, NPLOT=0
      REAL(KIND=wp), DIMENSION(nx) :: SIG, XSTOR, Y
      REAL(KIND=wp), DIMENSION(nx,mx) :: U, YPOSTOR
      REAL(KIND=wp), DIMENSION(mx,mx) :: V, XI
      EXTERNAL FUNCS
C
C*** End of declarations rewritten by SPAG
C
c#########################################################################
c
c     The variables input via namelist /fit/ are
c     (name   :: description (default value)
c
c     iguess  :: initial selection of exponents (2),
c                =0 exponents provided by the user,
c                =1 random selection,
c                =2 selection from an appropriate function
c     expo    :: list of initial exponents (needed if iguess=0)
c     noexp   :: number of exponents to be used (needed if iguess=1,2).
c                Must be larger or equal to the number of functions.
c     nplot   :: no functions (0) / all functions (1) on unit luplot (0)
c     lunumb  :: unit number for numerical functions input (13 for NUMBAS)
c     luplot  :: unit number for plotable functions output (if used) (17)
c     iswmol3 :: unit number for final exponents in SWMOL3 format (18)
c     ftol    :: convergence criterion usid in POWELL (1.d-9)
c     iprint  :: print flag (0)
c                =1 all iteration data,
c                =2 plus mesh.
c     rdlow   :: lower limit of the random selection (0.01)
c                (used if iguess=1)
c     rdup    :: upper limit of the random selection (0.49)
c                (used if iguess=1)
c     beta    :: even-tempered beta coefficient (0.016)
c                (used if iguess=2)
c     gamma   :: even-tempered gamma coefficient (1.39)
c                (used if iguess=2)
c
      NAMELIST /fit/ expo, iprint, noexp, ftol, luplot, nplot, 
     &   iguess, lunumb, rdlow, rdup, beta, gamma, iswmol3
c
c#########################################################################
c
      CALL PRINT_UKRMOL_HEADER(iout)
      WRITE(iout,10)
 10   FORMAT(/////11X,' Program GTOBAS',//)
c
      DO i=1, mx
         fix(i)=one
         expo(i)=error
      END DO
c
c *** READ IN DATA & SET-UP PARAMETERS ***
c
      READ(luinp,fit)
c
      IF(iguess.EQ.0)THEN
         DO i=1, mx
            IF(expo(i).LT.error)noexp=i
         END DO
      ELSE IF(iguess.EQ.1)THEN
         DO i=1, noexp
            CALL random_number(rand)
            expo(i)=rdlow+rand*rdup
         END DO
      ELSE IF(iguess.EQ.2)THEN
         DO i=1, noexp
            expo(i)=beta*gamma**i
         END DO
      END IF
c
c *** Read functions generated by NUMBAS ***
c
      CALL READNUM(lunumb,iprint-1,radmat,lval,nh,npt,xstor,ypostor,nx)
c
      nl=lval
      ma=noexp
c
c#########################################################################
c
c     npt is number of grid points
c     nx is the maximum size of the grid
c     nl is the angular momentum
c     nh is number of functions to be fitted
c     ma is number of exponents to be optimised
c
c#########################################################################
c
      IF(npt.GE.nx .OR. nh.GT.ma)THEN
         WRITE(iout,20)
 20      FORMAT(/'*** ERROR: grid too large or too few exponents ***'/)
         STOP
      END IF
c
c
      WRITE(iout,101)radmat
 101  FORMAT(/' Boundary radius                      = ',f4.1)
      WRITE(iout,102)npt
 102  FORMAT(' Number of radial mesh points         = ',i4)
      WRITE(iout,103)ma
 103  FORMAT(' Number of exponents to be optimized  = ',i4)
      WRITE(iout,104)
 104  FORMAT('--------------------------------------------')
c
      iqqq=0
c
      WRITE(iout,2)
 2    FORMAT(/' Initial exponents ')
      DO i=1, ma
         WRITE(iout,3)i, expo(i)
 3       FORMAT(i3,d14.6)
      END DO
c
c *** Store mesh in reverse odrder ***
c
      DO i=1, npt
         x(i)=xstor(npt+1-i)
      END DO
c
c *** Store functions in reverse order ***
c
      DO j=1, nh
         DO i=1, npt
            ypo(i,j)=ypostor(npt+1-i,j)
         END DO
      END DO
c
c *** SET LENGTH OF INITIAL 'SEARCH' VECTOR ***
c
      DO ii=1, ma
         DO jj=1, ma
            IF(ii.EQ.jj)THEN
               xi(ii,jj)=fix(ii)
            ELSE
               xi(ii,jj)=zero
            END IF
         END DO
      END DO
c
c *** CALL POWELL ROUTINES TO MINIMIZE FUNCTION FUNCP ***
c
      CALL POWELL(expo,xi,ma,ftol,iter,mx,iprint,fret)
c
      WRITE(iout,5)iter
 5    FORMAT(/' Total Number of iterations =',i3)
      WRITE(iout,55)fret
 55   FORMAT(/' Final minimisation function =',d14.6)
      WRITE(iout,6)
c
c *** Final exponents in decreasing order ***
c
 6    FORMAT(/' Final exponents ')
      DO i=1, ma-1
         j=i
         DO k=i+1, ma
            IF(expo(k).GT.expo(j))THEN
               j=k
            END IF
         END DO
         exch=expo(i)
         expo(i)=expo(j)
         expo(j)=exch
      END DO
      DO i=1, ma
         WRITE(iout,7)i, expo(i)
 7       FORMAT(i3,d14.6)
      END DO
c
      IF(nplot.EQ.1)THEN
c
         WRITE(iout,8)nh, luplot
 8       FORMAT(/' Writing ',i3,' functions to unit',i3,' for plotting')
         DO nnplot=1, nh
            DO i=1, npt
               sig(i)=small
               y(i)=ypo(i,nnplot)
            END DO
            DO i=1, ma
               a(i)=one
            END DO
            iqqq=1
c
            CALL SVDFIT(x,y,sig,npt,a,ma,u,v,w,nx,chisq,funcs,nl,expo)
c
c *** write functions & FITS TO luplot ***
c
            WRITE(luplot,9)nnplot
 9          FORMAT(' # Function ',i3)
            DO k=1, npt
               xp=x(k)
               finalfit=zero
               DO i=1, ma
                  finalfit=finalfit+a(i)*xp**nl*exp(-expo(i)*xp**2)
               END DO
               WRITE(luplot,97)x(k), finalfit, ypo(k,nnplot)
 97            FORMAT(f9.4,d14.6,d14.6)
            END DO
c
         END DO
      END IF
c
c *** Save exponents in SWMOL3 format (decreasing order) ***
c
      IF(iswmol3.GT.0)THEN
         DO i=1, ma-1
            j=i
            DO k=i+1, ma
               IF(expo(k).GT.expo(j))THEN
                  j=k
               END IF
            END DO
            exch=expo(i)
            expo(i)=expo(j)
            expo(j)=exch
         END DO
         WRITE(iout,98)iswmol3
 98      FORMAT(/' Saving exponents in SWMOL3 format to unit',i3/)
         DO i=1, ma
            WRITE(iswmol3,99)1, 1
            WRITE(iswmol3,999)expo(i), 1._wp
 99         FORMAT(i5,i5)
 999        FORMAT(f9.6,f6.0)
         END DO
      END IF
c
      STOP
      END PROGRAM GTOBAS
!*==funcs.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
c ***  FUNCS defines the Gaussian-type functions ***
c
      SUBROUTINE FUNCS(x,na,nl,ax,afunc)
      USE PRECISN, ONLY : WP                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: NA, NL
      REAL(KIND=wp) :: X
      REAL(KIND=wp), DIMENSION(na) :: AFUNC, AX
      INTENT (IN) AX, NA, NL, X
      INTENT (OUT) AFUNC
C
C Local variables
C
      INTEGER :: I
C
C*** End of declarations rewritten by SPAG
C
c
      DO i=1, na
         afunc(i)=x**nl*exp(-ax(i)*x**2)
      END DO
c
      RETURN
      END SUBROUTINE FUNCS
!*==funcp.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
c ***  FUNCTION TO MINIMIZE: see Nestmann et al, J Phys B, 23 (1990) L773 ***
c
      FUNCTION FUNCP(ax)
      USE GTOBAS_DATA, only:nx, mx, SMALL, CH, NH, NPT, MA, NL,
     &                      X, YPO, RADMAT, IQQQ
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, ONE=>XONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp), DIMENSION(mx) :: AX
      REAL(KIND=wp) :: FUNCP
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(mx) :: A, W
      REAL(KIND=wp) :: AAA, BBB, CHISQ, DDD, EEE, FFF, RR, SJ, SM1, SM2
      INTEGER :: H, I, J, K
      REAL(KIND=wp), DIMENSION(nx) :: SIG, Y
      REAL(KIND=wp), DIMENSION(nx,mx) :: U
      REAL(KIND=wp), DIMENSION(mx,mx) :: V
      EXTERNAL FUNCS
C
C*** End of declarations rewritten by SPAG
C
c
c
c
      fff=zero
      DO i=1, npt
         sig(i)=small
      END DO
c
c *** FIND COEFFICIENTS ***
c
      DO h=1, nh
         DO i=1, npt
            y(i)=ypo(i,h)
         END DO
         j=0
         IF(iqqq.EQ.0)THEN
            DO i=1, ma
               a(i)=one
               j=j+1
            END DO
            iqqq=1
         ELSE
            DO i=1, ma
               j=j+1
               a(i)=ch(h,i)
            END DO
         END IF
c
         CALL SVDFIT(x,y,sig,npt,a,ma,u,v,w,nx,chisq,FUNCS,nl,ax)
c
         DO i=1, ma
            ch(h,i)=a(i)
         END DO
         fff=fff+chisq
      END DO
c
      fff=zero
c
c *** TERM TO AVOID CONVERGENCE OF TWO EXPONENTS TO SAME VALUE ***
c
      eee=zero
      DO i=2, ma
         DO j=1, i-1
            eee=eee+exp(-radmat*abs(ax(i)/ax(j)-ax(j)/ax(i)))
         END DO
      END DO
c
      DO h=1, nh
         ddd=zero
         sm1=zero
         sm2=zero
c
c *** FORWARD fit (Standard) ***
c
         DO k=1, npt
c
c *** BACKWARD fit (worth trying sometimes) ***
c
c         do k=npt,1,-1
            rr=x(k)
            aaa=zero
            sj=ypo(k,h)
            DO i=1, ma
               bbb=ch(h,i)*rr**nl*exp(-ax(i)*rr**2)
               aaa=aaa+bbb
            END DO
            aaa=aaa-sj
            aaa=aaa*aaa
c
            ddd=sj*sj
            sm1=sm1+aaa
            sm2=sm2+ddd
         END DO
         fff=fff+(sm1/sm2)
      END DO
      funcp=fff+eee
C
      RETURN
      END FUNCTION FUNCP
!*==svdfit.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
 
c######################################################################
c
c     SVDFIT, SVBKSB and SVDCMP are Numerical Recipes subroutines used
c     to perform a least-square fit by use of a singular value
c     decomposition technique (ref [2] pp. 51 and 665).
c
c######################################################################
c
      SUBROUTINE SVDFIT(x,y,sig,ndata,a,ma,u,v,w,mp,chisq,FUNCS,nl,ax)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, ONE=>XONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: CHISQ
      INTEGER :: MA, MP, NDATA, NL
      REAL(KIND=wp), DIMENSION(ma) :: A, AX, W
      REAL(KIND=wp), DIMENSION(ndata) :: SIG, X, Y
      REAL(KIND=wp), DIMENSION(mp,ma) :: U
      REAL(KIND=wp), DIMENSION(ma,ma) :: V
      INTENT (IN) SIG, Y
      INTENT (INOUT) CHISQ, W
      EXTERNAL FUNCS
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(ma) :: AFUNC
      REAL(KIND=wp), DIMENSION(ndata) :: B
      INTEGER :: I, J
      REAL(KIND=wp) :: SUM, THRESH, TMP, TOL=1.E-15_wp, WMAX
C
C*** End of declarations rewritten by SPAG
C
c
      DO i=1, ndata
         CALL FUNCS(x(i),ma,nl,ax,afunc)
         tmp=one/sig(i)
         DO j=1, ma
            u(i,j)=afunc(j)*tmp
         END DO
         b(i)=y(i)*tmp
      END DO
c
      CALL SVDCMP(u,ndata,ma,mp,w,v)
c
      wmax=zero
      DO j=1, ma
         IF(w(j).GT.wmax)wmax=w(j)
      END DO
      thresh=tol*wmax
      DO j=1, ma
         IF(w(j).LT.thresh)w(j)=zero
      END DO
c
      CALL SVBKSB(u,w,v,ndata,ma,mp,b,a)
c
      chisq=zero
      DO i=1, ndata
         CALL FUNCS(x(i),ma,nl,ax,afunc)
         sum=zero
         DO j=1, ma
            sum=sum+a(j)*afunc(j)
         END DO
         chisq=chisq+((y(i)-sum)/sig(i))**2
      END DO
c
      RETURN
      END SUBROUTINE SVDFIT
!*==svbksb.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
      SUBROUTINE SVBKSB(u,w,v,m,n,mp,b,x)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: M, MP, N
      REAL(KIND=wp), DIMENSION(mp) :: B
      REAL(KIND=wp), DIMENSION(mp,n) :: U
      REAL(KIND=wp), DIMENSION(n,n) :: V
      REAL(KIND=wp), DIMENSION(n) :: W, X
      INTENT (IN) B, M, MP, N, U, V, W
      INTENT (OUT) X
C
C Local variables
C
      INTEGER :: I, J, JJ
      REAL(KIND=wp) :: S
      REAL(KIND=wp), DIMENSION(n) :: TMP
C
C*** End of declarations rewritten by SPAG
C
c
      DO j=1, n
         s=zero
         IF(w(j).NE.zero)THEN  ! jmc perhaps better for this comparison to be with tiny(1.d0) (and rephrased...)
            DO i=1, m
               s=s+u(i,j)*b(i)
            END DO
            s=s/w(j)
         END IF
         tmp(j)=s
      END DO
      DO j=1, n
         s=zero
         DO jj=1, n
            s=s+v(j,jj)*tmp(jj)
         END DO
         x(j)=s
      END DO
c
      RETURN
      END SUBROUTINE SVBKSB
!*==svdcmp.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
      SUBROUTINE SVDCMP(a,m,n,mp,w,v)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, ONE=>XONE, TWO=>XTWO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: M, MP, N
      REAL(KIND=wp), DIMENSION(mp,n) :: A
      REAL(KIND=wp), DIMENSION(n,n) :: V
      REAL(KIND=wp), DIMENSION(n) :: W
      INTENT (IN) M, MP, N
      INTENT (INOUT) A, V, W
C
C Local variables
C
      REAL(KIND=wp) :: ANORM, C, F, G, H, S, SCALE, X, Y, Z
      INTEGER :: I, ITS, J, JJ, K, L, NM
      REAL(KIND=wp), EXTERNAL :: PYTHAG
      REAL(KIND=wp), DIMENSION(n) :: RV1
C
C*** End of declarations rewritten by SPAG
C
c
      g=zero
      scale=zero
      anorm=zero
      DO i=1, n
         l=i+1
         rv1(i)=scale*g
         g=zero
         s=zero
         scale=zero
         IF(i.LE.m)THEN
            DO k=i, m
               scale=scale+abs(a(k,i))
            END DO
            IF(scale.NE.zero)THEN
               DO k=i, m
                  a(k,i)=a(k,i)/scale
                  s=s+a(k,i)*a(k,i)
               END DO
               f=a(i,i)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,i)=f-g
               DO j=l, n
                  s=zero
                  DO k=i, m
                     s=s+a(k,i)*a(k,j)
                  END DO
                  f=s/h
                  DO k=i, m
                     a(k,j)=a(k,j)+f*a(k,i)
                  END DO
               END DO
               DO k=i, m
                  a(k,i)=scale*a(k,i)
               END DO
            END IF
         END IF
         w(i)=scale*g
         g=zero
         s=zero
         scale=zero
         IF((i.LE.m) .AND. (i.NE.n))THEN
            DO k=l, n
               scale=scale+abs(a(i,k))
            END DO
            IF(scale.NE.zero)THEN
               DO k=l, n
                  a(i,k)=a(i,k)/scale
                  s=s+a(i,k)*a(i,k)
               END DO
               f=a(i,l)
               g=-sign(sqrt(s),f)
               h=f*g-s
               a(i,l)=f-g
               DO k=l, n
                  rv1(k)=a(i,k)/h
               END DO
               DO j=l, m
                  s=zero
                  DO k=l, n
                     s=s+a(j,k)*a(i,k)
                  END DO
                  DO k=l, n
                     a(j,k)=a(j,k)+s*rv1(k)
                  END DO
               END DO
               DO k=l, n
                  a(i,k)=scale*a(i,k)
               END DO
            END IF
         END IF
         anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      END DO
      DO i=n, 1, -1
         IF(i.LT.n)THEN
            IF(g.NE.zero)THEN
               DO j=l, n
                  v(j,i)=(a(i,j)/a(i,l))/g
               END DO
               DO j=l, n
                  s=zero
                  DO k=l, n
                     s=s+a(i,k)*v(k,j)
                  END DO
                  DO k=l, n
                     v(k,j)=v(k,j)+s*v(k,i)
                  END DO
               END DO
            END IF
            DO j=l, n
               v(i,j)=zero
               v(j,i)=zero
            END DO
         END IF
         v(i,i)=one
         g=rv1(i)
         l=i
      END DO
      DO i=min(m,n), 1, -1
         l=i+1
         g=w(i)
         DO j=l, n
            a(i,j)=zero
         END DO
         IF(g.NE.zero)THEN
            g=one/g
            DO j=l, n
               s=zero
               DO k=l, m
                  s=s+a(k,i)*a(k,j)
               END DO
               f=(s/a(i,i))*g
               DO k=i, m
                  a(k,j)=a(k,j)+f*a(k,i)
               END DO
            END DO
            DO j=i, m
               a(j,i)=a(j,i)*g
            END DO
         ELSE
            DO j=i, m
               a(j,i)=zero
            END DO
         END IF
         a(i,i)=a(i,i)+one
      END DO
      DO k=n, 1, -1
         DO its=1, 30
            DO l=k, 1, -1
               nm=l-1
               IF((abs(rv1(l))+anorm).EQ.anorm)GO TO 2
               IF((abs(w(nm))+anorm).EQ.anorm)EXIT
            END DO
            c=zero
            s=one
            DO i=l, k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               IF((abs(f)+anorm).EQ.anorm)EXIT
               g=w(i)
               h=pythag(f,g)
               w(i)=h
               h=one/h
               c=(g*h)
               s=-(f*h)
               DO j=1, m
                  y=a(j,nm)
                  z=a(j,i)
                  a(j,nm)=(y*c)+(z*s)
                  a(j,i)=-(y*s)+(z*c)
               END DO
            END DO
 2          z=w(k)
            IF(l.EQ.k)THEN
               IF(z.LT.zero)THEN
                  w(k)=-z
                  DO j=1, n
                     v(j,k)=-v(j,k)
                  END DO
               END IF
               GO TO 3
            END IF
            IF(its.EQ.30)STOP 'no convergence in svdcmp'
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(two*h*y)
            g=PYTHAG(f,one)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=one
            s=one
            DO j=l, nm
               i=j+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=PYTHAG(f,h)
               rv1(j)=z
               c=f/z
               s=h/z
               f=(x*c)+(g*s)
               g=-(x*s)+(g*c)
               h=y*s
               y=y*c
               DO jj=1, n
                  x=v(jj,j)
                  z=v(jj,i)
                  v(jj,j)=(x*c)+(z*s)
                  v(jj,i)=-(x*s)+(z*c)
               END DO
               z=pythag(f,h)
               w(j)=z
               IF(z.NE.zero)THEN
                  z=one/z
                  c=f*z
                  s=h*z
               END IF
               f=(c*g)+(s*y)
               x=-(s*g)+(c*y)
               DO jj=1, m
                  y=a(jj,j)
                  z=a(jj,i)
                  a(jj,j)=(y*c)+(z*s)
                  a(jj,i)=-(y*s)+(z*c)
               END DO
            END DO
            rv1(l)=zero
            rv1(k)=f
            w(k)=x
         END DO
 3       CONTINUE
      END DO
      RETURN
      END SUBROUTINE SVDCMP
!*==pythag.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
c *** Used by svdcmp ***
c
      FUNCTION PYTHAG(a,b)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, ONE=>XONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: A, B
      REAL(KIND=wp) :: PYTHAG
      INTENT (IN) A, B
C
C Local variables
C
      REAL(KIND=wp) :: ABSA, ABSB
C
C*** End of declarations rewritten by SPAG
C
      absa=abs(a)
      absb=abs(b)
      IF(absa.GT.absb)THEN
         pythag=absa*sqrt(one+(absb/absa)**2)
      ELSE
         IF(absb.EQ.zero)THEN
            pythag=zero
         ELSE
            pythag=absb*sqrt(one+(absa/absb)**2)
         END IF
      END IF
      RETURN
      END FUNCTION PYTHAG
!*==linmin.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
c######################################################################
c
c     POWELL, LINMIN and MNBRAK are Numerical Recipes subroutines used
c     to minimize a multi-dimensional function. The procedure is based
c     on the Powell's method which is the prototype of multidimensional
c     direction-set methods (ref. [2] p. 406).
c
c     LINMIN implements the one-dimensional line minimization used by
c     POWELL. MNBRAK is used by LINMIN to bracket the minimum.
c
c####################################################################
c
c     Line Minimization
c
      SUBROUTINE LINMIN(p,xi,n,fret)
      USE GTOBAS_DATA, only: NCOM, PCOM, XICOM
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, ONE=>XONE
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: FRET
      INTEGER :: N
      REAL(KIND=wp), DIMENSION(n) :: P, XI
      INTENT (IN) N
      INTENT (OUT) FRET
      INTENT (INOUT) P, XI
C
C Local variables
C
      REAL(KIND=wp) :: AX, BX, FA, FB, FX, TOL=1.E-5_wp, XMIN, XX
      REAL(KIND=wp), EXTERNAL :: BRENT, F1DIM
      INTEGER :: J
C
C*** End of declarations rewritten by SPAG
C
c
      ncom=n
      DO j=1, n
         pcom(j)=p(j)
         xicom(j)=xi(j)
      END DO
      ax=zero
      xx=one
      CALL mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,tol,xmin)
      DO j=1, n
         xi(j)=xmin*xi(j)
         p(j)=p(j)+xi(j)
      END DO
c
      RETURN
      END SUBROUTINE LINMIN
!*==f1dim.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
      FUNCTION F1DIM(x)
      USE GTOBAS_DATA, only: NCOM, PCOM, XICOM
      USE PRECISN, ONLY : WP                        
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: X
      REAL(KIND=wp) :: F1DIM
      INTENT (IN) X
C
C Local variables
C
      REAL(KIND=wp), DIMENSION(ncom) :: AX, XT
      REAL(KIND=wp), EXTERNAL :: FUNCP
      INTEGER :: J
C
C*** End of declarations rewritten by SPAG
C
c
c
      DO j=1, ncom
         xt(j)=pcom(j)+x*xicom(j)
         ax(j)=exp(xt(j))
      END DO
c
      f1dim=FUNCP(ax)
c
      RETURN
c
      END FUNCTION F1DIM
!*==brent.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
      FUNCTION BRENT(ax,bx,cx,f,tol,xmin)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, HALF=>XHALF, TWO=>XTWO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: AX, BX, CX, TOL, XMIN
      REAL(KIND=wp) :: BRENT
      REAL(KIND=wp), EXTERNAL :: F
      INTENT (IN) AX, BX, CX, TOL
      INTENT (OUT) XMIN
C
C Local variables
C
      REAL(KIND=wp) :: A, B, CGOLD=0.3819660_wp, D, E, ETEMP, FU, FV, 
     &                FW, FX, P, Q, R, TOL1, TOL2, U, 
     &                V, W, X, XM, ZEPS=1.0E-15_wp
      INTEGER :: ITER, ITMAX=1000
C
C*** End of declarations rewritten by SPAG
C
c
c
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=zero
      fx=f(x)
      fv=fx
      fw=fx
      DO iter=1, itmax
         xm=half*(a+b)
         tol1=tol*abs(x)+zeps
         tol2=two*tol1
         IF(abs(x-xm).LE.(tol2-half*(b-a)))GO TO 3
         IF(abs(e).GT.tol1)THEN
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=two*(q-r)
            IF(q.GT.zero)p=-p
            q=abs(q)
            etemp=e
            e=d
            IF(abs(p).GE.abs(half*q*etemp) .OR. p.LE.q*(a-x) .OR. 
     &         p.GE.q*(b-x))GO TO 1
            d=p/q
            u=x+d
            IF(u-a.LT.tol2 .OR. b-u.LT.tol2)d=sign(tol1,xm-x)
            GO TO 2
         END IF
 1       IF(x.GE.xm)THEN
            e=a-x
         ELSE
            e=b-x
         END IF
         d=cgold*e
 2       IF(abs(d).GE.tol1)THEN
            u=x+d
         ELSE
            u=x+sign(tol1,d)
         END IF
         fu=f(u)
         IF(fu.LE.fx)THEN
            IF(u.GE.x)THEN
               a=x
            ELSE
               b=x
            END IF
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
         ELSE
            IF(u.LT.x)THEN
               a=u
            ELSE
               b=u
            END IF
            IF(fu.LE.fw .OR. w.EQ.x)THEN
               v=w
               fv=fw
               w=u
               fw=fu
            ELSE IF(fu.LE.fv .OR. v.EQ.x .OR. v.EQ.w)THEN
               v=u
               fv=fu
            END IF
         END IF
      END DO
      STOP 'brent exceed max. iterations'
 3    xmin=x
      brent=fx
      RETURN
c
      END FUNCTION BRENT
!*==mnbrak.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
c *** BRACKET MINIMUM ***
c
      SUBROUTINE MNBRAK(ax,bx,cx,fa,fb,fc,func)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, TWO=>XTWO
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: AX, BX, CX, FA, FB, FC
      REAL(KIND=wp), EXTERNAL :: FUNC
      INTENT (INOUT) AX, BX, CX, FA, FB, FC
C
C Local variables
C
      REAL(KIND=wp) :: DUM, FU, GLIMIT=1000.0_wp, GOLD=1.618034_wp, Q,
     &                R, TINY=1.E-20_wp, U, ULIM
C
C*** End of declarations rewritten by SPAG
C
c
      fa=func(ax)
      fb=func(bx)
      IF(fb.GT.fa)THEN
         dum=ax
         ax=bx
         bx=dum
         dum=fb
         fb=fa
         fa=dum
      END IF
      cx=bx+gold*(bx-ax)
      fc=func(cx)
 1    IF(fb.GE.fc)THEN
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(two*sign(max(abs(q-r),tiny),q-r))
         ulim=bx+glimit*(cx-bx)
         IF((bx-u)*(u-cx).GT.ZERO)THEN
            fu=func(u)
            IF(fu.LT.fc)THEN
               ax=bx
               fa=fb
               bx=u
               fb=fu
               RETURN
            ELSE IF(fu.GT.fb)THEN
               cx=u
               fc=fu
               RETURN
            END IF
            u=cx+gold*(cx-bx)
            fu=func(u)
         ELSE IF((cx-u)*(u-ulim).GT.ZERO)THEN
            fu=func(u)
            IF(fu.LT.fc)THEN
               bx=cx
               cx=u
               u=cx+gold*(cx-bx)
               fb=fc
               fc=fu
               fu=func(u)
            END IF
         ELSE IF((u-ulim)*(ulim-cx).GE.ZERO)THEN
            u=ulim
            fu=func(u)
         ELSE
            u=cx+gold*(cx-bx)
            fu=func(u)
         END IF
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         GO TO 1
      END IF
c
      RETURN
      END SUBROUTINE MNBRAK
!*==powell.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
c
c *** POWELL ***
c
      SUBROUTINE POWELL(pps,xi,n,ftol,iter,mx,iprint,fret)
      USE PRECISN, ONLY : WP                        
      USE CONSTS, ONLY : ZERO=>XZERO, TWO=>XTWO
      USE GTOBAS_DATA, ONLY : OUT=>IOUT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL(KIND=wp) :: FRET, FTOL
      INTEGER :: IPRINT, ITER, MX, N
      REAL(KIND=wp), DIMENSION(n) :: PPS
      REAL(KIND=wp), DIMENSION(mx,mx) :: XI
      INTENT (IN) FTOL, IPRINT, MX
      INTENT (INOUT) FRET, ITER, PPS, XI
C
C Local variables
C
      REAL(KIND=wp) :: DEL, FP, FPTT, FUNCMINN=1.0E2_wp, T 
      REAL(KIND=wp), EXTERNAL :: FUNCP
      INTEGER :: I, IBIG, III, ITMAX=2000, J
      REAL(KIND=wp), DIMENSION(n) :: P, PT, PTT, XIT
C
C*** End of declarations rewritten by SPAG
C
c
c
      DO i=1, n
         p(i)=log(pps(i))
         pt(i)=p(i)
      END DO
c
      fret=FUNCP(pps)
c
      iter=0
 1    iter=iter+1
c
      IF(iprint.GT.0)THEN
         WRITE(out,7)iter, fret
 7       FORMAT(/' Starting iteration : ',i3,/' Func = ',d14.6)
         DO j=1, n
            WRITE(out,77)j, exp(p(j))
 77         FORMAT(i3,d14.6)
         END DO
      END IF
c
      IF(fret.LT.funcminn)funcminn=fret
c
      fp=fret
      ibig=0
      del=zero
      DO i=1, n
         DO j=1, n
            xit(j)=xi(j,i)
         END DO
         fptt=fret
         CALL LINMIN(p,xit,n,fret)
         IF(abs(fptt-fret).GT.del)THEN
            del=abs(fptt-fret)
            ibig=i
         END IF
      END DO
      IF(two*abs(fp-fret).LE.ftol*(abs(fp)+abs(fret)))RETURN
      IF(iter.EQ.itmax)THEN
         WRITE(out,*)'powell exceeding max. iterations'
         DO iii=1, n
            WRITE(out,*)p(iii)
         END DO
         STOP
      END IF
      DO j=1, n
         ptt(j)=two*p(j)-pt(j)
         xit(j)=p(j)-pt(j)
         pt(j)=p(j)
         pps(j)=exp(ptt(j))
      END DO
      fptt=FUNCP(pps)
      IF(fptt.GE.fp)GO TO 1
c
      t=two*(fp-two*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
c
      IF(t.GE.zero)GO TO 1
c
      CALL LINMIN(p,xit,n,fret)
      DO j=1, n
         xi(j,ibig)=xi(j,n)
         xi(j,n)=xit(j)
      END DO
c
      GO TO 1
c
      END SUBROUTINE POWELL
!*==readnum.spg  processed by SPAG 6.56Rc at 10:33 on 11 Mar 2010
      SUBROUTINE READNUM(NFTC,IFLOUT,ra,lval,inh,npt,rr,fun,nx)
C
C***********************************************************************
C
c     READNUM reads numerical continuum basis produced by NUMBAS
C
C***********************************************************************
C
      USE PRECISN, ONLY : WP                        
      USE GLOBAL_UTILS, ONLY : CWBOPN
      USE GTOBAS_DATA, ONLY : IWRITE=>IOUT
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER :: IFLOUT, INH, LVAL, NFTC, NPT, NX
      REAL(KIND=wp) :: RA
      REAL(KIND=wp), DIMENSION(nx,*) :: FUN
      REAL(KIND=wp), DIMENSION(nx) :: RR
      INTENT (IN) IFLOUT, NX
      INTENT (OUT) FUN, RA
      INTENT (INOUT) INH, LVAL, NPT, RR
C
C Local variables
C
      REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: AK, RX
C JMC      CHARACTER(LEN=8), DIMENSION(15) :: HEAD
C JMC changing the line above to the next line of code below to avoid having hard-wired 
C dimensions.  I'm assuming that the dim. 15 was because 15*8=120, 
C where 120 is the length of the title string written to the output unit of numcbas.f 
C (see subroutine WRHEAD in numcbas.f).  I have left the original line above in case 
C this assumption is incorrect!
      CHARACTER(LEN=120) :: HEAD
      INTEGER :: I, IR, IR1, IR2, J, NCF, NINTD, NINTP
      REAL(KIND=wp) :: R0
C
C*** End of declarations rewritten by SPAG
C
C
C     READ HEADER RECORDS AND PRINT
C
      CALL CWBOPN(nftc)
      READ(NFTC)HEAD
      WRITE(iwrite,2)
 2    FORMAT(///' Numerical basis read by READNUM'/)
      READ(NFTC)NCF, lval, npt
      WRITE(iwrite,3)Lval
 3    FORMAT(' Angular momentum LVAL =',i2)
C
      IF(npt.GT.nx)THEN
         WRITE(iwrite,4)npt, nx
 4       FORMAT(/' *** Too many mesh points npt=',i4,5x,'nx=',
     &          i4/' *** Increase nx in main routine')
         STOP
      END IF
C
C *** STORAGE ALLOCATION ***
C
      nintp=npt+1
      NINTD=NINTP*ncf
      ALLOCATE(rx(nintd),ak(ncf))
C
C *** READ mesh ***
C
      READ(NFTC)r0, (rr(IR),IR=1,npt)
      ra=rr(npt)
C
C *** READ functions ***
C
      IR2=0
      DO I=1, NCF
         IR1=IR2+1
         IR2=IR2+nintp
         READ(NFTC)(RX(IR),IR=IR1,IR2)
      END DO
c
      WRITE(iwrite,1111)
 1111 FORMAT(/' Selected poles '//'      Seqno  Lval  Energy')
      READ(NFTC)(AK(I),I=1,ncf)
      inh=0
      ir1=-npt
      DO J=1, ncf
         ir1=ir1+nintp
         inh=inh+1
         WRITE(iwrite,1112)inh, j, lval, ak(j)
 1112    FORMAT(3I5,f10.4)
         DO ir=1, npt
            fun(ir,inh)=rx(ir1+ir)
         END DO
      END DO
C
      IF(IFLOUT.GT.0)THEN
         WRITE(IWRITE,5454)
         DO J=1, npt
            WRITE(IWRITE,5353)J, rr(J)
         END DO
      END IF
      DEALLOCATE(rx,ak)
C
      RETURN
c
 5353 FORMAT(5X,I5,5X,F12.6)
 5454 FORMAT(//8X,'J',12X,'R'/)
c
      END SUBROUTINE READNUM
