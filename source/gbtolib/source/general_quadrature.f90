! Copyright 2019
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
module general_quadrature_gbl
use precisn_gbl
use utils_gbl, only: xermsg
use const_gbl, only: limit

 private

 public function_1d, function_2d, function_1d_many, n_7, x_7, w_7, n_10, x_10, w_10, gl_expand_A_B, get_gaussrule
 public quadrature_u_integral
 
 !DQAGS is the routine that performs the numerical quadrature
 public DQAGS

 public quad2d, gl2d, dqelg, quad1d

 !> \class <function_1d>
 !> This is a class that defines an abstract function of one variable, whose specific implementation is deferred.
 !> The purpose is to use this abstract function in some bspline-related routines which require a user-defined function as a parameters.
 type, abstract :: function_1d
    !> Number of function evaulation, number of subdivisions (applicable for adaptive quadratures).
    integer :: neval = 0, ndiv = 0
    !> Maximum number of sub-divisions of the area to integrate over.
    integer :: max_div = 2*limit
 contains
      !> This must be used in all routines using this object to evaluate the function. This symbol is resolved into one of wp_eval,
      !> ep_eval depending on the floating point type of X on input.
      generic, public :: eval => wp_eval, ep_eval
      !> \memberof function_1d
      procedure(wp_user_function_interface), deferred :: wp_eval
      !> \memberof function_1d
      procedure(ep_user_function_interface), deferred :: ep_eval
 end type function_1d
 abstract interface
      real(wp) function wp_user_function_interface(data,x)
         import :: function_1d, wp
         class(function_1d) :: data
         real(kind=wp), intent(in) :: x
      end function wp_user_function_interface
 end interface
 abstract interface
      real(ep1) function ep_user_function_interface(data,x)
         import :: function_1d, ep1
         class(function_1d) :: data
         real(kind=ep1), intent(in) :: x
      end function ep_user_function_interface
 end interface

 !> \class <function_1d_many>
 !> The same as function_1d but now the function evaluations are for n points at the same time.
 type, abstract, extends(function_1d) :: function_1d_many
   !no data components
 contains
      !> This must be used in all routines using this object to evaluate the function. This symbol is resolved into one of wp_eval, ep_eval depending on the floating point type of X on input.
      generic, public :: eval => wp_eval_many, ep_eval_many
      !> \memberof bound_user_function
      procedure(wp_user_function_interface_many), deferred :: wp_eval_many
      !> \memberof bound_user_function
      procedure(ep_user_function_interface_many), deferred :: ep_eval_many
 end type function_1d_many
 abstract interface
      function wp_user_function_interface_many(data,x,n)
         import :: function_1d_many, wp
         class(function_1d_many) :: data
         integer, intent(in) :: n
         real(kind=wp), intent(in) :: x(n)
         real(kind=wp) :: wp_user_function_interface_many(n)
      end function wp_user_function_interface_many
 end interface
 abstract interface
      function ep_user_function_interface_many(data,x,n)
         import :: function_1d_many, ep1
         class(function_1d_many) :: data
         integer, intent(in) :: n
         real(kind=ep1), intent(in) :: x(n)
         real(kind=ep1) :: ep_user_function_interface_many(n)
      end function ep_user_function_interface_many
 end interface

 !> \class <function_2d>
 !> This is a class that defines an abstract function of two variables, whose specific implementation is deferred.
 !> The purpose is to use this abstract function in some bspline-related routines which require a user-defined function as a parameters.
 type, abstract :: function_2d
    !> Number of function evaulation, number of subdivisions (applicable for adaptive quadratures).
    integer :: neval = 0, ndiv = 0
    !> Maximum number of sub-divisions of the area to integrate over: currently this is not used in quad2d.
    integer :: max_div = 2*limit
 contains
    procedure(fn2d_eval_interface), deferred :: eval
 end type function_2d
 abstract interface
      real(cfp) function fn2d_eval_interface(this,x,y)
         import :: function_2d, cfp
         class(function_2d) :: this
         real(kind=cfp), intent(in) :: x, y
      end function fn2d_eval_interface
 end interface

 !> Order of the Gauss-Legendre quadrature to which the x_7 and w_7 arrays correspond.
 integer, parameter :: n_7 = 7

 !> Weights for the Gauss-Legendre quadrature of order 7 on interval [0,1].
 real(kind=cfp), parameter :: w_7(2*n_7+1) = (/0.015376620998058634177314196788602209_cfp,&
                                           &0.035183023744054062354633708225333669_cfp,&
                                           &0.05357961023358596750593477334293465_cfp,&
                                           &0.06978533896307715722390239725551416_cfp,&
                                           &0.08313460290849696677660043024060441_cfp,&
                                           &0.09308050000778110551340028093321141_cfp,&
                                           &0.09921574266355578822805916322191966_cfp,&
                                           &0.10128912096278063644031009998375966_cfp,&
                                           &0.09921574266355578822805916322191966_cfp,&
                                           &0.09308050000778110551340028093321141_cfp,&
                                           &0.08313460290849696677660043024060441_cfp,&
                                           &0.06978533896307715722390239725551416_cfp,&
                                           &0.05357961023358596750593477334293465_cfp,&
                                           &0.035183023744054062354633708225333669_cfp,&
                                           &0.015376620998058634177314196788602209_cfp/)

 !> Abscissas for the Gauss-Legendre quadrature of order 7 on interval [0,1].
 real(kind=cfp), parameter :: x_7(2*n_7+1) = (/0.0060037409897572857552171407066937094_cfp,&
                                          &0.031363303799647047846120526144895264_cfp,&
                                          &0.075896708294786391899675839612891574_cfp,&
                                          &0.13779113431991497629190697269303100_cfp,&
                                          &0.21451391369573057623138663137304468_cfp,&
                                          &0.30292432646121831505139631450947727_cfp,&
                                          &0.39940295300128273884968584830270190_cfp,&
                                          &0.50000000000000000000000000000000000_cfp,&
                                          &0.60059704699871726115031415169729810_cfp,&
                                          &0.69707567353878168494860368549052273_cfp,&
                                          &0.78548608630426942376861336862695532_cfp,&
                                          &0.86220886568008502370809302730696900_cfp,&
                                          &0.92410329170521360810032416038710843_cfp,&
                                          &0.96863669620035295215387947385510474_cfp,&
                                          &0.99399625901024271424478285929330629_cfp/)

 !> Order of the Gauss-Legendre quadrature to which the x_10 and w_10 arrays correspond.
 integer, parameter :: n_10 = 10

 !> Abscissas for the Gauss-Legendre quadrature of order 10 on interval [0,1].
 real(kind=cfp), parameter :: x_10(2*n_10+1) = (/0.0031239146898052498698789820310295354_cfp,&
                                        &0.016386580716846852841688892546152419_cfp,&
                                        &0.039950332924799585604906433142515553_cfp,&
                                        &0.073318317708341358176374680706216165_cfp,&
                                        &0.11578001826216104569206107434688598_cfp,&
                                        &0.16643059790129384034701666500483042_cfp,&
                                        &0.22419058205639009647049060163784336_cfp,&
                                        &0.28782893989628060821316555572810597_cfp,&
                                        &0.35598934159879945169960374196769984_cfp,&
                                        &0.42721907291955245453148450883065683_cfp,&
                                        &0.50000000000000000000000000000000000_cfp,&
                                        &0.57278092708044754546851549116934317_cfp,&
                                        &0.64401065840120054830039625803230016_cfp,&
                                        &0.71217106010371939178683444427189403_cfp,&
                                        &0.77580941794360990352950939836215664_cfp,&
                                        &0.83356940209870615965298333499516958_cfp,&
                                        &0.88421998173783895430793892565311402_cfp,&
                                        &0.92668168229165864182362531929378384_cfp,&
                                        &0.96004966707520041439509356685748445_cfp,&
                                        &0.98361341928315314715831110745384758_cfp,&
                                        &0.99687608531019475013012101796897046_cfp/)

 !> Weights for the Gauss-Legendre quadrature of order 10 on interval [0,1].
 real(kind=cfp), parameter :: w_10(2*n_10+1) =  (/0.008008614128887166662112308429235508_cfp,&
                                           &0.018476894885426246899975334149664833_cfp,&
                                           &0.028567212713428604141817913236223979_cfp,&
                                           &0.038050056814189651008525826650091590_cfp,&
                                           &0.046722211728016930776644870556966044_cfp,&
                                           &0.05439864958357418883173728903505282_cfp,&
                                           &0.06091570802686426709768358856286680_cfp,&
                                           &0.06613446931666873089052628724838780_cfp,&
                                           &0.06994369739553657736106671193379156_cfp,&
                                           &0.07226220199498502953191358327687627_cfp,&
                                           &0.07304056682484521359599257384168559_cfp,&
                                           &0.07226220199498502953191358327687627_cfp,&
                                           &0.06994369739553657736106671193379156_cfp,&
                                           &0.06613446931666873089052628724838780_cfp,&
                                           &0.06091570802686426709768358856286680_cfp,&
                                           &0.05439864958357418883173728903505282_cfp,&
                                           &0.046722211728016930776644870556966044_cfp,&
                                           &0.038050056814189651008525826650091590_cfp,&
                                           &0.028567212713428604141817913236223979_cfp,&
                                           &0.018476894885426246899975334149664833_cfp,&
                                           &0.008008614128887166662112308429235508_cfp/)

contains

  !> Determines the n-point Gauss-Legendre quadrature rule for the interval [0,1]. Preferentially it uses the hard-coded rules for n = 15 and n = 21,
  !> otherwise it generates the rule using the routine gaussrule.
  subroutine get_gaussrule(n,x,w)
      implicit none
      integer, intent(in) ::  n
      real(kind=cfp), intent(out) :: x(n),w(n)

      real(kind=cfp), parameter :: quad_start = 0.0_cfp, quad_end = 1.0_cfp

        if (n .eq. 2*n_10+1) then
           x = x_10
           w = w_10
        elseif (n .eq. 2*n_7+1) then
           x = x_7
           w = w_7
        else
           call gaussrule(quad_start,quad_end,x,w,n)
        endif

  end subroutine get_gaussrule

  !> Determines the n-point Gauss-Legendre quadrature rule for the interval [r1,r2].
  subroutine gaussrule(r1,r2,x,w,n)
      use phys_const_gbl, only: pi
      implicit none
      integer, intent(in) ::  n
      real(kind=cfp), intent(in) :: r1,r2
      real(kind=cfp), intent(out) :: x(n),w(n)

      integer i,j,m
      real(kind=cfp) p1,p2,p3,pp,xl,xm,z,z1,eps

         eps = F1MACH(3,cfp_dummy)
         m=(n+1)/2
         xm=0.5_cfp*(r2+r1)
         xl=0.5_cfp*(r2-r1)
         do i=1,m
            z=cos(pi*(i-0.25_cfp)/(n+0.50_cfp))
            do
               p1=1.0_cfp
               p2=0.0_cfp
               do j=1,n
                  p3=p2
                  p2=p1
                  p1=((2.0_cfp*j-1.0_cfp)*z*p2-(j-1.0_cfp)*p3)/j
               enddo
               pp=n*(z*p1-p2)/(z*z-1.0_cfp)
               z1=z
               z=z1-p1/pp
               if (abs(z-z1).le.EPS) exit
            enddo
            x(i)=xm-xl*z
            x(n+1-i)=xm+xl*z
            w(i)=2.0_cfp*xl/((1.0_cfp-z*z)*pp*pp)
            w(n+1-i)=w(i)
         enddo !i

  end subroutine gaussrule

  !> Takes the Gauss-Legendre rule for the interval [0,1] and expands it for the given interval [A,B].
  subroutine gl_expand_A_B(x,w,n,x_AB,w_AB,A,B)
     implicit none
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: A, B
     real(kind=cfp), intent(in) :: x(2*n+1), w(2*n+1)
     real(kind=cfp), intent(out) :: x_AB(2*n+1), w_AB(2*n+1)

     integer :: i
     real(kind=cfp) :: delta

        delta = B-A
        do i=1,2*n+1
           x_AB(i) = x(i)*delta + A
           w_AB(i) = w(i)*delta
        enddo !i

  end subroutine gl_expand_A_B

 !> Adaptive 1D quadrature based on Gauss-Legendre rule.
 !> \param [in] f The 1D function to be integrated.
 !> \param [in] A Interval start.
 !> \param [in] B Interval end.
 !> \param [in] eps Required relative precision for the integral.
 !> \param [in] Qest Optional: an estimate of the integral over the interval as obtained by a call to gl1d. 
 recursive function quad1d(f,A,B,eps,Qest) result(I)
     implicit none
     class(function_1d) :: f
     real(kind=cfp), intent(in) :: A,B,eps
     real(kind=cfp), optional, intent(in) :: Qest

     real(kind=cfp) :: half, QP1, QP2, Q, I

        !Quit if too many sub-divisions have been performed: todo this should trigger an error/warning message!!
        !if (f%ndiv > f%max_div) return

        half = (A+B)*0.5_cfp

        !If the estimate of the quadrature on this rectangle using a previous call to gl2d is given then use it otherwise calculate it.
        if (present(Qest)) then
           Q = Qest
        else
           Q = gl1d(f,A,B)
        endif

        !Calculate quadratures on the four sub-rectangles
        QP1 = gl1d(f,A,half)
        QP2 = gl1d(f,half,B)

        I = QP1+QP2

        !Continue recursively on each rectangle if the desired precision has not been reached
        !todo instead of dividing all four try dividing the largest one, then the second, etc. and see if we converge without improving all at the same time.
        if (abs((Q-I)/I) > eps) then
           Q = quad1d(f,A,half,eps,QP1)
           I = Q
           Q = quad1d(f,half,B,eps,QP2)
           I = I + Q
        endif

        f%ndiv = f%ndiv + 1

 end function quad1d

 !> 1D Quadrature on rectangle using the Gauss-Legendre rule of order 7. The meaning of the input variables is identical to quad1d input parameters.
 function gl1d(f,A,B)
     implicit none
     class(function_1d) :: f
     real(kind=cfp), intent(in) :: A,B
     real(kind=cfp) :: gl1d

     integer :: i, j
     real(kind=cfp) :: x_AB(2*n_7+1), w_AB(2*n_7+1), val

        call gl_expand_A_B(x_7,w_7,n_7,x_AB,w_AB,A,B)

        gl1d = 0.0_cfp
        do i=1,2*n_7+1
           val = f%eval(x_AB(i))
           gl1d = gl1d + w_AB(i)*val
        enddo !i

        f%neval = f%neval + 2*n_7+1

 end function gl1d

 !> Adaptive 2D quadrature on rectangle based on Gauss-Kronrod rule. The algorithm is based on that of Romanowski published in Int. J. Q. Chem.
 !> \param [in] f The 2D function to be integrated.
 !> \param [in] Ax Rectangle X-coordinate start.
 !> \param [in] Bx Rectangle X-coordinate end.
 !> \param [in] Ay Rectangle Y-coordinate start.
 !> \param [in] By Rectangle Y-coordinate end.
 !> \param [in] eps Required relative precision for the integral.
 !> \param [in] Qest Optional: an estimate of the integral over the specified rectangle as obtained by a call to gl2d. 
 !>                  Note that the estimate must be obtained using gl2d for the algorithm to proceed correctly since it relies
 !>                  on the fixed-point G-K quadrature to calculate integrals on the sub-rectangles.
 recursive function quad2d(f,Ax,Bx,Ay,By,eps,Qest) result(I)
     implicit none
     class(function_2d) :: f
     real(kind=cfp), intent(in) :: Ax,Bx,Ay,By,eps
     real(kind=cfp), optional, intent(in) :: Qest

     real(kind=cfp) :: Xhalf, Yhalf, QP1, QP2, QP3, QP4, Q, I

        !Quit if too many sub-divisions have been performed: todo this should trigger an error/warning message!!
!        if (f%ndiv > f%max_div) return

        Xhalf = (Ax+Bx)*0.5_cfp
        Yhalf = (Ay+By)*0.5_cfp

        !If the estimate of the quadrature on this rectangle using a previous call to gl2d is given then use it otherwise calculate it.
        if (present(Qest)) then
           Q = Qest
        else
           Q = gl2d(f,Ax,Bx,Ay,By)
        endif

        !Calculate quadratures on the four sub-rectangles
        QP1 = gl2d(f,Ax,Xhalf,Ay,Yhalf)
        QP2 = gl2d(f,Ax,Xhalf,Yhalf,By)
        QP3 = gl2d(f,Xhalf,Bx,Yhalf,By)
        QP4 = gl2d(f,Xhalf,Bx,Ay,Yhalf)

        I = QP1+QP2+QP3+QP4

        !Continue recursively on each rectangle if the desired precision has not been reached
        !todo instead of dividing all four try dividing the largest one, then the second, etc. and see if we converge without improving all at the same time.
        if (abs((Q-I)/I) > eps) then
           Q = quad2d(f,Ax,Xhalf,Ay,Yhalf,eps,QP1)
           I = Q
           Q = quad2d(f,Ax,Xhalf,Yhalf,By,eps,QP2)
           I = I + Q
           Q = quad2d(f,Xhalf,Bx,Yhalf,By,eps,QP3)
           I = I + Q
           Q = quad2d(f,Xhalf,Bx,Ay,Yhalf,eps,QP4)
           I = I + Q
        endif

        f%ndiv = f%ndiv + 1

 end function quad2d

 !> 2D Quadrature on rectangle using the Gauss-Kronrod rule of order 8. The meaning of the input variables is identical to quad2d input parameters.
 function gl2d(f,Ax,Bx,Ay,By)
     implicit none
     class(function_2d) :: f
     real(kind=cfp), intent(in) :: Ax,Bx,Ay,By
     real(kind=cfp) :: gl2d

     !Abscissae and weights for the Gauss-Kronrod rule for interval [0,1] obtained using the Mathematica command:
     !n2d=8;
     !crule=NIntegrate`CartesianRuleData[{{"GaussKronrodRule","GaussPoints"->n2d},{"GaussKronrodRule","GaussPoints"->n2d}},36]
     integer, parameter :: n = 8
     real(kind=cfp), parameter :: x(2*n+1) = (/ &
        0.003310062059141922032055965490164602_cfp, &
        0.019855071751231884158219565715263505_cfp, &
        0.052939546576271789025819491230874338_cfp, &
        0.101666761293186630204223031762084782_cfp, &
        0.163822964527420661421844630953584475_cfp, &
        0.237233795041835507091130475405376825_cfp, &
        0.319649451035934021403725688851554280_cfp, &
        0.408282678752175097530261928819908010_cfp, &
        0.500000000000000000000000000000000000_cfp, &
        0.591717321247824902469738071180091990_cfp, &
        0.680350548964065978596274311148445720_cfp, &
        0.762766204958164492908869524594623175_cfp, &
        0.836177035472579338578155369046415525_cfp, &
        0.898333238706813369795776968237915218_cfp, &
        0.947060453423728210974180508769125662_cfp, &
        0.980144928248768115841780434284736495_cfp, &
        0.996689937940858077967944034509835398_cfp  &
     /)
     real(kind=cfp), parameter :: w(2*n+1) = (/ &
        0.00891119166035517757639348060137489490_cfp, &
        0.024719697501069654250181984723498447_cfp,   &
        0.0412411494656791653443125967228039477_cfp,  &
        0.055823185413419806611054079466970742_cfp,   &
        0.0681315546275861076311693726272531016_cfp,  &
        0.078326303084094200245124044243484369_cfp,   &
        0.0860353042776056559286474401019285433_cfp,  &
        0.09070001253401732153087426258627522_cfp,    &
        0.0922232028723458217644854778528214649_cfp,  &
        0.09070001253401732153087426258627522_cfp,    &
        0.0860353042776056559286474401019285433_cfp,  &
        0.078326303084094200245124044243484369_cfp,   &
        0.0681315546275861076311693726272531016_cfp,  &
        0.055823185413419806611054079466970742_cfp,   &
        0.0412411494656791653443125967228039477_cfp,  &
        0.024719697501069654250181984723498447_cfp,   &
        0.00891119166035517757639348060137489490_cfp  &
     /)

     integer :: i, j
     real(kind=cfp) :: Lx, Ly

        Lx = Bx-Ax
        Ly = By-Ay

        gl2d = 0.0_cfp
        do i=1,2*n+1
           do j=1,2*n+1
              gl2d = gl2d + w(i)*w(j)*f%eval(x(i)*Lx+Ax,x(j)*Ly+Ay)
           enddo !j
        enddo !i
        gl2d = gl2d*Lx*Ly

        f%neval = f%neval + (2*n+1)**2

 end function gl2d

!>***BEGIN PROLOGUE  DQAGS
!>***PURPOSE  The routine calculates an approximation result to a given
!>            Definite integral  I = Integral of F over (A,B),
!>            Hopefully satisfying following claim for accuracy
!>            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A1
!>***TYPE      real(kind=cfp) (QAGS-S, DQAGS-D)
!>***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
!>             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
!>             QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>        Computation of a definite integral
!>        Standard fortran subroutine
!>        Double precision version
!>
!>
!>        PARAMETERS
!>         ON ENTRY
!>            F      - class(function_1d)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>            A      - Double precision
!>                     Lower limit of integration
!>
!>            B      - Double precision
!>                     Upper limit of integration
!>
!>            EPSABS - Double precision
!>                     Absolute accuracy requested
!>            EPSREL - Double precision
!>                     Relative accuracy requested
!>                     If  EPSABS.LE.0
!>                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
!>                     The routine will end with IER = 6.
!>
!>         ON RETURN
!>            RESULT - Double precision
!>                     Approximation to the integral
!>
!>            ABSERR - Double precision
!>                     Estimate of the modulus of the absolute error,
!>                     which should equal or exceed ABS(I-RESULT)
!>
!>            NEVAL  - Integer
!>                     Number of integrand evaluations
!>
!>            IER    - Integer
!>                     IER = 0 Normal and reliable termination of the
!>                             routine. It is assumed that the requested
!>                             accuracy has been achieved.
!>                     IER.GT.0 Abnormal termination of the routine
!>                             The estimates for integral and error are
!>                             less reliable. It is assumed that the
!>                             requested accuracy has not been achieved.
!>            ERROR MESSAGES
!>                     IER = 1 Maximum number of subdivisions allowed
!>                             has been achieved. One can allow more sub-
!>                             divisions by increasing the value of LIMIT
!>                             (and taking the according dimension
!>                             adjustments into account. However, if
!>                             this yields no improvement it is advised
!>                             to analyze the integrand in order to
!>                             determine the integration difficulties. If
!>                             the position of a local difficulty can be
!>                             determined (E.G. SINGULARITY,
!>                             DISCONTINUITY WITHIN THE INTERVAL) one
!>                             will probably gain from splitting up the
!>                             interval at this point and calling the
!>                             integrator on the subranges. If possible,
!>                             an appropriate special-purpose integrator
!>                             should be used, which is designed for
!>                             handling the type of difficulty involved.
!>                         = 2 The occurrence of roundoff error is detec-
!>                             ted, which prevents the requested
!>                             tolerance from being achieved.
!>                             The error may be under-estimated.
!>                         = 3 Extremely bad integrand behaviour
!>                             occurs at some points of the integration
!>                             interval.
!>                         = 4 The algorithm does not converge.
!>                             Roundoff error is detected in the
!>                             Extrapolation table. It is presumed that
!>                             the requested tolerance cannot be
!>                             achieved, and that the returned result is
!>                             the best which can be obtained.
!>                         = 5 The integral is probably divergent, or
!>                             slowly convergent. It must be noted that
!>                             divergence can occur with any other value
!>                             of IER.
!>                         = 6 The input is invalid, because
!>                             (EPSABS.LE.0 AND
!>                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)
!>                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
!>                             RESULT, ABSERR, NEVAL, LAST are set to
!>                             zero.  Except when LIMIT or LENW is
!>                             invalid, IWORK(1), WORK(LIMIT*2+1) and
!>                             WORK(LIMIT*3+1) are set to zero, WORK(1)
!>                             is set to A and WORK(LIMIT+1) TO B.
!>
!>         DIMENSIONING PARAMETERS
!>            LIMIT - Integer
!>                    DIMENSIONING PARAMETER FOR IWORK
!>                    LIMIT determines the maximum number of subintervals
!>                    in the partition of the given integration interval
!>                    (A,B), LIMIT.GE.1.
!>                    IF LIMIT.LT.1, the routine will end with IER = 6.
!>
!>            LENW  - Integer
!>                    DIMENSIONING PARAMETER FOR WORK
!>                    LENW must be at least LIMIT*4.
!>                    If LENW.LT.LIMIT*4, the routine will end
!>                    with IER = 6.
!>
!>            LAST  - Integer
!>                    On return, LAST equals the number of subintervals
!>                    produced in the subdivision process, determines the
!>                    number of significant elements actually in the WORK
!>                    Arrays.
!>
!>         WORK ARRAYS
!>            IWORK - Integer
!>                    Vector of dimension at least LIMIT, the first K
!>                    elements of which contain pointers
!>                    to the error estimates over the subintervals
!>                    such that WORK(LIMIT*3+IWORK(1)),... ,
!>                    WORK(LIMIT*3+IWORK(K)) form a decreasing
!>                    sequence, with K = LAST IF LAST.LE.(LIMIT/2+2),
!>                    and K = LIMIT+1-LAST otherwise
!>
!>            WORK  - Double precision
!>                    Vector of dimension at least LENW
!>                    on return
!>                    WORK(1), ..., WORK(LAST) contain the left
!>                     end-points of the subintervals in the
!>                     partition of (A,B),
!>                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
!>                     the right end-points,
!>                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
!>                     the integral approximations over the subintervals,
!>                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
!>                     contain the error estimates.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  DQAGSE, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>***END PROLOGUE  DQAGS
!>
!>
      SUBROUTINE DQAGS (F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
      class(function_1d) :: F
      real(kind=cfp) A,ABSERR,B,EPSABS,EPSREL,RESULT,WORK!,F
      INTEGER IER,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
!
      DIMENSION IWORK(*),WORK(*)
!
      !DECLARE F AS: real(kind=cfp), EXTERNAL :: F
      !EXTERNAL F
!
!         CHECK VALIDITY OF LIMIT AND LENW.
!
!***FIRST EXECUTABLE STATEMENT  DQAGS
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0_cfp
      ABSERR = 0.0_cfp
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
!
!         PREPARE CALL FOR DQAGSE.
!
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
!
      CALL DQAGSE(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
!
!         CALL ERROR HANDLER IF NECESSARY.
!
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) THEN
         PRINT *,RESULT
         CALL XERMSG ('SLATEC', 'DQAGS', 'ABNORMAL RETURN', IER, LVL)
      ENDIF
      RETURN
      END SUBROUTINE DQAGS

!>***BEGIN PROLOGUE  DQAGSE
!>***PURPOSE  The routine calculates an approximation result to a given
!>            definite integral I = Integral of F over (A,B),
!>            hopefully satisfying following claim for accuracy
!>            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A1
!>***TYPE      real(kind=cfp) (QAGSE-S, DQAGSE-D)
!>***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
!>             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
!>             QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>        Computation of a definite integral
!>        Standard fortran subroutine
!>        Double precision version
!>
!>        PARAMETERS
!>         ON ENTRY
!>            F      - class(function_1d)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>            A      - Double precision
!>                     Lower limit of integration
!>
!>            B      - Double precision
!>                     Upper limit of integration
!>
!>            EPSABS - Double precision
!>                     Absolute accuracy requested
!>            EPSREL - Double precision
!>                     Relative accuracy requested
!>                     If  EPSABS.LE.0
!>                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
!>                     the routine will end with IER = 6.
!>
!>            LIMIT  - Integer
!>                     Gives an upper bound on the number of subintervals
!>                     in the partition of (A,B)
!>
!>         ON RETURN
!>            RESULT - Double precision
!>                     Approximation to the integral
!>
!>            ABSERR - Double precision
!>                     Estimate of the modulus of the absolute error,
!>                     which should equal or exceed ABS(I-RESULT)
!>
!>            NEVAL  - Integer
!>                     Number of integrand evaluations
!>
!>            IER    - Integer
!>                     IER = 0 Normal and reliable termination of the
!>                             routine. It is assumed that the requested
!>                             accuracy has been achieved.
!>                     IER.GT.0 Abnormal termination of the routine
!>                             the estimates for integral and error are
!>                             less reliable. It is assumed that the
!>                             requested accuracy has not been achieved.
!>            ERROR MESSAGES
!>                         = 1 Maximum number of subdivisions allowed
!>                             has been achieved. One can allow more sub-
!>                             divisions by increasing the value of LIMIT
!>                             (and taking the according dimension
!>                             adjustments into account). However, if
!>                             this yields no improvement it is advised
!>                             to analyze the integrand in order to
!>                             determine the integration difficulties. If
!>                             the position of a local difficulty can be
!>                             determined (e.g. singularity,
!>                             discontinuity within the interval) one
!>                             will probably gain from splitting up the
!>                             interval at this point and calling the
!>                             integrator on the subranges. If possible,
!>                             an appropriate special-purpose integrator
!>                             should be used, which is designed for
!>                             handling the type of difficulty involved.
!>                         = 2 The occurrence of roundoff error is detec-
!>                             ted, which prevents the requested
!>                             tolerance from being achieved.
!>                             The error may be under-estimated.
!>                         = 3 Extremely bad integrand behaviour
!>                             occurs at some points of the integration
!>                             interval.
!>                         = 4 The algorithm does not converge.
!>                             Roundoff error is detected in the
!>                             extrapolation table.
!>                             It is presumed that the requested
!>                             tolerance cannot be achieved, and that the
!>                             returned result is the best which can be
!>                             obtained.
!>                         = 5 The integral is probably divergent, or
!>                             slowly convergent. It must be noted that
!>                             divergence can occur with any other value
!>                             of IER.
!>                         = 6 The input is invalid, because
!>                             EPSABS.LE.0 and
!>                             EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28).
!>                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
!>                             IORD(1) and ELIST(1) are set to zero.
!>                             ALIST(1) and BLIST(1) are set to A and B
!>                             respectively.
!>
!>            ALIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the left end points
!>                     of the subintervals in the partition of the
!>                     given integration range (A,B)
!>
!>            BLIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the right end points
!>                     of the subintervals in the partition of the given
!>                     integration range (A,B)
!>
!>            RLIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the integral
!>                     approximations on the subintervals
!>
!>            ELIST  - Double precision
!>                     Vector of dimension at least LIMIT, the first
!>                      LAST  elements of which are the moduli of the
!>                     absolute error estimates on the subintervals
!>
!>            IORD   - Integer
!>                     Vector of dimension at least LIMIT, the first K
!>                     elements of which are pointers to the
!>                     error estimates over the subintervals,
!>                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
!>                     form a decreasing sequence, with K = LAST
!>                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
!>                     otherwise
!>
!>            LAST   - Integer
!>                     Number of subintervals actually produced in the
!>                     subdivision process
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, DQELG, DQK21, DQPSRT
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>***END PROLOGUE  DQAGSE
!>
!>            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
!>            LIMEXP IN SUBROUTINE DQELG (RLIST2 SHOULD BE OF DIMENSION
!>            (LIMEXP+2) AT LEAST).
!>
!>            LIST OF MAJOR VARIABLES
!>            -----------------------
!>
!>           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
!>                       CONSIDERED UP TO NOW
!>           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
!>                       CONSIDERED UP TO NOW
!>           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
!>                       (ALIST(I),BLIST(I))
!>           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2 CONTAINING
!>                       THE PART OF THE EPSILON TABLE WHICH IS STILL
!>                       NEEDED FOR FURTHER COMPUTATIONS
!>           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
!>           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
!>                       ESTIMATE
!>           ERRMAX    - ELIST(MAXERR)
!>           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
!>                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
!>           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
!>           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
!>           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
!>                       ABS(RESULT))
!>           *****1    - VARIABLE FOR THE LEFT INTERVAL
!>           *****2    - VARIABLE FOR THE RIGHT INTERVAL
!>           LAST      - INDEX FOR SUBDIVISION
!>           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
!>           NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
!>                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
!>                       INTEGRAL HAS BEEN OBTAINED IT IS PUT IN
!>                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
!>                       BY ONE.
!>           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
!>                       TO NOW, MULTIPLIED BY 1.5
!>           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
!>                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
!>           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS
!>                       ATTEMPTING TO PERFORM EXTRAPOLATION I.E. BEFORE
!>                       SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO
!>                       DECREASE THE VALUE OF ERLARG.
!>           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
!>                       IS NO LONGER ALLOWED (TRUE VALUE)
!>
!>            MACHINE DEPENDENT CONSTANTS
!>            ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!>           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
      SUBROUTINE DQAGSE (F, A, B, EPSABS, EPSREL, LIMIT, RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
      use const_gbl, only: max_epstab
      class(function_1d) :: F
      real(kind=cfp) A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,BLIST,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2, &
     &  DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,ERRMAX, &
     &  ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,OFLOW,RESABS,RESEPS,RESULT,&
     &  RES3LA,RLIST,RLIST2,SMALL,UFLOW
      real(kind=cfp) cfp_dummy
      INTEGER ID,IER,IERRO,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,KSGN,KTMIN,LAST,LIMIT,MAXERR,NEVAL,NRES,NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
!
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),RES3LA(3),RLIST(*),RLIST2(max_epstab)
!
!***FIRST EXECUTABLE STATEMENT  DQAGSE
      EPMACH = F1MACH(4,cfp_dummy)
!
!            TEST ON VALIDITY OF PARAMETERS
!            ------------------------------
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0_cfp
      ABSERR = 0.0_cfp
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0_cfp
      ELIST(1) = 0.0_cfp
      IF(EPSABS.LE.0.0_cfp.AND.EPSREL.LT.MAX(0.5E+02_cfp*EPMACH,0.5E-28_cfp)) IER = 6
      IF(IER.EQ.6) GO TO 999
!
!           FIRST APPROXIMATION TO THE INTEGRAL
!           -----------------------------------
!
      UFLOW = F1MACH(1,cfp_dummy)
      OFLOW = F1MACH(2,cfp_dummy)
      IERRO = 0
      CALL DQK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
!
!           TEST ON ACCURACY.
!
      DRES = ABS(RESULT)
      ERRBND = MAX(EPSABS,EPSREL*DRES)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      IF(ABSERR.LE.1.0E+02_cfp*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS).OR.ABSERR.EQ.0.0_cfp) GO TO 140
!
!           INITIALIZATION
!           --------------
!
      RLIST2(1) = RESULT
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      ABSERR = OFLOW
      NRMAX = 1
      NRES = 0
      NUMRL2 = 2
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KSGN = -1
      IF(DRES.GE.(0.1E+01_cfp-0.5E+02_cfp*EPMACH)*DEFABS) KSGN = 1
!
!           MAIN DO-LOOP
!           ------------
!
      DO 90 LAST = 2,LIMIT
!
!           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR
!           ESTIMATE.
!
        A1 = ALIST(MAXERR)
        B1 = 0.5_cfp*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        CALL DQK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
!
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!           AND ERROR AND TEST FOR ACCURACY.
!
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 15
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1E-04_cfp*ABS(AREA12).OR.ERRO12.LT.0.99_cfp*ERRMAX) GO TO 10
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   10   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   15   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
!
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
!
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
!
!           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
!           EQUALS LIMIT.
!
        IF(LAST.EQ.LIMIT) IER = 1
!
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT A POINT OF THE INTEGRATION RANGE.
!
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1E+01_cfp+0.1E+03_cfp*EPMACH)*(ABS(A2)+0.1E+04_cfp*UFLOW)) IER = 4
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
        IF(ERROR2.GT.ERROR1) GO TO 20
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 30
   20   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
!
!           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
!           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
!
   30   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
        IF(ERRSUM.LE.ERRBND) GO TO 115
! ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0) GO TO 100
        IF(LAST.EQ.2) GO TO 80
        IF(NOEXT) GO TO 90
        ERLARG = ERLARG-ERLAST
        IF(ABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 40
!
!           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
!           SMALLEST INTERVAL.
!
        IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
        EXTRAP = .TRUE.
        NRMAX = 2
   40   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 60
!
!           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
!           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
!           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
!
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 50 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
! ***JUMP OUT OF DO-LOOP
          IF(ABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
          NRMAX = NRMAX+1
   50   CONTINUE
!
!           PERFORM EXTRAPOLATION.
!
   60   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        CALL DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1E-02_cfp*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 70
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
! ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LE.ERTEST) GO TO 100
!
!           PREPARE BISECTION OF THE SMALLEST INTERVAL.
!
   70   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 100
        MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5_cfp
        ERLARG = ERRSUM
        GO TO 90
   80   SMALL = ABS(B-A)*0.375_cfp
        ERLARG = ERRSUM
        ERTEST = ERRBND
        RLIST2(2) = AREA
   90 CONTINUE
!
!           SET FINAL RESULT AND ERROR ESTIMATE.
!           ------------------------------------
!
  100 IF(ABSERR.EQ.OFLOW) GO TO 115
      IF(IER+IERRO.EQ.0) GO TO 110
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0_cfp.AND.AREA.NE.0.0_cfp) GO TO 105
      IF(ABSERR.GT.ERRSUM) GO TO 115
      IF(AREA.EQ.0.0_cfp) GO TO 130
      GO TO 110
  105 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA)) GO TO 115
!
!           TEST ON DIVERGENCE.
!
  110 IF(KSGN.EQ.(-1).AND.MAX(ABS(RESULT),ABS(AREA)).LE.DEFABS*0.1E-01_cfp) GO TO 130
      IF(0.1E-01_cfp.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1E+03_cfp.OR.ERRSUM.GT.ABS(AREA)) IER = 6
      GO TO 130
!
!           COMPUTE GLOBAL INTEGRAL SUM.
!
  115 RESULT = 0.0_cfp
      DO 120 K = 1,LAST
         RESULT = RESULT+RLIST(K)
  120 CONTINUE
      ABSERR = ERRSUM
  130 IF(IER.GT.2) IER = IER-1
  140 NEVAL = 42*LAST-21
  999 RETURN
      END SUBROUTINE DQAGSE

!>***BEGIN PROLOGUE  DQELG
!>***SUBSIDIARY
!>***PURPOSE  The routine determines the limit of a given sequence of
!>            approximations, by means of the Epsilon algorithm of
!>            P.Wynn. An estimate of the absolute error is also given.
!>            The condensed Epsilon table is computed. Only those
!>            elements needed for the computation of the next diagonal
!>            are preserved.
!>***LIBRARY   SLATEC
!>***TYPE      real(kind=cfp) (QELG-S, DQELG-D)
!>***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Epsilon algorithm
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS
!>              N      - Integer
!>                       EPSTAB(N) contains the new element in the
!>                       first column of the epsilon table.
!>
!>              EPSTAB - Double precision
!>                       Vector of dimension 52 containing the elements
!>                       of the two lower diagonals of the triangular
!>                       epsilon table. The elements are numbered
!>                       starting at the right-hand corner of the
!>                       triangle.
!>
!>              RESULT - Double precision
!>                       Resulting approximation to the integral
!>
!>              ABSERR - Double precision
!>                       Estimate of the absolute error computed from
!>                       RESULT and the 3 previous results
!>
!>              RES3LA - Double precision
!>                       Vector of dimension 3 containing the last 3
!>                       results
!>
!>              NRES   - Integer
!>                       Number of calls to the routine
!>                       (should be zero at first call)
!>
!>***SEE ALSO  DQAGIE, DQAGOE, DQAGPE, DQAGSE
!>***ROUTINES CALLED  F1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>***END PROLOGUE  DQELG
!>
!>           LIST OF MAJOR VARIABLES
!>           -----------------------
!>
!>           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
!>           E1       ELEMENT IN THE EPSILON TABLE IS BASED
!>           E2
!>           E3                 E0
!>                        E3    E1    NEW
!>                              E2
!>           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
!>                    DIAGONAL
!>           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
!>           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
!>                    OF ERROR
!>
!>           MACHINE DEPENDENT CONSTANTS
!>           ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
!>           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
!>           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
!>           DIAGONAL OF THE EPSILON TABLE IS DELETED.
!>
      SUBROUTINE DQELG (N, EPSTAB, RESULT, ABSERR, RES3LA, NRES)
      use const_gbl, only: max_epstab
      real(kind=cfp) ABSERR, DELTA1, DELTA2, DELTA3, EPMACH, EPSINF, EPSTAB, ERROR, ERR1, ERR2, ERR3, E0, E1, E1ABS
      real(kind=cfp) E2, E3, OFLOW, RES, RESULT, RES3LA, SS, TOL1, TOL2, TOL3, cfp_dummy
      INTEGER I,IB,IB2,IE,INDX,K1,K2,K3,LIMEXP,N,NEWELM,NRES,NUM
      DIMENSION EPSTAB(max_epstab),RES3LA(3)
!***FIRST EXECUTABLE STATEMENT  DQELG
      EPMACH = F1MACH(4,cfp_dummy)
      OFLOW = F1MACH(2,cfp_dummy)
      NRES = NRES+1
      ABSERR = OFLOW
      RESULT = EPSTAB(N)
      IF(N.LT.3) GO TO 100
      LIMEXP = max_epstab-2
      EPSTAB(N+2) = EPSTAB(N)
      NEWELM = (N-1)/2
      EPSTAB(N) = OFLOW
      NUM = N
      K1 = N
      DO 40 I = 1,NEWELM
        K2 = K1-1
        K3 = K1-2
        RES = EPSTAB(K1+2)
        E0 = EPSTAB(K3)
        E1 = EPSTAB(K2)
        E2 = RES
        E1ABS = ABS(E1)
        DELTA2 = E2-E1
        ERR2 = ABS(DELTA2)
        TOL2 = MAX(ABS(E2),E1ABS)*EPMACH
        DELTA3 = E1-E0
        ERR3 = ABS(DELTA3)
        TOL3 = MAX(E1ABS,ABS(E0))*EPMACH
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
!
!           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
!           ACCURACY, CONVERGENCE IS ASSUMED.
!           RESULT = E2
!           ABSERR = ABS(E1-E0)+ABS(E2-E1)
!
        RESULT = RES
        ABSERR = ERR2+ERR3
! ***JUMP OUT OF DO-LOOP
        GO TO 100
   10   E3 = EPSTAB(K1)
        EPSTAB(K1) = E1
        DELTA1 = E1-E3
        ERR1 = ABS(DELTA1)
        TOL1 = MAX(E1ABS,ABS(E3))*EPMACH
!
!           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
!           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
!
        IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
        SS = 0.1E+01_cfp/DELTA1+0.1E+01_cfp/DELTA2-0.1E+01_cfp/DELTA3
        EPSINF = ABS(SS*E1)
!
!           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
!           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
!           OF N.
!
        IF(EPSINF.GT.0.1E-03_cfp) GO TO 30
   20   N = I+I-1
! ***JUMP OUT OF DO-LOOP
        GO TO 50
!
!           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
!           THE VALUE OF RESULT.
!
   30   RES = E1+0.1E+01_cfp/SS
        EPSTAB(K1) = RES
        K1 = K1-2
        ERROR = ERR2+ABS(RES-E2)+ERR3
        IF(ERROR.GT.ABSERR) GO TO 40
        ABSERR = ERROR
        RESULT = RES
   40 CONTINUE
!
!           SHIFT THE TABLE.
!
   50 IF(N.EQ.LIMEXP) N = 2*(LIMEXP/2)-1
      IB = 1
      IF((NUM/2)*2.EQ.NUM) IB = 2
      IE = NEWELM+1
      DO 60 I=1,IE
        IB2 = IB+2
        EPSTAB(IB) = EPSTAB(IB2)
        IB = IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      INDX = NUM-N+1
      DO 70 I = 1,N
        EPSTAB(I)= EPSTAB(INDX)
        INDX = INDX+1
   70 CONTINUE
   80 IF(NRES.GE.4) GO TO 90
      RES3LA(NRES) = RESULT
      ABSERR = OFLOW
      GO TO 100
!
!           COMPUTE ERROR ESTIMATE
!
   90 ABSERR = ABS(RESULT-RES3LA(3))+ABS(RESULT-RES3LA(2))+ABS(RESULT-RES3LA(1))
      RES3LA(1) = RES3LA(2)
      RES3LA(2) = RES3LA(3)
      RES3LA(3) = RESULT
  100 ABSERR = MAX(ABSERR,0.5E+01_cfp*EPMACH*ABS(RESULT))
      RETURN
      END SUBROUTINE DQELG

!>***BEGIN PROLOGUE  DQK21
!>***PURPOSE  To compute I = Integral of F over (A,B), with error
!>                           estimate
!>                       J = Integral of ABS(F) over (A,B)
!>***LIBRARY   SLATEC (QUADPACK)
!>***CATEGORY  H2A1A2
!>***TYPE      real(kind=cfp) (QK21-S, DQK21-D)
!>***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Integration rules
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS
!>            ON ENTRY
!>            F      - class(function_1d)
!>                     Function whose method 'eval' defines the integrand
!>                     Function F(X).
!>
!>              A      - Double precision
!>                       Lower limit of integration
!>
!>              B      - Double precision
!>                       Upper limit of integration
!>
!>            ON RETURN
!>              RESULT - Double precision
!>                       Approximation to the integral I
!>                       RESULT is computed by applying the 21-POINT
!>                       KRONROD RULE (RESK) obtained by optimal addition
!>                       of abscissae to the 10-POINT GAUSS RULE (RESG).
!>
!>              ABSERR - Double precision
!>                       Estimate of the modulus of the absolute error,
!>                       which should not exceed ABS(I-RESULT)
!>
!>              RESABS - Double precision
!>                       Approximation to the integral J
!>
!>              RESASC - Double precision
!>                       Approximation to the integral of ABS(F-I/(B-A))
!>                       over (A,B)
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>***END PROLOGUE  DQK21
!>
!>
!>           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!>           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!>           CORRESPONDING WEIGHTS ARE GIVEN.
!>
!>           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
!>                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
!>                    GAUSS RULE
!>                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!>                    ADDED TO THE 10-POINT GAUSS RULE
!>
!>           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
!>
!>           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
!>
!>
!> GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
!> AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
!> BELL LABS, NOV. 1981.
!>
!>
!>           LIST OF MAJOR VARIABLES
!>           -----------------------
!>
!>           CENTR  - MID POINT OF THE INTERVAL
!>           HLGTH  - HALF-LENGTH OF THE INTERVAL
!>           ABSC   - ABSCISSA
!>           FVAL*  - FUNCTION VALUE
!>           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
!>           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
!>           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!>                    I.E. TO I/(B-A)
!>
!>
!>           MACHINE DEPENDENT CONSTANTS
!>           ---------------------------
!>
!>           EPMACH IS THE LARGEST RELATIVE SPACING.
!>           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!>
      SUBROUTINE DQK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
      class(function_1d) :: F
      real(kind=cfp) A, ABSC, ABSERR, B, CENTR, DHLGTH, EPMACH, FC, FSUM, FVAL1, FVAL2, FV1, FV2, HLGTH, RESABS, RESASC
      real(kind=cfp) RESG, RESK, RESKH, RESULT, UFLOW, WG, WGK, XGK, cfp_dummy, FV(21), POINTS(21)
      INTEGER J,JTW,JTWM1
!
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
!
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.066671344308688137593568809893332_cfp /
      DATA WG  (  2) / 0.149451349150580593145776339657697_cfp /
      DATA WG  (  3) / 0.219086362515982043995534934228163_cfp /
      DATA WG  (  4) / 0.269266719309996355091226921569469_cfp /
      DATA WG  (  5) / 0.295524224714752870173892994651338_cfp /
!
      DATA XGK (  1) / 0.995657163025808080735527280689003_cfp /
      DATA XGK (  2) / 0.973906528517171720077964012084452_cfp /
      DATA XGK (  3) / 0.930157491355708226001207180059508_cfp /
      DATA XGK (  4) / 0.865063366688984510732096688423493_cfp /
      DATA XGK (  5) / 0.780817726586416897063717578345042_cfp /
      DATA XGK (  6) / 0.679409568299024406234327365114874_cfp /
      DATA XGK (  7) / 0.562757134668604683339000099272694_cfp /
      DATA XGK (  8) / 0.433395394129247190799265943165784_cfp /
      DATA XGK (  9) / 0.294392862701460198131126603103866_cfp /
      DATA XGK ( 10) / 0.148874338981631210884826001129720_cfp /
      DATA XGK ( 11) / 0.000000000000000000000000000000000_cfp /
!
      DATA WGK (  1) / 0.011694638867371874278064396062192_cfp /
      DATA WGK (  2) / 0.032558162307964727478818972459390_cfp /
      DATA WGK (  3) / 0.054755896574351996031381300244580_cfp /
      DATA WGK (  4) / 0.075039674810919952767043140916190_cfp /
      DATA WGK (  5) / 0.093125454583697605535065465083366_cfp /
      DATA WGK (  6) / 0.109387158802297641899210590325805_cfp /
      DATA WGK (  7) / 0.123491976262065851077958109831074_cfp /
      DATA WGK (  8) / 0.134709217311473325928054001771707_cfp /
      DATA WGK (  9) / 0.142775938577060080797094273138717_cfp /
      DATA WGK ( 10) / 0.147739104901338491374841515972068_cfp /
      DATA WGK ( 11) / 0.149445554002916905664936468389821_cfp /
!
!***FIRST EXECUTABLE STATEMENT  DQK21
      EPMACH = F1MACH(4,cfp_dummy)
      UFLOW = F1MACH(1,cfp_dummy)
!
      CENTR = 0.5_cfp*(A+B)
      HLGTH = 0.5_cfp*(B-A)
      DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
!           FOR functions_1d EVALUATE ONE POINT AT A TIME.
!           FOR functions_1d_many FIRST FIND THE WHOLE
!           GRID AND THEN EVALUATE IN ALL GRIDPOINTS AT 
!           ONCE.
!

      SELECT TYPE (F)
         CLASS IS (function_1d)

            RESG = 0.0_cfp
            FC = F%eval(CENTR)
            RESK = WGK(11)*FC
            RESABS = ABS(RESK)
            DO 10 J=1,5
              JTW = 2*J
              ABSC = HLGTH*XGK(JTW)
              FVAL1 = F%eval(CENTR-ABSC)
              FVAL2 = F%eval(CENTR+ABSC)
              FV1(JTW) = FVAL1
              FV2(JTW) = FVAL2
              FSUM = FVAL1+FVAL2
              RESG = RESG+WG(J)*FSUM
              RESK = RESK+WGK(JTW)*FSUM
              RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
         10 CONTINUE
            DO 15 J = 1,5
              JTWM1 = 2*J-1
              ABSC = HLGTH*XGK(JTWM1)
              FVAL1 = F%eval(CENTR-ABSC)
              FVAL2 = F%eval(CENTR+ABSC)
              FV1(JTWM1) = FVAL1
              FV2(JTWM1) = FVAL2
              FSUM = FVAL1+FVAL2
              RESK = RESK+WGK(JTWM1)*FSUM
              RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
         15 CONTINUE
            RESKH = RESK*0.5_cfp
            RESASC = WGK(11)*ABS(FC-RESKH)
            DO 20 J=1,10
              RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
         20 CONTINUE
            RESULT = RESK*HLGTH
            RESABS = RESABS*DHLGTH
            RESASC = RESASC*DHLGTH
            ABSERR = ABS((RESK-RESG)*HLGTH)
            IF(RESASC.NE.0.0_cfp.AND.ABSERR.NE.0.0_cfp) ABSERR = RESASC*MIN(0.1E+01_cfp,(0.2E+03_cfp*ABSERR/RESASC)**1.5E+00_cfp)
            IF(RESABS.GT.UFLOW/(0.5E+02_cfp*EPMACH)) ABSERR = MAX((EPMACH*0.5E+02_cfp)*RESABS,ABSERR)
            RETURN

         CLASS IS (function_1d_many)
            
            K = 1
            POINTS(K) = CENTR
            DO J=1,5
              JTW = 2*J
              ABSC = HLGTH*XGK(JTW)
              K = K + 1
              POINTS(K) = CENTR-ABSC
              K = K + 1
              POINTS(K) = CENTR+ABSC
            ENDDO
            DO J = 1,5
              JTWM1 = 2*J-1
              ABSC = HLGTH*XGK(JTWM1)
              K = K + 1
              POINTS(K) = CENTR-ABSC
              K = K + 1
              POINTS(K) = CENTR+ABSC
            ENDDO
      !           NOW EVALUATE THE INTEGRAND AT THESE POINTS.
            FV = F%eval(POINTS,K)
      !
            RESG = 0.0_cfp
            K = 1
            FC = FV(K) !F%eval(CENTR)
            RESK = WGK(11)*FC
            RESABS = ABS(RESK)
            DO 30 J=1,5
              JTW = 2*J
              !ABSC = HLGTH*XGK(JTW)
              K = K + 1
              FVAL1 = FV(K) !F%eval(CENTR-ABSC)
              K = K + 1
              FVAL2 = FV(K) !F%eval(CENTR+ABSC)
              FV1(JTW) = FVAL1
              FV2(JTW) = FVAL2
              FSUM = FVAL1+FVAL2
              RESG = RESG+WG(J)*FSUM
              RESK = RESK+WGK(JTW)*FSUM
              RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
         30 CONTINUE
            DO 35 J = 1,5
              JTWM1 = 2*J-1
              !ABSC = HLGTH*XGK(JTWM1)
              K = K + 1
              FVAL1 = FV(K) !F%eval(CENTR-ABSC)
              K = K + 1
              FVAL2 = FV(K) !F%eval(CENTR+ABSC)
              FV1(JTWM1) = FVAL1
              FV2(JTWM1) = FVAL2
              FSUM = FVAL1+FVAL2
              RESK = RESK+WGK(JTWM1)*FSUM
              RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
         35 CONTINUE
            RESKH = RESK*0.5_cfp
            RESASC = WGK(11)*ABS(FC-RESKH)
            DO 40 J=1,10
              RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
         40 CONTINUE
            RESULT = RESK*HLGTH
            RESABS = RESABS*DHLGTH
            RESASC = RESASC*DHLGTH
            ABSERR = ABS((RESK-RESG)*HLGTH)
            IF(RESASC.NE.0.0_cfp.AND.ABSERR.NE.0.0_cfp) ABSERR = RESASC*MIN(0.1E+01_cfp,(0.2E+03_cfp*ABSERR/RESASC)**1.5E+00_cfp)
            IF(RESABS.GT.UFLOW/(0.5E+02_cfp*EPMACH)) ABSERR = MAX((EPMACH*0.5E+02_cfp)*RESABS,ABSERR)
            RETURN

         CLASS DEFAULT
            call xermsg('general_quadrature','DQK21','The type of function on input must be either function_1d or &
                        &function_1d_many.',1,2)
      END SELECT

      END SUBROUTINE DQK21

!>***BEGIN PROLOGUE  DQPSRT
!>***SUBSIDIARY
!>***PURPOSE  This routine maintains the descending ordering in the
!>            list of the local error estimated resulting from the
!>            interval subdivision process. At each call two error
!>            estimates are inserted using the sequential search
!>            method, top-down for the largest error estimate and
!>            bottom-up for the smallest error estimate.
!>***LIBRARY   SLATEC
!>***TYPE      real(kind=cfp) (QPSRT-S, DQPSRT-D)
!>***KEYWORDS  SEQUENTIAL SORTING
!>***AUTHOR  Piessens, Robert
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>           de Doncker, Elise
!>             Applied Mathematics and Programming Division
!>             K. U. Leuven
!>***DESCRIPTION
!>
!>           Ordering routine
!>           Standard fortran subroutine
!>           Double precision version
!>
!>           PARAMETERS (MEANING AT OUTPUT)
!>              LIMIT  - Integer
!>                       Maximum number of error estimates the list
!>                       can contain
!>
!>              LAST   - Integer
!>                       Number of error estimates currently in the list
!>
!>              MAXERR - Integer
!>                       MAXERR points to the NRMAX-th largest error
!>                       estimate currently in the list
!>
!>              ERMAX  - Double precision
!>                       NRMAX-th largest error estimate
!>                       ERMAX = ELIST(MAXERR)
!>
!>              ELIST  - Double precision
!>                       Vector of dimension LAST containing
!>                       the error estimates
!>
!>              IORD   - Integer
!>                       Vector of dimension LAST, the first K elements
!>                       of which contain pointers to the error
!>                       estimates, such that
!>                       ELIST(IORD(1)),...,  ELIST(IORD(K))
!>                       form a decreasing sequence, with
!>                       K = LAST if LAST.LE.(LIMIT/2+2), and
!>                       K = LIMIT+1-LAST otherwise
!>
!>              NRMAX  - Integer
!>                       MAXERR = IORD(NRMAX)
!>
!>***SEE ALSO  DQAGE, DQAGIE, DQAGPE, DQAWSE
!>***ROUTINES CALLED  (NONE)
!>***REVISION HISTORY  (YYMMDD)
!>   800101  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>***END PROLOGUE  DQPSRT
!>
      SUBROUTINE DQPSRT (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)
      real(kind=cfp) ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,NRMAX
      DIMENSION ELIST(*),IORD(*)
!
!           CHECK WHETHER THE LIST CONTAINS MORE THAN
!           TWO ERROR ESTIMATES.
!
!***FIRST EXECUTABLE STATEMENT  DQPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
!
!           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
!           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
!           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
!           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
!
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
!
!           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
!           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
!           SUBDIVISIONS STILL ALLOWED.
!
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
!
!           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
!           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
!
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
!
!           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
!
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
! ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
!
!           SET MAXERR AND ERMAX.
!
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END SUBROUTINE DQPSRT

!     Abscissas and weights of Gaussian Quadrature
!     produced by Mathematica code: AWGQ
!ccccc
!     for int_a^b{f(z(x))W(x)dx}
!     argument    z(x) = x
!     weight      W(x) = x**(2*L)*exp(-x**2)
!     lower limit a = 0
!     upper limit b = +Infinity
!     GQ order n==18,L=0,...,24
!ccccc
      subroutine quadrature_u_integral(x,w,n,L)
!     x: abscissas
!     w: weights
      implicit none
      integer, intent(in) :: n, L
      real(kind=cfp), intent(out) :: x(n),w(n)
!
      if (n .ne. 18 .or. L < 0 .or. L > 24) then
         call xermsg('general_quadrature','quadrature_u_integral', &
                     'On input n must be equal to 18 and L must be in the range [0;24].',1,1)
      endif

      if(n.eq.18 .and. L .eq.0) then
        x(1)=1.66490322202372327447714429806571665E-2_cfp
        w(1)=4.26142619814402343016279420043364307E-2_cfp
        x(2)=8.68590621084086644675129815579114681E-2_cfp
        w(2)=9.65666206000141056170006679582123331E-2_cfp
        x(3)=2.09868348130364142021838753038192800E-1_cfp
        w(3)=1.41519805800333141668255899625808090E-1_cfp
        x(4)=3.80820011068889100983758151725581925E-1_cfp
        w(4)=1.66987766912260043155372327141828977E-1_cfp
        x(5)=5.94030720753567083406433802547308049E-1_cfp
        w(5)=1.63315996982124497361219918422153983E-1_cfp
        x(6)=8.43863021635420192757674184557786877E-1_cfp
        w(6)=1.30699466898223965583228244414552703E-1_cfp
        x(7)=1.12530737639272964828726611490096105E0_cfp
        w(7)=8.33782286268739221855899400556967617E-2_cfp
        x(8)=1.43428771048618524292176673609752148E0_cfp
        w(8)=4.11112305586426925178613509285170154E-2_cfp
        x(9)=1.76777393564493681913999434323957822E0_cfp
        w(9)=1.51576539874531054894670516472746539E-2_cfp
        x(10)=2.12379869512188868535534629099934566E0_cfp
        w(10)=4.03326905747473179178466967905519779E-3_cfp
        x(11)=2.50145945877649262471348943058897348E0_cfp
        w(11)=7.44434183403738436876331846854961792E-4_cfp
        x(12)=2.90097223563890611834014365316071673E0_cfp
        w(12)=9.09462923551198320417749396628607376E-5_cfp
        x(13)=3.32384856894200723550010778630591449E0_cfp
        w(13)=6.93223179081282168585655849742917362E-6_cfp
        x(14)=3.77331513435005717325307716948541392E0_cfp
        w(14)=3.04417379930169850289910740424100627E-7_cfp
        x(15)=4.25524711266143974366696666598420946E0_cfp
        w(15)=6.85706226732795221220972953252013969E-9_cfp
        x(16)=4.78037806093042877680970682222312155E0_cfp
        w(16)=6.57362956131057685569607827427012282E-11_cfp
        x(17)=5.37053657981972462591916408220001534E0_cfp
        w(17)=1.89340654173428573497308066005032294E-13_cfp
        x(18)=6.08421686390342864341227897100851443E0_cfp
        w(18)=6.91219900671971931466795228079604552E-17_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.1) then
        x(1)=6.96597387614463863461268360119255008E-2_cfp
        w(1)=4.38010710869322772091139391683842559E-4_cfp
        x(2)=1.84250233164821135903889996522952205E-1_cfp
        w(2)=4.52345528757740956163124044016956708E-3_cfp
        x(3)=3.43833974822474353286235213086867212E-1_cfp
        w(3)=1.89586600144798759321653709264095944E-2_cfp
        x(4)=5.43529709713899590471053177209158257E-1_cfp
        w(4)=4.79414527587109422925230788080029104E-2_cfp
        x(5)=7.78312511818120447560133225962850912E-1_cfp
        w(5)=8.28796899353481576876922702502378390E-2_cfp
        x(6)=1.04358452637828314466385289587917067E0_cfp
        w(6)=1.02312458750870820684993378149978325E-1_cfp
        x(7)=1.33547731039452715326162323045169677E0_cfp
        w(7)=9.11511089627791459098751663933460109E-2_cfp
        x(8)=1.65098002448728840397045614887270270E0_cfp
        w(8)=5.82948234192263771924611382601405478E-2_cfp
        x(9)=1.98796998890432369458622168801317440E0_cfp
        w(9)=2.63703638151858236124550920981875174E-2_cfp
        x(10)=2.34521003143858924544943449453127818E0_cfp
        w(10)=8.25302204214347313662815866681467380E-3_cfp
        x(11)=2.72236326569406276810483548837218383E0_cfp
        w(11)=1.73458891460220161553366846867247856E-3_cfp
        x(12)=3.12007138598606330928003541692425606E0_cfp
        w(12)=2.35346858937890683159444095388262854E-4_cfp
        x(13)=3.54015705773118057606533036194218417E0_cfp
        w(13)=1.95386742765238714350691602807576062E-5_cfp
        x(14)=3.98606417260362170867298067718483602E0_cfp
        w(14)=9.20385065434960914678563186463837861E-7_cfp
        x(15)=4.46380357862327336127305250798099501E0_cfp
        w(15)=2.19744313630042844656879622820565791E-8_cfp
        x(16)=4.98416102038868502553107037163096837E0_cfp
        w(16)=2.21209650067245885354108615621711308E-10_cfp
        x(17)=5.56889625433947751475479086853450015E0_cfp
        w(17)=6.64342082704336917177622603512853062E-13_cfp
        x(18)=6.27611802050711333067389464978344722E0_cfp
        w(18)=2.51756748289203538984707382214850867E-16_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.2) then
        x(1)=1.39830776281895539207346772994635372E-1_cfp
        w(1)=4.96524059369746941590685321640428214E-5_cfp
        x(2)=2.91517154801946303139587008421406387E-1_cfp
        w(2)=1.13299820851504896343533596687429511E-3_cfp
        x(3)=4.80375215012040286364304863427899540E-1_cfp
        w(3)=8.71949062953273896750270171058388401E-3_cfp
        x(4)=7.02617699813400424746549993143255063E-1_cfp
        w(4)=3.53374155889078858953563992053835572E-2_cfp
        x(5)=9.54175840330874785838917084225570524E-1_cfp
        w(5)=8.83671717103415212392510016069602862E-2_cfp
        x(6)=1.23144185897714786792200437494762814E0_cfp
        w(6)=1.45912122360162008511180112143588900E-1_cfp
        x(7)=1.53148508324305159675255704652775485E0_cfp
        w(7)=1.63716709602761049879265989827980947E-1_cfp
        x(8)=1.85211523604885166794732758988439496E0_cfp
        w(8)=1.25867236374253090274130454924099594E-1_cfp
        x(9)=2.19188910545513490624414492812107034E0_cfp
        w(9)=6.60158541638975386761284598433547829E-2_cfp
        x(10)=2.55011123269857395362929013794292197E0_cfp
        w(10)=2.32882152629690711491813234020883799E-2_cfp
        x(11)=2.92686491530278752612211808986078928E0_cfp
        w(11)=5.39620109062547520909487680475306194E-3_cfp
        x(12)=3.32311023780885415745631983786424971E0_cfp
        w(12)=7.93180081667495451220204422504664259E-4_cfp
        x(13)=3.74090459657704715284183837211505453E0_cfp
        w(13)=7.03561109218502983930965669761217817E-5_cfp
        x(14)=4.18385676480146140147967378941692285E0_cfp
        w(14)=3.50200496756594740684915594677463356E-6_cfp
        x(15)=4.65807977278441997137849728100389946E0_cfp
        w(15)=8.75743034995545018345233201897895006E-8_cfp
        x(16)=5.17439420881499795520970671099213084E0_cfp
        w(16)=9.16955418698139349627415731603934786E-10_cfp
        x(17)=5.75449195515098567012157084046944392E0_cfp
        w(17)=2.84916339027130704927678552870742218E-12_cfp
        x(18)=6.45613676213529434601207520201523999E0_cfp
        w(18)=1.11343819444093965744569717678840604E-15_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.3) then
        x(1)=2.21042347831804064036320518232381608E-1_cfp
        w(1)=1.87609419327987681276349695328946994E-5_cfp
        x(2)=4.04193701670543745793589516907428835E-1_cfp
        w(2)=7.33564871855076838536543033962153189E-4_cfp
        x(3)=6.16981058332160358219285464701087495E-1_cfp
        w(3)=8.55973955428245142721291741557563812E-3_cfp
        x(4)=8.57335563732651242627645844516339083E-1_cfp
        w(4)=4.81989315332014451855534976054785023E-2_cfp
        x(5)=1.12227201921569280965738316087933105E0_cfp
        w(5)=1.56659724404253288374415944281306144E-1_cfp
        x(6)=1.40909356130027692636440991086671666E0_cfp
        w(6)=3.19212537935173358560041253528602716E-1_cfp
        x(7)=1.71564902445834074418493519758868781E0_cfp
        w(7)=4.24361848496319712679774492301188886E-1_cfp
        x(8)=2.04039712658213602491483587470992551E0_cfp
        w(8)=3.74371054844231675257896724323329369E-1_cfp
        x(9)=2.38241977947450243810423157044037132E0_cfp
        w(9)=2.19679672406939272309579722293110613E-1_cfp
        x(10)=2.74143668962930121488360134548162144E0_cfp
        w(10)=8.49769285355251988182208946626557946E-2_cfp
        x(11)=3.11785316134989576656191470347807762E0_cfp
        w(11)=2.12477564961835486571504469892567989E-2_cfp
        x(12)=3.51287431402958382783547883010099489E0_cfp
        w(12)=3.32717833134865085706578074723745530E-3_cfp
        x(13)=3.92873914838335006791675780695169743E0_cfp
        w(13)=3.11167075948141083957733300830860303E-4_cfp
        x(14)=4.36918423685037656817218582783504563E0_cfp
        w(14)=1.61946313115628473333630893021172538E-5_cfp
        x(15)=4.84040066503283612885304739145578601E0_cfp
        w(15)=4.20601500291350004761458941371027530E-7_cfp
        x(16)=5.35323134115378290975694805162742963E0_cfp
        w(16)=4.54936616662942787040459471318047931E-9_cfp
        x(17)=5.92930171342587487465762379872170168E0_cfp
        w(17)=1.45428034788464067314169991978149784E-11_cfp
        x(18)=6.62606355099873140541438478412559275E0_cfp
        w(18)=5.83246888602849600058345678336203556E-15_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.4) then
        x(1)=3.09471538715179342136033103154050383E-1_cfp
        w(1)=1.53458653374929850173173092640546892E-5_cfp
        x(2)=5.19517463604751889711400703171923575E-1_cfp
        w(2)=8.95278507904567343212247954008566991E-4_cfp
        x(3)=7.52187126765495241863409833107325824E-1_cfp
        w(3)=1.42120512952475269411542523280672770E-2_cfp
        x(4)=1.00738840777746244910824874775041893E0_cfp
        w(4)=1.02208524700270483282714559899668124E-1_cfp
        x(5)=1.28323079368934986392422578620691678E0_cfp
        w(5)=4.04497619737975202161836704670603298E-1_cfp
        x(6)=1.57783656997620409883611215749722719E0_cfp
        w(6)=9.66891203872603220857118006974446325E-1_cfp
        x(7)=1.88970880061291967156691621406573136E0_cfp
        w(7)=1.46424152879890959681799303572740593E0_cfp
        x(8)=2.21783094367793535295204131940000985E0_cfp
        w(8)=1.43752870128314616358809291006389463E0_cfp
        x(9)=2.56170170856233308071208250053171700E0_cfp
        w(9)=9.21380299845261504010795319167584861E-1_cfp
        x(10)=2.92136790317042818295355520389933935E0_cfp
        w(10)=3.83510755994203658276609717303272071E-1_cfp
        x(11)=3.29748830851791512730111615058381362E0_cfp
        w(11)=1.01945009748141840700391786652138105E-1_cfp
        x(12)=3.69146139048117626004400939931609426E0_cfp
        w(12)=1.68054011092531801924879312400658385E-2_cfp
        x(13)=4.10566984208105807727885824587710291E0_cfp
        w(13)=1.64143381676123143189776154629846962E-3_cfp
        x(14)=4.54395121779902132790442972448780936E0_cfp
        w(14)=8.86410622137526890854745407955806390E-5_cfp
        x(15)=5.01255704061073649518255271306032987E0_cfp
        w(15)=2.37614512235094845806633773055279765E-6_cfp
        x(16)=5.52234350492645748733936738343140766E0_cfp
        w(16)=2.64148343266786827818486556684649848E-8_cfp
        x(17)=6.09487148889392901153433061730503408E0_cfp
        w(17)=8.65028928187155457124692774024724292E-11_cfp
        x(18)=6.78730833083920894993749887101059934E0_cfp
        w(18)=3.54728478486354367640040422380968607E-14_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.5) then
        x(1)=4.02559198832693046295620595191141160E-1_cfp
        w(1)=2.17252054838271605225666825047729258E-5_cfp
        x(2)=6.35758658959463309525044659320851777E-1_cfp
        w(2)=1.73135307753104229871189432060905257E-3_cfp
        x(3)=8.85179697284295691158890107265576560E-1_cfp
        w(3)=3.48812249038488427030273834547604039E-2_cfp
        x(4)=1.15276228469813869074186616628736392E0_cfp
        w(4)=3.03301974473370496698609162664711404E-1_cfp
        x(5)=1.43766282947109014772836563278008782E0_cfp
        w(5)=1.39975917571440598480796349955722021E0_cfp
        x(6)=1.73873186938269578993492441687750372E0_cfp
        w(6)=3.79372708822801394570329146570094138E0_cfp
        x(7)=2.05502660100460404493770145583211827E0_cfp
        w(7)=6.37096779167225164360576507833994547E0_cfp
        x(8)=2.38596259445243322484818491842699321E0_cfp
        w(8)=6.81405163717928918803976053156485546E0_cfp
        x(9)=2.73137675224397987253426039192405215E0_cfp
        w(9)=4.69051031290060673653085563074706035E0_cfp
        x(10)=3.09157961314603964184855642887097985E0_cfp
        w(10)=2.07265739316353737691207066658545148E0_cfp
        x(11)=3.46743405235600865114956762199822179E0_cfp
        w(11)=5.79422826661811549333762132468936181E-1_cfp
        x(12)=3.86049438411670236565074016175388796E0_cfp
        w(12)=9.96814693993422510239633412384135378E-2_cfp
        x(13)=4.27325923217751028017186303669891504E0_cfp
        w(13)=1.00967415536390455547065893211175871E-2_cfp
        x(14)=4.70964731144690349049230544517512883E0_cfp
        w(14)=5.62510916215761167096163375511235204E-4_cfp
        x(15)=5.17595750045937237071249021484212002E0_cfp
        w(15)=1.54903431957645888406852293605405111E-5_cfp
        x(16)=5.68305272632472663410355630619419585E0_cfp
        w(16)=1.76294466637905058777301250742670715E-7_cfp
        x(17)=6.25243143639364445782030567578647230E0_cfp
        w(17)=5.89503520900470184971982827092992017E-10_cfp
        x(18)=6.94099991636645628977260340522825535E0_cfp
        w(18)=2.46475640028492684464193579996494132E-13_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.6) then
        x(1)=4.98546411434896622836481799360466948E-1_cfp
        w(1)=4.65264321280738789205483677113377000E-5_cfp
        x(2)=7.51824234905971771863831757981487510E-1_cfp
        w(2)=4.76487884330467804848989231499409552E-3_cfp
        x(3)=1.01553327516166491987625085565965052E0_cfp
        w(3)=1.16142816173783681481568657872378276E-1_cfp
        x(4)=1.29359512047128621005655370474506104E0_cfp
        w(4)=1.17556705578950302162884207368366232E0_cfp
        x(5)=1.58614018199074666543311199239528162E0_cfp
        w(5)=6.13850539934487902206826944539716910E0_cfp
        x(6)=1.89266372842049895437842777918730793E0_cfp
        w(6)=1.84152640613727774413585111940684885E1_cfp
        x(7)=2.21269654986486694506170454106047622E0_cfp
        w(7)=3.36416534971268941430930145703194435E1_cfp
        x(8)=2.54601562688668905126148918861748294E0_cfp
        w(8)=3.86002981637046062568773500866843875E1_cfp
        x(9)=2.89273738371050481804305003562514963E0_cfp
        w(9)=2.81855443607725426711698487152123851E1_cfp
        x(10)=3.25338914858160561058946596555098188E0_cfp
        w(10)=1.30908738460144337518388077538332369E1_cfp
        x(11)=3.62900129681848775455620860363405777E0_cfp
        w(11)=3.81772288164286903761159100430624142E0_cfp
        x(12)=4.02125606913624523062072434130867205E0_cfp
        w(12)=6.80927093394416463479555517833308272E-1_cfp
        x(13)=4.43274720581658295351180796900917766E0_cfp
        w(13)=7.11413627139417451721088276531099062E-2_cfp
        x(14)=4.86745962566255193264870312420023570E0_cfp
        w(14)=4.07088479262112340943691028614078324E-3_cfp
        x(15)=5.33172968466107915191317077071665320E0_cfp
        w(15)=1.14742020490390525188846675926289584E-4_cfp
        x(16)=5.83642225393084008368570491434584664E0_cfp
        w(16)=1.33284194208696173209213404635092921E-6_cfp
        x(17)=6.40297563127662440360281300056884830E0_cfp
        w(17)=4.53911628799967406986644222662870424E-9_cfp
        x(18)=7.08805514343297935309039481021615075E0_cfp
        w(18)=1.93062155686314764096066109595806282E-12_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.7) then
        x(1)=5.96207266641964871517582078605711308E-1_cfp
        w(1)=1.37997787644824499633794380459339630E-4_cfp
        x(2)=8.67023237135939796846754361081947937E-1_cfp
        w(2)=1.73717954160394412902947207589469307E-2_cfp
        x(3)=1.14305356623351161186529274050355105E0_cfp
        w(3)=4.94896011750160488836104897305023231E-1_cfp
        x(4)=1.43009960521982907440001443888054311E0_cfp
        w(4)=5.67300131462256592115073726456361404E0_cfp
        x(5)=1.72918534773021771725883161404705392E0_cfp
        w(5)=3.27884747049307128475193047399791114E1_cfp
        x(6)=2.04037683247552784212693030210469384E0_cfp
        w(6)=1.06973115662529862281783092675009218E2_cfp
        x(7)=2.36361315878201499317807386329143509E0_cfp
        w(7)=2.09584663139463733381744362624877768E2_cfp
        x(8)=2.69897767128574839106212733292505588E0_cfp
        w(8)=2.55032223267432455901819057270629049E2_cfp
        x(9)=3.04682146462536161082389992019685124E0_cfp
        w(9)=1.95705416038350975905565039165223704E2_cfp
        x(10)=3.40785301289563934695481655551146998E0_cfp
        w(10)=9.48161563847314703490127313321099358E1_cfp
        x(11)=3.78324240985502216031979500591813427E0_cfp
        w(11)=2.86674872993400470402419142644941641E1_cfp
        x(12)=4.17477855089679672113267929008187970E0_cfp
        w(12)=5.27411563631576951154330108573770194E0_cfp
        x(13)=4.58513435599228960965014363204530872E0_cfp
        w(13)=5.65979616305158692775897688491978437E-1_cfp
        x(14)=5.01834944845896666600079923913604569E0_cfp
        w(14)=3.31491483086272196050955373953465696E-2_cfp
        x(15)=5.48079008445081827732925802279969078E0_cfp
        w(15)=9.53564678576624389251768228616225834E-4_cfp
        x(16)=5.98331955673648172827012134239555544E0_cfp
        w(16)=1.12778775928527217693524866394213494E-5_cfp
        x(17)=6.54731833963892337525260901947064478E0_cfp
        w(17)=3.90359209407479385157307120488388120E-8_cfp
        x(18)=7.22922825300888738775497031691319761E0_cfp
        w(18)=1.68590127930001045732246161842619501E-11_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.8) then
        x(1)=6.94679979790960945676494559071147783E-1_cfp
        w(1)=5.33179638882069112031587122677201518E-4_cfp
        x(2)=9.80923345705274617547368399968106001E-1_cfp
        w(2)=7.97823213960694223560506996493518267E-2_cfp
        x(3)=1.26768389554706711966731135612683399E0_cfp
        w(3)=2.58868890561050726209382460144162902E0_cfp
        x(4)=1.56251946736982263551303173410533397E0_cfp
        w(4)=3.29228644524592689741332069436983355E1_cfp
        x(5)=1.86726840021948572171621524510529230E0_cfp
        w(5)=2.07171017698505066619726760837743273E2_cfp
        x(6)=2.18250317792039300073409541464279866E0_cfp
        w(6)=7.25317608879563872179966855636599297E2_cfp
        x(7)=2.50851786083585504061488171912000865E0_cfp
        w(7)=1.50764547438731855427310002995001071E3_cfp
        x(8)=2.84565856039101523388536167261287196E0_cfp
        w(8)=1.92854386939557011026312776654443964E3_cfp
        x(9)=3.19447604836858385063041579162643624E0_cfp
        w(9)=1.54411591819749076284113232106668780E3_cfp
        x(10)=3.55583229777799986614151998839630307E0_cfp
        w(10)=7.75761357670629663809178904137012202E2_cfp
        x(11)=3.93101558285051816533014137252520887E0_cfp
        w(11)=2.41986192163478254823865243370554984E2_cfp
        x(12)=4.32190492129501687199597037071726775E0_cfp
        w(12)=4.57365930915988559868701484903330295E1_cfp
        x(13)=4.73123998164370931956378571972754192E0_cfp
        w(13)=5.02447431248022414716886661229630944E0_cfp
        x(14)=5.16310610041642854984497739647618828E0_cfp
        w(14)=3.00368486270802426284366291961394289E-1_cfp
        x(15)=5.62389353760375240712947531066407849E0_cfp
        w(15)=8.79749825996402157673156628945389671E-3_cfp
        x(16)=6.12446140385537934255137943317329053E0_cfp
        w(16)=1.05729957392887994121581731394757287E-4_cfp
        x(17)=6.68613465602336547160743607050437706E0_cfp
        w(17)=3.71315455995752790791657762747349963E-7_cfp
        x(18)=7.36514689943095796438159403248851997E0_cfp
        w(18)=1.62591281848930305584844512815781628E-10_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.9) then
        x(1)=7.93355634773069815677984124048687220E-1_cfp
        w(1)=2.56706046255139890852241607045955363E-3_cfp
        x(2)=1.09326092240460102643993230480559953E0_cfp
        w(2)=4.44895661726104034808471671790464299E-1_cfp
        x(3)=1.38944888678260904819612478259579025E0_cfp
        w(3)=1.61178965198470667808778339079882238E1_cfp
        x(4)=1.69110532600263790166456543595347936E0_cfp
        w(4)=2.23865510129128101524958512764246520E2_cfp
        x(5)=2.00080881210550641975282013304918211E0_cfp
        w(5)=1.51431283817047890112954248011338923E3_cfp
        x(6)=2.31958277621661274120899214657166305E0_cfp
        w(6)=5.63085983140256361873523294029717183E3_cfp
        x(7)=2.64803227583117011051988123870469783E0_cfp
        w(7)=1.23132030208812749572820040002286780E4_cfp
        x(8)=2.98673113691212794201512282239723174E0_cfp
        w(8)=1.64435342989472551162903396537974491E4_cfp
        x(9)=3.33640202050077802203281863780538258E0_cfp
        w(9)=1.36589739492825853802094147556204999E4_cfp
        x(10)=3.69803855572132688870944438066546650E0_cfp
        w(10)=7.08262018839077924290209686271509579E3_cfp
        x(11)=4.07303019271583834296045416646288952E0_cfp
        w(11)=2.27048410136234135030452652567041180E3_cfp
        x(12)=4.46333299591960297142979696097897376E0_cfp
        w(12)=4.39435532255929148990494155724850315E2_cfp
        x(13)=4.87174352546312182891847357442278667E0_cfp
        w(13)=4.92855194602033692949926404645591247E1_cfp
        x(14)=5.30238569450684142499755784979893373E0_cfp
        w(14)=3.00044843057328841251737882506172797E0_cfp
        x(15)=5.76166911898272321353528607140229181E0_cfp
        w(15)=8.93066061434013291136477165713870509E-2_cfp
        x(16)=6.26044679758171621712368033894063623E0_cfp
        w(16)=1.08886698763659762830465564163028869E-3_cfp
        x(17)=6.81999041483146900501549779387091276E0_cfp
        w(17)=3.87450636108294075822423987584411557E-6_cfp
        x(18)=7.49633886443412713992102358127692080E0_cfp
        w(18)=1.71794762455863881220447320150205756E-9_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.10) then
        x(1)=8.91803144526040524849388438116912588E-1_cfp
        w(1)=1.48983976820956718869263629752327832E-2_cfp
        x(2)=1.20388346452342837096152843316411263E0_cfp
        w(2)=2.92981592176562792251421482192475048E0_cfp
        x(3)=1.50842019888827131089939581384870953E0_cfp
        w(3)=1.16678414955841321857211034068940868E2_cfp
        x(4)=1.81610177270244780128665073352136199E0_cfp
        w(4)=1.74801507151959463855550988245505335E3_cfp
        x(5)=2.13017949972116694293978727954321682E0_cfp
        w(5)=1.25840443114551281519531566154988624E4_cfp
        x(6)=2.45207999441698163547592840617317694E0_cfp
        w(6)=4.92938948898244631001197470691507978E4_cfp
        x(7)=2.78268277256283614821754114841314385E0_cfp
        w(7)=1.12643524195334957093995761174154264E5_cfp
        x(8)=3.12276081401334611509334280350605876E0_cfp
        w(8)=1.56180437237952946336598733268326576E5_cfp
        x(9)=3.47318626027830816826554691265476870E0_cfp
        w(9)=1.33979101299277829124703117925882545E5_cfp
        x(10)=3.83506685995555482339931273767085860E0_cfp
        w(10)=7.14316189298273364689261991003917653E4_cfp
        x(11)=4.20987959784466113575090450049743159E0_cfp
        w(11)=2.34584972352997529441300712116305672E4_cfp
        x(12)=4.59964708313748912238597030579460843E0_cfp
        w(12)=4.63689417785487931175728453145266341E3_cfp
        x(13)=5.00721483735540282709992121259114918E0_cfp
        w(13)=5.29760092412390117256822496901935576E2_cfp
        x(14)=5.43673960954410053721691919207138724E0_cfp
        w(14)=3.27816512598949242803942031155622642E1_cfp
        x(15)=5.89464666519635642917041742742789168E0_cfp
        w(15)=9.89980437503727458084656766564980702E-1_cfp
        x(16)=6.39178146104637574807135453497878383E0_cfp
        w(16)=1.22286057726028906715675333648569097E-2_cfp
        x(17)=6.94936456946136957346487107504005668E0_cfp
        w(17)=4.40352959026484079519758080865823672E-5_cfp
        x(18)=7.62325217607052814622257903801008187E0_cfp
        w(18)=1.97501765615089639560664865182373034E-8_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.11) then
        x(1)=9.89717799948204380578851301462278199E-1_cfp
        w(1)=1.01601248521203032918405442365057620E-1_cfp
        x(2)=1.31271216932867591939349827448433335E0_cfp
        w(2)=2.23003827737115322913121863922755517E1_cfp
        x(3)=1.62469550564953544167638853449098377E0_cfp
        w(3)=9.64141518295506542979522583638253421E2_cfp
        x(4)=1.93774085030386945590574129093730707E0_cfp
        w(4)=1.54263865916751819129808902285769796E4_cfp
        x(5)=2.25571166097383147277790166040021500E0_cfp
        w(5)=1.17246232676239353846730448560087698E5_cfp
        x(6)=2.58039659955890933554250409918985329E0_cfp
        w(6)=4.80654406912370963950331175711199462E5_cfp
        x(7)=2.91291901960953170336666877902423127E0_cfp
        w(7)=1.14159645624866116192992709767308891E6_cfp
        x(8)=3.25422750401078380593529402728316831E0_cfp
        w(8)=1.63599674519030674332940331091461733E6_cfp
        x(9)=3.60532535678082494342301307897022215E0_cfp
        w(9)=1.44397652446414382085491247152057540E6_cfp
        x(10)=3.96742016631719961017294517003696503E0_cfp
        w(10)=7.89109193590033624548370530568351172E5_cfp
        x(11)=4.34206535402448445724059713083662187E0_cfp
        w(11)=2.64786661614209295876571779129171079E5_cfp
        x(12)=4.73134152335398135763195522942354208E0_cfp
        w(12)=5.33355792557040774065807869350605366E4_cfp
        x(13)=5.13813669421571827583047156031825014E0_cfp
        w(13)=6.19567664729917673451038576696656040E3_cfp
        x(14)=5.56663577799913173208945653312896577E0_cfp
        w(14)=3.89080395970709814477431679696500962E2_cfp
        x(15)=6.02327670349868425879385207304948708E0_cfp
        w(15)=1.19055319767333484476789386418907289E1_cfp
        x(16)=6.51889632995562666850277482362361708E0_cfp
        w(16)=1.48818187573197086993005379784467770E-1_cfp
        x(17)=7.07466617981951614264459761956045497E0_cfp
        w(17)=5.41782502250855911852814231993656005E-4_cfp
        x(18)=7.74627046568063711030350230353596780E0_cfp
        w(18)=2.45570227023723681985070963439523842E-7_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.12) then
        x(1)=1.08688551983051780406795953648536361E0_cfp
        w(1)=7.97995745880380053298379394688088005E-1_cfp
        x(2)=1.41971723292536716546058552582006216E0_cfp
        w(2)=1.92875862830589987608195816394956962E2_cfp
        x(3)=1.73838554197713314808773869625240300E0_cfp
        w(3)=8.96147077021097131819243634584984816E3_cfp
        x(4)=2.05623915416397451116633348336946162E0_cfp
        w(4)=1.51905097208725633606196610140976195E5_cfp
        x(5)=2.37769965653726980754549489160944435E0_cfp
        w(5)=1.21095507087312941318845174315305415E6_cfp
        x(6)=2.70488224255496640962066447052331455E0_cfp
        w(6)=5.16765442668126444095350414945650729E6_cfp
        x(7)=3.03912823782714938536689518568968129E0_cfp
        w(7)=1.27002629993028913339438730549689883E7_cfp
        x(8)=3.38154219857358372770453346997693988E0_cfp
        w(8)=1.87425456744102362384803466100258868E7_cfp
        x(9)=3.73324343588202287819899211270714223E0_cfp
        w(9)=1.69680886663164140263236117418473619E7_cfp
        x(10)=4.09552760405116021423770558282351308E0_cfp
        w(10)=9.47998663621067592287975249655014514E6_cfp
        x(11)=4.47001541890602176333781001896007858E0_cfp
        w(11)=3.24315064971651599791006833995982868E6_cfp
        x(12)=4.85883843174432552134891597082753809E0_cfp
        w(12)=6.64475758289323601304448147497580685E5_cfp
        x(13)=5.26492183253348993259425874763681287E0_cfp
        w(13)=7.83589955141945342862473829161453748E4_cfp
        x(14)=5.69247484924329487090792693564130997E0_cfp
        w(14)=4.98720660221708850261496932064652795E3_cfp
        x(15)=6.14794564563810095899847093541949759E0_cfp
        w(15)=1.54448129854600898140704899486047993E2_cfp
        x(16)=6.64216171099825634795814784922729787E0_cfp
        w(16)=1.95171821308813889058019365607357703E0_cfp
        x(17)=7.19624747655279801085898509020078983E0_cfp
        w(17)=7.17720019691493310137498374292536309E-3_cfp
        x(18)=7.86572483448295023063513372140929763E0_cfp
        w(18)=3.28506880798168458249681507965100792E-6_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.13) then
        x(1)=1.18315771914641179785226912580633907E0_cfp
        w(1)=7.10332774424439533946070438053209172E0_cfp
        x(2)=1.52490135404964879450750680413586519E0_cfp
        w(2)=1.86968881366043144869539637417323621E3_cfp
        x(3)=1.84960612848650989241216598619545944E0_cfp
        w(3)=9.25796915813314648170073504214471996E4_cfp
        x(4)=2.17179696058613754765826417816369458E0_cfp
        w(4)=1.65156882694468260287507411471260520E6_cfp
        x(5)=2.49640557645019781706359054460207396E0_cfp
        w(5)=1.37350804304606721831292236129984226E7_cfp
        x(6)=2.82584292982539289786260829123375421E0_cfp
        w(6)=6.07435390750013858241963620328850252E7_cfp
        x(7)=3.16164630353158886917843009849876867E0_cfp
        w(7)=1.53907030813924638117935831184984764E8_cfp
        x(8)=3.50505970595200460338009440831540529E0_cfp
        w(8)=2.33175744694352005102419574745550634E8_cfp
        x(9)=3.85730578065823930170361327782202979E0_cfp
        w(9)=2.15969162643022731775806366743958164E8_cfp
        x(10)=4.21975843128053591305705425540667252E0_cfp
        w(10)=1.23088536811101327929736439113582808E8_cfp
        x(11)=4.59409805293187312057968210661295189E0_cfp
        w(11)=4.28524655383328033226023394319259900E7_cfp
        x(12)=4.98250127469258969368733491934255639E0_cfp
        w(12)=8.91656897095970860103800692699668652E6_cfp
        x(13)=5.38792601806462021138837892968231349E0_cfp
        w(13)=1.06601801019557493229822095263220869E6_cfp
        x(14)=5.81460263230535863997431692347486263E0_cfp
        w(14)=6.86836741309861551177032082071081029E4_cfp
        x(15)=6.26898752518943700897553517810217459E0_cfp
        w(15)=2.15063474788949180697339482176735803E3_cfp
        x(16)=6.76189825924512492661474911439551313E0_cfp
        w(16)=2.74509508579096161563436096004634987E1_cfp
        x(17)=7.31441402869316805940605787056905928E0_cfp
        w(17)=1.01891546004443769661837458484805278E-1_cfp
        x(18)=7.98190312727071118056273614605197205E0_cfp
        w(18)=4.70617541283569166640464602224870210E-5_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.14) then
        x(1)=1.27843345394585583592763342164570167E0_cfp
        w(1)=7.07285438629537333382571609406453953E1_cfp
        x(2)=1.62828860621501995581299276024448928E0_cfp
        w(2)=2.00859518361385222722376347157884942E4_cfp
        x(3)=1.95847330222687041070655510443731981E0_cfp
        w(3)=1.05259334005859944443787031884647944E6_cfp
        x(4)=2.28459846193505794443977706794576253E0_cfp
        w(4)=1.96529462678839804904373419419454379E7_cfp
        x(5)=2.61206334747972539161860001082870128E0_cfp
        w(5)=1.69740501815728145774403784602819740E8_cfp
        x(6)=2.94354790427561550502819291213877799E0_cfp
        w(6)=7.75074934092412751366018249335311408E8_cfp
        x(7)=3.28076650129193882099300494699847363E0_cfp
        w(7)=2.01835154063760872745255942916348574E9_cfp
        x(8)=3.62508857096955034739186739940299792E0_cfp
        w(8)=3.13114799069158928515215864158897448E9_cfp
        x(9)=3.97782938656912464718530234038855965E0_cfp
        w(9)=2.96050550036214632900220811419555799E9_cfp
        x(10)=4.34043283433191014118724842500355715E0_cfp
        w(10)=1.71805997042419588880573162451954295E9_cfp
        x(11)=4.71463258203063655076260467058931215E0_cfp
        w(11)=6.07731082079781491957349204597007774E8_cfp
        x(12)=5.10264539868372510779132879684107120E0_cfp
        w(12)=1.28250727539287544648783036100834174E8_cfp
        x(13)=5.50745820542601282603707009198426445E0_cfp
        w(13)=1.55268719991724888262933327306376780E7_cfp
        x(14)=5.93331979472586297744866280393687712E0_cfp
        w(14)=1.01172815864489272531258218656768471E6_cfp
        x(15)=6.38669317396651390120283901955649687E0_cfp
        w(15)=3.20032197583085428683151092649444342E4_cfp
        x(16)=6.87838558757702101502580904475267128E0_cfp
        w(16)=4.12304076440515744367903230609529508E2_cfp
        x(17)=7.42943274448259261403744672518364746E0_cfp
        w(17)=1.54367527081472015297006450415594990E0_cfp
        x(18)=8.09505725622077487612657836713415414E0_cfp
        w(18)=7.19058950422824847740100649019217927E-4_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.15) then
        x(1)=1.37264661215145504061322589730860928E0_cfp
        w(1)=7.79287247681844982377287569039727677E2_cfp
        x(2)=1.72991687076230457567145383015447561E0_cfp
        w(2)=2.36905725777899611556198864649944972E5_cfp
        x(3)=2.06510040182350858197091133721507918E0_cfp
        w(3)=1.30623065856894724056990403501707934E7_cfp
        x(4)=2.39481257666402050148100539257961328E0_cfp
        w(4)=2.54067983770664408610007214856926189E8_cfp
        x(5)=2.72488234722707746746226804049135401E0_cfp
        w(5)=2.27028823316658977264267053710128765E9_cfp
        x(6)=3.05823526360853833774557990476126947E0_cfp
        w(6)=1.06699680319018107750547059044329563E10_cfp
        x(7)=3.39674649641793836589264059721603780E0_cfp
        w(7)=2.84819872618528178992031678497430865E10_cfp
        x(8)=3.74189889481024255201440488194885700E0_cfp
        w(8)=4.51439237155272889636746493903229782E10_cfp
        x(9)=4.09509124476698834651696445563891890E0_cfp
        w(9)=4.34914314087190413426562916206456158E10_cfp
        x(10)=4.45783039072969036245815801294959006E0_cfp
        w(10)=2.56585935453514174323914310994079962E10_cfp
        x(11)=4.83189783465086452706750234269429187E0_cfp
        w(11)=9.20941918206178183259076473806368355E9_cfp
        x(12)=5.21954629523617674365213639293316126E0_cfp
        w(12)=1.96880588947078588328619797515659075E9_cfp
        x(13)=5.62378852893507425610665705287318715E0_cfp
        w(13)=2.41129880350868518410788562201025504E8_cfp
        x(14)=6.04888950905524756018659897086938881E0_cfp
        w(14)=1.58762446526194279693801182297562722E7_cfp
        x(15)=6.50131747548317093790427360149556700E0_cfp
        w(15)=5.06955627509398733019359957095272137E5_cfp
        x(16)=6.99186909082616791868189854778631558E0_cfp
        w(16)=6.58784965568888974445074571186067092E3_cfp
        x(17)=7.54153823166938638036154551161960984E0_cfp
        w(17)=2.48649874467079926013751503572223282E1_cfp
        x(18)=8.20540904166758442295324150932731590E0_cfp
        w(18)=1.16746415260404573115218647955973447E-2_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.16) then
        x(1)=1.46575664193453841834203160410697699E0_cfp
        w(1)=9.41529863367300755426078317838258415E3_cfp
        x(2)=1.82983266100276046726727238096528018E0_cfp
        w(2)=3.04360732586684399612709119442772277E6_cfp
        x(3)=2.16959638900322797367364426487935031E0_cfp
        w(3)=1.75689318035142960012874933787909644E8_cfp
        x(4)=2.50259402801648709305900340996896192E0_cfp
        w(4)=3.54583219650210227597266354730819579E9_cfp
        x(5)=2.83505054490680979296609970673167269E0_cfp
        w(5)=3.26751119370185784587343216064156627E10_cfp
        x(6)=3.17011657230397601100819471374220757E0_cfp
        w(6)=1.57636712042509739821115733804874213E11_cfp
        x(7)=3.50981394198127989013519926876478387E0_cfp
        w(7)=4.30369941540434036511990738428281169E11_cfp
        x(8)=3.85572856638077616609226848059914531E0_cfp
        w(8)=6.95615054180242823983381517390136532E11_cfp
        x(9)=4.20933491639282864038301222871391134E0_cfp
        w(9)=6.81735592353311438894317833536256389E11_cfp
        x(10)=4.57219677712974713729596485957022542E0_cfp
        w(10)=4.08324373364499905485773031649822793E11_cfp
        x(11)=4.94613883045445925078197070073099841E0_cfp
        w(11)=1.48532060704177988021381386537279922E11_cfp
        x(12)=5.33344616054572296156238160188958692E0_cfp
        w(12)=3.21346931710978856853802505420773156E10_cfp
        x(13)=5.73715465550608671384091217077509919E0_cfp
        w(13)=3.97802588373588374506408609919434451E9_cfp
        x(14)=6.16154354452148824962615473053084985E0_cfp
        w(14)=2.64455907679492186193141956562070123E8_cfp
        x(15)=6.61308515674246272271593946644555017E0_cfp
        w(15)=8.51887584380273028221130618234103179E6_cfp
        x(16)=7.10256540812010770634422279652505302E0_cfp
        w(16)=1.11597740680181930269767229082771118E5_cfp
        x(17)=7.65093790182775088314967356751913429E0_cfp
        w(17)=4.24407779280229774908443764628563909E2_cfp
        x(18)=8.31315491293884172449897759650393714E0_cfp
        w(18)=2.00760797477369823004155979998779070E-1_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.17) then
        x(1)=1.55774178583023084377574987355025662E0_cfp
        w(1)=1.23785811280532477841841003404642722E5_cfp
        x(2)=1.92808756948739143559888657173836472E0_cfp
        w(2)=4.23067226198706307357246226884654153E7_cfp
        x(3)=2.27206495362457697678855263779635126E0_cfp
        w(3)=2.54579850306779248163632638489319860E9_cfp
        x(4)=2.60808451692747837232826353317254000E0_cfp
        w(4)=5.31330022587970762742037945736972410E10_cfp
        x(5)=2.94273721223241192350900688122235669E0_cfp
        w(5)=5.03528368487617823697844369697526064E11_cfp
        x(6)=3.27938066812800149580896283141213627E0_cfp
        w(6)=2.48781920165626769971513662361924352E12_cfp
        x(7)=3.62017102700004937897363365487639458E0_cfp
        w(7)=6.93326415963842612863242226046831268E12_cfp
        x(8)=3.96678827695140530292263724238215334E0_cfp
        w(8)=1.14090453958203534503803809754521445E13_cfp
        x(9)=4.32077580423270690038609539961759504E0_cfp
        w(9)=1.13587742950927571453618423751709835E13_cfp
        x(10)=4.68374914141582385781678439946140252E0_cfp
        w(10)=6.89860642491628876444768770285175992E12_cfp
        x(11)=5.05757213763405713933130587342155624E0_cfp
        w(11)=2.54065331569244933969698046805840851E12_cfp
        x(12)=5.44455915445866169203739074650824252E0_cfp
        w(12)=5.55773387870467488415097234390766722E11_cfp
        x(13)=5.84776688541400496080141191308287677E0_cfp
        w(13)=6.94869699010671415133281355050203240E10_cfp
        x(14)=6.27148716693452770535271782621361411E0_cfp
        w(14)=4.66110042448852881414525575723526739E9_cfp
        x(15)=6.72219545640865109219545744750419777E0_cfp
        w(15)=1.51381034636829234354457131677406834E8_cfp
        x(16)=7.21066683516569275202460707627429405E0_cfp
        w(16)=1.99811923214295872800765974880983093E6_cfp
        x(17)=7.75781610308412900941171435632525185E0_cfp
        w(17)=7.65303770094612391331112596714584966E3_cfp
        x(18)=8.41846972425051619459004466234786610E0_cfp
        w(18)=3.64573161484917211951671851715905968E0_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.18) then
        x(1)=1.64859410695647857422983304400808786E0_cfp
        w(1)=1.75936502599057489504881684120033228E6_cfp
        x(2)=2.02473582763517093404620733611811342E0_cfp
        w(2)=6.32591481346530781204490608930847504E8_cfp
        x(3)=2.37260411375971597357381531329097327E0_cfp
        w(3)=3.95361404265599993932928356210201773E10_cfp
        x(4)=2.71141389197952770734221700175565625E0_cfp
        w(4)=8.50799526012119658466785749931745396E11_cfp
        x(5)=3.04809525545566513291702371009112749E0_cfp
        w(5)=8.27185045201900865810548457996867790E12_cfp
        x(6)=3.38619682097384263402281606058306834E0_cfp
        w(6)=4.17716526443243993787393707982887756E13_cfp
        x(7)=3.72799819542034481877011689608455142E0_cfp
        w(7)=1.18632017195928920348810260889996218E14_cfp
        x(8)=4.07526559212127785043932617674596164E0_cfp
        w(8)=1.98461323503012706135600718273341334E14_cfp
        x(9)=4.42960541964806719460934791792956239E0_cfp
        w(9)=2.00475848664145202638198563952723047E14_cfp
        x(10)=4.79268044589838071694273347963006249E0_cfp
        w(10)=1.23332290822058907679490854013008480E14_cfp
        x(11)=5.16639020443945247959000312659546689E0_cfp
        w(11)=4.59448996389014185475975213958742255E13_cfp
        x(12)=5.55307565597794615610472572747301234E0_cfp
        w(12)=1.01542627382346045919741800895227467E13_cfp
        x(13)=5.95581228516100485682785894873984661E0_cfp
        w(13)=1.28136046779253072520886001326374680E12_cfp
        x(14)=6.37890311522215736576780600815364813E0_cfp
        w(14)=8.66759790909421512652371451960210935E10_cfp
        x(15)=6.82882592018550598345550973074746093E0_cfp
        w(15)=2.83668098994254211501283510635449975E9_cfp
        x(16)=7.31634491866117176404009833612343203E0_cfp
        w(16)=3.77084172115233183888992040144299851E7_cfp
        x(17)=7.86233749391665671517442025841792973E0_cfp
        w(17)=1.45397485507688153634898779282594331E5_cfp
        x(18)=8.52150987718650211741824965353456313E0_cfp
        w(18)=6.97262612607715304047062575032052130E1_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.19) then
        x(1)=1.73831580863555587053264752950853715E0_cfp
        w(1)=2.68798503575077674969982294557582174E7_cfp
        x(2)=2.11983263394323055987874407680016813E0_cfp
        w(2)=1.01238883548835035183127565318186182E10_cfp
        x(3)=2.47130612566491311947351555943653687E0_cfp
        w(3)=6.55061052958007045423308878916348523E11_cfp
        x(4)=2.81270126470292497027093331678247578E0_cfp
        w(4)=1.44976622848318565572514423731580771E13_cfp
        x(5)=3.15126321929782072825630983457141978E0_cfp
        w(5)=1.44304091532518972527868448751521435E14_cfp
        x(6)=3.49071736838418542081047887037362954E0_cfp
        w(6)=7.43500979394481934124293557678514032E14_cfp
        x(7)=3.83345720983866998946941243446796659E0_cfp
        w(7)=2.14862907632167686321259464710322338E15_cfp
        x(8)=4.18132828583716222441380103125890491E0_cfp
        w(8)=3.64965535036717665079538976944487878E15_cfp
        x(9)=4.53599486616856885237947333795650399E0_cfp
        w(9)=3.73658312796767935802885053827619928E15_cfp
        x(10)=4.89916300939780047556311772281189437E0_cfp
        w(10)=2.32633406238698840277753270797265728E15_cfp
        x(11)=5.27276489184307720602605923043442158E0_cfp
        w(11)=8.75914198196007745565722372247730872E14_cfp
        x(12)=5.65916573641559815938931879123565293E0_cfp
        w(12)=1.95447430457594194586250429973886124E14_cfp
        x(13)=6.06145806444732722424750998233470162E0_cfp
        w(13)=2.48775553154572102794730659041168642E13_cfp
        x(14)=6.48395485537422385601809043749929956E0_cfp
        w(14)=1.69608638264828178054969347225909902E12_cfp
        x(15)=6.93313551163387087653760070892068684E0_cfp
        w(15)=5.59096410250073217238358755223391542E10_cfp
        x(16)=7.41975340767123887434602484663799630E0_cfp
        w(16)=7.48193004259346421816241862924063533E8_cfp
        x(17)=7.96464981872908334817000543336034023E0_cfp
        w(17)=2.90320405006134808842772908085945276E6_cfp
        x(18)=8.62241589504723133905399177979841171E0_cfp
        w(18)=1.40105404692725687491259210797619389E3_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.20) then
        x(1)=1.82691649546007614023760020885250770E0_cfp
        w(1)=4.39278048950106501274896349793303171E8_cfp
        x(2)=2.21343301686871364574794049707826098E0_cfp
        w(2)=1.72651595675611793172579568288277420E11_cfp
        x(3)=2.56825758411884518243064011742798318E0_cfp
        w(3)=1.15331964886124947765707667368038143E13_cfp
        x(4)=2.91205604510540811060016281976517366E0_cfp
        w(4)=2.61925358201389475390479721423539168E14_cfp
        x(5)=3.25236700956485112744285471241882288E0_cfp
        w(5)=2.66418744837271571818186696408022966E15_cfp
        x(6)=3.59307992609752922212715259471998387E0_cfp
        w(6)=1.39837472594879999856250485485884473E16_cfp
        x(7)=3.93669369307449638227125627681971525E0_cfp
        w(7)=4.10674731191389366806709347070077529E16_cfp
        x(8)=4.28512709122875362266918883700505694E0_cfp
        w(8)=7.07494099919537370951532535415888857E16_cfp
        x(9)=4.64009770629398855420261504165997246E0_cfp
        w(9)=7.33446678765074382981743520250877967E16_cfp
        x(10)=5.00335141932926004829754106249524825E0_cfp
        w(10)=4.61733944248671700028999841133393301E16_cfp
        x(11)=5.37685037790087615746295262086898341E0_cfp
        w(11)=1.75590521116981891322077520164820837E16_cfp
        x(12)=5.76298201663145553984720083511745935E0_cfp
        w(12)=3.95328222162677779964461822448904986E15_cfp
        x(13)=6.16485435717545578933701157024162355E0_cfp
        w(13)=5.07289074063911748171972386267629899E14_cfp
        x(14)=6.58678926364567158555615499858994267E0_cfp
        w(14)=3.48419669430611386737070298690509109E13_cfp
        x(15)=7.03526718119195468913442280726324505E0_cfp
        w(15)=1.15634345547250392077065490008149595E12_cfp
        x(16)=7.52103069496137818415936308266443130E0_cfp
        w(16)=1.55722574782417127688967172947316300E10_cfp
        x(17)=8.06488620782728935866119118481132018E0_cfp
        w(17)=6.07876257241579187785265415885458405E7_cfp
        x(18)=8.72131456031102382548463691447454303E0_cfp
        w(18)=2.95116355820101608419150385731840885E4_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.21) then
        x(1)=1.91441112506585041327652621296581681E0_cfp
        w(1)=7.64562546583423114102567349414598520E9_cfp
        x(2)=2.30559107165006188430068704122456252E0_cfp
        w(2)=3.12541993735253927510776436508829595E12_cfp
        x(3)=2.66353963575148215440486152939786193E0_cfp
        w(3)=2.15010339030362988604048398067067333E14_cfp
        x(4)=3.00957888771537397785418383938764709E0_cfp
        w(4)=5.00084072074911560889006549737875662E15_cfp
        x(5)=3.35152137591635791022024013346767624E0_cfp
        w(5)=5.18958081892649586695331085521321661E16_cfp
        x(6)=3.69340925171073427599408286476681599E0_cfp
        w(6)=2.77112788981367995420800945059392131E17_cfp
        x(7)=4.03783925040621443729364514149485002E0_cfp
        w(7)=8.26085485474067363057686387355379114E17_cfp
        x(8)=4.38679798651440337272685604916633562E0_cfp
        w(8)=1.44197062216053965605039989074294257E18_cfp
        x(9)=4.74205233819928479372443197628889879E0_cfp
        w(9)=1.51236622871323724600927511178275925E18_cfp
        x(10)=5.10538494378273552878780572881074122E0_cfp
        w(10)=9.62028028667489268520855392127745154E17_cfp
        x(11)=5.47878556341656131959492033974472185E0_cfp
        w(11)=3.69266798116552988191864411588535736E17_cfp
        x(12)=5.86466203500228230507941849159709844E0_cfp
        w(12)=8.38386389390851291136328912914788913E16_cfp
        x(13)=6.26613652840775683935100260209823830E0_cfp
        w(13)=1.08405354241938170903336205340814656E16_cfp
        x(14)=6.68753885504150983725430581694852930E0_cfp
        w(14)=7.49754600877542199402945743622991270E14_cfp
        x(15)=7.13535000274996438894026355178397897E0_cfp
        w(15)=2.50429024357282178988388109600130642E13_cfp
        x(16)=7.62030185043171498883701369861695445E0_cfp
        w(16)=3.39266447192097504525820951986455069E11_cfp
        x(17)=8.16316709623967627640773535645615480E0_cfp
        w(17)=1.33190042162443173887209729639514974E9_cfp
        x(18)=8.81832070114515585664826740831922151E0_cfp
        w(18)=6.50317620675318035863776916524571651E5_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.22) then
        x(1)=2.00081847038397631166531208272542330E0_cfp
        w(1)=1.41184157873625304356641248009201418E11_cfp
        x(2)=2.39635945963594176587054777262338896E0_cfp
        w(2)=5.98501781948215989495160069000292569E13_cfp
        x(3)=2.75722825528036277497601830490616544E0_cfp
        w(3)=4.23096086439603667355165620984190144E15_cfp
        x(4)=3.10536254714400216926907131869703722E0_cfp
        w(4)=1.00605237201768464230238687556766341E17_cfp
        x(5)=3.44883119070823808444808663039309741E0_cfp
        w(5)=1.06362137987329382478176835324245992E18_cfp
        x(6)=3.79181882374444000923747132640145098E0_cfp
        w(6)=5.77101686692524628686150301002271554E18_cfp
        x(7)=4.13701325256908545283160080722725649E0_cfp
        w(7)=1.74449298249548203814459730435283606E19_cfp
        x(8)=4.48646410722798158024142999563653086E0_cfp
        w(8)=3.08265931321394983831830531610775076E19_cfp
        x(9)=4.84198397971195431566660631821534536E0_cfp
        w(9)=3.26852486268626071075843633822184500E19_cfp
        x(10)=5.20538954337208563086003354681966688E0_cfp
        w(10)=2.09944617599145048143508877637046758E19_cfp
        x(11)=5.57869607840802524769397834542175274E0_cfp
        w(11)=8.12928094649690305832121967914402005E18_cfp
        x(12)=5.96433022345240496070848529225318130E0_cfp
        w(12)=1.86031341103288876255833801672925281E18_cfp
        x(13)=6.36542710113381080008164974506673686E0_cfp
        w(13)=2.42276241954206314983579423226248636E17_cfp
        x(14)=6.78632364657930430558073147925849280E0_cfp
        w(14)=1.68667901895188391866929942245687102E16_cfp
        x(15)=7.23350096231612685429518716753948171E0_cfp
        w(15)=5.66800671513533041222407936571321607E14_cfp
        x(16)=7.71768032580014564148646264948867436E0_cfp
        w(16)=7.72227022843801781220309829190962321E12_cfp
        x(17)=8.25960183474975386569029067710101569E0_cfp
        w(17)=3.04805057267319501134676313707922934E10_cfp
        x(18)=8.91353869391352499777152393201501128E0_cfp
        w(18)=1.49636312124452455381688470588362809E7_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.23) then
        x(1)=2.08615996172782600500258007788395278E0_cfp
        w(1)=2.75664940710350376957430245938608868E12_cfp
        x(2)=2.48578909225675378995162533846867224E0_cfp
        w(2)=1.20866703423439084116213660903811244E15_cfp
        x(3)=2.84939455237483106521942546106113306E0_cfp
        w(3)=8.76311919304763792551418382215774484E16_cfp
        x(4)=3.19949264691022555380026417645412837E0_cfp
        w(4)=2.12698132156764944597399147179225494E18_cfp
        x(5)=3.54439255460183613512917816706160966E0_cfp
        w(5)=2.28797106160449049956402277595353430E19_cfp
        x(6)=3.88841218602373491975620676302139190E0_cfp
        w(6)=1.26005232166761434161163078455844292E20_cfp
        x(7)=4.23432434242045536427221402755641203E0_cfp
        w(7)=3.85882800841056668537642232191275702E20_cfp
        x(8)=4.58423735582396001752867985834967746E0_cfp
        w(8)=6.89754452885708541237424045714733596E20_cfp
        x(9)=4.94000633510456208057656239880778529E0_cfp
        w(9)=7.38840783193647900433176228019712468E20_cfp
        x(10)=5.30347956016569232785008232111301042E0_cfp
        w(10)=4.78926851144317575429153533956546371E20_cfp
        x(11)=5.67669596648544975512625069117637202E0_cfp
        w(11)=1.86976614909834343613038783126767296E20_cfp
        x(12)=6.06209956704786336523777298634037765E0_cfp
        w(12)=4.31077039843678957780554256122160097E19_cfp
        x(13)=6.46283737576085919678533866649946679E0_cfp
        w(13)=5.65229725470626318838692349164790776E18_cfp
        x(14)=6.88325272497961109196169221209002964E0_cfp
        w(14)=3.95956160728687665543716967661210811E17_cfp
        x(15)=7.32982646469711525513821300073502979E0_cfp
        w(15)=1.33826372672839090481041580986231076E16_cfp
        x(16)=7.81326939238925599211280679883552660E0_cfp
        w(16)=1.83312938165449845064176179118975011E14_cfp
        x(17)=8.35429005060840040986257174897721629E0_cfp
        w(17)=7.27287032573564424959464740551347839E11_cfp
        x(18)=9.00706373425093158439350739425283940E0_cfp
        w(18)=3.58902422493339194799776576976936914E8_cfp
        return
      end if
!
      if(n.eq.18 .and. L .eq.24) then
        x(1)=2.17045881325478049196715569340253207E0_cfp
        w(1)=5.67387093251901772569994651016608164E13_cfp
        x(2)=2.57392894485200039966396109416931740E0_cfp
        w(2)=2.56704398340558554121511902587295100E16_cfp
        x(3)=2.94010508854556627651575677731482356E0_cfp
        w(3)=1.90550378343931648733164502320779843E18_cfp
        x(4)=3.29204836764636155984329372024509618E0_cfp
        w(4)=4.71449408984582244955831333586871323E19_cfp
        x(5)=3.63829375495560334431698892956246238E0_cfp
        w(5)=5.15399313519584989306596752604744908E20_cfp
        x(6)=3.98328409756303171354159431563205840E0_cfp
        w(6)=2.87829164021070899547727564830348090E21_cfp
        x(7)=4.32987171504845626947556023473820994E0_cfp
        w(7)=8.92264501153556143640414959458088189E21_cfp
        x(8)=4.68021976446296390693872212757414155E0_cfp
        w(8)=1.61215790504209265335421462439723814E22_cfp
        x(9)=5.03622300382885640038582528710655631E0_cfp
        w(9)=1.74352179297680798833365742976330791E22_cfp
        x(10)=5.39975914413938243897766738828255416E0_cfp
        w(10)=1.13993201817405723524107131344466860E22_cfp
        x(11)=5.77288910742801670861175654383607507E0_cfp
        w(11)=4.48502380653359712262085739744983158E21_cfp
        x(12)=6.15807300623361429103690873907539362E0_cfp
        w(12)=1.04132826908519427807084116310035364E21_cfp
        x(13)=6.55846879945624341257905798911136363E0_cfp
        w(13)=1.37418723448204464967240749502381088E20_cfp
        x(14)=6.97842557344470715029873212350026646E0_cfp
        w(14)=9.68343047875989668057918419851519494E18_cfp
        x(15)=7.42442361002177749675764415125484928E0_cfp
        w(15)=3.29076190024958611639335189311256406E17_cfp
        x(16)=7.90716336074223134605522351123728169E0_cfp
        w(16)=4.53078202141979094124044117724938746E15_cfp
        x(17)=8.44732280328958786613586303542128661E0_cfp
        w(17)=1.80642370870972189572678248632183301E13_cfp
        x(18)=9.09898291829757947222058693974974971E0_cfp
        w(18)=8.95877742782814566481519421010296156E9_cfp
        return
      end if
!
      end subroutine quadrature_u_integral

end module general_quadrature_gbl
