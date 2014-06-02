MODULE global_vars
! here we declare the arrays and parameters that will be accessible 
! by the main program and all subroutines

! Arakawa c-grid
integer, parameter :: jt = 600
integer, parameter :: kt = 120
integer, parameter :: kw = kt

! time
real, parameter :: dt = 0.1 ! timestep, s
real, parameter :: tmax = 7200. ! run time, s
integer, parameter :: ittmax = int( tmax / dt ) ! approx number of timesteps
integer ITT ! the timestep index

! space
real, parameter :: H = 6000. ! m, the height of the simulated domain
real, parameter :: L = 30000. ! the horizontal domain length

real, parameter :: dz = H / real( kt )
real, parameter :: dy = L / real( jt ) 

! vectors along z dimension
real, dimension(1:kt+2):: theta_0, theta_v0
real, dimension(1:kt+2):: qv_0
real, dimension(1:kw+2):: pi_0
real, dimension(1:kw+2):: pbar

! large 2D arrays
real, dimension(1:jt+2, 1:kt+2) :: theta_l, theta
real, dimension(1:jt+2, 1:kt+2) :: qw, qc, qv
real, dimension(1:jt+2, 1:kw+2) :: pi_1
real, dimension(1:jt+2, 1:kw+2) :: w
real, dimension(1:jt+2, 1:kw+2) :: v

real, dimension(1:jt+2, 1:kt+2, 1:2) :: ftheta_l
real, dimension(1:jt+2, 1:kt+2, 1:2) :: fqw
real, dimension(1:jt+2, 1:kw+2, 1:2) :: fpi_1
real, dimension(1:jt+2, 1:kw+2, 1:2) :: fw
real, dimension(1:jt+2, 1:kw+2, 1:2) :: fv

! physical constants
real, parameter :: g = 9.81 ! gravity, m/s2
real, parameter :: c_p = 1004. ! specific heat, J/kg/K
real, parameter :: R = 287. ! gas constant, J/kg/K
real, parameter :: L_f = 2.5104e+06 ! Latent heat of condensation, J/kg

real, parameter :: c_s = 50. ! speed of sound, m/s

real, parameter :: K_th = 50. ! eddy diffusivity of theta, m2/s
real, parameter :: K_v = 50. ! eddy diffusivity of v, m2/s
real, parameter :: K_w = 50. ! eddy diffusivity of w, m2/s

! parameters to be assigned later, useful to be global
real theta_top
real theta_bot
real A
real B

! frequently used indicies
integer J
integer K
integer N1
integer N2

END MODULE global_vars

PROGRAM qcom

USE global_vars

IMPLICIT NONE
!     pgf90 -o qcom qcom.f90 
!     ./qcom

integer, parameter rcount = 1 ! record counter for writing direct access

write(*,*) "horizontal resolution="
write(*,*) jt
write(*,*) "vertical resolution="
write(*,*) kt
write(*,*) "dz="
write(*,*) dz
write(*,*) "dy="
write(*,*) dy

! open files to write output, the 4 in recl signifies 4 bits per byte
open(unit = 50, file='gate2D_bub_602x122_6km_qc.dat', &
form='unformatted', access='direct', recl = (jt+2)*(kt+2)*4, action='write')
open(unit = 51, file='gate2D_bub_602x122_6km_qw.dat', &
form='unformatted', access='direct', recl = (jt+2)*(kt+2)*4, action='write')
open(unit = 52, file='gate2D_bub_602x122_6km_thetal.dat', &
form='unformatted', access='direct', recl = (jt+2)*(kt+2)*4, action='write')
open(unit = 53, file='gate2D_bub_602x122_6km_w.dat', &
form='unformatted', access='direct', recl = (jt+2)*(kt+2)*4, action='write')
open(unit = 54, file='gate2D_bub_602x122_6km_pi1.dat', &
form='unformatted', access='direct', recl = (jt+2)*(kt+2)*4, action='write')

open(unit = 55, file='gate2D_bub_602x122_6km_pbar.dat', &
form='unformatted', access='direct', recl = (kt+2)*4, action='write')


! begin
ITT = 1 ! itt is time step index

! initialize domain
CALL init

! initialization complete, call bound
CALL bound

! do any initial perturbation
CALL perturb

! write initial state of variables
write(50, rec = ITT) qc
write(51, rec = ITT) qw
write(52, rec = ITT) theta_l
write(53, rec = ITT) w
write(54, rec = ITT) pi_1
write(55, rec = ITT) pbar

! USE FORWARD SCHEME TO do first step
! ADAMS - BASHFORTH coefficients A and B
A = 1.
B = 0.
N1 = MOD ( ITT    , 2 ) + 1
N2 = MOD ( ITT - 1, 2 ) + 1

CALL rcalc ! calculate forcing terms
CALL AB ! update variables using a time scheme
CALL bound

! ADAMS - BASHFORTH TWO - LEVEL SCHEME

A =   3. / 2. 
B = - 1. / 2. 

DO ITT = 2, ittmax

if (MOD(ITT,100)==0) then
    write(*,*) (100*ITT/ittmax) ! displays simulation progress
end if

! make any perturbation
!if (MOD(ITT,6000)==0) then ! every 10 minutes
!CALL perturb
!end if

N1 = MOD ( ITT    , 2 ) + 1
N2 = MOD ( ITT - 1, 2 ) + 1

CALL rcalc ! calculate forcing terms
CALL AB ! update variables using a time scheme
CALL bound

! write variables every 60 seconds
if (MOD(ITT,600)==0) then
rcount = rcount + 1
! write initial state of variables
write(50, rec = rcount) qc
write(51, rec = rcount) qw
write(52, rec = rcount) theta_l
write(53, rec = rcount) w
write(54, rec = rcount) pi_1
write(55, rec = rcount) pbar

end if ! write

end do ! time loop

write(*,*) "Number of records written:"
write(*,*) rcount

contains

SUBROUTINE init()

USE global_vars

IMPLICIT NONE

real stab_fac ! factor to reduce stability
real RH ! relative humidity

! Define the initial atmospheric state
! these vectors must match the z grid, kt+2 elements long
open(unit = 23, file = 'theta0_gigaLES.txt')
open(unit = 24, file = 'qv0_gigaLES.txt')

do K=1,kt+2
     read(23,*) theta_0(K)
     read(24,*) qv_0(K)
end do
! correct units
qv_0 = qv_0 / 1000. ! kg/kg

! factor to reduce the stability (optional)
stab_fac = 0.
theta_0 = theta_0 - stab_fac * ( theta_0 - MINVAL( theta_0 ) )

! assign theta top and bottom values
theta_bot = ( theta_0(1) + theta_0(2) ) / 2.
theta_top = theta_0(kt+2)

! define pressure based on theta_v0
theta_v0 = theta_0 * ( 1 + 0.61 * qv_0 )
pi_0(1) = ( 1015. / 1000. ) ** ( R / c_p ) ! set base pressure (dz/2 underground) to be 1015mb
do K = 2, kt+2
     pi_0(K) = pi_0(K-1) - dz * g / ( c_p * ( theta_v0(K) + theta_v0(K-1) ) / 2. )
end do
pi_0(kt+2) = pi_0(kt+1)
pbar = 1000. * ( pi_0 ** ( c_p / R ) ) !mb

! set the RH to a constant value? redefine qv_0
!RH = 0.9
!do K = 1, kt + 2
!    qv_0(K) = RH * ( 0.622 * ES( theta_0(K) * pi_0(K) ) / ( 100 * pbar(K) - ES( theta_0(K) * pi_0(K) ) ) )
!end do
! now recalculate things that depend on qv_0
!theta_v0 = theta_0 * ( 1 + 0.61 * qv_0 )
!do K = 2, kt+2
!     pi_0(K) = pi_0(K-1) - dble(dz) * g / ( c_p * ( theta_v0(K) + theta_v0(K-1) ) / 2. )
!end do
!pi_0(kt+2) = pi_0(kt+1)
!pbar = 1000. * ( pi_0 ** ( c_p / R ) ) !mb


!     initialize the main part of the domain
DO J = 2, jt+1
     theta_l(J, 1:kt+2) = theta_0
     qv(J, 1:kt+2) = qv_0
END DO

theta = theta_l
qw = qv

v(1:jt+2, 1:kt+2) = 0.
w(1:jt+2, 1:kt+2) = 0.
pi_1(1:jt+2, 1:kt+2) = 0.
qc(1:jt+2, 1:kt+2) = 0.

fv(1:jt+2, 1:kt+2, 1:2) = 0.
fw(1:jt+2, 1:kt+2, 1:2) = 0.
ftheta_l(1:jt+2, 1:kt+2, 1:2) = 0.
fqw(1:jt+2, 1:kt+2, 1:2) = 0.
fpi_1(1:jt+2, 1:kt+2, 1:2) = 0.

END SUBROUTINE init

SUBROUTINE perturb()

use global_vars

IMPLICIT NONE
integer iyp
integer izp
real tbub
real esbub
real supersat
real rnd

! release a bubble

do iyp = 290, 310 ! that's 1km at 50m res
do izp = 10, 30 ! that's 1km at 50m z resolution 

theta_l(iyp, izp) = theta_l(iyp, izp) &
+ 1. * COS( ( 3.14159 / 2. ) * ( & 
0.5 * ( ( ( real(iyp) * dy - 15000. ) / 500. ) ** 2. ) &
+ 0.5 * ( ( ( real(izp) * dz - 1000. ) / 500. ) ** 2 ) ) ** 0.5 ) ** 2

tbub = theta_l(iyp,izp) * pi_0(izp)
esbub = ES( tbub )
! define the initial supersaturation of the bubble
supersat = 1.0
qw(iyp,izp) = supersat * ( 0.622 * esbub / ( 100 * pbar(izp) - esbub ) )

!qw(iyp, izp) = qw(iyp, izp) &
!+ 0.001 * COS( ( 3.14159 / 2. ) * ( & 
!0.5 * ( ( ( real(iyp) * dy - 15000. ) / 500. ) ** 2. ) &
!+ 0.5 * ( ( ( real(izp) * dz - 1000. ) / 500. ) ** 2 ) ) ** 0.5 ) ** 2

!w(iyp, izp) = w(iyp, izp) &
!+ 1. * COS( ( 3.14159 / 2. ) * ( & 
!0.5 * ( ( ( real(iyp) * dy - 15000. ) / 500. ) ** 2. ) &
!+ 0.5 * ( ( ( real(izp) * dz - 1000. ) / 500. ) ** 2 ) ) ** 0.5 ) ** 2

end do
end do


! randomly perturb a section
!CALL RANDOM_SEED ! initializes the PRNG
!izp = 15 ! z(15) = 1km altitude
!do iyp = 290, 310 ! 1 km wide in the center, dy = 50m
!    CALL RANDOM_NUMBER(rnd)
!    theta_l(iyp,izp) = theta_l(iyp,izp) + 1. * rnd
!    CALL RANDOM_NUMBER(rnd)
!    w(iyp,izp) = w(iyp,izp) + 1. * rnd
!    CALL RANDOM_NUMBER(rnd)
!    qw(iyp,izp) = qw(iyp,izp) + 0.001 * rnd
!end do

END SUBROUTINE perturb

SUBROUTINE bound()

use global_vars

IMPLICIT NONE

! conduction, constant boundary temp
! this is effectively a large scale forcing
!theta_l( 2:jt+1, 1) = 2. * theta_bot - theta_l( 2:jt+1, 2)
!theta_l( 2:jt+1, kt+2) = 2. * theta_top - theta_l( 2:jt+1, kt+1)

! roof layer equals level below
theta_l( 1:jt+2, kt+2) = theta_l( 1:jt+2, kt+1)
theta( 1:jt+2, kt+2) = theta( 1:jt+2, kt+1)
qv( 1:jt+2, kt+2) = qv( 1:jt+2, kt+1)
qc( 1:jt+2, kt+2) = qc( 1:jt+2, kt+1)
qw( 1:jt+2, kt+2) = qw( 1:jt+2, kt+1)
pi_1( 1:jt+2, kt+2) = pi_1( 1:jt+2, kt+1)

! free slip boundary
v(2:jt+1,1)=v(2:jt+1,2)
v(2:jt+1,kt+2)=v(2:jt+1,kt+1) 
! the term above used to =v(2:jt+1,kt+1)
! ground and roof condition for w
w(2:jt+1,1:2)=0. ! layer 1 is underground, layer 2 is ground level
w(2:jt+1,kt+2)=0.
        
! should find new qv0, theta_0 ??? otherwise isn't buoyancy forcing out of touch?
! for ik=1:kt+2
!     qv_0(ik)=mean(qv(2:jt+1,ik));
! end
        
! periodic boundaries
theta_l(1,1:kt+2)=theta_l(jt+1,1:kt+2)
theta_l(jt+2,1:kt+2)=theta_l(2,1:kt+2)
qw(1,1:kt+2)=qw(jt+1,1:kt+2)
qw(jt+2,1:kt+2)=qw(2,1:kt+2)
w(1,1:kt+2)=w(jt+1,1:kt+2)
w(jt+2,1:kt+2)=w(2,1:kt+2)
v(1,1:kt+2)=v(jt+1,1:kt+2)
v(jt+2,1:kt+2)=v(2,1:kt+2)
pi_1(1,1:kt+2)=pi_1(jt+1,1:kt+2)
pi_1(jt+2,1:kt+2)=pi_1(2,1:kt+2)
               
qv(1,1:kt+2)=qv(jt+1,1:kt+2)
qv(jt+2,1:kt+2)=qv(2,1:kt+2)
qc(1,1:kt+2)=qc(jt+1,1:kt+2)
qc(jt+2,1:kt+2)=qc(2,1:kt+2)
theta(1,1:kt+2)=theta(jt+1,1:kt+2)
theta(jt+2,1:kt+2)=theta(2,1:kt+2)

END SUBROUTINE bound

SUBROUTINE rcalc

use global_vars

IMPLICIT NONE

real d1,d2,d3,d4 ! dummy variables

!CALCULATES FORCING TERMS FOR V(J,K), ETC.; STORES THEM IN FV(J,K,N2)

! call adjust to see if any saturation might occur given the current situation

do K = 2, kt+1
do J = 2, jt+1
      ! use qc to get theta from theta_l
      theta(J,K) = theta_l(J,K) + ( L_f / ( c_p * pi_0(K) ) ) * qc(J,K)
      d1 = theta(J,K)
      d2 = qw(J,K) - qc(J,K)
      d3 = qc(J,K)
      d4 = 100. * pbar(K)
      CALL ADJUST(d1,d2,d3,d4)
      theta(J,K) = d1
      qv(J,K) = d2
      qc(J,K) = d3 
end do
end do

DO K = 2, kt+1
DO J = 2, jt+1
fv(J,K,N2) = -v(J,K)*((v(J+1,K)-v(J-1,K))/(2.*dy)) &
                 -((w(J-1,K)+w(J,K))/2.)*((v(J,K)-v(J,K-1))/dz)/2. &
                 -((w(J-1,K+1)+w(J,K+1))/2.)*((v(J,K+1)-v(J,K))/dz)/2. &
                 - c_p*theta_v0(K)*(pi_1(J,K)-pi_1(J-1,K))/dy &
                 + K_v*((v(J,K+1)-2.*v(J,K)+v(J,K-1))/(dz*dz) &
                 + (v(J+1,K)-2.*v(J,K)+v(J-1,K))/(dy*dy)) 
 
fw(J,K,N2) = -w(J,K)*((w(J,K+1)-w(J,K-1))/(2.*dz)) &
                 -((v(J,K-1)+v(J,K))/2.)*((w(J,K)-w(J-1,K))/dy)/2. &
                 -((v(J+1,K-1)+v(J+1,K))/2.)*((w(J+1,K)-w(J,K))/dy)/2. &
                 -c_p*((theta_v0(K-1)+theta_v0(K))/2.) &
                 * (pi_1(J,K)-pi_1(J,K-1))/(dz) &
                 + g*(((theta(J,K-1)+theta(J,K))/(theta_0(K-1) &
                 + theta_0(K))) - 1.+0.61*(((qv(J,K)+qv(J,K-1))/2.) &
                 - ((qv_0(K)+qv_0(K-1))/2.))-((qc(J,K)+qc(J,K-1))/2.)) &
                 + K_w*((w(J,K+1)-2.*w(J,K)+w(J,K-1))/(dz*dz) &
                 + (w(J+1,K)-2.*w(J,K)+w(J-1,K))/(dy*dy)) 
             
ftheta_l(J,K,N2) = -v(J,K)*((theta_l(J,K)-theta_l(J-1,K))/dy)/2. &
                 -v(J+1,K)*((theta_l(J+1,K)-theta_l(J,K))/dy)/2. &
                 -w(J,K)*((theta_l(J,K)-theta_l(J,K-1))/dz)/2. &
                 -w(J,K+1)*((theta_l(J,K+1)-theta_l(J,K))/dz)/2. &
       +K_th*((theta_l(J,K+1)-2.*theta_l(J,K)+theta_l(J,K-1))/(dz*dz) &
           +(theta_l(J+1,K)-2.*theta_l(J,K)+theta_l(J-1,K))/(dy*dy))
             
fqw(J,K,N2) = -v(J,K)*((qw(J,K)-qw(J-1,K))/dy)/2. &
                 -v(J+1,K)*((qw(J+1,K)-qw(J,K))/dy)/2. &
                 -w(J,K)*((qw(J,K)-qw(J,K-1))/dz)/2. &
                 -w(J,K+1)*((qw(J,K+1)-qw(J,K))/dz)/2. &
                 +K_th*((qw(J,K+1)-2.*qw(J,K)+qw(J,K-1))/(dz*dz)& 
                 +(qw(J+1,K)-2.*qw(J,K)+qw(J-1,K))/(dy*dy))
 
fpi_1(J,K,N2)= -((c_s**2)/(c_p*(theta_0(K)**2.))) &
                   * ( (v(J+1,K)*theta_v0(K) &
                   -  v(J,K)*theta_v0(K))/(dy) &
                   + (w(J,K+1)*((theta_v0(K+1)+theta_v0(K))/2.) &
                   -  w(J,K)*((theta_v0(K)+theta_v0(K-1))/2.))/(dz) )

END DO
END DO


END SUBROUTINE rcalc

SUBROUTINE AB

use global_vars

IMPLICIT NONE
integer ij, ik

!     THE FOLLOWING LOOP UPDATES V USING EITHER THE FORWARD OR THE ADAMS-BASHFORTH 
!     SCHEME DEPENDING ON THE VALUES OF A, B.  
!     SUBSCRIPT N2 OF FV ALWAYS REFERS TO THE MOST RECENTLY CALCULATED VALUES FOR FV.

DO ik = 2, kt+1
DO ij = 2, jt+1
theta_l(ij,ik)=theta_l(ij,ik)+dt*(A*ftheta_l(ij,ik,N2)+B*ftheta_l(ij,ik,N1))
w(ij,ik)=w(ij,ik)+dt*(A*fw(ij,ik,N2)+B*fw(ij,ik,N1))
v(ij,ik)=v(ij,ik)+dt*(A*fv(ij,ik,N2)+B*fv(ij,ik,N1))
pi_1(ij,ik)=pi_1(ij,ik)+dt*(A*fpi_1(ij,ik,N2)+B*fpi_1(ij,ik,N1))
qw(ij,ik)=qw(ij,ik)+dt*(A*fqw(ij,ik,N2)+B*fqw(ij,ik,N1))

END DO
END DO

END SUBROUTINE AB

SUBROUTINE ADJUST ( TH, QV, QC, PBAR_in )
IMPLICIT NONE
!     SUBROUTINE ADJUST ( TH, QV, QC, PBAR, qvs )
!
!     Performs an isobaric moist adiabatic adjustment.
!     The final state is either subsaturated with no liquid water
!     present, or exactly saturated with liquid water present.
!
!     This version iterates to obtain the adjustment, so accurate
!     guesses for TH and QV are not needed. For use in situations
!     where TH and QV will be successively updated, as in a parcel
!     model, set ITTMAX=1.
!
!     Units: SI (MKS)
!
!     Input --
!     TH: potential temperature, theta^* (K)
!     QV: mixing ratio of water vapor, q_v^* (kg/kg)
!     QC: mixing ratio of liquid water, q_c^* (kg/kg)
!     PBAR: pressure, p (Pa)
!
!     Output --
!     TH: adjusted potential temperature, theta^{n+1} (K)
!     QV: adjusted mixing ratio of water vapor, q_v^{n+1} (kg/kg)
!     QC: adjusted mixing ratio of liquid water, q_c^{n+1} (kg/kg)
!     QVS: saturation mixing ratio based on TH and PBAR
!
!     29 Nov 91 --
!     Exact Qsat now used in ALPHA and QSAT.
!     Corrected QVS1 for QC=0.
!     5 Feb 92 --
!     PIBAR calculated within subroutine now.
!     6 Feb 92 --
!     GAM corrected!
!

real, intent(IN) :: PBAR_in
real, intent(INOUT) :: TH, QV, QC
real THSTAR,TH1,QVSTAR,QV1,QCSTAR,esd
integer IT
real PIBAR,GAM,Tstar,es1,ALPHA,THFAC,QVSAT,QW1,QC1,QVS1,DTH,QVS
real, parameter :: HLF = 2.5104e+06 ! Latent heat of condensation, J/kg
real, parameter :: CP = 1004. ! specific hear, J/kg/K
real, parameter :: rgas = 287. ! gas constant, J/kg/K
real, parameter :: pzero = 1.0e+05 ! reference pressure, J/m3 = Pa
real, parameter :: DTCRIT = 0.001 ! tolerance, K
integer, parameter :: ITMAX = 2 ! max iterations

!     PBAR IS A ( HYDROSTATIC ) REFERENCE PRESSURE FIELD.
!
      IT = 1
      PIBAR = ( PBAR_in / pzero ) ** ( rgas / CP )      
GAM = HLF / ( CP * PIBAR )
!
      THSTAR = TH
      QVSTAR = QV
      QCSTAR = QC
!

   30 Tstar = THSTAR * PIBAR
      es1 = ES( Tstar )
!
      ALPHA = DESDT(Tstar)*0.622*PIBAR*PBAR_in/((PBAR_in-es1)**2)
      THFAC = GAM / ( 1. + GAM * ALPHA )
!
      QVSAT = 0.622 / ( PBAR_in - es1 ) * es1
      TH1 = THSTAR + THFAC * ( QVSTAR - QVSAT )
      QV1 = QVSAT + ALPHA * ( TH1 - THSTAR )
!     statement below gives same result as one above
!     QV1 = QVSTAR - ( TH1 - THSTAR ) / GAM
      QC1=QVSTAR-QV1+QCSTAR      
      !QW1 = QV + QC
      !QC1 = QW1 - QV1
!
      QVS1 = QV1
!
      IF ( QC1 .LT. 0. ) then
        QC1 = 0.
        QV1 = QVSTAR+QCSTAR
        TH1 = THSTAR + GAM * ( QVSTAR - QV1 )
        esd = ES( ( TH1 * PIBAR ) )
        QVS1 = 0.622 / ( PBAR_in -esd ) * esd
        !QVS1 = QVSAT + alpha * ( TH1 - THSTAR )
      ENDIF
!
      DTH = ( TH1 - THSTAR ) * PIBAR
!
      IF ( ABS( DTH ) .LT. DTCRIT .OR. IT .EQ. ITMAX ) then
        TH = TH1
        QV = QV1
        QC = QC1
        QVS = QVS1
      IF (ABS(DTH) .GT. DTCRIT) then
       write(*,*) "adjust fail"
      end if
        return
      endif
!
      THSTAR = TH1
      QVSTAR = QV1
      QCSTAR = QC1
!
      IT = IT + 1
!
      GO TO 30
!
END SUBROUTINE ADJUST

FUNCTION ES ( T )
IMPLICIT NONE
!
!     LOWE'S FORMULA FOR SATURATION VAPOR PRESSURE ( PA ).
!     T IS IN DEGREES KELVIN.
real, intent(IN) :: T
real TC, X
integer jj
real, dimension(1:7) :: C

C = ( /6.107800, &
       4.436519E-01, &
       1.428946E-02, &
       2.650648E-04, &
       3.031240E-06, &
       2.034081E-08, &
       6.136821E-11/ )
       
TC = T - 273.16
if ( TC .LT. - 50. ) then
      TC = - 50.
end if

X = C(7)
do jj = 1, 6
      X = X * TC + C(7-jj)
end do

ES = X * 100.
RETURN
      
END FUNCTION ES
      
FUNCTION DESDT ( T )
IMPLICIT NONE
!
!     LOWE'S FORMULA FOR THE DERIVATIVE OF
!     SATURATION VAPOR PRESSURE WITH RESPECT TO TEMPERATURE.
!     ES IS IN PASCALS. T IS IN DEGREES KELVIN.
!

real,intent(IN) :: T
real TC, X
integer jj
real, dimension(1:7) :: D

D = ( /4.438100E-01, &
       2.857003E-02, &
       7.938054E-04, &
       1.215215E-05, &
       1.036561E-07, &
       3.532422E-10, &
     - 7.090245E-13/ )

TC = T - 273.16
if ( TC .LT. - 50. ) then
TC = - 50.
end if

X = D(7)
do jj = 1, 6
      X = X * TC + D(7-jj)
end do

DESDT = X * 100.
RETURN

END FUNCTION DESDT
      
END PROGRAM qcom
