SUBROUTINE entrainment ()

use global_vars
 
implicit none

! the 4 main arrays

real interp_d ( jt+1, kt+1 )
real lr_edges ( jt+1, kt )
real tb_edges ( jt, kt+1 )

! make the interp_d grid
! must choose what the real number d is based off of
! such that d >= 0 is "convective" and d < 0 not
! for example, d cannot be set equal to qc because qc cannot have values less than zero
! after Dawe and Austin (2011) we define d as q_diff which is the excess of total water
! over the saturation specific humidity

! first find qsat(T,p)

qdiff = qw(J,K) - ( 0.622 * ES( theta_l(J, K) * pi_0(K) ) / ( 100 * pbar(K) - ES( theta_l(J, K) * pi_0(K) ) ) )


do J = 1, jt+1
  do K = 1, kt+1
    interp_d(J, K) = ( ( qw(J,K) - ( 0.622 * ES( theta_l(J, K) * pi_0(K) ) / ( 100 * pbar(K) - ES( theta_l(J, K) * pi_0(K) ) ) ) ) &
      + ( qw(J+1,K) - ( 0.622 * ES( theta_l(J+1, K) * pi_0(K) ) / ( 100 * pbar(K) - ES( theta_l(J+1, K) * pi_0(K) ) ) ) ) &
      + ( qw(J,K+1) - ( 0.622 * ES( theta_l(J, K+1) * pi_0(K+1) ) / ( 100 * pbar(K+1) - ES( theta_l(J, K+1 * pi_0(K+1) ) ) ) ) &
      + ( qw(J+1,K+1) - ( 0.622 * ES( theta_l(J+1, K+1) * pi_0(K+1) ) / ( 100 * pbar(K+1) - ES( theta_l(J+1, K+1) * pi_0(K+1) ) ) ) ) ) &
      / 4
  end do
end do

