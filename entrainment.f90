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
do J = 1, jt+1
  do K = 1, kt+1
    interp_d(J, K) = ( ( qw(J,K) - ( 0.622 * ES( theta_l(J, K) * pi_0(K) ) / ( 100 * pbar(K) - ES( theta_l(J, K) * pi_0(K) ) ) ) ) &
      + ( qw(J+1,K) - ( 0.622 * ES( theta_l(J+1, K) * pi_0(K) ) / ( 100 * pbar(K) - ES( theta_l(J+1, K) * pi_0(K) ) ) ) ) &
      + ( qw(J,K+1) - ( 0.622 * ES( theta_l(J, K+1) * pi_0(K+1) ) / ( 100 * pbar(K+1) - ES( theta_l(J, K+1 * pi_0(K+1) ) ) ) ) &
      + ( qw(J+1,K+1) - ( 0.622 * ES( theta_l(J+1, K+1) * pi_0(K+1) ) / ( 100 * pbar(K+1) - ES( theta_l(J+1, K+1) * pi_0(K+1) ) ) ) ) ) &
      / 4
  end do
end do

! define the possible categories
onec = [ 1, 2, 4, 8 ]
threec = [ 7, 11, 13, 14 ]
two_adj = [ 3, 6, 9, 12 ]
two_opp = [ 5, 10 ]

! now go over each box and check the vertices
do J = 1, jt
  do K = 1, kt
    a = interp_d(J, K+1) >= 0
    b = interp_d(J+1, K+1) >= 0
    c = interp_d(J+1, K) >= 0
    d = interp_d(J, K) >= 0
    
    mcode = 0
    if a then
      mcode = mcode + 8
    end if
    if b then
      mcode = mcode + 4
    end if
    if c then
      mcode = mcode + 2
    end if
    if d then
      mcode = mcode + 1
    end if
    
    if (mcode == 0) then
      cldarea(J, K, 2) = 0
    end if
    if (mcode == 15) then
      cldarea(J, K, 2) = 1
    end if
    
  end do
end do

