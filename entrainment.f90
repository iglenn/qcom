SUBROUTINE entrainment ()

!use global_vars ! in particular, cldarea(jt, kt) and udotw(jt, kt)

use global_vars, only : cldarea, udotw, J, K, qw, theta_l, pi_0, pbar, dy, dz, w, v
implicit none

real interp_d ( jt+1, kt+1 ) ! the convective condition
real lr_edges ( jt+1, kt ) ! left - right grid edges
real tb_edges ( jt, kt+1 ) ! top - bottom grid edges

integer, dimension(4,1) :: onec, threec, two_adj
integer, dimension(2,1) :: two_opp
integer mcode
logical a, b, c, d

! make the interp_d grid
! must choose what the real number d is such that d >= 0 is "convective" and d < 0 not
! for example, d cannot be set equal to qc because qc cannot have values less than zero
! after Dawe and Austin (2011) I define d as q_diff which is the excess of total water
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

! define the possible categories, besides all or nothing
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
    
    ! all
    if (mcode == 15) then
      cldarea(J, K, 2) = ( dy * dz )
    end if
    ! or nothing
    if (mcode == 0) then
      cldarea(J, K, 2) = 0
    end if

    ! one corner
    if (ANY(mcode == onec)) then
      ! go through each one
      if (mcode == 1) then ! it's d
        lr_edges(J,K) = - interp_d(J, K) * dz / ( interp_d(J, K+1) - interp_d(J, K)  )
        tb_edges(J,K) = - interp_d(J, K) * dy / ( interp_d(J+1, K) - interp_d(J, K)  )
        cldarea(J,K) = 0.5 * lr_edges(J,K) * tb_edges(J,K) 
        udotw(J,K) = lr_edges(J,K) * v(J+1,K+1) + tb_edges(J,K) * w(J+1, K+1)
        CYCLE
      end if
      if (mcode == 2) then ! it's c
        lr_edges(J+1,K) = - interp_d(J+1, K) * dz / ( interp_d(J+1, K+1) - interp_d(J+1, K)  )
        tb_edges(J,K) = - interp_d(J+1, K) * dy / ( interp_d(J, K) - interp_d(J+1, K)  )
        cldarea(J,K) = 0.5 * lr_edges(J+1,K) * tb_edges(J,K)
        udotw(J,K) = lr_edges(J+1,K) * v(J+2,K+1) + tb_edges(J,K) * w(J+1, K+1)
        CYCLE
      end if
      if (mcode == 4) then ! it's b
        lr_edges(J+1,K) = - interp_d(J+1, K+1) * dz / ( interp_d(J+1, K) - interp_d(J+1, K+1)  )
        tb_edges(J,K+1) = - interp_d(J+1, K+1) * dy / ( interp_d(J, K+1) - interp_d(J+1, K+1)  )
        cldarea(J,K) = 0.5 * lr_edges(J+1,K) * tb_edges(J,K+1)
        udotw(J,K) = lr_edges(J+1,K) * v(J+2,K+1) + tb_edges(J,K+1) * w(J+1, K+2)
        CYCLE
      end if
      if (mcode == 8) then ! it's a
        lr_edges(J,K) = - interp_d(J, K+1) * dz / ( interp_d(J, K) - interp_d(J, K+1)  )
        tb_edges(J,K+1) = - interp_d(J, K+1) * dy / ( interp_d(J+1, K+1) - interp_d(J, K+1)  )
        cldarea(J,K) = 0.5 * lr_edges(J,K) * tb_edges(J,K+1)
        udotw(J,K) = lr_edges(J,K) * v(J+1,K+1) + tb_edges(J,K+1) * w(J+1, K+2)
        CYCLE
      end if
    end if ! any onec
    
    ! three corners
    if (ANY(mcode == threec)) then ! same as before, just have to find inverse
      ! go through each one
      if (mcode == 14) then ! it's d
        lr_edges(J,K) = interp_d(J, K) * dz / ( interp_d(J, K+1) - interp_d(J, K)  )
        tb_edges(J,K) = interp_d(J, K) * dy / ( interp_d(J+1, K) - interp_d(J, K)  )
        cldarea(J,K) = ( dy * dz ) - 0.5 * lr_edges(J,K) * tb_edges(J,K) 
        udotw(J,K) = ( dz - lr_edges(J,K) ) * v(J+1,K+1) + ( dy - tb_edges(J,K) ) * w(J+1, K+1) + dz * v(J+2,K+1) + dy * w(J+1,K+2)
        CYCLE
      end if
      if (mcode == 13) then ! it's c
        lr_edges(J+1,K) = interp_d(J+1, K) * dz / ( interp_d(J+1, K+1) - interp_d(J+1, K)  )
        tb_edges(J,K) = interp_d(J+1, K) * dy / ( interp_d(J, K) - interp_d(J+1, K)  )
        cldarea(J,K) = ( dy * dz ) -  0.5 * lr_edges(J+1,K) * tb_edges(J,K)
        udotw(J,K) = ( dz - lr_edges(J+1,K) ) * v(J+2,K+1) + ( dy - tb_edges(J,K) ) * w(J+1, K+1) + dz * v(J+1,K+1) + dy * w(J+1,K+2)
        CYCLE
      end if
      if (mcode == 11) then ! it's b
        lr_edges(J+1,K) = interp_d(J+1, K+1) * dz / ( interp_d(J+1, K) - interp_d(J+1, K+1)  )
        tb_edges(J,K+1) = interp_d(J+1, K+1) * dy / ( interp_d(J, K+1) - interp_d(J+1, K+1)  )
        cldarea(J,K) = ( dy * dz ) - 0.5 * lr_edges(J+1,K) * tb_edges(J,K+1) 
        udotw(J,K) = ( dz - lr_edges(J+1,K) ) * v(J+2,K+1) + ( dy - tb_edges(J,K+1) ) * w(J+1, K+2) + dz * v(J+1,K+1) + dy * w(J+1,K+1)
        CYCLE
      end if
      if (mcode == 7) then ! it's a
        lr_edges(J,K) = interp_d(J, K+1) * dz / ( interp_d(J, K) - interp_d(J, K+1)  )
        tb_edges(J,K+1) = interp_d(J, K+1) * dy / ( interp_d(J+1, K+1) - interp_d(J, K+1)  )
        cldarea(J,K) =  ( dy * dz ) - 0.5 * lr_edges(J,K) * tb_edges(J,K+1)
        udotw(J,K) = ( dz - lr_edges(J,K) ) * v(J+1,K+1) + ( dy - tb_edges(J,K+1) ) * w(J+1, K+2) + dz * v(J+2,K+1) + dy * w(J+1,K+1)
        CYCLE
      end if
    end if ! any threec
    
    ! two corners, adjacent
    if (ANY(mcode == two_adj) then
      if (mcode == 3) then ! its c and d
        ! from d to a
        lr_edges(J,K) = - interp_d(J, K) * dz / ( interp_d(J,K) - interp_d(J, K+1) )
        ! from c to b
        lr_edges(J+1,K) =  - interp_d(J+1, K) * dz / ( interp_d(J+1,K+1) - interp_d(J+1, K) )
        cldarea(J,K) = 0.5 * dy * lr_edges(J,K) * lr_edges(J+1,K)
        udotw(J,K) = lr_edges(J,K) * v(J+1,K+1) + lr_edges(J+1,K) * v(J+2,K+1) + dy * w(J+1,K+1)
        CYCLE
      end if
      if (mcode == 12) then ! its a and b, just invert
        ! from d to a
        lr_edges(J,K) = dz - interp_d(J, K) * dz / ( interp_d(J,K) - interp_d(J, K+1) )
        ! from c to b
        lr_edges(J+1,K) = dz - interp_d(J+1, K) * dz / ( interp_d(J+1,K+1) - interp_d(J+1, K) )
        cldarea(J,K) = 0.5 * dy * lr_edges(J,K) * lr_edges(J+1,K)
        udotw(J,K) = lr_edges(J,K) * v(J+1,K+1) + lr_edges(J+1,K) * v(J+2,K+1) + dy * w(J+1,K+2)
        CYCLE
      end if
      if (mcode == 9) then ! its a and d
        ! from d to c
        tb_edges(J,K) = - interp_d(J, K) * dy / ( interp_d(J+1,K) - interp_d(J, K) )
        ! from a to b
        tb_edges(J,K+1) =  - interp_d(J, K+1) * dy / ( interp_d(J+1,K+1) - interp_d(J, K+1) )
        cldarea(J,K) = 0.5 * dz * tb_edges(J,K) * tb_edges(J,K+1)
        udotw(J,K) = tb_edges(J,K) * w(J+1,K+1) + tb_edges(J,K+1) * w(J+1,K+2) + dz * v(J+1,K+1)
        CYCLE
      end if
      if (mcode == 6) then ! its b and c, just invert
        ! from d to c
        tb_edges(J,K) = dy - interp_d(J, K) * dy / ( interp_d(J+1,K) - interp_d(J, K) )
        ! from a to b
        tb_edges(J,K+1) = dy - interp_d(J, K+1) * dy / ( interp_d(J+1,K+1) - interp_d(J, K+1) )
        cldarea(J,K) = 0.5 * dz * tb_edges(J,K) * tb_edges(J,K+1)
        udotw(J,K) = tb_edges(J,K) * w(J+1,K+1) + tb_edges(J,K+1) * w(J+1,K+2) + dz * v(J+2,K+1)
        CYCLE
      end if
    end if ! two adjacent
    
    ! two corners, opposite
    if (ANY(mcode == two_opp)) then
      if (mcode == 5) then ! b and d
        ! from b to c
        lr_edges(J+1,K) = - interp_d(J+1, K+1) * dz / ( interp_d(J+1, K) - interp_d(J+1, K+1)  )
        ! b to a
        tb_edges(J,K+1) = - interp_d(J+1, K+1) * dy / ( interp_d(J, K+1) - interp_d(J+1, K+1)  )
        ! d to a
        lr_edges(J,K) = - interp_d(J, K) * dz / ( interp_d(J, K+1) - interp_d(J, K)  )
        ! d to c
        tb_edges(J,K) = - interp_d(J, K) * dy / ( interp_d(J+1, K) - interp_d(J, K)  )
        cldarea(J,K) = 0.5 * ( lr_edges(J+1,K) * tb_edges(J,K+1) + lr_edges(J,K) * tb_edges(J,K) )
        udotw(J,K) = lr_edges(J,K) * v(J+1,K+1) + lr_edges(J+1,K) * v(J+2,K+1) + tb_edges(J,K) * w(J+1,K+1) + tb_edges(J,K+1) * w(J+1,K+2)
        CYCLE
      end if
      if (mcode == 10) then ! a and c
        ! from b to c
        lr_edges(J+1,K) = dz - interp_d(J+1, K+1) * dz / ( interp_d(J+1, K) - interp_d(J+1, K+1)  )
        ! b to a
        tb_edges(J,K+1) = dy - interp_d(J+1, K+1) * dy / ( interp_d(J, K+1) - interp_d(J+1, K+1)  )
        ! d to a
        lr_edges(J,K) = dz - interp_d(J, K) * dz / ( interp_d(J, K+1) - interp_d(J, K)  )
        ! d to c
        tb_edges(J,K) = dy - interp_d(J, K) * dy / ( interp_d(J+1, K) - interp_d(J, K)  )
        cldarea(J,K) = 0.5 * ( lr_edges(J+1,K) * tb_edges(J,K+1) + lr_edges(J,K) * tb_edges(J,K) )
        udotw(J,K) = lr_edges(J,K) * v(J+1,K+1) + lr_edges(J+1,K) * v(J+2,K+1) + tb_edges(J,K) * w(J+1,K+1) + tb_edges(J,K+1) * w(J+1,K+2)
        CYCLE
      end if
    end if ! any two opposite
    
  end do ! K
end do ! J


