subroutine entrainment ( jt_in, kt_in )

implicit none

! the 4 main arrays
integer, parameter :: jt = jt_in
integer, parameter :: kt = kt_in

real interp_d ( jt+1, kt+1 )
real lr_edges ( jt+1, kt )
real tb_edges ( jt, kt+1 )
real cldarea ( jt, kt, 2 )
