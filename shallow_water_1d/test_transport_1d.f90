program test_transport_1d
 
  use transport_mod
 
  implicit none
  
  integer, parameter :: nk = 1000, nt = 11000
  real :: u(nk), h(nk), dx(nk), dt, hn(nk)
  integer :: t

  u(1:nk) = 10.
  h(:) = 0. 
  h(2:11) = 1.
  dx(:) = 100.
  dt = 1.
  
  open(9,file="test_transport.out")
  write(9,*) h
  
  do t = 1,nt
    h(1)=h(nk-1)
    h(nk)=h(2)
    call transport_1d(u,h,dt,dx,hn)
    h(2:nk-1)=hn(2:nk-1)
  end do

  write(9,*) h

end program test_transport_1d





