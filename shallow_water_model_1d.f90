program main   
 
  use transport_mod
  use filter_mod

  implicit none
  
  integer, parameter :: nk = 1000, nt = 3000
  real :: u(nk), h(nk), dx(nk), dt, hn(nk), un(nk)
  integer :: t, i

  u(1:nk) = 0.
  h(:) = 10. 

  do i = 1,51

    h(i+nk/2) = sin((i-1.)/50.*3.14)*1. + 10. 
    
  end do

  dx(:) = 100.
  dt = 1.
  
  open(9,file="test_transport.out")
  
  do t = 1,nt

    h(1)=h(nk-1)
    h(nk)=h(2)
    u(1)=u(nk-1)
    u(nk)=u(2)
    call pressure_grad(u,h,dt,dx,un)
    u(2:nk-1) = un(2:nk-1)
    u(1)=u(nk-1)
    u(nk)=u(2)
    call transport_1d(u,h,dt,dx,hn)
    h(2:nk-1)=hn(2:nk-1)
    h(1)=h(nk-1)
    h(nk)=h(2)
    call shapiro_filter_1d(h)
    write(9,*) h
  end do

  contains

    subroutine pressure_grad(u,h,dt,dx,un)

      real, intent(in) :: u(:),h(:),dt,dx(:)
      real, intent(out) :: un(:)
      
      integer :: nk,k
      real, parameter :: g=9.81
      nk = size(u)
      
      do k=2,nk-1
        un(k) = u(k) - dt * g * (h(k+1) - h(k)) / dx(k)
      end do

    end subroutine pressure_grad  

end program main





