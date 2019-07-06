
module transport_mod
implicit none

 contains

  subroutine transport_1d(u,h,dt,dx,hn)

    real, intent(in) :: u(:),h(:),dt,dx(:)
    real, intent(out) :: hn(:)
   
    integer :: n,k
    
    n = size(h)

    do k = 2,n-1 
      hn(k) = h(k) - (dt/dx(k)) * &
      (0.5 * (u(k) + abs(u(k))) * h(k) + &
      0.5 * (u(k) - abs(u(k))) * h(k+1) - &
      0.5 * (u(k-1) + abs(u(k-1))) * h(k-1) - &
      0.5 * (u(k-1) - abs(u(k-1))) * h(k))

    end do
  end subroutine transport_1d
end module transport_mod

