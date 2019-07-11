
module transport_mod
implicit none

 contains

  subroutine transport_1d(u,h,dt,dx,hn,wd,wn)

    real, intent(in) :: u(:),h(:),dt,dx(:)
    real, intent(out) :: hn(:)
    real, intent(in), optional :: wd(:), wn(:)
   
    integer :: n,k
    real :: wd1(size(h)), wn1(size(u))
    
    wd1=1.0
    wn1=1.0

    if (present(wd)) wd1=wd
    if (present(wn)) wn1=wn

    n = size(h)

    do k = 2,n-1 
      hn(k) = h(k) - (dt/(dx(k)*wd1(k))) * &
        ( wn1(k) * (0.5*(u(k) + abs(u(k))) * h(k) + &
                   0.5*(u(k) - abs(u(k))) * h(k+1)) - &
        wn1(k-1) * (0.5*(u(k-1) + abs(u(k-1))) * h(k-1) + &
                    0.5*(u(k-1) - abs(u(k-1))) * h(k)))

    end do
  end subroutine transport_1d

  subroutine transport_2d(u,v,h,dt,dx,dy,hn,wdx,wnx,wdy,wny)

    real, intent(in), dimension(:,:) :: u,v,h,dx,dy
    real, intent(in), dimension(:,:), optional :: wdx, wdy, wnx, wny
    real, intent(out) :: hn(:,:)
    real, intent(in) :: dt
   
    integer :: nx,ny,i,j
    real, dimension(size(h,1),size(h,2)) :: hnx,hny,wdx1,wdy1,wnx1,wny1

    ny = size(h,2)
    nx = size(h,1)
    wdx1=1.0
    wdy1=1.0
    wnx1=1.0
    wny1=1.0

    if (present(wdx)) wdx1=wdx
    if (present(wdy)) wdy1=wdy
    if (present(wnx)) wnx1=wnx
    if (present(wny)) wny1=wny

    do j = 2,ny-1
      
      call transport_1d(u(:,j),h(:,j),dt,dx(:,j),hnx(:,j),wdx1(:,j),wnx1(:,j))
    
    end do

    do i = 2,nx-1
    
      call transport_1d(v(i,:),h(i,:),dt,dy(i,:),hny(i,:),wdy1(i,:),wny1(i,:))
    
    end do

    hn(:,:) = hnx(:,:) + hny(:,:) - h(:,:)

 end subroutine transport_2d

end module transport_mod

