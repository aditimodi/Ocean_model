program main   
 
  use transport_mod
  use filter_mod
  use grid_mod, only : define_grid, dxU, dxT, nlat, nlon, dyT, dyV
  use grid_mod, only : cos_latT, cos_latV
  use press_grad_mod
  implicit none
  
  integer :: nt=300
  real, allocatable, dimension(:,:) :: u, v, h, hn, un, vn
  real :: dt
  integer :: t, i, j

  call define_grid() 

  allocate(u(nlon,nlat))
  allocate(v(nlon,nlat))
  allocate(h(nlon,nlat))
  allocate(hn(nlon,nlat))
  allocate(un(nlon,nlat))
  allocate(vn(nlon,nlat))
  
  ! initial conditions

  u(1:nlon,1:nlat) = 0.
  v(1:nlon,1:nlat) = 0.
  h(:,:) = 10. 

  do i = 1,51
    do j = 1,11
    h(i+nlon/2,j+nlat/2) = sin((i-1.)/50.*3.14) * &
    sin((j-1.)/10.*3.14) + 10. 
  end do
  end do

  dt = 1800.0

  open(9,file="test_transport.out")
!  write(9,*) cos_latT(200,:)
!  write(9,*) dxT(200,:)
!  write(9,*) dyT(200,:)
  write(9,*) h(:,:)
  
  do t = 1,nt

    h(1,:)=h(nlon-1,:)
    h(nlon,:)=h(2,:)
    h(:,nlat) = h(:,nlat-1) 
    h(:,1) = h(:,2)
    u(1,:)=u(nlon-1,:)
    u(nlon,:)=u(2,:)
    v(1,:)=v(nlon-1,:)
    v(nlon,:)=v(2,:)

 !   call pressure_grad_2d(u,v,h,dt,dyV,dyV,un,vn)
    call pressure_grad_2d(u,v,h,dt,dxU,dyV,un,vn)
    u(2:nlon-1,2:nlat-1) = un(2:nlon-1,2:nlat-1)
    v(2:nlon-1,2:nlat-1) = vn(2:nlon-1,2:nlat-1)
    u(1,:)=u(nlon-1,:)
    u(nlon,:)=u(2,:)
    v(1,:)=v(nlon-1,:)
    v(nlon,:)=v(2,:)

  !  call transport_2d(u,v,h,dt,dyT,dyT,hn)
    call transport_2d(u,v,h,dt,dxT,dyT,hn,wdy=cos_latT,wny=cos_latV)
    h(2:nlon-1,2:nlat-1)=hn(2:nlon-1,2:nlat-1)
    h(1,:)=h(nlon-1,:)
    h(nlon,:)=h(2,:)
    h(:,nlat) = h(:,nlat-1) 
    h(:,1) = h(:,2)
    
    call shapiro_filter_2d(h,w=cos_latT)
    
    if (mod(t,20)==0) then 
      write(9,*) h(:,:)
      print *,t
    end if

  end do

end program main
