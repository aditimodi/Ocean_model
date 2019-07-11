module press_grad_mod

  implicit none

 contains

    subroutine pressure_grad_1d(u,h,dt,dx,un)

      real, intent(in) :: u(:),h(:),dt,dx(:)
      real, intent(out) :: un(:)
      
      integer :: nx,k
      real, parameter :: g=9.81
      nx = size(u)
      
      do k=2,nx-1
        un(k) = u(k) - dt * g * (h(k+1) - h(k)) / dx(k)
      end do

    end subroutine pressure_grad_1d  


    subroutine pressure_grad_2d(u,v,h,dt,dx,dy,un,vn)

      real, intent(in), dimension(:,:) :: u,v,h,dx,dy
      real, intent(in) :: dt
      real, intent(out), dimension(:,:) :: un,vn
      
      integer :: nx,k,ny,i,j
      real, parameter :: g=9.81
      nx = size(u,1)
      ny = size(u,2)
      
      do j = 2,ny-1
        call pressure_grad_1d(u(:,j),h(:,j),dt,dx(:,j),un(:,j))
      end do

      do i = 2,nx-1
        call pressure_grad_1d(v(i,:),h(i,:),dt,dy(i,:),vn(i,:))
      end do

    end subroutine pressure_grad_2d  

end module press_grad_mod






