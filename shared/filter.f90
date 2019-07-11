
module filter_mod

  implicit none

  real :: epsln=0.01

  contains

    subroutine shapiro_filter_1d(h,w)

      real, intent(inout) :: h(:)
      real, intent(in), optional :: w(:)
      
      real :: hn(size(h))
      real :: w1(size(h))
      integer :: nk, k
      nk=size(h)

      w1=1.0

      if (present(w)) w1=w

      do k=2,nk-1

        hn(k) = (1.-epsln)*h(k) +  epsln * &
          (w1(k-1)*h(k-1) + w1(k+1)* h(k+1)) / &
          (w1(k-1) + w1(k+1))

      end do

      h(2:nk-1) = hn(2:nk-1)

    end subroutine shapiro_filter_1d

    subroutine shapiro_filter_2d(h,w)

      real, intent(inout) :: h(:,:)
      real, intent(in), optional :: w(:,:)
      
      real :: hn(size(h,1),size(h,2))
      real :: w1(size(h,1),size(h,2))
      integer :: nx, ny, i, j
      nx = size(h,1)
      ny = size(h,2)

      w1(:,:)=1.0

      if (present(w)) w1=w

      do i=2,nx-1
        do j=2,ny-1

        hn(i,j) = (1.-epsln)*h(i,j) + epsln * &
          (w1(i-1,j)*h(i-1,j) + w1(i+1,j)*h(i+1,j) + &
          w1(i,j-1)*h(i,j-1) + w1(i,j+1)*h(i,j+1)) / &
          (w1(i-1,j) + w1(i+1,j) + &
          w1(i,j-1) + w1(i,j+1)) 
        
        end do
      end do

      h(2:nx-1,2:ny-1) = hn(2:nx-1,2:ny-1)

    end subroutine shapiro_filter_2d


end module filter_mod


