
module filter_mod

  implicit none

  real :: epsln=0.01

  contains

    subroutine shapiro_filter_1d(h)

      real, intent(inout) :: h(:)
      
      real :: hn(size(h))
      integer :: nk, k
      nk=size(h)

      do k=2,nk-1

        hn(k) = (1.-epsln)*h(k) + 0.5 * epsln * &
          (h(k-1) + h(k+1))

      end do

      h(2:nk-1) = hn(2:nk-1)

    end subroutine shapiro_filter_1d

end module filter_mod


