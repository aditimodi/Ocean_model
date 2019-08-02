
module coriolis_mod

  use constants_mod, only : OMEGA, DEG_TO_RAD

  implicit none

  real, allocatable, dimension(:,:) :: alphaU, alphaV, bnU, bnV, bdU, bdV, fU, fV
 
  contains

    subroutine init_coriolis(dt,latU,latV)

      real, intent(in) :: dt
      real, intent(in), dimension(:,:) :: latU, latV
      
      allocate(alphaU(1:size(latU,1),1:size(latU,2)))
      allocate(alphaV(1:size(latU,1),1:size(latU,2)))
      allocate(bnU(1:size(latU,1),1:size(latU,2)))
      allocate(bnV(1:size(latU,1),1:size(latU,2)))
      allocate(bdU(1:size(latU,1),1:size(latU,2)))
      allocate(bdV(1:size(latU,1),1:size(latU,2)))
      allocate(fU(1:size(latU,1),1:size(latU,2)))
      allocate(fV(1:size(latU,1),1:size(latU,2)))
  
      fU = 2 * omega * sin(latU*DEG_TO_RAD) 
      fV = 2 * omega * sin(latV*DEG_TO_RAD)
  
      alphaU = dt * fU
      alphaV = dt * fV
  
      bnU = 1.0 - (0.25 * alphaU**2)
      bnV = 1.0 - (0.25 * alphaV**2)
     
      bdU = 1.0 / (1.0 + 0.25 * alphaU**2)  
      bdV = 1.0 / (1.0 + 0.25 * alphaV**2)
  
    end subroutine init_coriolis
 
    subroutine calc_coriolis(u,v,un,vn) 

      real, intent(in), dimension(:,:) :: u, v 
      real, intent(out), dimension(:,:) :: un, vn

      real :: Vu(size(u,1),size(u,2))
      real :: Uv(size(u,1),size(u,2))

      integer :: i,j

      Vu(:,:) = 0.
      Uv(:,:) = 0.

      ! calculations are done considering first and last grid as boundary 

      do i=2,size(u,1)-1
        do j=2,size(u,2)-1
        
           Vu(i,j) = (v(i,j) + v(i+1,j) + v(i,j-1) + v(i+1,j-1)) * 0.25
           Uv(i,j) = (u(i-1,j) + u(i,j) + u(i-1,j+1) + u(i,j+1)) * 0.25
          
        end do
      end do

      un(:,:) =0.
      vn(:,:) =0.

      do i=2,size(u,1)-1
        do j=2,size(u,2)-1
       
          un(i,j) = (bnU(i,j) * u(i,j) + alphaU(i,j) * Vu(i,j)) * bdU(i,j) 
          vn(i,j) = (bnV(i,j) * v(i,j) - alphaV(i,j) * Uv(i,j)) * bdV(i,j)

        end do
      end do

    end subroutine calc_coriolis

end module coriolis_mod


