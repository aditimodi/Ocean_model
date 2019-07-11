
module grid_mod

  implicit none

  real, parameter :: RADIUS=637*10**4, pi=3.1415926536, deg2rad=pi/180.0

  real, allocatable, dimension(:,:) :: latT, latV, lonT, lonU, dxT, dxU, dyT, dyV, cos_latT, cos_latV
  real :: latS=-15, dlat=1, lonW=0, dlon=1
  integer ::  nlat=30, nlon=360
  logical :: cyclicX=.true.
  contains

    subroutine define_grid()

      integer :: i, j 
      
      namelist/grid_nml/latS,dlat,nlat,lonW,dlon,nlon,cyclicX
      open(9,file="input.nml",status="old")
      read(9, nml=grid_nml)
      write(6, nml=grid_nml)
      close(9)

      allocate(latT(nlon,nlat))
      allocate(latV(nlon,nlat))
      allocate(lonT(nlon,nlat))
      allocate(lonU(nlon,nlat))
      allocate(dxT(nlon,nlat))
      allocate(dxU(nlon,nlat))
      allocate(dyT(nlon,nlat))
      allocate(dyV(nlon,nlat))
      allocate(cos_latT(nlon,nlat))
      allocate(cos_latV(nlon,nlat))

      latV(:, 1)=latS
      latT(:, 1)=latS + dlat*0.5

      do j=2,nlat
        
        latV(:,j)=latV(:,j-1)+dlat
        latT(:,j)=latT(:,j-1)+dlat
        
      end do
      
      cos_latT(:,:) = cos(latT(:,:)*deg2rad)
      cos_latV(:,:) = cos(latV(:,:)*deg2rad)

      lonU(1,:)=lonW
      lonT(1,:)=lonW + dlon*0.5

      do i=2,nlon
        
        lonU(i,:)=lonU(i-1,:)+dlon
        lonT(i,:)=lonT(i-1,:)+dlon
        
      end do

      do j=1,nlat
        
        dxT(:,j) = RADIUS * cos_latT(:,j) * dlon * deg2rad
        
      end do

      dxU(:,:) = dxT(:,:)
      dyT(:,:) = RADIUS * dlat * deg2rad
      dyV(:,:) = dyT(:,:)

    end subroutine define_grid

end module grid_mod


