!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z

         ! Added for refined geometry
         integer :: n_jet
         real(WP), dimension(:), allocatable :: dy
         real(WP) :: Djet, scale_factor, jj, b, d, c, dy_min
         real(WP) :: n_jet_count ! debugging

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))

         ! For refined geometry 
         allocate(dy(ny))
         call param_read('D jet', Djet) ! Read in particle diameter
         call param_read('n jet', n_jet) ! Number of cells in jet domain


         ! Create simple rectilinear grid in X and Z
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         !> Grid in Y
         dy_min = Djet/real(n_jet,WP)

         ! Function Parameters
         scale_factor = 1.0_WP
         c = dy_min + scale_factor
         d = 0.85_WP

         do j = ny/2 + 2,ny+1
            jj = real(j - ny/2 -1, WP) ! Mask for index  
            dy(j-1) = scale_factor*tanh(d*jj - real(n_jet,WP)/(2.0_WP*d)) + c 
         end do

         ! Mirror other half of y grid
         do j= 1,ny/2
            dy(j) = dy(ny - j + 1)
         end do

         ! Scale grid cells so the sum = Ly
         scale_factor = (Ly - real(ny,WP)*dy_min) /sum(dy)
         dy = scale_factor * (dy - dy_min) + dy_min

         ! Build y array
         y(1) = -0.5_WP*Ly ! Initial point 
         do j = 2,ny+1
            y(j) = y(j-1) + dy(j-1)
         end do
         
         ! Double check number of elements in refined area - for debugging
         !n_jet_count = 0
         !do j = 2,ny+1
         !   if (abs(y(j) - 0) <= Djet/2.0_WP) then
         !      n_jet_count = n_jet_count + 1.0_WP
         !   end if
         !end do
         !print *, n_jet_count

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='vdjet')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Create masks for this config
      create_walls: block
         cfg%VF=1.0_WP
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
