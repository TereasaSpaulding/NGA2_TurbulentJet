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
         integer :: i,j,k,nx,ny,nz,n_jet, n_domain, jj
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z, dy_array, dy_norm
         real(WP) :: Djet, dy_max, dy_min, a, b, Norm

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         call param_read('D jet', Djet)
         !call param_read('n_jet', n_jet)
         !call param_read('dy min', dy_min)
         allocate(dy_array(ny))
         allocate(dy_norm(ny))

         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do

         ! Number of grid cells in jet domain
         n_jet = 32 ! -> Make input parameter

         ! Size/Number of large grid cells
         n_domain = ny - n_jet
         dy_max = Ly*1_WP/real(n_domain,WP)
 
         ! Scaling Factors 
         ! Make input param? to avoid complex numbers
         a = -440.9990
         b = -27.2211
         

         !> Start at y = 0, Create upper half of y grid, ny should be even
         do j=ny/2 +2,ny+1
            jj = j - (0.5*ny+1) ! Offset for index so tanh(x) is close to 0 near j = ny/2
            dy_array(j-1) = jj*jj + a*jj + b             
         end do

         !> Mirror other half of y grid
         do j=1,ny/2
            dy_array(j) = dy_array(ny + 1 - j)
         end do

         ! Normalize so that the sum of the cell sizes = Ly
         Norm = sum(dy_array)/Ly
         do j = 1,ny
            dy_norm(j) = dy_array(j)/Norm
         end do

         ! Create grid in y
         y(1) = -0.5_WP*Ly ! Initial point 
         do j = 2,ny+1
            y(j) = y(j-1) + dy_norm(j-1);
         end do

         
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
