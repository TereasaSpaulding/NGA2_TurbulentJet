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
         integer :: n_jet, n_large, jj, center_idx, ny_mask
         real(WP), dimension(:), allocatable :: dy, y_mask, dy_mask
         real(WP) :: Djet, y0, scaling, refined_area, avg
         real(WP) :: n_jet_count, exp ! debugging

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! For refined geometry 
         allocate(y_mask(ny+1))
         allocate(dy_mask(ny))
         call param_read('D jet', Djet) ! Read in particle diameter

         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         

         !> Grid in Y
         do j = ny/2 + 2,ny+1
            jj = (j - (ny/2+1))/ (ny/2+1) ! Mask for index  
            y_mask(j) = (cosh(jj) - 1)
         end do

         ! Mirror other half of y grid
         do j= 1,ny/2
            y_mask(j) = - y_mask(ny + 1 - j)
         end do

         scaling = Ly/sum(y_mask(2:ny+1) - y_mask(1:(ny)))
         y = scaling * y_mask
         dy = y(2:ny+1) - y(1:ny) ! find dy between elements
         center_idx = ny/2

         ! Average 5 center cells to avoid dy = 0
         avg = 1/5* sum(dy(center_idx -2: center_idx + 2)) 
         dy(center_idx-2:center_idx+2) = avg
         ! Average of center cell and +/- 3 cells
         avg = 1/2*(dy(center_idx)+ dy(center_idx - 3))
         ! Set new avg to centercel +/- 1 cell
         dy(center_idx-1) = avg
         dy(center_idx + 1) = avg
         ! Average of cell + 1idx from center and cell +3 idx from center
         avg =  1/2*(dy(center_idx - 1)+ dy(center_idx - 3))
         ! Set new avg to centercell +/- 2 cell
         dy(center_idx-2) = avg
         dy(center_idx + 2) = avg

         y(1) = -0.5_WP*Ly ! Initial point 
         do j = 2,ny+1
            y(j) = y(j-1) + dy(j-1)
         end do

         ! Double check number of elements in refined area - for debugging
         n_jet_count = 0
         do j = 2,ny+1
            if (abs(y(j) - 0) <= refined_area) then
               n_jet_count = n_jet_count + 1.0_WP
            end if
         end do
         print *, n_jet_count
         !print *, exp

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
