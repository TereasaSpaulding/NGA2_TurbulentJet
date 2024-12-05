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
         integer :: n_jet, n_large, jj, exp_count, min_idx, ny_mask
         real(WP), dimension(:), allocatable :: dy, y_mask, dy_mask
         real(WP) :: Djet, y0, scaling, refined_area
         real(WP) :: n_jet_count, exp ! debugging

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! For refined geometry 
         allocate(y_mask(ny+3))
         allocate(dy_mask(ny+1))
         call param_read('D jet', Djet) ! Read in particle diameter

         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         

         !> Grid in Y
         ny_mask = ny + 3 ! Add 3 since you lose 2 points due to difference + 1 will give you 0
         refined_area = Djet/2.0_WP
         n_jet = 15
         exp = 1.0_WP 

         ! Exponent must be odd
         exp_loop: do exp_count = 1,10 

            exp = exp + 2.0_WP

            do j=1,ny_mask
               y0 = real(j, WP)/real(ny_mask, WP)    ! Initial input values for power function
               y_mask(j) = (y0 - 0.5_WP)**exp ! Initial power function
            end do

            scaling = Ly/sum(y_mask(2:ny_mask-1) - y_mask(1:(ny_mask-2)) ) ! Factor to scale by to ensure domain length = Ly
            y_mask = scaling*y_mask;

            dy_mask = abs(y_mask(2:ny_mask-1) - y_mask(1:ny_mask -2))   ! Find grid size

            ! Find index of leftmost element within the radius
            min_idx = (ny_mask - 3) / 2 + 1 - n_jet / 2; 

            ! Check if elements within refined region are enough
            ! Allows +/- 3 cells in jet domain
            !if (abs(sum(dy_mask(min_idx:(min_idx + n_jet - 1))) - refined_area*2.0_WP) <= 3.0_WP) then
            !   exit exp_loop
            !end if

            ! Forces no less than n_jet cells in jet
            if (sum(dy_mask(min_idx:(min_idx + n_jet - 1))) <= refined_area * 2.0_WP) then
               exit exp_loop
            end if

         end do exp_loop

         dy = [dy_mask(1:((ny_mask - 3)/2 - 1)), dy_mask(((ny_mask - 3)/2 + 1):(ny_mask -2))];

         y(1) = -0.5_WP*Ly ! Initial point 
         do j = 2,ny+1
            y(j) = y(j-1) + dy(j-1);
         end do

         ! Double check number of elements in refined area - for debugging
         !n_jet_count = 0
         !do j = 2,ny+1
         !   if (abs(y(j) - 0) <= refined_area) then
         !      n_jet_count = n_jet_count + 1.0_WP
         !   end if
         !end do
         !print *, n_jet_count
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
