!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use hypre_str_class,   only: hypre_str
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use sgsmodel_class,    only: sgsmodel     ! SGS model for eddy viscosity
   ! Classes used in particle injection case
   use lpt_class,         only: lpt
   use partmesh_class,    only: partmesh
   implicit none
   private
   
   !> LPT Solver
   type(lpt),         public :: lp

   !> Single low Mach flow solver and scalar solver and corresponding time tracker
   type(hypre_str),   public :: ps
   type(ddadi),       public :: vs,ss
   type(lowmach),     public :: fs
   type(vdscalar),    public :: sc
   type(timetracker), public :: time
   type(sgsmodel),    public :: sgs   !< SGS model for eddy viscosity 
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   type(partmesh) :: pmesh
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile, lptfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resSC
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   
   !> Equation of state
   real(WP) :: rho0
   real(WP) :: Z_jet
   real(WP) :: D_jet
   real(WP) :: U_jet
   
   !> Integral of pressure residual
   real(WP) :: int_RP=0.0_WP

   !> Fluid definition
   real(WP) :: visc

   !> Mixing Layer Problem definition
   real(WP) :: delta,Udiff
   integer :: nwaveX,nwaveZ
   real(WP), dimension(:), allocatable :: wnumbX,wshiftX,wampX,wnumbZ,wshiftZ,wampZ

contains
   
   
   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator
   
   
   !> Function that localizes y+ boundary
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes z- boundary
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator
   
   
   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   

   !> Function that localizes jet at -x
   function jet(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: radius
      logical :: isIn
      isIn=.false.
      ! Jet in yz plane
      radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
      if (radius.le.0.5_WP*D_jet.and.i.eq.pg%imin) isIn=.true.
   end function jet
   
   !> Function that localizes jet at -x
   function jetsc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: radius
      logical :: isIn
      isIn=.false.
      ! Jet in yz plane
      radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
      if (radius.le.0.5_WP*D_jet.and.i.eq.pg%imin- 1) isIn=.true.
   end function jetsc

   !> Obtain density from equation of state based on Burke-Schumann
   subroutine get_rho()
      implicit none
      integer :: i,j,k
      do k=sc%cfg%kmino_,sc%cfg%kmaxo_
         do j=sc%cfg%jmino_,sc%cfg%jmaxo_
            do i=sc%cfg%imino_,sc%cfg%imaxo_
               sc%rho(i,j,k) = rho0    ! No density changes for nonreactive case
            end do
         end do
      end do
   end subroutine get_rho
  
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Read in the EOS info
      call param_read('rho0',rho0)      

      ! Read in inlet information
      call param_read('Z jet',Z_jet)
      call param_read('D jet',D_jet)
      call param_read('U jet',U_jet)
      
      ! Create a low-Mach flow solver with bconds
      create_velocity_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use lowmach_class,   only: dirichlet,clipped_neumann,slip
         real(WP) :: visc
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define jet boundary conditions
         call fs%add_bcond(name='jet'   ,type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=jet   )
         ! Use slip on the sides with correction 
         !call fs%add_bcond(name='yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         !call fs%add_bcond(name='ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         !call fs%add_bcond(name='zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         !call fs%add_bcond(name='zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Outflow on the right
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=18
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_velocity_solver
      
      
      ! Create a scalar solver
      create_scalar: block
         use vdscalar_class, only: dirichlet,neumann,quick
         real(WP) :: diffusivity
         ! Create scalar solver
         sc=vdscalar(cfg=cfg,scheme=quick,name='MixFrac')
         ! Define jet boundary conditions
         call sc%add_bcond(name='jet'   ,type=dirichlet,locator=jetsc   )
         ! Outflow on the right
         call sc%add_bcond(name='outflow',type=neumann,locator=xp_locator,dir='+x')
         ! Assign constant diffusivity
         call param_read('Dynamic diffusivity',diffusivity)
         sc%diff=diffusivity
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=13)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_scalar
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         ! Flow solver
         allocate(gradU(1:3,1:3,fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))   
         allocate(resU(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resV(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resW(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         ! Scalar solver
         allocate(resSC(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=fs%cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Initialize our mixture fraction field
      initialize_scalar: block
         use vdscalar_class, only: bcond
         integer :: n,i,j,k
         type(bcond), pointer :: mybc
         ! Zero initial field
         sc%SC=0.0_WP
         ! Apply BCs
         call sc%get_bcond('jet',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            sc%SC(i,j,k)=Z_jet
         end do
         ! Compute density
         call get_rho()
      end block initialize_scalar
      
      ! Initialize our velocity field
      initialize_velocity: block
         use lowmach_class, only: bcond
         integer :: n,i,j,k
         type(bcond), pointer :: mybc
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Set density from scalar
         fs%rho=sc%rho
         ! Form momentum
         call fs%rho_multiply
         ! Apply BCs
         call fs%apply_bcond(time%t,time%dt)
         call fs%get_bcond('jet',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)   =U_jet
            fs%rhoU(i,j,k)=U_jet*rho0   ! No Burke-Schumann for nonreactive case         
         end do
         ! Get cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         resSC=0.0_WP; call fs%get_div(drhodt=resSC)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity
      
      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs

      ! Initialize our LPT
      initialize_lpt: block
         use param,  only: param_read
         use random, only: random_uniform
         real(WP) :: dp
         integer :: i,np
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Get number of particles
         call param_read('Number of particles',np)
         ! Root process initializes np particles randomly
         if (lp%cfg%amRoot) then
            call lp%resize(np)
            do i=1,np
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! Assign random position in domain
               lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
               &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
               &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
               if (lp%cfg%nx.eq.1) lp%p(i)%pos(1)=0.0_WP
               if (lp%cfg%nz.eq.1) lp%p(i)%pos(3)=0.0_WP
               ! Give zero velocity
               lp%p(i)%vel=0.0_WP
               lp%p(i)%angVel=0.0_WP
               ! Zero out collision forces
               lp%p(i)%Acol=0.0_WP
               lp%p(i)%Tcol=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
         ! Distribute particles
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
      end block initialize_lpt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         pmesh=partmesh(nvar=0,nvec=0,name='lpt')
         call lp%update_partmesh(pmesh)
      end block create_pmesh
      
      ! Add Ensight output - vdjet
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='vdjet')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('density',sc%rho)
         call ens_out%add_scalar('mixfrac',sc%SC)
         ! Add particles and VF to output - particles
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl) ! from particle injection case
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(sc%SCmax,'Zmax')
         call mfile%add_column(sc%SCmin,'Zmin')
         call mfile%add_column(sc%rhomax,'RHOmax')
         call mfile%add_column(sc%rhomin,'RHOmin')
         call mfile%add_column(int_RP,'Int(RP)')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')      
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(sc%SCint,'SC integral')
         call consfile%add_column(sc%rhoint,'RHO integral')
         call consfile%add_column(sc%rhoSCint,'rhoSC integral')
         call consfile%write()
         
         ! Create LPT monitor
         call lp%get_max()
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()

      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old scalar
         sc%rhoold=sc%rho
         sc%SCold =sc%SC
         
         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! For SGS model
         ! Turbulence modeling
         sgs_modeling: block
            use sgsmodel_class, only: vreman
            resU=fs%rho
            call fs%get_gradu(gradU)
            call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
            fs%visc=visc+sgs%visc
         end block sgs_modeling
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ============= SCALAR SOLVER =======================
            ! Build mid-time scalar
            sc%SC=0.5_WP*(sc%SC+sc%SCold)
            
            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Assemble explicit residual
            resSC=time%dt*resSC-(2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)
            
            ! Form implicit residual
            call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            
            ! Advance scalar field
            sc%SC=2.0_WP*sc%SC-sc%SCold+resSC
            
            ! Apply all other boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            dirichlet_scalar: block
               use vdscalar_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               call sc%get_bcond('jet',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc%SC(i,j,k)=Z_jet
               end do
            end block dirichlet_scalar
            ! ===================================================
            
            ! ============ UPDATE PROPERTIES ====================
            ! Backup rhoSC
            !resSC=sc%rho*sc%SC
            ! Update density
            call get_rho()
            ! Rescale scalar for conservation
            !sc%SC=resSC/sc%rho
            ! UPDATE THE VISCOSITY
            ! UPDATE THE DIFFUSIVITY
            ! ===================================================
            
            ! ============ VELOCITY SOLVER ======================
            
            ! Build n+1 density
            fs%rho=0.5_WP*(sc%rho+sc%rhoold)
            
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Update momentum
            call fs%rho_multiply()

            ! Apply boundary conditions
            call fs%apply_bcond(time%tmid,time%dtmid)
            dirichlet_velocity: block
               use lowmach_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               call fs%get_bcond('jet',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k)   =U_jet
                  fs%rhoU(i,j,k)=U_jet*rho0 ! No Burke-Schumann for nonreactive case
               end do
            end block dirichlet_velocity
            
            ! Solve Poisson equation
            call sc%get_drhodt(dt=time%dt,drhodt=resSC)
            call fs%correct_mfr(drhodt=resSC)
            call fs%get_div(drhodt=resSC)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            call cfg%integrate(A=fs%psolv%rhs,integral=int_RP)
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide
            ! ===================================================
            
            ! Advance particles by dt
            call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc)

            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call sc%get_drhodt(dt=time%dt,drhodt=resSC)
         call fs%get_div(drhodt=resSC)
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         call lptfile%write()
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi, gradU)
      
   end subroutine simulation_final
   
   
end module simulation