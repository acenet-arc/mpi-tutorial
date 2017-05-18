       program diffuse
       implicit none
   
!
! Simulation parameters
!
       integer, parameter :: totpoints=1000
       real, parameter    :: xleft=-12., xright=+12.
       real, parameter    :: kappa=1.
       integer, parameter :: nsteps=100000
       integer, parameter :: plotsteps=50

!
! The calculated temperature, 
! and the known correct solution from theory
!
       real, allocatable :: x(:)
       real, allocatable, target :: temperature(:,:)
       real, allocatable :: theory(:)
       real, dimension(:), pointer :: old, new, tmp
       real :: fixedlefttemp, fixedrighttemp

       integer :: step
       integer :: i
       integer :: red, grey, white
       real :: time
       real :: dt, dx
       real :: error

!
! Parameters of the original temperature distribution
!
       real, parameter :: ao=1., sigmao = 1.
       real :: a, sigma

!
! Cell size and timestep size
!
       dx = (xright-xleft)/(totpoints-1)
       dt = dx**2 * kappa/10.

! 
! Allocate data, including ghost cells, for old and new timesteps.
! Theory doesn't need ghost cells, but we include them for simplicity.
!
       allocate(temperature(totpoints+2,2))
       allocate(theory(totpoints+2))
       allocate(x(totpoints+2))
!
! Set up initial conditions
!
       old => temperature(:,1)
       new => temperature(:,2)
       time = 0.
       x = xleft + [((i-1+0.5)*dx,i=1,totpoints+2)]
       old    = ao*exp(-(x)**2 / (2.*sigmao**2))
       theory= ao*exp(-(x)**2 / (2.*sigmao**2))

       fixedlefttemp = ao*exp(-(xleft-dx)**2 / (2.*sigmao**2))
       fixedrighttemp= ao*exp(-(xright+dx)**2 / (2.*sigmao**2))
!
! Set up plotting with PGPLOT
!
       call pgbeg(0, "/xwindow", 1, 1)
       call pgask(0)
       call pgenv(xleft, xright, 0., 1.5*ao, 0, 0)
       call pglab('x', 'Temperature', 'Diffusion Test')
       red = 2
       call pgscr(red,1.,0.,0.)
       grey = 3
       call pgscr(grey,.2,.2,.2)
       white = 4
       call pgscr(white,1.,1.,1.)
 
       call pgsls(1)
       call pgsci(white)
       call pgline(totpoints, x(2), theory(2))
       call pgsls(2)
       call pgsci(red)
       call pgline(totpoints, x(2), old(2))
!
! Evolve
!
       do step=1, nsteps
!
! Boundary conditions: keep endpoint temperatures fixed.
!
           old(1) = fixedlefttemp
           old(totpoints+2) = fixedrighttemp
!
! Update solution
!
           forall (i=2:totpoints+1)
               new(i) = old(i) + dt*kappa/(dx**2) * &
                         (old(i+1) - 2*old(i) + old(i-1))
           end forall
           time = time + dt

           if (mod(step, plotsteps) == 0) then
                call pgbbuf()
                call pgeras
                call pgsls(2)   !dashed
                call pgslw(12)  !thick
                call pgsci(red)
                call pgline(totpoints, x(2), new(2))
            endif
! 
! Update theoretical (correct) solution
!
           sigma = sqrt(2.*kappa*time + sigmao**2)
           a = ao*sigmao/sigma
           theory = a*exp(-(x)**2 / (2.*sigma**2))

           if (mod(step, plotsteps) == 0) then
               call pgsls(1)   
               call pgslw(6)   
               call pgsci(white)
               call pgline(totpoints, x(2), theory(2))
               call pgebuf()
           endif
           error = sqrt(sum((theory(1:totpoints+1) - new(1:totpoints+1))**2))
           print *, 'Step = ', step, 'Time = ', time, ' Err = ', error

           tmp => old
           old => new
           new => tmp
       enddo

       call pgend
       deallocate(temperature)
       deallocate(theory)
       deallocate(x)
       end
