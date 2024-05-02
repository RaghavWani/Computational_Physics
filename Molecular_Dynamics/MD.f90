!----------------------Modules------------------------------------
module all_parameters
    implicit none
    !Problem fixed parameters:
    integer, parameter :: lx=20,ly=20,lz=20 !box dimension
    integer, parameter :: n_part=1000 !number of particles
    integer, parameter :: KbT=1, mass=1 !constants
    integer, parameter :: n_iter_eqlbm=50000 !number of iterations considered for sys to reach equilibrium
    integer, parameter :: n_iter=20000 !number of iterations
    real*8,  parameter :: r_c=2.50d0, sigma=1.0d0, r_s=4.50d0*sigma !cutoff radius, particle radius, neighbour list radius
    real*8,  parameter :: dr=0.1*sigma, dv=0.05d0 !g(r) bin size, ...
    real*8             :: post(3*n_part), velocity(3*n_part), force(3*n_part), old_force(3*n_part),acc(3*n_part), new_acc(3*n_part) !main ingredients
    integer            :: no_neigh(n_part), neigh_list(n_part,n_part)
    real*8             :: force_lj,KE,potential_energy,potential_lj
    real*8             :: gr(nint(2*lx/(2*dr)))

    !other essential params:
    real*8, parameter :: sig_6=sigma**6, sig_12=sigma**12
    real*8, parameter :: eps=4.0d0 !=4*epsilon; potential prefactor or potential depth
    real*8, parameter :: Fc=eps*((12*sig_12/(r_c**13)) - (6*sig_6/(r_c**7))) !constants in Verlet algo
    real*8, parameter :: Vc=Fc*r_c + eps*((sigma/r_c)**(12) - (sigma/r_c)**(6)) !constants in Verlet algo

    !some variables:
    real*8, parameter :: dt=0.0025 !time interval 
    real*8, parameter :: dt2by2 = (0.5*(dt**2))/mass, dtby2 = (0.5*(dt))/mass
    real*8,parameter  :: avg_rho=n_part/dfloat(lx*ly*lz) !average density of particles 
    real*8            :: dist(nint(lx/(2*dv))),vel_arr(n_part) !speed distribution
end module all_parameters

module all_subroutines
    
 contains
    include 'post_init.f90'   !to initalize position
    include 'vel_init.f90'    !to initalize velocity
    include 'calc_force.f90'  !to calculate force
    include 'update_post.f90' !to update position
    include 'update_vel.f90'  !to update velocity
    include 'thermostat.f90'  !to vary energy at each iteration
    include 'neigh_list_finder.f90' !to find neighbour list
    include 'calc_force_neighb.f90' !to calculate force using neigh_list
    include 'calc_gr.f90'     !to calculate pair correlation function
    include 'speed_dist.f90'  !to calculate speed 
    
end module all_subroutines

!----------------------Main Program------------------------------------
program molecular_dynamics
    use all_parameters
    use all_subroutines
    implicit none 
    real :: caller,avg_vx, avg_vy, avg_vz     
    integer :: time,i

    !initialization:
    time=0
    CALL post_init
    CALL vel_init
    CALL calc_force
    write(*,*) potential_energy, KE, (KE + potential_energy)
    write(*,*) 'Initization complete'

    !--------------------Equilibrating the system and storing physically quantities:
    open(15, file='energy_dt(0.0025).dat',action='write') !file saving all energy (pe,ke,tot) after each iteration
    ! write(15,*) avg_vx, avg_vy, avg_vz
    ! write(15,*) time, potential_energy, kinetic_energy, total_energy

    do time=1,n_iter_eqlbm
        CALL update_post !updating position 
        old_force = force !store force at time t
        CALL calc_force !calculate force at time t+dt
        CALL update_vel !update velocity using forces at t & t+dt

        !store physical quantities:
        if (mod(time,10)==0) then
            write(15,*) avg_vx, avg_vy, avg_vz !for momentum conservation
            write(15,*) time, potential_energy, KE, (KE + potential_energy) !for energy conservation
        end if
        if (mod(time,1000)==0) CALL thermostat  !modifying energy to classically expected value until equilibrium is reached
        if (mod(time,1000)==0) write(*,*) time, potential_energy, KE, (KE + potential_energy)
    end do
    write(*,*) "Equilibrium reached..."

    !--------------------Faster: Equilibrating the system using neighbour list 
    !Neighbour list: considering contribution only from nearby particles, that lie within radius r_s:

    open(16, file='neighbour_list.dat',action='write') !file saving neighbour list
    do time=1,n_iter_eqlbm
        CALL update_post !updating position 
        old_force = force !store force at time t
        CALL calc_force_neighb !calculate force at time t+dt, but considering only neighbouring contributions
        CALL update_vel !update velocity using forces at t & t+dt

        if (mod(time,40)==0) call neigh_list_finder !neighbour list
        if (mod(time,1000)==0) CALL thermostat !thermostat
        if (mod(time,1000)==0) write(*,*) time, potential_energy, KE, (KE + potential_energy)
    end do

    !--------------------Post eqmb, no thermostat
    gr=0
    dist=0
    do time=1,n_iter
        CALL update_post
        old_force = force
        CALL calc_force_neighb
        CALL update_vel 

        if (mod(time,40)==0) call neigh_list_finder !neighbour list
        if (mod(time,100)==0) call calc_gr !pair correlation function
        if (mod(time,100)==0) call speed_dist !speed distribution
        !store physical quantities:
        ! if (mod(time,10)==0) then
        !     write(15,*) time, potential_energy/n_part, KE/n_part, (KE + potential_energy)/n_part
        ! end if
        if (mod(time,1000)==0) write(*,*) time, potential_energy, KE, (KE + potential_energy)
    end do

    !--------------------pair correlation function: g(r) data 
    caller=n_iter/100
    open(unit=18,file='gr.dat')
    do i=1,size(gr)
        gr(i) = gr(i)/(4*3.1415*dr*i*dr*i*dr*avg_rho*caller*n_part) 
        write(18,*) dr*i, gr(i)
    enddo

    !--------------------speed distribution:
    dist=dist/(n_part*dv)
    open(unit=72,file='vel_dist.dat')
    do i=1,size(dist)
        write(72,*) dist(i)
    enddo
    
end program molecular_dynamics