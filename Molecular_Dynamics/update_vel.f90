subroutine update_vel
    use all_parameters
    implicit none
    real*8 :: avg_vx, avg_vy, avg_vz
    integer :: i

    !Using Verlet algorithm:
    ! v(t+dt)=v(t) + 0.5*dt*( f(t) + f(t+dt) )/mass
    
    avg_vx=0.0d0; avg_vy=0.0d0; avg_vz=0.0d0
    KE= 0.0d0

    do i=1,n_part
        velocity(3*i-2) = velocity(3*i-2) + dtby2*(old_force(3*i-2) + force(3*i-2))
        velocity(3*i-1) = velocity(3*i-1) + dtby2*(old_force(3*i-1) + force(3*i-1))
        velocity(3*i)   = velocity(3*i)   + dtby2*(old_force(3*i)   + force(3*i))

        !Average velocity:
        avg_vx = avg_vx + velocity(3*i-2)
        avg_vy = avg_vy + velocity(3*i-1)
        avg_vz = avg_vz + velocity(3*i)

        KE = KE + ((velocity(3*i-2))**2 + (velocity(3*i-1))**2 + (velocity(3*i))**2) !Kinetic energy
    end do
    
    KE = 0.50d0*dfloat(mass)*KE
end subroutine update_vel