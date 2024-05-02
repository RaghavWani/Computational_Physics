subroutine vel_init
    use all_parameters
    implicit none
    real*8 :: avg_vx, avg_vy, avg_vz
    real*8 :: r,A
    integer :: i

    KE=0.0d0
    open(1, file='initial_velocity.dat',action='write') !file saving initial velocity
    write(1,*) 'vel_x, vel_y, vel_z'

    ! velocity = sqrt(12*KbT/mass)*(random number between -0.5 0.5) along each direction
    A=dsqrt(12.0d0)*(KbT/mass)

    do i=1,n_part
        CALL random_number(r)
        velocity(3*i-2) = A*(r - 0.5d0) !along x

        CALL random_number(r)
        velocity(3*i-1) = A*(r - 0.5d0) !along y

        CALL random_number(r)
        velocity(3*i)   = A*(r - 0.5d0) !along z
    end do
    
    !average vel along each direction
    avg_vx=0.0d0; avg_vy=0.0d0; avg_vz=0.0d0
    
    do i=1,n_part
        avg_vx=velocity(3*i-2) + avg_vx 
        avg_vy=velocity(3*i-1) + avg_vy
        avg_vz=velocity(3*i) + avg_vz
    end do

    avg_vx=avg_vx/dfloat(n_part)
    avg_vy=avg_vy/dfloat(n_part)
    avg_vz=avg_vz/dfloat(n_part)

    write(*,*) avg_vx, avg_vy, avg_vz

    !velocity correction - frame shifting - all velocity wiyh respect to Center of Mass frame of the system
    do i=1,n_part
        velocity(3*i-2) = avg_vx - velocity(3*i-2) !along x
        velocity(3*i-1) = avg_vy - velocity(3*i-1) !along y
        velocity(3*i)   = avg_vz - velocity(3*i)   !along z

        KE = KE + ((velocity(3*i-2))**2 + (velocity(3*i-1))**2 + (velocity(3*i))**2)
    end do
    KE = 0.50d0*dfloat(mass)*KE
    
    do i=1,n_part
        write(1,*) velocity(3*i-2), velocity(3*i-1), velocity(3*i)
    end do
    close(1)

end subroutine vel_init