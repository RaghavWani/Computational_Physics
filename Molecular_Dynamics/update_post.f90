subroutine update_post
    use all_parameters
    implicit none
    integer :: i
    !

    !Using Verlet algorithm:
    ! x(t+dt)=x(t) + v(t)*dt + (0.5*dt*f(t))/mass
    ! post=modulo(post,dfloat(lx))

    do i=1,n_part
        post(3*i - 2) = post(3*i - 2) + velocity(3*i - 2)*dt + dt2by2*force(3*i - 2)
        post(3*i - 1) = post(3*i - 1) + velocity(3*i - 1)*dt + dt2by2*force(3*i - 1)
        post(3*i)     = post(3*i)     + velocity(3*i)*dt     + dt2by2*force(3*i)

        !Periodic Boundary Conditions:
        post(3*i - 2) = modulo(post(3*i - 2),dfloat(lx))
        post(3*i - 1) = modulo(post(3*i - 1),dfloat(ly))
        post(3*i)     = modulo(post(3*i),dfloat(lz))
    end do
    
end subroutine update_post