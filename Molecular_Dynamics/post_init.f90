subroutine post_init
    !initializing position of particles randomly but uniformly over box size
    
    use all_parameters
    implicit none
        real*8 :: r,x1,y1,z1,dx,dy,dz,x2,y2,z2
        integer :: i,j!,k,N_L,q
    
    open(2, file='initial_position.dat',action='write') !file saving initial velocity
    write(2,*) 'x, y, z'

    do i=1,n_part
25      CALL random_number(r)
        post(3*i-2) = r*lx

        CALL random_number(r)
        post(3*i-1) = r*ly

        CALL random_number(r)
        post(3*i)   = r*lz

        !To ensure that particles don't overlap, we consider particles that are atleast 1 sigma far
        if (i>=2) then 
            do j=i-1,1,-1
                x1=post(3*i-2); y1=post(3*i-1); z1=post(3*i)
                x2=post(3*j-2); y2=post(3*j-1); z2=post(3*j)
                dx=(x2-x1); dy=(y2-y1); dz=(z2-z1)

                dx = dx - nint(dx/lx)*lx
                dy = dy - nint(dy/ly)*ly
                dz = dz - nint(dz/lz)*lz

                r=dsqrt((dx)**2 + (dy)**2 + (dz)**2)
                if (r.le.sigma) then
                    goto 25
                end if
            end do
        end if
    end do
    
    do i=1,n_part
        write(2,*) post(3*i-2), post(3*i-1), post(3*i)
    end do
    close(2)
    
end subroutine post_init