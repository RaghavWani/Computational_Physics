subroutine calc_gr
    use all_parameters
    implicit none
    real*8 :: x1,y1,z1,dx,dy,dz,r_12
    integer :: jj,kk,loc

    do jj=1, n_part
        x1 = post(3*jj - 2); y1 = post(3*jj - 1); z1 = post(3*jj) !position of (jj)th particle 

        do kk= 1, n_part !positions of particle other than (jj)th one 
            if (jj/=kk) then
                !As we are considering PBC, we are accounting for minimum image convention here (if statements):
                dx=x1-post(3*kk - 2) !dist along x
                dy=y1-post(3*kk - 1) !dist along y
                dz=z1-post(3*kk) !dist along z

                if (abs(dx).ge.(dfloat(lx)/2.0d0)) dx = (dfloat(lx) -abs(dx))*(-1.0d0*dx/abs(dx))
                if (abs(dy).ge.(dfloat(ly)/2.0d0)) dy = (dfloat(ly) -abs(dy))*(-1.0d0*dy/abs(dy))
                if (abs(dz).ge.(dfloat(lz)/2.0d0)) dz = (dfloat(lz) -abs(dz))*(-1.0d0*dz/abs(dz))

                r_12 = dsqrt(dx**2 + dy**2 + dz**2) !distance between (jj)th and (kk)th particle 
                
                !pair correlation:
                loc=int(r_12/dr) + 1
                gr(loc) = gr(loc) + 1
            end if
        end do
    end do
    
end subroutine calc_gr