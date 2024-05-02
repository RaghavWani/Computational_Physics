subroutine calc_force_neighb
    use all_parameters
    implicit none
    real*8 :: x1,y1,z1,dx,dy,dz,r_12
    integer :: jj,kk,p

    potential_energy = 0.0d0
    force=0.0d0

    do jj=1, (n_part - 1)
        x1 = post(3*jj - 2); y1 = post(3*jj - 1); z1 = post(3*jj) !position of (jj)th particle 

        do kk= 1, no_neigh(jj) !positions of particle other than (jj)th one 
            p=neigh_list(jj,kk) !(kk)th neighbour of particle (jj)

            !As we are considering PBC, we are accounting for minimum image convention here (if statements):
            dx=x1-post(3*p - 2) !dist along x
            dy=y1-post(3*p - 1) !dist along y
            dz=z1-post(3*p) !dist along z

            if (abs(dx).ge.(dfloat(lx)/2.0d0)) dx = (dfloat(lx) -abs(dx))*(-1.0d0*dx/abs(dx))
            if (abs(dy).ge.(dfloat(ly)/2.0d0)) dy = (dfloat(ly) -abs(dy))*(-1.0d0*dy/abs(dy))
            if (abs(dz).ge.(dfloat(lz)/2.0d0)) dz = (dfloat(lz) -abs(dz))*(-1.0d0*dz/abs(dz))

            r_12 = dsqrt(dx**2 + dy**2 + dz**2) !distance between (jj)th and (kk)th particle 
            
            !consider only those particles for force calculation that lie at distances within r_c along each direction.
            !as force and potential is zero beyond r_c
            if ((r_12.le.r_c)) then
                ! write(*,*) r_12
                force_lj = eps*((12.0d0*(sig_12/(r_12**(13)))) - (6.0d0*(sig_6/(r_12**(7))))) - Fc !Lennard-Jones force
                
                !force on (jj)th particle along each direction:
                force(3*jj -2) = force(3*jj -2) + force_lj*(dx/(r_12)) !along x
                force(3*jj -1) = force(3*jj -1) + force_lj*(dy/(r_12)) !along y
                force(3*jj)    = force(3*jj)    + force_lj*(dz/(r_12)) !along z
                
                !force on other particle in opposite direction:
                force(3*p -2) = force(3*p -2) - force_lj*(dx/(r_12)) !along x
                force(3*p -1) = force(3*p -1) - force_lj*(dy/(r_12)) !along y
                force(3*p)    = force(3*p)    - force_lj*(dz/(r_12)) !along z

                potential_lj = eps*((sigma/r_12)**12 - (sigma/r_12)**6) + Fc*r_12 - Vc !Lennard-Jones potential 
                potential_energy = potential_energy + potential_lj
                                                    
            end if
        end do
    end do
end subroutine calc_force_neighb