subroutine speed_dist
    use all_parameters
    implicit none
    integer::counter,i,j
    real*8::vx,shell_v !speed along x-direction only
    vel_arr=0
    counter=0
    do i=1,n_part
        counter=counter+1
        vx=velocity(3*i-2)
        vel_arr(counter)=vx
    enddo
    
    shell_v=0
    do i=1,size(dist)
        shell_v=shell_v + dv !speed of particles between v and v+dv
        do j=1,size(vel_arr)
            if ((shell_v-dv<vel_arr(j)).and.(vel_arr(j)<shell_v)) dist(i)=dist(i)+ 1
        enddo
    enddo    
end subroutine speed_dist