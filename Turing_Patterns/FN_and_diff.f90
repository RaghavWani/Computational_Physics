!Fitzhugh-Nagumo + diffusion equation
!dA/dt= D_a*d^2(A) + R1; where R1=(A - A^3 - B + alpha)
!dB/dt= D_b*d^2(B) + R2; where R2=beta*(A - B)

program diffusion3
    implicit none
    integer, parameter :: lx=60, ly=60! BOUNDARIES at x,y =0 and x,y=33. 
    real*8 :: old_tempA(1:lx, 1:ly), tempA(1:lx, 1:ly)
    real*8 :: old_tempB(1:lx, 1:ly), tempB(1:lx, 1:ly)
    integer :: jj, ii
    real*8 ::  prefactor, densityA, densityB, limit!,init_temp_high,init_temp_low,  
    integer :: xp, xn, yp, yn !dummy
    character (len=30) :: filename,filenamel
    integer :: iunit, iunitB, test, counter
    
    !parametrs for simulations
    real*8 :: dx=1.0d0, dy=1.0d0, dt=0.002d0 
    !real*8 :: toler = 0.0001d0
    real*8 :: Diff_a=1.0d0, Diff_b=100.0d0 !diffusion constants
    real*8 :: alpha=0.01, beta=10.0d0 !Fitzhugh-Nagumo parametrs
    integer :: n_snapshots=1000, tot_niter=40000, niter
    
    old_tempA= 0.0d0
    old_tempB= 0.0d0

    ! old_tempA(1,1) = 0.1d0
    ! old_tempA(1,ly) = 0.1d0
    ! old_tempA(lx,ly) = 0.1d0
    ! old_tempA(lx,1) = 0.1d0
    
    ! old_tempA(lx,ly/4) = 0.1d0
    ! old_tempA(lx,ly/2) = 0.1d0
    ! old_tempA(lx,3*ly/4) = 0.1d0

    call random_number(old_tempA)
    call random_number(old_tempB)
    old_tempA=old_tempA - 0.50d0; 
    old_tempB=old_tempB - 0.50d0
          
    limit = 0.0001d0 !limit of CONVERGENCE
    
    ! iunit=71; ci=0
    ! write(filename, '("initialize_", i0,".dat")') ci
    open(unit=1, file='int_tempA.dat')
    open(unit=2, file='int_tempB.dat')
    ! WRITE DOWN THE BOUNDARY CONDITIONS at ZERoth ITERATION> 
    do ii=1, lx
        do jj=1, ly
            write(1, *) ii, jj, old_tempA(ii,jj)
            write(2, *) ii, jj, old_tempB(ii,jj)
        enddo
    enddo
    close(1)
    close(2)
        
    !dx=0.05d0; dy=0.05d0
    test=0; counter =0; prefactor = (0.5d0*dx*dx*dy*dy)/(dx*dx + dy*dy)
    do niter =1,tot_niter!LOOP OVER ITERATIONS.
        !counter = counter +1
        test=0
        
        do jj=1,ly ! UPDATE STEP: going across a row
            xp = jj+1; xn= jj-1
            
            do ii=1, lx !going down a column
                if(jj.eq.1) xn=ly
                if(jj.eq.ly) xp=1
                
                yp = ii + 1; yn = ii - 1
                
                if(ii.eq.1) yn =lx
                if(ii.eq.lx) yp=1
                
                !dA/dt= D_a*d^2(A) + R1; where R1=(A - A^3 - B + alpha)
                tempA(ii, jj) = old_tempA(ii,jj) + dt*diff_a*(old_tempA(yn,jj) + old_tempA(yp,jj) + &
                & old_tempA(ii,xp) + old_tempA(ii,xn) - 4.0d0*old_tempA(ii,jj))
                tempA(ii,jj) = tempA(ii,jj) +dt*(old_tempA(ii,jj) - (old_tempA(ii,jj))**3 + alpha - old_tempB(ii,jj)) 
                
                !dB/dt= D_b*d^2(B) + R2; where R2=beta*(A - B)
                tempB(ii, jj) = old_tempB(ii,jj) + dt*diff_b*(old_tempB(yn,jj) + old_tempB(yp,jj) + old_tempB(ii,xp) + &
            & old_tempB(ii,xn) - 4.0d0*old_tempB(ii,jj))
                tempB(ii,jj) = tempB(ii,jj) +dt*beta*(old_tempA(ii,jj) - old_tempB(ii,jj))
                
            enddo
            
        end do
        
        ! do jj=2,ly-1! CHECKING FOR CONVERGENCE AT EACH LATTICE SITE 
        !     do ii=2, lx-1
        !         if((abs (tempA(jj,ii) - old_tempA(jj,ii))).gt.limit) test=1
        !         if((abs (tempB(jj,ii) - old_tempB(jj,ii))).gt.limit) test=1
        !     enddo 
        ! enddo
        
        ! if(test.eq.0) exit! EXIT CONDITION. 
        old_tempA= tempA! AFTER TEST CONDITION
        old_tempB= tempB
    
    ! write(*,*) 'counter', counter
       
    !WRITE THE evolving equation RESULT: SOLUTION TO THE EON with PBC.
    if(mod(niter,n_snapshots).eq.0) then
        iunit=71; iunitB=72
        
        write(filename, '("initialize_", i0,".dat")') niter
        write(filenamel,'("initializeB_", i0, ".dat")') niter
        
        densityA=0.0d0; densityB=0.0d0
        open(newunit=iunit, file=filename) ! 
        open(newunit=iunitB,file=filenamel) !
        do ii=1,lx
            do jj=1,ly
                write(iunit,*) ii,jj, tempA(ii,jj)
                write(iunitB,*) ii,jj, tempB(ii,jj) 
                densityA = densityA + tempA(ii,jj) 
                densityB = densityB + tempB(ii,jj)
            enddo
        enddo
        close(iunit)
        close(iunitB)
        write(*,*) niter, densityA, densityB
    endif
    enddo ! do niter =1,10000
    
End Program diffusion3  