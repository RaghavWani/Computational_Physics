program ising_model
    implicit none
      integer,allocatable,dimension(:,:,:) :: spin !for 3d 
      integer :: int_spin, T_temp, n_equil !initial value of spin (up or down)
      integer :: L, niter !lattice size, number of iterations
      real :: Eng,M,mag !energy, Magnetization, Magnetization per spin
      real*8 :: T, J_ising=1.0d0 !Temperature, Ising hamiltonian constant
      real*8 :: cv,chi, av_M,av_E !Specific heat, Chi, avg magnetization, avg energy
      real*8 :: av_M2,av_E2, av_M_N,av_E_N !  (avg magnetization)^2 , (avg energy)^2

      !dummy variable
      integer :: i,j,k,N,a,b,c,d,e,f,time,mm,nn,oo 
      real*8::Ei,Ef,dE,h,u
      real :: r
    
    ! Spin initialization:
    print*,'Enter the lattice size'
    read*, L
    print*,'Initial spin of entire lattice (+1 or -1)'
    read*, int_spin
    
    allocate(spin(L,L,L))
    Eng=0.0d0
    M=0.0d0
    N=L*L*L !for 3d
  
    !Initialize the Lattice
    do k=1,L !z-axis
      do j=1,L !y axis
        do i=1,L !x axis
          spin(i,j,k)=int_spin !if initially all spin up or down 
        end do
      end do
    end do
  
    !Calculate initial magnetization and energy
    do k=1,L
      do i=1,L
        do j=1,L
  
          !some constants
          a=i+1; b=i-1; 
          c=j+1; d=j-1; 
          e=k+1; f=k-1
          
          !Periodic Boundary Conditions:
          if(i==L) a=1; if(i==1) b=L !along x
          if(j==1) d=L; if(j==L) c=1 !along y
          if(k==L) e=1; if(k==1) f=L !along z
          
          M = M + spin(i,j,k) !magnetization = sum of all the spins
          !energy = J*S1*(sum of neighbouring spins of S1)
          Eng = (Eng) - J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+spin(i,j,f)))
        end do
      end do
    end do
  
    mag = M/(float(N)) !magnetization per spin
    Eng= real(Eng*0.5d0) !energy correction for twice calculation of energy
  
    print*,'initial energy E and energy per spin =', Eng, Eng/float(N)
    print*,'initial magnetization M and magnetization per spin =', M, mag
  
    !Initialization Complete..

    !Now evolve it reach equillibrium, for different temperature values
    open(1,file='q3_.dat')
    write(1,*) 'MCS step','Magnetization per spin','Energy per spin' 
    open(2,file='q7_10_L9.dat')
    write(2,*) 'T', 'av_m', 'av_e', 'cv', 'chi' 

    do T_temp=470,380,-2 !temperature loop

        !T=dfloat(T_temp)/10.0d0 !fix T
        T=T_temp/100.0d0 !fix T
        av_M = 0.0d0; av_E= 0.0d0
        av_M2= 0.0d0; av_E2= 0.0d0 
        av_M_N= 0.0d0; av_E_N= 0.0d0

        do time=1,niter !loop over no of Monte Carlo Steps
            do oo=1,L
                do mm=1,L
                    do nn=1,L
                
                        !random number between 0 to L+1
                        call random_number(r); i=int(r*float(L))+1
                        call random_number(r); j=int(r*float(L))+1
                        call random_number(r); k=int(r*float(L))+1

                        a=i+1; b=i-1; 
                        c=j+1; d=j-1; 
                        e=k+1; f=k-1

                        !!Periodic Boundary Conditions:
                        if(i==L) a=1; if(i==1) b=L !along x
                        if(j==1) d=L; if(j==L) c=1 !along y
                        if(k==L) e=1; if(k==1) f=L !along z
                        
                        !initial energy:
                        Ei=-J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+spin(i,j,f)))
                        
                        !flip the spin
                        spin(i,j,k)=-spin(i,j,k)

                        !energy after spin flip:
                        Ef=-J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+spin(i,j,f)))
                        
                        dE=Ef-Ei !energy difference after changing spin
                        !accept spin flip if energy is reduced, orelse accept it with probability exp(-dE/T)
                        if(dE<=0.0d0)then
                            E=E+dE
                            M=M+(2.0d0*float(spin(i,j,k)))
                        else
                            u=exp(-dE/T)
                            call random_number(h)
                            if(h<u)then
                                E=E+dE
                                M=M+(2.0d0*float(spin(i,j,k)))
                            else
                                spin(i,j,k)=-spin(i,j,k)
                            end if
                        end if

                    end do !do nn=1,L
                end do !do mm=1,L
            end do !do oo=1,L

            write(1,*) time,M/float(N),E/float(N)
        end do !do time=1,niter

        print*, 'doing cv, chi for', T
        av_m = av_m/dfloat(niter-n_equil) !dividing by number of microstates over which we have calculated the avg
        av_e = av_e/dfloat(niter-n_equil)
        
        cv=(av_e2 - av_e_n*av_e_n)/(T*T) !Cv= del(E)/del(T) 
        cv = cv/dfloat(niter-n_equil)

        chi=(av_m2 - av_m_n*av_m_n)/(T) !Cv= del(M)/del(B)
        chi = chi/dfloat(niter-n_equil)

        write(2,*)T, av_m, av_e, cv, chi
    end do !do T_temp=470,380,-2

    close(1)
    deallocate(spin)
    
  end program ising_model
  