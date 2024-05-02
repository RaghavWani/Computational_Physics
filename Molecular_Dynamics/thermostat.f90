subroutine thermostat
    use all_parameters
    implicit none
        real*8  :: KE_actual,Sf
        integer :: i
    KE =0.0d0
    do i=1,n_part
        KE = KE + ((velocity(3*i-2))**2 + (velocity(3*i-1))**2 + (velocity(3*i))**2) !Kinetic energy
    end do
    KE = 0.50d0*dfloat(mass)*KE

    !rescaling kinetic energy
    KE_actual = 1.50d0*KbT*dfloat(n_part)
    Sf = dsqrt(KE_actual/ KE) !scaling factor

    do i=1,n_part
        velocity(3*i - 2) = velocity(3*i - 2)*Sf
        velocity(3*i - 1) = velocity(3*i - 1)*Sf
        velocity(3*i)     = velocity(3*i)*Sf

        KE = KE + ((velocity(3*i-2))**2 + (velocity(3*i-1))**2 + (velocity(3*i))**2) !Kinetic energy
    end do
    KE = 0.50d0*dfloat(mass)*KE
end subroutine thermostat