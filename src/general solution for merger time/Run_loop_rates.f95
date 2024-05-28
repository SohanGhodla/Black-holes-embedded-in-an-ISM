program run_loop

    implicit none
    real(kind=8) rate_at_t_Z(301, 13), rate_at_z(301), W(81,301,10,13), sample_rates(41,301,10,13), SFRD(13427, 13), &
            rate_at_t_Z_p(301,13,10)

    real(kind=8) Lookbacktime(13427), Redshifts(13427)
    real(kind=8) t, t_z, t_z_0, dt, del_t_old, del_t_new
    real(kind=8) z_merger, z_form, T_p, a, t_new, Mtotal, Mtotal_new_l, Mtotal_new_r
    integer itime_new, iMtot_l, iMtot_r

    integer, parameter :: num_times = 301 ! number of elements in elaytime array
    real(kind=8), dimension(num_times) :: delaytime
    real(kind=8) :: time_val = 6.0d0

    character(len=12) :: t_z_string, t_z_0_string, z_string
    integer redshift, itime, Z, iratio, stat, stat2, redshift_BBH_form
    integer i, j ,k ,l, iMtot, iMtot_new

    ! ---------------- ISM and BBHs -----------------
    real(kind=8) :: c, msun, G, rand_num
    real(kind=8) :: m_H, c_s, yr, rho, alpha, n_H, delta_t, Mtot, Mtot_new, m1, m2, mass_ratio


    if (command_argument_count().ne.3) then
        ! write(*,*) "ERROR: 2 - STOPPING"
        stop
    end if

    call get_command_argument(1,t_z_string)
    call get_command_argument(2,t_z_0_string)
    call get_command_argument(3,z_string)

    read(t_z_string, *) t_z
    read(t_z_0_string, *) t_z_0
    read(z_string,*) z_merger


    do i = 1, num_times
        delaytime(i) = time_val
        time_val = time_val + 0.1d0
    end do


    ! ---------------- ISM and BBHs -----------------
    ! ---- Constants ----
    c = 3.0d8
    G = 6.67d-11
    yr = 365.25d0 * 24.0d0 * 3600.0d0
    msun = 1.99d30
    m_H = 1.671d-27  ! kg
    n_H = 0.01d0

    ! ---- Baryonic matter type ----
    c_s = sqrt(1.0d0 * 1.38d-23 * 100.0d0 / m_H)
    rho = n_H * m_H / (1.0d-2)**3  ! kg/m^3
    alpha = 4.0d0 * (4 * atan(1.0d0)) * 1.12d0 * G**2d0 * rho / c_s**3d0


    sample_rates(1:41,1:301,1:10, 1:13) = 0d0
    rate_at_t_Z(1:301, 1:13)            = 0d0
    rate_at_z(1:301)                    = 0d0
    W(1:81,1:301,1:10, 1:13)            = 0d0
    SFRD(1:13427, 1:13)                 = 0d0
    rate_at_t_Z_p(1:301,1:13,1:10)      = 0d0

    stat = 0
    stat2 = 0
    open(8,FILE='Lookbacktime.txt',STATUS='old',ACCESS='SEQUENTIAL',ACTION='read')
    read(8,*,iostat=stat) (Lookbacktime(i), i=1,13427)
    close(8)

    open(8,FILE='Redshifts.txt',STATUS='old',ACCESS='SEQUENTIAL',ACTION='read')
    read(8,*,iostat=stat) (Redshifts(i), i=1,13427)
    close(8)

    open(9,FILE='sample_rates.txt',STATUS='old',ACCESS='SEQUENTIAL',ACTION='read')
    do i = 1,13
        do j = 1,10
            do k = 1,301
                read(9,*, iostat = stat2) (sample_rates(l, k ,j ,i), l=1,41)
            enddo
        enddo
    enddo
    close(9)


    open(9,FILE='SFRD.txt',STATUS='old',ACCESS='SEQUENTIAL',ACTION='read')
    stat2 = 0
    do i = 1,13
        read(9,*, iostat = stat2) (SFRD(l,i), l=1,13427)
    enddo
    close(9)


    t = 200d6
    dt = 5d6
    ! print*, t_z, t_z_0
    do while (t < t_z)  ! t_z is the age of universe at merger redshift z.

        redshift  = minloc(abs(Lookbacktime - ((t_z_0 - t)/1d9) ), 1)  !find_nearest(Lookbacktime, (t_z_0 - t)/1e9)

        redshift_BBH_form  = minloc(abs(Lookbacktime - ((t_z_0 - t + 1d7)/1d9) ), 1)   ! --- Both redshift and redshift_at_BBHform are same
        z_form = Redshifts(redshift_BBH_form)

        ! ------------ ISM and BBHs ------------
        delta_t = t_z - t - 1d7  ! The latter term to account for BBH formation

        rho = n_H * m_H / (1.0d-2)**3 * (1+z_form)**3d0 !/ Mtot  ! kg/m^3
        alpha = 4.0d0 * (4 * atan(1.0d0)) * 1.12d0 * G**2d0 * rho / c_s**3d0

        do itime = 1,200 !301  ! these are old merger time bin index - need to associate new time with these bin index.
            ! t_new = 10**(delaytime(itime)) ! just initially

            do iMtot = 10,35 !1,40 ! not considering the other masses as they are either too low or high (saves time)
                do iratio = 1,10
                    T_p =   10**(delaytime(itime))  ! in years
                    Mtot = 10**((iMtot - 10) / 10d0)

                    mass_ratio = 0.1*iratio
                    m1 = Mtot/(mass_ratio + 1d0)
                    m2 = Mtot - m1


                    call random_number(rand_num)
                    rand_num = FLOOR(100*rand_num + 1) ! gives num b/w 1-10
                    if (rand_num <= 3 ) then
                        call calc_T(t_new, T_p, m1, m2, alpha)  ! to update t to new t merger (coupled BHs)
                        itime_new = 1+nint(10d0*(log10(t_new)-6d0))  ! restructuring the delaytime bins so that the max value is again in the 42 bin (age of the universe)

                        ! To put bounds on time -- don't want to over calculate
                        if(itime_new <= 1) itime_new=1
                        if(itime_new>=301) itime_new=301

                        Mtot_new = m1 * (-1d0 / (m1 * msun * alpha * t_new * yr - 1d0))  + &  ! mass growing seperatly
                                m2 * (-1d0 / (m2 * msun * alpha * t_new * yr - 1d0))
                        iMtot_new = nint(log10(Mtot_new)*10d0)+10
                        if(iMtot_new<=1) iMtot_new = 1
                        if(iMtot_new>=80) iMtot_new = 80

                       
                    else
                        t_new = T_p

                        itime_new = itime  ! restructuring the delaytime bins so that the max value is again in the 42 bin (age of the universe)
                        Mtot_new = Mtot
                    endif

                       
                    if(t_new <= delta_t) then
                        if((delta_t) <= 10**(delaytime(itime_new)).and.itime_new == 1) then ! do not change the 1d6 as it very small anyway

                            !                        do Z = 1,13
                            W(iMtot_new + 1, 1, iratio, 1:13) =  W(iMtot_new + 1, 1, iratio, 1:13) &
                                    + sample_rates(iMtot + 1, itime, iratio, 1:13) * SFRD(redshift, 1:13) * dt
                            !W[Z][p].iloc[1, 2:41] += sample_rates[Z][p].iloc[1, 2:41] *  SFRD[Z][redshift] * dt
                            !                            rate_at_t_Z_p(1,Z, iratio)= sum(W(2:81,1,iratio, Z))
                            !                            rate_at_t_Z(1, Z) = sum(rate_at_t_Z_p(1, Z, 1:10))
                            !                        enddo
                            !                        rate_at_z(1) = sum(rate_at_t_Z(1, 1:13))


                        else if ( ((delta_t) >= 10**(delaytime(itime_new) - 0.05d0)) &
                                .and. ((delta_t) < 10**(delaytime(itime_new) + 0.05d0)) ) then

                            del_t_old = (((10**(5.95d0 + 0.1d0 *(itime)))*(10**(0.1d0) - 1)))
                            del_t_new = (((10**(5.95d0 + 0.1d0 *(itime_new)))*(10**(0.1d0) - 1)))

                            !                        do Z = 1,13
                            W(iMtot_new+1, itime_new, iratio, 1:13)=W(iMtot_new+1, itime_new, iratio, 1:13)&
                                    + sample_rates(iMtot + 1, itime, iratio, 1:13)&
                                            * SFRD(redshift, 1:13) * dt * (del_t_old/del_t_new)
                            !W[Z][p].iloc[itime, 2:41] += sample_rates[Z][p].iloc[itime, 2:41] *  SFRD[Z][redshift] * dt
                            !                            rate_at_t_Z_p(itime_new,Z, iratio) = sum(W(2:81,itime_new,iratio, Z))
                            !                            rate_at_t_Z(itime_new, Z) = sum(rate_at_t_Z_p(itime_new, Z, 1:10))
                            !                        enddo
                            !                        rate_at_z(itime_new) = sum(rate_at_t_Z(itime_new, 1:13))  !!! the eventually left over itime bins remain zero as their matter did not end up in mergering at t_z. Note
                            ! that the W array does store the data from the previous instance of the loop even though the rate_at_z array does not.
                        endif
                    endif
                                                                                                
                enddo
            enddo
        enddo

        ! endif
        !    print*,t, t_z
        t  = t + dt !10**(dt - (1e-2)*index)
    enddo


    open(10,FILE='./files/W_'//trim(z_string)//".txt",STATUS='unknown',ACCESS='SEQUENTIAL',ACTION='write')
    do itime = 1,301
        do Z=1,13
            do iratio=1,10
                rate_at_t_Z_p(itime,Z, iratio) = sum(W(2:81,itime,iratio, Z))

                write(10,*) (W(iMtot, itime, iratio,  Z), iMtot=1,81)
            end do
            rate_at_t_Z(itime, Z) = sum(rate_at_t_Z_p(itime, Z, 1:10))
        end do
        rate_at_z(itime) = sum(rate_at_t_Z(itime, 1:13))
    end do
    close(10)

    open(11,FILE='./files/rate_at_z_'//trim(z_string)//".txt",STATUS='unknown',ACCESS='SEQUENTIAL',ACTION='write')
    write(11, *) (rate_at_z(i), i=1,301)
    close(11)
    !
    !open(10,FILE='./files/W_'//trim(z_string)//".txt",STATUS='unknown',ACCESS='SEQUENTIAL',ACTION='write')
    !do Z=1,13
    !    do iratio = 1,10
    !        do itime = 1,301
    !            write(10,*) (W(iMtot, itime, iratio,  Z), iMtot=1,81)
    !        enddo
    !    enddo
    !enddo
    !close(10)

end program run_loop