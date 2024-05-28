subroutine calc_T(T_new, T_p, m1, m2, alpha)  !T_p should be in years, retuns T in years
  implicit none

  real(kind=8) :: T_p, l, r, m, T_new, m1, m2, alpha, msun, sec_yr
  integer :: i

  msun = 1.99d30
  l = 1d0
  sec_yr = 365.25d0*24d0*3600d0
  r = 13.8d9 * sec_yr
  ! T_p = 1d7
  ! z = 10d0

T_p = T_p * sec_yr
if(m2/m1 <= 0.8) then
  do i = 1, 30  ! performing a binary search
    m = (l+r)/2d0
    if (isnan(equation(m,  T_p, m1, m2, alpha)) .or. (l > 0.95*r) .or. (r < 1.05*l))  then
      T_new = 1d0 / ((m1+m2) * msun /2d0 * alpha) * (- ( 1d0 / (14d0 * (m1+ &
            m2) * msun /2d0 * alpha * T_p + 1d0)) ** (1d0 / 14d0) + 1d0)
      T_new = T_new / sec_yr
      exit
    
    else if (equation(m,  T_p, m1, m2, alpha) > 0d0) then
      r = m
    else
      l = m
    end if
    T_new = l / sec_yr   ! returns T in years
  !  write(*,*) T_new
  end do

else
  T_new = 1d0 / ((m1+m2) * msun /2d0 * alpha) * (- ( 1d0 / (14d0 * (m1+ &
          m2) * msun /2d0 * alpha * T_p + 1d0)) ** (1d0 / 14d0) + 1d0)
  T_new = T_new / sec_yr
end if


contains
real(kind=8) function equation(T,  T_p, m1, m2, alpha)
real(kind=8) :: T, T_p, m1, m2, alpha, msun, c, d, trunc_log_cT, trunc_log_dT, temp

msun = 1.99d30
c = m1 * msun * alpha
d = m2 * msun * alpha

temp = -300d0
if (isnan(log(1-c*T)) .or.  (log(1d0-c*T) < temp)) then
  trunc_log_cT = log(10**temp)
else
  trunc_log_cT = log(1d0-c*T)
endif

if (isnan(log(1-d*T)) .or.  (log(1d0-d*T) < temp)) then
  trunc_log_dT = log(10**temp)
else
  trunc_log_dT = log(1d0-d*T)
endif
!
equation = -T_p - ((((c - d)*T*(126d0*c**16d0*T**3d0*(-1d0 + d*T)**9d0 - 42d0*c**15d0*T**2d0*(-1d0 + d*T)**9d0* &
          (12d0 + 49d0*d*T) + 42d0*c**14d0*T*(-1d0 + d*T)**9d0*(18d0 + 196d0*d*T + &
          419d0*d**2d0*T**2d0) - 42d0*c**13d0*(-1d0 + d*T)**9d0*(12d0 + 294d0*d*T + 1676d0*d**2d0*T**2d0 + &
          3013d0*d**3d0*T**3d0) - 56d0*d**13d0*(9d0 - 36d0*d*T + 84d0*d**2d0*T**2d0 - 126d0*d**3d0*T**3d0 + &
          126d0*d**4d0*T**4d0 - 84d0*d**5d0*T**5d0 + 36d0*d**6d0*T**6d0 - 9d0*d**7d0*T**7d0 + d**8d0*T**8d0) + &
          7d0*c*d**12d0*(936d0 - 3636d0*d*T + 8004d0*d**2d0*T**2d0 - 11046d0*d**3d0*T**3d0 + &
                  9702d0*d**4d0*T**4d0 - 5124d0*d**5d0*T**5d0 + 1236d0*d**6d0*T**6d0 + 171d0*d**7d0*T**7d0 - &
                  179d0*d**8d0*T**8d0 + 32d0*d**9d0*T**9d0) + c**2d0*d**11d0*(-39312d0 + 147420d0*d*T - &
          300804d0*d**2d0*T**2d0 + 367206d0*d**3d0*T**3d0 - 255150d0*d**4d0*T**4d0 + 65436d0*d**5d0*T**5d0 + &
          39324d0*d**6d0*T**6d0 - 37971d0*d**7d0*T**7d0 + 10519d0*d**8d0*T**8d0 - 28d0*d**9d0*T**9d0 - &
          336d0*d**10d0*T**10d0) + c**3d0*d**10d0*(144144d0 - 517608d0*d*T + 952224d0*d**2d0*T**2d0 - &
          947730d0*d**3d0*T**3d0 + 346122d0*d**4d0*T**4d0 + 277536d0*d**5d0*T**5d0 - 391848d0*d**6d0*T**6d0 + &
          169965d0*d**7d0*T**7d0 - 9545d0*d**8d0*T**8d0 - 13726d0*d**9d0*T**9d0 + 2562d0*d**10d0*T**10d0 + &
          224d0*d**11d0*T**11d0) + c**9d0*d**5d0*T*(2594592d0 - 13837824d0*d*T + &
          31423392d0*d**2d0*T**2d0 - 36756720d0*d**3d0*T**3d0 + 20396376d0*d**4d0*T**4d0 - &
          2234232d0*d**5d0*T**5d0 + 332046d0*d**6d0*T**6d0 - 5255250d0*d**7d0*T**7d0 + &
          4499066d0*d**8d0*T**8d0 - 991146d0*d**9d0*T**9d0 - 267241d0*d**10d0*T**10d0 + &
          109261d0*d**11d0*T**11d0) + c**12d0*d*(-6552d0 - 32760d0*d*T + 194376d0*d**2d0*T**2d0 + &
          805896d0*d**3d0*T**3d0 - 7999992d0*d**4d0*T**4d0 + 27230112d0*d**5d0*T**5d0 - &
          54224352d0*d**6d0*T**6d0 + 71149806d0*d**7d0*T**7d0 - 63818118d0*d**8d0*T**8d0 + &
          39009516d0*d**9d0*T**9d0 - 15615756d0*d**10d0*T**10d0 + 3703011d0*d**11d0*T**11d0 - &
          395243d0*d**12d0*T**12d0) + c**10d0*d**3d0*(-144144d0 + 1297296d0*d*T - &
          9081072d0*d**2d0*T**2d0 + 37189152d0*d**3d0*T**3d0 - 89080992d0*d**4d0*T**4d0 + &
          130666536d0*d**5d0*T**5d0 - 118126008d0*d**6d0*T**6d0 + 61675614d0*d**7d0*T**7d0 - &
          13960518d0*d**8d0*T**8d0 - 1164306d0*d**9d0*T**9d0 + 306150d0*d**10d0*T**10d0 + &
          597623d0*d**11d0*T**11d0 - 179027d0*d**12d0*T**12d0) + &
          c**8d0*d**5d0*(-648648d0 + 1297296d0*d*T + 1873872d0*d**2d0*T**2d0 - 8738730d0*d**3d0*T**3d0 + &
                  8594586d0*d**4d0*T**4d0 + 3531528d0*d**5d0*T**5d0 - 13168584d0*d**6d0*T**6d0 + &
                  8822385d0*d**7d0*T**7d0 + 103675d0*d**8d0*T**8d0 - 2636062d0*d**9d0*T**9d0 + &
                  954798d0*d**10d0*T**10d0 + 39065d0*d**11d0*T**11d0 - 52901d0*d**12d0*T**12d0) - &
          c**4d0*d**9d0*(360360d0 - 1225224d0*d*T + 1937208d0*d**2d0*T**2d0 - 1244334d0*d**3d0*T**3d0 - &
                  660114d0*d**4d0*T**4d0 + 1803984d0*d**5d0*T**5d0 - 1248000d0*d**6d0*T**6d0 + &
                  216099d0*d**7d0*T**7d0 + 161681d0*d**8d0*T**8d0 - 80210d0*d**9d0*T**9d0 + 3666d0*d**10d0*T**10d0 + &
                  2548d0*d**11d0*T**11d0 + 56d0*d**12d0*T**12d0) + c**5d0*d**8d0*(648648d0 - 2054052d0*d*T + &
          2534532d0*d**2d0*T**2d0 + 51870d0*d**3d0*T**3d0 - 3926286d0*d**4d0*T**4d0 + &
          4594044d0*d**5d0*T**5d0 - 1729884d0*d**6d0*T**6d0 - 678015d0*d**7d0*T**7d0 + &
          811291d0*d**8d0*T**8d0 - 188890d0*d**9d0*T**9d0 - 32214d0*d**10d0*T**10d0 + &
          12545d0*d**11d0*T**11d0 + 763d0*d**12d0*T**12d0) - c**6d0*d**7d0*(864864d0 - 2486484d0*d*T + &
          1837836d0*d**2d0*T**2d0 + 3525522d0*d**3d0*T**3d0 - 8594586d0*d**4d0*T**4d0 + &
          6306300d0*d**5d0*T**5d0 + 471900d0*d**6d0*T**6d0 - 3489057d0*d**7d0*T**7d0 + &
          1842269d0*d**8d0*T**8d0 - 66638d0*d**9d0*T**9d0 - 198042d0*d**10d0*T**10d0 + &
          33007d0*d**11d0*T**11d0 + 4853d0*d**12d0*T**12d0) + c**7d0*d**6d0*(864864d0 - 2162160d0*d*T - &
          72072d0*d**2d0*T**2d0 + 7405398d0*d**3d0*T**3d0 - 11081070d0*d**4d0*T**4d0 + &
          3531528d0*d**5d0*T**5d0 + 6507072d0*d**6d0*T**6d0 - 7321743d0*d**7d0*T**7d0 + &
          2049619d0*d**8d0*T**8d0 + 823394d0*d**9d0*T**9d0 - 558714d0*d**10d0*T**10d0 + &
          39065d0*d**11d0*T**11d0 + 19171d0*d**12d0*T**12d0) + &
          c**11d0*d**2d0*(39312d0 + 58968d0*d*T - 2299752d0*d**2d0*T**2d0 + 14152320d0*d**3d0*T**3d0 - &
                  47882016d0*d**4d0*T**4d0 + 103128480d0*d**5d0*T**5d0 - 147590352d0*d**6d0*T**6d0 + &
                  140987574d0*d**7d0*T**7d0 - 87184734d0*d**8d0*T**8d0 + 31641324d0*d**9d0*T**9d0 - &
                  4604184d0*d**10d0*T**10d0 - 699673d0*d**11d0*T**11d0 + 253405d0*d**12d0*T**12d0)))/&
          ((-1d0 + c*T)**4d0*(-1d0 + d*T)**9d0) + 360360d0*c**9d0*d**4d0*trunc_log_cT - &
          360360d0*c**9d0*d**4d0*trunc_log_dT)/(504d0*(c - d)**14d0))

  end function equation

end subroutine 