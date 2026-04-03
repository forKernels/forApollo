! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_time.f90 — Precision time system conversions
!
! Time scales: JD, MJD, Unix, UTC, TAI, TT, TDB, GPS, GMST
! Calendar ↔ JD conversions, leap second table.
!
! Future: migrates to forTime when that library ships.
! Until then, these routines are self-contained.
!
! Error codes: info = 0 ok, 3 invalid input.
! All JD values are c_double (Julian Date, days since -4712 Jan 1.5).

! ==========================================================================
! fa_time_jd_to_mjd — Julian Date to Modified Julian Date
! ==========================================================================
subroutine fa_time_jd_to_mjd(jd, mjd, info) &
    bind(C, name="fa_time_jd_to_mjd")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: jd
    real(c_double), intent(out)       :: mjd
    integer(c_int), intent(out)       :: info

    info = 0
    mjd = jd - 2400000.5d0
end subroutine

! ==========================================================================
! fa_time_mjd_to_jd — Modified Julian Date to Julian Date
! ==========================================================================
subroutine fa_time_mjd_to_jd(mjd, jd, info) &
    bind(C, name="fa_time_mjd_to_jd")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: mjd
    real(c_double), intent(out)       :: jd
    integer(c_int), intent(out)       :: info

    info = 0
    jd = mjd + 2400000.5d0
end subroutine

! ==========================================================================
! fa_time_unix_to_jd — Unix epoch seconds to Julian Date
! ==========================================================================
subroutine fa_time_unix_to_jd(unix_sec, jd, info) &
    bind(C, name="fa_time_unix_to_jd")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: unix_sec
    real(c_double), intent(out)       :: jd
    integer(c_int), intent(out)       :: info

    info = 0
    ! Unix epoch = JD 2440587.5 (1970-01-01 00:00:00 UTC)
    jd = unix_sec / 86400.0d0 + 2440587.5d0
end subroutine

! ==========================================================================
! fa_time_jd_to_unix — Julian Date to Unix epoch seconds
! ==========================================================================
subroutine fa_time_jd_to_unix(jd, unix_sec, info) &
    bind(C, name="fa_time_jd_to_unix")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: jd
    real(c_double), intent(out)       :: unix_sec
    integer(c_int), intent(out)       :: info

    info = 0
    unix_sec = (jd - 2440587.5d0) * 86400.0d0
end subroutine

! ==========================================================================
! fa_time_utc_to_tai — UTC JD to TAI JD
! ==========================================================================
subroutine fa_time_utc_to_tai(utc_jd, tai_jd, info) &
    bind(C, name="fa_time_utc_to_tai")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: utc_jd
    real(c_double), intent(out)       :: tai_jd
    integer(c_int), intent(out)       :: info

    real(c_double) :: dt_ls

    info = 0
    call get_leap_seconds(utc_jd, dt_ls)
    tai_jd = utc_jd + dt_ls / 86400.0d0

contains

    subroutine get_leap_seconds(jd, dt)
        implicit none
        real(c_double), intent(in)  :: jd
        real(c_double), intent(out) :: dt

        ! Leap second table: JD boundaries and cumulative leap seconds
        ! 28 entries from 1972-01-01 to 2017-01-01
        integer, parameter :: N_LS = 28
        real(c_double) :: ls_jd(N_LS), ls_val(N_LS)
        integer :: i

        ! JD values for leap second introduction dates
        ls_jd( 1) = 2441317.5d0;  ls_val( 1) = 10.0d0  ! 1972-01-01
        ls_jd( 2) = 2441499.5d0;  ls_val( 2) = 11.0d0  ! 1972-07-01
        ls_jd( 3) = 2441683.5d0;  ls_val( 3) = 12.0d0  ! 1973-01-01
        ls_jd( 4) = 2442048.5d0;  ls_val( 4) = 13.0d0  ! 1974-01-01
        ls_jd( 5) = 2442413.5d0;  ls_val( 5) = 14.0d0  ! 1975-01-01
        ls_jd( 6) = 2442778.5d0;  ls_val( 6) = 15.0d0  ! 1976-01-01
        ls_jd( 7) = 2443144.5d0;  ls_val( 7) = 16.0d0  ! 1977-01-01
        ls_jd( 8) = 2443509.5d0;  ls_val( 8) = 17.0d0  ! 1978-01-01
        ls_jd( 9) = 2443874.5d0;  ls_val( 9) = 18.0d0  ! 1979-01-01
        ls_jd(10) = 2444239.5d0;  ls_val(10) = 19.0d0  ! 1980-01-01
        ls_jd(11) = 2444786.5d0;  ls_val(11) = 20.0d0  ! 1981-07-01
        ls_jd(12) = 2445151.5d0;  ls_val(12) = 21.0d0  ! 1982-07-01
        ls_jd(13) = 2445516.5d0;  ls_val(13) = 22.0d0  ! 1983-07-01
        ls_jd(14) = 2446247.5d0;  ls_val(14) = 23.0d0  ! 1985-07-01
        ls_jd(15) = 2447161.5d0;  ls_val(15) = 24.0d0  ! 1988-01-01
        ls_jd(16) = 2447892.5d0;  ls_val(16) = 25.0d0  ! 1990-01-01
        ls_jd(17) = 2448257.5d0;  ls_val(17) = 26.0d0  ! 1991-01-01
        ls_jd(18) = 2448804.5d0;  ls_val(18) = 27.0d0  ! 1992-07-01
        ls_jd(19) = 2449169.5d0;  ls_val(19) = 28.0d0  ! 1993-07-01
        ls_jd(20) = 2449534.5d0;  ls_val(20) = 29.0d0  ! 1994-07-01
        ls_jd(21) = 2450083.5d0;  ls_val(21) = 30.0d0  ! 1996-01-01
        ls_jd(22) = 2450630.5d0;  ls_val(22) = 31.0d0  ! 1997-07-01
        ls_jd(23) = 2451179.5d0;  ls_val(23) = 32.0d0  ! 1999-01-01
        ls_jd(24) = 2453736.5d0;  ls_val(24) = 33.0d0  ! 2006-01-01
        ls_jd(25) = 2454832.5d0;  ls_val(25) = 34.0d0  ! 2009-01-01
        ls_jd(26) = 2456109.5d0;  ls_val(26) = 35.0d0  ! 2012-07-01
        ls_jd(27) = 2457204.5d0;  ls_val(27) = 36.0d0  ! 2015-07-01
        ls_jd(28) = 2457754.5d0;  ls_val(28) = 37.0d0  ! 2017-01-01

        ! Before 1972: no leap seconds in the modern sense
        if (jd < ls_jd(1)) then
            dt = 0.0d0
            return
        end if

        ! Find applicable leap second count (last entry <= jd)
        dt = ls_val(N_LS)
        do i = N_LS, 2, -1
            if (jd >= ls_jd(i)) then
                dt = ls_val(i)
                return
            end if
        end do
        dt = ls_val(1)
    end subroutine

end subroutine

! ==========================================================================
! fa_time_tai_to_utc — TAI JD to UTC JD
! ==========================================================================
subroutine fa_time_tai_to_utc(tai_jd, utc_jd, info) &
    bind(C, name="fa_time_tai_to_utc")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: tai_jd
    real(c_double), intent(out)       :: utc_jd
    integer(c_int), intent(out)       :: info

    real(c_double) :: dt_ls

    info = 0
    ! Approximate: subtract leap seconds based on TAI estimate of UTC
    ! This is iterative in principle but one step is sufficient for sub-ms
    call get_leap_seconds_approx(tai_jd, dt_ls)
    utc_jd = tai_jd - dt_ls / 86400.0d0

contains

    subroutine get_leap_seconds_approx(jd, dt)
        implicit none
        real(c_double), intent(in)  :: jd
        real(c_double), intent(out) :: dt

        integer, parameter :: N_LS = 28
        real(c_double) :: ls_jd(N_LS), ls_val(N_LS)
        integer :: i

        ls_jd( 1) = 2441317.5d0;  ls_val( 1) = 10.0d0
        ls_jd( 2) = 2441499.5d0;  ls_val( 2) = 11.0d0
        ls_jd( 3) = 2441683.5d0;  ls_val( 3) = 12.0d0
        ls_jd( 4) = 2442048.5d0;  ls_val( 4) = 13.0d0
        ls_jd( 5) = 2442413.5d0;  ls_val( 5) = 14.0d0
        ls_jd( 6) = 2442778.5d0;  ls_val( 6) = 15.0d0
        ls_jd( 7) = 2443144.5d0;  ls_val( 7) = 16.0d0
        ls_jd( 8) = 2443509.5d0;  ls_val( 8) = 17.0d0
        ls_jd( 9) = 2443874.5d0;  ls_val( 9) = 18.0d0
        ls_jd(10) = 2444239.5d0;  ls_val(10) = 19.0d0
        ls_jd(11) = 2444786.5d0;  ls_val(11) = 20.0d0
        ls_jd(12) = 2445151.5d0;  ls_val(12) = 21.0d0
        ls_jd(13) = 2445516.5d0;  ls_val(13) = 22.0d0
        ls_jd(14) = 2446247.5d0;  ls_val(14) = 23.0d0
        ls_jd(15) = 2447161.5d0;  ls_val(15) = 24.0d0
        ls_jd(16) = 2447892.5d0;  ls_val(16) = 25.0d0
        ls_jd(17) = 2448257.5d0;  ls_val(17) = 26.0d0
        ls_jd(18) = 2448804.5d0;  ls_val(18) = 27.0d0
        ls_jd(19) = 2449169.5d0;  ls_val(19) = 28.0d0
        ls_jd(20) = 2449534.5d0;  ls_val(20) = 29.0d0
        ls_jd(21) = 2450083.5d0;  ls_val(21) = 30.0d0
        ls_jd(22) = 2450630.5d0;  ls_val(22) = 31.0d0
        ls_jd(23) = 2451179.5d0;  ls_val(23) = 32.0d0
        ls_jd(24) = 2453736.5d0;  ls_val(24) = 33.0d0
        ls_jd(25) = 2454832.5d0;  ls_val(25) = 34.0d0
        ls_jd(26) = 2456109.5d0;  ls_val(26) = 35.0d0
        ls_jd(27) = 2457204.5d0;  ls_val(27) = 36.0d0
        ls_jd(28) = 2457754.5d0;  ls_val(28) = 37.0d0

        if (jd < ls_jd(1)) then
            dt = 0.0d0
            return
        end if

        dt = ls_val(N_LS)
        do i = N_LS, 2, -1
            if (jd >= ls_jd(i)) then
                dt = ls_val(i)
                return
            end if
        end do
        dt = ls_val(1)
    end subroutine

end subroutine

! ==========================================================================
! fa_time_tai_to_tt — TAI JD to Terrestrial Time JD
! ==========================================================================
subroutine fa_time_tai_to_tt(tai_jd, tt_jd, info) &
    bind(C, name="fa_time_tai_to_tt")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: tai_jd
    real(c_double), intent(out)       :: tt_jd
    integer(c_int), intent(out)       :: info

    info = 0
    ! TT = TAI + 32.184 seconds
    tt_jd = tai_jd + 32.184d0 / 86400.0d0
end subroutine

! ==========================================================================
! fa_time_tt_to_tai — TT JD to TAI JD
! ==========================================================================
subroutine fa_time_tt_to_tai(tt_jd, tai_jd, info) &
    bind(C, name="fa_time_tt_to_tai")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: tt_jd
    real(c_double), intent(out)       :: tai_jd
    integer(c_int), intent(out)       :: info

    info = 0
    tai_jd = tt_jd - 32.184d0 / 86400.0d0
end subroutine

! ==========================================================================
! fa_time_tt_to_tdb — TT JD to Barycentric Dynamical Time JD
! ==========================================================================
subroutine fa_time_tt_to_tdb(tt_jd, tdb_jd, info) &
    bind(C, name="fa_time_tt_to_tdb")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: tt_jd
    real(c_double), intent(out)       :: tdb_jd
    integer(c_int), intent(out)       :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: g

    info = 0
    ! Mean anomaly of Earth (degrees → radians)
    g = 2.0d0 * PI * (357.53d0 + 0.9856003d0 * (tt_jd - 2451545.0d0)) / 360.0d0
    ! Fairhead & Bretagnon approximation (dominant term)
    tdb_jd = tt_jd + 0.001658d0 * sin(g + 0.0167d0 * sin(g)) / 86400.0d0
end subroutine

! ==========================================================================
! fa_time_utc_to_gps — UTC JD to GPS time (seconds since GPS epoch)
! ==========================================================================
subroutine fa_time_utc_to_gps(utc_jd, gps_sec, info) &
    bind(C, name="fa_time_utc_to_gps")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: utc_jd
    real(c_double), intent(out)       :: gps_sec
    integer(c_int), intent(out)       :: info

    interface
        subroutine fa_time_utc_to_tai(u_jd, t_jd, inf) &
            bind(C, name="fa_time_utc_to_tai")
            use iso_c_binding
            real(c_double), value, intent(in) :: u_jd
            real(c_double), intent(out)       :: t_jd
            integer(c_int), intent(out)       :: inf
        end subroutine
    end interface

    real(c_double) :: tai_jd
    real(c_double), parameter :: GPS_EPOCH_JD = 2444244.5d0

    info = 0
    call fa_time_utc_to_tai(utc_jd, tai_jd, info)
    if (info /= 0) return
    gps_sec = (tai_jd - GPS_EPOCH_JD) * 86400.0d0 - 19.0d0
end subroutine

! ==========================================================================
! fa_time_gps_to_utc — GPS time (seconds) to UTC JD
! ==========================================================================
subroutine fa_time_gps_to_utc(gps_sec, utc_jd, info) &
    bind(C, name="fa_time_gps_to_utc")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: gps_sec
    real(c_double), intent(out)       :: utc_jd
    integer(c_int), intent(out)       :: info

    interface
        subroutine fa_time_tai_to_utc(t_jd, u_jd, inf) &
            bind(C, name="fa_time_tai_to_utc")
            use iso_c_binding
            real(c_double), value, intent(in) :: t_jd
            real(c_double), intent(out)       :: u_jd
            integer(c_int), intent(out)       :: inf
        end subroutine
    end interface

    real(c_double) :: tai_jd
    real(c_double), parameter :: GPS_EPOCH_JD = 2444244.5d0

    info = 0
    ! TAI = GPS + 19 seconds, then convert to JD
    tai_jd = GPS_EPOCH_JD + (gps_sec + 19.0d0) / 86400.0d0
    call fa_time_tai_to_utc(tai_jd, utc_jd, info)
end subroutine

! ==========================================================================
! fa_time_gmst — Greenwich Mean Sidereal Time (radians)
! ==========================================================================
subroutine fa_time_gmst(ut1_jd, gmst_rad, info) &
    bind(C, name="fa_time_gmst")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: ut1_jd
    real(c_double), intent(out)       :: gmst_rad
    integer(c_int), intent(out)       :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double), parameter :: TWOPI = 2.0d0 * PI
    real(c_double) :: T, gmst_sec

    info = 0
    ! Julian centuries from J2000.0
    T = (ut1_jd - 2451545.0d0) / 36525.0d0

    ! IAU 1982 GMST model (seconds of time)
    gmst_sec = 67310.54841d0 &
             + (876600.0d0 * 3600.0d0 + 8640184.812866d0) * T &
             + 0.093104d0 * T * T &
             - 6.2d-6 * T * T * T

    ! Convert seconds to radians (1 full rotation = 86400 sidereal seconds)
    gmst_rad = mod(gmst_sec * TWOPI / 86400.0d0, TWOPI)
    if (gmst_rad < 0.0d0) gmst_rad = gmst_rad + TWOPI
end subroutine

! ==========================================================================
! fa_time_cal_to_jd — Calendar date to Julian Date
! ==========================================================================
subroutine fa_time_cal_to_jd(year, month, day, hour, minute, second, jd, info) &
    bind(C, name="fa_time_cal_to_jd")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: year, month, day, hour, minute
    real(c_double), value, intent(in) :: second
    real(c_double), intent(out)       :: jd
    integer(c_int), intent(out)       :: info

    integer :: y, m, A, B

    info = 0

    ! Adjust for Jan/Feb (treat as months 13/14 of previous year)
    if (month <= 2) then
        y = year - 1
        m = month + 12
    else
        y = year
        m = month
    end if

    ! Gregorian calendar correction
    A = y / 100
    B = 2 - A + A / 4

    jd = int(365.25d0 * (y + 4716)) &
       + int(30.6001d0 * (m + 1)) &
       + day + B - 1524.5d0 &
       + (hour + minute / 60.0d0 + second / 3600.0d0) / 24.0d0
end subroutine

! ==========================================================================
! fa_time_jd_to_cal — Julian Date to calendar date
! ==========================================================================
subroutine fa_time_jd_to_cal(jd, year, month, day, hour, minute, second, info) &
    bind(C, name="fa_time_jd_to_cal")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: jd
    integer(c_int), intent(out)       :: year, month, day, hour, minute
    real(c_double), intent(out)       :: second
    integer(c_int), intent(out)       :: info

    integer :: Z, A, alpha, B, C, D, E
    real(c_double) :: F, day_frac

    info = 0

    ! Integer and fractional parts
    Z = int(jd + 0.5d0)
    F = (jd + 0.5d0) - Z

    if (Z < 2299161) then
        A = Z
    else
        alpha = int((Z - 1867216.25d0) / 36524.25d0)
        A = Z + 1 + alpha - alpha / 4
    end if

    B = A + 1524
    C = int((B - 122.1d0) / 365.25d0)
    D = int(365.25d0 * C)
    E = int((B - D) / 30.6001d0)

    day_frac = B - D - int(30.6001d0 * E) + F

    day = int(day_frac)
    day_frac = (day_frac - day) * 24.0d0
    hour = int(day_frac)
    day_frac = (day_frac - hour) * 60.0d0
    minute = int(day_frac)
    second = (day_frac - minute) * 60.0d0

    if (E < 14) then
        month = E - 1
    else
        month = E - 13
    end if

    if (month > 2) then
        year = C - 4716
    else
        year = C - 4715
    end if
end subroutine

! ==========================================================================
! fa_time_leap_seconds — Get leap second count for UTC JD
! ==========================================================================
subroutine fa_time_leap_seconds(utc_jd, dt_ls, info) &
    bind(C, name="fa_time_leap_seconds")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in) :: utc_jd
    real(c_double), intent(out)       :: dt_ls
    integer(c_int), intent(out)       :: info

    integer, parameter :: N_LS = 28
    real(c_double) :: ls_jd(N_LS), ls_val(N_LS)
    integer :: i

    info = 0

    ls_jd( 1) = 2441317.5d0;  ls_val( 1) = 10.0d0  ! 1972-01-01
    ls_jd( 2) = 2441499.5d0;  ls_val( 2) = 11.0d0  ! 1972-07-01
    ls_jd( 3) = 2441683.5d0;  ls_val( 3) = 12.0d0  ! 1973-01-01
    ls_jd( 4) = 2442048.5d0;  ls_val( 4) = 13.0d0  ! 1974-01-01
    ls_jd( 5) = 2442413.5d0;  ls_val( 5) = 14.0d0  ! 1975-01-01
    ls_jd( 6) = 2442778.5d0;  ls_val( 6) = 15.0d0  ! 1976-01-01
    ls_jd( 7) = 2443144.5d0;  ls_val( 7) = 16.0d0  ! 1977-01-01
    ls_jd( 8) = 2443509.5d0;  ls_val( 8) = 17.0d0  ! 1978-01-01
    ls_jd( 9) = 2443874.5d0;  ls_val( 9) = 18.0d0  ! 1979-01-01
    ls_jd(10) = 2444239.5d0;  ls_val(10) = 19.0d0  ! 1980-01-01
    ls_jd(11) = 2444786.5d0;  ls_val(11) = 20.0d0  ! 1981-07-01
    ls_jd(12) = 2445151.5d0;  ls_val(12) = 21.0d0  ! 1982-07-01
    ls_jd(13) = 2445516.5d0;  ls_val(13) = 22.0d0  ! 1983-07-01
    ls_jd(14) = 2446247.5d0;  ls_val(14) = 23.0d0  ! 1985-07-01
    ls_jd(15) = 2447161.5d0;  ls_val(15) = 24.0d0  ! 1988-01-01
    ls_jd(16) = 2447892.5d0;  ls_val(16) = 25.0d0  ! 1990-01-01
    ls_jd(17) = 2448257.5d0;  ls_val(17) = 26.0d0  ! 1991-01-01
    ls_jd(18) = 2448804.5d0;  ls_val(18) = 27.0d0  ! 1992-07-01
    ls_jd(19) = 2449169.5d0;  ls_val(19) = 28.0d0  ! 1993-07-01
    ls_jd(20) = 2449534.5d0;  ls_val(20) = 29.0d0  ! 1994-07-01
    ls_jd(21) = 2450083.5d0;  ls_val(21) = 30.0d0  ! 1996-01-01
    ls_jd(22) = 2450630.5d0;  ls_val(22) = 31.0d0  ! 1997-07-01
    ls_jd(23) = 2451179.5d0;  ls_val(23) = 32.0d0  ! 1999-01-01
    ls_jd(24) = 2453736.5d0;  ls_val(24) = 33.0d0  ! 2006-01-01
    ls_jd(25) = 2454832.5d0;  ls_val(25) = 34.0d0  ! 2009-01-01
    ls_jd(26) = 2456109.5d0;  ls_val(26) = 35.0d0  ! 2012-07-01
    ls_jd(27) = 2457204.5d0;  ls_val(27) = 36.0d0  ! 2015-07-01
    ls_jd(28) = 2457754.5d0;  ls_val(28) = 37.0d0  ! 2017-01-01

    if (utc_jd < ls_jd(1)) then
        dt_ls = 0.0d0
        return
    end if

    dt_ls = ls_val(N_LS)
    do i = N_LS, 2, -1
        if (utc_jd >= ls_jd(i)) then
            dt_ls = ls_val(i)
            return
        end if
    end do
    dt_ls = ls_val(1)
end subroutine
