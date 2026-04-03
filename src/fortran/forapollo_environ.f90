! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_environ.f90 — Environment models
!
! Atmosphere, gravity, solar radiation pressure, geodesics, magnetic field.
! All bind(C, name="fa_*"). Stateless, reentrant.
!
! Error codes: info = 0 ok, 3 invalid input, 4 convergence failure.

! ==========================================================================
! fa_environ_atmosphere_exp — Exponential atmosphere model
! ==========================================================================
subroutine fa_environ_atmosphere_exp(h, rho0, h_scale, rho, info) &
    bind(C, name="fa_environ_atmosphere_exp")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: h, rho0, h_scale
    real(c_double), intent(out)        :: rho
    integer(c_int), intent(out)        :: info

    info = 0
    if (h_scale <= 0.0d0) then
        info = 3
        rho = 0.0d0
        return
    end if
    rho = rho0 * exp(-h / h_scale)
end subroutine

! ==========================================================================
! fa_environ_atmosphere_us76 — US Standard Atmosphere 1976
! ==========================================================================
subroutine fa_environ_atmosphere_us76(h, rho, T_atm, p_atm, info) &
    bind(C, name="fa_environ_atmosphere_us76")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: h
    real(c_double), intent(out)        :: rho, T_atm, p_atm
    integer(c_int), intent(out)        :: info

    ! Constants
    real(c_double), parameter :: g0 = 9.80665d0       ! m/s^2
    real(c_double), parameter :: M_air = 0.0289644d0  ! kg/mol
    real(c_double), parameter :: R_gas = 8.31447d0    ! J/(mol·K)

    ! Layer base altitudes (km), temperatures (K), lapse rates (K/km), pressures (Pa)
    integer, parameter :: NL = 7
    real(c_double) :: h_base(NL), T_base(NL), lapse(NL), p_base(NL)
    real(c_double) :: h_m, dh
    integer :: layer, i

    info = 0

    h_base = (/ 0.0d0, 11.0d0, 20.0d0, 32.0d0, 47.0d0, 51.0d0, 71.0d0 /)
    T_base = (/ 288.15d0, 216.65d0, 216.65d0, 228.65d0, 270.65d0, 270.65d0, 214.65d0 /)
    lapse  = (/ -6.5d0, 0.0d0, 1.0d0, 2.8d0, 0.0d0, -2.8d0, -2.0d0 /)

    if (h < 0.0d0 .or. h > 86.0d0) then
        info = 3
        rho = 0.0d0
        T_atm = 0.0d0
        p_atm = 0.0d0
        return
    end if

    ! Compute base pressures layer by layer
    p_base(1) = 101325.0d0  ! sea level Pa
    do i = 2, NL
        dh = (h_base(i) - h_base(i-1)) * 1000.0d0  ! convert km to m
        if (abs(lapse(i-1)) < 1.0d-10) then
            ! Isothermal layer
            p_base(i) = p_base(i-1) * exp(-g0 * M_air * dh / (R_gas * T_base(i-1)))
        else
            ! Gradient layer
            p_base(i) = p_base(i-1) * (T_base(i) / T_base(i-1)) &
                       ** (-g0 * M_air / (R_gas * lapse(i-1) * 1.0d-3))
        end if
    end do

    ! Find layer
    layer = NL
    do i = NL, 2, -1
        if (h >= h_base(i)) then
            layer = i
            exit
        end if
        if (i == 2) layer = 1
    end do

    ! Temperature at altitude
    dh = h - h_base(layer)
    T_atm = T_base(layer) + lapse(layer) * dh

    ! Pressure at altitude
    h_m = dh * 1000.0d0  ! km to m
    if (abs(lapse(layer)) < 1.0d-10) then
        p_atm = p_base(layer) * exp(-g0 * M_air * h_m / (R_gas * T_base(layer)))
    else
        p_atm = p_base(layer) * (T_atm / T_base(layer)) &
               ** (-g0 * M_air / (R_gas * lapse(layer) * 1.0d-3))
    end if

    ! Density from ideal gas law
    rho = p_atm * M_air / (R_gas * T_atm)
end subroutine

! ==========================================================================
! fa_environ_gravity_pointmass — Point-mass gravity
! ==========================================================================
subroutine fa_environ_gravity_pointmass(r_vec, mu, g_vec, info) &
    bind(C, name="fa_environ_gravity_pointmass")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r_vec(3)
    real(c_double), value, intent(in)  :: mu
    real(c_double), intent(out)        :: g_vec(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: r_mag, r3

    info = 0
    r_mag = sqrt(r_vec(1)**2 + r_vec(2)**2 + r_vec(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        g_vec = 0.0d0
        return
    end if
    r3 = r_mag * r_mag * r_mag
    g_vec(1) = -mu * r_vec(1) / r3
    g_vec(2) = -mu * r_vec(2) / r3
    g_vec(3) = -mu * r_vec(3) / r3
end subroutine

! ==========================================================================
! fa_environ_gravity_j2 — J2 oblateness gravity
! ==========================================================================
subroutine fa_environ_gravity_j2(r_vec, mu, J2, R_eq, g_vec, info) &
    bind(C, name="fa_environ_gravity_j2")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r_vec(3)
    real(c_double), value, intent(in)  :: mu, J2, R_eq
    real(c_double), intent(out)        :: g_vec(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: x, y, z, r_mag, r2, r5, r7, Re2, fac, z2_r2

    info = 0
    x = r_vec(1); y = r_vec(2); z = r_vec(3)
    r_mag = sqrt(x*x + y*y + z*z)
    if (r_mag < 1.0d-12) then
        info = 3
        g_vec = 0.0d0
        return
    end if

    r2 = r_mag * r_mag
    r5 = r2 * r2 * r_mag
    Re2 = R_eq * R_eq
    z2_r2 = z * z / r2
    fac = -1.5d0 * J2 * mu * Re2 / r5

    ! Point-mass + J2 perturbation
    g_vec(1) = -mu * x / (r2 * r_mag) + fac * x * (1.0d0 - 5.0d0 * z2_r2)
    g_vec(2) = -mu * y / (r2 * r_mag) + fac * y * (1.0d0 - 5.0d0 * z2_r2)
    g_vec(3) = -mu * z / (r2 * r_mag) + fac * z * (3.0d0 - 5.0d0 * z2_r2)
end subroutine

! ==========================================================================
! fa_environ_gravity_j4 — J2 + J4 gravity
! ==========================================================================
subroutine fa_environ_gravity_j4(r_vec, mu, J2, J4, R_eq, g_vec, info) &
    bind(C, name="fa_environ_gravity_j4")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r_vec(3)
    real(c_double), value, intent(in)  :: mu, J2, J4, R_eq
    real(c_double), intent(out)        :: g_vec(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: x, y, z, r_mag, r2, r4, r7, r9, Re2, Re4
    real(c_double) :: z2, z4, z2_r2, z4_r4, fac2, fac4

    info = 0
    x = r_vec(1); y = r_vec(2); z = r_vec(3)
    r_mag = sqrt(x*x + y*y + z*z)
    if (r_mag < 1.0d-12) then
        info = 3
        g_vec = 0.0d0
        return
    end if

    r2 = r_mag * r_mag
    r4 = r2 * r2
    Re2 = R_eq * R_eq
    Re4 = Re2 * Re2
    z2 = z * z
    z4 = z2 * z2
    z2_r2 = z2 / r2
    z4_r4 = z4 / r4

    ! J2 contribution
    fac2 = -1.5d0 * J2 * mu * Re2 / (r4 * r_mag)

    ! J4 contribution
    fac4 = (15.0d0 / 8.0d0) * J4 * mu * Re4 / (r4 * r4 * r_mag)

    ! Point-mass
    g_vec(1) = -mu * x / (r2 * r_mag)
    g_vec(2) = -mu * y / (r2 * r_mag)
    g_vec(3) = -mu * z / (r2 * r_mag)

    ! J2
    g_vec(1) = g_vec(1) + fac2 * x * (1.0d0 - 5.0d0 * z2_r2)
    g_vec(2) = g_vec(2) + fac2 * y * (1.0d0 - 5.0d0 * z2_r2)
    g_vec(3) = g_vec(3) + fac2 * z * (3.0d0 - 5.0d0 * z2_r2)

    ! J4
    g_vec(1) = g_vec(1) + fac4 * x * (3.0d0 - 42.0d0 * z2_r2 + 63.0d0 * z4_r4)
    g_vec(2) = g_vec(2) + fac4 * y * (3.0d0 - 42.0d0 * z2_r2 + 63.0d0 * z4_r4)
    g_vec(3) = g_vec(3) + fac4 * z * (15.0d0 - 70.0d0 * z2_r2 + 63.0d0 * z4_r4)
end subroutine

! ==========================================================================
! fa_environ_srp — Solar radiation pressure acceleration
! ==========================================================================
subroutine fa_environ_srp(r_sat, r_sun, A_over_m, Cr, a_srp, info) &
    bind(C, name="fa_environ_srp")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r_sat(3), r_sun(3)
    real(c_double), value, intent(in)  :: A_over_m, Cr
    real(c_double), intent(out)        :: a_srp(3)
    integer(c_int), intent(out)        :: info

    ! Solar radiation pressure at 1 AU (N/m^2)
    real(c_double), parameter :: P_sr = 4.56d-6
    ! 1 AU in km
    real(c_double), parameter :: AU_km = 1.495978707d8

    real(c_double) :: d(3), d_mag, d_hat(3), r_au2

    info = 0

    ! Vector from satellite to sun
    d(1) = r_sun(1) - r_sat(1)
    d(2) = r_sun(2) - r_sat(2)
    d(3) = r_sun(3) - r_sat(3)
    d_mag = sqrt(d(1)**2 + d(2)**2 + d(3)**2)

    if (d_mag < 1.0d-12) then
        info = 3
        a_srp = 0.0d0
        return
    end if

    d_hat(1) = d(1) / d_mag
    d_hat(2) = d(2) / d_mag
    d_hat(3) = d(3) / d_mag

    ! Scale pressure by inverse square of distance (in AU)
    r_au2 = (d_mag / AU_km) ** 2

    ! Acceleration away from sun (opposite to sun direction)
    ! a = -P_sr * Cr * (A/m) * (1 AU / r)^2 * d_hat
    a_srp(1) = -P_sr * Cr * A_over_m / r_au2 * d_hat(1)
    a_srp(2) = -P_sr * Cr * A_over_m / r_au2 * d_hat(2)
    a_srp(3) = -P_sr * Cr * A_over_m / r_au2 * d_hat(3)
end subroutine

! ==========================================================================
! fa_environ_geodesic_vincenty — Vincenty inverse geodesic
! ==========================================================================
subroutine fa_environ_geodesic_vincenty(lat1, lon1, lat2, lon2, a_body, f_body, &
                                         dist, az1, az2, info) &
    bind(C, name="fa_environ_geodesic_vincenty")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: lat1, lon1, lat2, lon2, a_body, f_body
    real(c_double), intent(out)        :: dist, az1, az2
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: b, U1, U2, L, lambda, lambda_old
    real(c_double) :: sin_U1, cos_U1, sin_U2, cos_U2
    real(c_double) :: sin_lam, cos_lam, sin_sigma, cos_sigma, sigma
    real(c_double) :: sin_alpha, cos2_alpha, cos_2sigma_m, C_vin
    real(c_double) :: u2_vin, A_vin, B_vin, delta_sigma
    integer :: iter

    info = 0
    b = a_body * (1.0d0 - f_body)

    U1 = atan((1.0d0 - f_body) * tan(lat1))
    U2 = atan((1.0d0 - f_body) * tan(lat2))
    sin_U1 = sin(U1); cos_U1 = cos(U1)
    sin_U2 = sin(U2); cos_U2 = cos(U2)

    L = lon2 - lon1
    lambda = L

    do iter = 1, 200
        sin_lam = sin(lambda)
        cos_lam = cos(lambda)

        sin_sigma = sqrt((cos_U2 * sin_lam)**2 + &
                        (cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_lam)**2)

        if (sin_sigma < 1.0d-14) then
            ! Co-incident points
            dist = 0.0d0
            az1 = 0.0d0
            az2 = 0.0d0
            return
        end if

        cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lam
        sigma = atan2(sin_sigma, cos_sigma)

        sin_alpha = cos_U1 * cos_U2 * sin_lam / sin_sigma
        cos2_alpha = 1.0d0 - sin_alpha * sin_alpha

        if (abs(cos2_alpha) < 1.0d-14) then
            cos_2sigma_m = 0.0d0
        else
            cos_2sigma_m = cos_sigma - 2.0d0 * sin_U1 * sin_U2 / cos2_alpha
        end if

        C_vin = f_body / 16.0d0 * cos2_alpha * (4.0d0 + f_body * (4.0d0 - 3.0d0 * cos2_alpha))

        lambda_old = lambda
        lambda = L + (1.0d0 - C_vin) * f_body * sin_alpha * &
                (sigma + C_vin * sin_sigma * (cos_2sigma_m + C_vin * cos_sigma * &
                (-1.0d0 + 2.0d0 * cos_2sigma_m * cos_2sigma_m)))

        if (abs(lambda - lambda_old) < 1.0d-12) exit

        if (iter == 200) then
            info = 4
            dist = 0.0d0
            az1 = 0.0d0
            az2 = 0.0d0
            return
        end if
    end do

    u2_vin = cos2_alpha * (a_body * a_body - b * b) / (b * b)
    A_vin = 1.0d0 + u2_vin / 16384.0d0 * (4096.0d0 + u2_vin * &
           (-768.0d0 + u2_vin * (320.0d0 - 175.0d0 * u2_vin)))
    B_vin = u2_vin / 1024.0d0 * (256.0d0 + u2_vin * &
           (-128.0d0 + u2_vin * (74.0d0 - 47.0d0 * u2_vin)))

    delta_sigma = B_vin * sin_sigma * (cos_2sigma_m + B_vin / 4.0d0 * &
                 (cos_sigma * (-1.0d0 + 2.0d0 * cos_2sigma_m**2) - &
                  B_vin / 6.0d0 * cos_2sigma_m * (-3.0d0 + 4.0d0 * sin_sigma**2) * &
                  (-3.0d0 + 4.0d0 * cos_2sigma_m**2)))

    dist = b * A_vin * (sigma - delta_sigma)

    az1 = atan2(cos_U2 * sin(lambda), &
               cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos(lambda))
    az2 = atan2(cos_U1 * sin(lambda), &
               -sin_U1 * cos_U2 + cos_U1 * sin_U2 * cos(lambda))
end subroutine

! ==========================================================================
! fa_environ_geodesic_haversine — Great-circle distance on sphere
! ==========================================================================
subroutine fa_environ_geodesic_haversine(lat1, lon1, lat2, lon2, R_body, dist, info) &
    bind(C, name="fa_environ_geodesic_haversine")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: lat1, lon1, lat2, lon2, R_body
    real(c_double), intent(out)        :: dist
    integer(c_int), intent(out)        :: info

    real(c_double) :: dlat, dlon, a_hav, c_hav

    info = 0
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a_hav = sin(dlat / 2.0d0)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2.0d0)**2
    c_hav = 2.0d0 * atan2(sqrt(a_hav), sqrt(1.0d0 - a_hav))
    dist = R_body * c_hav
end subroutine

! ==========================================================================
! fa_environ_magnetic_dipole — Dipole magnetic field model
! ==========================================================================
subroutine fa_environ_magnetic_dipole(r_vec, m_dipole, B_vec, info) &
    bind(C, name="fa_environ_magnetic_dipole")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r_vec(3), m_dipole(3)
    real(c_double), intent(out)        :: B_vec(3)
    integer(c_int), intent(out)        :: info

    ! mu_0 / (4 * pi) in SI: 1e-7 T·m/A
    real(c_double), parameter :: mu0_4pi = 1.0d-7
    real(c_double) :: r_mag, r2, r3, r5, m_dot_r, r_hat(3)

    info = 0
    r_mag = sqrt(r_vec(1)**2 + r_vec(2)**2 + r_vec(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        B_vec = 0.0d0
        return
    end if

    r2 = r_mag * r_mag
    r3 = r2 * r_mag
    r5 = r2 * r3

    m_dot_r = m_dipole(1) * r_vec(1) + m_dipole(2) * r_vec(2) + m_dipole(3) * r_vec(3)

    ! B = (mu0/4pi) * (3*(m·r)*r/r^5 - m/r^3)
    B_vec(1) = mu0_4pi * (3.0d0 * m_dot_r * r_vec(1) / r5 - m_dipole(1) / r3)
    B_vec(2) = mu0_4pi * (3.0d0 * m_dot_r * r_vec(2) / r5 - m_dipole(2) / r3)
    B_vec(3) = mu0_4pi * (3.0d0 * m_dot_r * r_vec(3) / r5 - m_dipole(3) / r3)
end subroutine
