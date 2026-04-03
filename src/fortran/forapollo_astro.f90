! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_astro.f90 — Astrodynamics utilities
!
! Kepler solver, orbital mechanics, transfers, eclipse geometry,
! planetary constants, Stumpff functions, universal Kepler propagation.
!
! All bind(C, name="fa_*"). Stateless, reentrant.
! Error codes: info = 0 ok, 3 invalid input, 4 convergence failure.

! ==========================================================================
! fa_astro_kepler_solve — Solve Kepler's equation M = E - e*sin(E)
! ==========================================================================
subroutine fa_astro_kepler_solve(mean_anom, ecc, ecc_anom, info) &
    bind(C, name="fa_astro_kepler_solve")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: mean_anom, ecc
    real(c_double), intent(out)        :: ecc_anom
    integer(c_int), intent(out)        :: info

    real(c_double) :: d_ea, f_val, fp_val
    integer :: iter

    info = 0
    if (ecc < 0.0d0 .or. ecc >= 1.0d0) then
        info = 3
        ecc_anom = 0.0d0
        return
    end if

    ! Initial guess
    ecc_anom = mean_anom + 0.85d0 * sign(1.0d0, sin(mean_anom)) * ecc

    do iter = 1, 50
        f_val = ecc_anom - ecc * sin(ecc_anom) - mean_anom
        fp_val = 1.0d0 - ecc * cos(ecc_anom)
        if (abs(fp_val) < 1.0d-15) then
            info = 4
            return
        end if
        d_ea = f_val / fp_val
        ecc_anom = ecc_anom - d_ea
        if (abs(d_ea) < 1.0d-14) return
    end do

    info = 4  ! did not converge
end subroutine

! ==========================================================================
! fa_astro_true_anomaly — Eccentric anomaly to true anomaly
! ==========================================================================
subroutine fa_astro_true_anomaly(ecc_anom, ecc, nu, info) &
    bind(C, name="fa_astro_true_anomaly")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: ecc_anom, ecc
    real(c_double), intent(out)        :: nu
    integer(c_int), intent(out)        :: info

    info = 0
    if (ecc < 0.0d0 .or. ecc >= 1.0d0) then
        info = 3
        nu = 0.0d0
        return
    end if
    nu = 2.0d0 * atan2(sqrt(1.0d0 + ecc) * sin(ecc_anom / 2.0d0), &
                        sqrt(1.0d0 - ecc) * cos(ecc_anom / 2.0d0))
end subroutine

! ==========================================================================
! fa_astro_mean_anomaly — Eccentric anomaly to mean anomaly
! ==========================================================================
subroutine fa_astro_mean_anomaly(ecc_anom, ecc, mean_anom, info) &
    bind(C, name="fa_astro_mean_anomaly")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: ecc_anom, ecc
    real(c_double), intent(out)        :: mean_anom
    integer(c_int), intent(out)        :: info

    info = 0
    mean_anom = ecc_anom - ecc * sin(ecc_anom)
end subroutine

! ==========================================================================
! fa_astro_vis_viva — Vis-viva equation: v = sqrt(mu*(2/r - 1/a))
! ==========================================================================
subroutine fa_astro_vis_viva(r, a, mu, v, info) &
    bind(C, name="fa_astro_vis_viva")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: r, a, mu
    real(c_double), intent(out)        :: v
    integer(c_int), intent(out)        :: info

    real(c_double) :: val

    info = 0
    if (r <= 0.0d0 .or. abs(a) < 1.0d-12 .or. mu <= 0.0d0) then
        info = 3
        v = 0.0d0
        return
    end if
    val = mu * (2.0d0 / r - 1.0d0 / a)
    if (val < 0.0d0) then
        info = 3
        v = 0.0d0
        return
    end if
    v = sqrt(val)
end subroutine

! ==========================================================================
! fa_astro_period — Orbital period: T = 2*pi*sqrt(a^3/mu)
! ==========================================================================
subroutine fa_astro_period(a, mu, T, info) &
    bind(C, name="fa_astro_period")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: a, mu
    real(c_double), intent(out)        :: T
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0

    info = 0
    if (a <= 0.0d0 .or. mu <= 0.0d0) then
        info = 3
        T = 0.0d0
        return
    end if
    T = 2.0d0 * PI * sqrt(a**3 / mu)
end subroutine

! ==========================================================================
! fa_astro_soi — Sphere of influence: r_soi = a * (m_body/m_central)^(2/5)
! ==========================================================================
subroutine fa_astro_soi(a_orbit, m_body, m_central, r_soi, info) &
    bind(C, name="fa_astro_soi")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: a_orbit, m_body, m_central
    real(c_double), intent(out)        :: r_soi
    integer(c_int), intent(out)        :: info

    info = 0
    if (a_orbit <= 0.0d0 .or. m_body <= 0.0d0 .or. m_central <= 0.0d0) then
        info = 3
        r_soi = 0.0d0
        return
    end if
    r_soi = a_orbit * (m_body / m_central) ** 0.4d0
end subroutine

! ==========================================================================
! fa_astro_hohmann — Hohmann transfer delta-v and time of flight
! ==========================================================================
subroutine fa_astro_hohmann(r1, r2, mu, dv1, dv2, tof, info) &
    bind(C, name="fa_astro_hohmann")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: r1, r2, mu
    real(c_double), intent(out)        :: dv1, dv2, tof
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: a_t, v_c1, v_c2, v_t1, v_t2

    info = 0
    if (r1 <= 0.0d0 .or. r2 <= 0.0d0 .or. mu <= 0.0d0) then
        info = 3
        dv1 = 0.0d0; dv2 = 0.0d0; tof = 0.0d0
        return
    end if

    ! Transfer orbit semi-major axis
    a_t = (r1 + r2) / 2.0d0

    ! Circular velocities
    v_c1 = sqrt(mu / r1)
    v_c2 = sqrt(mu / r2)

    ! Transfer orbit velocities at periapsis and apoapsis
    v_t1 = sqrt(mu * (2.0d0 / r1 - 1.0d0 / a_t))
    v_t2 = sqrt(mu * (2.0d0 / r2 - 1.0d0 / a_t))

    ! Delta-v (magnitudes, signed for direction)
    dv1 = v_t1 - v_c1
    dv2 = v_c2 - v_t2

    ! Time of flight = half the transfer orbit period
    tof = PI * sqrt(a_t**3 / mu)
end subroutine

! ==========================================================================
! fa_astro_bielliptic — Bi-elliptic transfer
! ==========================================================================
subroutine fa_astro_bielliptic(r1, r2, r_int, mu, dv1, dv2, dv3, tof, info) &
    bind(C, name="fa_astro_bielliptic")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: r1, r2, r_int, mu
    real(c_double), intent(out)        :: dv1, dv2, dv3, tof
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: a1, a2, v_c1, v_c2, v1_dep, v1_arr, v2_dep, v2_arr

    info = 0
    if (r1 <= 0.0d0 .or. r2 <= 0.0d0 .or. r_int <= 0.0d0 .or. mu <= 0.0d0) then
        info = 3
        dv1 = 0.0d0; dv2 = 0.0d0; dv3 = 0.0d0; tof = 0.0d0
        return
    end if

    ! First transfer ellipse: r1 → r_int
    a1 = (r1 + r_int) / 2.0d0
    ! Second transfer ellipse: r_int → r2
    a2 = (r_int + r2) / 2.0d0

    v_c1 = sqrt(mu / r1)
    v_c2 = sqrt(mu / r2)

    ! First burn: depart r1
    v1_dep = sqrt(mu * (2.0d0 / r1 - 1.0d0 / a1))
    dv1 = v1_dep - v_c1

    ! Second burn: at r_int, transfer between ellipses
    v1_arr = sqrt(mu * (2.0d0 / r_int - 1.0d0 / a1))
    v2_dep = sqrt(mu * (2.0d0 / r_int - 1.0d0 / a2))
    dv2 = v2_dep - v1_arr

    ! Third burn: circularize at r2
    v2_arr = sqrt(mu * (2.0d0 / r2 - 1.0d0 / a2))
    dv3 = v_c2 - v2_arr

    ! Total time of flight
    tof = PI * sqrt(a1**3 / mu) + PI * sqrt(a2**3 / mu)
end subroutine

! ==========================================================================
! fa_astro_eclipse_conical — Conical shadow model
! ==========================================================================
subroutine fa_astro_eclipse_conical(pos_sat, pos_sun, pos_body, rad_body, shadow, info) &
    bind(C, name="fa_astro_eclipse_conical")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: pos_sat(3), pos_sun(3), pos_body(3)
    real(c_double), value, intent(in)  :: rad_body
    real(c_double), intent(out)        :: shadow
    integer(c_int), intent(out)        :: info

    ! Sun radius (km)
    real(c_double), parameter :: RAD_SUN = 696000.0d0

    real(c_double) :: s2b(3), s2sun(3), d_body, d_sun
    real(c_double) :: theta_body, theta_sun, theta_sep

    info = 0

    ! Vector from satellite to body center
    s2b(1) = pos_body(1) - pos_sat(1)
    s2b(2) = pos_body(2) - pos_sat(2)
    s2b(3) = pos_body(3) - pos_sat(3)

    ! Vector from satellite to sun
    s2sun(1) = pos_sun(1) - pos_sat(1)
    s2sun(2) = pos_sun(2) - pos_sat(2)
    s2sun(3) = pos_sun(3) - pos_sat(3)

    d_body = sqrt(s2b(1)**2 + s2b(2)**2 + s2b(3)**2)
    d_sun = sqrt(s2sun(1)**2 + s2sun(2)**2 + s2sun(3)**2)

    if (d_body < 1.0d-12 .or. d_sun < 1.0d-12) then
        info = 3
        shadow = 0.0d0
        return
    end if

    ! Angular radii as seen from satellite
    theta_body = asin(min(rad_body / d_body, 1.0d0))
    theta_sun = asin(min(RAD_SUN / d_sun, 1.0d0))

    ! Angular separation between body center and sun
    theta_sep = acos(max(-1.0d0, min(1.0d0, &
                (s2b(1)*s2sun(1) + s2b(2)*s2sun(2) + s2b(3)*s2sun(3)) / (d_body * d_sun))))

    ! Shadow determination
    if (theta_sep >= theta_body + theta_sun) then
        ! Full sunlight
        shadow = 0.0d0
    else if (theta_sep <= theta_body - theta_sun) then
        ! Full shadow (body completely occults sun)
        shadow = 1.0d0
    else
        ! Penumbra — linear interpolation
        shadow = 1.0d0 - (theta_sep - (theta_body - theta_sun)) / (2.0d0 * theta_sun)
        shadow = max(0.0d0, min(1.0d0, shadow))
    end if
end subroutine

! ==========================================================================
! fa_astro_ground_track — Sub-satellite point (lat, lon) from ECI position
! ==========================================================================
subroutine fa_astro_ground_track(r_eci, gmst, lat, lon, info) &
    bind(C, name="fa_astro_ground_track")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r_eci(3)
    real(c_double), value, intent(in)  :: gmst
    real(c_double), intent(out)        :: lat, lon
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: x_ecef, y_ecef, r_xy

    info = 0

    ! Rotate ECI to ECEF via GMST
    x_ecef = r_eci(1) * cos(gmst) + r_eci(2) * sin(gmst)
    y_ecef = -r_eci(1) * sin(gmst) + r_eci(2) * cos(gmst)

    r_xy = sqrt(x_ecef**2 + y_ecef**2)

    ! Geocentric latitude
    lat = atan2(r_eci(3), r_xy)

    ! Longitude
    lon = atan2(y_ecef, x_ecef)
end subroutine

! ==========================================================================
! fa_astro_planetary_mu — Gravitational parameter (km^3/s^2)
! ==========================================================================
subroutine fa_astro_planetary_mu(body_id, mu, info) &
    bind(C, name="fa_astro_planetary_mu")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: body_id
    real(c_double), intent(out)        :: mu
    integer(c_int), intent(out)        :: info

    info = 0
    select case (body_id)
    case (1);  mu = 22032.09d0          ! Mercury
    case (2);  mu = 324858.63d0         ! Venus
    case (3);  mu = 398600.4418d0       ! Earth
    case (4);  mu = 42828.375d0         ! Mars
    case (5);  mu = 126686534.9d0       ! Jupiter
    case (6);  mu = 37931187.0d0        ! Saturn
    case (7);  mu = 5793939.0d0         ! Uranus
    case (8);  mu = 6836529.0d0         ! Neptune
    case (9);  mu = 871.0d0             ! Pluto
    case (10); mu = 4902.800066d0       ! Moon
    case (11); mu = 132712440041.94d0   ! Sun
    case default
        info = 3
        mu = 0.0d0
    end select
end subroutine

! ==========================================================================
! fa_astro_planetary_radius — Mean equatorial radius (km)
! ==========================================================================
subroutine fa_astro_planetary_radius(body_id, radius, info) &
    bind(C, name="fa_astro_planetary_radius")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: body_id
    real(c_double), intent(out)        :: radius
    integer(c_int), intent(out)        :: info

    info = 0
    select case (body_id)
    case (1);  radius = 2439.7d0        ! Mercury
    case (2);  radius = 6051.8d0        ! Venus
    case (3);  radius = 6378.137d0      ! Earth
    case (4);  radius = 3396.2d0        ! Mars
    case (5);  radius = 71492.0d0       ! Jupiter
    case (6);  radius = 60268.0d0       ! Saturn
    case (7);  radius = 25559.0d0       ! Uranus
    case (8);  radius = 24764.0d0       ! Neptune
    case (9);  radius = 1188.3d0        ! Pluto
    case (10); radius = 1737.4d0        ! Moon
    case (11); radius = 696000.0d0      ! Sun
    case default
        info = 3
        radius = 0.0d0
    end select
end subroutine

! ==========================================================================
! fa_astro_stumpff_c2 — Stumpff function c2(psi)
! ==========================================================================
subroutine fa_astro_stumpff_c2(psi, c2, info) &
    bind(C, name="fa_astro_stumpff_c2")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: psi
    real(c_double), intent(out)        :: c2
    integer(c_int), intent(out)        :: info

    real(c_double) :: sp

    info = 0
    if (psi > 1.0d-6) then
        sp = sqrt(psi)
        c2 = (1.0d0 - cos(sp)) / psi
    else if (psi < -1.0d-6) then
        sp = sqrt(-psi)
        c2 = (cosh(sp) - 1.0d0) / (-psi)
    else
        c2 = 0.5d0
    end if
end subroutine

! ==========================================================================
! fa_astro_stumpff_c3 — Stumpff function c3(psi)
! ==========================================================================
subroutine fa_astro_stumpff_c3(psi, c3, info) &
    bind(C, name="fa_astro_stumpff_c3")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: psi
    real(c_double), intent(out)        :: c3
    integer(c_int), intent(out)        :: info

    real(c_double) :: sp

    info = 0
    if (psi > 1.0d-6) then
        sp = sqrt(psi)
        c3 = (sp - sin(sp)) / (psi * sp)
    else if (psi < -1.0d-6) then
        sp = sqrt(-psi)
        c3 = (sinh(sp) - sp) / ((-psi) * sp)
    else
        c3 = 1.0d0 / 6.0d0
    end if
end subroutine

! ==========================================================================
! fa_astro_universal_kepler — Universal variable Kepler propagation
! ==========================================================================
subroutine fa_astro_universal_kepler(r0, vr0, alpha, dt, mu, chi, info) &
    bind(C, name="fa_astro_universal_kepler")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: r0, vr0, alpha, dt, mu
    real(c_double), intent(out)        :: chi
    integer(c_int), intent(out)        :: info

    real(c_double) :: chi_n, psi, c2, c3, r, f_val, fp_val, dchi
    integer :: iter, sub_info

    info = 0
    if (r0 <= 0.0d0 .or. mu <= 0.0d0) then
        info = 3
        chi = 0.0d0
        return
    end if

    ! Initial guess
    chi_n = sqrt(mu) * abs(alpha) * dt

    do iter = 1, 50
        psi = chi_n * chi_n * alpha

        call stumpff_c2_local(psi, c2)
        call stumpff_c3_local(psi, c3)

        r = chi_n**2 * c2 + vr0 / sqrt(mu) * chi_n * (1.0d0 - psi * c3) + r0 * (1.0d0 - psi * c2)

        if (abs(r) < 1.0d-15) then
            info = 4
            chi = chi_n
            return
        end if

        f_val = r0 * vr0 / sqrt(mu) * chi_n**2 * c2 &
              + (1.0d0 - r0 * alpha) * chi_n**3 * c3 &
              + r0 * chi_n - sqrt(mu) * dt

        fp_val = r

        dchi = f_val / fp_val
        chi_n = chi_n - dchi

        if (abs(dchi) < 1.0d-12) then
            chi = chi_n
            return
        end if
    end do

    info = 4
    chi = chi_n

contains

    subroutine stumpff_c2_local(ps, cc2)
        implicit none
        real(c_double), intent(in)  :: ps
        real(c_double), intent(out) :: cc2
        real(c_double) :: sps
        if (ps > 1.0d-6) then
            sps = sqrt(ps)
            cc2 = (1.0d0 - cos(sps)) / ps
        else if (ps < -1.0d-6) then
            sps = sqrt(-ps)
            cc2 = (cosh(sps) - 1.0d0) / (-ps)
        else
            cc2 = 0.5d0
        end if
    end subroutine

    subroutine stumpff_c3_local(ps, cc3)
        implicit none
        real(c_double), intent(in)  :: ps
        real(c_double), intent(out) :: cc3
        real(c_double) :: sps
        if (ps > 1.0d-6) then
            sps = sqrt(ps)
            cc3 = (sps - sin(sps)) / (ps * sps)
        else if (ps < -1.0d-6) then
            sps = sqrt(-ps)
            cc3 = (sinh(sps) - sps) / ((-ps) * sps)
        else
            cc3 = 1.0d0 / 6.0d0
        end if
    end subroutine

end subroutine
