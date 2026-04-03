! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_coords.f90 — Coordinate frame transformations
!
! Inertial, rotating, local, orbital, geodetic frame conversions.
! Orbital element conversions (classical, equinoctial, Cartesian).
! All bind(C, name="fa_*"). Stateless, reentrant.
!
! Frame IDs:
!   1 = ECI_J2000    10 = ECEF        20 = NED    30 = PQW
!   2 = ICRF         11 = Moon-fixed   21 = ENU    31 = RSW
!                    12 = Body-fixed   22 = LVLH   32 = VNC
!   40 = WGS84 geodetic   41 = Generic ellipsoid
!   50 = Topocentric      51 = Pinhole camera
!
! Error codes: info = 0 ok, 3 invalid input, 4 convergence failure.
! All rotation matrices R(9) are flat 3x3 row-major.

! ==========================================================================
! fa_coords_eci_to_ecef — ECI to ECEF via GMST rotation
! ==========================================================================
subroutine fa_coords_eci_to_ecef(x_eci, x_ecef, gmst, info) &
    bind(C, name="fa_coords_eci_to_ecef")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_eci(3)
    real(c_double), intent(out)        :: x_ecef(3)
    real(c_double), value, intent(in)  :: gmst
    integer(c_int), intent(out)        :: info

    real(c_double) :: cg, sg

    info = 0
    cg = cos(gmst)
    sg = sin(gmst)

    ! R_z(gmst) rotation
    x_ecef(1) =  cg * x_eci(1) + sg * x_eci(2)
    x_ecef(2) = -sg * x_eci(1) + cg * x_eci(2)
    x_ecef(3) = x_eci(3)
end subroutine

! ==========================================================================
! fa_coords_ecef_to_eci — ECEF to ECI via GMST rotation
! ==========================================================================
subroutine fa_coords_ecef_to_eci(x_ecef, x_eci, gmst, info) &
    bind(C, name="fa_coords_ecef_to_eci")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_ecef(3)
    real(c_double), intent(out)        :: x_eci(3)
    real(c_double), value, intent(in)  :: gmst
    integer(c_int), intent(out)        :: info

    real(c_double) :: cg, sg

    info = 0
    cg = cos(gmst)
    sg = sin(gmst)

    ! R_z(-gmst) = R_z(gmst)^T
    x_eci(1) = cg * x_ecef(1) - sg * x_ecef(2)
    x_eci(2) = sg * x_ecef(1) + cg * x_ecef(2)
    x_eci(3) = x_ecef(3)
end subroutine

! ==========================================================================
! fa_coords_ecef_to_geodetic — ECEF to geodetic (Bowring's method)
! ==========================================================================
subroutine fa_coords_ecef_to_geodetic(x_ecef, lat, lon, alt, a_body, f_body, info) &
    bind(C, name="fa_coords_ecef_to_geodetic")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_ecef(3)
    real(c_double), intent(out)        :: lat, lon, alt
    real(c_double), value, intent(in)  :: a_body, f_body
    integer(c_int), intent(out)        :: info

    real(c_double) :: x, y, z, e2, b, ep2, p, theta, N_ell
    integer :: iter

    info = 0
    x = x_ecef(1); y = x_ecef(2); z = x_ecef(3)

    e2 = 2.0d0 * f_body - f_body * f_body
    b = a_body * (1.0d0 - f_body)
    ep2 = (a_body * a_body - b * b) / (b * b)
    p = sqrt(x * x + y * y)

    ! Longitude
    lon = atan2(y, x)

    ! Bowring's iterative method for latitude
    theta = atan2(z * a_body, p * b)

    do iter = 1, 5
        lat = atan2(z + ep2 * b * sin(theta)**3, &
                    p - e2 * a_body * cos(theta)**3)
        theta = atan2(z, p * (1.0d0 - e2 * a_body / &
                sqrt(a_body**2 * cos(lat)**2 + b**2 * sin(lat)**2)))
    end do

    ! Altitude
    N_ell = a_body / sqrt(1.0d0 - e2 * sin(lat)**2)
    if (abs(cos(lat)) > 1.0d-10) then
        alt = p / cos(lat) - N_ell
    else
        alt = abs(z) / sin(lat) - N_ell * (1.0d0 - e2)
    end if
end subroutine

! ==========================================================================
! fa_coords_geodetic_to_ecef — Geodetic to ECEF
! ==========================================================================
subroutine fa_coords_geodetic_to_ecef(lat, lon, alt, x_ecef, a_body, f_body, info) &
    bind(C, name="fa_coords_geodetic_to_ecef")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: lat, lon, alt, a_body, f_body
    real(c_double), intent(out)        :: x_ecef(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: e2, N_ell, clat, slat, clon, slon

    info = 0
    e2 = 2.0d0 * f_body - f_body * f_body
    clat = cos(lat); slat = sin(lat)
    clon = cos(lon); slon = sin(lon)

    N_ell = a_body / sqrt(1.0d0 - e2 * slat * slat)

    x_ecef(1) = (N_ell + alt) * clat * clon
    x_ecef(2) = (N_ell + alt) * clat * slon
    x_ecef(3) = (N_ell * (1.0d0 - e2) + alt) * slat
end subroutine

! ==========================================================================
! fa_coords_ecef_to_ned — ECEF to NED relative to reference point
! ==========================================================================
subroutine fa_coords_ecef_to_ned(x_ecef, x_ned, lat_ref, lon_ref, alt_ref, &
                                  a_body, f_body, info) &
    bind(C, name="fa_coords_ecef_to_ned")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_ecef(3)
    real(c_double), intent(out)        :: x_ned(3)
    real(c_double), value, intent(in)  :: lat_ref, lon_ref, alt_ref, a_body, f_body
    integer(c_int), intent(out)        :: info

    real(c_double) :: ref_ecef(3), dx(3)
    real(c_double) :: clat, slat, clon, slon, e2_loc, N_loc

    info = 0

    ! Inline geodetic→ECEF for reference point
    e2_loc = 2.0d0 * f_body - f_body * f_body
    clat = cos(lat_ref); slat = sin(lat_ref)
    clon = cos(lon_ref); slon = sin(lon_ref)
    N_loc = a_body / sqrt(1.0d0 - e2_loc * slat * slat)
    ref_ecef(1) = (N_loc + alt_ref) * clat * clon
    ref_ecef(2) = (N_loc + alt_ref) * clat * slon
    ref_ecef(3) = (N_loc * (1.0d0 - e2_loc) + alt_ref) * slat

    dx(1) = x_ecef(1) - ref_ecef(1)
    dx(2) = x_ecef(2) - ref_ecef(2)
    dx(3) = x_ecef(3) - ref_ecef(3)

    clat = cos(lat_ref); slat = sin(lat_ref)
    clon = cos(lon_ref); slon = sin(lon_ref)

    ! ECEF→NED rotation
    x_ned(1) = -slat * clon * dx(1) - slat * slon * dx(2) + clat * dx(3)
    x_ned(2) = -slon * dx(1) + clon * dx(2)
    x_ned(3) = -clat * clon * dx(1) - clat * slon * dx(2) - slat * dx(3)
end subroutine

! ==========================================================================
! fa_coords_ned_to_ecef — NED to ECEF
! ==========================================================================
subroutine fa_coords_ned_to_ecef(x_ned, x_ecef, lat_ref, lon_ref, alt_ref, &
                                  a_body, f_body, info) &
    bind(C, name="fa_coords_ned_to_ecef")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_ned(3)
    real(c_double), intent(out)        :: x_ecef(3)
    real(c_double), value, intent(in)  :: lat_ref, lon_ref, alt_ref, a_body, f_body
    integer(c_int), intent(out)        :: info

    real(c_double) :: ref_ecef(3), dx(3)
    real(c_double) :: clat, slat, clon, slon, e2_loc, N_loc

    info = 0

    ! Inline geodetic→ECEF for reference point
    e2_loc = 2.0d0 * f_body - f_body * f_body
    clat = cos(lat_ref); slat = sin(lat_ref)
    clon = cos(lon_ref); slon = sin(lon_ref)
    N_loc = a_body / sqrt(1.0d0 - e2_loc * slat * slat)
    ref_ecef(1) = (N_loc + alt_ref) * clat * clon
    ref_ecef(2) = (N_loc + alt_ref) * clat * slon
    ref_ecef(3) = (N_loc * (1.0d0 - e2_loc) + alt_ref) * slat

    ! NED→ECEF rotation (transpose of ECEF→NED)
    dx(1) = -slat * clon * x_ned(1) - slon * x_ned(2) - clat * clon * x_ned(3)
    dx(2) = -slat * slon * x_ned(1) + clon * x_ned(2) - clat * slon * x_ned(3)
    dx(3) =  clat * x_ned(1)                           - slat * x_ned(3)

    x_ecef(1) = ref_ecef(1) + dx(1)
    x_ecef(2) = ref_ecef(2) + dx(2)
    x_ecef(3) = ref_ecef(3) + dx(3)

end subroutine

! ==========================================================================
! fa_coords_ecef_to_enu — ECEF to ENU
! ==========================================================================
subroutine fa_coords_ecef_to_enu(x_ecef, x_enu, lat_ref, lon_ref, alt_ref, &
                                  a_body, f_body, info) &
    bind(C, name="fa_coords_ecef_to_enu")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_ecef(3)
    real(c_double), intent(out)        :: x_enu(3)
    real(c_double), value, intent(in)  :: lat_ref, lon_ref, alt_ref, a_body, f_body
    integer(c_int), intent(out)        :: info

    real(c_double) :: ref_ecef(3), dx(3)
    real(c_double) :: clat, slat, clon, slon, e2_loc, N_loc

    info = 0

    ! Inline geodetic→ECEF for reference point
    e2_loc = 2.0d0 * f_body - f_body * f_body
    clat = cos(lat_ref); slat = sin(lat_ref)
    clon = cos(lon_ref); slon = sin(lon_ref)
    N_loc = a_body / sqrt(1.0d0 - e2_loc * slat * slat)
    ref_ecef(1) = (N_loc + alt_ref) * clat * clon
    ref_ecef(2) = (N_loc + alt_ref) * clat * slon
    ref_ecef(3) = (N_loc * (1.0d0 - e2_loc) + alt_ref) * slat

    dx(1) = x_ecef(1) - ref_ecef(1)
    dx(2) = x_ecef(2) - ref_ecef(2)
    dx(3) = x_ecef(3) - ref_ecef(3)

    ! ECEF→ENU rotation
    x_enu(1) = -slon * dx(1) + clon * dx(2)
    x_enu(2) = -slat * clon * dx(1) - slat * slon * dx(2) + clat * dx(3)
    x_enu(3) =  clat * clon * dx(1) + clat * slon * dx(2) + slat * dx(3)
end subroutine

! ==========================================================================
! fa_coords_enu_to_ecef — ENU to ECEF
! ==========================================================================
subroutine fa_coords_enu_to_ecef(x_enu, x_ecef, lat_ref, lon_ref, alt_ref, &
                                  a_body, f_body, info) &
    bind(C, name="fa_coords_enu_to_ecef")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_enu(3)
    real(c_double), intent(out)        :: x_ecef(3)
    real(c_double), value, intent(in)  :: lat_ref, lon_ref, alt_ref, a_body, f_body
    integer(c_int), intent(out)        :: info

    real(c_double) :: ref_ecef(3), dx(3)
    real(c_double) :: clat, slat, clon, slon, e2_loc, N_loc

    info = 0

    ! Inline geodetic→ECEF for reference point
    e2_loc = 2.0d0 * f_body - f_body * f_body
    clat = cos(lat_ref); slat = sin(lat_ref)
    clon = cos(lon_ref); slon = sin(lon_ref)
    N_loc = a_body / sqrt(1.0d0 - e2_loc * slat * slat)
    ref_ecef(1) = (N_loc + alt_ref) * clat * clon
    ref_ecef(2) = (N_loc + alt_ref) * clat * slon
    ref_ecef(3) = (N_loc * (1.0d0 - e2_loc) + alt_ref) * slat

    ! ENU→ECEF rotation (transpose)
    dx(1) = -slon * x_enu(1) - slat * clon * x_enu(2) + clat * clon * x_enu(3)
    dx(2) =  clon * x_enu(1) - slat * slon * x_enu(2) + clat * slon * x_enu(3)
    dx(3) =                     clat * x_enu(2)        + slat * x_enu(3)

    x_ecef(1) = ref_ecef(1) + dx(1)
    x_ecef(2) = ref_ecef(2) + dx(2)
    x_ecef(3) = ref_ecef(3) + dx(3)
end subroutine

! ==========================================================================
! fa_coords_eci_to_lvlh — ECI to LVLH (RSW convention)
! ==========================================================================
subroutine fa_coords_eci_to_lvlh(x_eci, x_lvlh, r_ref, v_ref, info) &
    bind(C, name="fa_coords_eci_to_lvlh")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_eci(3), r_ref(3), v_ref(3)
    real(c_double), intent(out)        :: x_lvlh(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: R_hat(3), W_hat(3), S_hat(3)
    real(c_double) :: r_mag, h(3), h_mag

    info = 0

    r_mag = sqrt(r_ref(1)**2 + r_ref(2)**2 + r_ref(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        x_lvlh = 0.0d0
        return
    end if

    ! R = radial
    R_hat(1) = r_ref(1) / r_mag
    R_hat(2) = r_ref(2) / r_mag
    R_hat(3) = r_ref(3) / r_mag

    ! h = r x v (angular momentum)
    h(1) = r_ref(2) * v_ref(3) - r_ref(3) * v_ref(2)
    h(2) = r_ref(3) * v_ref(1) - r_ref(1) * v_ref(3)
    h(3) = r_ref(1) * v_ref(2) - r_ref(2) * v_ref(1)
    h_mag = sqrt(h(1)**2 + h(2)**2 + h(3)**2)

    if (h_mag < 1.0d-12) then
        info = 3
        x_lvlh = 0.0d0
        return
    end if

    ! W = cross-track
    W_hat(1) = h(1) / h_mag
    W_hat(2) = h(2) / h_mag
    W_hat(3) = h(3) / h_mag

    ! S = along-track (W x R)
    S_hat(1) = W_hat(2) * R_hat(3) - W_hat(3) * R_hat(2)
    S_hat(2) = W_hat(3) * R_hat(1) - W_hat(1) * R_hat(3)
    S_hat(3) = W_hat(1) * R_hat(2) - W_hat(2) * R_hat(1)

    ! Project into LVLH
    x_lvlh(1) = R_hat(1) * x_eci(1) + R_hat(2) * x_eci(2) + R_hat(3) * x_eci(3)
    x_lvlh(2) = S_hat(1) * x_eci(1) + S_hat(2) * x_eci(2) + S_hat(3) * x_eci(3)
    x_lvlh(3) = W_hat(1) * x_eci(1) + W_hat(2) * x_eci(2) + W_hat(3) * x_eci(3)
end subroutine

! ==========================================================================
! fa_coords_lvlh_to_eci — LVLH to ECI
! ==========================================================================
subroutine fa_coords_lvlh_to_eci(x_lvlh, x_eci, r_ref, v_ref, info) &
    bind(C, name="fa_coords_lvlh_to_eci")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_lvlh(3), r_ref(3), v_ref(3)
    real(c_double), intent(out)        :: x_eci(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: R_hat(3), W_hat(3), S_hat(3)
    real(c_double) :: r_mag, h(3), h_mag

    info = 0

    r_mag = sqrt(r_ref(1)**2 + r_ref(2)**2 + r_ref(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        x_eci = 0.0d0
        return
    end if

    R_hat = r_ref / r_mag

    h(1) = r_ref(2) * v_ref(3) - r_ref(3) * v_ref(2)
    h(2) = r_ref(3) * v_ref(1) - r_ref(1) * v_ref(3)
    h(3) = r_ref(1) * v_ref(2) - r_ref(2) * v_ref(1)
    h_mag = sqrt(h(1)**2 + h(2)**2 + h(3)**2)

    if (h_mag < 1.0d-12) then
        info = 3
        x_eci = 0.0d0
        return
    end if

    W_hat = h / h_mag

    S_hat(1) = W_hat(2) * R_hat(3) - W_hat(3) * R_hat(2)
    S_hat(2) = W_hat(3) * R_hat(1) - W_hat(1) * R_hat(3)
    S_hat(3) = W_hat(1) * R_hat(2) - W_hat(2) * R_hat(1)

    ! Transpose: ECI = R^T * LVLH
    x_eci(1) = R_hat(1) * x_lvlh(1) + S_hat(1) * x_lvlh(2) + W_hat(1) * x_lvlh(3)
    x_eci(2) = R_hat(2) * x_lvlh(1) + S_hat(2) * x_lvlh(2) + W_hat(2) * x_lvlh(3)
    x_eci(3) = R_hat(3) * x_lvlh(1) + S_hat(3) * x_lvlh(2) + W_hat(3) * x_lvlh(3)
end subroutine

! ==========================================================================
! fa_coords_cart_to_keplerian — Cartesian to classical orbital elements
! ==========================================================================
subroutine fa_coords_cart_to_keplerian(rv, mu, oe, info) &
    bind(C, name="fa_coords_cart_to_keplerian")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: rv(6)
    real(c_double), value, intent(in)  :: mu
    real(c_double), intent(out)        :: oe(6)
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: r(3), v(3), h(3), n(3), e_vec(3)
    real(c_double) :: r_mag, v_mag, h_mag, n_mag, e_mag, energy
    real(c_double) :: a, ecc, inc, raan, omega, nu

    info = 0
    r = rv(1:3); v = rv(4:6)

    r_mag = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    v_mag = sqrt(v(1)**2 + v(2)**2 + v(3)**2)

    if (r_mag < 1.0d-12 .or. mu <= 0.0d0) then
        info = 3
        oe = 0.0d0
        return
    end if

    ! Angular momentum h = r x v
    h(1) = r(2) * v(3) - r(3) * v(2)
    h(2) = r(3) * v(1) - r(1) * v(3)
    h(3) = r(1) * v(2) - r(2) * v(1)
    h_mag = sqrt(h(1)**2 + h(2)**2 + h(3)**2)

    ! Node vector n = k x h
    n(1) = -h(2)
    n(2) =  h(1)
    n(3) = 0.0d0
    n_mag = sqrt(n(1)**2 + n(2)**2)

    ! Eccentricity vector e = (v x h)/mu - r/|r|
    e_vec(1) = (v(2) * h(3) - v(3) * h(2)) / mu - r(1) / r_mag
    e_vec(2) = (v(3) * h(1) - v(1) * h(3)) / mu - r(2) / r_mag
    e_vec(3) = (v(1) * h(2) - v(2) * h(1)) / mu - r(3) / r_mag
    e_mag = sqrt(e_vec(1)**2 + e_vec(2)**2 + e_vec(3)**2)

    ! Specific energy
    energy = v_mag**2 / 2.0d0 - mu / r_mag

    ! Semi-major axis
    if (abs(1.0d0 - e_mag) > 1.0d-10) then
        a = -mu / (2.0d0 * energy)
    else
        a = h_mag**2 / mu  ! parabolic: use semi-latus rectum
    end if

    ! Eccentricity
    ecc = e_mag

    ! Inclination
    inc = acos(max(-1.0d0, min(1.0d0, h(3) / h_mag)))

    ! RAAN
    if (n_mag > 1.0d-12) then
        raan = acos(max(-1.0d0, min(1.0d0, n(1) / n_mag)))
        if (n(2) < 0.0d0) raan = 2.0d0 * PI - raan
    else
        raan = 0.0d0
    end if

    ! Argument of periapsis
    if (n_mag > 1.0d-12 .and. e_mag > 1.0d-12) then
        omega = acos(max(-1.0d0, min(1.0d0, &
                (n(1) * e_vec(1) + n(2) * e_vec(2)) / (n_mag * e_mag))))
        if (e_vec(3) < 0.0d0) omega = 2.0d0 * PI - omega
    else
        omega = 0.0d0
    end if

    ! True anomaly
    if (e_mag > 1.0d-12) then
        nu = acos(max(-1.0d0, min(1.0d0, &
             (e_vec(1) * r(1) + e_vec(2) * r(2) + e_vec(3) * r(3)) / (e_mag * r_mag))))
        if (r(1) * v(1) + r(2) * v(2) + r(3) * v(3) < 0.0d0) nu = 2.0d0 * PI - nu
    else
        nu = 0.0d0
    end if

    oe(1) = a
    oe(2) = ecc
    oe(3) = inc
    oe(4) = raan
    oe(5) = omega
    oe(6) = nu
end subroutine

! ==========================================================================
! fa_coords_keplerian_to_cart — Classical orbital elements to Cartesian
! ==========================================================================
subroutine fa_coords_keplerian_to_cart(oe, mu, rv, info) &
    bind(C, name="fa_coords_keplerian_to_cart")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: oe(6)
    real(c_double), value, intent(in)  :: mu
    real(c_double), intent(out)        :: rv(6)
    integer(c_int), intent(out)        :: info

    real(c_double) :: a, ecc, inc, raan, omega, nu
    real(c_double) :: p, r_mag, r_pqw(3), v_pqw(3)
    real(c_double) :: cos_raan, sin_raan, cos_aop, sin_aop, cos_inc, sin_inc, cos_nu, sin_nu

    info = 0
    a = oe(1); ecc = oe(2); inc = oe(3)
    raan = oe(4); omega = oe(5); nu = oe(6)

    if (a <= 0.0d0 .or. ecc < 0.0d0 .or. ecc >= 1.0d0 .or. mu <= 0.0d0) then
        info = 3
        rv = 0.0d0
        return
    end if

    ! Semi-latus rectum
    p = a * (1.0d0 - ecc * ecc)
    cos_nu = cos(nu); sin_nu = sin(nu)
    r_mag = p / (1.0d0 + ecc * cos_nu)

    ! Position and velocity in perifocal frame
    r_pqw(1) = r_mag * cos_nu
    r_pqw(2) = r_mag * sin_nu
    r_pqw(3) = 0.0d0

    v_pqw(1) = -sqrt(mu / p) * sin_nu
    v_pqw(2) =  sqrt(mu / p) * (ecc + cos_nu)
    v_pqw(3) = 0.0d0

    ! Rotation PQW → ECI
    cos_raan = cos(raan); sin_raan = sin(raan)
    cos_aop = cos(omega); sin_aop = sin(omega)
    cos_inc = cos(inc);   sin_inc = sin(inc)

    ! Position
    rv(1) = (cos_raan*cos_aop - sin_raan*sin_aop*cos_inc) * r_pqw(1) &
          + (-cos_raan*sin_aop - sin_raan*cos_aop*cos_inc) * r_pqw(2)
    rv(2) = (sin_raan*cos_aop + cos_raan*sin_aop*cos_inc) * r_pqw(1) &
          + (-sin_raan*sin_aop + cos_raan*cos_aop*cos_inc) * r_pqw(2)
    rv(3) = (sin_aop*sin_inc) * r_pqw(1) + (cos_aop*sin_inc) * r_pqw(2)

    ! Velocity
    rv(4) = (cos_raan*cos_aop - sin_raan*sin_aop*cos_inc) * v_pqw(1) &
          + (-cos_raan*sin_aop - sin_raan*cos_aop*cos_inc) * v_pqw(2)
    rv(5) = (sin_raan*cos_aop + cos_raan*sin_aop*cos_inc) * v_pqw(1) &
          + (-sin_raan*sin_aop + cos_raan*cos_aop*cos_inc) * v_pqw(2)
    rv(6) = (sin_aop*sin_inc) * v_pqw(1) + (cos_aop*sin_inc) * v_pqw(2)
end subroutine

! ==========================================================================
! fa_coords_cart_to_equinoctial — Cartesian to equinoctial elements
! ==========================================================================
subroutine fa_coords_cart_to_equinoctial(rv, mu, eq, info) &
    bind(C, name="fa_coords_cart_to_equinoctial")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: rv(6)
    real(c_double), value, intent(in)  :: mu
    real(c_double), intent(out)        :: eq(6)
    integer(c_int), intent(out)        :: info

    interface
        subroutine fa_coords_cart_to_keplerian(rv, mu, oe, info) &
            bind(C, name="fa_coords_cart_to_keplerian")
            use iso_c_binding
            real(c_double), intent(in)         :: rv(6)
            real(c_double), value, intent(in)  :: mu
            real(c_double), intent(out)        :: oe(6)
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: oe(6), a, ecc, inc, raan, omega, nu, Lg

    ! First convert to classical elements, then to equinoctial
    call fa_coords_cart_to_keplerian(rv, mu, oe, info)
    if (info /= 0) return

    a = oe(1); ecc = oe(2); inc = oe(3)
    raan = oe(4); omega = oe(5); nu = oe(6)

    ! p = semi-latus rectum
    eq(1) = a * (1.0d0 - ecc * ecc)
    ! f = e*cos(omega + RAAN)
    eq(2) = ecc * cos(omega + raan)
    ! g = e*sin(omega + RAAN)
    eq(3) = ecc * sin(omega + raan)
    ! h = tan(i/2)*cos(RAAN)
    eq(4) = tan(inc / 2.0d0) * cos(raan)
    ! k = tan(i/2)*sin(RAAN)
    eq(5) = tan(inc / 2.0d0) * sin(raan)
    ! L = true longitude = RAAN + omega + nu
    Lg = raan + omega + nu
    eq(6) = mod(Lg, 2.0d0 * PI)
    if (eq(6) < 0.0d0) eq(6) = eq(6) + 2.0d0 * PI
end subroutine

! ==========================================================================
! fa_coords_equinoctial_to_cart — Equinoctial elements to Cartesian
! ==========================================================================
subroutine fa_coords_equinoctial_to_cart(eq, mu, rv, info) &
    bind(C, name="fa_coords_equinoctial_to_cart")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: eq(6)
    real(c_double), value, intent(in)  :: mu
    real(c_double), intent(out)        :: rv(6)
    integer(c_int), intent(out)        :: info

    real(c_double) :: p, f_eq, g_eq, h_eq, k_eq, Lg
    real(c_double) :: alpha2, s2, w, r_mag, cL, sL
    real(c_double) :: r1, r2, v1, v2
    real(c_double) :: fx, fy, fz, gx, gy, gz

    info = 0
    p = eq(1); f_eq = eq(2); g_eq = eq(3)
    h_eq = eq(4); k_eq = eq(5); Lg = eq(6)

    if (p <= 0.0d0 .or. mu <= 0.0d0) then
        info = 3
        rv = 0.0d0
        return
    end if

    alpha2 = h_eq**2 - k_eq**2
    s2 = 1.0d0 + h_eq**2 + k_eq**2
    cL = cos(Lg); sL = sin(Lg)
    w = 1.0d0 + f_eq * cL + g_eq * sL
    r_mag = p / w

    r1 = r_mag * cL
    r2 = r_mag * sL

    ! Unit vectors f, g in equinoctial frame
    fx = (1.0d0 + alpha2) / s2
    fy = 2.0d0 * h_eq * k_eq / s2
    fz = -2.0d0 * k_eq / s2

    gx = 2.0d0 * h_eq * k_eq / s2
    gy = (1.0d0 - alpha2) / s2
    gz = 2.0d0 * h_eq / s2

    ! Position
    rv(1) = r1 * fx + r2 * gx
    rv(2) = r1 * fy + r2 * gy
    rv(3) = r1 * fz + r2 * gz

    ! Velocity components in equinoctial frame
    v1 = -sqrt(mu / p) * (sL + g_eq)
    v2 =  sqrt(mu / p) * (cL + f_eq)

    rv(4) = v1 * fx + v2 * gx
    rv(5) = v1 * fy + v2 * gy
    rv(6) = v1 * fz + v2 * gz
end subroutine

! ==========================================================================
! fa_coords_rsw_matrix — RSW rotation matrix from r, v
! ==========================================================================
subroutine fa_coords_rsw_matrix(r, v, R_rsw, info) &
    bind(C, name="fa_coords_rsw_matrix")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r(3), v(3)
    real(c_double), intent(out)        :: R_rsw(9)
    integer(c_int), intent(out)        :: info

    real(c_double) :: R_hat(3), W_hat(3), S_hat(3)
    real(c_double) :: r_mag, h(3), h_mag

    info = 0
    r_mag = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        R_rsw = 0.0d0
        return
    end if

    R_hat = r / r_mag

    h(1) = r(2) * v(3) - r(3) * v(2)
    h(2) = r(3) * v(1) - r(1) * v(3)
    h(3) = r(1) * v(2) - r(2) * v(1)
    h_mag = sqrt(h(1)**2 + h(2)**2 + h(3)**2)

    if (h_mag < 1.0d-12) then
        info = 3
        R_rsw = 0.0d0
        return
    end if

    W_hat = h / h_mag

    S_hat(1) = W_hat(2) * R_hat(3) - W_hat(3) * R_hat(2)
    S_hat(2) = W_hat(3) * R_hat(1) - W_hat(1) * R_hat(3)
    S_hat(3) = W_hat(1) * R_hat(2) - W_hat(2) * R_hat(1)

    ! Row-major 3x3: row 0 = R, row 1 = S, row 2 = W
    R_rsw(1) = R_hat(1); R_rsw(2) = R_hat(2); R_rsw(3) = R_hat(3)
    R_rsw(4) = S_hat(1); R_rsw(5) = S_hat(2); R_rsw(6) = S_hat(3)
    R_rsw(7) = W_hat(1); R_rsw(8) = W_hat(2); R_rsw(9) = W_hat(3)
end subroutine

! ==========================================================================
! fa_coords_vnc_matrix — VNC rotation matrix from r, v
! ==========================================================================
subroutine fa_coords_vnc_matrix(r, v, R_vnc, info) &
    bind(C, name="fa_coords_vnc_matrix")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r(3), v(3)
    real(c_double), intent(out)        :: R_vnc(9)
    integer(c_int), intent(out)        :: info

    real(c_double) :: V_hat(3), N_hat(3), C_hat(3)
    real(c_double) :: v_mag, h(3), h_mag

    info = 0
    v_mag = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
    if (v_mag < 1.0d-12) then
        info = 3
        R_vnc = 0.0d0
        return
    end if

    V_hat = v / v_mag

    h(1) = r(2) * v(3) - r(3) * v(2)
    h(2) = r(3) * v(1) - r(1) * v(3)
    h(3) = r(1) * v(2) - r(2) * v(1)
    h_mag = sqrt(h(1)**2 + h(2)**2 + h(3)**2)

    if (h_mag < 1.0d-12) then
        info = 3
        R_vnc = 0.0d0
        return
    end if

    N_hat = h / h_mag

    ! C = V x N
    C_hat(1) = V_hat(2) * N_hat(3) - V_hat(3) * N_hat(2)
    C_hat(2) = V_hat(3) * N_hat(1) - V_hat(1) * N_hat(3)
    C_hat(3) = V_hat(1) * N_hat(2) - V_hat(2) * N_hat(1)

    ! Row-major 3x3: row 0 = V, row 1 = N, row 2 = C
    R_vnc(1) = V_hat(1); R_vnc(2) = V_hat(2); R_vnc(3) = V_hat(3)
    R_vnc(4) = N_hat(1); R_vnc(5) = N_hat(2); R_vnc(6) = N_hat(3)
    R_vnc(7) = C_hat(1); R_vnc(8) = C_hat(2); R_vnc(9) = C_hat(3)
end subroutine

! ==========================================================================
! fa_coords_rotation — Get 3x3 rotation matrix between frame IDs
! ==========================================================================
subroutine fa_coords_rotation(from_id, to_id, t, params, np, R, info) &
    bind(C, name="fa_coords_rotation")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: from_id, to_id, np
    real(c_double), value, intent(in)  :: t
    real(c_double), intent(in)         :: params(np)
    real(c_double), intent(out)        :: R(9)
    integer(c_int), intent(out)        :: info

    real(c_double) :: cg, sg

    info = 0
    R = 0.0d0

    ! ECI (1) ↔ ECEF (10): params(1) = gmst or compute from t
    if ((from_id == 1 .and. to_id == 10) .or. (from_id == 10 .and. to_id == 1)) then
        if (np >= 1) then
            cg = cos(params(1))
            sg = sin(params(1))
        else
            info = 3
            return
        end if

        if (from_id == 1 .and. to_id == 10) then
            ! Rz(gmst)
            R(1) =  cg; R(2) = sg; R(3) = 0.0d0
            R(4) = -sg; R(5) = cg; R(6) = 0.0d0
            R(7) = 0.0d0; R(8) = 0.0d0; R(9) = 1.0d0
        else
            ! Rz(-gmst)
            R(1) = cg; R(2) = -sg; R(3) = 0.0d0
            R(4) = sg; R(5) =  cg; R(6) = 0.0d0
            R(7) = 0.0d0; R(8) = 0.0d0; R(9) = 1.0d0
        end if
    else if (from_id == to_id) then
        ! Identity
        R(1) = 1.0d0; R(5) = 1.0d0; R(9) = 1.0d0
    else
        info = 3  ! unsupported frame pair
    end if
end subroutine

! ==========================================================================
! fa_coords_transform — Master transform dispatch
! ==========================================================================
subroutine fa_coords_transform(from_id, to_id, n, x_in, x_out, t, params, np, info) &
    bind(C, name="fa_coords_transform")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: from_id, to_id, n, np
    real(c_double), intent(in)         :: x_in(n), params(np)
    real(c_double), intent(out)        :: x_out(n)
    real(c_double), value, intent(in)  :: t
    integer(c_int), intent(out)        :: info

    interface
        subroutine fa_coords_rotation(fid, tid, tt, par, npp, rot, inf) &
            bind(C, name="fa_coords_rotation")
            use iso_c_binding
            integer(c_int), value, intent(in) :: fid, tid, npp
            real(c_double), value, intent(in) :: tt
            real(c_double), intent(in)        :: par(npp)
            real(c_double), intent(out)       :: rot(9)
            integer(c_int), intent(out)       :: inf
        end subroutine
    end interface

    real(c_double) :: R(9)
    integer :: sub_info

    info = 0

    if (n < 3) then
        info = 3
        x_out(1:n) = 0.0d0
        return
    end if

    if (from_id == to_id) then
        x_out(1:n) = x_in(1:n)
        return
    end if

    ! Get rotation matrix
    call fa_coords_rotation(from_id, to_id, t, params, np, R, sub_info)
    if (sub_info /= 0) then
        info = sub_info
        x_out(1:n) = 0.0d0
        return
    end if

    ! Apply rotation to position (first 3 components)
    x_out(1) = R(1)*x_in(1) + R(2)*x_in(2) + R(3)*x_in(3)
    x_out(2) = R(4)*x_in(1) + R(5)*x_in(2) + R(6)*x_in(3)
    x_out(3) = R(7)*x_in(1) + R(8)*x_in(2) + R(9)*x_in(3)

    ! If velocity present (n >= 6), rotate velocity too
    if (n >= 6) then
        x_out(4) = R(1)*x_in(4) + R(2)*x_in(5) + R(3)*x_in(6)
        x_out(5) = R(4)*x_in(4) + R(5)*x_in(5) + R(6)*x_in(6)
        x_out(6) = R(7)*x_in(4) + R(8)*x_in(5) + R(9)*x_in(6)
    end if

    ! Copy remaining components unchanged
    if (n > 6) then
        x_out(7:n) = x_in(7:n)
    end if
end subroutine
