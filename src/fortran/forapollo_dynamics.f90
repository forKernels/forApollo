! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_dynamics.f90 — Tracking dynamics models with analytic Jacobians
!
! Built-in dynamics catalog for state estimation. Each model provides:
!   - f(x): continuous-time state derivative x_dot = f(x, u, t)
!   - F(x): analytic Jacobian df/dx (flat row-major n×n)
!
! Model IDs (tracking group):
!   40 = Constant velocity   — state [pos(d), vel(d)], n = 2*d
!   41 = Constant acceleration — state [pos(d), vel(d), acc(d)], n = 3*d
!   42 = Constant turn rate  — state [x, y, vx, vy, omega], n = 5
!
! Model IDs (orbital group):
!    1 = Kepler two-body      — state [r(3), v(3)], n=6, params=[mu]
!    2 = J2 oblateness         — state [r(3), v(3)], n=6, params=[mu, J2, R_eq]
!    3 = CR3BP rotating frame  — state [r(3), v(3)], n=6, params=[mu_ratio]
!    4 = Atmospheric drag      — state [r(3), v(3), beta], n=7, params=[mu, rho0, h_scale, R_body]
!
! All bind(C) routines use value for scalar args, assumed-size for arrays.
! Internal model subroutines do NOT have bind(C) and do NOT use value.
! Matrices are flat 1D row-major: element (row, col) = F((row-1)*n + col).
!
! Error codes: info = 0 ok, 3 invalid input.

subroutine fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info) &
    bind(C, name="fa_dynamics_dispatch")
    use iso_c_binding
    implicit none

    integer(c_int), value  :: model_id, n, nu, np
    real(c_double), value  :: t
    real(c_double), intent(in)  :: x(n)
    real(c_double), intent(in)  :: u(nu)
    real(c_double), intent(in)  :: params(np)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info

    info = 0
    x_dot(1:n) = 0.0d0

    select case (model_id)
    case (1)
        call dynamics_kepler(n, x, params, np, x_dot, info)
    case (2)
        call dynamics_j2(n, x, params, np, x_dot, info)
    case (3)
        call dynamics_cr3bp(n, x, params, np, x_dot, info)
    case (4)
        call dynamics_drag(n, x, params, np, x_dot, info)
    case (40)
        call dynamics_const_vel(n, x, x_dot, info)
    case (41)
        call dynamics_const_accel(n, x, x_dot, info)
    case (42)
        call dynamics_const_turn(n, x, x_dot, info)
    case default
        info = 3
        return
    end select

contains

    ! =========================================================================
    ! Constant Velocity (model_id = 40)
    ! State: [pos_1, ..., pos_d, vel_1, ..., vel_d], n = 2*d
    ! x_dot = [vel_1, ..., vel_d, 0, ..., 0]
    ! =========================================================================
    subroutine dynamics_const_vel(n, x, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        integer :: d

        ! n must be even (n = 2*d)
        if (mod(n, 2) /= 0 .or. n < 2) then
            info = 3
            return
        end if

        info = 0
        d = n / 2

        ! Position derivatives = velocities
        x_dot(1:d) = x(d+1 : n)

        ! Velocity derivatives = 0 (constant velocity assumption)
        x_dot(d+1 : n) = 0.0d0

    end subroutine dynamics_const_vel

    ! =========================================================================
    ! Constant Acceleration (model_id = 41)
    ! State: [pos_1..d, vel_1..d, acc_1..d], n = 3*d
    ! x_dot = [vel_1..d, acc_1..d, 0..0]
    ! =========================================================================
    subroutine dynamics_const_accel(n, x, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        integer :: d

        ! n must be divisible by 3 (n = 3*d)
        if (mod(n, 3) /= 0 .or. n < 3) then
            info = 3
            return
        end if

        info = 0
        d = n / 3

        ! Position derivatives = velocities
        x_dot(1:d) = x(d+1 : 2*d)

        ! Velocity derivatives = accelerations
        x_dot(d+1 : 2*d) = x(2*d+1 : n)

        ! Acceleration derivatives = 0 (constant acceleration assumption)
        x_dot(2*d+1 : n) = 0.0d0

    end subroutine dynamics_const_accel

    ! =========================================================================
    ! Constant Turn Rate (model_id = 42)
    ! State: [x, y, vx, vy, omega], n = 5 exactly
    ! x_dot = [vx, vy, -omega*vy, omega*vx, 0]
    !
    ! This is the coordinated turn model in 2D where the object moves with
    ! constant speed and constant angular rate omega.
    ! =========================================================================
    subroutine dynamics_const_turn(n, x, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: vx, vy, omega

        ! n must be exactly 5
        if (n /= 5) then
            info = 3
            return
        end if

        info = 0

        vx    = x(3)
        vy    = x(4)
        omega = x(5)

        x_dot(1) = vx                ! dx/dt = vx
        x_dot(2) = vy                ! dy/dt = vy
        x_dot(3) = -omega * vy       ! dvx/dt = -omega * vy
        x_dot(4) =  omega * vx       ! dvy/dt =  omega * vx
        x_dot(5) = 0.0d0             ! domega/dt = 0 (constant turn rate)

    end subroutine dynamics_const_turn

    ! =========================================================================
    ! Kepler Two-Body (model_id = 1)
    ! State: [rx, ry, rz, vx, vy, vz], n=6
    ! Params: [mu], np >= 1
    ! x_dot = [v; -mu/|r|^3 * r]
    ! info=1 if |r| < 1e-10 (degenerate), info=3 if n/=6 or np<1
    ! =========================================================================
    subroutine dynamics_kepler(n, x, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, rx, ry, rz, r, r3

        if (n /= 6 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        mu = params(1)
        rx = x(1); ry = x(2); rz = x(3)
        r = sqrt(rx*rx + ry*ry + rz*rz)

        if (r < 1.0d-10) then
            info = 1
            x_dot(1:6) = 0.0d0
            return
        end if

        r3 = r * r * r

        ! Position derivatives = velocity
        x_dot(1) = x(4)
        x_dot(2) = x(5)
        x_dot(3) = x(6)

        ! Velocity derivatives = gravitational acceleration
        x_dot(4) = -mu * rx / r3
        x_dot(5) = -mu * ry / r3
        x_dot(6) = -mu * rz / r3

    end subroutine dynamics_kepler

    ! =========================================================================
    ! J2 Oblateness Perturbation (model_id = 2)
    ! State: [rx, ry, rz, vx, vy, vz], n=6
    ! Params: [mu, J2, R_eq], np >= 3
    ! Kepler + J2 perturbation acceleration
    ! =========================================================================
    subroutine dynamics_j2(n, x, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, J2, R_eq, rx, ry, rz, r, r2, r5
        real(c_double) :: z2_over_r2, coeff

        if (n /= 6 .or. np < 3) then
            info = 3
            return
        end if

        ! Start with Kepler base dynamics
        call dynamics_kepler(n, x, params, np, x_dot, info)
        if (info /= 0) return

        mu   = params(1)
        J2   = params(2)
        R_eq = params(3)

        rx = x(1); ry = x(2); rz = x(3)
        r2 = rx*rx + ry*ry + rz*rz
        r  = sqrt(r2)
        r5 = r2 * r2 * r

        z2_over_r2 = rz * rz / r2
        coeff = -1.5d0 * J2 * mu * R_eq * R_eq / r5

        ! J2 perturbation accelerations added to Kepler
        x_dot(4) = x_dot(4) + coeff * rx * (1.0d0 - 5.0d0 * z2_over_r2)
        x_dot(5) = x_dot(5) + coeff * ry * (1.0d0 - 5.0d0 * z2_over_r2)
        x_dot(6) = x_dot(6) + coeff * rz * (3.0d0 - 5.0d0 * z2_over_r2)

    end subroutine dynamics_j2

    ! =========================================================================
    ! Circular Restricted Three-Body Problem (model_id = 3)
    ! State: [x, y, z, vx, vy, vz], n=6 in rotating frame
    ! Params: [mu_ratio] where mu = m2/(m1+m2), np >= 1
    !
    ! Primary m1 at (-mu, 0, 0), secondary m2 at (1-mu, 0, 0)
    ! Equations in rotating frame with non-dimensional units:
    !   vx_dot = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3
    !   vy_dot = -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3
    !   vz_dot = -(1-mu)*z/r1^3 - mu*z/r2^3
    ! =========================================================================
    subroutine dynamics_cr3bp(n, x, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, mu1, px, py, pz, pvx, pvy, pvz
        real(c_double) :: d1x, d2x, r1_2, r2_2, r1_3, r2_3

        if (n /= 6 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        mu  = params(1)         ! mass ratio m2/(m1+m2)
        mu1 = 1.0d0 - mu        ! 1 - mu

        px  = x(1); py  = x(2); pz  = x(3)
        pvx = x(4); pvy = x(5); pvz = x(6)

        ! Distances to primaries
        d1x = px + mu            ! x - (-mu) = x + mu
        d2x = px - mu1           ! x - (1-mu)

        r1_2 = d1x*d1x + py*py + pz*pz
        r2_2 = d2x*d2x + py*py + pz*pz

        if (r1_2 < 1.0d-20 .or. r2_2 < 1.0d-20) then
            info = 1
            x_dot(1:6) = 0.0d0
            return
        end if

        r1_3 = r1_2 * sqrt(r1_2)   ! r1^3
        r2_3 = r2_2 * sqrt(r2_2)   ! r2^3

        ! Position derivatives = velocity
        x_dot(1) = pvx
        x_dot(2) = pvy
        x_dot(3) = pvz

        ! Velocity derivatives in rotating frame
        x_dot(4) = 2.0d0*pvy + px - mu1*d1x/r1_3 - mu*d2x/r2_3
        x_dot(5) = -2.0d0*pvx + py - mu1*py/r1_3 - mu*py/r2_3
        x_dot(6) = -mu1*pz/r1_3 - mu*pz/r2_3

    end subroutine dynamics_cr3bp

    ! =========================================================================
    ! Atmospheric Drag (model_id = 4)
    ! State: [rx, ry, rz, vx, vy, vz, beta], n=7
    !   beta = Cd*A/m (ballistic coefficient, estimated as state)
    ! Params: [mu, rho0, h_scale, R_body], np >= 4
    !
    ! Kepler gravity + exponential atmospheric drag:
    !   rho = rho0 * exp(-(|r| - R_body) / h_scale)
    !   a_drag = -0.5 * rho * beta * |v| * v
    !   beta_dot = 0 (constant ballistic coefficient)
    ! =========================================================================
    subroutine dynamics_drag(n, x, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, rho0, h_scale, R_body, beta
        real(c_double) :: rx, ry, rz, vx, vy, vz
        real(c_double) :: r, r3, v_mag, rho, drag_coeff

        if (n /= 7 .or. np < 4) then
            info = 3
            return
        end if

        info = 0
        mu      = params(1)
        rho0    = params(2)
        h_scale = params(3)
        R_body  = params(4)

        rx = x(1); ry = x(2); rz = x(3)
        vx = x(4); vy = x(5); vz = x(6)
        beta = x(7)

        r = sqrt(rx*rx + ry*ry + rz*rz)

        if (r < 1.0d-10) then
            info = 1
            x_dot(1:7) = 0.0d0
            return
        end if

        r3 = r * r * r
        v_mag = sqrt(vx*vx + vy*vy + vz*vz)

        ! Atmospheric density (exponential model)
        rho = rho0 * exp(-(r - R_body) / h_scale)

        ! Drag coefficient: -0.5 * rho * beta * |v|
        drag_coeff = -0.5d0 * rho * beta * v_mag

        ! Position derivatives = velocity
        x_dot(1) = vx
        x_dot(2) = vy
        x_dot(3) = vz

        ! Velocity derivatives = Kepler gravity + drag
        x_dot(4) = -mu * rx / r3 + drag_coeff * vx
        x_dot(5) = -mu * ry / r3 + drag_coeff * vy
        x_dot(6) = -mu * rz / r3 + drag_coeff * vz

        ! Ballistic coefficient is constant
        x_dot(7) = 0.0d0

    end subroutine dynamics_drag

end subroutine fa_dynamics_dispatch


! =============================================================================
! Jacobian dispatch — analytic df/dx for each built-in dynamics model
! F is flat row-major n×n: element (i,j) stored at F((i-1)*n + j)
! =============================================================================
subroutine fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np, F, info) &
    bind(C, name="fa_dynamics_jacobian")
    use iso_c_binding
    implicit none

    integer(c_int), value  :: model_id, n, nu, np
    real(c_double), value  :: t
    real(c_double), intent(in)  :: x(n)
    real(c_double), intent(in)  :: u(nu)
    real(c_double), intent(in)  :: params(np)
    real(c_double), intent(out) :: F(n*n)
    integer(c_int), intent(out) :: info

    info = 0
    F(1:n*n) = 0.0d0

    select case (model_id)
    case (1)
        call jacobian_kepler(n, x, params, np, F, info)
    case (2)
        call jacobian_j2(n, x, params, np, F, info)
    case (3)
        call jacobian_cr3bp(n, x, params, np, F, info)
    case (4)
        call jacobian_drag(n, x, params, np, F, info)
    case (40)
        call jacobian_const_vel(n, x, F, info)
    case (41)
        call jacobian_const_accel(n, x, F, info)
    case (42)
        call jacobian_const_turn(n, x, F, info)
    case default
        info = 3
        return
    end select

contains

    ! =========================================================================
    ! Jacobian for Constant Velocity
    !
    !       [ 0  I ]
    !   F = [ 0  0 ]
    !
    ! where each block is d×d. The identity block maps velocity to position.
    ! Row-major: F((i-1)*n + j)
    ! =========================================================================
    subroutine jacobian_const_vel(n, x, F, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        integer :: d, i

        if (mod(n, 2) /= 0 .or. n < 2) then
            info = 3
            return
        end if

        info = 0
        d = n / 2

        F(1:n*n) = 0.0d0

        ! Upper-right identity block: dpos_i/dvel_i = 1
        ! Row i (1..d), Col d+i → index (i-1)*n + (d+i)
        do i = 1, d
            F((i-1)*n + d + i) = 1.0d0
        end do

    end subroutine jacobian_const_vel

    ! =========================================================================
    ! Jacobian for Constant Acceleration
    !
    !       [ 0  I  0 ]
    !   F = [ 0  0  I ]
    !       [ 0  0  0 ]
    !
    ! Each block is d×d. Identity blocks link pos→vel and vel→acc.
    ! =========================================================================
    subroutine jacobian_const_accel(n, x, F, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        integer :: d, i

        if (mod(n, 3) /= 0 .or. n < 3) then
            info = 3
            return
        end if

        info = 0
        d = n / 3

        F(1:n*n) = 0.0d0

        ! Block (1,2): dpos_i/dvel_i = 1
        ! Row i (1..d), Col d+i
        do i = 1, d
            F((i-1)*n + d + i) = 1.0d0
        end do

        ! Block (2,3): dvel_i/dacc_i = 1
        ! Row d+i (d+1..2d), Col 2d+i
        do i = 1, d
            F((d+i-1)*n + 2*d + i) = 1.0d0
        end do

    end subroutine jacobian_const_accel

    ! =========================================================================
    ! Jacobian for Constant Turn Rate
    !
    ! State: [x, y, vx, vy, omega]
    ! f = [vx, vy, -omega*vy, omega*vx, 0]
    !
    !       [ 0  0   1       0      0     ]
    !   F = [ 0  0   0       1      0     ]
    !       [ 0  0   0      -omega  -vy   ]
    !       [ 0  0   omega   0       vx   ]
    !       [ 0  0   0       0       0    ]
    !
    ! Partials:
    !   df3/dvx = 0, df3/dvy = -omega, df3/domega = -vy
    !   df4/dvx = omega, df4/dvy = 0, df4/domega = vx
    ! =========================================================================
    subroutine jacobian_const_turn(n, x, F, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: vx, vy, omega

        if (n /= 5) then
            info = 3
            return
        end if

        info = 0

        vx    = x(3)
        vy    = x(4)
        omega = x(5)

        F(1:n*n) = 0.0d0

        ! Row 1: dx/dt = vx → df1/dvx = 1 → F(1,3) = (0)*5+3 = 3
        F(3) = 1.0d0

        ! Row 2: dy/dt = vy → df2/dvy = 1 → F(2,4) = (1)*5+4 = 9
        F(9) = 1.0d0

        ! Row 3: dvx/dt = -omega*vy
        !   df3/dvy    = -omega → F(3,4) = (2)*5+4 = 14
        !   df3/domega = -vy    → F(3,5) = (2)*5+5 = 15
        F(14) = -omega
        F(15) = -vy

        ! Row 4: dvy/dt = omega*vx
        !   df4/dvx    =  omega → F(4,3) = (3)*5+3 = 18
        !   df4/domega =  vx    → F(4,5) = (3)*5+5 = 20
        F(18) = omega
        F(20) = vx

        ! Row 5: domega/dt = 0 → all zeros (already set)

    end subroutine jacobian_const_turn

    ! =========================================================================
    ! Jacobian for Kepler Two-Body (model_id = 1)
    !
    ! State: [rx, ry, rz, vx, vy, vz], n=6
    ! F is 6x6 row-major:
    !   Top-left  3x3 = 0 (dr_dot/dr = 0)
    !   Top-right 3x3 = I (dr_dot/dv = I)
    !   Bottom-left 3x3 = gravity gradient tensor
    !   Bottom-right 3x3 = 0 (dv_dot/dv = 0)
    !
    ! Gravity gradient: d(-mu*xi/r^3)/dxj = -mu*(delta_ij*r^2 - 3*xi*xj)/r^5
    ! =========================================================================
    subroutine jacobian_kepler(n, x, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, rx, ry, rz, r2, r, r5
        real(c_double) :: rr(3)
        integer :: i, j, row, col

        if (n /= 6 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        mu = params(1)
        rx = x(1); ry = x(2); rz = x(3)
        rr(1) = rx; rr(2) = ry; rr(3) = rz

        r2 = rx*rx + ry*ry + rz*rz
        r  = sqrt(r2)

        if (r < 1.0d-10) then
            info = 1
            F(1:n*n) = 0.0d0
            return
        end if

        r5 = r2 * r2 * r

        F(1:n*n) = 0.0d0

        ! Top-right 3x3 identity: dr_dot/dv = I
        ! Row i (1..3), Col 3+i → F((i-1)*6 + 3+i)
        do i = 1, 3
            F((i-1)*6 + 3 + i) = 1.0d0
        end do

        ! Bottom-left 3x3: gravity gradient tensor
        ! d(a_i)/d(r_j) = -mu*(delta_ij*r^2 - 3*r_i*r_j)/r^5
        ! Row 3+i, Col j → F((3+i-1)*6 + j)
        do i = 1, 3
            do j = 1, 3
                row = 3 + i
                col = j
                if (i == j) then
                    F((row-1)*6 + col) = -mu * (r2 - 3.0d0*rr(i)*rr(j)) / r5
                else
                    F((row-1)*6 + col) = -mu * (-3.0d0*rr(i)*rr(j)) / r5
                end if
            end do
        end do

    end subroutine jacobian_kepler

    ! =========================================================================
    ! Jacobian for J2 Oblateness (model_id = 2)
    !
    ! Kepler Jacobian + J2 perturbation partial derivatives
    ! The J2 acceleration components:
    !   a_J2_x = C * x * (1 - 5*z^2/r^2)
    !   a_J2_y = C * y * (1 - 5*z^2/r^2)
    !   a_J2_z = C * z * (3 - 5*z^2/r^2)
    ! where C = -1.5*J2*mu*R_eq^2/r^5
    !
    ! Partial derivatives are computed analytically.
    ! =========================================================================
    subroutine jacobian_j2(n, x, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, J2, R_eq, rx, ry, rz
        real(c_double) :: r2, r, r5, r7, z2
        real(c_double) :: C0, t1, t2
        integer :: row, col

        if (n /= 6 .or. np < 3) then
            info = 3
            return
        end if

        ! Start with Kepler Jacobian
        call jacobian_kepler(n, x, params, np, F, info)
        if (info /= 0) return

        mu   = params(1)
        J2   = params(2)
        R_eq = params(3)

        rx = x(1); ry = x(2); rz = x(3)
        r2 = rx*rx + ry*ry + rz*rz
        r  = sqrt(r2)
        r5 = r2 * r2 * r
        r7 = r5 * r2
        z2 = rz * rz

        C0 = -1.5d0 * J2 * mu * R_eq * R_eq

        ! J2 partial derivatives added to bottom-left 3x3 block
        ! These are d(a_J2_i)/d(r_j), added to F(3+i, j)
        !
        ! Using the chain rule on a_J2_x = C0/r^5 * x * (1 - 5*z^2/r^2):
        ! d(a_J2_x)/dx = C0 * [(1-5*z2/r2)/r5 + x*d/dx((1-5*z2/r2)/r5)]
        ! d/dx(1/r5) = -5*x/r7
        ! d/dx(z2/r2) = -2*z2*x/r4
        ! So d/dx((1-5*z2/r2)/r5) = -5*x/r7*(1-5*z2/r2) + (-5)*(-2*z2*x/r4)/r5
        !                          = -5*x/r7*(1-5*z2/r2) + 10*z2*x/(r4*r5)
        !                          = -5*x/r7 + 25*z2*x/r9 + 10*z2*x/r9
        !                          = -5*x/r7 + 35*z2*x/r9
        ! d(a_J2_x)/dx = C0*[(1-5*z2/r2)/r5 + x*(-5*x/r7 + 35*z2*x/r9)]
        !              = C0/r5*[(1-5*z2/r2) - 5*x^2/r2 + 35*z2*x^2/r4]
        !              = C0/r5*[1 - 5*(z2+x^2)/r2 + 35*z2*x^2/r4]
        !
        ! General pattern for d(a_J2_i)/d(r_j):
        ! We compute these term by term.

        ! --- d(a_J2_x)/drx  (row=4, col=1)
        row = 4; col = 1
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * (1.0d0 - 5.0d0*(z2 + rx*rx)/r2 + 35.0d0*z2*rx*rx/(r2*r2))

        ! --- d(a_J2_x)/dry  (row=4, col=2)
        row = 4; col = 2
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * rx * (-5.0d0*ry/r2 + 35.0d0*z2*ry/(r2*r2))

        ! --- d(a_J2_x)/drz  (row=4, col=3)
        row = 4; col = 3
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * rx * (-5.0d0*rz/r2 - 10.0d0*rz/r2 + 35.0d0*z2*rz/(r2*r2))

        ! --- d(a_J2_y)/drx  (row=5, col=1) — symmetric with d(a_J2_x)/dry swapping x↔y
        row = 5; col = 1
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * ry * (-5.0d0*rx/r2 + 35.0d0*z2*rx/(r2*r2))

        ! --- d(a_J2_y)/dry  (row=5, col=2)
        row = 5; col = 2
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * (1.0d0 - 5.0d0*(z2 + ry*ry)/r2 + 35.0d0*z2*ry*ry/(r2*r2))

        ! --- d(a_J2_y)/drz  (row=5, col=3)
        row = 5; col = 3
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * ry * (-5.0d0*rz/r2 - 10.0d0*rz/r2 + 35.0d0*z2*rz/(r2*r2))

        ! --- d(a_J2_z)/drx  (row=6, col=1)
        ! a_J2_z = C0/r5 * z * (3 - 5*z^2/r^2)
        ! d/dx[C0*z*(3-5*z2/r2)/r5] = C0*z*d/dx[(3-5*z2/r2)/r5]
        !   = C0*z*[-5*x*(3-5*z2/r2)/r7 + 10*z2*x/r9]
        !   = C0*z/r5*[-5*x*(3-5*z2/r2)/r2 + 10*z2*x/r4]
        !   = C0*z*x/r5*[-15/r2 + 25*z2/r4 + 10*z2/r4]
        !   = C0*z*x/r5*[-15/r2 + 35*z2/r4]
        row = 6; col = 1
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * rz * rx * (-15.0d0/r2 + 35.0d0*z2/(r2*r2))

        ! --- d(a_J2_z)/dry  (row=6, col=2) — same form as drx with ry
        row = 6; col = 2
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * rz * ry * (-15.0d0/r2 + 35.0d0*z2/(r2*r2))

        ! --- d(a_J2_z)/drz  (row=6, col=3)
        ! d/dz[C0*z*(3-5*z2/r2)/r5]
        !   = C0/r5*[(3-5*z2/r2) + z*d/dz(3-5*z2/r2) + z*(3-5*z2/r2)*(-5*z/r2)]
        ! d/dz(3-5*z2/r2) = -10*z/r2 + 10*z3/r4 = -10*z*(r2-z2)/r4
        ! Wait, let me be more careful:
        ! d/dz(-5*z2/r2) = -5*(2*z*r2 - z2*2*z)/r4 = -5*(2*z*r2 - 2*z3)/r4 = -10*z*(r2-z2)/r4
        ! Full: C0*{(3-5*z2/r2)/r5 + z*[-10*z*(r2-z2)/(r4*r5)] + z*(3-5*z2/r2)*(-5*z/r7)}
        !     = C0/r5*{(3-5*z2/r2) - 10*z2*(r2-z2)/r4 - 5*z2*(3-5*z2/r2)/r2}
        !     = C0/r5*{3 - 5*z2/r2 - 10*z2/r2 + 10*z4/r4 - 15*z2/r2 + 25*z4/r4}
        !     = C0/r5*{3 - 30*z2/r2 + 35*z4/r4}
        row = 6; col = 3
        F((row-1)*6+col) = F((row-1)*6+col) + &
            C0/r5 * (3.0d0 - 30.0d0*z2/r2 + 35.0d0*z2*z2/(r2*r2))

    end subroutine jacobian_j2

    ! =========================================================================
    ! Jacobian for CR3BP (model_id = 3)
    !
    ! State: [x, y, z, vx, vy, vz], n=6 in rotating frame
    ! Params: [mu_ratio], np >= 1
    !
    ! The Jacobian has structure:
    !   [ 0    I    ]
    !   [ Uxx  Cori ]
    !
    ! where Uxx is the gravity+centrifugal gradient and Cori is the Coriolis
    ! matrix [[0, 2, 0], [-2, 0, 0], [0, 0, 0]].
    ! =========================================================================
    subroutine jacobian_cr3bp(n, x, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, mu1, px, py, pz
        real(c_double) :: d1x, d2x
        real(c_double) :: r1_2, r2_2, r1, r2, r1_3, r2_3, r1_5, r2_5
        integer :: i

        if (n /= 6 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        mu  = params(1)
        mu1 = 1.0d0 - mu

        px = x(1); py = x(2); pz = x(3)

        d1x = px + mu       ! x - (-mu)
        d2x = px - mu1      ! x - (1-mu)

        r1_2 = d1x*d1x + py*py + pz*pz
        r2_2 = d2x*d2x + py*py + pz*pz

        if (r1_2 < 1.0d-20 .or. r2_2 < 1.0d-20) then
            info = 1
            F(1:n*n) = 0.0d0
            return
        end if

        r1 = sqrt(r1_2)
        r2 = sqrt(r2_2)
        r1_3 = r1_2 * r1
        r2_3 = r2_2 * r2
        r1_5 = r1_3 * r1_2
        r2_5 = r2_3 * r2_2

        F(1:n*n) = 0.0d0

        ! Top-right 3x3 identity: position dot depends on velocity
        do i = 1, 3
            F((i-1)*6 + 3 + i) = 1.0d0
        end do

        ! Bottom-left 3x3: partial of acceleration w.r.t. position
        ! Uxx = d^2 Omega / dxi dxj where Omega is the effective potential
        ! The pseudo-potential: Omega = 0.5*(x^2+y^2) + (1-mu)/r1 + mu/r2
        !
        ! Omega_xx = 1 - (1-mu)/r1^3 + 3*(1-mu)*d1x^2/r1^5
        !              - mu/r2^3     + 3*mu*d2x^2/r2^5
        ! Omega_yy = 1 - (1-mu)/r1^3 + 3*(1-mu)*y^2/r1^5
        !              - mu/r2^3     + 3*mu*y^2/r2^5
        ! Omega_zz = -(1-mu)/r1^3 + 3*(1-mu)*z^2/r1^5
        !            - mu/r2^3    + 3*mu*z^2/r2^5
        ! Omega_xy = 3*(1-mu)*d1x*y/r1^5 + 3*mu*d2x*y/r2^5
        ! Omega_xz = 3*(1-mu)*d1x*z/r1^5 + 3*mu*d2x*z/r2^5
        ! Omega_yz = 3*(1-mu)*y*z/r1^5   + 3*mu*y*z/r2^5

        ! Row 4, Col 1: d(vx_dot)/dx = Omega_xx
        F(3*6+1) = 1.0d0 - mu1/r1_3 + 3.0d0*mu1*d1x*d1x/r1_5 &
                          - mu/r2_3  + 3.0d0*mu*d2x*d2x/r2_5

        ! Row 4, Col 2: d(vx_dot)/dy = Omega_xy
        F(3*6+2) = 3.0d0*mu1*d1x*py/r1_5 + 3.0d0*mu*d2x*py/r2_5

        ! Row 4, Col 3: d(vx_dot)/dz = Omega_xz
        F(3*6+3) = 3.0d0*mu1*d1x*pz/r1_5 + 3.0d0*mu*d2x*pz/r2_5

        ! Row 5, Col 1: d(vy_dot)/dx = Omega_xy (symmetric)
        F(4*6+1) = F(3*6+2)

        ! Row 5, Col 2: d(vy_dot)/dy = Omega_yy
        F(4*6+2) = 1.0d0 - mu1/r1_3 + 3.0d0*mu1*py*py/r1_5 &
                          - mu/r2_3  + 3.0d0*mu*py*py/r2_5

        ! Row 5, Col 3: d(vy_dot)/dz = Omega_yz
        F(4*6+3) = 3.0d0*mu1*py*pz/r1_5 + 3.0d0*mu*py*pz/r2_5

        ! Row 6, Col 1: d(vz_dot)/dx = Omega_xz (symmetric)
        F(5*6+1) = F(3*6+3)

        ! Row 6, Col 2: d(vz_dot)/dy = Omega_yz (symmetric)
        F(5*6+2) = F(4*6+3)

        ! Row 6, Col 3: d(vz_dot)/dz = Omega_zz
        F(5*6+3) = -mu1/r1_3 + 3.0d0*mu1*pz*pz/r1_5 &
                   -mu/r2_3  + 3.0d0*mu*pz*pz/r2_5

        ! Bottom-right 3x3: Coriolis terms
        ! d(vx_dot)/dvy = 2
        F(3*6+5) = 2.0d0
        ! d(vy_dot)/dvx = -2
        F(4*6+4) = -2.0d0

    end subroutine jacobian_cr3bp

    ! =========================================================================
    ! Jacobian for Atmospheric Drag (model_id = 4)
    !
    ! State: [rx, ry, rz, vx, vy, vz, beta], n=7
    ! Params: [mu, rho0, h_scale, R_body], np >= 4
    !
    ! Jacobian is 7x7:
    !   Rows 1-3: dr/dt = v → same identity block as Kepler (row i, col 3+i)
    !   Rows 4-6: dv/dt = gravity + drag → partials w.r.t. r, v, and beta
    !   Row 7: dbeta/dt = 0 → all zeros
    ! =========================================================================
    subroutine jacobian_drag(n, x, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mu, rho0, h_scale, R_body, beta
        real(c_double) :: rx, ry, rz, vx, vy, vz
        real(c_double) :: r2, r, r5, v2, v_mag
        real(c_double) :: rho, drho_dr_over_r, half_rho_beta
        real(c_double) :: rr(3), vv(3)
        integer :: i, j, row, col

        if (n /= 7 .or. np < 4) then
            info = 3
            return
        end if

        info = 0
        mu      = params(1)
        rho0    = params(2)
        h_scale = params(3)
        R_body  = params(4)

        rx = x(1); ry = x(2); rz = x(3)
        vx = x(4); vy = x(5); vz = x(6)
        beta = x(7)
        rr = [rx, ry, rz]
        vv = [vx, vy, vz]

        r2 = rx*rx + ry*ry + rz*rz
        r  = sqrt(r2)

        if (r < 1.0d-10) then
            info = 1
            F(1:n*n) = 0.0d0
            return
        end if

        r5 = r2 * r2 * r
        v2 = vx*vx + vy*vy + vz*vz
        v_mag = sqrt(v2)

        rho = rho0 * exp(-(r - R_body) / h_scale)

        F(1:n*n) = 0.0d0

        ! Top-right 3x3 identity: dr_dot/dv = I
        do i = 1, 3
            F((i-1)*7 + 3 + i) = 1.0d0
        end do

        ! Gravity gradient (same as Kepler, but 7x7 indexing)
        ! d(a_grav_i)/d(r_j) = -mu*(delta_ij*r2 - 3*ri*rj)/r5
        do i = 1, 3
            do j = 1, 3
                row = 3 + i  ! rows 4,5,6
                col = j      ! cols 1,2,3
                if (i == j) then
                    F((row-1)*7 + col) = -mu * (r2 - 3.0d0*rr(i)*rr(j)) / r5
                else
                    F((row-1)*7 + col) = 3.0d0 * mu * rr(i) * rr(j) / r5
                end if
            end do
        end do

        ! Drag partials
        ! a_drag_i = -0.5 * rho * beta * |v| * v_i
        !
        ! d(a_drag_i)/d(r_j):
        !   drho/dr_j = rho * (-1/h_scale) * r_j/r
        !   d(a_drag_i)/d(r_j) = -0.5 * drho/dr_j * beta * |v| * v_i
        !                       = 0.5 * rho * beta * |v| * v_i / (h_scale * r) * r_j
        if (v_mag > 1.0d-30) then
            half_rho_beta = 0.5d0 * rho * beta
            drho_dr_over_r = -rho / (h_scale * r)  ! drho/d|r| * 1/r = (d/dr_j factor)

            ! d(a_drag_i)/d(r_j) = -0.5*beta*|v|*v_i * drho/d|r| * r_j/r
            !                     = half_rho_beta * |v| * v_i * r_j / (h_scale * r)
            ! (because drho/dr_j = rho*(-1/h_scale)*r_j/r, and a_drag has -0.5*rho*beta
            !  so the signs: -0.5*(-rho/h_scale*rj/r)*beta*|v|*vi = 0.5*rho*beta*|v|*vi*rj/(h_s*r))
            do i = 1, 3
                do j = 1, 3
                    row = 3 + i
                    col = j
                    F((row-1)*7+col) = F((row-1)*7+col) + &
                        half_rho_beta * v_mag * vv(i) * rr(j) / (h_scale * r)
                end do
            end do

            ! d(a_drag_i)/d(v_j):
            !   a_drag_i = -0.5*rho*beta*|v|*v_i
            !   d/dv_j(-0.5*rho*beta*|v|*v_i) = -0.5*rho*beta*(v_j*v_i/|v| + |v|*delta_ij)
            !                                  = -half_rho_beta*(v_i*v_j/|v| + |v|*delta_ij)
            do i = 1, 3
                do j = 1, 3
                    row = 3 + i
                    col = 3 + j
                    F((row-1)*7+col) = -half_rho_beta * (vv(i)*vv(j)/v_mag + v_mag * merge(1.0d0, 0.0d0, i==j))
                end do
            end do

            ! d(a_drag_i)/d(beta):
            !   a_drag_i = -0.5*rho*beta*|v|*v_i
            !   d/d(beta) = -0.5*rho*|v|*v_i
            do i = 1, 3
                row = 3 + i
                col = 7
                F((row-1)*7+col) = -0.5d0 * rho * v_mag * vv(i)
            end do
        end if

        ! Row 7: dbeta/dt = 0 → all zeros (already set)

    end subroutine jacobian_drag

end subroutine fa_dynamics_jacobian
