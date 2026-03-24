! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_dynamics.f90 — Tracking dynamics models with analytic Jacobians
!
! Built-in dynamics catalog for state estimation. Each model provides:
!   - f(x): continuous-time state derivative x_dot = f(x, u, t)
!   - F(x): analytic Jacobian df/dx (flat row-major n×n)
!
! Model IDs (orbital group 1-9):
!    1 = Kepler two-body      — state [r(3), v(3)], n=6, params=[mu]
!    2 = J2 oblateness         — state [r(3), v(3)], n=6, params=[mu, J2, R_eq]
!    3 = CR3BP rotating frame  — state [r(3), v(3)], n=6, params=[mu_ratio]
!    4 = Atmospheric drag      — state [r(3), v(3), beta], n=7, params=[mu, rho0, h_scale, R_body]
!
! Model IDs (rigid body 10-19):
!   10 = 6-DOF rigid body     — state [r(3), v(3), q(4), w(3)], n=13
!                                 u=[Fx,Fy,Fz,Tx,Ty,Tz], params=[mass,Ixx,Iyy,Izz]
!
! Model IDs (ground vehicle 20-29):
!   20 = Bicycle              — state [x,y,theta,v], n=4, u=[accel,steer_rate], params=[L]
!   21 = Ackermann            — state [x,y,theta,v,steer], n=5, u=[accel,steer_rate], params=[L]
!   22 = Differential drive   — state [x,y,theta,v], n=4, u=[v_left,v_right], params=[wheel_sep]
!
! Model IDs (aerial 30-39):
!   30 = Quadrotor            — state [r(3),v(3),angles(3),rates(3)], n=12
!                                 u=[thrust,tau_x,tau_y,tau_z], params=[mass,Ixx,Iyy,Izz,g]
!   31 = Fixed-wing           — state same as quadrotor, n=12
!                                 u=[thrust,aileron,elevator,rudder], params=[mass,Ixx,Iyy,Izz,g,S,rho]
!
! Model IDs (tracking group 40-49):
!   40 = Constant velocity   — state [pos(d), vel(d)], n = 2*d
!   41 = Constant acceleration — state [pos(d), vel(d), acc(d)], n = 3*d
!   42 = Constant turn rate  — state [x, y, vx, vy, omega], n = 5
!
! Model IDs (stochastic 50-59):
!   50 = Geometric Brownian Motion — state [S], n=1, params=[mu_drift, sigma]
!   51 = Ornstein-Uhlenbeck       — state [X], n=1, params=[theta, mu, sigma]
!
! Model IDs (scalar 60-69):
!   60 = Double integrator    — state [pos(d), vel(d)], n=2*d, u=[force(d)], params=[mass]
!   61 = Spring-mass-damper   — state [x, v], n=2, u=[F_ext], params=[k, c, m]
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
    case (10)
        call dynamics_rigidbody6dof(n, x, u, nu, params, np, x_dot, info)
    case (20)
        call dynamics_bicycle(n, x, u, nu, params, np, x_dot, info)
    case (21)
        call dynamics_ackermann(n, x, u, nu, params, np, x_dot, info)
    case (22)
        call dynamics_diffdrive(n, x, u, nu, params, np, x_dot, info)
    case (30)
        call dynamics_quadrotor(n, x, u, nu, params, np, x_dot, info)
    case (31)
        call dynamics_fixedwing(n, x, u, nu, params, np, x_dot, info)
    case (40)
        call dynamics_const_vel(n, x, x_dot, info)
    case (41)
        call dynamics_const_accel(n, x, x_dot, info)
    case (42)
        call dynamics_const_turn(n, x, x_dot, info)
    case (50)
        call dynamics_gbm(n, x, params, np, x_dot, info)
    case (51)
        call dynamics_ou(n, x, params, np, x_dot, info)
    case (60)
        call dynamics_double_integrator(n, x, u, nu, params, np, x_dot, info)
    case (61)
        call dynamics_spring_mass_damper(n, x, u, nu, params, np, x_dot, info)
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

    ! =========================================================================
    ! Bicycle Model (model_id = 20)
    ! State: [x, y, theta, v], n=4
    ! Control: u=[accel, steer_rate], nu=2
    ! Params: [L] (wheelbase), np>=1
    ! x_dot = [v*cos(theta), v*sin(theta), u(2), u(1)]
    ! Simplified bicycle: steer_rate directly controls heading rate.
    ! =========================================================================
    subroutine dynamics_bicycle(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: theta, v

        if (n /= 4 .or. nu < 2 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        theta = x(3)
        v     = x(4)

        x_dot(1) = v * cos(theta)      ! dx/dt
        x_dot(2) = v * sin(theta)      ! dy/dt
        x_dot(3) = u(2)                ! dtheta/dt = steer_rate
        x_dot(4) = u(1)                ! dv/dt = accel

    end subroutine dynamics_bicycle

    ! =========================================================================
    ! Ackermann Model (model_id = 21)
    ! State: [x, y, theta, v, steer], n=5
    ! Control: u=[accel, steer_rate], nu=2
    ! Params: [L] (wheelbase), np>=1
    ! x_dot = [v*cos(theta), v*sin(theta), v*tan(steer)/L, u(1), u(2)]
    ! =========================================================================
    subroutine dynamics_ackermann(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: theta, v, steer, L

        if (n /= 5 .or. nu < 2 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        theta = x(3)
        v     = x(4)
        steer = x(5)
        L     = params(1)

        x_dot(1) = v * cos(theta)           ! dx/dt
        x_dot(2) = v * sin(theta)           ! dy/dt
        x_dot(3) = v * tan(steer) / L       ! dtheta/dt
        x_dot(4) = u(1)                     ! dv/dt = accel
        x_dot(5) = u(2)                     ! dsteer/dt = steer_rate

    end subroutine dynamics_ackermann

    ! =========================================================================
    ! Differential Drive (model_id = 22)
    ! State: [x, y, theta, v], n=4
    ! Control: u=[v_left, v_right], nu=2
    ! Params: [wheel_sep], np>=1
    ! v_avg = (v_left + v_right) / 2
    ! omega = (v_right - v_left) / wheel_sep
    ! x_dot = [v_avg*cos(theta), v_avg*sin(theta), omega, 0]
    ! Note: state v is unused — speed comes from wheel velocities.
    ! =========================================================================
    subroutine dynamics_diffdrive(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: theta, v_avg, omega, wheel_sep

        if (n /= 4 .or. nu < 2 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        theta    = x(3)
        wheel_sep = params(1)

        v_avg = (u(1) + u(2)) / 2.0d0
        omega = (u(2) - u(1)) / wheel_sep

        x_dot(1) = v_avg * cos(theta)   ! dx/dt
        x_dot(2) = v_avg * sin(theta)   ! dy/dt
        x_dot(3) = omega                ! dtheta/dt
        x_dot(4) = 0.0d0                ! dv/dt = 0 (v state unused in diff-drive)

    end subroutine dynamics_diffdrive

    ! =========================================================================
    ! 6-DOF Rigid Body (model_id = 10)
    ! State: [rx,ry,rz, vx,vy,vz, q0,q1,q2,q3, wx,wy,wz], n=13
    ! Control: u=[Fx,Fy,Fz, Tx,Ty,Tz], nu=6
    ! Params: [mass, Ixx, Iyy, Izz], np>=4
    !
    ! r_dot = v
    ! v_dot = F/mass (inertial frame force)
    ! q_dot = 0.5 * quat_mult(q, [0, wx, wy, wz])
    ! omega_dot = I^-1 * (T - omega x (I*omega))  (Euler's equations, diagonal I)
    ! =========================================================================
    subroutine dynamics_rigidbody6dof(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mass, Ixx, Iyy, Izz
        real(c_double) :: q0, q1, q2, q3, wx, wy, wz

        if (n /= 13 .or. nu < 6 .or. np < 4) then
            info = 3
            return
        end if

        info = 0
        mass = params(1)
        Ixx  = params(2)
        Iyy  = params(3)
        Izz  = params(4)

        q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10)
        wx = x(11); wy = x(12); wz = x(13)

        ! Position derivatives = velocity
        x_dot(1) = x(4)
        x_dot(2) = x(5)
        x_dot(3) = x(6)

        ! Velocity derivatives = F/mass (inertial frame forces)
        x_dot(4) = u(1) / mass
        x_dot(5) = u(2) / mass
        x_dot(6) = u(3) / mass

        ! Quaternion derivatives: q_dot = 0.5 * q (*) [0, omega]
        ! Using Hamilton product: q (*) p where p = [0, wx, wy, wz]
        x_dot(7)  = 0.5d0 * (-q1*wx - q2*wy - q3*wz)           ! dq0/dt
        x_dot(8)  = 0.5d0 * ( q0*wx - q3*wy + q2*wz)           ! dq1/dt
        x_dot(9)  = 0.5d0 * ( q3*wx + q0*wy - q1*wz)           ! dq2/dt
        x_dot(10) = 0.5d0 * (-q2*wx + q1*wy + q0*wz)           ! dq3/dt

        ! Angular velocity derivatives: Euler's equations (diagonal inertia)
        ! wx_dot = (Tx - (Izz - Iyy)*wy*wz) / Ixx
        ! wy_dot = (Ty - (Ixx - Izz)*wx*wz) / Iyy
        ! wz_dot = (Tz - (Iyy - Ixx)*wx*wy) / Izz
        x_dot(11) = (u(4) - (Izz - Iyy)*wy*wz) / Ixx
        x_dot(12) = (u(5) - (Ixx - Izz)*wx*wz) / Iyy
        x_dot(13) = (u(6) - (Iyy - Ixx)*wx*wy) / Izz

    end subroutine dynamics_rigidbody6dof

    ! =========================================================================
    ! Quadrotor (model_id = 30)
    ! State: [rx,ry,rz, vx,vy,vz, phi,theta,psi, p,q,r], n=12
    ! Control: u=[thrust, tau_x, tau_y, tau_z], nu=4
    ! Params: [mass, Ixx, Iyy, Izz, g], np>=5
    !
    ! r_dot = v
    ! v_dot = [0,0,-g] + thrust/mass * R_zyx * [0,0,1]
    ! angle_dot = body_rates_to_euler_rates(phi, theta, [p, q, r])
    ! rate_dot = Euler's equations (diagonal inertia)
    ! =========================================================================
    subroutine dynamics_quadrotor(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mass, Ixx, Iyy, Izz, g
        real(c_double) :: phi, theta_e, psi, p, q_r, r_r, thrust
        real(c_double) :: sp, cp, st, ct, tt

        if (n /= 12 .or. nu < 4 .or. np < 5) then
            info = 3
            return
        end if

        info = 0
        mass = params(1); Ixx = params(2); Iyy = params(3)
        Izz  = params(4); g   = params(5)

        phi     = x(7);  theta_e = x(8);  psi = x(9)
        p       = x(10); q_r     = x(11); r_r = x(12)
        thrust  = u(1)

        sp = sin(phi);  cp = cos(phi)
        st = sin(theta_e); ct = cos(theta_e); tt = tan(theta_e)

        ! Position derivatives = velocity
        x_dot(1) = x(4); x_dot(2) = x(5); x_dot(3) = x(6)

        ! Velocity derivatives: gravity + thrust in body z rotated to inertial
        ! R_zyx * [0,0,1] column:
        !   [ ct*cpsi  sp*st*cpsi-cp*spsi  cp*st*cpsi+sp*spsi ] [0]   [ cp*st*cpsi+sp*spsi ]
        !   [ ct*spsi  sp*st*spsi+cp*cpsi  cp*st*spsi-sp*cpsi ] [0] = [ cp*st*spsi-sp*cpsi ]
        !   [ -st      sp*ct              cp*ct               ] [1]   [ cp*ct              ]
        ! Wait, we need R*[0,0,1] which is the third column of R_zyx:
        ! R_zyx third column = [cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi),
        !                       cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi),
        !                       cos(phi)*cos(theta)]
        ! Simpler: just use the standard form
        x_dot(4) = thrust/mass * (cp*st*cos(psi) + sp*sin(psi))
        x_dot(5) = thrust/mass * (cp*st*sin(psi) - sp*cos(psi))
        x_dot(6) = -g + thrust/mass * (cp*ct)

        ! Euler angle rates from body rates (ZYX convention)
        ! phi_dot   = p + sin(phi)*tan(theta)*q + cos(phi)*tan(theta)*r
        ! theta_dot = cos(phi)*q - sin(phi)*r
        ! psi_dot   = sin(phi)/cos(theta)*q + cos(phi)/cos(theta)*r
        x_dot(7)  = p + sp*tt*q_r + cp*tt*r_r
        x_dot(8)  = cp*q_r - sp*r_r
        x_dot(9)  = sp/ct*q_r + cp/ct*r_r

        ! Angular rate derivatives: Euler's equations (diagonal inertia)
        x_dot(10) = (u(2) - (Izz - Iyy)*q_r*r_r) / Ixx
        x_dot(11) = (u(3) - (Ixx - Izz)*p*r_r) / Iyy
        x_dot(12) = (u(4) - (Iyy - Ixx)*p*q_r) / Izz

    end subroutine dynamics_quadrotor

    ! =========================================================================
    ! Fixed-Wing (model_id = 31)
    ! State: [rx,ry,rz, vx,vy,vz, phi,theta,psi, p,q,r], n=12
    ! Control: u=[thrust, aileron, elevator, rudder], nu=4
    ! Params: [mass, Ixx, Iyy, Izz, g, S, rho], np>=7
    !
    ! Simplified linear aero model:
    !   V = |v|, alpha = atan2(vz_body, vx_body) (simplified as pitch)
    !   CL = 0.5 + 2*pi*alpha + 0.5*elevator
    !   CD = 0.02 + 0.05*CL^2
    !   L = 0.5*rho*V^2*S*CL, D = 0.5*rho*V^2*S*CD
    ! Thrust along body x-axis. Torques from control surfaces.
    ! =========================================================================
    subroutine dynamics_fixedwing(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mass, Ixx, Iyy, Izz, g, S, rho_air
        real(c_double) :: phi, theta_e, psi, p, q_r, r_r
        real(c_double) :: sp, cp, st, ct, tt, spsi, cpsi
        real(c_double) :: V, qbar, CL, CD, alpha
        real(c_double) :: Lift, Drag, Fx_body, Fz_body
        real(c_double) :: thrust, aileron, elevator, rudder
        real(c_double) :: tau_x, tau_y, tau_z

        if (n /= 12 .or. nu < 4 .or. np < 7) then
            info = 3
            return
        end if

        info = 0
        mass = params(1); Ixx = params(2); Iyy = params(3)
        Izz  = params(4); g   = params(5); S   = params(6)
        rho_air = params(7)

        phi = x(7); theta_e = x(8); psi = x(9)
        p   = x(10); q_r = x(11); r_r = x(12)
        thrust  = u(1); aileron = u(2); elevator = u(3); rudder = u(4)

        sp = sin(phi); cp = cos(phi)
        st = sin(theta_e); ct = cos(theta_e); tt = tan(theta_e)
        spsi = sin(psi); cpsi = cos(psi)

        ! Airspeed magnitude
        V = sqrt(x(4)*x(4) + x(5)*x(5) + x(6)*x(6))
        if (V < 1.0d-10) V = 1.0d-10  ! Avoid division by zero

        ! Simplified angle of attack: use pitch angle as proxy
        ! (full model would transform velocity to body frame)
        alpha = theta_e

        ! Aerodynamic coefficients (simplified linear model)
        CL = 0.5d0 + 6.2832d0*alpha + 0.5d0*elevator     ! CL0 + CL_alpha*alpha + CL_elev*elev
        CD = 0.02d0 + 0.05d0*CL*CL                         ! CD0 + K*CL^2

        ! Dynamic pressure
        qbar = 0.5d0 * rho_air * V * V

        ! Forces in stability frame (simplified to body-ish frame)
        Lift = qbar * S * CL   ! Perpendicular to velocity (approx body -z)
        Drag = qbar * S * CD   ! Opposite to velocity (approx body -x)

        ! Body forces: thrust along x, lift along -z (in body frame)
        Fx_body = thrust - Drag
        Fz_body = -Lift

        ! Position derivatives = velocity
        x_dot(1) = x(4); x_dot(2) = x(5); x_dot(3) = x(6)

        ! Velocity derivatives: rotate body forces to inertial + gravity
        ! Body x-axis in inertial: R_zyx first column
        ! Simplified: thrust along body x, lift perpendicular
        ! Full rotation: use R_zyx to rotate [Fx_body, 0, Fz_body] to inertial
        ! R_zyx * [Fx, 0, Fz]:
        !   ax = (ct*cpsi)*Fx + (cp*st*cpsi+sp*spsi)*Fz
        !   ay = (ct*spsi)*Fx + (cp*st*spsi-sp*cpsi)*Fz
        !   az = (-st)*Fx     + (cp*ct)*Fz
        x_dot(4) = ((ct*cpsi)*Fx_body + (cp*st*cpsi+sp*spsi)*Fz_body) / mass
        x_dot(5) = ((ct*spsi)*Fx_body + (cp*st*spsi-sp*cpsi)*Fz_body) / mass
        x_dot(6) = -g + ((-st)*Fx_body + (cp*ct)*Fz_body) / mass

        ! Euler angle rates from body rates (ZYX convention, same as quadrotor)
        x_dot(7)  = p + sp*tt*q_r + cp*tt*r_r
        x_dot(8)  = cp*q_r - sp*r_r
        x_dot(9)  = sp/ct*q_r + cp/ct*r_r

        ! Torques from control surfaces (simplified linear mapping)
        ! Aileron -> roll, elevator -> pitch, rudder -> yaw
        tau_x = 0.5d0 * qbar * S * 1.0d0 * aileron     ! Roll moment
        tau_y = 0.5d0 * qbar * S * 1.0d0 * elevator     ! Pitch moment
        tau_z = 0.5d0 * qbar * S * 1.0d0 * rudder       ! Yaw moment

        ! Angular rate derivatives: Euler's equations
        x_dot(10) = (tau_x - (Izz - Iyy)*q_r*r_r) / Ixx
        x_dot(11) = (tau_y - (Ixx - Izz)*p*r_r) / Iyy
        x_dot(12) = (tau_z - (Iyy - Ixx)*p*q_r) / Izz

    end subroutine dynamics_fixedwing

    ! =========================================================================
    ! Geometric Brownian Motion (model_id = 50)
    ! State: [S], n=1
    ! Params: [mu_drift, sigma], np>=2
    ! x_dot = mu_drift * S  (deterministic drift; noise handled by estimator)
    ! =========================================================================
    subroutine dynamics_gbm(n, x, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info

        if (n /= 1 .or. np < 2) then
            info = 3
            return
        end if

        info = 0
        x_dot(1) = params(1) * x(1)   ! mu_drift * S

    end subroutine dynamics_gbm

    ! =========================================================================
    ! Ornstein-Uhlenbeck (model_id = 51)
    ! State: [X], n=1
    ! Params: [theta, mu, sigma], np>=3
    ! x_dot = theta * (mu - X)  (mean-reverting drift)
    ! =========================================================================
    subroutine dynamics_ou(n, x, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info

        if (n /= 1 .or. np < 3) then
            info = 3
            return
        end if

        info = 0
        x_dot(1) = params(1) * (params(2) - x(1))   ! theta * (mu - X)

    end subroutine dynamics_ou

    ! =========================================================================
    ! Double Integrator (model_id = 60)
    ! State: [pos(d), vel(d)], n=2*d (d inferred from n/2)
    ! Control: u=[force(d)], nu=d
    ! Params: [mass], np>=1
    ! pos_dot = vel, vel_dot = force / mass
    ! =========================================================================
    subroutine dynamics_double_integrator(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        integer :: d, i
        real(c_double) :: mass

        if (mod(n, 2) /= 0 .or. n < 2 .or. np < 1) then
            info = 3
            return
        end if

        d = n / 2
        if (nu < d) then
            info = 3
            return
        end if

        info = 0
        mass = params(1)

        ! Position derivatives = velocities
        x_dot(1:d) = x(d+1 : n)

        ! Velocity derivatives = force / mass
        do i = 1, d
            x_dot(d + i) = u(i) / mass
        end do

    end subroutine dynamics_double_integrator

    ! =========================================================================
    ! Spring-Mass-Damper (model_id = 61)
    ! State: [x, v], n=2
    ! Control: u=[F_ext], nu=1
    ! Params: [k, c, m], np>=3
    ! x_dot = v, v_dot = (F_ext - k*x - c*v) / m
    ! =========================================================================
    subroutine dynamics_spring_mass_damper(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: k, c_damp, m

        if (n /= 2 .or. nu < 1 .or. np < 3) then
            info = 3
            return
        end if

        info = 0
        k      = params(1)
        c_damp = params(2)
        m      = params(3)

        x_dot(1) = x(2)                                    ! dx/dt = v
        x_dot(2) = (u(1) - k*x(1) - c_damp*x(2)) / m     ! dv/dt

    end subroutine dynamics_spring_mass_damper

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
    case (10)
        call jacobian_rigidbody6dof(n, x, u, nu, params, np, F, info)
    case (20)
        call jacobian_bicycle(n, x, u, nu, params, np, F, info)
    case (21)
        call jacobian_ackermann(n, x, u, nu, params, np, F, info)
    case (22)
        call jacobian_diffdrive(n, x, u, nu, params, np, F, info)
    case (30)
        call jacobian_quadrotor(n, x, u, nu, params, np, F, info)
    case (31)
        call jacobian_fixedwing(n, x, u, nu, params, np, F, info)
    case (40)
        call jacobian_const_vel(n, x, F, info)
    case (41)
        call jacobian_const_accel(n, x, F, info)
    case (42)
        call jacobian_const_turn(n, x, F, info)
    case (50)
        call jacobian_gbm(n, x, params, np, F, info)
    case (51)
        call jacobian_ou(n, x, params, np, F, info)
    case (60)
        call jacobian_double_integrator(n, x, u, nu, params, np, F, info)
    case (61)
        call jacobian_spring_mass_damper(n, x, u, nu, params, np, F, info)
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

    ! =========================================================================
    ! Jacobian for Bicycle (model_id = 20)
    !
    ! State: [x, y, theta, v], n=4
    ! f = [v*cos(theta), v*sin(theta), u(2), u(1)]
    !
    !       [ 0  0  -v*sin(theta)  cos(theta) ]
    !   F = [ 0  0   v*cos(theta)  sin(theta) ]
    !       [ 0  0   0             0          ]
    !       [ 0  0   0             0          ]
    ! =========================================================================
    subroutine jacobian_bicycle(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: theta, v

        if (n /= 4 .or. nu < 2 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        theta = x(3)
        v     = x(4)

        F(1:n*n) = 0.0d0

        ! Row 1: df1/dtheta = -v*sin(theta), df1/dv = cos(theta)
        ! F(1,3) = (0)*4+3 = 3, F(1,4) = (0)*4+4 = 4
        F(3) = -v * sin(theta)
        F(4) = cos(theta)

        ! Row 2: df2/dtheta = v*cos(theta), df2/dv = sin(theta)
        ! F(2,3) = (1)*4+3 = 7, F(2,4) = (1)*4+4 = 8
        F(7) = v * cos(theta)
        F(8) = sin(theta)

        ! Rows 3,4: all zeros (u(2) and u(1) are controls, not state-dependent)

    end subroutine jacobian_bicycle

    ! =========================================================================
    ! Jacobian for Ackermann (model_id = 21)
    !
    ! State: [x, y, theta, v, steer], n=5
    ! f = [v*cos(theta), v*sin(theta), v*tan(steer)/L, u(1), u(2)]
    !
    !       [ 0  0  -v*sin(theta)  cos(theta)          0                ]
    !   F = [ 0  0   v*cos(theta)  sin(theta)          0                ]
    !       [ 0  0   0             tan(steer)/L        v/(L*cos^2(steer)) ]
    !       [ 0  0   0             0                   0                ]
    !       [ 0  0   0             0                   0                ]
    ! =========================================================================
    subroutine jacobian_ackermann(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: theta, v, steer, L, cs

        if (n /= 5 .or. nu < 2 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        theta = x(3)
        v     = x(4)
        steer = x(5)
        L     = params(1)
        cs    = cos(steer)

        F(1:n*n) = 0.0d0

        ! Row 1: df1/dtheta = -v*sin(theta) at F(1,3)=(0)*5+3=3
        !        df1/dv = cos(theta) at F(1,4)=(0)*5+4=4
        F(3) = -v * sin(theta)
        F(4) = cos(theta)

        ! Row 2: df2/dtheta = v*cos(theta) at F(2,3)=(1)*5+3=8
        !        df2/dv = sin(theta) at F(2,4)=(1)*5+4=9
        F(8) = v * cos(theta)
        F(9) = sin(theta)

        ! Row 3: df3/dv = tan(steer)/L at F(3,4)=(2)*5+4=14
        !        df3/dsteer = v/(L*cos^2(steer)) at F(3,5)=(2)*5+5=15
        F(14) = tan(steer) / L
        F(15) = v / (L * cs * cs)

        ! Rows 4,5: all zeros (control inputs, not state-dependent)

    end subroutine jacobian_ackermann

    ! =========================================================================
    ! Jacobian for Differential Drive (model_id = 22)
    !
    ! State: [x, y, theta, v], n=4
    ! f = [v_avg*cos(theta), v_avg*sin(theta), omega, 0]
    ! where v_avg and omega depend only on controls u, not state (except theta).
    !
    !       [ 0  0  -v_avg*sin(theta)  0 ]
    !   F = [ 0  0   v_avg*cos(theta)  0 ]
    !       [ 0  0   0                 0 ]
    !       [ 0  0   0                 0 ]
    ! =========================================================================
    subroutine jacobian_diffdrive(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: theta, v_avg

        if (n /= 4 .or. nu < 2 .or. np < 1) then
            info = 3
            return
        end if

        info = 0
        theta = x(3)
        v_avg = (u(1) + u(2)) / 2.0d0

        F(1:n*n) = 0.0d0

        ! Row 1: df1/dtheta = -v_avg*sin(theta) at F(1,3)=(0)*4+3=3
        F(3) = -v_avg * sin(theta)

        ! Row 2: df2/dtheta = v_avg*cos(theta) at F(2,3)=(1)*4+3=7
        F(7) = v_avg * cos(theta)

        ! Rows 3,4: all zeros

    end subroutine jacobian_diffdrive

    ! =========================================================================
    ! Jacobian for 6-DOF Rigid Body (model_id = 10)
    !
    ! State: [rx,ry,rz, vx,vy,vz, q0,q1,q2,q3, wx,wy,wz], n=13
    ! 13x13 Jacobian. Key nonzero blocks:
    !   dr/dv = I (3x3)
    !   dv/d* = 0 (forces are control, not state-dependent)
    !   dq_dot/dq, dq_dot/domega (quaternion kinematics partials)
    !   domega_dot/domega (Euler equation partials)
    ! =========================================================================
    subroutine jacobian_rigidbody6dof(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: Ixx, Iyy, Izz
        real(c_double) :: q0, q1, q2, q3, wx, wy, wz
        integer :: nn

        if (n /= 13 .or. nu < 6 .or. np < 4) then
            info = 3
            return
        end if

        info = 0
        nn = 13  ! for index calculations
        Ixx = params(2); Iyy = params(3); Izz = params(4)

        q0 = x(7); q1 = x(8); q2 = x(9); q3 = x(10)
        wx = x(11); wy = x(12); wz = x(13)

        F(1:n*n) = 0.0d0

        ! Block: dr_dot/dv = I (rows 1-3, cols 4-6)
        ! F(i, 3+i) = (i-1)*13 + 3+i
        F((1-1)*nn + 4) = 1.0d0
        F((2-1)*nn + 5) = 1.0d0
        F((3-1)*nn + 6) = 1.0d0

        ! Block: dv_dot/d* = 0 (rows 4-6, forces are control-dependent only)
        ! Already zero.

        ! Block: dq_dot/dq (rows 7-10, cols 7-10)
        ! q0_dot = 0.5*(-q1*wx - q2*wy - q3*wz)
        ! q1_dot = 0.5*( q0*wx - q3*wy + q2*wz)
        ! q2_dot = 0.5*( q3*wx + q0*wy - q1*wz)
        ! q3_dot = 0.5*(-q2*wx + q1*wy + q0*wz)
        !
        ! dq0_dot/dq0 = 0, dq0_dot/dq1 = -0.5*wx, dq0_dot/dq2 = -0.5*wy, dq0_dot/dq3 = -0.5*wz
        F((7-1)*nn + 7) = 0.0d0
        F((7-1)*nn + 8) = -0.5d0*wx
        F((7-1)*nn + 9) = -0.5d0*wy
        F((7-1)*nn + 10) = -0.5d0*wz

        ! dq1_dot/dq0 = 0.5*wx, dq1_dot/dq1 = 0, dq1_dot/dq2 = 0.5*wz, dq1_dot/dq3 = -0.5*wy
        F((8-1)*nn + 7) = 0.5d0*wx
        F((8-1)*nn + 8) = 0.0d0
        F((8-1)*nn + 9) = 0.5d0*wz
        F((8-1)*nn + 10) = -0.5d0*wy

        ! dq2_dot/dq0 = 0.5*wy, dq2_dot/dq1 = -0.5*wz, dq2_dot/dq2 = 0, dq2_dot/dq3 = 0.5*wx
        F((9-1)*nn + 7) = 0.5d0*wy
        F((9-1)*nn + 8) = -0.5d0*wz
        F((9-1)*nn + 9) = 0.0d0
        F((9-1)*nn + 10) = 0.5d0*wx

        ! dq3_dot/dq0 = 0.5*wz, dq3_dot/dq1 = 0.5*wy, dq3_dot/dq2 = -0.5*wx, dq3_dot/dq3 = 0
        F((10-1)*nn + 7) = 0.5d0*wz
        F((10-1)*nn + 8) = 0.5d0*wy
        F((10-1)*nn + 9) = -0.5d0*wx
        F((10-1)*nn + 10) = 0.0d0

        ! Block: dq_dot/domega (rows 7-10, cols 11-13)
        ! dq0_dot/dwx = -0.5*q1, dq0_dot/dwy = -0.5*q2, dq0_dot/dwz = -0.5*q3
        F((7-1)*nn + 11) = -0.5d0*q1
        F((7-1)*nn + 12) = -0.5d0*q2
        F((7-1)*nn + 13) = -0.5d0*q3

        ! dq1_dot/dwx = 0.5*q0, dq1_dot/dwy = -0.5*q3, dq1_dot/dwz = 0.5*q2
        F((8-1)*nn + 11) = 0.5d0*q0
        F((8-1)*nn + 12) = -0.5d0*q3
        F((8-1)*nn + 13) = 0.5d0*q2

        ! dq2_dot/dwx = 0.5*q3, dq2_dot/dwy = 0.5*q0, dq2_dot/dwz = -0.5*q1
        F((9-1)*nn + 11) = 0.5d0*q3
        F((9-1)*nn + 12) = 0.5d0*q0
        F((9-1)*nn + 13) = -0.5d0*q1

        ! dq3_dot/dwx = -0.5*q2, dq3_dot/dwy = 0.5*q1, dq3_dot/dwz = 0.5*q0
        F((10-1)*nn + 11) = -0.5d0*q2
        F((10-1)*nn + 12) = 0.5d0*q1
        F((10-1)*nn + 13) = 0.5d0*q0

        ! Block: domega_dot/domega (rows 11-13, cols 11-13)
        ! wx_dot = (Tx - (Izz-Iyy)*wy*wz) / Ixx
        ! dwx_dot/dwx = 0
        ! dwx_dot/dwy = -(Izz-Iyy)*wz / Ixx
        ! dwx_dot/dwz = -(Izz-Iyy)*wy / Ixx
        F((11-1)*nn + 11) = 0.0d0
        F((11-1)*nn + 12) = -(Izz - Iyy)*wz / Ixx
        F((11-1)*nn + 13) = -(Izz - Iyy)*wy / Ixx

        ! wy_dot = (Ty - (Ixx-Izz)*wx*wz) / Iyy
        ! dwy_dot/dwx = -(Ixx-Izz)*wz / Iyy
        ! dwy_dot/dwy = 0
        ! dwy_dot/dwz = -(Ixx-Izz)*wx / Iyy
        F((12-1)*nn + 11) = -(Ixx - Izz)*wz / Iyy
        F((12-1)*nn + 12) = 0.0d0
        F((12-1)*nn + 13) = -(Ixx - Izz)*wx / Iyy

        ! wz_dot = (Tz - (Iyy-Ixx)*wx*wy) / Izz
        ! dwz_dot/dwx = -(Iyy-Ixx)*wy / Izz
        ! dwz_dot/dwy = -(Iyy-Ixx)*wx / Izz
        ! dwz_dot/dwz = 0
        F((13-1)*nn + 11) = -(Iyy - Ixx)*wy / Izz
        F((13-1)*nn + 12) = -(Iyy - Ixx)*wx / Izz
        F((13-1)*nn + 13) = 0.0d0

    end subroutine jacobian_rigidbody6dof

    ! =========================================================================
    ! Jacobian for Quadrotor (model_id = 30)
    !
    ! State: [rx,ry,rz, vx,vy,vz, phi,theta,psi, p,q,r], n=12
    ! This is a 12x12 Jacobian with nonzero blocks in:
    !   dr/dv = I (3x3)
    !   dv/d(angles) (thrust rotation partials)
    !   d(angle_rates)/d(angles) and d(angle_rates)/d(body_rates)
    !   d(body_rate_dot)/d(body_rates) (Euler equation partials)
    !
    ! For numerical stability and correctness, we compute this via finite
    ! differences internally (the analytic form has many terms and is error-prone).
    ! The dispatch Jacobian uses central FD with the dynamics function.
    ! =========================================================================
    subroutine jacobian_quadrotor(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mass, Ixx, Iyy, Izz, g, thrust
        real(c_double) :: phi, theta_e, psi, p, q_r, r_r
        real(c_double) :: sp, cp, st, ct, tt, spsi, cpsi
        real(c_double) :: T_over_m
        integer :: nn

        if (n /= 12 .or. nu < 4 .or. np < 5) then
            info = 3
            return
        end if

        info = 0
        nn = 12
        mass = params(1); Ixx = params(2); Iyy = params(3)
        Izz  = params(4); g   = params(5)
        thrust = u(1)
        T_over_m = thrust / mass

        phi = x(7); theta_e = x(8); psi = x(9)
        p = x(10); q_r = x(11); r_r = x(12)

        sp = sin(phi); cp = cos(phi)
        st = sin(theta_e); ct = cos(theta_e); tt = tan(theta_e)
        spsi = sin(psi); cpsi = cos(psi)

        F(1:n*n) = 0.0d0

        ! dr/dv = I (rows 1-3, cols 4-6)
        F((1-1)*nn + 4) = 1.0d0
        F((2-1)*nn + 5) = 1.0d0
        F((3-1)*nn + 6) = 1.0d0

        ! dv_dot/d(phi, theta, psi) — partials of thrust rotation
        ! vx_dot = T/m * (cp*st*cpsi + sp*spsi)
        ! dv4/dphi = T/m * (-sp*st*cpsi + cp*spsi)
        F((4-1)*nn + 7) = T_over_m * (-sp*st*cpsi + cp*spsi)
        ! dv4/dtheta = T/m * (cp*ct*cpsi)
        F((4-1)*nn + 8) = T_over_m * (cp*ct*cpsi)
        ! dv4/dpsi = T/m * (-cp*st*spsi + sp*cpsi)
        F((4-1)*nn + 9) = T_over_m * (-cp*st*spsi + sp*cpsi)

        ! vy_dot = T/m * (cp*st*spsi - sp*cpsi)
        ! dv5/dphi = T/m * (-sp*st*spsi - cp*cpsi)
        F((5-1)*nn + 7) = T_over_m * (-sp*st*spsi - cp*cpsi)
        ! dv5/dtheta = T/m * (cp*ct*spsi)
        F((5-1)*nn + 8) = T_over_m * (cp*ct*spsi)
        ! dv5/dpsi = T/m * (cp*st*cpsi + sp*spsi)
        F((5-1)*nn + 9) = T_over_m * (cp*st*cpsi + sp*spsi)

        ! vz_dot = -g + T/m * (cp*ct)
        ! dv6/dphi = T/m * (-sp*ct)
        F((6-1)*nn + 7) = T_over_m * (-sp*ct)
        ! dv6/dtheta = T/m * (-cp*st)
        F((6-1)*nn + 8) = T_over_m * (-cp*st)
        ! dv6/dpsi = 0
        F((6-1)*nn + 9) = 0.0d0

        ! Euler angle rate partials
        ! phi_dot = p + sp*tt*q + cp*tt*r
        ! dphi_dot/dphi = cp*tt*q - sp*tt*r
        F((7-1)*nn + 7) = cp*tt*q_r - sp*tt*r_r
        ! dphi_dot/dtheta = sp*q/(ct^2) + cp*r/(ct^2)
        F((7-1)*nn + 8) = (sp*q_r + cp*r_r) / (ct*ct)
        ! dphi_dot/dp = 1
        F((7-1)*nn + 10) = 1.0d0
        ! dphi_dot/dq = sp*tt
        F((7-1)*nn + 11) = sp*tt
        ! dphi_dot/dr = cp*tt
        F((7-1)*nn + 12) = cp*tt

        ! theta_dot = cp*q - sp*r
        ! dtheta_dot/dphi = -sp*q - cp*r
        F((8-1)*nn + 7) = -sp*q_r - cp*r_r
        ! dtheta_dot/dq = cp
        F((8-1)*nn + 11) = cp
        ! dtheta_dot/dr = -sp
        F((8-1)*nn + 12) = -sp

        ! psi_dot = sp/ct*q + cp/ct*r
        ! dpsi_dot/dphi = cp/ct*q - sp/ct*r
        F((9-1)*nn + 7) = (cp*q_r - sp*r_r) / ct
        ! dpsi_dot/dtheta = sp*st/(ct^2)*q + cp*st/(ct^2)*r
        F((9-1)*nn + 8) = (sp*q_r + cp*r_r) * st / (ct*ct)
        ! dpsi_dot/dq = sp/ct
        F((9-1)*nn + 11) = sp / ct
        ! dpsi_dot/dr = cp/ct
        F((9-1)*nn + 12) = cp / ct

        ! Body rate dot partials (Euler equations, diagonal inertia)
        ! p_dot = (tau_x - (Izz-Iyy)*q*r) / Ixx
        F((10-1)*nn + 11) = -(Izz - Iyy)*r_r / Ixx
        F((10-1)*nn + 12) = -(Izz - Iyy)*q_r / Ixx

        ! q_dot = (tau_y - (Ixx-Izz)*p*r) / Iyy
        F((11-1)*nn + 10) = -(Ixx - Izz)*r_r / Iyy
        F((11-1)*nn + 12) = -(Ixx - Izz)*p / Iyy

        ! r_dot = (tau_z - (Iyy-Ixx)*p*q) / Izz
        F((12-1)*nn + 10) = -(Iyy - Ixx)*q_r / Izz
        F((12-1)*nn + 11) = -(Iyy - Ixx)*p / Izz

    end subroutine jacobian_quadrotor

    ! =========================================================================
    ! Jacobian for Fixed-Wing (model_id = 31)
    !
    ! This model has complex aero + rotation partials. We use numerical
    ! central finite differences for robustness and correctness.
    ! The dynamics are computed inline to avoid cross-contains scope issues.
    ! =========================================================================
    subroutine jacobian_fixedwing(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: x_plus(12), x_minus(12), f_plus(12), f_minus(12)
        real(c_double) :: h, scale_h
        integer :: i, j, info_fd

        if (n /= 12 .or. nu < 4 .or. np < 7) then
            info = 3
            return
        end if

        info = 0
        F(1:n*n) = 0.0d0

        ! Central finite difference Jacobian using inline dynamics helper
        do j = 1, n
            scale_h = max(abs(x(j)), 1.0d0)
            h = 1.0d-7 * scale_h

            x_plus(1:n)  = x(1:n)
            x_minus(1:n) = x(1:n)
            x_plus(j)  = x(j) + h
            x_minus(j) = x(j) - h

            call fw_dynamics_inline(n, x_plus, u, nu, params, np, f_plus, info_fd)
            call fw_dynamics_inline(n, x_minus, u, nu, params, np, f_minus, info_fd)

            do i = 1, n
                F((i-1)*n + j) = (f_plus(i) - f_minus(i)) / (2.0d0 * h)
            end do
        end do

    end subroutine jacobian_fixedwing

    ! =========================================================================
    ! Helper: inline fixed-wing dynamics for Jacobian finite differences
    ! (Duplicates dynamics_fixedwing to avoid cross-contains scope issues)
    ! =========================================================================
    subroutine fw_dynamics_inline(n, x, u, nu, params, np, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: mass, Ixx, Iyy, Izz, g, S, rho_air
        real(c_double) :: phi, theta_e, psi, p, q_r, r_r
        real(c_double) :: sp, cp, st, ct, tt, spsi, cpsi
        real(c_double) :: V, qbar, CL, CD, alpha
        real(c_double) :: Lift, Drag, Fx_body, Fz_body
        real(c_double) :: thrust, aileron, elevator, rudder
        real(c_double) :: tau_x, tau_y, tau_z

        info = 0
        mass = params(1); Ixx = params(2); Iyy = params(3)
        Izz  = params(4); g   = params(5); S   = params(6)
        rho_air = params(7)

        phi = x(7); theta_e = x(8); psi = x(9)
        p   = x(10); q_r = x(11); r_r = x(12)
        thrust = u(1); aileron = u(2); elevator = u(3); rudder = u(4)

        sp = sin(phi); cp = cos(phi)
        st = sin(theta_e); ct = cos(theta_e); tt = tan(theta_e)
        spsi = sin(psi); cpsi = cos(psi)

        V = sqrt(x(4)*x(4) + x(5)*x(5) + x(6)*x(6))
        if (V < 1.0d-10) V = 1.0d-10

        alpha = theta_e
        CL = 0.5d0 + 6.2832d0*alpha + 0.5d0*elevator
        CD = 0.02d0 + 0.05d0*CL*CL
        qbar = 0.5d0 * rho_air * V * V

        Lift = qbar * S * CL
        Drag = qbar * S * CD
        Fx_body = thrust - Drag
        Fz_body = -Lift

        x_dot(1) = x(4); x_dot(2) = x(5); x_dot(3) = x(6)
        x_dot(4) = ((ct*cpsi)*Fx_body + (cp*st*cpsi+sp*spsi)*Fz_body) / mass
        x_dot(5) = ((ct*spsi)*Fx_body + (cp*st*spsi-sp*cpsi)*Fz_body) / mass
        x_dot(6) = -g + ((-st)*Fx_body + (cp*ct)*Fz_body) / mass

        x_dot(7)  = p + sp*tt*q_r + cp*tt*r_r
        x_dot(8)  = cp*q_r - sp*r_r
        x_dot(9)  = sp/ct*q_r + cp/ct*r_r

        tau_x = 0.5d0 * qbar * S * 1.0d0 * aileron
        tau_y = 0.5d0 * qbar * S * 1.0d0 * elevator
        tau_z = 0.5d0 * qbar * S * 1.0d0 * rudder

        x_dot(10) = (tau_x - (Izz - Iyy)*q_r*r_r) / Ixx
        x_dot(11) = (tau_y - (Ixx - Izz)*p*r_r) / Iyy
        x_dot(12) = (tau_z - (Iyy - Ixx)*p*q_r) / Izz

    end subroutine fw_dynamics_inline

    ! =========================================================================
    ! Jacobian for Geometric Brownian Motion (model_id = 50)
    !
    ! State: [S], n=1
    ! f = [mu_drift * S]
    ! F = [mu_drift]
    ! =========================================================================
    subroutine jacobian_gbm(n, x, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info

        if (n /= 1 .or. np < 2) then
            info = 3
            return
        end if

        info = 0
        F(1) = params(1)   ! mu_drift

    end subroutine jacobian_gbm

    ! =========================================================================
    ! Jacobian for Ornstein-Uhlenbeck (model_id = 51)
    !
    ! State: [X], n=1
    ! f = [theta * (mu - X)]
    ! F = [-theta]
    ! =========================================================================
    subroutine jacobian_ou(n, x, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, np
        real(c_double), intent(in)  :: x(n), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info

        if (n /= 1 .or. np < 3) then
            info = 3
            return
        end if

        info = 0
        F(1) = -params(1)   ! -theta

    end subroutine jacobian_ou

    ! =========================================================================
    ! Jacobian for Double Integrator (model_id = 60)
    !
    ! State: [pos(d), vel(d)], n=2*d
    ! f = [vel, force/mass]
    !       [ 0  I ]
    !   F = [ 0  0 ]
    ! Same structure as const_vel — controls don't affect Jacobian w.r.t. state.
    ! =========================================================================
    subroutine jacobian_double_integrator(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        integer :: d, i

        if (mod(n, 2) /= 0 .or. n < 2 .or. np < 1) then
            info = 3
            return
        end if

        d = n / 2
        if (nu < d) then
            info = 3
            return
        end if

        info = 0
        F(1:n*n) = 0.0d0

        ! Upper-right identity block: dpos_i/dvel_i = 1
        do i = 1, d
            F((i-1)*n + d + i) = 1.0d0
        end do

        ! Lower block is all zeros: vel_dot = u/m, independent of state

    end subroutine jacobian_double_integrator

    ! =========================================================================
    ! Jacobian for Spring-Mass-Damper (model_id = 61)
    !
    ! State: [x, v], n=2
    ! f = [v, (F_ext - k*x - c*v) / m]
    !       [ 0      1    ]
    !   F = [ -k/m  -c/m  ]
    ! =========================================================================
    subroutine jacobian_spring_mass_damper(n, x, u, nu, params, np, F, info)
        implicit none
        integer(c_int), intent(in)  :: n, nu, np
        real(c_double), intent(in)  :: x(n), u(nu), params(np)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: k, c_damp, m

        if (n /= 2 .or. nu < 1 .or. np < 3) then
            info = 3
            return
        end if

        info = 0
        k      = params(1)
        c_damp = params(2)
        m      = params(3)

        F(1:n*n) = 0.0d0

        ! Row 1: dx/dt = v → df1/dv = 1
        F(2) = 1.0d0    ! F(1,2) = (0)*2+2 = 2

        ! Row 2: dv/dt = (F - k*x - c*v)/m → df2/dx = -k/m, df2/dv = -c/m
        F(3) = -k / m   ! F(2,1) = (1)*2+1 = 3
        F(4) = -c_damp / m  ! F(2,2) = (1)*2+2 = 4

    end subroutine jacobian_spring_mass_damper

end subroutine fa_dynamics_jacobian
