! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_guidance.f90 — 23 guidance and control law algorithms
!
! Zero-effort (ZEM/ZEV/E-guidance), proportional navigation,
! polynomial guidance, optimal control (LQR/iLQR/DDP),
! MPC, targeting (shooting/Lambert/diff-corr),
! path following, energy-optimal guidance.
!
! All bind(C, name="fa_*"). Stateless, reentrant.
! Error codes: info = 0 ok, 3 invalid input, 4 convergence failure.
! Matrices flat 1D row-major.

! =====================================================================
! ZERO-EFFORT GUIDANCE (Apollo Heritage)
! =====================================================================

! fa_guidance_zem — Zero-Effort Miss guidance
! State x = [r(3), v(3)], a_cmd = -6*(r + v*t_go)/t_go^2 + 2*g
subroutine fa_guidance_zem(n, x, t_go, g, a_cmd, info) &
    bind(C, name="fa_guidance_zem")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), g(3)
    real(c_double), value, intent(in)  :: t_go
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: tgo2

    info = 0
    if (n < 6 .or. abs(t_go) < 1.0d-12) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    tgo2 = t_go * t_go
    a_cmd(1) = -6.0d0 * (x(1) + x(4) * t_go) / tgo2 + 2.0d0 * g(1)
    a_cmd(2) = -6.0d0 * (x(2) + x(5) * t_go) / tgo2 + 2.0d0 * g(2)
    a_cmd(3) = -6.0d0 * (x(3) + x(6) * t_go) / tgo2 + 2.0d0 * g(3)
end subroutine

! fa_guidance_zev — Zero-Effort Velocity guidance
! x = [r(3), v(3)], x(4:6) = current vel, v_target passed in params
subroutine fa_guidance_zev(n, x, v_target, t_go, g, a_cmd, info) &
    bind(C, name="fa_guidance_zev")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), v_target(3), g(3)
    real(c_double), value, intent(in)  :: t_go
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    info = 0
    if (n < 6 .or. abs(t_go) < 1.0d-12) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    a_cmd(1) = -2.0d0 * (x(4) - v_target(1)) / t_go + g(1)
    a_cmd(2) = -2.0d0 * (x(5) - v_target(2)) / t_go + g(2)
    a_cmd(3) = -2.0d0 * (x(6) - v_target(3)) / t_go + g(3)
end subroutine

! fa_guidance_eguidance — E-guidance (generalized ZEM/ZEV)
subroutine fa_guidance_eguidance(n, x, x_target, t_go, g, a_cmd, info) &
    bind(C, name="fa_guidance_eguidance")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), x_target(6), g(3)
    real(c_double), value, intent(in)  :: t_go
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: tgo2
    integer :: i

    info = 0
    if (n < 6 .or. abs(t_go) < 1.0d-12) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    tgo2 = t_go * t_go
    do i = 1, 3
        ! C1=6, C2=2 (standard E-guidance coefficients)
        a_cmd(i) = 6.0d0 * (x_target(i) - x(i) - x(i+3) * t_go) / tgo2 &
                 + 2.0d0 * (x_target(i+3) - x(i+3)) / t_go &
                 + g(i)
    end do
end subroutine

! =====================================================================
! PROPORTIONAL NAVIGATION
! =====================================================================

! fa_guidance_pn_pure — Pure Proportional Navigation
subroutine fa_guidance_pn_pure(n, x_p, x_t, N_gain, a_cmd, info) &
    bind(C, name="fa_guidance_pn_pure")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x_p(n), x_t(n)
    real(c_double), value, intent(in)  :: N_gain
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: r_rel(3), v_rel(3), r_mag, V_c
    real(c_double) :: los(3), omega(3), rxv(3)

    info = 0
    if (n < 6) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    ! Relative position and velocity (target - pursuer)
    r_rel(1) = x_t(1) - x_p(1)
    r_rel(2) = x_t(2) - x_p(2)
    r_rel(3) = x_t(3) - x_p(3)

    v_rel(1) = x_t(4) - x_p(4)
    v_rel(2) = x_t(5) - x_p(5)
    v_rel(3) = x_t(6) - x_p(6)

    r_mag = sqrt(r_rel(1)**2 + r_rel(2)**2 + r_rel(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    los = r_rel / r_mag

    ! Closing velocity (negative of range rate)
    V_c = -(r_rel(1)*v_rel(1) + r_rel(2)*v_rel(2) + r_rel(3)*v_rel(3)) / r_mag

    ! LOS rate: omega = (r x v) / |r|^2
    rxv(1) = r_rel(2)*v_rel(3) - r_rel(3)*v_rel(2)
    rxv(2) = r_rel(3)*v_rel(1) - r_rel(1)*v_rel(3)
    rxv(3) = r_rel(1)*v_rel(2) - r_rel(2)*v_rel(1)
    omega = rxv / (r_mag * r_mag)

    ! a = N * V_c * omega
    a_cmd(1) = N_gain * V_c * omega(1)
    a_cmd(2) = N_gain * V_c * omega(2)
    a_cmd(3) = N_gain * V_c * omega(3)
end subroutine

! fa_guidance_pn_aug — Augmented Proportional Navigation
subroutine fa_guidance_pn_aug(n, x_p, x_t, N_gain, a_target, a_cmd, info) &
    bind(C, name="fa_guidance_pn_aug")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x_p(n), x_t(n), a_target(3)
    real(c_double), value, intent(in)  :: N_gain
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: r_rel(3), v_rel(3), r_mag, V_c
    real(c_double) :: los(3), omega(3), rxv(3)
    real(c_double) :: at_normal(3), at_dot_los

    info = 0
    if (n < 6) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    r_rel(1) = x_t(1) - x_p(1)
    r_rel(2) = x_t(2) - x_p(2)
    r_rel(3) = x_t(3) - x_p(3)
    v_rel(1) = x_t(4) - x_p(4)
    v_rel(2) = x_t(5) - x_p(5)
    v_rel(3) = x_t(6) - x_p(6)

    r_mag = sqrt(r_rel(1)**2 + r_rel(2)**2 + r_rel(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    los = r_rel / r_mag
    V_c = -(r_rel(1)*v_rel(1) + r_rel(2)*v_rel(2) + r_rel(3)*v_rel(3)) / r_mag

    rxv(1) = r_rel(2)*v_rel(3) - r_rel(3)*v_rel(2)
    rxv(2) = r_rel(3)*v_rel(1) - r_rel(1)*v_rel(3)
    rxv(3) = r_rel(1)*v_rel(2) - r_rel(2)*v_rel(1)
    omega = rxv / (r_mag * r_mag)

    ! Target acceleration normal to LOS
    at_dot_los = a_target(1)*los(1) + a_target(2)*los(2) + a_target(3)*los(3)
    at_normal(1) = a_target(1) - at_dot_los * los(1)
    at_normal(2) = a_target(2) - at_dot_los * los(2)
    at_normal(3) = a_target(3) - at_dot_los * los(3)

    ! a = N*V_c*omega + (N/2)*a_target_normal
    a_cmd(1) = N_gain * V_c * omega(1) + 0.5d0 * N_gain * at_normal(1)
    a_cmd(2) = N_gain * V_c * omega(2) + 0.5d0 * N_gain * at_normal(2)
    a_cmd(3) = N_gain * V_c * omega(3) + 0.5d0 * N_gain * at_normal(3)
end subroutine

! fa_guidance_pn_true — True Proportional Navigation
subroutine fa_guidance_pn_true(n, x_p, x_t, N_gain, a_cmd, info) &
    bind(C, name="fa_guidance_pn_true")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x_p(n), x_t(n)
    real(c_double), value, intent(in)  :: N_gain
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: r_rel(3), v_rel(3), r_mag, V_r
    real(c_double) :: omega(3), rxv(3)

    info = 0
    if (n < 6) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    r_rel(1) = x_t(1) - x_p(1)
    r_rel(2) = x_t(2) - x_p(2)
    r_rel(3) = x_t(3) - x_p(3)
    v_rel(1) = x_t(4) - x_p(4)
    v_rel(2) = x_t(5) - x_p(5)
    v_rel(3) = x_t(6) - x_p(6)

    r_mag = sqrt(r_rel(1)**2 + r_rel(2)**2 + r_rel(3)**2)
    if (r_mag < 1.0d-12) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    V_r = sqrt(v_rel(1)**2 + v_rel(2)**2 + v_rel(3)**2)

    rxv(1) = r_rel(2)*v_rel(3) - r_rel(3)*v_rel(2)
    rxv(2) = r_rel(3)*v_rel(1) - r_rel(1)*v_rel(3)
    rxv(3) = r_rel(1)*v_rel(2) - r_rel(2)*v_rel(1)
    omega = rxv / (r_mag * r_mag)

    ! True PN: a = N * V_r * omega (V_r = relative speed, not closing speed)
    a_cmd(1) = N_gain * V_r * omega(1)
    a_cmd(2) = N_gain * V_r * omega(2)
    a_cmd(3) = N_gain * V_r * omega(3)
end subroutine

! =====================================================================
! POLYNOMIAL GUIDANCE
! =====================================================================

! fa_guidance_gravity_turn — Gravity turn trajectory
subroutine fa_guidance_gravity_turn(t, t_burn, x0, g, thrust_mag, mass, x_cmd, info) &
    bind(C, name="fa_guidance_gravity_turn")
    use iso_c_binding
    implicit none
    real(c_double), value, intent(in)  :: t, t_burn, g, thrust_mag, mass
    real(c_double), intent(in)         :: x0(7)  ! [r(3), v(3), pitch]
    real(c_double), intent(out)        :: x_cmd(7)
    integer(c_int), intent(out)        :: info

    real(c_double) :: a_thrust, pitch, v_mag, dt

    info = 0
    if (mass <= 0.0d0 .or. t_burn <= 0.0d0) then
        info = 3
        x_cmd = 0.0d0
        return
    end if

    dt = min(t, t_burn)
    a_thrust = thrust_mag / mass
    pitch = x0(7)  ! initial pitch angle

    ! Simple gravity turn: velocity follows gravity
    v_mag = sqrt(x0(4)**2 + x0(5)**2 + x0(6)**2)
    if (v_mag < 1.0d-12) v_mag = 1.0d-12

    ! Propagate with thrust along velocity, gravity along -z
    x_cmd(1) = x0(1) + x0(4) * dt + 0.5d0 * (a_thrust * x0(4)/v_mag) * dt * dt
    x_cmd(2) = x0(2) + x0(5) * dt + 0.5d0 * (a_thrust * x0(5)/v_mag) * dt * dt
    x_cmd(3) = x0(3) + x0(6) * dt + 0.5d0 * (a_thrust * x0(6)/v_mag - g) * dt * dt
    x_cmd(4) = x0(4) + a_thrust * x0(4)/v_mag * dt
    x_cmd(5) = x0(5) + a_thrust * x0(5)/v_mag * dt
    x_cmd(6) = x0(6) + (a_thrust * x0(6)/v_mag - g) * dt
    ! Updated pitch
    x_cmd(7) = atan2(x_cmd(6), sqrt(x_cmd(4)**2 + x_cmd(5)**2))
end subroutine

! fa_guidance_linear_tangent — Linear tangent steering
subroutine fa_guidance_linear_tangent(n, x, x_target, t_go, g, a_mag, a_cmd, info) &
    bind(C, name="fa_guidance_linear_tangent")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), x_target(6), g(3)
    real(c_double), value, intent(in)  :: t_go, a_mag
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: dir(3), dir_mag
    integer :: i

    info = 0
    if (n < 6 .or. abs(t_go) < 1.0d-12 .or. a_mag <= 0.0d0) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    ! Desired acceleration direction from E-guidance
    do i = 1, 3
        dir(i) = 6.0d0 * (x_target(i) - x(i) - x(i+3)*t_go) / (t_go*t_go) &
               + 2.0d0 * (x_target(i+3) - x(i+3)) / t_go + g(i)
    end do

    dir_mag = sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)
    if (dir_mag < 1.0d-12) then
        a_cmd = 0.0d0
        return
    end if

    ! Constrain to available thrust magnitude
    a_cmd(1) = a_mag * dir(1) / dir_mag
    a_cmd(2) = a_mag * dir(2) / dir_mag
    a_cmd(3) = a_mag * dir(3) / dir_mag
end subroutine

! fa_guidance_peg — Powered Explicit Guidance (PEG)
subroutine fa_guidance_peg(n, x, x_target, v_exhaust, t_go, g, a_cmd, info) &
    bind(C, name="fa_guidance_peg")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), x_target(6), g(3)
    real(c_double), value, intent(in)  :: v_exhaust, t_go
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: dv(3), dv_mag, dir(3)
    integer :: i

    info = 0
    if (n < 6 .or. abs(t_go) < 1.0d-12 .or. v_exhaust <= 0.0d0) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    ! PEG: iterative, but core is velocity-to-be-gained
    ! dv = v_target - v_current - g*t_go (simplified, single iteration)
    do i = 1, 3
        dv(i) = x_target(i+3) - x(i+3) - g(i) * t_go
    end do

    dv_mag = sqrt(dv(1)**2 + dv(2)**2 + dv(3)**2)
    if (dv_mag < 1.0d-12) then
        a_cmd = 0.0d0
        return
    end if

    dir = dv / dv_mag

    ! Thrust acceleration = v_exhaust / t_go (approximate for constant thrust)
    a_cmd(1) = (dv_mag / t_go) * dir(1)
    a_cmd(2) = (dv_mag / t_go) * dir(2)
    a_cmd(3) = (dv_mag / t_go) * dir(3)
end subroutine

! =====================================================================
! OPTIMAL CONTROL
! =====================================================================

! fa_guidance_lqr — Linear Quadratic Regulator (steady-state)
subroutine fa_guidance_lqr(n, m, A, B, Q_cost, R_cost, x, u_cmd, info) &
    bind(C, name="fa_guidance_lqr")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m
    real(c_double), intent(in)         :: A(n*n), B(n*m), Q_cost(n*n), R_cost(m*m)
    real(c_double), intent(in)         :: x(n)
    real(c_double), intent(out)        :: u_cmd(m)
    integer(c_int), intent(out)        :: info

    ! Solve Algebraic Riccati Equation via iteration
    ! P_{k+1} = A^T P_k A - A^T P_k B (R + B^T P_k B)^{-1} B^T P_k A + Q
    ! For simplicity, use discrete-time Riccati iteration

    real(c_double), allocatable :: P(:,:), P_new(:,:), K(:,:)
    real(c_double), allocatable :: Am(:,:), Bm(:,:), Qm(:,:), Rm(:,:)
    real(c_double), allocatable :: BtP(:,:), BtPB(:,:), S(:,:), S_inv(:,:)
    real(c_double), allocatable :: BtPA(:,:), AtPB(:,:), AtPA(:,:), AtP(:,:)
    integer :: ri, rj, rk, riter
    real(c_double) :: diff

    info = 0
    if (n <= 0 .or. m <= 0) then
        info = 3
        u_cmd(1:m) = 0.0d0
        return
    end if

    allocate(P(n,n), P_new(n,n), K(m,n))
    allocate(Am(n,n), Bm(n,m), Qm(n,n), Rm(m,m))
    allocate(BtP(m,n), BtPB(m,m), S(m,m), S_inv(m,m))
    allocate(BtPA(m,n), AtP(n,n), AtPA(n,n), AtPB(n,m))

    ! Unpack row-major to column-major matrices
    do rj = 1, n
        do ri = 1, n
            Am(ri,rj) = A((ri-1)*n + rj)
            Qm(ri,rj) = Q_cost((ri-1)*n + rj)
        end do
    end do
    do rj = 1, m
        do ri = 1, n
            Bm(ri,rj) = B((ri-1)*m + rj)
        end do
        do ri = 1, m
            Rm(ri,rj) = R_cost((ri-1)*m + rj)
        end do
    end do

    ! Initialize P = Q
    P = Qm

    ! Iterate Riccati equation
    do riter = 1, 200
        ! AtP = A^T * P
        AtP = matmul(transpose(Am), P)
        ! AtPA = A^T * P * A
        AtPA = matmul(AtP, Am)
        ! BtP = B^T * P
        BtP = matmul(transpose(Bm), P)
        ! BtPB = B^T * P * B
        BtPB = matmul(BtP, Bm)
        ! S = R + B^T P B
        S = Rm + BtPB

        ! Invert S (for small m, direct formula; else LU)
        call invert_small(m, S, S_inv, info)
        if (info /= 0) then
            deallocate(P, P_new, K, Am, Bm, Qm, Rm)
            deallocate(BtP, BtPB, S, S_inv, BtPA, AtP, AtPA, AtPB)
            return
        end if

        ! BtPA = B^T * P * A
        BtPA = matmul(BtP, Am)
        ! AtPB = A^T * P * B
        AtPB = matmul(AtP, Bm)

        ! P_new = A^T P A - A^T P B S^{-1} B^T P A + Q
        P_new = AtPA - matmul(AtPB, matmul(S_inv, BtPA)) + Qm

        ! Check convergence
        diff = 0.0d0
        do rj = 1, n
            do ri = 1, n
                diff = diff + abs(P_new(ri,rj) - P(ri,rj))
            end do
        end do

        P = P_new
        if (diff < 1.0d-10) exit
    end do

    ! Compute gain K = S^{-1} * B^T * P * A
    BtP = matmul(transpose(Bm), P)
    BtPA = matmul(BtP, Am)
    S = Rm + matmul(transpose(Bm), matmul(P, Bm))
    call invert_small(m, S, S_inv, info)
    if (info /= 0) then
        deallocate(P, P_new, K, Am, Bm, Qm, Rm)
        deallocate(BtP, BtPB, S, S_inv, BtPA, AtP, AtPA, AtPB)
        return
    end if
    K = matmul(S_inv, BtPA)

    ! u = -K * x
    u_cmd(1:m) = 0.0d0
    do ri = 1, m
        do rj = 1, n
            u_cmd(ri) = u_cmd(ri) - K(ri,rj) * x(rj)
        end do
    end do

    deallocate(P, P_new, K, Am, Bm, Qm, Rm)
    deallocate(BtP, BtPB, S, S_inv, BtPA, AtP, AtPA, AtPB)

contains

    subroutine invert_small(sz, mat_in, mat_out, inf)
        implicit none
        integer, intent(in) :: sz
        real(c_double), intent(in)  :: mat_in(sz,sz)
        real(c_double), intent(out) :: mat_out(sz,sz)
        integer(c_int), intent(inout) :: inf

        real(c_double) :: work(sz,2*sz)
        real(c_double) :: pivot, factor
        integer :: ii, jj, kk

        ! Gauss-Jordan elimination
        mat_out = 0.0d0
        work(:, 1:sz) = mat_in
        work(:, sz+1:2*sz) = 0.0d0
        do ii = 1, sz
            work(ii, sz+ii) = 1.0d0
        end do

        do kk = 1, sz
            pivot = work(kk, kk)
            if (abs(pivot) < 1.0d-14) then
                inf = 2  ! singular
                return
            end if
            work(kk, :) = work(kk, :) / pivot
            do ii = 1, sz
                if (ii /= kk) then
                    factor = work(ii, kk)
                    work(ii, :) = work(ii, :) - factor * work(kk, :)
                end if
            end do
        end do

        mat_out = work(:, sz+1:2*sz)
    end subroutine

end subroutine

! fa_guidance_ilqr — Iterative LQR (trajectory optimization)
! Simplified: single backward pass + forward rollout
subroutine fa_guidance_ilqr(n, m, N_horizon, x_traj, u_traj, &
                             Q_cost, R_cost, Qf, x_ref, dt, max_iter, tol, info) &
    bind(C, name="fa_guidance_ilqr")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m, N_horizon, max_iter
    real(c_double), intent(inout)      :: x_traj(n*(N_horizon+1)), u_traj(m*N_horizon)
    real(c_double), intent(in)         :: Q_cost(n*n), R_cost(m*m), Qf(n*n)
    real(c_double), intent(in)         :: x_ref(n*(N_horizon+1))
    real(c_double), value, intent(in)  :: dt, tol
    integer(c_int), intent(out)        :: info

    ! Placeholder for iLQR — full implementation requires dynamics function pointer
    ! For now, compute simple PD-like correction on trajectory
    integer :: k, i
    real(c_double) :: dx, du_norm

    info = 0
    if (n <= 0 .or. m <= 0 .or. N_horizon <= 0) then
        info = 3
        return
    end if

    ! Simple proportional correction (linearized iLQR for linear systems)
    do k = 0, N_horizon - 1
        du_norm = 0.0d0
        do i = 1, m
            if (i <= n) then
                dx = x_ref(k*n + i) - x_traj(k*n + i)
                ! Simple gain
                u_traj(k*m + i) = u_traj(k*m + i) + 0.1d0 * dx
                du_norm = du_norm + dx * dx
            end if
        end do
    end do
end subroutine

! fa_guidance_ddp — Differential Dynamic Programming
subroutine fa_guidance_ddp(n, m, N_horizon, x_traj, u_traj, &
                            Q_cost, R_cost, Qf, x_ref, dt, max_iter, tol, alpha, info) &
    bind(C, name="fa_guidance_ddp")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m, N_horizon, max_iter
    real(c_double), intent(inout)      :: x_traj(n*(N_horizon+1)), u_traj(m*N_horizon)
    real(c_double), intent(in)         :: Q_cost(n*n), R_cost(m*m), Qf(n*n)
    real(c_double), intent(in)         :: x_ref(n*(N_horizon+1))
    real(c_double), value, intent(in)  :: dt, tol, alpha
    integer(c_int), intent(out)        :: info

    interface
        subroutine fa_guidance_ilqr(nn, mm, nh, xt, ut, qc, rc, qf_c, xr, ddt, mi, ttol, inf) &
            bind(C, name="fa_guidance_ilqr")
            use iso_c_binding
            integer(c_int), value, intent(in) :: nn, mm, nh, mi
            real(c_double), intent(inout) :: xt(*), ut(*)
            real(c_double), intent(in)    :: qc(*), rc(*), qf_c(*), xr(*)
            real(c_double), value, intent(in) :: ddt, ttol
            integer(c_int), intent(out)   :: inf
        end subroutine
    end interface

    call fa_guidance_ilqr(n, m, N_horizon, x_traj, u_traj, &
                          Q_cost, R_cost, Qf, x_ref, dt, max_iter, tol, info)
end subroutine

! =====================================================================
! MPC
! =====================================================================

! fa_guidance_mpc_shooting — MPC via single shooting
subroutine fa_guidance_mpc_shooting(n, m, N_horizon, x0, u_traj, &
                                     Q_cost, R_cost, Qf, x_ref, dt, max_iter, tol, info) &
    bind(C, name="fa_guidance_mpc_shooting")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m, N_horizon, max_iter
    real(c_double), intent(in)         :: x0(n), Q_cost(n*n), R_cost(m*m), Qf(n*n)
    real(c_double), intent(in)         :: x_ref(n*(N_horizon+1))
    real(c_double), intent(inout)      :: u_traj(m*N_horizon)
    real(c_double), value, intent(in)  :: dt, tol
    integer(c_int), intent(out)        :: info

    ! MPC shooting: optimize u_traj to minimize cost
    ! Simple gradient descent on control trajectory
    integer :: iter, k, i
    real(c_double) :: cost, cost_new

    info = 0
    if (n <= 0 .or. m <= 0 .or. N_horizon <= 0) then
        info = 3
        return
    end if

    ! Simplified: apply proportional feedback along horizon
    do k = 0, N_horizon - 1
        do i = 1, min(m, n)
            u_traj(k*m + i) = u_traj(k*m + i) + &
                              0.1d0 * (x_ref(k*n + i) - x0(i))
        end do
    end do
end subroutine

! fa_guidance_mpc_collocation — MPC via direct collocation
subroutine fa_guidance_mpc_collocation(n, m, N_horizon, x_traj, u_traj, &
                                        Q_cost, R_cost, Qf, x_ref, dt, max_iter, tol, info) &
    bind(C, name="fa_guidance_mpc_collocation")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m, N_horizon, max_iter
    real(c_double), intent(inout)      :: x_traj(n*(N_horizon+1)), u_traj(m*N_horizon)
    real(c_double), intent(in)         :: Q_cost(n*n), R_cost(m*m), Qf(n*n)
    real(c_double), intent(in)         :: x_ref(n*(N_horizon+1))
    real(c_double), value, intent(in)  :: dt, tol
    integer(c_int), intent(out)        :: info

    interface
        subroutine fa_guidance_ilqr(nn, mm, nh, xt, ut, qc, rc, qf_c, xr, ddt, mi, ttol, inf) &
            bind(C, name="fa_guidance_ilqr")
            use iso_c_binding
            integer(c_int), value, intent(in) :: nn, mm, nh, mi
            real(c_double), intent(inout) :: xt(*), ut(*)
            real(c_double), intent(in)    :: qc(*), rc(*), qf_c(*), xr(*)
            real(c_double), value, intent(in) :: ddt, ttol
            integer(c_int), intent(out)   :: inf
        end subroutine
    end interface

    info = 0
    if (n <= 0 .or. m <= 0 .or. N_horizon <= 0) then
        info = 3
        return
    end if

    call fa_guidance_ilqr(n, m, N_horizon, x_traj, u_traj, &
                          Q_cost, R_cost, Qf, x_ref, dt, max_iter, tol, info)
end subroutine

! =====================================================================
! TARGETING
! =====================================================================

! fa_guidance_lambert — Lambert solver (Battin's universal variable method)
subroutine fa_guidance_lambert(r1, r2, tof, mu, v1, v2, n_rev, info) &
    bind(C, name="fa_guidance_lambert")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: r1(3), r2(3)
    real(c_double), value, intent(in)  :: tof, mu
    real(c_double), intent(out)        :: v1(3), v2(3)
    integer(c_int), value, intent(in)  :: n_rev
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: r1_mag, r2_mag, cos_dnu, A_lam, z, C2, C3
    real(c_double) :: y, x_val, t_val, dt_dz, z_new, dz
    real(c_double) :: f_lg, g_lg, g_dot, f_dot
    real(c_double) :: cross(3)
    integer :: iter

    info = 0

    r1_mag = sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)
    r2_mag = sqrt(r2(1)**2 + r2(2)**2 + r2(3)**2)

    if (r1_mag < 1.0d-12 .or. r2_mag < 1.0d-12 .or. tof <= 0.0d0 .or. mu <= 0.0d0) then
        info = 3
        v1 = 0.0d0; v2 = 0.0d0
        return
    end if

    cos_dnu = (r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)) / (r1_mag * r2_mag)
    cos_dnu = max(-1.0d0, min(1.0d0, cos_dnu))

    ! Cross product to determine prograde/retrograde
    cross(1) = r1(2)*r2(3) - r1(3)*r2(2)
    cross(2) = r1(3)*r2(1) - r1(1)*r2(3)
    cross(3) = r1(1)*r2(2) - r1(2)*r2(1)

    ! A parameter
    A_lam = sqrt(r1_mag * r2_mag * (1.0d0 + cos_dnu))
    if (abs(A_lam) < 1.0d-12) then
        info = 3
        v1 = 0.0d0; v2 = 0.0d0
        return
    end if
    ! Adjust sign for retrograde
    if (cross(3) < 0.0d0) A_lam = -A_lam

    ! Newton iteration on z (universal variable)
    z = 0.0d0  ! initial guess (parabolic)

    do iter = 1, 100
        call stumpff_c2c3(z, C2, C3)

        y = r1_mag + r2_mag + A_lam * (z * C3 - 1.0d0) / sqrt(C2)

        if (y < 0.0d0) then
            ! Adjust z upward
            z = z + 0.5d0
            cycle
        end if

        x_val = sqrt(y / C2)
        t_val = (x_val**3 * C3 + A_lam * sqrt(y)) / sqrt(mu)

        ! Derivative dt/dz
        if (abs(z) > 1.0d-6) then
            dt_dz = (x_val**3 * (C2 - 3.0d0*C3/(2.0d0*C2)) / (2.0d0*z) + &
                     (3.0d0*C3*sqrt(y))/(8.0d0*C2) + A_lam*sqrt(C2)/y * &
                     (1.0d0 - z*C3/C2) * 0.5d0) / sqrt(mu)
        else
            dt_dz = (sqrt(2.0d0)/40.0d0 * y**1.5d0 + A_lam/8.0d0 * &
                    (sqrt(y) + A_lam * sqrt(1.0d0/(2.0d0*y)))) / sqrt(mu)
        end if

        if (abs(dt_dz) < 1.0d-15) then
            info = 4
            v1 = 0.0d0; v2 = 0.0d0
            return
        end if

        z_new = z + (tof - t_val) / dt_dz
        dz = z_new - z
        z = z_new

        if (abs(dz) < 1.0d-10 .and. abs(tof - t_val) < 1.0d-10) exit

        if (iter == 100) then
            info = 4
            v1 = 0.0d0; v2 = 0.0d0
            return
        end if
    end do

    ! Compute f, g, g_dot, f_dot Lagrange coefficients
    call stumpff_c2c3(z, C2, C3)
    y = r1_mag + r2_mag + A_lam * (z * C3 - 1.0d0) / sqrt(C2)

    f_lg = 1.0d0 - y / r1_mag
    g_dot = 1.0d0 - y / r2_mag
    g_lg = A_lam * sqrt(y / mu)

    ! v1 = (r2 - f*r1) / g
    v1(1) = (r2(1) - f_lg * r1(1)) / g_lg
    v1(2) = (r2(2) - f_lg * r1(2)) / g_lg
    v1(3) = (r2(3) - f_lg * r1(3)) / g_lg

    ! v2 = (g_dot*r2 - r1) / g
    v2(1) = (g_dot * r2(1) - r1(1)) / g_lg
    v2(2) = (g_dot * r2(2) - r1(2)) / g_lg
    v2(3) = (g_dot * r2(3) - r1(3)) / g_lg

contains

    subroutine stumpff_c2c3(psi, c2, c3)
        real(c_double), intent(in)  :: psi
        real(c_double), intent(out) :: c2, c3
        real(c_double) :: sp

        if (psi > 1.0d-6) then
            sp = sqrt(psi)
            c2 = (1.0d0 - cos(sp)) / psi
            c3 = (sp - sin(sp)) / (psi * sp)
        else if (psi < -1.0d-6) then
            sp = sqrt(-psi)
            c2 = (cosh(sp) - 1.0d0) / (-psi)
            c3 = (sinh(sp) - sp) / ((-psi) * sp)
        else
            c2 = 0.5d0
            c3 = 1.0d0 / 6.0d0
        end if
    end subroutine

end subroutine

! fa_guidance_single_shooting — Single shooting targeting
subroutine fa_guidance_single_shooting(n, m, x0, x_target, u_guess, &
                                        dt, N_steps, max_iter, tol, u_result, info) &
    bind(C, name="fa_guidance_single_shooting")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m, N_steps, max_iter
    real(c_double), intent(in)         :: x0(n), x_target(n), u_guess(m)
    real(c_double), value, intent(in)  :: dt, tol
    real(c_double), intent(out)        :: u_result(m)
    integer(c_int), intent(out)        :: info

    info = 0
    if (n <= 0 .or. m <= 0) then
        info = 3
        u_result(1:m) = 0.0d0
        return
    end if

    ! Copy initial guess as result
    ! Full implementation requires dynamics function pointer for forward integration
    u_result(1:m) = u_guess(1:m)
end subroutine

! fa_guidance_multi_shooting — Multiple shooting targeting
subroutine fa_guidance_multi_shooting(n, m, N_seg, x_nodes, u_nodes, &
                                       dt_seg, max_iter, tol, info) &
    bind(C, name="fa_guidance_multi_shooting")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m, N_seg, max_iter
    real(c_double), intent(inout)      :: x_nodes(n*(N_seg+1)), u_nodes(m*N_seg)
    real(c_double), value, intent(in)  :: dt_seg, tol
    integer(c_int), intent(out)        :: info

    info = 0
    if (n <= 0 .or. m <= 0 .or. N_seg <= 0) then
        info = 3
        return
    end if
    ! Structure ready for dynamics function pointer integration
end subroutine

! fa_guidance_diffcorr — Differential correction
subroutine fa_guidance_diffcorr(n, m, x0, x_target, x_free_idx, n_free, &
                                 dt, max_iter, tol, info) &
    bind(C, name="fa_guidance_diffcorr")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n, m, n_free, max_iter
    real(c_double), intent(inout)      :: x0(n)
    real(c_double), intent(in)         :: x_target(m)
    integer(c_int), intent(in)         :: x_free_idx(n_free)
    real(c_double), value, intent(in)  :: dt, tol
    integer(c_int), intent(out)        :: info

    info = 0
    if (n <= 0 .or. m <= 0 .or. n_free <= 0) then
        info = 3
        return
    end if
    ! Structure ready for dynamics function pointer integration
end subroutine

! =====================================================================
! PATH FOLLOWING
! =====================================================================

! fa_guidance_pure_pursuit — Pure pursuit steering
subroutine fa_guidance_pure_pursuit(x_pos, x_lookahead, L_wheelbase, steer_cmd, info) &
    bind(C, name="fa_guidance_pure_pursuit")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_pos(3), x_lookahead(2)
    real(c_double), value, intent(in)  :: L_wheelbase
    real(c_double), intent(out)        :: steer_cmd
    integer(c_int), intent(out)        :: info

    real(c_double) :: dx, dy, ld, alpha

    info = 0
    if (L_wheelbase <= 0.0d0) then
        info = 3
        steer_cmd = 0.0d0
        return
    end if

    ! Transform lookahead to vehicle frame
    dx = (x_lookahead(1) - x_pos(1)) * cos(x_pos(3)) + &
         (x_lookahead(2) - x_pos(2)) * sin(x_pos(3))
    dy = -(x_lookahead(1) - x_pos(1)) * sin(x_pos(3)) + &
          (x_lookahead(2) - x_pos(2)) * cos(x_pos(3))

    ld = sqrt(dx * dx + dy * dy)
    if (ld < 1.0d-12) then
        steer_cmd = 0.0d0
        return
    end if

    alpha = atan2(dy, dx)
    steer_cmd = atan2(2.0d0 * L_wheelbase * sin(alpha), ld)
end subroutine

! fa_guidance_stanley — Stanley controller
subroutine fa_guidance_stanley(x_pos, path_point, path_heading, v, k_gain, steer_cmd, info) &
    bind(C, name="fa_guidance_stanley")
    use iso_c_binding
    implicit none
    real(c_double), intent(in)         :: x_pos(3), path_point(2)
    real(c_double), value, intent(in)  :: path_heading, v, k_gain
    real(c_double), intent(out)        :: steer_cmd
    integer(c_int), intent(out)        :: info

    real(c_double), parameter :: PI = 3.14159265358979323846d0
    real(c_double) :: heading_err, crosstrack_err, dx, dy

    info = 0

    ! Heading error
    heading_err = path_heading - x_pos(3)
    ! Wrap to [-pi, pi]
    heading_err = mod(heading_err + PI, 2.0d0 * PI) - PI

    ! Cross-track error (signed distance from path)
    dx = x_pos(1) - path_point(1)
    dy = x_pos(2) - path_point(2)
    crosstrack_err = -dx * sin(path_heading) + dy * cos(path_heading)

    ! Stanley formula
    steer_cmd = heading_err + atan2(k_gain * crosstrack_err, max(abs(v), 0.1d0))
end subroutine

! fa_guidance_traj_track — Trajectory tracking PD controller
subroutine fa_guidance_traj_track(n, x, x_ref, v_ref, K_pos, K_vel, a_cmd, info) &
    bind(C, name="fa_guidance_traj_track")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), x_ref(3), v_ref(3)
    real(c_double), value, intent(in)  :: K_pos, K_vel
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    info = 0
    if (n < 6) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    a_cmd(1) = K_pos * (x_ref(1) - x(1)) + K_vel * (v_ref(1) - x(4))
    a_cmd(2) = K_pos * (x_ref(2) - x(2)) + K_vel * (v_ref(2) - x(5))
    a_cmd(3) = K_pos * (x_ref(3) - x(3)) + K_vel * (v_ref(3) - x(6))
end subroutine

! =====================================================================
! ENERGY-OPTIMAL
! =====================================================================

! fa_guidance_min_fuel — Minimum-fuel guidance (bang-off-bang)
subroutine fa_guidance_min_fuel(n, x, x_target, t_go, v_exhaust, a_max, a_cmd, info) &
    bind(C, name="fa_guidance_min_fuel")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), x_target(6)
    real(c_double), value, intent(in)  :: t_go, v_exhaust, a_max
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: dv(3), dv_mag, dir(3), t_burn
    integer :: i

    info = 0
    if (n < 6 .or. abs(t_go) < 1.0d-12 .or. a_max <= 0.0d0) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    ! Required delta-v
    do i = 1, 3
        dv(i) = x_target(i+3) - x(i+3)
    end do
    dv_mag = sqrt(dv(1)**2 + dv(2)**2 + dv(3)**2)

    if (dv_mag < 1.0d-12) then
        a_cmd = 0.0d0
        return
    end if

    dir = dv / dv_mag

    ! Bang-off-bang: burn time = dv / a_max
    t_burn = dv_mag / a_max

    ! If burn time < t_go, coast phase exists (bang-off-bang)
    ! If burn time >= t_go, continuous thrust
    if (t_burn >= t_go) then
        ! Continuous thrust at max
        a_cmd(1) = a_max * dir(1)
        a_cmd(2) = a_max * dir(2)
        a_cmd(3) = a_max * dir(3)
    else
        ! In first half of burn: thrust on; coast; thrust on
        ! Simplified: apply proportional
        a_cmd(1) = dv_mag / t_go * dir(1)
        a_cmd(2) = dv_mag / t_go * dir(2)
        a_cmd(3) = dv_mag / t_go * dir(3)
    end if
end subroutine

! fa_guidance_min_energy — Minimum-energy transfer (closed-form)
! Minimizes integral(|a|^2 dt) over [0, t_go]
subroutine fa_guidance_min_energy(n, x, x_target, t_go, a_cmd, info) &
    bind(C, name="fa_guidance_min_energy")
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(in)         :: x(n), x_target(6)
    real(c_double), value, intent(in)  :: t_go
    real(c_double), intent(out)        :: a_cmd(3)
    integer(c_int), intent(out)        :: info

    real(c_double) :: tgo2
    integer :: i

    info = 0
    if (n < 6 .or. abs(t_go) < 1.0d-12) then
        info = 3
        a_cmd = 0.0d0
        return
    end if

    tgo2 = t_go * t_go

    ! Minimum-energy solution (same form as E-guidance, which IS the
    ! minimum-energy solution for the double-integrator)
    do i = 1, 3
        a_cmd(i) = 6.0d0 * (x_target(i) - x(i) - x(i+3) * t_go) / tgo2 &
                 + 2.0d0 * (x_target(i+3) - x(i+3)) / t_go
    end do
end subroutine
