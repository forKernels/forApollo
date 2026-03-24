! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_estimate.f90 — State estimation: KF, EKF, IEKF, UKF, ESKF, SR variants
!
! Core state estimation routines with ternary measurement gating.
! All matrices are flat 1D row-major at the bind(C) boundary.
!
! Entry points:
!   fa_kf_predict      — Linear Kalman predict: x = F*x, P = F*P*F^T + Q
!   fa_kf_update       — Linear Kalman update with ternary gating
!   fa_ekf_predict     — Extended Kalman predict via propagator + STM
!   fa_ekf_update      — Extended Kalman update with ternary gating
!   fa_iekf_update     — Iterated EKF update (re-linearization loop)
!   fa_ukf_predict     — Unscented Kalman predict (sigma points)
!   fa_ukf_update      — Unscented Kalman update (sigma points)
!   fa_eskf_predict    — Error-State Kalman predict (Apollo heritage)
!   fa_eskf_update     — Error-State Kalman update on error state
!   fa_eskf_inject     — Inject error state into nominal
!   fa_srekf_predict   — Square-Root EKF predict (Cholesky factor)
!   fa_srekf_update    — Square-Root EKF update (Cholesky factor)
!   fa_srukf_predict   — Square-Root UKF predict (Cholesky factor)
!   fa_srukf_update    — Square-Root UKF update (Cholesky factor)
!
! Ternary measurement gating:
!   For each measurement i, compute chi2 = y(i)^2 / S(i,i)
!     chi2 > 16.0 (4-sigma) → validity(i) = -1, rejected
!     chi2 > 9.0  (3-sigma) → validity(i) =  0, uncertain (skip)
!     else                   → validity(i) = +1, fused
!
! Internal matrix helpers (no bind(C), no external deps):
!   mat_multiply, mat_transpose, mat_add, mat_subtract, mat_identity,
!   mat_solve_symm, mat_cholesky, mat_outer, mat_scale, vec_norm
!
! Joseph form for P update: P = (I-KH)*P*(I-KH)^T + K*R*K^T
! This is numerically superior to simple (I-KH)*P.
!
! Error codes: info = 0 ok, 2 not positive definite, 3 invalid input,
!              4 propagation/observation failure or max iterations reached.
!
! Row-major convention: element (row, col) = M((row-1)*ncol + col)


! ==========================================================================
! fa_kf_predict — Linear Kalman predict step
!
! In-place update:
!   x_new = F * x
!   P_new = F * P * F^T + Q
! ==========================================================================
subroutine fa_kf_predict(n, x, P, F, Q, info) bind(C, name="fa_kf_predict")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: P(n*n)
    real(c_double), intent(in)         :: F(n*n)
    real(c_double), intent(in)         :: Q(n*n)
    integer(c_int), intent(out)        :: info

    real(c_double), allocatable :: x_new(:), FP(:), FPFT(:), FT(:)
    integer :: nn

    info = 0
    nn = n * n

    ! --- Input validation ---
    if (n <= 0) then
        info = 3
        return
    end if

    allocate(x_new(n), FP(nn), FPFT(nn), FT(nn))

    ! x_new = F * x (matrix-vector: treat x as n x 1)
    call mat_vec_multiply(F, x, x_new, n, n)
    x(1:n) = x_new(1:n)

    ! P_new = F * P * F^T + Q
    call mat_multiply(F, P, FP, n, n, n)
    call mat_transpose(F, FT, n, n)
    call mat_multiply(FP, FT, FPFT, n, n, n)
    call mat_add(FPFT, Q, P, n, n)

    deallocate(x_new, FP, FPFT, FT)

end subroutine fa_kf_predict


! ==========================================================================
! fa_kf_update — Linear Kalman update step with ternary gating
!
! Innovation: y = z - H*x
! Innovation covariance: S = H*P*H^T + R
! Ternary gating per measurement on chi2 = y(i)^2 / S(i,i)
! Kalman gain: K = P * H^T * S^-1
! Joseph form: P = (I-KH)*P*(I-KH)^T + K*R*K^T
! ==========================================================================
subroutine fa_kf_update(n, m, x, P, z, H, R, validity, info) bind(C, name="fa_kf_update")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, m
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: P(n*n)
    real(c_double), intent(in)         :: z(m)
    real(c_double), intent(in)         :: H(m*n)
    real(c_double), intent(in)         :: R(m*m)
    integer(c_int), intent(out)        :: validity(m)
    integer(c_int), intent(out)        :: info

    real(c_double), allocatable :: z_pred(:), y(:)
    real(c_double), allocatable :: HT(:), PH(:), S(:), S_inv(:)
    real(c_double), allocatable :: K(:), IKH(:), eye_n(:)
    real(c_double), allocatable :: IKHP(:), IKHT(:), IKHP_IKHT(:)
    real(c_double), allocatable :: KR(:), KT(:), KRKT(:)
    real(c_double) :: chi2
    integer :: i, nn, mm, solve_info
    integer :: n_valid
    ! Filtered versions for gating
    real(c_double), allocatable :: H_filt(:), R_filt(:), y_filt(:)
    integer :: j, row_out

    info = 0
    nn = n * n
    mm = m * m

    ! --- Input validation ---
    if (n <= 0 .or. m <= 0) then
        info = 3
        return
    end if

    allocate(z_pred(m), y(m))
    allocate(HT(n*m), PH(n*m), S(mm), S_inv(mm))
    allocate(eye_n(nn))

    ! z_pred = H * x
    call mat_vec_multiply(H, x, z_pred, m, n)

    ! Innovation: y = z - z_pred
    y(1:m) = z(1:m) - z_pred(1:m)

    ! S = H * P * H^T + R
    call mat_transpose(H, HT, m, n)           ! HT is n x m
    call mat_multiply(P, HT, PH, n, n, m)     ! PH is n x m
    call mat_multiply(H, PH, S, m, n, m)      ! S is m x m
    call mat_add(S, R, S, m, m)

    ! --- Ternary gating ---
    validity(1:m) = 1   ! default: fuse
    do i = 1, m
        ! S(i,i) is at position (i-1)*m + i in flat row-major
        if (S((i-1)*m + i) > 0.0d0) then
            chi2 = y(i) * y(i) / S((i-1)*m + i)
        else
            ! Degenerate covariance — reject
            chi2 = 999.0d0
        end if

        if (chi2 > 16.0d0) then
            validity(i) = -1   ! rejected (> 4-sigma)
        else if (chi2 > 9.0d0) then
            validity(i) = 0    ! uncertain (> 3-sigma), skip
        end if
    end do

    ! Count valid measurements
    n_valid = 0
    do i = 1, m
        if (validity(i) == 1) n_valid = n_valid + 1
    end do

    ! If no valid measurements, skip update
    if (n_valid == 0) then
        deallocate(z_pred, y, HT, PH, S, S_inv, eye_n)
        return
    end if

    ! Build filtered H, R, y with only valid measurements
    allocate(H_filt(n_valid * n), R_filt(n_valid * n_valid), y_filt(n_valid))

    row_out = 0
    do i = 1, m
        if (validity(i) == 1) then
            row_out = row_out + 1
            y_filt(row_out) = y(i)
            ! Copy row i of H to row row_out of H_filt
            do j = 1, n
                H_filt((row_out - 1) * n + j) = H((i - 1) * n + j)
            end do
        end if
    end do

    ! Build filtered R (extract valid rows and columns)
    R_filt(1:n_valid*n_valid) = 0.0d0
    row_out = 0
    do i = 1, m
        if (validity(i) /= 1) cycle
        row_out = row_out + 1
        j = 0
        block
            integer :: k, col_out
            col_out = 0
            do k = 1, m
                if (validity(k) /= 1) cycle
                col_out = col_out + 1
                R_filt((row_out - 1) * n_valid + col_out) = R((i - 1) * m + k)
            end do
        end block
    end do

    ! Recompute S_filt = H_filt * P * H_filt^T + R_filt
    deallocate(HT, PH, S, S_inv)
    allocate(HT(n * n_valid), PH(n * n_valid))
    allocate(S(n_valid * n_valid), S_inv(n_valid * n_valid))

    call mat_transpose(H_filt, HT, n_valid, n)       ! HT is n x n_valid
    call mat_multiply(P, HT, PH, n, n, n_valid)      ! PH is n x n_valid
    call mat_multiply(H_filt, PH, S, n_valid, n, n_valid)  ! S is n_valid x n_valid
    call mat_add(S, R_filt, S, n_valid, n_valid)

    ! K = P * H_filt^T * S^-1
    ! Solve S * K^T = (P * H_filt^T)^T  →  actually compute S_inv then K = PH * S_inv
    ! For numerical stability, solve S * X = I to get S_inv
    call mat_identity(S_inv, n_valid)
    call mat_solve_symm(S, S_inv, n_valid, n_valid, solve_info)
    if (solve_info /= 0) then
        info = 4
        deallocate(z_pred, y, HT, PH, S, S_inv, eye_n)
        deallocate(H_filt, R_filt, y_filt)
        return
    end if

    ! K = PH * S_inv  (n x n_valid) * (n_valid x n_valid) = n x n_valid
    allocate(K(n * n_valid))
    call mat_multiply(PH, S_inv, K, n, n_valid, n_valid)

    ! x = x + K * y_filt
    block
        real(c_double) :: dx(n)
        call mat_vec_multiply(K, y_filt, dx, n, n_valid)
        x(1:n) = x(1:n) + dx(1:n)
    end block

    ! Joseph form: P = (I - K*H_filt) * P * (I - K*H_filt)^T + K * R_filt * K^T
    allocate(IKH(nn), IKHP(nn), IKHT(nn), IKHP_IKHT(nn))
    allocate(KR(n * n_valid), KT(n_valid * n), KRKT(nn))

    ! IKH = I - K * H_filt
    call mat_identity(eye_n, n)
    block
        real(c_double), allocatable :: KH(:)
        allocate(KH(nn))
        call mat_multiply(K, H_filt, KH, n, n_valid, n)   ! KH is n x n
        IKH(1:nn) = eye_n(1:nn) - KH(1:nn)
        deallocate(KH)
    end block

    ! IKHP = IKH * P
    call mat_multiply(IKH, P, IKHP, n, n, n)

    ! IKHT = IKH^T
    call mat_transpose(IKH, IKHT, n, n)

    ! IKHP_IKHT = IKHP * IKHT
    call mat_multiply(IKHP, IKHT, IKHP_IKHT, n, n, n)

    ! KR = K * R_filt  (n x n_valid) * (n_valid x n_valid) = n x n_valid
    call mat_multiply(K, R_filt, KR, n, n_valid, n_valid)

    ! KT = K^T  (n_valid x n)
    call mat_transpose(K, KT, n, n_valid)

    ! KRKT = KR * KT  (n x n_valid) * (n_valid x n) = n x n
    call mat_multiply(KR, KT, KRKT, n, n_valid, n)

    ! P = IKHP_IKHT + KRKT
    call mat_add(IKHP_IKHT, KRKT, P, n, n)

    deallocate(z_pred, y, HT, PH, S, S_inv, eye_n)
    deallocate(H_filt, R_filt, y_filt)
    deallocate(K, IKH, IKHP, IKHT, IKHP_IKHT, KR, KT, KRKT)

end subroutine fa_kf_update


! ==========================================================================
! fa_ekf_predict — Extended Kalman predict step
!
! Propagates state using fa_propagate (or f_ptr if non-null).
! Propagates covariance via STM: P = Phi * P * Phi^T + Q
!   where Phi comes from fa_propagate_stm for built-in models.
! For custom f_ptr without df_ptr: uses identity Phi, sets info=4.
! ==========================================================================
subroutine fa_ekf_predict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, n_steps, info) &
    bind(C, name="fa_ekf_predict")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, model_id, np, n_steps
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: P(n*n)
    type(c_funptr), value              :: f_ptr
    type(c_funptr), value              :: df_ptr
    real(c_double), intent(in)         :: Q(n*n)
    real(c_double), value, intent(in)  :: dt
    real(c_double), intent(in)         :: params(np)
    integer(c_int), intent(out)        :: info

    ! Interface for fa_propagate
    interface
        subroutine fa_propagate(n, x, u, nu, f_ptr, model_id, params, np, dt, n_steps, info) &
            bind(C, name="fa_propagate")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, nu, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(in)         :: u(nu), params(np)
            real(c_double), value, intent(in)  :: dt
            type(c_funptr), value              :: f_ptr
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_propagate_stm(n, x, phi, u, nu, f_ptr, df_ptr, model_id, params, np, dt, n_steps, info) &
            bind(C, name="fa_propagate_stm")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, nu, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n), phi(n*n)
            real(c_double), intent(in)         :: u(nu), params(np)
            real(c_double), value, intent(in)  :: dt
            type(c_funptr), value              :: f_ptr, df_ptr
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double), allocatable :: Phi(:), PhiP(:), PhiT(:), PhiPPhiT(:)
    real(c_double) :: u_dummy(1)
    integer :: nn, prop_info
    logical :: use_custom

    info = 0
    nn = n * n

    ! --- Input validation ---
    if (n <= 0) then
        info = 3
        return
    end if

    if (dt == 0.0d0) return

    use_custom = c_associated(f_ptr)

    allocate(Phi(nn), PhiP(nn), PhiT(nn), PhiPPhiT(nn))

    ! Initialize Phi to identity
    call mat_identity(Phi, n)

    ! Zero-length control input
    u_dummy(1) = 0.0d0

    if (use_custom .and. .not. c_associated(df_ptr)) then
        ! Custom dynamics without Jacobian — propagate state only, Phi = I
        ! This is a fallback; covariance propagation will be approximate
        call fa_propagate(n, x, u_dummy, 0, f_ptr, model_id, params, np, dt, n_steps, prop_info)
        if (prop_info /= 0) then
            info = 4
            deallocate(Phi, PhiP, PhiT, PhiPPhiT)
            return
        end if
        ! Phi stays identity — flag as approximate
        info = 4
    else
        ! Use STM propagation (built-in model or custom with Jacobian)
        call fa_propagate_stm(n, x, Phi, u_dummy, 0, f_ptr, df_ptr, model_id, &
                              params, np, dt, n_steps, prop_info)
        if (prop_info /= 0) then
            info = 4
            deallocate(Phi, PhiP, PhiT, PhiPPhiT)
            return
        end if
    end if

    ! P_new = Phi * P * Phi^T + Q
    call mat_multiply(Phi, P, PhiP, n, n, n)
    call mat_transpose(Phi, PhiT, n, n)
    call mat_multiply(PhiP, PhiT, PhiPPhiT, n, n, n)
    call mat_add(PhiPPhiT, Q, P, n, n)

    deallocate(Phi, PhiP, PhiT, PhiPPhiT)

end subroutine fa_ekf_predict


! ==========================================================================
! fa_ekf_update — Extended Kalman update step with ternary gating
!
! Compute predicted measurement z_pred = h(x) via fa_observe_dispatch or h_ptr
! Compute Jacobian H = dh/dx via fa_observe_jacobian or dh_ptr
! Then standard Kalman update with ternary gating (same as fa_kf_update)
! ==========================================================================
subroutine fa_ekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info) &
    bind(C, name="fa_ekf_update")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, m, obs_id, nop
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: P(n*n)
    real(c_double), intent(in)         :: z(m)
    type(c_funptr), value              :: h_ptr
    type(c_funptr), value              :: dh_ptr
    real(c_double), intent(in)         :: R(m*m)
    real(c_double), intent(in)         :: obs_params(nop)
    integer(c_int), intent(out)        :: validity(m)
    integer(c_int), intent(out)        :: info

    ! Interfaces for observation dispatch
    interface
        subroutine fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info) &
            bind(C, name="fa_observe_dispatch")
            use iso_c_binding
            integer(c_int), value  :: obs_id, n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine

        subroutine fa_observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H, info) &
            bind(C, name="fa_observe_jacobian")
            use iso_c_binding
            integer(c_int), value  :: obs_id, n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: H(m*n)
            integer(c_int), intent(out) :: info
        end subroutine

        ! Interface for fa_kf_update (reuse linear update)
        subroutine fa_kf_update(n, m, x, P, z_inn, H, R, validity, info) &
            bind(C, name="fa_kf_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            real(c_double), intent(in)         :: z_inn(m)
            real(c_double), intent(in)         :: H(m*n)
            real(c_double), intent(in)         :: R(m*m)
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    ! Abstract interfaces for user-supplied function pointers
    abstract interface
        subroutine custom_observe_t(n, x, m, t, obs_params, nop, z_pred, info) bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine

        subroutine custom_obs_jacobian_t(n, x, m, t, obs_params, nop, H, info) bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: H(m*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    procedure(custom_observe_t), pointer     :: h_custom
    procedure(custom_obs_jacobian_t), pointer :: dh_custom
    logical :: use_custom_h, use_custom_dh

    real(c_double), allocatable :: z_pred(:), H(:), z_eff(:)
    integer :: obs_info

    info = 0

    ! --- Input validation ---
    if (n <= 0 .or. m <= 0) then
        info = 3
        return
    end if

    use_custom_h  = c_associated(h_ptr)
    use_custom_dh = c_associated(dh_ptr)
    if (use_custom_h)  call c_f_procpointer(h_ptr, h_custom)
    if (use_custom_dh) call c_f_procpointer(dh_ptr, dh_custom)

    allocate(z_pred(m), H(m * n), z_eff(m))

    ! Compute predicted measurement z_pred = h(x)
    if (use_custom_h) then
        call h_custom(n, x, m, 0.0d0, obs_params, nop, z_pred, obs_info)
    else
        call fa_observe_dispatch(obs_id, n, x, m, 0.0d0, obs_params, nop, z_pred, obs_info)
    end if
    if (obs_info /= 0) then
        info = 4
        deallocate(z_pred, H, z_eff)
        return
    end if

    ! Compute Jacobian H = dh/dx
    if (use_custom_dh) then
        call dh_custom(n, x, m, 0.0d0, obs_params, nop, H, obs_info)
    else
        call fa_observe_jacobian(obs_id, n, x, m, 0.0d0, obs_params, nop, H, obs_info)
    end if
    if (obs_info /= 0) then
        info = 4
        deallocate(z_pred, H, z_eff)
        return
    end if

    ! Build effective measurement: z_eff = z - z_pred + H*x
    ! This way fa_kf_update computes innovation as z_eff - H*x = z - z_pred
    ! which is the correct nonlinear innovation
    block
        real(c_double), allocatable :: Hx(:)
        integer :: i
        allocate(Hx(m))
        call mat_vec_multiply(H, x, Hx, m, n)
        do i = 1, m
            z_eff(i) = z(i) - z_pred(i) + Hx(i)
        end do
        deallocate(Hx)
    end block

    ! Delegate to linear update with the linearized quantities
    call fa_kf_update(n, m, x, P, z_eff, H, R, validity, info)

    deallocate(z_pred, H, z_eff)

end subroutine fa_ekf_update


! ==========================================================================
! Internal matrix helper subroutines
! No bind(C), no external dependencies. Will be replaced by forMath later.
! All operate on flat 1D row-major arrays.
! ==========================================================================

! --------------------------------------------------------------------------
! mat_vec_multiply — y = A * x, A is (m x n) flat row-major, x is (n), y is (m)
! --------------------------------------------------------------------------
subroutine mat_vec_multiply(A, x, y, m, n)
    use iso_c_binding
    implicit none
    integer, intent(in) :: m, n
    real(c_double), intent(in)  :: A(m*n), x(n)
    real(c_double), intent(out) :: y(m)
    integer :: i, j
    real(c_double) :: s

    do i = 1, m
        s = 0.0d0
        do j = 1, n
            s = s + A((i-1)*n + j) * x(j)
        end do
        y(i) = s
    end do
end subroutine mat_vec_multiply

! --------------------------------------------------------------------------
! mat_multiply — C = A * B
! A is (m x k) flat row-major, B is (k x n_col), C is (m x n_col)
! --------------------------------------------------------------------------
subroutine mat_multiply(A, B, C, m, k, n_col)
    use iso_c_binding
    implicit none
    integer, intent(in) :: m, k, n_col
    real(c_double), intent(in)  :: A(m*k), B(k*n_col)
    real(c_double), intent(out) :: C(m*n_col)
    integer :: i, j, l
    real(c_double) :: s

    do i = 1, m
        do j = 1, n_col
            s = 0.0d0
            do l = 1, k
                s = s + A((i-1)*k + l) * B((l-1)*n_col + j)
            end do
            C((i-1)*n_col + j) = s
        end do
    end do
end subroutine mat_multiply

! --------------------------------------------------------------------------
! mat_transpose — AT = A^T
! A is (m x n) flat row-major, AT is (n x m)
! --------------------------------------------------------------------------
subroutine mat_transpose(A, AT, m, n)
    use iso_c_binding
    implicit none
    integer, intent(in) :: m, n
    real(c_double), intent(in)  :: A(m*n)
    real(c_double), intent(out) :: AT(n*m)
    integer :: i, j

    do i = 1, m
        do j = 1, n
            AT((j-1)*m + i) = A((i-1)*n + j)
        end do
    end do
end subroutine mat_transpose

! --------------------------------------------------------------------------
! mat_add — C = A + B, all (m x n) flat row-major
! C may alias A or B.
! --------------------------------------------------------------------------
subroutine mat_add(A, B, C, m, n)
    use iso_c_binding
    implicit none
    integer, intent(in) :: m, n
    real(c_double), intent(in)  :: A(m*n), B(m*n)
    real(c_double), intent(out) :: C(m*n)
    integer :: i

    do i = 1, m * n
        C(i) = A(i) + B(i)
    end do
end subroutine mat_add

! --------------------------------------------------------------------------
! mat_identity — I = eye(n), flat row-major (n x n)
! --------------------------------------------------------------------------
subroutine mat_identity(eye, n)
    use iso_c_binding
    implicit none
    integer, intent(in) :: n
    real(c_double), intent(out) :: eye(n*n)
    integer :: i

    eye(1:n*n) = 0.0d0
    do i = 1, n
        eye((i-1)*n + i) = 1.0d0
    end do
end subroutine mat_identity

! --------------------------------------------------------------------------
! mat_solve_symm — Solve A * X = B in-place (X overwrites B)
!
! A is (n x n) symmetric positive definite, flat row-major
! B is (n x m_rhs) flat row-major, overwritten with solution X
!
! Uses Cholesky decomposition: A = L * L^T
! Then solve L * Y = B (forward substitution)
! Then solve L^T * X = Y (back substitution)
!
! For small systems this is efficient enough. The Zig layer will
! dispatch to forMath BLAS (dpotrf/dpotrs) for large systems.
!
! info = 0 success, 1 not positive definite
! --------------------------------------------------------------------------
subroutine mat_solve_symm(A, B, n, m_rhs, info)
    use iso_c_binding
    implicit none
    integer, intent(in)    :: n, m_rhs
    real(c_double), intent(in)    :: A(n*n)
    real(c_double), intent(inout) :: B(n*m_rhs)
    integer, intent(out)   :: info

    real(c_double), allocatable :: L(:)   ! n x n lower triangular
    real(c_double) :: s
    integer :: i, j, k

    info = 0

    allocate(L(n*n))
    L(1:n*n) = 0.0d0

    ! Cholesky factorization: A = L * L^T (row-major)
    do i = 1, n
        do j = 1, i
            s = A((i-1)*n + j)
            do k = 1, j - 1
                s = s - L((i-1)*n + k) * L((j-1)*n + k)
            end do
            if (i == j) then
                if (s <= 0.0d0) then
                    info = 1   ! Not positive definite
                    deallocate(L)
                    return
                end if
                L((i-1)*n + j) = sqrt(s)
            else
                L((i-1)*n + j) = s / L((j-1)*n + j)
            end if
        end do
    end do

    ! Forward substitution: solve L * Y = B (column by column of B)
    do k = 1, m_rhs
        do i = 1, n
            s = B((i-1)*m_rhs + k)
            do j = 1, i - 1
                s = s - L((i-1)*n + j) * B((j-1)*m_rhs + k)
            end do
            B((i-1)*m_rhs + k) = s / L((i-1)*n + i)
        end do
    end do

    ! Back substitution: solve L^T * X = Y (column by column of B)
    do k = 1, m_rhs
        do i = n, 1, -1
            s = B((i-1)*m_rhs + k)
            do j = i + 1, n
                s = s - L((j-1)*n + i) * B((j-1)*m_rhs + k)
            end do
            B((i-1)*m_rhs + k) = s / L((i-1)*n + i)
        end do
    end do

    deallocate(L)

end subroutine mat_solve_symm


! --------------------------------------------------------------------------
! mat_cholesky — Cholesky decomposition A = L * L^T
!
! A(n*n) input symmetric positive definite (flat row-major)
! L(n*n) output lower triangular (flat row-major)
! info = 0 success, 2 not positive definite
! --------------------------------------------------------------------------
subroutine mat_cholesky(A, L, n, info)
    use iso_c_binding
    implicit none
    integer, intent(in)  :: n
    real(c_double), intent(in)  :: A(n*n)
    real(c_double), intent(out) :: L(n*n)
    integer, intent(out) :: info
    integer :: i, j, k
    real(c_double) :: s

    info = 0
    L(1:n*n) = 0.0d0

    do i = 1, n
        do j = 1, i
            s = A((i-1)*n + j)
            do k = 1, j - 1
                s = s - L((i-1)*n + k) * L((j-1)*n + k)
            end do
            if (i == j) then
                if (s <= 0.0d0) then
                    info = 2  ! Not positive definite
                    return
                end if
                L((i-1)*n + j) = sqrt(s)
            else
                L((i-1)*n + j) = s / L((j-1)*n + j)
            end if
        end do
    end do

end subroutine mat_cholesky


! --------------------------------------------------------------------------
! mat_subtract — C = A - B, all (m x n) flat row-major
! C may alias A or B.
! --------------------------------------------------------------------------
subroutine mat_subtract(A, B, C, m, n)
    use iso_c_binding
    implicit none
    integer, intent(in) :: m, n
    real(c_double), intent(in)  :: A(m*n), B(m*n)
    real(c_double), intent(out) :: C(m*n)
    integer :: i

    do i = 1, m * n
        C(i) = A(i) - B(i)
    end do
end subroutine mat_subtract


! --------------------------------------------------------------------------
! mat_outer — C = a * b^T, a is (n), b is (m), C is (n x m) flat row-major
! --------------------------------------------------------------------------
subroutine mat_outer(a, b, C, n, m)
    use iso_c_binding
    implicit none
    integer, intent(in) :: n, m
    real(c_double), intent(in)  :: a(n), b(m)
    real(c_double), intent(out) :: C(n*m)
    integer :: i, j

    do i = 1, n
        do j = 1, m
            C((i-1)*m + j) = a(i) * b(j)
        end do
    end do
end subroutine mat_outer


! --------------------------------------------------------------------------
! mat_scale — A = alpha * A, flat array of length len
! --------------------------------------------------------------------------
subroutine mat_scale(A, alpha, len)
    use iso_c_binding
    implicit none
    integer, intent(in) :: len
    real(c_double), intent(inout) :: A(len)
    real(c_double), intent(in)    :: alpha
    integer :: i

    do i = 1, len
        A(i) = alpha * A(i)
    end do
end subroutine mat_scale


! --------------------------------------------------------------------------
! vec_norm — Euclidean norm of a vector
! --------------------------------------------------------------------------
function vec_norm(v, n) result(nrm)
    use iso_c_binding
    implicit none
    integer, intent(in) :: n
    real(c_double), intent(in) :: v(n)
    real(c_double) :: nrm
    integer :: i

    nrm = 0.0d0
    do i = 1, n
        nrm = nrm + v(i) * v(i)
    end do
    nrm = sqrt(nrm)
end function vec_norm


! ==========================================================================
! fa_iekf_update — Iterated Extended Kalman Filter update
!
! Re-linearizes at updated state estimate up to max_iter times.
! Converges when ||x_new - x_old|| < tol.
! Each iteration: re-evaluate h(x_i) and H(x_i), compute modified gain.
! info=4 if max_iter reached without convergence.
! ==========================================================================
subroutine fa_iekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, &
                           max_iter, tol, validity, info) &
    bind(C, name="fa_iekf_update")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, m, obs_id, nop, max_iter
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: P(n*n)
    real(c_double), intent(in)         :: z(m)
    type(c_funptr), value              :: h_ptr
    type(c_funptr), value              :: dh_ptr
    real(c_double), intent(in)         :: R(m*m)
    real(c_double), intent(in)         :: obs_params(nop)
    real(c_double), value, intent(in)  :: tol
    integer(c_int), intent(out)        :: validity(m)
    integer(c_int), intent(out)        :: info

    ! Interfaces for observation dispatch
    interface
        subroutine fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info) &
            bind(C, name="fa_observe_dispatch")
            use iso_c_binding
            integer(c_int), value  :: obs_id, n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine

        subroutine fa_observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H, info) &
            bind(C, name="fa_observe_jacobian")
            use iso_c_binding
            integer(c_int), value  :: obs_id, n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: H(m*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    abstract interface
        subroutine custom_observe_iekf_t(n, x, m, t, obs_params, nop, z_pred, info) bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine

        subroutine custom_obs_jac_iekf_t(n, x, m, t, obs_params, nop, H, info) bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: H(m*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    procedure(custom_observe_iekf_t), pointer     :: h_custom
    procedure(custom_obs_jac_iekf_t), pointer      :: dh_custom
    logical :: use_custom_h, use_custom_dh

    real(c_double), allocatable :: x0(:), xi(:), x_new(:), dx(:)
    real(c_double), allocatable :: Hi(:), z_pred_i(:), innov(:)
    real(c_double), allocatable :: HT(:), PHT(:), S(:), S_inv(:), K(:)
    real(c_double), allocatable :: eye_n(:), IKH(:), IKHP(:), IKHT(:), IKHP_IKHT(:)
    real(c_double), allocatable :: KR(:), KT(:), KRKT(:)
    real(c_double) :: chi2, diff_norm
    real(c_double) :: vec_norm
    integer :: iter, i, nn, mm, obs_info, solve_info

    info = 0
    nn = n * n
    mm = m * m

    if (n <= 0 .or. m <= 0) then
        info = 3
        return
    end if

    use_custom_h  = c_associated(h_ptr)
    use_custom_dh = c_associated(dh_ptr)
    if (use_custom_h)  call c_f_procpointer(h_ptr, h_custom)
    if (use_custom_dh) call c_f_procpointer(dh_ptr, dh_custom)

    allocate(x0(n), xi(n), x_new(n), dx(n))
    allocate(Hi(m*n), z_pred_i(m), innov(m))
    allocate(HT(n*m), PHT(n*m), S(mm), S_inv(mm), K(n*m))
    allocate(eye_n(nn), IKH(nn), IKHP(nn), IKHT(nn), IKHP_IKHT(nn))
    allocate(KR(n*m), KT(m*n), KRKT(nn))

    ! Save initial state as x0 (linearization reference)
    x0(1:n) = x(1:n)
    xi(1:n) = x(1:n)

    ! Iterative re-linearization loop
    do iter = 1, max_iter

        ! Evaluate h(x_i) at current iterate
        if (use_custom_h) then
            call h_custom(n, xi, m, 0.0d0, obs_params, nop, z_pred_i, obs_info)
        else
            call fa_observe_dispatch(obs_id, n, xi, m, 0.0d0, obs_params, nop, z_pred_i, obs_info)
        end if
        if (obs_info /= 0) then
            info = 4
            goto 900
        end if

        ! Evaluate H(x_i) at current iterate
        if (use_custom_dh) then
            call dh_custom(n, xi, m, 0.0d0, obs_params, nop, Hi, obs_info)
        else
            call fa_observe_jacobian(obs_id, n, xi, m, 0.0d0, obs_params, nop, Hi, obs_info)
        end if
        if (obs_info /= 0) then
            info = 4
            goto 900
        end if

        ! S = H_i * P * H_i^T + R
        call mat_transpose(Hi, HT, m, n)
        call mat_multiply(P, HT, PHT, n, n, m)
        call mat_multiply(Hi, PHT, S, m, n, m)
        call mat_add(S, R, S, m, m)

        ! K_i = P * H_i^T * S^-1
        call mat_identity(S_inv, m)
        call mat_solve_symm(S, S_inv, m, m, solve_info)
        if (solve_info /= 0) then
            info = 4
            goto 900
        end if
        call mat_multiply(PHT, S_inv, K, n, m, m)

        ! Innovation: z - h(x_i) - H_i * (x0 - x_i)
        ! First compute H_i * (x0 - x_i)
        dx(1:n) = x0(1:n) - xi(1:n)
        call mat_vec_multiply(Hi, dx, innov, m, n)
        ! innov = z - z_pred_i - H_i*(x0 - x_i)
        do i = 1, m
            innov(i) = z(i) - z_pred_i(i) - innov(i)
        end do

        ! x_new = x0 + K_i * innov
        call mat_vec_multiply(K, innov, dx, n, m)
        x_new(1:n) = x0(1:n) + dx(1:n)

        ! Check convergence
        dx(1:n) = x_new(1:n) - xi(1:n)
        diff_norm = vec_norm(dx, n)

        xi(1:n) = x_new(1:n)

        if (diff_norm < tol) exit
    end do

    ! If we exhausted iterations without converging, set info=4
    if (iter > max_iter) then
        info = 4
    end if

    ! Ternary gating on final innovation
    ! Recompute z_pred and H at final xi
    if (use_custom_h) then
        call h_custom(n, xi, m, 0.0d0, obs_params, nop, z_pred_i, obs_info)
    else
        call fa_observe_dispatch(obs_id, n, xi, m, 0.0d0, obs_params, nop, z_pred_i, obs_info)
    end if
    if (use_custom_dh) then
        call dh_custom(n, xi, m, 0.0d0, obs_params, nop, Hi, obs_info)
    else
        call fa_observe_jacobian(obs_id, n, xi, m, 0.0d0, obs_params, nop, Hi, obs_info)
    end if

    ! S_final = H * P * H^T + R
    call mat_transpose(Hi, HT, m, n)
    call mat_multiply(P, HT, PHT, n, n, m)
    call mat_multiply(Hi, PHT, S, m, n, m)
    call mat_add(S, R, S, m, m)

    ! Ternary gating
    validity(1:m) = 1
    do i = 1, m
        if (S((i-1)*m + i) > 0.0d0) then
            chi2 = (z(i) - z_pred_i(i))**2 / S((i-1)*m + i)
        else
            chi2 = 999.0d0
        end if
        if (chi2 > 16.0d0) then
            validity(i) = -1
        else if (chi2 > 9.0d0) then
            validity(i) = 0
        end if
    end do

    ! Update state
    x(1:n) = xi(1:n)

    ! Update P using final H: Joseph form P = (I-KH)*P*(I-KH)^T + K*R*K^T
    ! Recompute K with final H
    call mat_identity(S_inv, m)
    call mat_solve_symm(S, S_inv, m, m, solve_info)
    if (solve_info /= 0) then
        info = 4
        goto 900
    end if
    call mat_multiply(PHT, S_inv, K, n, m, m)

    ! IKH = I - K * H
    call mat_identity(eye_n, n)
    block
        real(c_double), allocatable :: KH(:)
        allocate(KH(nn))
        call mat_multiply(K, Hi, KH, n, m, n)
        IKH(1:nn) = eye_n(1:nn) - KH(1:nn)
        deallocate(KH)
    end block

    call mat_multiply(IKH, P, IKHP, n, n, n)
    call mat_transpose(IKH, IKHT, n, n)
    call mat_multiply(IKHP, IKHT, IKHP_IKHT, n, n, n)
    call mat_multiply(K, R, KR, n, m, m)
    call mat_transpose(K, KT, n, m)
    call mat_multiply(KR, KT, KRKT, n, m, n)
    call mat_add(IKHP_IKHT, KRKT, P, n, n)

900 continue
    deallocate(x0, xi, x_new, dx)
    deallocate(Hi, z_pred_i, innov)
    deallocate(HT, PHT, S, S_inv, K)
    deallocate(eye_n, IKH, IKHP, IKHT, IKHP_IKHT)
    deallocate(KR, KT, KRKT)

end subroutine fa_iekf_update


! ==========================================================================
! fa_ukf_predict — Unscented Kalman Filter predict step
!
! Generates 2n+1 sigma points, propagates through dynamics, reconstructs
! predicted mean and covariance.
! lambda = alpha^2*(n+kappa) - n
! Default UKF params if passed as 0: alpha=1e-3, beta=2.0, kappa=0.0
! ==========================================================================
subroutine fa_ukf_predict(n, x, P, f_ptr, model_id, params, np, Q, dt, &
                           alpha, beta_ukf, kappa, info) &
    bind(C, name="fa_ukf_predict")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, model_id, np
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: P(n*n)
    type(c_funptr), value              :: f_ptr
    real(c_double), intent(in)         :: params(np)
    real(c_double), intent(in)         :: Q(n*n)
    real(c_double), value, intent(in)  :: dt
    real(c_double), value, intent(in)  :: alpha, beta_ukf, kappa
    integer(c_int), intent(out)        :: info

    ! Interface for fa_propagate
    interface
        subroutine fa_propagate(n, x, u, nu, f_ptr, model_id, params, np, dt, n_steps, info) &
            bind(C, name="fa_propagate")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, nu, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(in)         :: u(nu), params(np)
            real(c_double), value, intent(in)  :: dt
            type(c_funptr), value              :: f_ptr
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double) :: a_eff, b_eff, k_eff, lam, scale_factor
    real(c_double) :: wm0, wc0, wmi
    integer :: n_sigma, nn, i, j, prop_info, chol_info
    real(c_double), allocatable :: L(:), Xi(:,:), Xi_prop(:,:)
    real(c_double), allocatable :: x_pred(:), diff(:), P_new(:), outer_tmp(:)
    real(c_double) :: u_dummy(1)

    info = 0
    nn = n * n

    if (n <= 0) then
        info = 3
        return
    end if

    if (dt == 0.0d0) return

    ! Default UKF parameters
    a_eff = alpha;  if (a_eff == 0.0d0) a_eff = 1.0d-3
    b_eff = beta_ukf; if (b_eff == 0.0d0) b_eff = 2.0d0
    k_eff = kappa  ! kappa=0 is valid default

    lam = a_eff * a_eff * (dble(n) + k_eff) - dble(n)
    scale_factor = dble(n) + lam
    n_sigma = 2 * n + 1

    ! Weights
    wm0 = lam / scale_factor
    wc0 = wm0 + (1.0d0 - a_eff*a_eff + b_eff)
    wmi = 1.0d0 / (2.0d0 * scale_factor)

    ! Cholesky of (n+lambda)*P
    allocate(L(nn))
    block
        real(c_double), allocatable :: scaled_P(:)
        allocate(scaled_P(nn))
        scaled_P(1:nn) = scale_factor * P(1:nn)
        call mat_cholesky(scaled_P, L, n, chol_info)
        deallocate(scaled_P)
    end block
    if (chol_info /= 0) then
        info = 2
        deallocate(L)
        return
    end if

    ! Generate sigma points (each column of Xi is a sigma point)
    allocate(Xi(n, n_sigma), Xi_prop(n, n_sigma))

    ! X_0 = x
    Xi(1:n, 1) = x(1:n)

    ! X_i = x + L_col_i, X_{n+i} = x - L_col_i
    ! L is row-major: column j of L = L((i-1)*n + j) for row i
    do i = 1, n
        do j = 1, n
            ! L column i: element at row j is L((j-1)*n + i)
            Xi(j, 1+i)   = x(j) + L((j-1)*n + i)
            Xi(j, 1+n+i) = x(j) - L((j-1)*n + i)
        end do
    end do

    ! Propagate each sigma point through dynamics
    u_dummy(1) = 0.0d0
    do i = 1, n_sigma
        Xi_prop(1:n, i) = Xi(1:n, i)
        call fa_propagate(n, Xi_prop(1:n, i), u_dummy, 0, f_ptr, model_id, &
                          params, np, dt, 10, prop_info)
        if (prop_info /= 0) then
            info = 4
            deallocate(L, Xi, Xi_prop)
            return
        end if
    end do

    ! Reconstruct predicted mean: x_pred = sum(W_mi * X_i)
    allocate(x_pred(n), diff(n), P_new(nn), outer_tmp(nn))
    x_pred(1:n) = wm0 * Xi_prop(1:n, 1)
    do i = 2, n_sigma
        x_pred(1:n) = x_pred(1:n) + wmi * Xi_prop(1:n, i)
    end do

    ! Reconstruct predicted covariance: P_pred = sum(W_ci * (Xi - x_pred)(Xi - x_pred)^T) + Q
    P_new(1:nn) = 0.0d0

    ! Sigma point 0
    diff(1:n) = Xi_prop(1:n, 1) - x_pred(1:n)
    call mat_outer(diff, diff, outer_tmp, n, n)
    call mat_scale(outer_tmp, wc0, nn)
    call mat_add(P_new, outer_tmp, P_new, n, n)

    ! Sigma points 1..2n
    do i = 2, n_sigma
        diff(1:n) = Xi_prop(1:n, i) - x_pred(1:n)
        call mat_outer(diff, diff, outer_tmp, n, n)
        call mat_scale(outer_tmp, wmi, nn)
        call mat_add(P_new, outer_tmp, P_new, n, n)
    end do

    ! Add process noise
    call mat_add(P_new, Q, P, n, n)
    x(1:n) = x_pred(1:n)

    deallocate(L, Xi, Xi_prop, x_pred, diff, P_new, outer_tmp)

end subroutine fa_ukf_predict


! ==========================================================================
! fa_ukf_update — Unscented Kalman Filter update step
!
! Generate sigma points, pass through measurement model, compute cross-
! covariance, Kalman gain, and update state/covariance.
! ==========================================================================
subroutine fa_ukf_update(n, m, x, P, z, h_ptr, obs_id, R, obs_params, nop, &
                          alpha, beta_ukf, kappa, validity, info) &
    bind(C, name="fa_ukf_update")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, m, obs_id, nop
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: P(n*n)
    real(c_double), intent(in)         :: z(m)
    type(c_funptr), value              :: h_ptr
    real(c_double), intent(in)         :: R(m*m)
    real(c_double), intent(in)         :: obs_params(nop)
    real(c_double), value, intent(in)  :: alpha, beta_ukf, kappa
    integer(c_int), intent(out)        :: validity(m)
    integer(c_int), intent(out)        :: info

    ! Interface for observation dispatch
    interface
        subroutine fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info) &
            bind(C, name="fa_observe_dispatch")
            use iso_c_binding
            integer(c_int), value  :: obs_id, n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    abstract interface
        subroutine custom_observe_ukf_t(n, x, m, t, obs_params, nop, z_pred, info) bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    procedure(custom_observe_ukf_t), pointer :: h_custom
    logical :: use_custom_h

    real(c_double) :: a_eff, b_eff, k_eff, lam, scale_factor
    real(c_double) :: wm0, wc0, wmi
    integer :: n_sigma, nn, mm, i, j, chol_info, obs_info, solve_info
    real(c_double), allocatable :: L(:)
    real(c_double), allocatable :: Xi(:,:), Zi(:,:)
    real(c_double), allocatable :: z_pred(:), diff_x(:), diff_z(:)
    real(c_double), allocatable :: Pzz(:), Pxz(:), Pzz_inv(:), K(:)
    real(c_double), allocatable :: outer_xz(:), outer_zz(:)
    real(c_double), allocatable :: KPzz(:), KT(:), KPzzKT(:)
    real(c_double) :: chi2

    info = 0
    nn = n * n
    mm = m * m

    if (n <= 0 .or. m <= 0) then
        info = 3
        return
    end if

    use_custom_h = c_associated(h_ptr)
    if (use_custom_h) call c_f_procpointer(h_ptr, h_custom)

    ! Default UKF parameters
    a_eff = alpha;  if (a_eff == 0.0d0) a_eff = 1.0d-3
    b_eff = beta_ukf; if (b_eff == 0.0d0) b_eff = 2.0d0
    k_eff = kappa

    lam = a_eff * a_eff * (dble(n) + k_eff) - dble(n)
    scale_factor = dble(n) + lam
    n_sigma = 2 * n + 1

    wm0 = lam / scale_factor
    wc0 = wm0 + (1.0d0 - a_eff*a_eff + b_eff)
    wmi = 1.0d0 / (2.0d0 * scale_factor)

    ! Cholesky of (n+lambda)*P
    allocate(L(nn))
    block
        real(c_double), allocatable :: scaled_P(:)
        allocate(scaled_P(nn))
        scaled_P(1:nn) = scale_factor * P(1:nn)
        call mat_cholesky(scaled_P, L, n, chol_info)
        deallocate(scaled_P)
    end block
    if (chol_info /= 0) then
        info = 2
        deallocate(L)
        return
    end if

    ! Generate sigma points
    allocate(Xi(n, n_sigma), Zi(m, n_sigma))
    Xi(1:n, 1) = x(1:n)
    do i = 1, n
        do j = 1, n
            Xi(j, 1+i)   = x(j) + L((j-1)*n + i)
            Xi(j, 1+n+i) = x(j) - L((j-1)*n + i)
        end do
    end do

    ! Pass each sigma point through measurement model
    do i = 1, n_sigma
        if (use_custom_h) then
            call h_custom(n, Xi(1:n, i), m, 0.0d0, obs_params, nop, Zi(1:m, i), obs_info)
        else
            call fa_observe_dispatch(obs_id, n, Xi(1:n, i), m, 0.0d0, obs_params, nop, Zi(1:m, i), obs_info)
        end if
        if (obs_info /= 0) then
            info = 4
            deallocate(L, Xi, Zi)
            return
        end if
    end do

    ! Predicted measurement mean: z_pred = sum(W_mi * Z_i)
    allocate(z_pred(m), diff_x(n), diff_z(m))
    allocate(Pzz(mm), Pxz(n*m), Pzz_inv(mm), K(n*m))
    allocate(outer_xz(n*m), outer_zz(mm))
    allocate(KPzz(n*m), KT(m*n), KPzzKT(nn))

    z_pred(1:m) = wm0 * Zi(1:m, 1)
    do i = 2, n_sigma
        z_pred(1:m) = z_pred(1:m) + wmi * Zi(1:m, i)
    end do

    ! Pzz = sum(W_ci * (Z_i - z_pred)(Z_i - z_pred)^T) + R
    Pzz(1:mm) = 0.0d0
    Pxz(1:n*m) = 0.0d0

    ! Sigma point 0
    diff_z(1:m) = Zi(1:m, 1) - z_pred(1:m)
    diff_x(1:n) = Xi(1:n, 1) - x(1:n)
    call mat_outer(diff_z, diff_z, outer_zz, m, m)
    call mat_scale(outer_zz, wc0, mm)
    call mat_add(Pzz, outer_zz, Pzz, m, m)
    call mat_outer(diff_x, diff_z, outer_xz, n, m)
    call mat_scale(outer_xz, wc0, n*m)
    call mat_add(Pxz, outer_xz, Pxz, n, m)

    ! Sigma points 1..2n
    do i = 2, n_sigma
        diff_z(1:m) = Zi(1:m, i) - z_pred(1:m)
        diff_x(1:n) = Xi(1:n, i) - x(1:n)
        call mat_outer(diff_z, diff_z, outer_zz, m, m)
        call mat_scale(outer_zz, wmi, mm)
        call mat_add(Pzz, outer_zz, Pzz, m, m)
        call mat_outer(diff_x, diff_z, outer_xz, n, m)
        call mat_scale(outer_xz, wmi, n*m)
        call mat_add(Pxz, outer_xz, Pxz, n, m)
    end do

    ! Add measurement noise
    call mat_add(Pzz, R, Pzz, m, m)

    ! Ternary gating on innovation
    validity(1:m) = 1
    do i = 1, m
        if (Pzz((i-1)*m + i) > 0.0d0) then
            chi2 = (z(i) - z_pred(i))**2 / Pzz((i-1)*m + i)
        else
            chi2 = 999.0d0
        end if
        if (chi2 > 16.0d0) then
            validity(i) = -1
        else if (chi2 > 9.0d0) then
            validity(i) = 0
        end if
    end do

    ! K = Pxz * Pzz^-1
    Pzz_inv(1:mm) = 0.0d0
    call mat_identity(Pzz_inv, m)
    call mat_solve_symm(Pzz, Pzz_inv, m, m, solve_info)
    if (solve_info /= 0) then
        info = 4
        deallocate(L, Xi, Zi, z_pred, diff_x, diff_z)
        deallocate(Pzz, Pxz, Pzz_inv, K, outer_xz, outer_zz)
        deallocate(KPzz, KT, KPzzKT)
        return
    end if
    call mat_multiply(Pxz, Pzz_inv, K, n, m, m)

    ! x += K * (z - z_pred)
    diff_z(1:m) = z(1:m) - z_pred(1:m)
    block
        real(c_double) :: dx(n)
        call mat_vec_multiply(K, diff_z, dx, n, m)
        x(1:n) = x(1:n) + dx(1:n)
    end block

    ! P -= K * Pzz * K^T
    ! Recompute Pzz (it was modified by solve) — use Pxz*K^T shortcut: P = P - K * Pzz_orig * K^T
    ! Actually Pzz was overwritten by solve. Use: P = P - Pxz * K^T (equivalent since K = Pxz * Pzz^-1)
    ! P_new = P - K * Pzz * K^T. Since Pzz was modified, compute P = P - Pxz * Pzz_inv^T * Pxz^T
    ! Simplest correct form: P = P - K * S * K^T where S = Pzz before solve
    ! But Pzz was modified. Use P = P - Pxz * K^T (since K*Pzz = Pxz, so Pxz*K^T = K*Pzz*K^T when Pzz is symmetric)
    ! Actually: K = Pxz * Pzz^-1, so K * Pzz = Pxz. Then K * Pzz * K^T = Pxz * K^T. Correct.
    call mat_transpose(K, KT, n, m)
    call mat_multiply(Pxz, KT, KPzzKT, n, m, n)
    call mat_subtract(P, KPzzKT, P, n, n)

    deallocate(L, Xi, Zi, z_pred, diff_x, diff_z)
    deallocate(Pzz, Pxz, Pzz_inv, K, outer_xz, outer_zz)
    deallocate(KPzz, KT, KPzzKT)

end subroutine fa_ukf_update


! ==========================================================================
! fa_eskf_predict — Error-State Kalman Filter predict
!
! Apollo heritage: propagates nominal state via dynamics, error state stays
! zero, covariance propagated via STM: P = Phi * P * Phi^T + Q
! ==========================================================================
subroutine fa_eskf_predict(n, x_nom, dx, P, f_ptr, model_id, params, np, Q, dt, n_steps, info) &
    bind(C, name="fa_eskf_predict")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, model_id, np, n_steps
    real(c_double), intent(inout)      :: x_nom(n)
    real(c_double), intent(inout)      :: dx(n)
    real(c_double), intent(inout)      :: P(n*n)
    type(c_funptr), value              :: f_ptr
    real(c_double), intent(in)         :: params(np)
    real(c_double), intent(in)         :: Q(n*n)
    real(c_double), value, intent(in)  :: dt
    integer(c_int), intent(out)        :: info

    ! Interface for fa_propagate and fa_propagate_stm
    interface
        subroutine fa_propagate(n, x, u, nu, f_ptr, model_id, params, np, dt, n_steps, info) &
            bind(C, name="fa_propagate")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, nu, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(in)         :: u(nu), params(np)
            real(c_double), value, intent(in)  :: dt
            type(c_funptr), value              :: f_ptr
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_propagate_stm(n, x, phi, u, nu, f_ptr, df_ptr, model_id, params, np, dt, n_steps, info) &
            bind(C, name="fa_propagate_stm")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, nu, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n), phi(n*n)
            real(c_double), intent(in)         :: u(nu), params(np)
            real(c_double), value, intent(in)  :: dt
            type(c_funptr), value              :: f_ptr, df_ptr
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double), allocatable :: Phi(:), PhiP(:), PhiT(:), PhiPPhiT(:)
    real(c_double) :: u_dummy(1)
    integer :: nn, prop_info

    info = 0
    nn = n * n

    if (n <= 0) then
        info = 3
        return
    end if

    if (dt == 0.0d0) return

    u_dummy(1) = 0.0d0

    allocate(Phi(nn), PhiP(nn), PhiT(nn), PhiPPhiT(nn))

    ! Initialize Phi to identity
    call mat_identity(Phi, n)

    ! Propagate nominal state and get STM
    ! We need both state propagation and STM. Use fa_propagate_stm with a copy.
    block
        real(c_double), allocatable :: x_copy(:)
        allocate(x_copy(n))
        x_copy(1:n) = x_nom(1:n)
        call fa_propagate_stm(n, x_copy, Phi, u_dummy, 0, f_ptr, c_null_funptr, model_id, &
                              params, np, dt, n_steps, prop_info)
        x_nom(1:n) = x_copy(1:n)
        deallocate(x_copy)
    end block
    if (prop_info /= 0) then
        info = 4
        deallocate(Phi, PhiP, PhiT, PhiPPhiT)
        return
    end if

    ! Error state stays zero after prediction (reset after injection)
    dx(1:n) = 0.0d0

    ! Propagate covariance: P = Phi * P * Phi^T + Q
    call mat_multiply(Phi, P, PhiP, n, n, n)
    call mat_transpose(Phi, PhiT, n, n)
    call mat_multiply(PhiP, PhiT, PhiPPhiT, n, n, n)
    call mat_add(PhiPPhiT, Q, P, n, n)

    deallocate(Phi, PhiP, PhiT, PhiPPhiT)

end subroutine fa_eskf_predict


! ==========================================================================
! fa_eskf_update — Error-State Kalman Filter update
!
! Innovation: y = z - h(x_nom)
! H evaluated at x_nom
! Standard Kalman update on error state dx: dx = K * y, P = (I-KH)*P
! ==========================================================================
subroutine fa_eskf_update(n, m, x_nom, dx, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, &
                           validity, info) &
    bind(C, name="fa_eskf_update")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, m, obs_id, nop
    real(c_double), intent(in)         :: x_nom(n)
    real(c_double), intent(inout)      :: dx(n)
    real(c_double), intent(inout)      :: P(n*n)
    real(c_double), intent(in)         :: z(m)
    type(c_funptr), value              :: h_ptr
    type(c_funptr), value              :: dh_ptr
    real(c_double), intent(in)         :: R(m*m)
    real(c_double), intent(in)         :: obs_params(nop)
    integer(c_int), intent(out)        :: validity(m)
    integer(c_int), intent(out)        :: info

    ! Interfaces for observation dispatch
    interface
        subroutine fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info) &
            bind(C, name="fa_observe_dispatch")
            use iso_c_binding
            integer(c_int), value  :: obs_id, n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine

        subroutine fa_observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H, info) &
            bind(C, name="fa_observe_jacobian")
            use iso_c_binding
            integer(c_int), value  :: obs_id, n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: H(m*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    abstract interface
        subroutine custom_observe_eskf_t(n, x, m, t, obs_params, nop, z_pred, info) bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: z_pred(m)
            integer(c_int), intent(out) :: info
        end subroutine

        subroutine custom_obs_jac_eskf_t(n, x, m, t, obs_params, nop, H, info) bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, m, nop
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: obs_params(nop)
            real(c_double), intent(out) :: H(m*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    procedure(custom_observe_eskf_t), pointer     :: h_custom
    procedure(custom_obs_jac_eskf_t), pointer      :: dh_custom
    logical :: use_custom_h, use_custom_dh

    real(c_double), allocatable :: z_pred(:), y(:), H(:)
    real(c_double), allocatable :: HT(:), PHT(:), S(:), S_inv(:), K(:)
    real(c_double), allocatable :: eye_n(:), IKH(:), IKHP(:), IKHT(:), IKHP_IKHT(:)
    real(c_double), allocatable :: KR(:), KT_mat(:), KRKT(:)
    real(c_double) :: chi2
    integer :: i, nn, mm, obs_info, solve_info

    info = 0
    nn = n * n
    mm = m * m

    if (n <= 0 .or. m <= 0) then
        info = 3
        return
    end if

    use_custom_h  = c_associated(h_ptr)
    use_custom_dh = c_associated(dh_ptr)
    if (use_custom_h)  call c_f_procpointer(h_ptr, h_custom)
    if (use_custom_dh) call c_f_procpointer(dh_ptr, dh_custom)

    allocate(z_pred(m), y(m), H(m*n))
    allocate(HT(n*m), PHT(n*m), S(mm), S_inv(mm), K(n*m))
    allocate(eye_n(nn), IKH(nn), IKHP(nn), IKHT(nn), IKHP_IKHT(nn))
    allocate(KR(n*m), KT_mat(m*n), KRKT(nn))

    ! Compute z_pred = h(x_nom)
    if (use_custom_h) then
        call h_custom(n, x_nom, m, 0.0d0, obs_params, nop, z_pred, obs_info)
    else
        call fa_observe_dispatch(obs_id, n, x_nom, m, 0.0d0, obs_params, nop, z_pred, obs_info)
    end if
    if (obs_info /= 0) then
        info = 4
        goto 910
    end if

    ! Compute H at x_nom
    if (use_custom_dh) then
        call dh_custom(n, x_nom, m, 0.0d0, obs_params, nop, H, obs_info)
    else
        call fa_observe_jacobian(obs_id, n, x_nom, m, 0.0d0, obs_params, nop, H, obs_info)
    end if
    if (obs_info /= 0) then
        info = 4
        goto 910
    end if

    ! Innovation: y = z - h(x_nom)
    y(1:m) = z(1:m) - z_pred(1:m)

    ! S = H * P * H^T + R
    call mat_transpose(H, HT, m, n)
    call mat_multiply(P, HT, PHT, n, n, m)
    call mat_multiply(H, PHT, S, m, n, m)
    call mat_add(S, R, S, m, m)

    ! Ternary gating
    validity(1:m) = 1
    do i = 1, m
        if (S((i-1)*m + i) > 0.0d0) then
            chi2 = y(i) * y(i) / S((i-1)*m + i)
        else
            chi2 = 999.0d0
        end if
        if (chi2 > 16.0d0) then
            validity(i) = -1
        else if (chi2 > 9.0d0) then
            validity(i) = 0
        end if
    end do

    ! K = P * H^T * S^-1
    call mat_identity(S_inv, m)
    call mat_solve_symm(S, S_inv, m, m, solve_info)
    if (solve_info /= 0) then
        info = 4
        goto 910
    end if
    call mat_multiply(PHT, S_inv, K, n, m, m)

    ! dx = K * y (update error state)
    call mat_vec_multiply(K, y, dx, n, m)

    ! Joseph form: P = (I-KH)*P*(I-KH)^T + K*R*K^T
    call mat_identity(eye_n, n)
    block
        real(c_double), allocatable :: KH(:)
        allocate(KH(nn))
        call mat_multiply(K, H, KH, n, m, n)
        IKH(1:nn) = eye_n(1:nn) - KH(1:nn)
        deallocate(KH)
    end block

    call mat_multiply(IKH, P, IKHP, n, n, n)
    call mat_transpose(IKH, IKHT, n, n)
    call mat_multiply(IKHP, IKHT, IKHP_IKHT, n, n, n)
    call mat_multiply(K, R, KR, n, m, m)
    call mat_transpose(K, KT_mat, n, m)
    call mat_multiply(KR, KT_mat, KRKT, n, m, n)
    call mat_add(IKHP_IKHT, KRKT, P, n, n)

910 continue
    deallocate(z_pred, y, H)
    deallocate(HT, PHT, S, S_inv, K)
    deallocate(eye_n, IKH, IKHP, IKHT, IKHP_IKHT)
    deallocate(KR, KT_mat, KRKT)

end subroutine fa_eskf_update


! ==========================================================================
! fa_eskf_inject — Inject error state into nominal state
!
! x_nom = x_nom + dx
! dx = 0 (reset error state after injection)
! P stays unchanged (already reflects uncertainty)
! ==========================================================================
subroutine fa_eskf_inject(n, x_nom, dx, P, info) bind(C, name="fa_eskf_inject")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n
    real(c_double), intent(inout)      :: x_nom(n)
    real(c_double), intent(inout)      :: dx(n)
    real(c_double), intent(inout)      :: P(n*n)
    integer(c_int), intent(out)        :: info

    integer :: i

    info = 0

    if (n <= 0) then
        info = 3
        return
    end if

    ! Inject error into nominal
    ! For general states, additive injection: x_nom = x_nom + dx
    ! TODO: For quaternion states, use multiplicative correction:
    !   q_nom = delta_q(dx(quat_indices)) * q_nom
    do i = 1, n
        x_nom(i) = x_nom(i) + dx(i)
    end do

    ! Reset error state
    dx(1:n) = 0.0d0

    ! P stays — it already reflects the updated uncertainty

end subroutine fa_eskf_inject


! ==========================================================================
! fa_srekf_predict — Square-Root EKF predict step
!
! S(n*n) — lower Cholesky factor of P (P = S*S^T), flat row-major
! Sq(n*n) — lower Cholesky factor of Q
!
! Simplified implementation: convert S→P, run regular EKF predict, convert P→S
! TODO: Proper QR-based implementation: form [Phi*S, Sq], QR factorize
!       to get new S directly without forming full P. This avoids the
!       numerical loss of squaring then taking sqrt.
! ==========================================================================
subroutine fa_srekf_predict(n, x, S, f_ptr, df_ptr, model_id, params, np, Sq, dt, n_steps, info) &
    bind(C, name="fa_srekf_predict")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, model_id, np, n_steps
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: S(n*n)
    type(c_funptr), value              :: f_ptr
    type(c_funptr), value              :: df_ptr
    real(c_double), intent(in)         :: Sq(n*n)
    real(c_double), value, intent(in)  :: dt
    real(c_double), intent(in)         :: params(np)
    integer(c_int), intent(out)        :: info

    ! Interface for fa_ekf_predict
    interface
        subroutine fa_ekf_predict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, n_steps, info) &
            bind(C, name="fa_ekf_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            type(c_funptr), value              :: f_ptr, df_ptr
            real(c_double), intent(in)         :: Q(n*n)
            real(c_double), value, intent(in)  :: dt
            real(c_double), intent(in)         :: params(np)
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double), allocatable :: P(:), Q(:), ST(:), SqT(:)
    integer :: nn, chol_info

    info = 0
    nn = n * n

    if (n <= 0) then
        info = 3
        return
    end if

    allocate(P(nn), Q(nn), ST(nn), SqT(nn))

    ! Convert S → P = S * S^T
    call mat_transpose(S, ST, n, n)
    call mat_multiply(S, ST, P, n, n, n)

    ! Convert Sq → Q = Sq * Sq^T
    call mat_transpose(Sq, SqT, n, n)
    call mat_multiply(Sq, SqT, Q, n, n, n)

    ! Run regular EKF predict
    call fa_ekf_predict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, n_steps, info)

    if (info /= 0) then
        deallocate(P, Q, ST, SqT)
        return
    end if

    ! Convert P → S via Cholesky
    call mat_cholesky(P, S, n, chol_info)
    if (chol_info /= 0) then
        info = 2
    end if

    deallocate(P, Q, ST, SqT)

end subroutine fa_srekf_predict


! ==========================================================================
! fa_srekf_update — Square-Root EKF update step
!
! Simplified: convert S→P, run regular EKF update, Cholesky P→S.
! TODO: Proper rank-1 Cholesky update/downdate for efficiency.
! ==========================================================================
subroutine fa_srekf_update(n, m, x, S, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info) &
    bind(C, name="fa_srekf_update")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, m, obs_id, nop
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: S(n*n)
    real(c_double), intent(in)         :: z(m)
    type(c_funptr), value              :: h_ptr
    type(c_funptr), value              :: dh_ptr
    real(c_double), intent(in)         :: R(m*m)
    real(c_double), intent(in)         :: obs_params(nop)
    integer(c_int), intent(out)        :: validity(m)
    integer(c_int), intent(out)        :: info

    ! Interface for fa_ekf_update
    interface
        subroutine fa_ekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info) &
            bind(C, name="fa_ekf_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m, obs_id, nop
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            real(c_double), intent(in)         :: z(m)
            type(c_funptr), value              :: h_ptr, dh_ptr
            real(c_double), intent(in)         :: R(m*m)
            real(c_double), intent(in)         :: obs_params(nop)
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double), allocatable :: P(:), ST(:)
    integer :: nn, chol_info

    info = 0
    nn = n * n

    if (n <= 0 .or. m <= 0) then
        info = 3
        return
    end if

    allocate(P(nn), ST(nn))

    ! Convert S → P = S * S^T
    call mat_transpose(S, ST, n, n)
    call mat_multiply(S, ST, P, n, n, n)

    ! Run regular EKF update
    call fa_ekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info)

    if (info /= 0) then
        deallocate(P, ST)
        return
    end if

    ! Convert P → S via Cholesky
    call mat_cholesky(P, S, n, chol_info)
    if (chol_info /= 0) then
        info = 2
    end if

    deallocate(P, ST)

end subroutine fa_srekf_update


! ==========================================================================
! fa_srukf_predict — Square-Root UKF predict step
!
! Simplified: convert S→P, run regular UKF predict, Cholesky P→S.
! TODO: Proper SR-UKF using QR decomposition on sigma point residuals
!       for direct square-root factor propagation.
! ==========================================================================
subroutine fa_srukf_predict(n, x, S, f_ptr, model_id, params, np, Sq, dt, &
                             alpha, beta_ukf, kappa, info) &
    bind(C, name="fa_srukf_predict")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, model_id, np
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: S(n*n)
    type(c_funptr), value              :: f_ptr
    real(c_double), intent(in)         :: params(np)
    real(c_double), intent(in)         :: Sq(n*n)
    real(c_double), value, intent(in)  :: dt
    real(c_double), value, intent(in)  :: alpha, beta_ukf, kappa
    integer(c_int), intent(out)        :: info

    ! Interface for fa_ukf_predict
    interface
        subroutine fa_ukf_predict(n, x, P, f_ptr, model_id, params, np, Q, dt, &
                                   alpha, beta_ukf, kappa, info) &
            bind(C, name="fa_ukf_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, model_id, np
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            type(c_funptr), value              :: f_ptr
            real(c_double), intent(in)         :: params(np)
            real(c_double), intent(in)         :: Q(n*n)
            real(c_double), value, intent(in)  :: dt
            real(c_double), value, intent(in)  :: alpha, beta_ukf, kappa
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double), allocatable :: P(:), Q(:), ST(:), SqT(:)
    integer :: nn, chol_info

    info = 0
    nn = n * n

    if (n <= 0) then
        info = 3
        return
    end if

    allocate(P(nn), Q(nn), ST(nn), SqT(nn))

    ! Convert S → P = S * S^T
    call mat_transpose(S, ST, n, n)
    call mat_multiply(S, ST, P, n, n, n)

    ! Convert Sq → Q = Sq * Sq^T
    call mat_transpose(Sq, SqT, n, n)
    call mat_multiply(Sq, SqT, Q, n, n, n)

    ! Run regular UKF predict
    call fa_ukf_predict(n, x, P, f_ptr, model_id, params, np, Q, dt, alpha, beta_ukf, kappa, info)

    if (info /= 0) then
        deallocate(P, Q, ST, SqT)
        return
    end if

    ! Convert P → S via Cholesky
    call mat_cholesky(P, S, n, chol_info)
    if (chol_info /= 0) then
        info = 2
    end if

    deallocate(P, Q, ST, SqT)

end subroutine fa_srukf_predict


! ==========================================================================
! fa_srukf_update — Square-Root UKF update step
!
! Simplified: convert S→P, run regular UKF update, Cholesky P→S.
! TODO: Proper SR-UKF update using QR decomposition on measurement
!       sigma point residuals for direct square-root factor update.
! ==========================================================================
subroutine fa_srukf_update(n, m, x, S, z, h_ptr, obs_id, R, obs_params, nop, &
                            alpha, beta_ukf, kappa, validity, info) &
    bind(C, name="fa_srukf_update")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, m, obs_id, nop
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(inout)      :: S(n*n)
    real(c_double), intent(in)         :: z(m)
    type(c_funptr), value              :: h_ptr
    real(c_double), intent(in)         :: R(m*m)
    real(c_double), intent(in)         :: obs_params(nop)
    real(c_double), value, intent(in)  :: alpha, beta_ukf, kappa
    integer(c_int), intent(out)        :: validity(m)
    integer(c_int), intent(out)        :: info

    ! Interface for fa_ukf_update
    interface
        subroutine fa_ukf_update(n, m, x, P, z, h_ptr, obs_id, R, obs_params, nop, &
                                  alpha, beta_ukf, kappa, validity, info) &
            bind(C, name="fa_ukf_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m, obs_id, nop
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            real(c_double), intent(in)         :: z(m)
            type(c_funptr), value              :: h_ptr
            real(c_double), intent(in)         :: R(m*m)
            real(c_double), intent(in)         :: obs_params(nop)
            real(c_double), value, intent(in)  :: alpha, beta_ukf, kappa
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine
    end interface

    real(c_double), allocatable :: P(:), ST(:)
    integer :: nn, chol_info

    info = 0
    nn = n * n

    if (n <= 0 .or. m <= 0) then
        info = 3
        return
    end if

    allocate(P(nn), ST(nn))

    ! Convert S → P = S * S^T
    call mat_transpose(S, ST, n, n)
    call mat_multiply(S, ST, P, n, n, n)

    ! Run regular UKF update
    call fa_ukf_update(n, m, x, P, z, h_ptr, obs_id, R, obs_params, nop, alpha, beta_ukf, kappa, validity, info)

    if (info /= 0) then
        deallocate(P, ST)
        return
    end if

    ! Convert P → S via Cholesky
    call mat_cholesky(P, S, n, chol_info)
    if (chol_info /= 0) then
        info = 2
    end if

    deallocate(P, ST)

end subroutine fa_srukf_update
