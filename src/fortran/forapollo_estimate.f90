! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_estimate.f90 — Linear Kalman Filter and Extended Kalman Filter
!
! Core state estimation routines with ternary measurement gating.
! All matrices are flat 1D row-major at the bind(C) boundary.
!
! Entry points:
!   fa_kf_predict   — Linear Kalman predict: x = F*x, P = F*P*F^T + Q
!   fa_kf_update    — Linear Kalman update with ternary gating
!   fa_ekf_predict  — Extended Kalman predict via propagator + STM
!   fa_ekf_update   — Extended Kalman update with ternary gating
!
! Ternary measurement gating:
!   For each measurement i, compute chi2 = y(i)^2 / S(i,i)
!     chi2 > 16.0 (4-sigma) → validity(i) = -1, rejected
!     chi2 > 9.0  (3-sigma) → validity(i) =  0, uncertain (skip)
!     else                   → validity(i) = +1, fused
!
! Internal matrix helpers (no bind(C), no external deps):
!   mat_multiply, mat_transpose, mat_add, mat_identity, mat_solve_symm
!
! Joseph form for P update: P = (I-KH)*P*(I-KH)^T + K*R*K^T
! This is numerically superior to simple (I-KH)*P.
!
! Error codes: info = 0 ok, 3 invalid input, 4 propagation/observation failure.
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
