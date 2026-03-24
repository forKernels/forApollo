! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_propagate.f90 — State propagation via RK4 integration
!
! Three bind(C) entry points:
!   fa_propagate       — single state, fixed-step RK4
!   fa_propagate_stm   — state + state transition matrix (augmented RK4)
!   fa_propagate_batch — OpenMP parallel batch propagation
!
! The propagator calls fa_dynamics_dispatch / fa_dynamics_jacobian for
! built-in model IDs, or user-supplied C function pointers for custom
! dynamics. This keeps the propagation engine model-agnostic.
!
! RK4 (classical 4th-order Runge-Kutta) with fixed sub-stepping:
!   h = dt / n_steps
!   For each sub-step:
!     k1 = h * f(t, x)
!     k2 = h * f(t + h/2, x + k1/2)
!     k3 = h * f(t + h/2, x + k2/2)
!     k4 = h * f(t + h, x + k3)
!     x  = x + (k1 + 2*k2 + 2*k3 + k4) / 6
!     t  = t + h
!
! Error codes: info = 0 ok, 3 invalid input, 4 dynamics failure.
!
! Matrices are flat row-major: element (row, col) = M((row-1)*n + col).

! ==========================================================================
! fa_propagate — single state RK4 propagation
! ==========================================================================
subroutine fa_propagate(n, x, u, nu, f_ptr, model_id, params, np, dt, n_steps, info) &
    bind(C, name="fa_propagate")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, nu, model_id, np, n_steps
    real(c_double), intent(inout)      :: x(n)
    real(c_double), intent(in)         :: u(nu), params(np)
    real(c_double), value, intent(in)  :: dt
    type(c_funptr), value              :: f_ptr
    integer(c_int), intent(out)        :: info

    ! Interface for fa_dynamics_dispatch (external bind(C) routine)
    interface
        subroutine fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info) &
            bind(C, name="fa_dynamics_dispatch")
            use iso_c_binding
            integer(c_int), value  :: model_id, n, nu, np
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: u(nu)
            real(c_double), intent(in)  :: params(np)
            real(c_double), intent(out) :: x_dot(n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    ! Interface for user-supplied dynamics function pointer
    abstract interface
        subroutine custom_dynamics_t(n, x, u, nu, t, params, np, x_dot, info) &
            bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, nu, np
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: u(nu)
            real(c_double), intent(in)  :: params(np)
            real(c_double), intent(out) :: x_dot(n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    procedure(custom_dynamics_t), pointer :: f_custom
    real(c_double) :: h, t
    real(c_double), allocatable :: k1(:), k2(:), k3(:), k4(:), x_tmp(:)
    integer :: steps_actual, step, dyn_info
    logical :: use_custom

    info = 0

    ! --- Input validation ---
    if (n <= 0) then
        info = 3
        return
    end if

    ! dt = 0 means no propagation
    if (dt == 0.0d0) return

    ! Clamp n_steps
    if (n_steps <= 0) then
        steps_actual = 1
    else
        steps_actual = n_steps
    end if

    h = dt / dble(steps_actual)
    t = 0.0d0

    ! Determine whether to use custom function pointer
    use_custom = c_associated(f_ptr)
    if (use_custom) then
        call c_f_procpointer(f_ptr, f_custom)
    end if

    allocate(k1(n), k2(n), k3(n), k4(n), x_tmp(n))

    ! --- RK4 sub-stepping loop ---
    do step = 1, steps_actual

        ! k1 = h * f(t, x)
        if (use_custom) then
            call f_custom(n, x, u, nu, t, params, np, k1, dyn_info)
        else
            call fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, k1, dyn_info)
        end if
        if (dyn_info /= 0) then
            info = 4
            deallocate(k1, k2, k3, k4, x_tmp)
            return
        end if
        k1 = h * k1

        ! k2 = h * f(t + h/2, x + k1/2)
        x_tmp = x + 0.5d0 * k1
        if (use_custom) then
            call f_custom(n, x_tmp, u, nu, t + 0.5d0*h, params, np, k2, dyn_info)
        else
            call fa_dynamics_dispatch(model_id, n, x_tmp, u, nu, t + 0.5d0*h, params, np, k2, dyn_info)
        end if
        if (dyn_info /= 0) then
            info = 4
            deallocate(k1, k2, k3, k4, x_tmp)
            return
        end if
        k2 = h * k2

        ! k3 = h * f(t + h/2, x + k2/2)
        x_tmp = x + 0.5d0 * k2
        if (use_custom) then
            call f_custom(n, x_tmp, u, nu, t + 0.5d0*h, params, np, k3, dyn_info)
        else
            call fa_dynamics_dispatch(model_id, n, x_tmp, u, nu, t + 0.5d0*h, params, np, k3, dyn_info)
        end if
        if (dyn_info /= 0) then
            info = 4
            deallocate(k1, k2, k3, k4, x_tmp)
            return
        end if
        k3 = h * k3

        ! k4 = h * f(t + h, x + k3)
        x_tmp = x + k3
        if (use_custom) then
            call f_custom(n, x_tmp, u, nu, t + h, params, np, k4, dyn_info)
        else
            call fa_dynamics_dispatch(model_id, n, x_tmp, u, nu, t + h, params, np, k4, dyn_info)
        end if
        if (dyn_info /= 0) then
            info = 4
            deallocate(k1, k2, k3, k4, x_tmp)
            return
        end if
        k4 = h * k4

        ! Update state
        x = x + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4) / 6.0d0
        t = t + h

    end do

    deallocate(k1, k2, k3, k4, x_tmp)

end subroutine fa_propagate


! ==========================================================================
! fa_propagate_stm — state + state transition matrix propagation
!
! Propagates both x(n) and Phi(n*n) simultaneously using augmented RK4.
! Phi_dot = F(x) * Phi, where F is the Jacobian df/dx.
! Phi should be initialized to identity before calling.
! ==========================================================================
subroutine fa_propagate_stm(n, x, phi, u, nu, f_ptr, df_ptr, model_id, params, np, dt, n_steps, info) &
    bind(C, name="fa_propagate_stm")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, nu, model_id, np, n_steps
    real(c_double), intent(inout)      :: x(n), phi(n*n)
    real(c_double), intent(in)         :: u(nu), params(np)
    real(c_double), value, intent(in)  :: dt
    type(c_funptr), value              :: f_ptr    ! custom dynamics (null → model_id)
    type(c_funptr), value              :: df_ptr   ! custom Jacobian (null → model_id)
    integer(c_int), intent(out)        :: info

    ! Interface for fa_dynamics_dispatch
    interface
        subroutine fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info) &
            bind(C, name="fa_dynamics_dispatch")
            use iso_c_binding
            integer(c_int), value  :: model_id, n, nu, np
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: u(nu)
            real(c_double), intent(in)  :: params(np)
            real(c_double), intent(out) :: x_dot(n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    ! Interface for fa_dynamics_jacobian
    interface
        subroutine fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np, F, info) &
            bind(C, name="fa_dynamics_jacobian")
            use iso_c_binding
            integer(c_int), value  :: model_id, n, nu, np
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: u(nu)
            real(c_double), intent(in)  :: params(np)
            real(c_double), intent(out) :: F(n*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    ! Abstract interfaces for user-supplied function pointers
    abstract interface
        subroutine custom_dynamics_t(n, x, u, nu, t, params, np, x_dot, info) &
            bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, nu, np
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: u(nu)
            real(c_double), intent(in)  :: params(np)
            real(c_double), intent(out) :: x_dot(n)
            integer(c_int), intent(out) :: info
        end subroutine

        subroutine custom_jacobian_t(n, x, u, nu, t, params, np, F, info) &
            bind(C)
            use iso_c_binding
            integer(c_int), value  :: n, nu, np
            real(c_double), value  :: t
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(in)  :: u(nu)
            real(c_double), intent(in)  :: params(np)
            real(c_double), intent(out) :: F(n*n)
            integer(c_int), intent(out) :: info
        end subroutine
    end interface

    procedure(custom_dynamics_t), pointer :: f_custom
    procedure(custom_jacobian_t), pointer :: df_custom
    logical :: use_custom_f, use_custom_df
    real(c_double) :: h, t
    integer :: steps_actual, step, dyn_info, nn
    integer :: i, j, k_idx

    ! RK4 stages for state x
    real(c_double), allocatable :: kx1(:), kx2(:), kx3(:), kx4(:), x_tmp(:)
    ! RK4 stages for Phi (flat n*n)
    real(c_double), allocatable :: kp1(:), kp2(:), kp3(:), kp4(:), phi_tmp(:)
    ! Jacobian workspace
    real(c_double), allocatable :: F_jac(:), phi_dot_tmp(:)

    info = 0
    nn = n * n

    ! --- Input validation ---
    if (n <= 0) then
        info = 3
        return
    end if

    if (dt == 0.0d0) return

    if (n_steps <= 0) then
        steps_actual = 1
    else
        steps_actual = n_steps
    end if

    h = dt / dble(steps_actual)
    t = 0.0d0

    ! Resolve function pointers
    use_custom_f  = c_associated(f_ptr)
    use_custom_df = c_associated(df_ptr)
    if (use_custom_f)  call c_f_procpointer(f_ptr, f_custom)
    if (use_custom_df) call c_f_procpointer(df_ptr, df_custom)

    allocate(kx1(n), kx2(n), kx3(n), kx4(n), x_tmp(n))
    allocate(kp1(nn), kp2(nn), kp3(nn), kp4(nn), phi_tmp(nn))
    allocate(F_jac(nn), phi_dot_tmp(nn))

    ! --- RK4 sub-stepping loop (augmented: x and Phi together) ---
    do step = 1, steps_actual

        ! --- Stage 1: evaluate at (t, x, phi) ---
        call eval_dynamics(n, x, u, nu, t, params, np, kx1, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call eval_jacobian(n, x, u, nu, t, params, np, F_jac, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call mat_mul_flat(n, F_jac, phi, kp1)
        kx1 = h * kx1
        kp1 = h * kp1

        ! --- Stage 2: evaluate at (t + h/2, x + kx1/2, phi + kp1/2) ---
        x_tmp   = x   + 0.5d0 * kx1
        phi_tmp = phi  + 0.5d0 * kp1
        call eval_dynamics(n, x_tmp, u, nu, t + 0.5d0*h, params, np, kx2, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call eval_jacobian(n, x_tmp, u, nu, t + 0.5d0*h, params, np, F_jac, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call mat_mul_flat(n, F_jac, phi_tmp, kp2)
        kx2 = h * kx2
        kp2 = h * kp2

        ! --- Stage 3: evaluate at (t + h/2, x + kx2/2, phi + kp2/2) ---
        x_tmp   = x   + 0.5d0 * kx2
        phi_tmp = phi  + 0.5d0 * kp2
        call eval_dynamics(n, x_tmp, u, nu, t + 0.5d0*h, params, np, kx3, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call eval_jacobian(n, x_tmp, u, nu, t + 0.5d0*h, params, np, F_jac, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call mat_mul_flat(n, F_jac, phi_tmp, kp3)
        kx3 = h * kx3
        kp3 = h * kp3

        ! --- Stage 4: evaluate at (t + h, x + kx3, phi + kp3) ---
        x_tmp   = x   + kx3
        phi_tmp = phi  + kp3
        call eval_dynamics(n, x_tmp, u, nu, t + h, params, np, kx4, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call eval_jacobian(n, x_tmp, u, nu, t + h, params, np, F_jac, dyn_info)
        if (dyn_info /= 0) then; info = 4; goto 999; end if
        call mat_mul_flat(n, F_jac, phi_tmp, kp4)
        kx4 = h * kx4
        kp4 = h * kp4

        ! --- Update state and STM ---
        x   = x   + (kx1 + 2.0d0*kx2 + 2.0d0*kx3 + kx4) / 6.0d0
        phi = phi  + (kp1 + 2.0d0*kp2 + 2.0d0*kp3 + kp4) / 6.0d0
        t   = t + h

    end do

999 continue
    deallocate(kx1, kx2, kx3, kx4, x_tmp)
    deallocate(kp1, kp2, kp3, kp4, phi_tmp)
    deallocate(F_jac, phi_dot_tmp)
    return

contains

    ! Evaluate dynamics: custom or built-in
    subroutine eval_dynamics(n, xv, uv, nu, tv, pv, np, xdot, ierr)
        integer, intent(in) :: n, nu, np
        real(c_double), intent(in) :: xv(n), uv(nu), pv(np)
        real(c_double), intent(in) :: tv
        real(c_double), intent(out) :: xdot(n)
        integer, intent(out) :: ierr
        if (use_custom_f) then
            call f_custom(n, xv, uv, nu, tv, pv, np, xdot, ierr)
        else
            call fa_dynamics_dispatch(model_id, n, xv, uv, nu, tv, pv, np, xdot, ierr)
        end if
    end subroutine

    ! Evaluate Jacobian: custom or built-in
    subroutine eval_jacobian(n, xv, uv, nu, tv, pv, np, Fv, ierr)
        integer, intent(in) :: n, nu, np
        real(c_double), intent(in) :: xv(n), uv(nu), pv(np)
        real(c_double), intent(in) :: tv
        real(c_double), intent(out) :: Fv(n*n)
        integer, intent(out) :: ierr
        if (use_custom_df) then
            call df_custom(n, xv, uv, nu, tv, pv, np, Fv, ierr)
        else
            call fa_dynamics_jacobian(model_id, n, xv, uv, nu, tv, pv, np, Fv, ierr)
        end if
    end subroutine

    ! Flat row-major matrix multiply: C = A * B, all n*n flat arrays
    ! C(i,j) = sum_k A(i,k) * B(k,j)
    ! Row-major: element (i,j) at index (i-1)*n + j
    subroutine mat_mul_flat(n, A, B, C_out)
        integer, intent(in) :: n
        real(c_double), intent(in)  :: A(n*n), B(n*n)
        real(c_double), intent(out) :: C_out(n*n)
        integer :: ii, jj, kk, idx
        real(c_double) :: s

        do ii = 1, n
            do jj = 1, n
                s = 0.0d0
                do kk = 1, n
                    s = s + A((ii-1)*n + kk) * B((kk-1)*n + jj)
                end do
                C_out((ii-1)*n + jj) = s
            end do
        end do
    end subroutine

end subroutine fa_propagate_stm


! ==========================================================================
! fa_propagate_batch — OpenMP parallel batch propagation
!
! Propagates n_states independent state vectors using the same dynamics
! model and control input. Useful for particle filters, Monte Carlo, etc.
!
! x_batch(n * n_states) — flat array, state i at x_batch((i-1)*n + 1 : i*n)
! info_batch(n_states)  — per-state error codes
!
! OpenMP parallelization when n_states >= 100.
! ==========================================================================
subroutine fa_propagate_batch(n, n_states, x_batch, u, nu, f_ptr, model_id, params, np, dt, n_steps, info_batch) &
    bind(C, name="fa_propagate_batch")
    use iso_c_binding
    implicit none

    integer(c_int), value, intent(in)  :: n, n_states, nu, model_id, np, n_steps
    real(c_double), intent(inout)      :: x_batch(n * n_states)
    real(c_double), intent(in)         :: u(nu), params(np)
    real(c_double), value, intent(in)  :: dt
    type(c_funptr), value              :: f_ptr
    integer(c_int), intent(out)        :: info_batch(n_states)

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

    integer :: i, offset

    ! Validate
    if (n <= 0 .or. n_states <= 0) then
        info_batch(1:n_states) = 3
        return
    end if

    ! Propagate each state independently
    ! Use OpenMP for large batches (>= 100 states)
    if (n_states >= 100) then
        !$omp parallel do default(none) &
        !$omp shared(n, n_states, x_batch, u, nu, f_ptr, model_id, params, np, dt, n_steps, info_batch) &
        !$omp private(i, offset)
        do i = 1, n_states
            offset = (i - 1) * n
            call fa_propagate(n, x_batch(offset+1:offset+n), u, nu, f_ptr, model_id, &
                              params, np, dt, n_steps, info_batch(i))
        end do
        !$omp end parallel do
    else
        do i = 1, n_states
            offset = (i - 1) * n
            call fa_propagate(n, x_batch(offset+1:offset+n), u, nu, f_ptr, model_id, &
                              params, np, dt, n_steps, info_batch(i))
        end do
    end if

end subroutine fa_propagate_batch
