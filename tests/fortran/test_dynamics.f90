! Copyright The Fantastic Planet — By David Clabaugh
!
! test_dynamics.f90 — Unit tests for forapollo_dynamics tracking models
!
! Tests constant-velocity, constant-acceleration, and constant-turn-rate
! dynamics models including analytic Jacobian verification via central
! finite differencing.
!
! Compile:
!   gfortran -O3 -std=f2008 -fall-intrinsics tests/fortran/test_dynamics.f90 \
!            build/obj/forapollo_dynamics.o -o build/test_dynamics -fopenmp
!
! Uses error stop on any failure.

program test_dynamics
    use iso_c_binding
    implicit none

    ! -------------------------------------------------------------------------
    ! Interface blocks for bind(C) subroutines.
    ! CRITICAL: without these, gfortran passes value args by reference → segfault.
    ! -------------------------------------------------------------------------
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

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    write(*,'(A)') '=== forApollo dynamics tests ==='
    write(*,*)

    ! --- Constant velocity tests ---
    call test_const_vel_basic()
    call test_const_vel_jacobian()

    ! --- Constant acceleration tests ---
    call test_const_accel_basic()
    call test_const_accel_jacobian()

    ! --- Constant turn rate tests ---
    call test_const_turn_basic()
    call test_const_turn_jacobian()

    ! --- Edge case / error tests ---
    call test_invalid_model_id()
    call test_const_vel_odd_n()
    call test_const_accel_bad_n()
    call test_const_turn_bad_n()

    ! --- Vehicle dynamics tests ---
    call test_bicycle_straight()
    call test_bicycle_jacobian()
    call test_ackermann_jacobian()
    call test_diffdrive_jacobian()

    ! --- Rigid body dynamics tests ---
    call test_rigidbody_torquefree()
    call test_rigidbody_jacobian()

    ! --- Aerial dynamics tests ---
    call test_quadrotor_hover()
    call test_quadrotor_jacobian()

    ! --- Stochastic dynamics tests ---
    call test_gbm_drift()
    call test_ou_mean_reversion()

    ! --- Scalar dynamics tests ---
    call test_spring_mass()
    call test_spring_mass_jacobian()

    ! --- Orbital dynamics tests ---
    call test_kepler_circular_orbit()
    call test_kepler_jacobian()
    call test_j2_perturbation()
    call test_j2_jacobian()
    call test_cr3bp_l1()
    call test_cr3bp_jacobian()
    call test_drag_model()
    call test_drag_jacobian()
    call test_orbital_edge_cases()

    ! --- Summary ---
    write(*,*)
    write(*,'(A,I0,A,I0,A)') '=== ', n_passed, ' passed, ', n_failed, ' failed ==='

    if (n_failed > 0) then
        error stop 'TESTS FAILED'
    end if

contains

    ! =========================================================================
    ! Helper: check if two vectors are approximately equal
    ! =========================================================================
    subroutine check_vec(test_name, n, actual, expected, tol, pass)
        character(len=*), intent(in)  :: test_name
        integer, intent(in)           :: n
        real(c_double), intent(in)    :: actual(n), expected(n)
        real(c_double), intent(in)    :: tol
        logical, intent(out)          :: pass
        integer :: i
        real(c_double) :: diff

        pass = .true.
        do i = 1, n
            diff = abs(actual(i) - expected(i))
            if (diff > tol) then
                write(*,'(A,A,A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') &
                    '  FAIL [', test_name, '] element ', i, &
                    ': got ', actual(i), ' expected ', expected(i), ' diff ', diff
                pass = .false.
            end if
        end do
    end subroutine check_vec

    ! =========================================================================
    ! Helper: report pass/fail
    ! =========================================================================
    subroutine report(test_name, pass)
        character(len=*), intent(in) :: test_name
        logical, intent(in)          :: pass

        if (pass) then
            write(*,'(A,A)') '  PASS: ', test_name
            n_passed = n_passed + 1
        else
            write(*,'(A,A)') '  FAIL: ', test_name
            n_failed = n_failed + 1
        end if
    end subroutine report

    ! =========================================================================
    ! Helper: verify analytic Jacobian against central finite difference
    !
    ! For each state element j, compute:
    !   F_fd(:, j) = (f(x + h*e_j) - f(x - h*e_j)) / (2h)
    ! then compare against analytic F.
    ! =========================================================================
    subroutine verify_jacobian(test_name, model_id_val, n_val, x, tol, pass)
        character(len=*), intent(in)  :: test_name
        integer, intent(in)           :: model_id_val, n_val
        real(c_double), intent(in)    :: x(n_val)
        real(c_double), intent(in)    :: tol
        logical, intent(out)          :: pass

        real(c_double) :: F_analytic(n_val * n_val)
        real(c_double) :: x_plus(n_val), x_minus(n_val)
        real(c_double) :: f_plus(n_val), f_minus(n_val)
        real(c_double) :: F_fd_col(n_val)
        real(c_double) :: dummy_u(1), dummy_p(1)
        real(c_double) :: h, diff
        integer(c_int) :: info
        integer :: i, j

        h = 1.0d-7

        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        ! Get analytic Jacobian
        call fa_dynamics_jacobian(model_id_val, n_val, x, dummy_u, 1, 0.0d0, dummy_p, 1, &
                                  F_analytic, info)
        if (info /= 0) then
            write(*,'(A,A,A,I0)') '  FAIL [', test_name, '] analytic Jacobian returned info=', info
            pass = .false.
            return
        end if

        pass = .true.

        ! Central finite difference for each column j
        do j = 1, n_val
            x_plus(1:n_val)  = x(1:n_val)
            x_minus(1:n_val) = x(1:n_val)
            x_plus(j)  = x(j) + h
            x_minus(j) = x(j) - h

            call fa_dynamics_dispatch(model_id_val, n_val, x_plus, dummy_u, 1, 0.0d0, &
                                      dummy_p, 1, f_plus, info)
            call fa_dynamics_dispatch(model_id_val, n_val, x_minus, dummy_u, 1, 0.0d0, &
                                      dummy_p, 1, f_minus, info)

            ! F_fd column j = (f_plus - f_minus) / (2h)
            F_fd_col(1:n_val) = (f_plus(1:n_val) - f_minus(1:n_val)) / (2.0d0 * h)

            ! Compare row i, col j: analytic F((i-1)*n + j) vs F_fd_col(i)
            do i = 1, n_val
                diff = abs(F_analytic((i-1)*n_val + j) - F_fd_col(i))
                if (diff > tol) then
                    write(*,'(A,A,A,I0,A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') &
                        '  FAIL [', test_name, '] F(', i, ',', j, &
                        '): analytic=', F_analytic((i-1)*n_val + j), &
                        ' fd=', F_fd_col(i), ' diff=', diff
                    pass = .false.
                end if
            end do
        end do

    end subroutine verify_jacobian

    ! =========================================================================
    ! Test: Constant velocity basic
    ! State: [pos_x, pos_y, vel_x, vel_y] = [0, 0, 1, 2]
    ! Expected x_dot = [1, 2, 0, 0]
    ! =========================================================================
    subroutine test_const_vel_basic()
        real(c_double) :: x(4), x_dot(4), expected(4)
        real(c_double) :: dummy_u(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = [0.0d0, 0.0d0, 1.0d0, 2.0d0]
        expected = [1.0d0, 2.0d0, 0.0d0, 0.0d0]
        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_dynamics_dispatch(40, 4, x, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot, info)

        pass = (info == 0)
        if (pass) call check_vec('const_vel_basic', 4, x_dot, expected, 1.0d-14, pass)
        call report('const_vel_basic', pass)
    end subroutine test_const_vel_basic

    ! =========================================================================
    ! Test: Constant velocity Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_const_vel_jacobian()
        real(c_double) :: x(4)
        logical :: pass

        x = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        call verify_jacobian('const_vel_jacobian', 40, 4, x, 1.0d-5, pass)
        call report('const_vel_jacobian', pass)
    end subroutine test_const_vel_jacobian

    ! =========================================================================
    ! Test: Constant acceleration basic (2D, n=6)
    ! State: [pos_x, pos_y, vel_x, vel_y, acc_x, acc_y] = [0,0,1,0,0,0]
    ! Expected x_dot = [1, 0, 0, 0, 0, 0]
    ! =========================================================================
    subroutine test_const_accel_basic()
        real(c_double) :: x(6), x_dot(6), expected(6)
        real(c_double) :: dummy_u(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0]
        expected = [1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_dynamics_dispatch(41, 6, x, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot, info)

        pass = (info == 0)
        if (pass) call check_vec('const_accel_basic', 6, x_dot, expected, 1.0d-14, pass)
        call report('const_accel_basic', pass)
    end subroutine test_const_accel_basic

    ! =========================================================================
    ! Test: Constant acceleration Jacobian
    ! =========================================================================
    subroutine test_const_accel_jacobian()
        real(c_double) :: x(6)
        logical :: pass

        x = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0]
        call verify_jacobian('const_accel_jacobian', 41, 6, x, 1.0d-5, pass)
        call report('const_accel_jacobian', pass)
    end subroutine test_const_accel_jacobian

    ! =========================================================================
    ! Test: Constant turn rate basic
    ! State: [x, y, vx, vy, omega] = [0, 0, 10, 0, 0.1]
    ! Expected x_dot = [10, 0, -0.1*0, 0.1*10, 0] = [10, 0, 0, 1, 0]
    ! =========================================================================
    subroutine test_const_turn_basic()
        real(c_double) :: x(5), x_dot(5), expected(5)
        real(c_double) :: dummy_u(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = [0.0d0, 0.0d0, 10.0d0, 0.0d0, 0.1d0]
        expected = [10.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0]
        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_dynamics_dispatch(42, 5, x, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot, info)

        pass = (info == 0)
        if (pass) call check_vec('const_turn_basic', 5, x_dot, expected, 1.0d-14, pass)
        call report('const_turn_basic', pass)
    end subroutine test_const_turn_basic

    ! =========================================================================
    ! Test: Constant turn rate Jacobian
    ! Use a state with nonzero vx, vy, omega to exercise all partial derivatives
    ! =========================================================================
    subroutine test_const_turn_jacobian()
        real(c_double) :: x(5)
        logical :: pass

        x = [1.0d0, 2.0d0, 5.0d0, 3.0d0, 0.2d0]
        call verify_jacobian('const_turn_jacobian', 42, 5, x, 1.0d-5, pass)
        call report('const_turn_jacobian', pass)
    end subroutine test_const_turn_jacobian

    ! =========================================================================
    ! Edge case: invalid model_id → info = 3
    ! =========================================================================
    subroutine test_invalid_model_id()
        real(c_double) :: x(4), x_dot(4)
        real(c_double) :: dummy_u(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_dynamics_dispatch(999, 4, x, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot, info)

        pass = (info == 3)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [invalid_model_id] expected info=3, got info=', info
        end if
        call report('invalid_model_id', pass)
    end subroutine test_invalid_model_id

    ! =========================================================================
    ! Edge case: odd n for constant velocity → info = 3
    ! =========================================================================
    subroutine test_const_vel_odd_n()
        real(c_double) :: x(3), x_dot(3)
        real(c_double) :: dummy_u(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = [1.0d0, 2.0d0, 3.0d0]
        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_dynamics_dispatch(40, 3, x, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot, info)

        pass = (info == 3)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [const_vel_odd_n] expected info=3, got info=', info
        end if
        call report('const_vel_odd_n', pass)
    end subroutine test_const_vel_odd_n

    ! =========================================================================
    ! Edge case: n not divisible by 3 for constant acceleration → info = 3
    ! =========================================================================
    subroutine test_const_accel_bad_n()
        real(c_double) :: x(4), x_dot(4)
        real(c_double) :: dummy_u(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_dynamics_dispatch(41, 4, x, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot, info)

        pass = (info == 3)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [const_accel_bad_n] expected info=3, got info=', info
        end if
        call report('const_accel_bad_n', pass)
    end subroutine test_const_accel_bad_n

    ! =========================================================================
    ! Edge case: n /= 5 for constant turn → info = 3
    ! =========================================================================
    subroutine test_const_turn_bad_n()
        real(c_double) :: x(4), x_dot(4)
        real(c_double) :: dummy_u(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        dummy_u(1) = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_dynamics_dispatch(42, 4, x, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot, info)

        pass = (info == 3)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [const_turn_bad_n] expected info=3, got info=', info
        end if
        call report('const_turn_bad_n', pass)
    end subroutine test_const_turn_bad_n

    ! =========================================================================
    ! Helper: verify analytic Jacobian against central finite difference
    ! with actual params array (for orbital models that need real parameters)
    ! =========================================================================
    subroutine verify_jacobian_params(test_name, model_id_val, n_val, x, &
                                       params_val, np_val, tol, pass)
        character(len=*), intent(in)  :: test_name
        integer, intent(in)           :: model_id_val, n_val, np_val
        real(c_double), intent(in)    :: x(n_val), params_val(np_val)
        real(c_double), intent(in)    :: tol
        logical, intent(out)          :: pass

        real(c_double) :: F_analytic(n_val * n_val)
        real(c_double) :: x_plus(n_val), x_minus(n_val)
        real(c_double) :: f_plus(n_val), f_minus(n_val)
        real(c_double) :: F_fd_col(n_val)
        real(c_double) :: dummy_u(1)
        real(c_double) :: h, diff, scale_h
        integer(c_int) :: info
        integer :: i, j

        dummy_u(1) = 0.0d0

        ! Get analytic Jacobian
        call fa_dynamics_jacobian(model_id_val, n_val, x, dummy_u, 1, 0.0d0, &
                                   params_val, np_val, F_analytic, info)
        if (info /= 0) then
            write(*,'(A,A,A,I0)') '  FAIL [', test_name, &
                '] analytic Jacobian returned info=', info
            pass = .false.
            return
        end if

        pass = .true.

        ! Central finite difference for each column j
        do j = 1, n_val
            ! Scale h relative to state magnitude for numerical stability
            scale_h = max(abs(x(j)), 1.0d0)
            h = 1.0d-7 * scale_h

            x_plus(1:n_val)  = x(1:n_val)
            x_minus(1:n_val) = x(1:n_val)
            x_plus(j)  = x(j) + h
            x_minus(j) = x(j) - h

            call fa_dynamics_dispatch(model_id_val, n_val, x_plus, dummy_u, 1, &
                                       0.0d0, params_val, np_val, f_plus, info)
            call fa_dynamics_dispatch(model_id_val, n_val, x_minus, dummy_u, 1, &
                                       0.0d0, params_val, np_val, f_minus, info)

            F_fd_col(1:n_val) = (f_plus(1:n_val) - f_minus(1:n_val)) / (2.0d0 * h)

            do i = 1, n_val
                diff = abs(F_analytic((i-1)*n_val + j) - F_fd_col(i))
                ! Use relative tolerance for large values
                if (diff > tol * max(abs(F_analytic((i-1)*n_val+j)), abs(F_fd_col(i)), 1.0d0)) then
                    write(*,'(A,A,A,I0,A,I0,A,ES14.7,A,ES14.7,A,ES10.3)') &
                        '  FAIL [', test_name, '] F(', i, ',', j, &
                        '): analytic=', F_analytic((i-1)*n_val + j), &
                        ' fd=', F_fd_col(i), ' diff=', diff
                    pass = .false.
                end if
            end do
        end do

    end subroutine verify_jacobian_params

    ! =========================================================================
    ! Helper: verify analytic Jacobian against central finite difference
    ! with actual params AND control arrays (for models that use u)
    ! =========================================================================
    subroutine verify_jacobian_params_ctrl(test_name, model_id_val, n_val, x, &
                                            u_val, nu_val, params_val, np_val, tol, pass)
        character(len=*), intent(in)  :: test_name
        integer, intent(in)           :: model_id_val, n_val, nu_val, np_val
        real(c_double), intent(in)    :: x(n_val), u_val(nu_val), params_val(np_val)
        real(c_double), intent(in)    :: tol
        logical, intent(out)          :: pass

        real(c_double) :: F_analytic(n_val * n_val)
        real(c_double) :: x_plus(n_val), x_minus(n_val)
        real(c_double) :: f_plus(n_val), f_minus(n_val)
        real(c_double) :: F_fd_col(n_val)
        real(c_double) :: h, diff, scale_h
        integer(c_int) :: info
        integer :: i, j

        ! Get analytic Jacobian
        call fa_dynamics_jacobian(model_id_val, n_val, x, u_val, nu_val, 0.0d0, &
                                   params_val, np_val, F_analytic, info)
        if (info /= 0) then
            write(*,'(A,A,A,I0)') '  FAIL [', test_name, &
                '] analytic Jacobian returned info=', info
            pass = .false.
            return
        end if

        pass = .true.

        ! Central finite difference for each column j
        do j = 1, n_val
            scale_h = max(abs(x(j)), 1.0d0)
            h = 1.0d-7 * scale_h

            x_plus(1:n_val)  = x(1:n_val)
            x_minus(1:n_val) = x(1:n_val)
            x_plus(j)  = x(j) + h
            x_minus(j) = x(j) - h

            call fa_dynamics_dispatch(model_id_val, n_val, x_plus, u_val, nu_val, &
                                       0.0d0, params_val, np_val, f_plus, info)
            call fa_dynamics_dispatch(model_id_val, n_val, x_minus, u_val, nu_val, &
                                       0.0d0, params_val, np_val, f_minus, info)

            F_fd_col(1:n_val) = (f_plus(1:n_val) - f_minus(1:n_val)) / (2.0d0 * h)

            do i = 1, n_val
                diff = abs(F_analytic((i-1)*n_val + j) - F_fd_col(i))
                if (diff > tol * max(abs(F_analytic((i-1)*n_val+j)), abs(F_fd_col(i)), 1.0d0)) then
                    write(*,'(A,A,A,I0,A,I0,A,ES14.7,A,ES14.7,A,ES10.3)') &
                        '  FAIL [', test_name, '] F(', i, ',', j, &
                        '): analytic=', F_analytic((i-1)*n_val + j), &
                        ' fd=', F_fd_col(i), ' diff=', diff
                    pass = .false.
                end if
            end do
        end do

    end subroutine verify_jacobian_params_ctrl

    ! =========================================================================
    ! Test: Bicycle straight line
    ! v=10, theta=0, steer_rate=0, accel=0 → x_dot=[10,0,0,0]
    ! =========================================================================
    subroutine test_bicycle_straight()
        real(c_double) :: x(4), x_dot(4), expected(4)
        real(c_double) :: u(2), params(1)
        integer(c_int) :: info
        logical :: pass

        x = [0.0d0, 0.0d0, 0.0d0, 10.0d0]    ! [x, y, theta, v]
        u = [0.0d0, 0.0d0]                     ! [accel, steer_rate]
        params(1) = 2.5d0                       ! wheelbase L
        expected = [10.0d0, 0.0d0, 0.0d0, 0.0d0]

        call fa_dynamics_dispatch(20, 4, x, u, 2, 0.0d0, params, 1, x_dot, info)

        pass = (info == 0)
        if (pass) call check_vec('bicycle_straight', 4, x_dot, expected, 1.0d-14, pass)
        call report('bicycle_straight', pass)
    end subroutine test_bicycle_straight

    ! =========================================================================
    ! Test: Bicycle Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_bicycle_jacobian()
        real(c_double) :: x(4), u(2), params(1)
        logical :: pass

        x = [1.0d0, 2.0d0, 0.3d0, 8.0d0]     ! nonzero theta and v
        u = [1.0d0, 0.2d0]
        params(1) = 2.5d0

        call verify_jacobian_params_ctrl('bicycle_jacobian', 20, 4, x, u, 2, params, 1, 1.0d-5, pass)
        call report('bicycle_jacobian', pass)
    end subroutine test_bicycle_jacobian

    ! =========================================================================
    ! Test: Ackermann Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_ackermann_jacobian()
        real(c_double) :: x(5), u(2), params(1)
        logical :: pass

        x = [1.0d0, 2.0d0, 0.3d0, 5.0d0, 0.1d0]  ! nonzero steer angle
        u = [0.5d0, 0.1d0]
        params(1) = 2.7d0    ! wheelbase

        call verify_jacobian_params_ctrl('ackermann_jacobian', 21, 5, x, u, 2, params, 1, 1.0d-5, pass)
        call report('ackermann_jacobian', pass)
    end subroutine test_ackermann_jacobian

    ! =========================================================================
    ! Test: Differential drive Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_diffdrive_jacobian()
        real(c_double) :: x(4), u(2), params(1)
        logical :: pass

        x = [1.0d0, 2.0d0, 0.5d0, 3.0d0]
        u = [1.0d0, 1.5d0]        ! left/right wheel velocities
        params(1) = 0.5d0          ! wheel separation

        call verify_jacobian_params_ctrl('diffdrive_jacobian', 22, 4, x, u, 2, params, 1, 1.0d-5, pass)
        call report('diffdrive_jacobian', pass)
    end subroutine test_diffdrive_jacobian

    ! =========================================================================
    ! Test: 6-DOF rigid body torque-free spin
    ! omega=[0,0,1], no forces/torques → omega_dot=0, quat rotates
    ! =========================================================================
    subroutine test_rigidbody_torquefree()
        real(c_double) :: x(13), x_dot(13), u(6), params(4)
        integer(c_int) :: info
        logical :: pass

        ! State: at origin, no velocity, identity quaternion, spinning about z
        x = [0.0d0, 0.0d0, 0.0d0, &      ! position
             0.0d0, 0.0d0, 0.0d0, &       ! velocity
             1.0d0, 0.0d0, 0.0d0, 0.0d0, & ! quaternion (identity)
             0.0d0, 0.0d0, 1.0d0]          ! omega = [0,0,1]

        u = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]  ! no forces/torques
        params = [10.0d0, 1.0d0, 1.0d0, 1.0d0]  ! mass=10, Ixx=Iyy=Izz=1 (sphere)

        call fa_dynamics_dispatch(10, 13, x, u, 6, 0.0d0, params, 4, x_dot, info)

        pass = (info == 0)

        ! omega_dot should be zero (torque-free, symmetric body)
        if (pass .and. abs(x_dot(11)) > 1.0d-14) pass = .false.
        if (pass .and. abs(x_dot(12)) > 1.0d-14) pass = .false.
        if (pass .and. abs(x_dot(13)) > 1.0d-14) pass = .false.

        ! Quaternion should be rotating: q_dot = 0.5*q(*)[0,0,0,1]
        ! q=[1,0,0,0], omega=[0,0,1]
        ! q0_dot = -0.5*(0+0+0) = 0
        ! q3_dot = 0.5*(0+0+1*1) = 0.5
        if (pass .and. abs(x_dot(7) - 0.0d0) > 1.0d-14) pass = .false.
        if (pass .and. abs(x_dot(10) - 0.5d0) > 1.0d-14) pass = .false.

        call report('rigidbody_torquefree', pass)
    end subroutine test_rigidbody_torquefree

    ! =========================================================================
    ! Test: 6-DOF rigid body Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_rigidbody_jacobian()
        real(c_double) :: x(13), u(6), params(4)
        logical :: pass

        ! General state with nonzero quaternion and angular velocity
        x = [1.0d0, 2.0d0, 3.0d0, &        ! position
             0.1d0, -0.2d0, 0.3d0, &        ! velocity
             0.5d0, 0.5d0, 0.5d0, 0.5d0, &  ! quaternion (normalized)
             0.1d0, 0.2d0, 0.3d0]           ! omega

        u = [1.0d0, 0.5d0, -0.3d0, 0.1d0, -0.1d0, 0.2d0]  ! forces and torques
        params = [10.0d0, 2.0d0, 3.0d0, 4.0d0]  ! mass, Ixx, Iyy, Izz

        call verify_jacobian_params_ctrl('rigidbody_jacobian', 10, 13, x, u, 6, params, 4, 1.0d-5, pass)
        call report('rigidbody_jacobian', pass)
    end subroutine test_rigidbody_jacobian

    ! =========================================================================
    ! Test: Quadrotor hover
    ! thrust=mass*g, zero angles/rates → v_dot ≈ [0,0,0]
    ! =========================================================================
    subroutine test_quadrotor_hover()
        real(c_double) :: x(12), x_dot(12), u(4), params(5)
        integer(c_int) :: info
        logical :: pass
        real(c_double) :: mass, g

        mass = 1.5d0
        g    = 9.81d0

        ! Hovering: zero position, zero velocity, zero angles, zero body rates
        x = [0.0d0, 0.0d0, 1.0d0, &   ! position (1m altitude)
             0.0d0, 0.0d0, 0.0d0, &   ! velocity = 0
             0.0d0, 0.0d0, 0.0d0, &   ! phi, theta, psi = 0
             0.0d0, 0.0d0, 0.0d0]     ! p, q, r = 0

        u = [mass*g, 0.0d0, 0.0d0, 0.0d0]  ! thrust = weight, no torques
        params = [mass, 0.01d0, 0.01d0, 0.02d0, g]  ! mass, Ixx, Iyy, Izz, g

        call fa_dynamics_dispatch(30, 12, x, u, 4, 0.0d0, params, 5, x_dot, info)

        pass = (info == 0)

        ! velocity derivatives should be approximately zero (hover)
        if (pass .and. abs(x_dot(4)) > 1.0d-10) then
            write(*,'(A,ES14.7)') '  FAIL [quadrotor_hover] vx_dot=', x_dot(4)
            pass = .false.
        end if
        if (pass .and. abs(x_dot(5)) > 1.0d-10) then
            write(*,'(A,ES14.7)') '  FAIL [quadrotor_hover] vy_dot=', x_dot(5)
            pass = .false.
        end if
        if (pass .and. abs(x_dot(6)) > 1.0d-10) then
            write(*,'(A,ES14.7)') '  FAIL [quadrotor_hover] vz_dot=', x_dot(6)
            pass = .false.
        end if

        call report('quadrotor_hover', pass)
    end subroutine test_quadrotor_hover

    ! =========================================================================
    ! Test: Quadrotor Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_quadrotor_jacobian()
        real(c_double) :: x(12), u(4), params(5)
        logical :: pass

        ! General state with small angles and rates to avoid singularity
        x = [1.0d0, 2.0d0, 3.0d0, &      ! position
             0.5d0, -0.3d0, 0.1d0, &      ! velocity
             0.1d0, 0.05d0, 0.2d0, &      ! phi, theta, psi (small angles)
             0.01d0, -0.02d0, 0.03d0]     ! p, q, r

        u = [15.0d0, 0.1d0, -0.05d0, 0.02d0]
        params = [1.5d0, 0.01d0, 0.01d0, 0.02d0, 9.81d0]

        call verify_jacobian_params_ctrl('quadrotor_jacobian', 30, 12, x, u, 4, params, 5, 1.0d-4, pass)
        call report('quadrotor_jacobian', pass)
    end subroutine test_quadrotor_jacobian

    ! =========================================================================
    ! Test: GBM drift
    ! S=100, mu=0.05 → x_dot = 5.0
    ! =========================================================================
    subroutine test_gbm_drift()
        real(c_double) :: x(1), x_dot(1), params(2)
        real(c_double) :: dummy_u(1)
        integer(c_int) :: info
        logical :: pass

        x(1) = 100.0d0
        params = [0.05d0, 0.2d0]     ! mu_drift=0.05, sigma=0.2
        dummy_u(1) = 0.0d0

        call fa_dynamics_dispatch(50, 1, x, dummy_u, 1, 0.0d0, params, 2, x_dot, info)

        pass = (info == 0)
        if (pass .and. abs(x_dot(1) - 5.0d0) > 1.0d-14) then
            write(*,'(A,ES14.7)') '  FAIL [gbm_drift] x_dot=', x_dot(1)
            pass = .false.
        end if
        call report('gbm_drift', pass)
    end subroutine test_gbm_drift

    ! =========================================================================
    ! Test: Ornstein-Uhlenbeck mean reversion
    ! X=110, theta=0.5, mu=100 → x_dot = 0.5*(100-110) = -5.0
    ! =========================================================================
    subroutine test_ou_mean_reversion()
        real(c_double) :: x(1), x_dot(1), params(3)
        real(c_double) :: dummy_u(1)
        integer(c_int) :: info
        logical :: pass

        x(1) = 110.0d0
        params = [0.5d0, 100.0d0, 0.3d0]   ! theta=0.5, mu=100, sigma=0.3
        dummy_u(1) = 0.0d0

        call fa_dynamics_dispatch(51, 1, x, dummy_u, 1, 0.0d0, params, 3, x_dot, info)

        pass = (info == 0)
        if (pass .and. abs(x_dot(1) - (-5.0d0)) > 1.0d-14) then
            write(*,'(A,ES14.7)') '  FAIL [ou_mean_reversion] x_dot=', x_dot(1)
            pass = .false.
        end if
        call report('ou_mean_reversion', pass)
    end subroutine test_ou_mean_reversion

    ! =========================================================================
    ! Test: Spring-mass-damper
    ! x=1, v=0, k=10, c=1, m=1, F=0 → x_dot=[0, -10]
    ! =========================================================================
    subroutine test_spring_mass()
        real(c_double) :: x(2), x_dot(2), expected(2)
        real(c_double) :: u(1), params(3)
        integer(c_int) :: info
        logical :: pass

        x = [1.0d0, 0.0d0]        ! x=1, v=0
        u(1) = 0.0d0              ! no external force
        params = [10.0d0, 1.0d0, 1.0d0]   ! k=10, c=1, m=1
        expected = [0.0d0, -10.0d0]  ! v_dot = (0 - 10*1 - 1*0)/1 = -10

        call fa_dynamics_dispatch(61, 2, x, u, 1, 0.0d0, params, 3, x_dot, info)

        pass = (info == 0)
        if (pass) call check_vec('spring_mass', 2, x_dot, expected, 1.0d-14, pass)
        call report('spring_mass', pass)
    end subroutine test_spring_mass

    ! =========================================================================
    ! Test: Spring-mass-damper Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_spring_mass_jacobian()
        real(c_double) :: x(2), u(1), params(3)
        logical :: pass

        x = [0.5d0, 1.0d0]
        u(1) = 2.0d0
        params = [10.0d0, 1.0d0, 2.0d0]   ! k=10, c=1, m=2

        call verify_jacobian_params_ctrl('spring_mass_jacobian', 61, 2, x, u, 1, params, 3, 1.0d-5, pass)
        call report('spring_mass_jacobian', pass)
    end subroutine test_spring_mass_jacobian

    ! =========================================================================
    ! Test: Kepler circular orbit
    ! x=[6778e3, 0, 0, 0, 7669, 0], mu=3.986004418e14
    ! Verify velocity propagation and radial gravity
    ! =========================================================================
    subroutine test_kepler_circular_orbit()
        real(c_double) :: x(6), x_dot(6), params(1)
        real(c_double) :: dummy_u(1)
        integer(c_int) :: info
        logical :: pass
        real(c_double) :: r, expected_grav

        x = [6778.0d3, 0.0d0, 0.0d0, 0.0d0, 7669.0d0, 0.0d0]
        params(1) = 3.986004418d14  ! mu for Earth
        dummy_u(1) = 0.0d0

        call fa_dynamics_dispatch(1, 6, x, dummy_u, 1, 0.0d0, params, 1, x_dot, info)

        pass = (info == 0)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [kepler_circular] info=', info
            call report('kepler_circular_orbit', pass)
            return
        end if

        ! x_dot(1:3) should be velocity = [0, 7669, 0]
        if (abs(x_dot(1) - 0.0d0) > 1.0d-10) pass = .false.
        if (abs(x_dot(2) - 7669.0d0) > 1.0d-10) pass = .false.
        if (abs(x_dot(3) - 0.0d0) > 1.0d-10) pass = .false.

        ! x_dot(4) should be -mu/r^2 (radial gravity)
        r = 6778.0d3
        expected_grav = -params(1) / (r * r)   ! approx -8.674
        if (abs(x_dot(4) - expected_grav) > 1.0d-3) then
            write(*,'(A,ES14.7,A,ES14.7)') '  FAIL [kepler_circular] ax=', &
                x_dot(4), ' expected=', expected_grav
            pass = .false.
        end if

        ! x_dot(5) and x_dot(6) should be 0 (no tangential/out-of-plane gravity)
        if (abs(x_dot(5)) > 1.0d-10) pass = .false.
        if (abs(x_dot(6)) > 1.0d-10) pass = .false.

        call report('kepler_circular_orbit', pass)
    end subroutine test_kepler_circular_orbit

    ! =========================================================================
    ! Test: Kepler Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_kepler_jacobian()
        real(c_double) :: x(6), params(1)
        logical :: pass

        x = [6778.0d3, 0.0d0, 0.0d0, 0.0d0, 7669.0d0, 0.0d0]
        params(1) = 3.986004418d14

        call verify_jacobian_params('kepler_jacobian', 1, 6, x, params, 1, 1.0d-5, pass)
        call report('kepler_jacobian', pass)
    end subroutine test_kepler_jacobian

    ! =========================================================================
    ! Test: J2 perturbation
    ! Same orbit + J2=1.08263e-3, R_eq=6378137
    ! At equatorial orbit (z=0), J2 should add non-zero x-acceleration
    ! =========================================================================
    subroutine test_j2_perturbation()
        real(c_double) :: x(6), x_dot_kep(6), x_dot_j2(6), params_kep(1), params_j2(3)
        real(c_double) :: dummy_u(1)
        integer(c_int) :: info
        logical :: pass
        real(c_double) :: delta_ax

        x = [6778.0d3, 0.0d0, 0.0d0, 0.0d0, 7669.0d0, 0.0d0]
        params_kep(1) = 3.986004418d14
        params_j2 = [3.986004418d14, 1.08263d-3, 6378137.0d0]
        dummy_u(1) = 0.0d0

        ! Get pure Kepler
        call fa_dynamics_dispatch(1, 6, x, dummy_u, 1, 0.0d0, params_kep, 1, x_dot_kep, info)
        ! Get J2
        call fa_dynamics_dispatch(2, 6, x, dummy_u, 1, 0.0d0, params_j2, 3, x_dot_j2, info)

        pass = (info == 0)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [j2_perturbation] info=', info
            call report('j2_perturbation', pass)
            return
        end if

        ! J2 should add a non-zero perturbation to x-acceleration
        delta_ax = x_dot_j2(4) - x_dot_kep(4)
        if (abs(delta_ax) < 1.0d-10) then
            write(*,'(A,ES14.7)') '  FAIL [j2_perturbation] delta_ax too small: ', delta_ax
            pass = .false.
        end if

        ! Velocity components should be identical (J2 only affects acceleration)
        if (abs(x_dot_j2(1) - x_dot_kep(1)) > 1.0d-14) pass = .false.
        if (abs(x_dot_j2(2) - x_dot_kep(2)) > 1.0d-14) pass = .false.
        if (abs(x_dot_j2(3) - x_dot_kep(3)) > 1.0d-14) pass = .false.

        call report('j2_perturbation', pass)
    end subroutine test_j2_perturbation

    ! =========================================================================
    ! Test: J2 Jacobian (finite-difference verification)
    ! Use an inclined orbit state to exercise z-dependent terms
    ! =========================================================================
    subroutine test_j2_jacobian()
        real(c_double) :: x(6), params(3)
        logical :: pass

        ! Inclined orbit: nonzero z component to exercise all J2 partials
        x = [4000.0d3, 3000.0d3, 2000.0d3, -1000.0d0, 5000.0d0, 3000.0d0]
        params = [3.986004418d14, 1.08263d-3, 6378137.0d0]

        call verify_jacobian_params('j2_jacobian', 2, 6, x, params, 3, 1.0d-4, pass)
        call report('j2_jacobian', pass)
    end subroutine test_j2_jacobian

    ! =========================================================================
    ! Test: CR3BP near L1
    ! Earth-Moon system, mu_ratio ~ 0.01215
    ! State near L1: [0.8369, 0, 0, 0, 0, 0]
    ! Accelerations should be small (near equilibrium)
    ! =========================================================================
    subroutine test_cr3bp_l1()
        real(c_double) :: x(6), x_dot(6), params(1)
        real(c_double) :: dummy_u(1)
        integer(c_int) :: info
        logical :: pass
        real(c_double) :: accel_mag

        x = [0.8369d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
        params(1) = 0.01215d0  ! Earth-Moon mass ratio
        dummy_u(1) = 0.0d0

        call fa_dynamics_dispatch(3, 6, x, dummy_u, 1, 0.0d0, params, 1, x_dot, info)

        pass = (info == 0)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [cr3bp_l1] info=', info
            call report('cr3bp_l1_equilibrium', pass)
            return
        end if

        ! Velocity derivatives should be zero (zero velocity in)
        if (abs(x_dot(1)) > 1.0d-14) pass = .false.
        if (abs(x_dot(2)) > 1.0d-14) pass = .false.
        if (abs(x_dot(3)) > 1.0d-14) pass = .false.

        ! Accelerations should be small near L1 (not exact, so allow tolerance)
        ! L1 is approximate, so just check magnitude is reasonably small
        accel_mag = sqrt(x_dot(4)**2 + x_dot(5)**2 + x_dot(6)**2)
        if (accel_mag > 0.5d0) then
            write(*,'(A,ES14.7)') '  FAIL [cr3bp_l1] accel magnitude too large: ', accel_mag
            pass = .false.
        end if

        call report('cr3bp_l1_equilibrium', pass)
    end subroutine test_cr3bp_l1

    ! =========================================================================
    ! Test: CR3BP Jacobian (finite-difference verification)
    ! Use a general state with nonzero velocity for full partial coverage
    ! =========================================================================
    subroutine test_cr3bp_jacobian()
        real(c_double) :: x(6), params(1)
        logical :: pass

        x = [0.5d0, 0.1d0, 0.05d0, 0.01d0, -0.02d0, 0.005d0]
        params(1) = 0.01215d0

        call verify_jacobian_params('cr3bp_jacobian', 3, 6, x, params, 1, 1.0d-5, pass)
        call report('cr3bp_jacobian', pass)
    end subroutine test_cr3bp_jacobian

    ! =========================================================================
    ! Test: Drag model
    ! LEO orbit with drag. Verify drag acceleration opposes velocity direction.
    ! =========================================================================
    subroutine test_drag_model()
        real(c_double) :: x_kep(6), x_drag(7), x_dot_kep(6), x_dot_drag(7)
        real(c_double) :: params_kep(1), params_drag(4)
        real(c_double) :: dummy_u(1)
        integer(c_int) :: info
        logical :: pass
        real(c_double) :: drag_ax, drag_ay, drag_az, dot_product_val

        ! Kepler reference (no drag)
        x_kep = [6778.0d3, 0.0d0, 0.0d0, 0.0d0, 7669.0d0, 0.0d0]
        params_kep(1) = 3.986004418d14
        dummy_u(1) = 0.0d0

        call fa_dynamics_dispatch(1, 6, x_kep, dummy_u, 1, 0.0d0, params_kep, 1, x_dot_kep, info)

        ! Drag model: same orbit with ballistic coefficient
        x_drag = [6778.0d3, 0.0d0, 0.0d0, 0.0d0, 7669.0d0, 0.0d0, 0.01d0]
        ! params: mu, rho0, h_scale, R_body
        params_drag = [3.986004418d14, 1.225d0, 8500.0d0, 6371.0d3]

        call fa_dynamics_dispatch(4, 7, x_drag, dummy_u, 1, 0.0d0, params_drag, 4, x_dot_drag, info)

        pass = (info == 0)
        if (.not. pass) then
            write(*,'(A,I0)') '  FAIL [drag_model] info=', info
            call report('drag_model', pass)
            return
        end if

        ! Extract drag-only acceleration (subtract Kepler part)
        drag_ax = x_dot_drag(4) - x_dot_kep(4)
        drag_ay = x_dot_drag(5) - x_dot_kep(5)
        drag_az = x_dot_drag(6) - x_dot_kep(6)

        ! Drag should oppose velocity: dot(a_drag, v) < 0
        dot_product_val = drag_ax * x_drag(4) + drag_ay * x_drag(5) + drag_az * x_drag(6)
        if (dot_product_val >= 0.0d0) then
            write(*,'(A,ES14.7)') '  FAIL [drag_model] drag not opposing velocity, dot=', &
                dot_product_val
            pass = .false.
        end if

        ! beta_dot should be zero
        if (abs(x_dot_drag(7)) > 1.0d-14) then
            write(*,'(A,ES14.7)') '  FAIL [drag_model] beta_dot not zero: ', x_dot_drag(7)
            pass = .false.
        end if

        call report('drag_model', pass)
    end subroutine test_drag_model

    ! =========================================================================
    ! Test: Drag Jacobian (finite-difference verification)
    ! =========================================================================
    subroutine test_drag_jacobian()
        real(c_double) :: x(7), params(4)
        logical :: pass

        x = [6778.0d3, 100.0d3, 50.0d3, -200.0d0, 7669.0d0, 100.0d0, 0.01d0]
        params = [3.986004418d14, 1.225d0, 8500.0d0, 6371.0d3]

        call verify_jacobian_params('drag_jacobian', 4, 7, x, params, 4, 1.0d-4, pass)
        call report('drag_jacobian', pass)
    end subroutine test_drag_jacobian

    ! =========================================================================
    ! Test: Orbital edge cases
    ! n!=6 for Kepler -> info=3, np<1 for Kepler -> info=3
    ! =========================================================================
    subroutine test_orbital_edge_cases()
        real(c_double) :: x4(4), x_dot4(4), x6(6), x_dot6(6)
        real(c_double) :: dummy_u(1), dummy_p(1)
        real(c_double) :: empty_p(1)
        integer(c_int) :: info
        logical :: pass

        dummy_u(1) = 0.0d0
        dummy_p(1) = 3.986004418d14

        pass = .true.

        ! Kepler with n=4 (wrong state dimension) -> info=3
        x4 = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        call fa_dynamics_dispatch(1, 4, x4, dummy_u, 1, 0.0d0, dummy_p, 1, x_dot4, info)
        if (info /= 3) then
            write(*,'(A,I0)') '  FAIL [orbital_edge_cases] Kepler n=4 expected info=3, got ', info
            pass = .false.
        end if

        ! Kepler with np=0 -> info=3
        ! We pass np=0 with a dummy array; the subroutine checks np<1
        empty_p(1) = 0.0d0
        x6 = [6778.0d3, 0.0d0, 0.0d0, 0.0d0, 7669.0d0, 0.0d0]
        call fa_dynamics_dispatch(1, 6, x6, dummy_u, 1, 0.0d0, empty_p, 0, x_dot6, info)
        if (info /= 3) then
            write(*,'(A,I0)') '  FAIL [orbital_edge_cases] Kepler np=0 expected info=3, got ', info
            pass = .false.
        end if

        call report('orbital_edge_cases', pass)
    end subroutine test_orbital_edge_cases

end program test_dynamics
