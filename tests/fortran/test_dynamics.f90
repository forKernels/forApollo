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

end program test_dynamics
