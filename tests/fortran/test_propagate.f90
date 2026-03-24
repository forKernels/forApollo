! Copyright The Fantastic Planet — By David Clabaugh
!
! test_propagate.f90 — Unit tests for forapollo_propagate (RK4, STM, batch)
!
! Compile:
!   gfortran -O3 -std=f2008 -fall-intrinsics tests/fortran/test_propagate.f90 \
!            build/obj/forapollo_dynamics.o build/obj/forapollo_propagate.o \
!            -o build/test_propagate -fopenmp
!
! Uses error stop on any failure.

program test_propagate
    use iso_c_binding
    implicit none

    ! -------------------------------------------------------------------------
    ! Interface blocks for bind(C) subroutines
    ! -------------------------------------------------------------------------
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

        subroutine fa_propagate_batch(n, n_states, x_batch, u, nu, f_ptr, model_id, params, np, dt, n_steps, info_batch) &
            bind(C, name="fa_propagate_batch")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, n_states, nu, model_id, np, n_steps
            real(c_double), intent(inout)      :: x_batch(n * n_states)
            real(c_double), intent(in)         :: u(nu), params(np)
            real(c_double), value, intent(in)  :: dt
            type(c_funptr), value              :: f_ptr
            integer(c_int), intent(out)        :: info_batch(n_states)
        end subroutine

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

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    write(*,'(A)') '=== forApollo propagation tests ==='
    write(*,*)

    ! --- Single state propagation tests ---
    call test_const_vel_single_step()
    call test_const_vel_multi_step()
    call test_kepler_quarter_orbit()

    ! --- STM propagation tests ---
    call test_stm_const_vel()

    ! --- Batch propagation tests ---
    call test_batch_const_vel()

    ! --- Edge case tests ---
    call test_invalid_n()
    call test_dt_zero()

    ! --- Summary ---
    write(*,*)
    write(*,'(A)') '==================================='
    write(*,'(A,I0,A,I0,A)') 'Results: ', n_passed, ' passed, ', n_failed, ' failed'
    if (n_failed > 0) then
        write(*,'(A)') 'SOME TESTS FAILED'
        error stop 1
    else
        write(*,'(A)') 'ALL TESTS PASSED'
    end if

contains

    ! =====================================================================
    ! Test 1: Const-vel propagation, single step
    ! x=[0,0,1,0], dt=10, n_steps=1 → x=[10,0,1,0] (exact for linear)
    ! =====================================================================
    subroutine test_const_vel_single_step()
        real(c_double) :: x(4), u(1), params(1)
        integer(c_int) :: info
        type(c_funptr) :: null_ptr

        write(*,'(A)', advance='no') '  const_vel single step ... '

        x = [0.0d0, 0.0d0, 1.0d0, 0.0d0]
        u(1) = 0.0d0
        params(1) = 0.0d0
        null_ptr = c_null_funptr

        call fa_propagate(4, x, u, 1, null_ptr, 40, params, 1, 10.0d0, 1, info)

        if (info /= 0) then
            write(*,'(A,I0)') 'FAIL (info=', info
            n_failed = n_failed + 1
            return
        end if

        ! RK4 on linear dynamics is exact
        if (abs(x(1) - 10.0d0) > 1.0d-10 .or. abs(x(2)) > 1.0d-10 .or. &
            abs(x(3) - 1.0d0) > 1.0d-10 .or. abs(x(4)) > 1.0d-10) then
            write(*,'(A)') 'FAIL (wrong state)'
            write(*,'(A,4F12.6)') '    got: ', x
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 2: Const-vel propagation, multi-step (should be identical)
    ! =====================================================================
    subroutine test_const_vel_multi_step()
        real(c_double) :: x(4), u(1), params(1)
        integer(c_int) :: info
        type(c_funptr) :: null_ptr

        write(*,'(A)', advance='no') '  const_vel multi step  ... '

        x = [0.0d0, 0.0d0, 1.0d0, 0.0d0]
        u(1) = 0.0d0
        params(1) = 0.0d0
        null_ptr = c_null_funptr

        call fa_propagate(4, x, u, 1, null_ptr, 40, params, 1, 10.0d0, 10, info)

        if (info /= 0) then
            write(*,'(A,I0)') 'FAIL (info=', info
            n_failed = n_failed + 1
            return
        end if

        if (abs(x(1) - 10.0d0) > 1.0d-10 .or. abs(x(2)) > 1.0d-10 .or. &
            abs(x(3) - 1.0d0) > 1.0d-10 .or. abs(x(4)) > 1.0d-10) then
            write(*,'(A)') 'FAIL (wrong state)'
            write(*,'(A,4F12.6)') '    got: ', x
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 3: Kepler circular orbit — quarter orbit propagation
    ! Initial: r=[6778e3, 0, 0], v=[0, 7669, 0] (LEO circular)
    ! mu = 3.986004418e14
    ! After quarter orbit (~T/4): r≈[0, 6778e3, 0], v≈[-7669, 0, 0]
    ! Tolerance: ~1% for RK4 with 100 steps
    ! =====================================================================
    subroutine test_kepler_quarter_orbit()
        real(c_double) :: x(6), u(1), params(1)
        real(c_double) :: mu, r0, v0, period, dt
        integer(c_int) :: info
        type(c_funptr) :: null_ptr
        real(c_double) :: tol

        write(*,'(A)', advance='no') '  kepler quarter orbit  ... '

        r0 = 6778.0d3        ! meters
        mu = 3.986004418d14  ! m^3/s^2
        v0 = sqrt(mu / r0)   ! circular orbit speed

        ! Orbital period T = 2*pi*sqrt(r0^3/mu)
        period = 2.0d0 * acos(-1.0d0) * sqrt(r0**3 / mu)
        dt = period / 4.0d0   ! quarter orbit

        x = [r0, 0.0d0, 0.0d0, 0.0d0, v0, 0.0d0]
        u(1) = 0.0d0
        params(1) = mu
        null_ptr = c_null_funptr

        call fa_propagate(6, x, u, 1, null_ptr, 1, params, 1, dt, 100, info)

        if (info /= 0) then
            write(*,'(A,I0)') 'FAIL (info=', info
            n_failed = n_failed + 1
            return
        end if

        ! After quarter orbit: should be near [0, r0, 0, -v0, 0, 0]
        ! Use 1% tolerance (RK4 with 100 steps on nonlinear Kepler)
        tol = 0.01d0

        if (abs(x(1)) / r0 > tol .or. &
            abs(x(2) - r0) / r0 > tol .or. &
            abs(x(3)) / r0 > tol .or. &
            abs(x(4) + v0) / v0 > tol .or. &
            abs(x(5)) / v0 > tol .or. &
            abs(x(6)) / v0 > tol) then
            write(*,'(A)') 'FAIL (orbit error > 1%)'
            write(*,'(A,3ES14.6)') '    r: ', x(1:3)
            write(*,'(A,3ES14.6)') '    v: ', x(4:6)
            write(*,'(A,3ES14.6)') '    expect r: ', 0.0d0, r0, 0.0d0
            write(*,'(A,3ES14.6)') '    expect v: ', -v0, 0.0d0, 0.0d0
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 4: STM propagation — const-vel, dt=10
    ! For linear dynamics, Phi should be exact:
    !   [1 0 dt 0]
    !   [0 1 0  dt]
    !   [0 0 1  0 ]
    !   [0 0 0  1 ]
    ! (row-major flat)
    ! =====================================================================
    subroutine test_stm_const_vel()
        real(c_double) :: x(4), phi(16), u(1), params(1)
        integer(c_int) :: info
        type(c_funptr) :: null_f, null_df
        real(c_double) :: phi_expected(16)
        real(c_double) :: dt_val
        integer :: i

        write(*,'(A)', advance='no') '  STM const_vel         ... '

        dt_val = 10.0d0
        x = [0.0d0, 0.0d0, 1.0d0, 0.0d0]
        u(1) = 0.0d0
        params(1) = 0.0d0
        null_f  = c_null_funptr
        null_df = c_null_funptr

        ! Initialize Phi = Identity (row-major flat)
        phi = 0.0d0
        phi(1)  = 1.0d0   ! (1,1)
        phi(6)  = 1.0d0   ! (2,2)
        phi(11) = 1.0d0   ! (3,3)
        phi(16) = 1.0d0   ! (4,4)

        call fa_propagate_stm(4, x, phi, u, 1, null_f, null_df, 40, params, 1, dt_val, 1, info)

        if (info /= 0) then
            write(*,'(A,I0)') 'FAIL (info=', info
            n_failed = n_failed + 1
            return
        end if

        ! Check state (same as single-step test)
        if (abs(x(1) - 10.0d0) > 1.0d-10 .or. abs(x(3) - 1.0d0) > 1.0d-10) then
            write(*,'(A)') 'FAIL (wrong state)'
            n_failed = n_failed + 1
            return
        end if

        ! Expected Phi (row-major): [1,0,dt,0, 0,1,0,dt, 0,0,1,0, 0,0,0,1]
        phi_expected = [1.0d0, 0.0d0, dt_val, 0.0d0, &
                        0.0d0, 1.0d0, 0.0d0,  dt_val, &
                        0.0d0, 0.0d0, 1.0d0,  0.0d0, &
                        0.0d0, 0.0d0, 0.0d0,  1.0d0]

        do i = 1, 16
            if (abs(phi(i) - phi_expected(i)) > 1.0d-10) then
                write(*,'(A,I0,A,F12.6,A,F12.6)') 'FAIL (phi(', i, ')=', phi(i), &
                    ' expected ', phi_expected(i)
                n_failed = n_failed + 1
                return
            end if
        end do

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 5: Batch propagation — 5 states with const-vel, dt=1
    ! Each state has different initial position but same velocity [1,0]
    ! =====================================================================
    subroutine test_batch_const_vel()
        integer, parameter :: ns = 5, n = 4
        real(c_double) :: x_batch(n * ns), u(1), params(1)
        integer(c_int) :: info_batch(ns)
        type(c_funptr) :: null_ptr
        integer :: i
        real(c_double) :: expected_x1

        write(*,'(A)', advance='no') '  batch const_vel       ... '

        u(1) = 0.0d0
        params(1) = 0.0d0
        null_ptr = c_null_funptr

        ! Initialize 5 states: pos_x = i, pos_y = 0, vx = 1, vy = 0
        do i = 1, ns
            x_batch((i-1)*n + 1) = dble(i)    ! pos_x
            x_batch((i-1)*n + 2) = 0.0d0      ! pos_y
            x_batch((i-1)*n + 3) = 1.0d0      ! vx
            x_batch((i-1)*n + 4) = 0.0d0      ! vy
        end do

        call fa_propagate_batch(n, ns, x_batch, u, 1, null_ptr, 40, params, 1, 1.0d0, 1, info_batch)

        ! Check all info codes are 0
        do i = 1, ns
            if (info_batch(i) /= 0) then
                write(*,'(A,I0,A,I0)') 'FAIL (info_batch(', i, ')=', info_batch(i)
                n_failed = n_failed + 1
                return
            end if
        end do

        ! Check each state propagated correctly: pos_x should be i + 1
        do i = 1, ns
            expected_x1 = dble(i) + 1.0d0
            if (abs(x_batch((i-1)*n + 1) - expected_x1) > 1.0d-10 .or. &
                abs(x_batch((i-1)*n + 3) - 1.0d0) > 1.0d-10) then
                write(*,'(A,I0,A)') 'FAIL (state ', i, ' wrong)'
                n_failed = n_failed + 1
                return
            end if
        end do

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 6a: Edge case — n <= 0 returns info=3
    ! =====================================================================
    subroutine test_invalid_n()
        real(c_double) :: x(1), u(1), params(1)
        integer(c_int) :: info
        type(c_funptr) :: null_ptr

        write(*,'(A)', advance='no') '  invalid n (n=0)       ... '

        x(1) = 0.0d0
        u(1) = 0.0d0
        params(1) = 0.0d0
        null_ptr = c_null_funptr

        call fa_propagate(0, x, u, 1, null_ptr, 40, params, 1, 1.0d0, 1, info)

        if (info /= 3) then
            write(*,'(A,I0)') 'FAIL (expected info=3, got ', info
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 6b: Edge case — dt=0 leaves x unchanged
    ! =====================================================================
    subroutine test_dt_zero()
        real(c_double) :: x(4), u(1), params(1)
        integer(c_int) :: info
        type(c_funptr) :: null_ptr

        write(*,'(A)', advance='no') '  dt=0 no-op            ... '

        x = [5.0d0, 3.0d0, 1.0d0, 2.0d0]
        u(1) = 0.0d0
        params(1) = 0.0d0
        null_ptr = c_null_funptr

        call fa_propagate(4, x, u, 1, null_ptr, 40, params, 1, 0.0d0, 10, info)

        if (info /= 0) then
            write(*,'(A,I0)') 'FAIL (info=', info
            n_failed = n_failed + 1
            return
        end if

        if (abs(x(1) - 5.0d0) > 1.0d-15 .or. abs(x(2) - 3.0d0) > 1.0d-15 .or. &
            abs(x(3) - 1.0d0) > 1.0d-15 .or. abs(x(4) - 2.0d0) > 1.0d-15) then
            write(*,'(A)') 'FAIL (state changed)'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

end program test_propagate
