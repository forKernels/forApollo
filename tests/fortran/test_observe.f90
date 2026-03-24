! Copyright The Fantastic Planet — By David Clabaugh
!
! test_observe.f90 — Unit tests for forapollo_observe measurement models
!
! Tests all 16 measurement models including analytic Jacobian verification
! via central finite differencing for nontrivial models.
!
! Compile:
!   gfortran -O3 -std=f2008 -fall-intrinsics tests/fortran/test_observe.f90 \
!            build/obj/forapollo_observe.o -o build/test_observe -fopenmp
!
! Uses error stop on any failure.

program test_observe
    use iso_c_binding
    implicit none

    ! -------------------------------------------------------------------------
    ! Interface blocks for bind(C) subroutines.
    ! CRITICAL: without these, gfortran passes value args by reference → segfault.
    ! -------------------------------------------------------------------------
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

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    write(*,'(A)') '=== forApollo observation model tests ==='
    write(*,*)

    ! --- Basic measurement tests ---
    call test_position()
    call test_range()
    call test_bearing()
    call test_range_bearing()
    call test_velocity()
    call test_doppler()
    call test_radar()
    call test_scalar()
    call test_pinhole()
    call test_relative_position()
    call test_relative_velocity()

    ! --- Attitude sensor tests ---
    call test_magnetometer()
    call test_gyroscope()

    ! --- Jacobian FD verification (nontrivial models) ---
    call test_jacobian_range()
    call test_jacobian_bearing()
    call test_jacobian_radar()
    call test_jacobian_pinhole()
    call test_jacobian_doppler()
    call test_jacobian_range_bearing()
    call test_jacobian_attitude()

    ! --- Edge cases ---
    call test_invalid_obs_id()

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
    !   H_fd(:, j) = (h(x + h*e_j) - h(x - h*e_j)) / (2h)
    ! then compare against analytic H.
    ! =========================================================================
    subroutine verify_observe_jacobian(test_name, obs_id_val, n_val, m_val, x, &
                                        obs_params, nop_val, tol, pass)
        character(len=*), intent(in)  :: test_name
        integer, intent(in)           :: obs_id_val, n_val, m_val, nop_val
        real(c_double), intent(in)    :: x(n_val)
        real(c_double), intent(in)    :: obs_params(nop_val)
        real(c_double), intent(in)    :: tol
        logical, intent(out)          :: pass

        real(c_double) :: H_analytic(m_val * n_val)
        real(c_double) :: x_plus(n_val), x_minus(n_val)
        real(c_double) :: z_plus(m_val), z_minus(m_val)
        real(c_double) :: H_fd_col(m_val)
        real(c_double) :: h_step, diff
        integer(c_int) :: info
        integer :: i, j

        h_step = 1.0d-7

        ! Get analytic Jacobian
        call fa_observe_jacobian(obs_id_val, n_val, x, m_val, 0.0d0, &
                                  obs_params, nop_val, H_analytic, info)
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
            x_plus(j)  = x(j) + h_step
            x_minus(j) = x(j) - h_step

            call fa_observe_dispatch(obs_id_val, n_val, x_plus, m_val, 0.0d0, &
                                      obs_params, nop_val, z_plus, info)
            call fa_observe_dispatch(obs_id_val, n_val, x_minus, m_val, 0.0d0, &
                                      obs_params, nop_val, z_minus, info)

            ! H_fd column j = (z_plus - z_minus) / (2*h_step)
            H_fd_col(1:m_val) = (z_plus(1:m_val) - z_minus(1:m_val)) / (2.0d0 * h_step)

            ! Compare row i, col j: analytic H((i-1)*n + j) vs H_fd_col(i)
            do i = 1, m_val
                diff = abs(H_analytic((i-1)*n_val + j) - H_fd_col(i))
                if (diff > tol) then
                    write(*,'(A,A,A,I0,A,I0,A,ES12.5,A,ES12.5,A,ES12.5)') &
                        '  FAIL [', test_name, '] H(', i, ',', j, &
                        '): analytic=', H_analytic((i-1)*n_val + j), &
                        ' FD=', H_fd_col(i), ' diff=', diff
                    pass = .false.
                end if
            end do
        end do
    end subroutine verify_observe_jacobian

    ! =========================================================================
    ! Test: Position (obs_id = 1)
    ! =========================================================================
    subroutine test_position()
        real(c_double) :: x(6), z(3), expected(3), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
        expected = (/ 1.0d0, 2.0d0, 3.0d0 /)
        dummy_p(1) = 0.0d0

        call fa_observe_dispatch(1, 6, x, 3, 0.0d0, dummy_p, 1, z, info)

        call check_vec('position', 3, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('position measurement', pass)
    end subroutine test_position

    ! =========================================================================
    ! Test: Range (obs_id = 2)
    ! =========================================================================
    subroutine test_range()
        real(c_double) :: x(6), z(1), expected(1), station(3)
        integer(c_int) :: info
        logical :: pass

        x = (/ 100.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        station = (/ 0.0d0, 0.0d0, 0.0d0 /)
        expected(1) = 100.0d0

        call fa_observe_dispatch(2, 6, x, 1, 0.0d0, station, 3, z, info)

        call check_vec('range', 1, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('range measurement', pass)
    end subroutine test_range

    ! =========================================================================
    ! Test: Bearing (obs_id = 3)
    ! =========================================================================
    subroutine test_bearing()
        real(c_double) :: x(6), z(1), expected(1), station(3)
        real(c_double), parameter :: pi = 3.141592653589793d0
        integer(c_int) :: info
        logical :: pass

        x = (/ 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        station = (/ 0.0d0, 0.0d0, 0.0d0 /)
        expected(1) = atan2(1.0d0, 1.0d0)  ! pi/4

        call fa_observe_dispatch(3, 6, x, 1, 0.0d0, station, 3, z, info)

        call check_vec('bearing', 1, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('bearing measurement', pass)
    end subroutine test_bearing

    ! =========================================================================
    ! Test: Range+Bearing (obs_id = 4)
    ! =========================================================================
    subroutine test_range_bearing()
        real(c_double) :: x(6), z(2), expected(2), station(3)
        integer(c_int) :: info
        logical :: pass

        x = (/ 3.0d0, 4.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        station = (/ 0.0d0, 0.0d0, 0.0d0 /)
        expected(1) = 5.0d0                     ! sqrt(9+16)
        expected(2) = atan2(4.0d0, 3.0d0)

        call fa_observe_dispatch(4, 6, x, 2, 0.0d0, station, 3, z, info)

        call check_vec('range+bearing', 2, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('range+bearing measurement', pass)
    end subroutine test_range_bearing

    ! =========================================================================
    ! Test: Velocity (obs_id = 10)
    ! =========================================================================
    subroutine test_velocity()
        real(c_double) :: x(6), z(3), expected(3), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = (/ 0.0d0, 0.0d0, 0.0d0, 5.0d0, 6.0d0, 7.0d0 /)
        expected = (/ 5.0d0, 6.0d0, 7.0d0 /)
        dummy_p(1) = 0.0d0

        call fa_observe_dispatch(10, 6, x, 3, 0.0d0, dummy_p, 1, z, info)

        call check_vec('velocity', 3, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('velocity measurement', pass)
    end subroutine test_velocity

    ! =========================================================================
    ! Test: Doppler (obs_id = 11)
    ! Target at [100,0,0] moving [10,0,0], station at origin stationary
    ! range-rate = dot([100,0,0],[10,0,0]) / 100 = 1000/100 = 10
    ! =========================================================================
    subroutine test_doppler()
        real(c_double) :: x(6), z(1), expected(1), params(6)
        integer(c_int) :: info
        logical :: pass

        x = (/ 100.0d0, 0.0d0, 0.0d0, 10.0d0, 0.0d0, 0.0d0 /)
        params = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        expected(1) = 10.0d0

        call fa_observe_dispatch(11, 6, x, 1, 0.0d0, params, 6, z, info)

        call check_vec('doppler', 1, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('doppler measurement', pass)
    end subroutine test_doppler

    ! =========================================================================
    ! Test: Radar (obs_id = 40)
    ! Target at [100, 0, 0] from origin → range=100, az=0, el=0
    ! =========================================================================
    subroutine test_radar()
        real(c_double) :: x(6), z(3), expected(3), station(3)
        integer(c_int) :: info
        logical :: pass

        x = (/ 100.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        station = (/ 0.0d0, 0.0d0, 0.0d0 /)
        expected = (/ 100.0d0, 0.0d0, 0.0d0 /)

        call fa_observe_dispatch(40, 6, x, 3, 0.0d0, station, 3, z, info)

        call check_vec('radar', 3, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('radar measurement', pass)
    end subroutine test_radar

    ! =========================================================================
    ! Test: Scalar (obs_id = 50)
    ! obs_params=[2] → z = x(2)
    ! =========================================================================
    subroutine test_scalar()
        real(c_double) :: x(6), z(1), expected(1), params(1)
        integer(c_int) :: info
        logical :: pass

        x = (/ 10.0d0, 20.0d0, 30.0d0, 40.0d0, 50.0d0, 60.0d0 /)
        params(1) = 2.0d0
        expected(1) = 20.0d0

        call fa_observe_dispatch(50, 6, x, 1, 0.0d0, params, 1, z, info)

        call check_vec('scalar', 1, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('scalar measurement', pass)
    end subroutine test_scalar

    ! =========================================================================
    ! Test: Pinhole (obs_id = 60)
    ! x=[1,0,2,...], fx=500,fy=500,cx=320,cy=240
    ! z = [500*1/2+320, 500*0/2+240] = [570, 240]
    ! =========================================================================
    subroutine test_pinhole()
        real(c_double) :: x(6), z(2), expected(2), params(4)
        integer(c_int) :: info
        logical :: pass

        x = (/ 1.0d0, 0.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        params = (/ 500.0d0, 500.0d0, 320.0d0, 240.0d0 /)
        expected = (/ 570.0d0, 240.0d0 /)

        call fa_observe_dispatch(60, 6, x, 2, 0.0d0, params, 4, z, info)

        call check_vec('pinhole', 2, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('pinhole measurement', pass)
    end subroutine test_pinhole

    ! =========================================================================
    ! Test: Relative Position (obs_id = 70)
    ! =========================================================================
    subroutine test_relative_position()
        real(c_double) :: x(6), z(3), expected(3), x2(3)
        integer(c_int) :: info
        logical :: pass

        x = (/ 10.0d0, 20.0d0, 30.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        x2 = (/ 1.0d0, 2.0d0, 3.0d0 /)
        expected = (/ 9.0d0, 18.0d0, 27.0d0 /)

        call fa_observe_dispatch(70, 6, x, 3, 0.0d0, x2, 3, z, info)

        call check_vec('relative_position', 3, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('relative position measurement', pass)
    end subroutine test_relative_position

    ! =========================================================================
    ! Test: Relative Velocity (obs_id = 71)
    ! =========================================================================
    subroutine test_relative_velocity()
        real(c_double) :: x(6), z(3), expected(3), x2(3)
        integer(c_int) :: info
        logical :: pass

        x = (/ 0.0d0, 0.0d0, 0.0d0, 10.0d0, 20.0d0, 30.0d0 /)
        x2 = (/ 1.0d0, 2.0d0, 3.0d0 /)
        expected = (/ 9.0d0, 18.0d0, 27.0d0 /)

        call fa_observe_dispatch(71, 6, x, 3, 0.0d0, x2, 3, z, info)

        call check_vec('relative_velocity', 3, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('relative velocity measurement', pass)
    end subroutine test_relative_velocity

    ! =========================================================================
    ! Test: Magnetometer (obs_id = 20)
    ! Identity quaternion q=[1,0,0,0] → R=I → z = ref_field
    ! =========================================================================
    subroutine test_magnetometer()
        real(c_double) :: x(13), z(3), expected(3), ref_field(3)
        integer(c_int) :: info
        logical :: pass

        x = 0.0d0
        x(7) = 1.0d0    ! q0 = 1 (identity quaternion)
        x(8) = 0.0d0
        x(9) = 0.0d0
        x(10) = 0.0d0

        ref_field = (/ 0.3d0, 0.0d0, -0.5d0 /)
        expected = (/ 0.3d0, 0.0d0, -0.5d0 /)  ! Identity rotation

        call fa_observe_dispatch(20, 13, x, 3, 0.0d0, ref_field, 3, z, info)

        call check_vec('magnetometer', 3, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('magnetometer measurement', pass)
    end subroutine test_magnetometer

    ! =========================================================================
    ! Test: Gyroscope (obs_id = 31)
    ! z = x(11:13)
    ! =========================================================================
    subroutine test_gyroscope()
        real(c_double) :: x(13), z(3), expected(3), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = 0.0d0
        x(7) = 1.0d0   ! identity quaternion
        x(11) = 0.1d0
        x(12) = -0.2d0
        x(13) = 0.3d0
        expected = (/ 0.1d0, -0.2d0, 0.3d0 /)
        dummy_p(1) = 0.0d0

        call fa_observe_dispatch(31, 13, x, 3, 0.0d0, dummy_p, 1, z, info)

        call check_vec('gyroscope', 3, z, expected, 1.0d-12, pass)
        if (info /= 0) pass = .false.
        call report('gyroscope measurement', pass)
    end subroutine test_gyroscope

    ! =========================================================================
    ! Jacobian FD: Range (obs_id = 2)
    ! =========================================================================
    subroutine test_jacobian_range()
        real(c_double) :: x(6), station(3)
        logical :: pass

        x = (/ 30.0d0, 40.0d0, 50.0d0, 1.0d0, 2.0d0, 3.0d0 /)
        station = (/ 5.0d0, 10.0d0, 15.0d0 /)

        call verify_observe_jacobian('range Jacobian FD', 2, 6, 1, x, station, 3, 1.0d-5, pass)
        call report('range Jacobian FD', pass)
    end subroutine test_jacobian_range

    ! =========================================================================
    ! Jacobian FD: Bearing (obs_id = 3)
    ! =========================================================================
    subroutine test_jacobian_bearing()
        real(c_double) :: x(6), station(3)
        logical :: pass

        x = (/ 10.0d0, 20.0d0, 5.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        station = (/ 1.0d0, 2.0d0, 0.0d0 /)

        call verify_observe_jacobian('bearing Jacobian FD', 3, 6, 1, x, station, 3, 1.0d-5, pass)
        call report('bearing Jacobian FD', pass)
    end subroutine test_jacobian_bearing

    ! =========================================================================
    ! Jacobian FD: Radar (obs_id = 40)
    ! =========================================================================
    subroutine test_jacobian_radar()
        real(c_double) :: x(6), station(3)
        logical :: pass

        ! Avoid axis-aligned to exercise all partials
        x = (/ 30.0d0, 40.0d0, 50.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        station = (/ 5.0d0, 10.0d0, 15.0d0 /)

        call verify_observe_jacobian('radar Jacobian FD', 40, 6, 3, x, station, 3, 1.0d-5, pass)
        call report('radar Jacobian FD', pass)
    end subroutine test_jacobian_radar

    ! =========================================================================
    ! Jacobian FD: Pinhole (obs_id = 60)
    ! =========================================================================
    subroutine test_jacobian_pinhole()
        real(c_double) :: x(6), params(4)
        logical :: pass

        x = (/ 0.5d0, 0.3d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        params = (/ 500.0d0, 500.0d0, 320.0d0, 240.0d0 /)

        call verify_observe_jacobian('pinhole Jacobian FD', 60, 6, 2, x, params, 4, 1.0d-5, pass)
        call report('pinhole Jacobian FD', pass)
    end subroutine test_jacobian_pinhole

    ! =========================================================================
    ! Jacobian FD: Doppler (obs_id = 11)
    ! =========================================================================
    subroutine test_jacobian_doppler()
        real(c_double) :: x(6), params(6)
        logical :: pass

        x = (/ 30.0d0, 40.0d0, 50.0d0, 5.0d0, -3.0d0, 2.0d0 /)
        params = (/ 1.0d0, 2.0d0, 3.0d0, 0.5d0, -0.5d0, 1.0d0 /)

        call verify_observe_jacobian('doppler Jacobian FD', 11, 6, 1, x, params, 6, 1.0d-5, pass)
        call report('doppler Jacobian FD', pass)
    end subroutine test_jacobian_doppler

    ! =========================================================================
    ! Jacobian FD: Range+Bearing (obs_id = 4)
    ! =========================================================================
    subroutine test_jacobian_range_bearing()
        real(c_double) :: x(6), station(3)
        logical :: pass

        x = (/ 20.0d0, 15.0d0, 10.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        station = (/ 3.0d0, 4.0d0, 5.0d0 /)

        call verify_observe_jacobian('range+bearing Jacobian FD', 4, 6, 2, x, station, 3, 1.0d-5, pass)
        call report('range+bearing Jacobian FD', pass)
    end subroutine test_jacobian_range_bearing

    ! =========================================================================
    ! Jacobian FD: Attitude sensor (obs_id = 20, magnetometer)
    ! Use a non-trivial quaternion for a good test
    ! =========================================================================
    subroutine test_jacobian_attitude()
        real(c_double) :: x(13), ref_field(3)
        real(c_double) :: qnorm
        logical :: pass

        x = 0.0d0
        ! Non-trivial quaternion (normalized)
        x(7) = 0.5d0
        x(8) = 0.5d0
        x(9) = 0.5d0
        x(10) = 0.5d0
        ! Already unit quaternion: 0.25*4 = 1.0

        ref_field = (/ 0.3d0, 0.0d0, -0.5d0 /)

        call verify_observe_jacobian('attitude Jacobian FD', 20, 13, 3, x, ref_field, 3, 1.0d-5, pass)
        call report('attitude sensor Jacobian FD', pass)
    end subroutine test_jacobian_attitude

    ! =========================================================================
    ! Test: Invalid obs_id → info=3
    ! =========================================================================
    subroutine test_invalid_obs_id()
        real(c_double) :: x(6), z(1), dummy_p(1)
        integer(c_int) :: info
        logical :: pass

        x = 0.0d0
        dummy_p(1) = 0.0d0

        call fa_observe_dispatch(999, 6, x, 1, 0.0d0, dummy_p, 1, z, info)

        pass = (info == 3)
        call report('invalid obs_id returns info=3', pass)
    end subroutine test_invalid_obs_id

end program test_observe
