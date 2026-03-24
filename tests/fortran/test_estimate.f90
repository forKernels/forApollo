! Copyright The Fantastic Planet — By David Clabaugh
!
! test_estimate.f90 — Unit tests for forapollo_estimate (KF, EKF)
!
! Compile:
!   gfortran -O3 -std=f2008 -fall-intrinsics tests/fortran/test_estimate.f90 \
!            build/obj/forapollo_dynamics.o build/obj/forapollo_observe.o \
!            build/obj/forapollo_propagate.o build/obj/forapollo_estimate.o \
!            -o build/test_estimate -fopenmp
!
! Uses error stop on any failure.

program test_estimate
    use iso_c_binding
    implicit none

    ! -------------------------------------------------------------------------
    ! Interface blocks for bind(C) subroutines
    ! -------------------------------------------------------------------------
    interface
        subroutine fa_kf_predict(n, x, P, F, Q, info) bind(C, name="fa_kf_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            real(c_double), intent(in)         :: F(n*n)
            real(c_double), intent(in)         :: Q(n*n)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_kf_update(n, m, x, P, z, H, R, validity, info) bind(C, name="fa_kf_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            real(c_double), intent(in)         :: z(m)
            real(c_double), intent(in)         :: H(m*n)
            real(c_double), intent(in)         :: R(m*m)
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_ekf_predict(n, x, P, f_ptr, df_ptr, Q, dt, model_id, params, np, n_steps, info) &
            bind(C, name="fa_ekf_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            type(c_funptr), value              :: f_ptr
            type(c_funptr), value              :: df_ptr
            real(c_double), intent(in)         :: Q(n*n)
            real(c_double), value, intent(in)  :: dt
            real(c_double), intent(in)         :: params(np)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_ekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info) &
            bind(C, name="fa_ekf_update")
            use iso_c_binding
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
        end subroutine
    end interface

    integer :: n_passed, n_failed

    n_passed = 0
    n_failed = 0

    write(*,'(A)') '=== forApollo estimation tests ==='
    write(*,*)

    ! --- Linear Kalman filter tests ---
    call test_kf_predict_constvel()
    call test_kf_update_position()
    call test_kf_convergence()

    ! --- Extended Kalman filter tests ---
    call test_ekf_predict_constvel()
    call test_ekf_update_range()
    call test_ekf_tracking_10steps()
    call test_ternary_gating()

    ! --- Edge cases ---
    call test_edge_cases()

    ! --- Summary ---
    write(*,*)
    write(*,'(A,I3,A,I3,A)') '=== Results: ', n_passed, ' passed, ', n_failed, ' failed ==='
    if (n_failed > 0) error stop 'ESTIMATION TESTS FAILED'
    write(*,'(A)') 'All estimation tests passed.'

contains

    ! =====================================================================
    ! Helper: build const-vel F matrix for 2D (n=4)
    ! State: [px, py, vx, vy], F = [I dt*I; 0 I] row-major
    ! =====================================================================
    subroutine build_constvel_F(F, dt)
        real(c_double), intent(out) :: F(16)
        real(c_double), intent(in)  :: dt

        F = 0.0d0
        ! Row 1: [1, 0, dt, 0]
        F(1)  = 1.0d0; F(3)  = dt
        ! Row 2: [0, 1, 0, dt]
        F(6)  = 1.0d0; F(8)  = dt
        ! Row 3: [0, 0, 1, 0]
        F(11) = 1.0d0
        ! Row 4: [0, 0, 0, 1]
        F(16) = 1.0d0
    end subroutine

    ! =====================================================================
    ! Helper: build identity matrix (flat row-major)
    ! =====================================================================
    subroutine build_identity(I_mat, n)
        integer, intent(in) :: n
        real(c_double), intent(out) :: I_mat(n*n)
        integer :: i
        I_mat = 0.0d0
        do i = 1, n
            I_mat((i-1)*n + i) = 1.0d0
        end do
    end subroutine

    ! =====================================================================
    ! Test 1: KF predict — const-vel 2D
    ! =====================================================================
    subroutine test_kf_predict_constvel()
        real(c_double) :: x(4), P(16), F(16), Q(16)
        integer(c_int) :: info

        write(*,'(A)', advance='no') '  KF predict const-vel... '

        ! Initial state: pos=(0,0), vel=(1,0)
        x = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P, 4)
        call build_constvel_F(F, 1.0d0)
        Q = 0.0d0  ! No process noise

        call fa_kf_predict(4, x, P, F, Q, info)

        ! After predict: x should be (1, 0, 1, 0)
        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        if (abs(x(1) - 1.0d0) > 1.0d-12 .or. abs(x(2)) > 1.0d-12 .or. &
            abs(x(3) - 1.0d0) > 1.0d-12 .or. abs(x(4)) > 1.0d-12) then
            write(*,'(A)') 'FAIL (state mismatch)'
            write(*,'(A,4F10.4)') '    x = ', x
            n_failed = n_failed + 1
            return
        end if

        ! P should have grown (off-diagonal terms from F*P*F^T)
        ! P_new = F*I*F^T = F*F^T. Check P(1,1) = 1 + dt^2 = 2
        if (abs(P(1) - 2.0d0) > 1.0d-12) then
            write(*,'(A,F10.4)') 'FAIL (P(1,1)=', P(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 2: KF update — position observation
    ! =====================================================================
    subroutine test_kf_update_position()
        real(c_double) :: x(4), P(16), z(2), H(8), R(4)
        integer(c_int) :: validity(2), info

        write(*,'(A)', advance='no') '  KF update position obs... '

        ! State: pos=(0,0), vel=(1,0), P = I
        x = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P, 4)

        ! Observe position: H = [1 0 0 0; 0 1 0 0]
        H = 0.0d0
        H(1) = 1.0d0   ! H(1,1)
        H(6) = 1.0d0   ! H(2,2)

        ! Measurement: z = (1.5, 0.1), R = 0.5*I
        z = (/ 1.5d0, 0.1d0 /)
        R = 0.0d0
        R(1) = 0.5d0   ! R(1,1)
        R(4) = 0.5d0   ! R(2,2)

        call fa_kf_update(4, 2, x, P, z, H, R, validity, info)

        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! All measurements should be valid (innovation is small)
        if (validity(1) /= 1 .or. validity(2) /= 1) then
            write(*,'(A,2I3)') 'FAIL (validity=', validity
            n_failed = n_failed + 1
            return
        end if

        ! x should move toward z. With P=I, R=0.5*I:
        ! K = P*H^T*(H*P*H^T+R)^-1 = I2*1/(1+0.5) = 2/3 for position components
        ! x_pos = 0 + 2/3 * 1.5 = 1.0
        ! x_py  = 0 + 2/3 * 0.1 = 0.0667
        if (abs(x(1) - 1.0d0) > 1.0d-10 .or. abs(x(2) - 0.1d0*2.0d0/3.0d0) > 1.0d-10) then
            write(*,'(A)') 'FAIL (state did not move toward measurement)'
            write(*,'(A,4F10.4)') '    x = ', x
            n_failed = n_failed + 1
            return
        end if

        ! P should decrease (we gained information)
        if (P(1) >= 1.0d0) then
            write(*,'(A,F10.4)') 'FAIL (P(1,1) did not decrease: ', P(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 3: KF convergence — 20 predict-update cycles
    ! =====================================================================
    subroutine test_kf_convergence()
        real(c_double) :: x(4), P(16), F(16), Q(16), z(2), H(8), R(4)
        integer(c_int) :: validity(2), info
        real(c_double) :: P_diag_initial, P_diag_final
        integer :: step

        write(*,'(A)', advance='no') '  KF convergence 20 steps... '

        ! Initial state: pos=(0,0), vel=(1,0)
        x = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P, 4)
        P = P * 10.0d0  ! Large initial uncertainty

        call build_constvel_F(F, 1.0d0)

        ! Small process noise
        Q = 0.0d0
        Q(1)  = 0.01d0; Q(6)  = 0.01d0
        Q(11) = 0.01d0; Q(16) = 0.01d0

        ! Position observation
        H = 0.0d0
        H(1) = 1.0d0; H(6) = 1.0d0
        R = 0.0d0
        R(1) = 1.0d0; R(4) = 1.0d0

        P_diag_initial = P(1)

        do step = 1, 20
            call fa_kf_predict(4, x, P, F, Q, info)
            if (info /= 0) then
                write(*,'(A,I3,A,I3)') 'FAIL (predict info=', info, ' step=', step; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if

            ! Measurement near true position (truth: pos = step * vel)
            z(1) = dble(step) * 1.0d0 + 0.1d0  ! Slightly noisy
            z(2) = 0.05d0

            call fa_kf_update(4, 2, x, P, z, H, R, validity, info)
            if (info /= 0) then
                write(*,'(A,I3,A,I3)') 'FAIL (update info=', info, ' step=', step; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if
        end do

        P_diag_final = P(1)

        ! Covariance should have decreased significantly
        if (P_diag_final >= P_diag_initial) then
            write(*,'(A,F10.4,A,F10.4)') 'FAIL (P did not converge: ', P_diag_initial, ' -> ', P_diag_final
            n_failed = n_failed + 1
            return
        end if

        ! Position uncertainty should be small after 20 updates
        if (P_diag_final > 2.0d0) then
            write(*,'(A,F10.4)') 'FAIL (P(1,1) still large: ', P_diag_final; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A,F8.4,A,F8.4,A)') 'PASS (P(1,1): ', P_diag_initial, ' -> ', P_diag_final, ')'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 4: EKF predict — const-vel using model_id=40
    ! =====================================================================
    subroutine test_ekf_predict_constvel()
        real(c_double) :: x(4), P(16), Q(16), params(1)
        integer(c_int) :: info

        write(*,'(A)', advance='no') '  EKF predict const-vel (model 40)... '

        ! const-vel 2D: model_id=40, state=[px,py,vx,vy], no params needed
        x = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P, 4)
        Q = 0.0d0
        params(1) = 0.0d0

        call fa_ekf_predict(4, x, P, c_null_funptr, c_null_funptr, Q, 1.0d0, 40, params, 1, 10, info)

        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! After 1s of const-vel: pos=(1,0), vel=(1,0)
        if (abs(x(1) - 1.0d0) > 1.0d-8 .or. abs(x(2)) > 1.0d-8 .or. &
            abs(x(3) - 1.0d0) > 1.0d-8 .or. abs(x(4)) > 1.0d-8) then
            write(*,'(A)') 'FAIL (state mismatch)'
            write(*,'(A,4F10.6)') '    x = ', x
            n_failed = n_failed + 1
            return
        end if

        ! Covariance should have grown (P = Phi*P*Phi^T, Phi has off-diag)
        if (P(1) <= 1.0d0) then
            write(*,'(A,F10.4)') 'FAIL (P(1,1) did not grow: ', P(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 5: EKF update — range observation
    ! =====================================================================
    subroutine test_ekf_update_range()
        real(c_double) :: x(6), P(36), z(1), R(1), obs_params(3)
        integer(c_int) :: validity(1), info
        real(c_double) :: x_before

        write(*,'(A)', advance='no') '  EKF update range obs... '

        ! 6-state: [px,py,pz,vx,vy,vz]
        x = (/ 100.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        call build_identity(P, 6)
        P = P * 100.0d0  ! Large initial uncertainty

        ! Range observation from station at origin
        ! obs_id = 2 (range), obs_params = station position
        obs_params = (/ 0.0d0, 0.0d0, 0.0d0 /)

        ! Measured range = 105 (truth is 100, measurement says farther)
        z(1) = 105.0d0
        R(1) = 1.0d0

        x_before = x(1)

        call fa_ekf_update(6, 1, x, P, z, c_null_funptr, c_null_funptr, 2, R, obs_params, 3, validity, info)

        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        if (validity(1) /= 1) then
            write(*,'(A,I3)') 'FAIL (validity=', validity(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! x(1) should move toward 105 (measurement says farther)
        if (x(1) <= x_before) then
            write(*,'(A,F10.4)') 'FAIL (x(1) did not increase: ', x(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! P should decrease (gained information from measurement)
        if (P(1) >= 100.0d0) then
            write(*,'(A,F10.4)') 'FAIL (P(1,1) did not decrease: ', P(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A,F8.2,A)') 'PASS (x(1) moved to ', x(1), ')'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 6: EKF tracking — 10 steps const-vel + position obs
    ! =====================================================================
    subroutine test_ekf_tracking_10steps()
        real(c_double) :: x(4), P(16), Q(16), params(1)
        real(c_double) :: z(2), R(4), obs_params(1)
        integer(c_int) :: validity(2), info
        real(c_double) :: truth_px, truth_py, err_initial, err_final
        integer :: step

        write(*,'(A)', advance='no') '  EKF tracking 10 steps... '

        ! Initial estimate: off by 5 in position
        x = (/ 5.0d0, 5.0d0, 1.0d0, 0.5d0 /)
        call build_identity(P, 4)
        P = P * 25.0d0  ! Large uncertainty

        ! Process noise
        Q = 0.0d0
        Q(1) = 0.01d0; Q(6) = 0.01d0; Q(11) = 0.01d0; Q(16) = 0.01d0

        ! Observation: position (obs_id=1)
        R = 0.0d0
        R(1) = 1.0d0; R(4) = 1.0d0
        obs_params(1) = 0.0d0
        params(1) = 0.0d0

        ! Truth: starts at (0,0) with vel (1, 0.5)
        truth_px = 0.0d0
        truth_py = 0.0d0

        err_initial = sqrt((x(1) - truth_px)**2 + (x(2) - truth_py)**2)

        do step = 1, 10
            ! Predict
            call fa_ekf_predict(4, x, P, c_null_funptr, c_null_funptr, Q, 1.0d0, 40, params, 1, 10, info)
            if (info /= 0) then
                write(*,'(A,I3,A,I3)') 'FAIL (predict info=', info, ' step=', step; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if

            ! Truth propagation
            truth_px = truth_px + 1.0d0
            truth_py = truth_py + 0.5d0

            ! Noisy measurement of truth position
            z(1) = truth_px + 0.1d0
            z(2) = truth_py - 0.05d0

            ! Update with position observation (obs_id=1, m=3 for 3D but we need m=2)
            ! Use linear KF update since position observation is linear
            block
                real(c_double) :: H(8)
                H = 0.0d0
                H(1) = 1.0d0   ! H(1,1)
                H(6) = 1.0d0   ! H(2,2)
                call fa_kf_update(4, 2, x, P, z, H, R, validity, info)
            end block
            if (info /= 0) then
                write(*,'(A,I3,A,I3)') 'FAIL (update info=', info, ' step=', step; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if
        end do

        err_final = sqrt((x(1) - truth_px)**2 + (x(2) - truth_py)**2)

        ! Error should decrease significantly
        if (err_final >= err_initial) then
            write(*,'(A,F10.4,A,F10.4)') 'FAIL (error grew: ', err_initial, ' -> ', err_final
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A,F8.4,A,F8.4,A)') 'PASS (err: ', err_initial, ' -> ', err_final, ')'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 7: Ternary gating — reject outlier, accept normal
    ! =====================================================================
    subroutine test_ternary_gating()
        real(c_double) :: x(2), P(4), z(2), H(4), R(4)
        integer(c_int) :: validity(2), info
        real(c_double) :: x_save(2)

        write(*,'(A)', advance='no') '  Ternary gating test... '

        ! Simple 2-state system
        x = (/ 0.0d0, 0.0d0 /)
        call build_identity(P, 2)

        ! H = I(2)
        call build_identity(H, 2)

        ! R = I(2)
        call build_identity(R, 2)

        ! z(1) = huge outlier (innovation >> 4 sigma), z(2) = normal
        ! S(i,i) = P(i,i) + R(i,i) = 1 + 1 = 2
        ! For rejection: chi2 = y^2/S > 16 → |y| > sqrt(32) ≈ 5.66
        ! For uncertain: chi2 = y^2/S > 9  → |y| > sqrt(18) ≈ 4.24
        z(1) = 100.0d0   ! Huge outlier: chi2 = 10000/2 = 5000 >> 16
        z(2) = 0.5d0     ! Normal: chi2 = 0.25/2 = 0.125 < 9

        x_save = x

        call fa_kf_update(2, 2, x, P, z, H, R, validity, info)

        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! Measurement 1 should be rejected
        if (validity(1) /= -1) then
            write(*,'(A,I3)') 'FAIL (outlier not rejected, validity(1)=', validity(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! Measurement 2 should be accepted
        if (validity(2) /= 1) then
            write(*,'(A,I3)') 'FAIL (normal not accepted, validity(2)=', validity(2); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! x(2) should have moved toward z(2)=0.5 but x(1) should be less affected
        ! Since only valid measurements are fused, x(1) should stay near 0
        ! and x(2) should move toward 0.5
        if (abs(x(2)) < 0.1d0) then
            write(*,'(A,F10.4)') 'FAIL (x(2) did not move toward measurement: ', x(2); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A,2I3,A)') 'PASS (validity=', validity, ')'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 8: Edge cases — n=0, m=0
    ! =====================================================================
    subroutine test_edge_cases()
        real(c_double) :: x(1), P(1), F(1), Q(1), z(1), H(1), R(1)
        integer(c_int) :: validity(1), info

        write(*,'(A)', advance='no') '  Edge cases (n=0, m=0)... '

        ! KF predict with n=0
        call fa_kf_predict(0, x, P, F, Q, info)
        if (info /= 3) then
            write(*,'(A,I3)') 'FAIL (kf_predict n=0 info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! KF update with n=0
        call fa_kf_update(0, 1, x, P, z, H, R, validity, info)
        if (info /= 3) then
            write(*,'(A,I3)') 'FAIL (kf_update n=0 info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! KF update with m=0
        call fa_kf_update(1, 0, x, P, z, H, R, validity, info)
        if (info /= 3) then
            write(*,'(A,I3)') 'FAIL (kf_update m=0 info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! EKF predict with n=0
        call fa_ekf_predict(0, x, P, c_null_funptr, c_null_funptr, Q, 1.0d0, 40, F, 1, 1, info)
        if (info /= 3) then
            write(*,'(A,I3)') 'FAIL (ekf_predict n=0 info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! EKF update with m=0
        call fa_ekf_update(1, 0, x, P, z, c_null_funptr, c_null_funptr, 1, R, Q, 1, validity, info)
        if (info /= 3) then
            write(*,'(A,I3)') 'FAIL (ekf_update m=0 info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

end program test_estimate
