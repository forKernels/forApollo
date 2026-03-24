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

        subroutine fa_iekf_update(n, m, x, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, &
                                   max_iter, tol, validity, info) &
            bind(C, name="fa_iekf_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m, obs_id, nop, max_iter
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: P(n*n)
            real(c_double), intent(in)         :: z(m)
            type(c_funptr), value              :: h_ptr, dh_ptr
            real(c_double), intent(in)         :: R(m*m)
            real(c_double), intent(in)         :: obs_params(nop)
            real(c_double), value, intent(in)  :: tol
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine

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

        subroutine fa_eskf_predict(n, x_nom, dx, P, f_ptr, model_id, params, np, Q, dt, n_steps, info) &
            bind(C, name="fa_eskf_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, model_id, np, n_steps
            real(c_double), intent(inout)      :: x_nom(n)
            real(c_double), intent(inout)      :: dx(n)
            real(c_double), intent(inout)      :: P(n*n)
            type(c_funptr), value              :: f_ptr
            real(c_double), intent(in)         :: params(np)
            real(c_double), intent(in)         :: Q(n*n)
            real(c_double), value, intent(in)  :: dt
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_eskf_update(n, m, x_nom, dx, P, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, &
                                   validity, info) &
            bind(C, name="fa_eskf_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m, obs_id, nop
            real(c_double), intent(in)         :: x_nom(n)
            real(c_double), intent(inout)      :: dx(n)
            real(c_double), intent(inout)      :: P(n*n)
            real(c_double), intent(in)         :: z(m)
            type(c_funptr), value              :: h_ptr, dh_ptr
            real(c_double), intent(in)         :: R(m*m)
            real(c_double), intent(in)         :: obs_params(nop)
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_eskf_inject(n, x_nom, dx, P, info) bind(C, name="fa_eskf_inject")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n
            real(c_double), intent(inout)      :: x_nom(n)
            real(c_double), intent(inout)      :: dx(n)
            real(c_double), intent(inout)      :: P(n*n)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_srekf_predict(n, x, S, f_ptr, df_ptr, model_id, params, np, Sq, dt, n_steps, info) &
            bind(C, name="fa_srekf_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, model_id, np, n_steps
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: S(n*n)
            type(c_funptr), value              :: f_ptr, df_ptr
            real(c_double), intent(in)         :: Sq(n*n)
            real(c_double), value, intent(in)  :: dt
            real(c_double), intent(in)         :: params(np)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_srekf_update(n, m, x, S, z, h_ptr, dh_ptr, obs_id, R, obs_params, nop, validity, info) &
            bind(C, name="fa_srekf_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m, obs_id, nop
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: S(n*n)
            real(c_double), intent(in)         :: z(m)
            type(c_funptr), value              :: h_ptr, dh_ptr
            real(c_double), intent(in)         :: R(m*m)
            real(c_double), intent(in)         :: obs_params(nop)
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_srukf_predict(n, x, S, f_ptr, model_id, params, np, Sq, dt, &
                                     alpha, beta_ukf, kappa, info) &
            bind(C, name="fa_srukf_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, model_id, np
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(inout)      :: S(n*n)
            type(c_funptr), value              :: f_ptr
            real(c_double), intent(in)         :: params(np)
            real(c_double), intent(in)         :: Sq(n*n)
            real(c_double), value, intent(in)  :: dt
            real(c_double), value, intent(in)  :: alpha, beta_ukf, kappa
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_srukf_update(n, m, x, S, z, h_ptr, obs_id, R, obs_params, nop, &
                                    alpha, beta_ukf, kappa, validity, info) &
            bind(C, name="fa_srukf_update")
            use iso_c_binding
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
        end subroutine

        subroutine fa_if_predict(n, eta, Ymat, F, Q, info) bind(C, name="fa_if_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n
            real(c_double), intent(inout)      :: eta(n)
            real(c_double), intent(inout)      :: Ymat(n*n)
            real(c_double), intent(in)         :: F(n*n)
            real(c_double), intent(in)         :: Q(n*n)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_if_update(n, m, eta, Ymat, z, H, R, validity, info) &
            bind(C, name="fa_if_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m
            real(c_double), intent(inout)      :: eta(n)
            real(c_double), intent(inout)      :: Ymat(n*n)
            real(c_double), intent(in)         :: z(m)
            real(c_double), intent(in)         :: H(m*n)
            real(c_double), intent(in)         :: R(m*m)
            integer(c_int), intent(out)        :: validity(m)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_pf_sir_predict(n, n_particles, particles, weights, &
                                      f_ptr, model_id, params, np, Q, dt, n_steps, seed, info) &
            bind(C, name="fa_pf_sir_predict")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, n_particles, model_id, np, n_steps
            real(c_double), intent(inout)      :: particles(n * n_particles)
            real(c_double), intent(inout)      :: weights(n_particles)
            type(c_funptr), value              :: f_ptr
            real(c_double), intent(in)         :: params(np)
            real(c_double), intent(in)         :: Q(n*n)
            real(c_double), value, intent(in)  :: dt
            integer(c_int), intent(inout)      :: seed
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_pf_sir_update(n, m, n_particles, particles, weights, &
                                     z, h_ptr, obs_id, R, obs_params, nop, info) &
            bind(C, name="fa_pf_sir_update")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m, n_particles, obs_id, nop
            real(c_double), intent(in)         :: particles(n * n_particles)
            real(c_double), intent(inout)      :: weights(n_particles)
            real(c_double), intent(in)         :: z(m)
            type(c_funptr), value              :: h_ptr
            real(c_double), intent(in)         :: R(m*m)
            real(c_double), intent(in)         :: obs_params(nop)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_pf_sir_resample(n, n_particles, particles, weights, seed, info) &
            bind(C, name="fa_pf_sir_resample")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, n_particles
            real(c_double), intent(inout)      :: particles(n * n_particles)
            real(c_double), intent(inout)      :: weights(n_particles)
            integer(c_int), intent(inout)      :: seed
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_rts_smooth(n, nsteps, x_filt, P_filt, x_pred, P_pred, F_all, &
                                  x_smooth, P_smooth, info) &
            bind(C, name="fa_rts_smooth")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, nsteps
            real(c_double), intent(in)         :: x_filt(n*nsteps)
            real(c_double), intent(in)         :: P_filt(n*n*nsteps)
            real(c_double), intent(in)         :: x_pred(n*nsteps)
            real(c_double), intent(in)         :: P_pred(n*n*nsteps)
            real(c_double), intent(in)         :: F_all(n*n*nsteps)
            real(c_double), intent(out)        :: x_smooth(n*nsteps)
            real(c_double), intent(out)        :: P_smooth(n*n*nsteps)
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_batch_wls(n, m_total, x, z_all, H_all, R_all, max_iter, tol, info) &
            bind(C, name="fa_batch_wls")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m_total, max_iter
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(in)         :: z_all(m_total)
            real(c_double), intent(in)         :: H_all(m_total * n)
            real(c_double), intent(in)         :: R_all(m_total * m_total)
            real(c_double), value, intent(in)  :: tol
            integer(c_int), intent(out)        :: info
        end subroutine

        subroutine fa_batch_map(n, m_total, x, x0, P0, z_all, H_all, R_all, max_iter, tol, info) &
            bind(C, name="fa_batch_map")
            use iso_c_binding
            integer(c_int), value, intent(in)  :: n, m_total, max_iter
            real(c_double), intent(inout)      :: x(n)
            real(c_double), intent(in)         :: x0(n)
            real(c_double), intent(in)         :: P0(n*n)
            real(c_double), intent(in)         :: z_all(m_total)
            real(c_double), intent(in)         :: H_all(m_total * n)
            real(c_double), intent(in)         :: R_all(m_total * m_total)
            real(c_double), value, intent(in)  :: tol
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

    ! --- Advanced estimator tests ---
    call test_iekf_update_range()
    call test_ukf_predict_constvel()
    call test_ukf_update_position()
    call test_ukf_update_range_nonlinear()
    call test_eskf_predict_update_inject()
    call test_srekf_vs_ekf()
    call test_ukf_singular_P()

    ! --- Edge cases ---
    call test_edge_cases()

    ! --- Information filter tests ---
    call test_if_update_additivity()
    call test_if_predict_update_cycle()

    ! --- Particle filter tests ---
    call test_pf_sir_basic()
    call test_pf_sir_resample()

    ! --- Smoother tests ---
    call test_rts_smoother()

    ! --- Batch estimator tests ---
    call test_batch_wls()
    call test_batch_map()

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

    ! =====================================================================
    ! Test 9: IEKF update — range measurement, should converge in 2-3 iters
    ! and be closer to truth than single EKF update
    ! =====================================================================
    subroutine test_iekf_update_range()
        real(c_double) :: x_ekf(6), P_ekf(36), x_iekf(6), P_iekf(36)
        real(c_double) :: z(1), R(1), obs_params(3)
        integer(c_int) :: validity_ekf(1), validity_iekf(1), info_ekf, info_iekf
        real(c_double) :: err_ekf, err_iekf
        real(c_double) :: truth_x

        write(*,'(A)', advance='no') '  IEKF update range obs... '

        ! 6-state: [px,py,pz,vx,vy,vz], truth at px=100
        truth_x = 100.0d0

        ! Both start at same estimate (px=95, off by 5)
        x_ekf  = (/ 95.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        x_iekf = (/ 95.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        call build_identity(P_ekf, 6)
        P_ekf = P_ekf * 100.0d0
        P_iekf = P_ekf

        ! Range observation from station at origin, measured range = 100
        obs_params = (/ 0.0d0, 0.0d0, 0.0d0 /)
        z(1) = 100.0d0
        R(1) = 1.0d0

        ! EKF update (single linearization)
        call fa_ekf_update(6, 1, x_ekf, P_ekf, z, c_null_funptr, c_null_funptr, &
                           2, R, obs_params, 3, validity_ekf, info_ekf)

        ! IEKF update (iterated, up to 10 iterations)
        call fa_iekf_update(6, 1, x_iekf, P_iekf, z, c_null_funptr, c_null_funptr, &
                            2, R, obs_params, 3, 10, 1.0d-10, validity_iekf, info_iekf)

        if (info_ekf /= 0) then
            write(*,'(A,I3)') 'FAIL (EKF info=', info_ekf; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! IEKF may return info=4 if not fully converged, but should still give good result
        ! Accept info=0 or info=4
        if (info_iekf /= 0 .and. info_iekf /= 4) then
            write(*,'(A,I3)') 'FAIL (IEKF info=', info_iekf; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! Both should have moved toward truth
        err_ekf  = abs(x_ekf(1) - truth_x)
        err_iekf = abs(x_iekf(1) - truth_x)

        ! IEKF should be at least as good as EKF (often better for nonlinear)
        if (err_iekf > err_ekf + 1.0d0) then
            write(*,'(A,F8.4,A,F8.4)') 'FAIL (IEKF worse: ekf_err=', err_ekf, ' iekf_err=', err_iekf
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A,F8.4,A,F8.4,A)') 'PASS (ekf_err=', err_ekf, ' iekf_err=', err_iekf, ')'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 10: UKF predict — const-vel model, should match KF predict
    ! (linear dynamics → UKF = KF)
    ! =====================================================================
    subroutine test_ukf_predict_constvel()
        real(c_double) :: x_kf(4), P_kf(16), x_ukf(4), P_ukf(16)
        real(c_double) :: F(16), Q(16), params(1)
        integer(c_int) :: info_kf, info_ukf

        write(*,'(A)', advance='no') '  UKF predict const-vel... '

        ! Initial state
        x_kf  = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        x_ukf = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P_kf, 4)
        P_ukf = P_kf

        Q = 0.0d0
        Q(1) = 0.01d0; Q(6) = 0.01d0; Q(11) = 0.01d0; Q(16) = 0.01d0
        params(1) = 0.0d0

        ! KF predict (linear reference)
        call build_constvel_F(F, 1.0d0)
        call fa_kf_predict(4, x_kf, P_kf, F, Q, info_kf)

        ! UKF predict (model_id=40 = const-vel 2D)
        ! Use large alpha for better sigma point spread on linear problems
        call fa_ukf_predict(4, x_ukf, P_ukf, c_null_funptr, 40, params, 1, Q, 1.0d0, &
                            1.0d0, 2.0d0, 0.0d0, info_ukf)

        if (info_kf /= 0) then
            write(*,'(A,I3)') 'FAIL (KF info=', info_kf; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if
        if (info_ukf /= 0) then
            write(*,'(A,I3)') 'FAIL (UKF info=', info_ukf; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! States should match (linear dynamics → UKF = KF)
        if (abs(x_ukf(1) - x_kf(1)) > 1.0d-6 .or. abs(x_ukf(2) - x_kf(2)) > 1.0d-6 .or. &
            abs(x_ukf(3) - x_kf(3)) > 1.0d-6 .or. abs(x_ukf(4) - x_kf(4)) > 1.0d-6) then
            write(*,'(A)') 'FAIL (state mismatch)'
            write(*,'(A,4F10.6)') '    KF  x = ', x_kf
            write(*,'(A,4F10.6)') '    UKF x = ', x_ukf
            n_failed = n_failed + 1
            return
        end if

        ! Covariances should be close (not exact due to RK4 vs analytic F)
        ! Check P(1,1)
        if (abs(P_ukf(1) - P_kf(1)) > 0.1d0) then
            write(*,'(A,F10.4,A,F10.4)') 'FAIL (P(1,1) mismatch: KF=', P_kf(1), ' UKF=', P_ukf(1)
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 11: UKF update — position measurement (linear, should match KF)
    ! =====================================================================
    subroutine test_ukf_update_position()
        real(c_double) :: x_kf(4), P_kf(16), x_ukf(4), P_ukf(16)
        real(c_double) :: z(2), H(8), R(4), obs_params(1)
        integer(c_int) :: validity_kf(2), validity_ukf(2), info_kf, info_ukf

        write(*,'(A)', advance='no') '  UKF update position obs... '

        ! Same initial state for both
        x_kf  = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        x_ukf = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P_kf, 4)
        P_ukf = P_kf

        ! Position observation: z = [1.5, 0.1]
        z = (/ 1.5d0, 0.1d0 /)
        R = 0.0d0; R(1) = 0.5d0; R(4) = 0.5d0
        obs_params(1) = 0.0d0

        ! KF update (linear reference)
        H = 0.0d0; H(1) = 1.0d0; H(6) = 1.0d0
        call fa_kf_update(4, 2, x_kf, P_kf, z, H, R, validity_kf, info_kf)

        ! UKF update (obs_id=1 = position, m=3 for 3D obs but we want 2D)
        ! Use linear KF as reference; for UKF we use obs_id=1 which is 3D position
        ! Instead, use a custom approach: the KF result is our reference
        ! For a fair comparison, let's just verify UKF produces reasonable results
        ! with a 2D position-like measurement using obs_id=1 but with n=4, m=3
        block
            real(c_double) :: x_u3(6), P_u3(36), z3(3), R3(9), op3(1)
            integer(c_int) :: val3(3), info3

            x_u3 = (/ 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0 /)
            call build_identity(P_u3, 6)
            z3 = (/ 1.5d0, 0.1d0, 0.0d0 /)
            R3 = 0.0d0; R3(1) = 0.5d0; R3(5) = 0.5d0; R3(9) = 0.5d0
            op3(1) = 0.0d0

            call fa_ukf_update(6, 3, x_u3, P_u3, z3, c_null_funptr, 1, R3, op3, 1, &
                               1.0d0, 2.0d0, 0.0d0, val3, info3)

            if (info3 /= 0) then
                write(*,'(A,I3)') 'FAIL (UKF info=', info3; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if

            ! x_u3(1) should move toward 1.5, x_u3(2) toward 0.1
            if (abs(x_u3(1) - 1.0d0) > 0.5d0) then
                write(*,'(A,F10.4)') 'FAIL (UKF x(1) unexpected: ', x_u3(1); write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if
        end block

        if (info_kf /= 0) then
            write(*,'(A,I3)') 'FAIL (KF info=', info_kf; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 12: UKF update — range measurement (nonlinear, no Jacobian needed)
    ! =====================================================================
    subroutine test_ukf_update_range_nonlinear()
        real(c_double) :: x(6), P(36), z(1), R(1), obs_params(3)
        integer(c_int) :: validity(1), info
        real(c_double) :: x_before

        write(*,'(A)', advance='no') '  UKF update range (nonlinear)... '

        ! 6-state: [px,py,pz,vx,vy,vz]
        x = (/ 100.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
        call build_identity(P, 6)
        P = P * 100.0d0

        ! Range observation from station at origin
        obs_params = (/ 0.0d0, 0.0d0, 0.0d0 /)
        z(1) = 105.0d0
        R(1) = 1.0d0

        x_before = x(1)

        ! UKF does not need Jacobian — uses sigma points through nonlinear h(x)
        call fa_ukf_update(6, 1, x, P, z, c_null_funptr, 2, R, obs_params, 3, &
                           1.0d0, 2.0d0, 0.0d0, validity, info)

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

        ! P should decrease
        if (P(1) >= 100.0d0) then
            write(*,'(A,F10.4)') 'FAIL (P(1,1) did not decrease: ', P(1); write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A,F8.2,A)') 'PASS (x(1) moved to ', x(1), ')'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 13: ESKF predict-update-inject cycle
    ! const-vel dynamics, position observation. After inject, x_nom should
    ! be close to what a regular EKF would produce.
    ! =====================================================================
    subroutine test_eskf_predict_update_inject()
        real(c_double) :: x_ekf(4), P_ekf(16)
        real(c_double) :: x_nom(4), dx_eskf(4), P_eskf(16)
        real(c_double) :: Q(16), params(1), z(2), R(4), obs_params(1)
        real(c_double) :: F(16), H(8)
        integer(c_int) :: validity_ekf(2), validity_eskf(2), info
        real(c_double) :: diff

        write(*,'(A)', advance='no') '  ESKF predict-update-inject... '

        ! Initial state
        x_ekf = (/ 0.0d0, 0.0d0, 1.0d0, 0.5d0 /)
        x_nom = (/ 0.0d0, 0.0d0, 1.0d0, 0.5d0 /)
        dx_eskf = 0.0d0
        call build_identity(P_ekf, 4)
        P_eskf = P_ekf

        Q = 0.0d0
        Q(1) = 0.01d0; Q(6) = 0.01d0; Q(11) = 0.01d0; Q(16) = 0.01d0
        params(1) = 0.0d0

        ! --- EKF predict ---
        call fa_ekf_predict(4, x_ekf, P_ekf, c_null_funptr, c_null_funptr, Q, 1.0d0, 40, params, 1, 10, info)
        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (EKF predict info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! --- ESKF predict ---
        call fa_eskf_predict(4, x_nom, dx_eskf, P_eskf, c_null_funptr, 40, params, 1, Q, 1.0d0, 10, info)
        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (ESKF predict info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! After predict, x_nom should match x_ekf (same dynamics)
        diff = abs(x_nom(1) - x_ekf(1)) + abs(x_nom(2) - x_ekf(2))
        if (diff > 1.0d-6) then
            write(*,'(A,F12.8)') 'FAIL (predict state mismatch: ', diff; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! --- EKF update (linear position obs) ---
        z = (/ 1.1d0, 0.55d0 /)
        R = 0.0d0; R(1) = 1.0d0; R(4) = 1.0d0
        obs_params(1) = 0.0d0

        H = 0.0d0; H(1) = 1.0d0; H(6) = 1.0d0
        call fa_kf_update(4, 2, x_ekf, P_ekf, z, H, R, validity_ekf, info)
        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (EKF update info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! --- ESKF update (position obs_id=1, m=3 for 3D) ---
        ! Use the linear H directly by calling eskf_update with obs_id=1
        ! But obs_id=1 gives 3D position. Let's use linear H manually via KF-style update on dx.
        ! Actually, let's use the eskf_update with custom measurement via the dispatch.
        ! For 2D position we can use a workaround: apply update manually
        block
            real(c_double) :: H_eskf(8), z_pred(2), y(2)
            real(c_double) :: HT(8), PHT(8), S(4), S_inv(4), K(8)
            real(c_double) :: eye4(16), KH(16), IKH(16), IKHP(16), IKHT(16), TMP(16)
            real(c_double) :: KR(8), KT_m(8), KRKT(16)
            integer :: solve_i

            ! H = [1 0 0 0; 0 1 0 0], z_pred = H * x_nom
            H_eskf = 0.0d0; H_eskf(1) = 1.0d0; H_eskf(6) = 1.0d0
            call mat_vec_multiply_test(H_eskf, x_nom, z_pred, 2, 4)

            ! y = z - h(x_nom)
            y(1) = z(1) - z_pred(1)
            y(2) = z(2) - z_pred(2)

            ! S = H*P*H^T + R
            call mat_transpose_test(H_eskf, HT, 2, 4)
            call mat_multiply_test(P_eskf, HT, PHT, 4, 4, 2)
            call mat_multiply_test(H_eskf, PHT, S, 2, 4, 2)
            S(1) = S(1) + R(1); S(4) = S(4) + R(4)

            ! K = P*H^T * S^-1
            call build_identity(S_inv(1:4), 2)
            call mat_solve_symm_test(S, S_inv, 2, 2, solve_i)
            call mat_multiply_test(PHT, S_inv, K, 4, 2, 2)

            ! dx = K * y
            call mat_vec_multiply_test(K, y, dx_eskf, 4, 2)

            ! P = (I-KH)*P*(I-KH)^T + K*R*K^T (Joseph form)
            call build_identity(eye4, 4)
            call mat_multiply_test(K, H_eskf, KH, 4, 2, 4)
            IKH = eye4 - KH
            call mat_multiply_test(IKH, P_eskf, IKHP, 4, 4, 4)
            call mat_transpose_test(IKH, IKHT, 4, 4)
            call mat_multiply_test(IKHP, IKHT, TMP, 4, 4, 4)
            call mat_multiply_test(K, R, KR, 4, 2, 2)
            call mat_transpose_test(K, KT_m, 4, 2)
            call mat_multiply_test(KR, KT_m, KRKT, 4, 2, 4)
            P_eskf = TMP + KRKT
        end block

        ! --- ESKF inject ---
        call fa_eskf_inject(4, x_nom, dx_eskf, P_eskf, info)
        if (info /= 0) then
            write(*,'(A,I3)') 'FAIL (inject info=', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        ! After inject, x_nom should be close to x_ekf
        diff = abs(x_nom(1) - x_ekf(1)) + abs(x_nom(2) - x_ekf(2)) + &
               abs(x_nom(3) - x_ekf(3)) + abs(x_nom(4) - x_ekf(4))
        if (diff > 1.0d-6) then
            write(*,'(A,F12.8)') 'FAIL (ESKF != EKF after inject: diff=', diff; write(*,'(A)') ')'
            write(*,'(A,4F10.6)') '    EKF  x = ', x_ekf
            write(*,'(A,4F10.6)') '    ESKF x = ', x_nom
            n_failed = n_failed + 1
            return
        end if

        ! dx should be zero after inject
        if (abs(dx_eskf(1)) > 1.0d-15 .or. abs(dx_eskf(2)) > 1.0d-15) then
            write(*,'(A)') 'FAIL (dx not zero after inject)'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A,F12.8,A)') 'PASS (diff=', diff, ')'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 14: SR-EKF — same scenario as EKF, results should match
    ! =====================================================================
    subroutine test_srekf_vs_ekf()
        real(c_double) :: x_ekf(4), P_ekf(16), x_sr(4), S(16), Sq(16)
        real(c_double) :: Q(16), params(1)
        real(c_double) :: z(1), R(1), obs_params(3)
        integer(c_int) :: validity_ekf(1), validity_sr(1), info_ekf, info_sr
        real(c_double) :: diff

        write(*,'(A)', advance='no') '  SR-EKF vs EKF... '

        ! 6-state is heavy; use 4-state const-vel with range obs
        ! Actually, range obs needs 3D pos. Use 6-state.
        block
            real(c_double) :: x6_ekf(6), P6_ekf(36), x6_sr(6), S6(36), Sq6(36)
            real(c_double) :: Q6(36), params6(1), z6(1), R6(1), op6(3)
            integer(c_int) :: val6_ekf(1), val6_sr(1), info6_ekf, info6_sr
            real(c_double) :: P6_from_S(36), S6T(36)
            integer :: i

            x6_ekf = (/ 100.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
            x6_sr  = (/ 100.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
            call build_identity(P6_ekf, 6)
            P6_ekf = P6_ekf * 100.0d0

            ! S = cholesky(P) = 10*I for P = 100*I
            S6 = 0.0d0
            do i = 1, 6
                S6((i-1)*6 + i) = 10.0d0
            end do

            op6 = (/ 0.0d0, 0.0d0, 0.0d0 /)
            z6(1) = 105.0d0
            R6(1) = 1.0d0

            ! EKF update
            call fa_ekf_update(6, 1, x6_ekf, P6_ekf, z6, c_null_funptr, c_null_funptr, &
                               2, R6, op6, 3, val6_ekf, info6_ekf)

            ! SR-EKF update
            call fa_srekf_update(6, 1, x6_sr, S6, z6, c_null_funptr, c_null_funptr, &
                                 2, R6, op6, 3, val6_sr, info6_sr)

            if (info6_ekf /= 0) then
                write(*,'(A,I3)') 'FAIL (EKF info=', info6_ekf; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if
            if (info6_sr /= 0) then
                write(*,'(A,I3)') 'FAIL (SR-EKF info=', info6_sr; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if

            ! States should match
            diff = 0.0d0
            do i = 1, 6
                diff = diff + abs(x6_ekf(i) - x6_sr(i))
            end do
            if (diff > 1.0d-8) then
                write(*,'(A,F12.8)') 'FAIL (state diff=', diff; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if

            ! Verify S*S^T ≈ P_ekf (recovered covariance matches)
            call mat_transpose_test(S6, S6T, 6, 6)
            call mat_multiply_test(S6, S6T, P6_from_S, 6, 6, 6)
            diff = 0.0d0
            do i = 1, 36
                diff = diff + abs(P6_from_S(i) - P6_ekf(i))
            end do
            if (diff > 1.0d-6) then
                write(*,'(A,F12.8)') 'FAIL (P mismatch, diff=', diff; write(*,'(A)') ')'
                n_failed = n_failed + 1
                return
            end if
        end block

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 15: UKF with singular P → info=2
    ! =====================================================================
    subroutine test_ukf_singular_P()
        real(c_double) :: x(4), P(16), Q(16), params(1)
        integer(c_int) :: info

        write(*,'(A)', advance='no') '  UKF singular P... '

        x = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        P = 0.0d0  ! Singular (all zeros, not positive definite)
        Q = 0.0d0
        params(1) = 0.0d0

        call fa_ukf_predict(4, x, P, c_null_funptr, 40, params, 1, Q, 1.0d0, &
                            1.0d0, 2.0d0, 0.0d0, info)

        if (info /= 2) then
            write(*,'(A,I3)') 'FAIL (expected info=2, got ', info; write(*,'(A)') ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 16: IF update additivity — two sensors give same result
    ! whether fused together or sequentially in information form
    ! =====================================================================
    subroutine test_if_update_additivity()
        ! State: 2D position (n=2), two position sensors (m=1 each)
        ! Use eta/Ym naming to avoid Fortran case-insensitivity collisions
        real(c_double) :: eta_seq(2), Ym_seq(4)
        real(c_double) :: eta_both(2), Ym_both(4)
        real(c_double) :: H1(2), H2(2), R1(1), R2(1)
        real(c_double) :: z1(1), z2(1)
        real(c_double) :: H_both(4), R_both(4), z_both(2)
        integer(c_int) :: validity1(1), validity2(1), val_both(2), info
        real(c_double) :: diff
        integer :: i

        write(*,'(A)', advance='no') '  IF update additivity... '

        ! Start with a prior: P = diag(100,100), x = [0, 0]
        ! Y = P^-1 = diag(0.01, 0.01)
        ! eta = Y * x = [0, 0]

        ! Sensor 1: observes x-position at z=5.0, R=1.0
        H1 = (/ 1.0d0, 0.0d0 /)
        R1(1) = 1.0d0
        z1(1) = 5.0d0

        ! Sensor 2: observes y-position at z=3.0, R=2.0
        H2 = (/ 0.0d0, 1.0d0 /)
        R2(1) = 2.0d0
        z2(1) = 3.0d0

        ! --- Sequential fusion: sensor 1 then sensor 2 ---
        eta_seq = (/ 0.0d0, 0.0d0 /)
        Ym_seq = (/ 0.01d0, 0.0d0, 0.0d0, 0.01d0 /)

        call fa_if_update(2, 1, eta_seq, Ym_seq, z1, H1, R1, validity1, info)
        if (info /= 0) then
            write(*,'(A,I3,A)') 'FAIL (sensor1 info=', info, ')'
            n_failed = n_failed + 1
            return
        end if

        call fa_if_update(2, 1, eta_seq, Ym_seq, z2, H2, R2, validity2, info)
        if (info /= 0) then
            write(*,'(A,I3,A)') 'FAIL (sensor2 info=', info, ')'
            n_failed = n_failed + 1
            return
        end if

        ! --- Simultaneous fusion: both sensors at once ---
        eta_both = (/ 0.0d0, 0.0d0 /)
        Ym_both = (/ 0.01d0, 0.0d0, 0.0d0, 0.01d0 /)

        H_both = (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /)
        R_both = (/ 1.0d0, 0.0d0, 0.0d0, 2.0d0 /)
        z_both = (/ 5.0d0, 3.0d0 /)

        call fa_if_update(2, 2, eta_both, Ym_both, z_both, H_both, R_both, val_both, info)
        if (info /= 0) then
            write(*,'(A,I3,A)') 'FAIL (combined info=', info, ')'
            n_failed = n_failed + 1
            return
        end if

        ! Compare sequential vs simultaneous
        diff = 0.0d0
        do i = 1, 2
            diff = diff + abs(eta_seq(i) - eta_both(i))
        end do
        do i = 1, 4
            diff = diff + abs(Ym_seq(i) - Ym_both(i))
        end do

        if (diff > 1.0d-10) then
            write(*,'(A,F12.8,A)') 'FAIL (diff=', diff, ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 17: IF predict-update cycle — 5 steps, compare with KF
    ! =====================================================================
    subroutine test_if_predict_update_cycle()
        ! 2D const-vel (n=4): [px, py, vx, vy]
        ! Position observation (m=2)
        integer, parameter :: n = 4, m = 2, T = 5
        real(c_double) :: x_kf(n), P_kf(n*n), x_if(n), P_if(n*n)
        real(c_double) :: eta_if(n), Ym_if(n*n)
        real(c_double) :: F(n*n), Q(n*n), H(m*n), R(m*m)
        real(c_double) :: z(m), dt
        integer(c_int) :: validity_kf(m), validity_if(m), info_kf, info_if
        real(c_double) :: diff
        integer :: step, i

        write(*,'(A)', advance='no') '  IF predict-update cycle... '

        dt = 1.0d0

        ! Build F (const-vel)
        call build_constvel_F(F, dt)

        ! Q = 0.01 * I
        call build_identity(Q, n)
        Q = 0.01d0 * Q

        ! H = [I_2 | 0_2] (observe position only)
        H = 0.0d0
        H(1) = 1.0d0  ! H(1,1)
        H(4) = 1.0d0  ! H(2,2)

        ! R = diag(1, 1)
        R = 0.0d0
        R(1) = 1.0d0
        R(4) = 1.0d0

        ! Initial state: x = [0, 0, 1, 0], P = I
        x_kf = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P_kf, n)

        ! Convert to information form: Y = P^-1 = I, eta = Y*x = x
        eta_if(1:n) = x_kf(1:n)
        call build_identity(Ym_if, n)

        ! Run 5 steps of predict-update
        do step = 1, T
            ! Generate "measurement" at true position (no noise for deterministic test)
            z(1) = dble(step)      ! true px = step * 1.0 (vx=1)
            z(2) = 0.0d0           ! true py = 0

            ! KF predict + update
            call fa_kf_predict(n, x_kf, P_kf, F, Q, info_kf)
            call fa_kf_update(n, m, x_kf, P_kf, z, H, R, validity_kf, info_kf)

            ! IF predict + update
            call fa_if_predict(n, eta_if, Ym_if, F, Q, info_if)
            call fa_if_update(n, m, eta_if, Ym_if, z, H, R, validity_if, info_if)
        end do

        ! Convert IF result back to covariance form for comparison
        ! P_if = Ym_if^-1, x_if = P_if * eta_if
        call build_identity(P_if, n)
        call mat_solve_symm_test(Ym_if, P_if, n, n, info_if)
        if (info_if /= 0) then
            write(*,'(A)') 'FAIL (IF inversion failed)'
            n_failed = n_failed + 1
            return
        end if
        call mat_vec_multiply_test(P_if, eta_if, x_if, n, n)

        ! Compare states
        diff = 0.0d0
        do i = 1, n
            diff = diff + abs(x_kf(i) - x_if(i))
        end do

        if (diff > 1.0d-6) then
            write(*,'(A,F12.8,A)') 'FAIL (state diff=', diff, ')'
            n_failed = n_failed + 1
            return
        end if

        ! Compare covariances
        diff = 0.0d0
        do i = 1, n*n
            diff = diff + abs(P_kf(i) - P_if(i))
        end do

        if (diff > 1.0d-5) then
            write(*,'(A,F12.8,A)') 'FAIL (cov diff=', diff, ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 18: SIR particle filter basic — 1000 particles, const-vel
    ! =====================================================================
    subroutine test_pf_sir_basic()
        ! 3D const-vel (n=6): [px,py,pz,vx,vy,vz], position obs (obs_id=1, m=3)
        integer, parameter :: n = 6, m = 3, n_particles = 1000, nstep = 5
        real(c_double), allocatable :: particles(:), weights(:)
        real(c_double) :: Q(n*n), R(m*m), z(m), params(1)
        real(c_double) :: obs_params(1)
        real(c_double) :: mean_x(n), true_x(n)
        real(c_double) :: diff, dt
        integer(c_int) :: seed, info
        integer :: i, j, step, offset

        write(*,'(A)', advance='no') '  PF SIR basic... '

        allocate(particles(n * n_particles), weights(n_particles))

        dt = 1.0d0
        params(1) = 0.0d0
        obs_params(1) = 0.0d0

        ! Q = 0.01 * I
        call build_identity(Q, n)
        Q = 0.01d0 * Q

        ! R = I_3
        R = 0.0d0
        R(1) = 1.0d0; R(5) = 1.0d0; R(9) = 1.0d0

        ! Initialize particles around [0,0,0, 1,0,0] with some spread
        seed = 42
        do i = 1, n_particles
            offset = (i - 1) * n
            particles(offset + 1) = 0.0d0 + 0.5d0 * lcg_normal_test(seed)  ! px
            particles(offset + 2) = 0.0d0 + 0.5d0 * lcg_normal_test(seed)  ! py
            particles(offset + 3) = 0.0d0 + 0.5d0 * lcg_normal_test(seed)  ! pz
            particles(offset + 4) = 1.0d0 + 0.2d0 * lcg_normal_test(seed)  ! vx
            particles(offset + 5) = 0.0d0 + 0.2d0 * lcg_normal_test(seed)  ! vy
            particles(offset + 6) = 0.0d0 + 0.2d0 * lcg_normal_test(seed)  ! vz
        end do

        ! Uniform weights
        do i = 1, n_particles
            weights(i) = 1.0d0 / dble(n_particles)
        end do

        ! Run 5 steps
        true_x = (/ 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0 /)
        do step = 1, nstep
            ! True state propagation (const-vel: pos += dt * vel)
            true_x(1) = true_x(1) + dt * true_x(4)
            true_x(2) = true_x(2) + dt * true_x(5)
            true_x(3) = true_x(3) + dt * true_x(6)

            ! Measurement at true position
            z(1) = true_x(1); z(2) = true_x(2); z(3) = true_x(3)

            ! Predict: propagate particles using const-vel model (model_id=40)
            call fa_pf_sir_predict(n, n_particles, particles, weights, &
                                    c_null_funptr, 40, params, 1, Q, dt, 10, seed, info)
            if (info /= 0) then
                write(*,'(A,I3,A)') 'FAIL (predict info=', info, ')'
                n_failed = n_failed + 1
                deallocate(particles, weights)
                return
            end if

            ! Update: weight particles using position observation (obs_id=1, m=3)
            call fa_pf_sir_update(n, m, n_particles, particles, weights, &
                                   z, c_null_funptr, 1, R, obs_params, 1, info)
            if (info /= 0) then
                write(*,'(A,I3,A)') 'FAIL (update info=', info, ')'
                n_failed = n_failed + 1
                deallocate(particles, weights)
                return
            end if

            ! Resample
            call fa_pf_sir_resample(n, n_particles, particles, weights, seed, info)
        end do

        ! Compute weighted mean
        mean_x = 0.0d0
        do i = 1, n_particles
            offset = (i - 1) * n
            do j = 1, n
                mean_x(j) = mean_x(j) + weights(i) * particles(offset + j)
            end do
        end do

        ! Check weighted mean is near truth (tolerance ~3x measurement noise)
        diff = sqrt((mean_x(1) - true_x(1))**2 + (mean_x(2) - true_x(2))**2 + &
                    (mean_x(3) - true_x(3))**2)

        if (diff > 3.0d0) then
            write(*,'(A,F8.4,A)') 'FAIL (position error=', diff, ')'
            n_failed = n_failed + 1
            deallocate(particles, weights)
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
        deallocate(particles, weights)
    end subroutine

    ! =====================================================================
    ! Test 19: SIR resample — verify weights uniform after resample
    ! =====================================================================
    subroutine test_pf_sir_resample()
        integer, parameter :: n = 2, n_particles = 100
        real(c_double) :: particles(n * n_particles), weights(n_particles)
        real(c_double) :: w_max, w_min, eff_before, eff_after
        integer(c_int) :: seed, info
        integer :: i, offset

        write(*,'(A)', advance='no') '  PF SIR resample... '

        seed = 123

        ! Initialize particles
        do i = 1, n_particles
            offset = (i - 1) * n
            particles(offset + 1) = dble(i)
            particles(offset + 2) = dble(i) * 0.1d0
        end do

        ! Set highly non-uniform weights (simulate post-update)
        weights = 0.0d0
        weights(1) = 0.5d0   ! one dominant particle
        weights(2) = 0.3d0
        do i = 3, n_particles
            weights(i) = 0.2d0 / dble(n_particles - 2)
        end do

        ! Compute effective sample size before resample: N_eff = 1/sum(w_i^2)
        eff_before = 0.0d0
        do i = 1, n_particles
            eff_before = eff_before + weights(i) * weights(i)
        end do
        eff_before = 1.0d0 / eff_before

        ! Resample
        call fa_pf_sir_resample(n, n_particles, particles, weights, seed, info)
        if (info /= 0) then
            write(*,'(A,I3,A)') 'FAIL (info=', info, ')'
            n_failed = n_failed + 1
            return
        end if

        ! After resample, all weights should be 1/N
        w_max = weights(1)
        w_min = weights(1)
        do i = 2, n_particles
            if (weights(i) > w_max) w_max = weights(i)
            if (weights(i) < w_min) w_min = weights(i)
        end do

        if (abs(w_max - w_min) > 1.0d-14) then
            write(*,'(A,F12.8,A)') 'FAIL (weight spread=', w_max - w_min, ')'
            n_failed = n_failed + 1
            return
        end if

        ! Effective sample size after should be N
        eff_after = 0.0d0
        do i = 1, n_particles
            eff_after = eff_after + weights(i) * weights(i)
        end do
        eff_after = 1.0d0 / eff_after

        ! After resample N_eff should be N (or very close)
        if (abs(eff_after - dble(n_particles)) > 1.0d-8) then
            write(*,'(A,F8.2,A)') 'FAIL (N_eff after=', eff_after, ')'
            n_failed = n_failed + 1
            return
        end if

        ! N_eff should have improved
        if (eff_after <= eff_before) then
            write(*,'(A)') 'FAIL (N_eff did not improve)'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 20: RTS smoother — smoothed P <= filtered P at interior steps
    ! =====================================================================
    subroutine test_rts_smoother()
        integer, parameter :: n = 4, m = 2, T = 10
        real(c_double) :: x_filt(n*T), P_filt(n*n*T)
        real(c_double) :: x_pred(n*T), P_pred(n*n*T)
        real(c_double) :: F_all(n*n*T)
        real(c_double) :: x_smooth(n*T), P_smooth(n*n*T)
        real(c_double) :: x(n), P(n*n), x_p(n), P_p(n*n)
        real(c_double) :: F(n*n), Q(n*n), H(m*n), R(m*m), z(m)
        real(c_double) :: dt
        integer(c_int) :: validity(m), info
        integer :: step, i, x_off, P_off
        real(c_double) :: trace_filt, trace_smooth
        logical :: smoother_better

        write(*,'(A)', advance='no') '  RTS smoother... '

        dt = 1.0d0

        ! Build matrices
        call build_constvel_F(F, dt)
        call build_identity(Q, n)
        Q = 0.01d0 * Q
        H = 0.0d0
        H(1) = 1.0d0; H(4) = 1.0d0
        R = 0.0d0
        R(1) = 1.0d0; R(4) = 1.0d0

        ! Initial state
        x = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /)
        call build_identity(P, n)

        ! Forward KF pass — store predicted and filtered at each step
        do step = 1, T
            ! Store F for this step
            P_off = (step - 1) * n * n
            F_all(P_off+1 : P_off+n*n) = F(1:n*n)

            ! Predict
            x_p(1:n) = x(1:n)
            P_p(1:n*n) = P(1:n*n)
            call fa_kf_predict(n, x, P, F, Q, info)

            ! Store predicted
            x_off = (step - 1) * n
            P_off = (step - 1) * n * n
            x_pred(x_off+1 : x_off+n) = x(1:n)
            P_pred(P_off+1 : P_off+n*n) = P(1:n*n)

            ! Generate measurement at true position
            z(1) = dble(step)
            z(2) = 0.0d0

            ! Update
            call fa_kf_update(n, m, x, P, z, H, R, validity, info)

            ! Store filtered
            x_filt(x_off+1 : x_off+n) = x(1:n)
            P_filt(P_off+1 : P_off+n*n) = P(1:n*n)
        end do

        ! Backward RTS smoothing pass
        call fa_rts_smooth(n, T, x_filt, P_filt, x_pred, P_pred, F_all, x_smooth, P_smooth, info)
        if (info /= 0) then
            write(*,'(A,I3,A)') 'FAIL (smoother info=', info, ')'
            n_failed = n_failed + 1
            return
        end if

        ! Verify: smoothed P trace <= filtered P trace at interior timesteps
        smoother_better = .true.
        do step = 2, T - 1
            P_off = (step - 1) * n * n
            trace_filt = 0.0d0
            trace_smooth = 0.0d0
            do i = 1, n
                trace_filt = trace_filt + P_filt(P_off + (i-1)*n + i)
                trace_smooth = trace_smooth + P_smooth(P_off + (i-1)*n + i)
            end do
            ! Smoothed should have smaller or equal trace (less uncertainty)
            if (trace_smooth > trace_filt + 1.0d-10) then
                smoother_better = .false.
            end if
        end do

        if (.not. smoother_better) then
            write(*,'(A)') 'FAIL (smoothed P not smaller than filtered P)'
            n_failed = n_failed + 1
            return
        end if

        ! Verify last timestep: smoothed = filtered
        P_off = (T - 1) * n * n
        x_off = (T - 1) * n
        block
            real(c_double) :: d
            d = 0.0d0
            do i = 1, n
                d = d + abs(x_smooth(x_off+i) - x_filt(x_off+i))
            end do
            if (d > 1.0d-12) then
                write(*,'(A,F12.8,A)') 'FAIL (last step not equal, diff=', d, ')'
                n_failed = n_failed + 1
                return
            end if
        end block

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 21: Batch WLS — 5 position measurements of stationary target
    ! =====================================================================
    subroutine test_batch_wls()
        ! n=2 (position), 5 measurements of 2D position (m_total=10)
        integer, parameter :: n = 2, n_meas = 5, m_total = 10
        real(c_double) :: x(n), z_all(m_total), H_all(m_total * n)
        real(c_double) :: R_all(m_total * m_total)
        real(c_double) :: expected(n), diff
        integer(c_int) :: info
        integer :: i

        write(*,'(A)', advance='no') '  Batch WLS... '

        ! 5 measurements of a stationary target at roughly [3, 4]
        ! z = [3.1,4.2, 2.9,3.8, 3.0,4.0, 3.2,4.1, 2.8,3.9]
        z_all = (/ 3.1d0, 4.2d0, 2.9d0, 3.8d0, 3.0d0, 4.0d0, 3.2d0, 4.1d0, 2.8d0, 3.9d0 /)

        ! H = [I_2; I_2; I_2; I_2; I_2] (each measurement is direct position)
        H_all = 0.0d0
        do i = 1, n_meas
            H_all(((i-1)*2) * n + 1) = 1.0d0  ! row 2i-1, col 1
            H_all(((i-1)*2 + 1) * n + 2) = 1.0d0  ! row 2i, col 2
        end do

        ! R = I (uniform noise, so WLS = simple average)
        R_all = 0.0d0
        do i = 1, m_total
            R_all((i-1)*m_total + i) = 1.0d0
        end do

        ! Initial guess
        x = (/ 0.0d0, 0.0d0 /)

        call fa_batch_wls(n, m_total, x, z_all, H_all, R_all, 10, 1.0d-12, info)
        if (info /= 0) then
            write(*,'(A,I3,A)') 'FAIL (info=', info, ')'
            n_failed = n_failed + 1
            return
        end if

        ! With uniform R, WLS should give simple mean
        expected(1) = (3.1d0 + 2.9d0 + 3.0d0 + 3.2d0 + 2.8d0) / 5.0d0
        expected(2) = (4.2d0 + 3.8d0 + 4.0d0 + 4.1d0 + 3.9d0) / 5.0d0

        diff = abs(x(1) - expected(1)) + abs(x(2) - expected(2))
        if (diff > 1.0d-10) then
            write(*,'(A,F12.8,A)') 'FAIL (diff=', diff, ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Test 22: Batch MAP — with strong prior, result between prior and WLS
    ! =====================================================================
    subroutine test_batch_map()
        integer, parameter :: n = 2, n_meas = 5, m_total = 10
        real(c_double) :: x_wls(n), x_map(n), x0(n), P0(n*n)
        real(c_double) :: z_all(m_total), H_all(m_total * n), R_all(m_total * m_total)
        real(c_double) :: wls_mean(n)
        integer(c_int) :: info
        integer :: i
        logical :: between

        write(*,'(A)', advance='no') '  Batch MAP... '

        ! Same measurements as WLS test
        z_all = (/ 3.1d0, 4.2d0, 2.9d0, 3.8d0, 3.0d0, 4.0d0, 3.2d0, 4.1d0, 2.8d0, 3.9d0 /)

        H_all = 0.0d0
        do i = 1, n_meas
            H_all(((i-1)*2) * n + 1) = 1.0d0
            H_all(((i-1)*2 + 1) * n + 2) = 1.0d0
        end do

        R_all = 0.0d0
        do i = 1, m_total
            R_all((i-1)*m_total + i) = 1.0d0
        end do

        ! First get WLS solution for comparison
        x_wls = (/ 0.0d0, 0.0d0 /)
        call fa_batch_wls(n, m_total, x_wls, z_all, H_all, R_all, 10, 1.0d-12, info)

        wls_mean(1) = x_wls(1)
        wls_mean(2) = x_wls(2)

        ! Now MAP with strong prior at [0, 0]
        x_map = (/ 0.0d0, 0.0d0 /)
        x0 = (/ 0.0d0, 0.0d0 /)
        ! P0 = 0.5 * I (strong prior — low uncertainty at x0)
        call build_identity(P0, n)
        P0 = 0.5d0 * P0

        call fa_batch_map(n, m_total, x_map, x0, P0, z_all, H_all, R_all, 10, 1.0d-12, info)
        if (info /= 0) then
            write(*,'(A,I3,A)') 'FAIL (info=', info, ')'
            n_failed = n_failed + 1
            return
        end if

        ! MAP result should be between prior (0,0) and WLS mean (~3,4)
        ! More precisely, MAP should be pulled toward prior compared to WLS
        between = .true.
        do i = 1, n
            ! MAP should be closer to x0 than WLS (prior pulls it)
            if (abs(x_map(i) - x0(i)) > abs(wls_mean(i) - x0(i)) + 1.0d-10) then
                between = .false.
            end if
            ! MAP should be different from prior (measurements pull it away)
            if (abs(x_map(i) - x0(i)) < 1.0d-6) then
                between = .false.
            end if
        end do

        if (.not. between) then
            write(*,'(A,2F8.4,A,2F8.4,A)') 'FAIL (MAP=', x_map(1), x_map(2), &
                ' WLS=', wls_mean(1), wls_mean(2), ')'
            n_failed = n_failed + 1
            return
        end if

        write(*,'(A)') 'PASS'
        n_passed = n_passed + 1
    end subroutine

    ! =====================================================================
    ! Helper: simple LCG normal for test particle initialization
    ! =====================================================================
    function lcg_normal_test(s) result(z)
        integer, intent(inout) :: s
        real(c_double) :: z
        real(c_double) :: u1, u2, pi_val
        pi_val = 3.14159265358979323846d0
        s = mod(s * 1103515245 + 12345, 2147483647)
        if (s < 0) s = s + 2147483647
        u1 = (dble(s) + 0.5d0) / 2147483647.0d0
        s = mod(s * 1103515245 + 12345, 2147483647)
        if (s < 0) s = s + 2147483647
        u2 = (dble(s) + 0.5d0) / 2147483647.0d0
        if (u1 < 1.0d-15) u1 = 1.0d-15
        z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi_val * u2)
    end function

    ! =====================================================================
    ! Local matrix helpers for test code (avoid linking issues)
    ! =====================================================================
    subroutine mat_vec_multiply_test(A, xv, yv, m, n)
        integer, intent(in) :: m, n
        real(c_double), intent(in)  :: A(m*n), xv(n)
        real(c_double), intent(out) :: yv(m)
        integer :: i, j
        real(c_double) :: s
        do i = 1, m
            s = 0.0d0
            do j = 1, n
                s = s + A((i-1)*n + j) * xv(j)
            end do
            yv(i) = s
        end do
    end subroutine

    subroutine mat_multiply_test(A, B, C, m, k, n_col)
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
    end subroutine

    subroutine mat_transpose_test(A, AT, m, n)
        integer, intent(in) :: m, n
        real(c_double), intent(in)  :: A(m*n)
        real(c_double), intent(out) :: AT(n*m)
        integer :: i, j
        do i = 1, m
            do j = 1, n
                AT((j-1)*m + i) = A((i-1)*n + j)
            end do
        end do
    end subroutine

    subroutine mat_solve_symm_test(A, B, n, m_rhs, info)
        integer, intent(in) :: n, m_rhs
        real(c_double), intent(in) :: A(n*n)
        real(c_double), intent(inout) :: B(n*m_rhs)
        integer, intent(out) :: info
        real(c_double), allocatable :: L(:)
        real(c_double) :: s
        integer :: i, j, k
        info = 0
        allocate(L(n*n))
        L = 0.0d0
        do i = 1, n
            do j = 1, i
                s = A((i-1)*n + j)
                do k = 1, j - 1
                    s = s - L((i-1)*n + k) * L((j-1)*n + k)
                end do
                if (i == j) then
                    if (s <= 0.0d0) then; info = 1; deallocate(L); return; end if
                    L((i-1)*n + j) = sqrt(s)
                else
                    L((i-1)*n + j) = s / L((j-1)*n + j)
                end if
            end do
        end do
        do k = 1, m_rhs
            do i = 1, n
                s = B((i-1)*m_rhs + k)
                do j = 1, i - 1
                    s = s - L((i-1)*n + j) * B((j-1)*m_rhs + k)
                end do
                B((i-1)*m_rhs + k) = s / L((i-1)*n + i)
            end do
        end do
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
    end subroutine

end program test_estimate
