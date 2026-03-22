! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_dynamics.f90 — Tracking dynamics models with analytic Jacobians
!
! Built-in dynamics catalog for state estimation. Each model provides:
!   - f(x): continuous-time state derivative x_dot = f(x, u, t)
!   - F(x): analytic Jacobian df/dx (flat row-major n×n)
!
! Model IDs (tracking group):
!   40 = Constant velocity   — state [pos(d), vel(d)], n = 2*d
!   41 = Constant acceleration — state [pos(d), vel(d), acc(d)], n = 3*d
!   42 = Constant turn rate  — state [x, y, vx, vy, omega], n = 5
!
! All bind(C) routines use value for scalar args, assumed-size for arrays.
! Internal model subroutines do NOT have bind(C) and do NOT use value.
! Matrices are flat 1D row-major: element (row, col) = F((row-1)*n + col).
!
! Error codes: info = 0 ok, 3 invalid input.

subroutine fa_dynamics_dispatch(model_id, n, x, u, nu, t, params, np, x_dot, info) &
    bind(C, name="fa_dynamics_dispatch")
    use iso_c_binding
    implicit none

    integer(c_int), value  :: model_id, n, nu, np
    real(c_double), value  :: t
    real(c_double), intent(in)  :: x(n)
    real(c_double), intent(in)  :: u(nu)
    real(c_double), intent(in)  :: params(np)
    real(c_double), intent(out) :: x_dot(n)
    integer(c_int), intent(out) :: info

    info = 0
    x_dot(1:n) = 0.0d0

    select case (model_id)
    case (40)
        call dynamics_const_vel(n, x, x_dot, info)
    case (41)
        call dynamics_const_accel(n, x, x_dot, info)
    case (42)
        call dynamics_const_turn(n, x, x_dot, info)
    case default
        info = 3
        return
    end select

contains

    ! =========================================================================
    ! Constant Velocity (model_id = 40)
    ! State: [pos_1, ..., pos_d, vel_1, ..., vel_d], n = 2*d
    ! x_dot = [vel_1, ..., vel_d, 0, ..., 0]
    ! =========================================================================
    subroutine dynamics_const_vel(n, x, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        integer :: d

        ! n must be even (n = 2*d)
        if (mod(n, 2) /= 0 .or. n < 2) then
            info = 3
            return
        end if

        info = 0
        d = n / 2

        ! Position derivatives = velocities
        x_dot(1:d) = x(d+1 : n)

        ! Velocity derivatives = 0 (constant velocity assumption)
        x_dot(d+1 : n) = 0.0d0

    end subroutine dynamics_const_vel

    ! =========================================================================
    ! Constant Acceleration (model_id = 41)
    ! State: [pos_1..d, vel_1..d, acc_1..d], n = 3*d
    ! x_dot = [vel_1..d, acc_1..d, 0..0]
    ! =========================================================================
    subroutine dynamics_const_accel(n, x, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        integer :: d

        ! n must be divisible by 3 (n = 3*d)
        if (mod(n, 3) /= 0 .or. n < 3) then
            info = 3
            return
        end if

        info = 0
        d = n / 3

        ! Position derivatives = velocities
        x_dot(1:d) = x(d+1 : 2*d)

        ! Velocity derivatives = accelerations
        x_dot(d+1 : 2*d) = x(2*d+1 : n)

        ! Acceleration derivatives = 0 (constant acceleration assumption)
        x_dot(2*d+1 : n) = 0.0d0

    end subroutine dynamics_const_accel

    ! =========================================================================
    ! Constant Turn Rate (model_id = 42)
    ! State: [x, y, vx, vy, omega], n = 5 exactly
    ! x_dot = [vx, vy, -omega*vy, omega*vx, 0]
    !
    ! This is the coordinated turn model in 2D where the object moves with
    ! constant speed and constant angular rate omega.
    ! =========================================================================
    subroutine dynamics_const_turn(n, x, x_dot, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: x_dot(n)
        integer(c_int), intent(out) :: info
        real(c_double) :: vx, vy, omega

        ! n must be exactly 5
        if (n /= 5) then
            info = 3
            return
        end if

        info = 0

        vx    = x(3)
        vy    = x(4)
        omega = x(5)

        x_dot(1) = vx                ! dx/dt = vx
        x_dot(2) = vy                ! dy/dt = vy
        x_dot(3) = -omega * vy       ! dvx/dt = -omega * vy
        x_dot(4) =  omega * vx       ! dvy/dt =  omega * vx
        x_dot(5) = 0.0d0             ! domega/dt = 0 (constant turn rate)

    end subroutine dynamics_const_turn

end subroutine fa_dynamics_dispatch


! =============================================================================
! Jacobian dispatch — analytic df/dx for each built-in dynamics model
! F is flat row-major n×n: element (i,j) stored at F((i-1)*n + j)
! =============================================================================
subroutine fa_dynamics_jacobian(model_id, n, x, u, nu, t, params, np, F, info) &
    bind(C, name="fa_dynamics_jacobian")
    use iso_c_binding
    implicit none

    integer(c_int), value  :: model_id, n, nu, np
    real(c_double), value  :: t
    real(c_double), intent(in)  :: x(n)
    real(c_double), intent(in)  :: u(nu)
    real(c_double), intent(in)  :: params(np)
    real(c_double), intent(out) :: F(n*n)
    integer(c_int), intent(out) :: info

    info = 0
    F(1:n*n) = 0.0d0

    select case (model_id)
    case (40)
        call jacobian_const_vel(n, x, F, info)
    case (41)
        call jacobian_const_accel(n, x, F, info)
    case (42)
        call jacobian_const_turn(n, x, F, info)
    case default
        info = 3
        return
    end select

contains

    ! =========================================================================
    ! Jacobian for Constant Velocity
    !
    !       [ 0  I ]
    !   F = [ 0  0 ]
    !
    ! where each block is d×d. The identity block maps velocity to position.
    ! Row-major: F((i-1)*n + j)
    ! =========================================================================
    subroutine jacobian_const_vel(n, x, F, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        integer :: d, i

        if (mod(n, 2) /= 0 .or. n < 2) then
            info = 3
            return
        end if

        info = 0
        d = n / 2

        F(1:n*n) = 0.0d0

        ! Upper-right identity block: dpos_i/dvel_i = 1
        ! Row i (1..d), Col d+i → index (i-1)*n + (d+i)
        do i = 1, d
            F((i-1)*n + d + i) = 1.0d0
        end do

    end subroutine jacobian_const_vel

    ! =========================================================================
    ! Jacobian for Constant Acceleration
    !
    !       [ 0  I  0 ]
    !   F = [ 0  0  I ]
    !       [ 0  0  0 ]
    !
    ! Each block is d×d. Identity blocks link pos→vel and vel→acc.
    ! =========================================================================
    subroutine jacobian_const_accel(n, x, F, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        integer :: d, i

        if (mod(n, 3) /= 0 .or. n < 3) then
            info = 3
            return
        end if

        info = 0
        d = n / 3

        F(1:n*n) = 0.0d0

        ! Block (1,2): dpos_i/dvel_i = 1
        ! Row i (1..d), Col d+i
        do i = 1, d
            F((i-1)*n + d + i) = 1.0d0
        end do

        ! Block (2,3): dvel_i/dacc_i = 1
        ! Row d+i (d+1..2d), Col 2d+i
        do i = 1, d
            F((d+i-1)*n + 2*d + i) = 1.0d0
        end do

    end subroutine jacobian_const_accel

    ! =========================================================================
    ! Jacobian for Constant Turn Rate
    !
    ! State: [x, y, vx, vy, omega]
    ! f = [vx, vy, -omega*vy, omega*vx, 0]
    !
    !       [ 0  0   1       0      0     ]
    !   F = [ 0  0   0       1      0     ]
    !       [ 0  0   0      -omega  -vy   ]
    !       [ 0  0   omega   0       vx   ]
    !       [ 0  0   0       0       0    ]
    !
    ! Partials:
    !   df3/dvx = 0, df3/dvy = -omega, df3/domega = -vy
    !   df4/dvx = omega, df4/dvy = 0, df4/domega = vx
    ! =========================================================================
    subroutine jacobian_const_turn(n, x, F, info)
        implicit none
        integer(c_int), intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: F(n*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: vx, vy, omega

        if (n /= 5) then
            info = 3
            return
        end if

        info = 0

        vx    = x(3)
        vy    = x(4)
        omega = x(5)

        F(1:n*n) = 0.0d0

        ! Row 1: dx/dt = vx → df1/dvx = 1 → F(1,3) = (0)*5+3 = 3
        F(3) = 1.0d0

        ! Row 2: dy/dt = vy → df2/dvy = 1 → F(2,4) = (1)*5+4 = 9
        F(9) = 1.0d0

        ! Row 3: dvx/dt = -omega*vy
        !   df3/dvy    = -omega → F(3,4) = (2)*5+4 = 14
        !   df3/domega = -vy    → F(3,5) = (2)*5+5 = 15
        F(14) = -omega
        F(15) = -vy

        ! Row 4: dvy/dt = omega*vx
        !   df4/dvx    =  omega → F(4,3) = (3)*5+3 = 18
        !   df4/domega =  vx    → F(4,5) = (3)*5+5 = 20
        F(18) = omega
        F(20) = vx

        ! Row 5: domega/dt = 0 → all zeros (already set)

    end subroutine jacobian_const_turn

end subroutine fa_dynamics_jacobian
