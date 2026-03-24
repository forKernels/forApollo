! Copyright The Fantastic Planet — By David Clabaugh
!
! forapollo_observe.f90 — Measurement model catalog with analytic Jacobians
!
! Built-in observation models for state estimation. Each model provides:
!   - h(x): predicted measurement z_pred = h(x, t, obs_params)
!   - H(x): analytic Jacobian dh/dx (m×n, flat row-major)
!
! Model IDs (position/range group 1-9):
!    1 = Position          — m=3, z = x(1:3)
!    2 = Range             — m=1, z = ||x(1:3) - station||, params=station(3)
!    3 = Bearing           — m=1, z = atan2(x(2)-s(2), x(1)-s(1)), params=station(3)
!    4 = Range+Bearing     — m=2, z = [range, bearing], params=station(3)
!
! Model IDs (velocity group 10-19):
!   10 = Velocity          — m=3, z = x(4:6)
!   11 = Doppler           — m=1, range-rate, params=station_pos(3)+station_vel(3)
!
! Model IDs (attitude sensors 20-29):
!   20 = Magnetometer      — m=3, z = R_body^T * ref_field, params=ref_field(3)
!   21 = Star tracker      — m=3, z = R_body^T * star_vec, params=star_vec(3)
!   22 = Sun sensor         — m=3, z = R_body^T * sun_vec, params=sun_vec(3)
!
! Model IDs (IMU sensors 30-39):
!   30 = Accelerometer     — m=3, z = R_body^T * (-gravity), params=gravity(3)
!   31 = Gyroscope         — m=3, z = x(11:13) angular rate from state
!
! Model IDs (radar 40-49):
!   40 = Radar             — m=3, z = [range, azimuth, elevation], params=station(3)
!
! Model IDs (scalar 50-59):
!   50 = Scalar            — m=1, z = x(int(obs_params(1)))
!
! Model IDs (camera 60-69):
!   60 = Pinhole           — m=2, z = [fx*x(1)/x(3)+cx, fy*x(2)/x(3)+cy], params=[fx,fy,cx,cy]
!
! Model IDs (relative 70-79):
!   70 = Relative position — m=3, z = x(1:3) - x2(1:3), params=x2(3)
!   71 = Relative velocity — m=3, z = x(4:6) - x2(1:3), params=x2(3)
!
! All bind(C) routines use value for scalar args, assumed-size for arrays.
! Internal model subroutines do NOT have bind(C) and do NOT use value.
! Jacobian H is flat 1D row-major: element (row, col) = H((row-1)*n + col).
!
! Error codes: info = 0 ok, 3 invalid input.

subroutine fa_observe_dispatch(obs_id, n, x, m, t, obs_params, nop, z_pred, info) &
    bind(C, name="fa_observe_dispatch")
    use iso_c_binding
    implicit none

    integer(c_int), value  :: obs_id, n, m, nop
    real(c_double), value  :: t
    real(c_double), intent(in)  :: x(n)
    real(c_double), intent(in)  :: obs_params(nop)
    real(c_double), intent(out) :: z_pred(m)
    integer(c_int), intent(out) :: info

    info = 0
    z_pred(1:m) = 0.0d0

    select case (obs_id)
    case (1)
        call observe_position(n, x, m, z_pred, info)
    case (2)
        call observe_range(n, x, m, obs_params, nop, z_pred, info)
    case (3)
        call observe_bearing(n, x, m, obs_params, nop, z_pred, info)
    case (4)
        call observe_range_bearing(n, x, m, obs_params, nop, z_pred, info)
    case (10)
        call observe_velocity(n, x, m, z_pred, info)
    case (11)
        call observe_doppler(n, x, m, obs_params, nop, z_pred, info)
    case (20)
        call observe_magnetometer(n, x, m, obs_params, nop, z_pred, info)
    case (21)
        call observe_star_tracker(n, x, m, obs_params, nop, z_pred, info)
    case (22)
        call observe_sun_sensor(n, x, m, obs_params, nop, z_pred, info)
    case (30)
        call observe_accelerometer(n, x, m, obs_params, nop, z_pred, info)
    case (31)
        call observe_gyroscope(n, x, m, z_pred, info)
    case (40)
        call observe_radar(n, x, m, obs_params, nop, z_pred, info)
    case (50)
        call observe_scalar(n, x, m, obs_params, nop, z_pred, info)
    case (60)
        call observe_pinhole(n, x, m, obs_params, nop, z_pred, info)
    case (70)
        call observe_relative_position(n, x, m, obs_params, nop, z_pred, info)
    case (71)
        call observe_relative_velocity(n, x, m, obs_params, nop, z_pred, info)
    case default
        info = 3
        return
    end select

contains

    ! =========================================================================
    ! Position (obs_id = 1)
    ! z = x(1:3), m = 3, n >= 3
    ! =========================================================================
    subroutine observe_position(n, x, m, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info

        if (n < 3 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        z_pred(1) = x(1)
        z_pred(2) = x(2)
        z_pred(3) = x(3)

    end subroutine observe_position

    ! =========================================================================
    ! Range (obs_id = 2)
    ! z = ||x(1:3) - station||, m = 1, params = station(3)
    ! =========================================================================
    subroutine observe_range(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy, dz, r

        if (n < 3 .or. m /= 1 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        dz = x(3) - obs_params(3)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        z_pred(1) = r

    end subroutine observe_range

    ! =========================================================================
    ! Bearing (obs_id = 3)
    ! z = atan2(x(2)-s(2), x(1)-s(1)), m = 1, params = station(3)
    ! =========================================================================
    subroutine observe_bearing(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy

        if (n < 3 .or. m /= 1 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        z_pred(1) = atan2(dy, dx)

    end subroutine observe_bearing

    ! =========================================================================
    ! Range + Bearing (obs_id = 4)
    ! z = [range, bearing], m = 2, params = station(3)
    ! =========================================================================
    subroutine observe_range_bearing(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy, dz, r

        if (n < 3 .or. m /= 2 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        dz = x(3) - obs_params(3)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        z_pred(1) = r
        z_pred(2) = atan2(dy, dx)

    end subroutine observe_range_bearing

    ! =========================================================================
    ! Velocity (obs_id = 10)
    ! z = x(4:6), m = 3, n >= 6
    ! =========================================================================
    subroutine observe_velocity(n, x, m, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info

        if (n < 6 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        z_pred(1) = x(4)
        z_pred(2) = x(5)
        z_pred(3) = x(6)

    end subroutine observe_velocity

    ! =========================================================================
    ! Doppler / Range-Rate (obs_id = 11)
    ! z = dot(dx, dv) / ||dx||
    ! where dx = x(1:3) - station_pos, dv = x(4:6) - station_vel
    ! params = [station_pos(3), station_vel(3)], nop >= 6
    ! =========================================================================
    subroutine observe_doppler(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx(3), dv(3), r, dot_dxdv

        if (n < 6 .or. m /= 1 .or. nop < 6) then
            info = 3
            return
        end if

        info = 0
        dx(1) = x(1) - obs_params(1)
        dx(2) = x(2) - obs_params(2)
        dx(3) = x(3) - obs_params(3)
        dv(1) = x(4) - obs_params(4)
        dv(2) = x(5) - obs_params(5)
        dv(3) = x(6) - obs_params(6)

        r = sqrt(dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))

        if (r < 1.0d-15) then
            ! Degenerate case: station at same position as target
            info = 3
            return
        end if

        dot_dxdv = dx(1)*dv(1) + dx(2)*dv(2) + dx(3)*dv(3)
        z_pred(1) = dot_dxdv / r

    end subroutine observe_doppler

    ! =========================================================================
    ! Helper: build rotation matrix from quaternion x(7:10) = [q0,q1,q2,q3]
    ! R rotates from inertial to body frame: v_body = R * v_inertial
    ! =========================================================================
    subroutine quat_to_dcm(x, n, R)
        implicit none
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out) :: R(3,3)
        real(c_double) :: q0, q1, q2, q3

        q0 = x(7)
        q1 = x(8)
        q2 = x(9)
        q3 = x(10)

        R(1,1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
        R(1,2) = 2.0d0*(q1*q2 + q0*q3)
        R(1,3) = 2.0d0*(q1*q3 - q0*q2)

        R(2,1) = 2.0d0*(q1*q2 - q0*q3)
        R(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3
        R(2,3) = 2.0d0*(q2*q3 + q0*q1)

        R(3,1) = 2.0d0*(q1*q3 + q0*q2)
        R(3,2) = 2.0d0*(q2*q3 - q0*q1)
        R(3,3) = q0*q0 - q1*q1 - q2*q2 + q3*q3

    end subroutine quat_to_dcm

    ! =========================================================================
    ! Magnetometer (obs_id = 20)
    ! z = R * ref_field_inertial, m = 3, n >= 10
    ! State quaternion at x(7:10), params = ref_field(3)
    ! =========================================================================
    subroutine observe_magnetometer(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: R(3,3), v(3)

        if (n < 10 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        v(1) = obs_params(1)
        v(2) = obs_params(2)
        v(3) = obs_params(3)

        call quat_to_dcm(x, n, R)

        ! z = R * v (matrix-vector product)
        z_pred(1) = R(1,1)*v(1) + R(1,2)*v(2) + R(1,3)*v(3)
        z_pred(2) = R(2,1)*v(1) + R(2,2)*v(2) + R(2,3)*v(3)
        z_pred(3) = R(3,1)*v(1) + R(3,2)*v(2) + R(3,3)*v(3)

    end subroutine observe_magnetometer

    ! =========================================================================
    ! Star Tracker (obs_id = 21)
    ! z = R * star_inertial, m = 3, n >= 10
    ! Same pattern as magnetometer with star vector
    ! =========================================================================
    subroutine observe_star_tracker(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: R(3,3), v(3)

        if (n < 10 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        v(1) = obs_params(1)
        v(2) = obs_params(2)
        v(3) = obs_params(3)

        call quat_to_dcm(x, n, R)

        z_pred(1) = R(1,1)*v(1) + R(1,2)*v(2) + R(1,3)*v(3)
        z_pred(2) = R(2,1)*v(1) + R(2,2)*v(2) + R(2,3)*v(3)
        z_pred(3) = R(3,1)*v(1) + R(3,2)*v(2) + R(3,3)*v(3)

    end subroutine observe_star_tracker

    ! =========================================================================
    ! Sun Sensor (obs_id = 22)
    ! z = R * sun_inertial, m = 3, n >= 10
    ! Same pattern as magnetometer with sun vector
    ! =========================================================================
    subroutine observe_sun_sensor(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: R(3,3), v(3)

        if (n < 10 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        v(1) = obs_params(1)
        v(2) = obs_params(2)
        v(3) = obs_params(3)

        call quat_to_dcm(x, n, R)

        z_pred(1) = R(1,1)*v(1) + R(1,2)*v(2) + R(1,3)*v(3)
        z_pred(2) = R(2,1)*v(1) + R(2,2)*v(2) + R(2,3)*v(3)
        z_pred(3) = R(3,1)*v(1) + R(3,2)*v(2) + R(3,3)*v(3)

    end subroutine observe_sun_sensor

    ! =========================================================================
    ! Accelerometer (obs_id = 30)
    ! Static case: z = R * (-gravity_inertial), m = 3, n >= 10
    ! params = gravity_inertial(3) (e.g., [0, 0, -9.81])
    ! =========================================================================
    subroutine observe_accelerometer(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: R(3,3), neg_g(3)

        if (n < 10 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        ! Negate gravity: accelerometer measures reaction to gravity
        neg_g(1) = -obs_params(1)
        neg_g(2) = -obs_params(2)
        neg_g(3) = -obs_params(3)

        call quat_to_dcm(x, n, R)

        z_pred(1) = R(1,1)*neg_g(1) + R(1,2)*neg_g(2) + R(1,3)*neg_g(3)
        z_pred(2) = R(2,1)*neg_g(1) + R(2,2)*neg_g(2) + R(2,3)*neg_g(3)
        z_pred(3) = R(3,1)*neg_g(1) + R(3,2)*neg_g(2) + R(3,3)*neg_g(3)

    end subroutine observe_accelerometer

    ! =========================================================================
    ! Gyroscope (obs_id = 31)
    ! z = x(11:13) angular rate from state, m = 3, n >= 13
    ! Assumes state layout: [r(3), v(3), q(4), omega(3)]
    ! =========================================================================
    subroutine observe_gyroscope(n, x, m, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info

        if (n < 13 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        z_pred(1) = x(11)
        z_pred(2) = x(12)
        z_pred(3) = x(13)

    end subroutine observe_gyroscope

    ! =========================================================================
    ! Radar (obs_id = 40)
    ! z = [range, azimuth, elevation], m = 3, params = station(3)
    ! azimuth = atan2(dy, dx), elevation = atan2(dz, sqrt(dx^2+dy^2))
    ! =========================================================================
    subroutine observe_radar(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy, dz, r, rxy

        if (n < 3 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        dz = x(3) - obs_params(3)

        r = sqrt(dx*dx + dy*dy + dz*dz)
        rxy = sqrt(dx*dx + dy*dy)

        z_pred(1) = r                      ! range
        z_pred(2) = atan2(dy, dx)           ! azimuth
        z_pred(3) = atan2(dz, rxy)          ! elevation

    end subroutine observe_radar

    ! =========================================================================
    ! Scalar (obs_id = 50)
    ! z = x(int(obs_params(1))), m = 1
    ! Observe a single state element by index
    ! =========================================================================
    subroutine observe_scalar(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        integer :: idx

        if (m /= 1 .or. nop < 1) then
            info = 3
            return
        end if

        idx = nint(obs_params(1))
        if (idx < 1 .or. idx > n) then
            info = 3
            return
        end if

        info = 0
        z_pred(1) = x(idx)

    end subroutine observe_scalar

    ! =========================================================================
    ! Pinhole Camera (obs_id = 60)
    ! z = [fx*x(1)/x(3) + cx, fy*x(2)/x(3) + cy], m = 2
    ! State is point in camera frame [X, Y, Z, ...]
    ! params = [fx, fy, cx, cy], nop >= 4
    ! =========================================================================
    subroutine observe_pinhole(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info
        real(c_double) :: fx, fy, cx, cy, inv_z

        if (n < 3 .or. m /= 2 .or. nop < 4) then
            info = 3
            return
        end if

        ! Avoid division by zero (point behind or on camera)
        if (abs(x(3)) < 1.0d-15) then
            info = 3
            return
        end if

        info = 0
        fx = obs_params(1)
        fy = obs_params(2)
        cx = obs_params(3)
        cy = obs_params(4)
        inv_z = 1.0d0 / x(3)

        z_pred(1) = fx * x(1) * inv_z + cx
        z_pred(2) = fy * x(2) * inv_z + cy

    end subroutine observe_pinhole

    ! =========================================================================
    ! Relative Position (obs_id = 70)
    ! z = x(1:3) - x2(1:3), m = 3, params = x2(3)
    ! =========================================================================
    subroutine observe_relative_position(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info

        if (n < 3 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        z_pred(1) = x(1) - obs_params(1)
        z_pred(2) = x(2) - obs_params(2)
        z_pred(3) = x(3) - obs_params(3)

    end subroutine observe_relative_position

    ! =========================================================================
    ! Relative Velocity (obs_id = 71)
    ! z = x(4:6) - x2(1:3), m = 3, params = x2(3)
    ! =========================================================================
    subroutine observe_relative_velocity(n, x, m, obs_params, nop, z_pred, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: z_pred(m)
        integer(c_int), intent(out) :: info

        if (n < 6 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        z_pred(1) = x(4) - obs_params(1)
        z_pred(2) = x(5) - obs_params(2)
        z_pred(3) = x(6) - obs_params(3)

    end subroutine observe_relative_velocity

end subroutine fa_observe_dispatch


! =============================================================================
! Jacobian dispatch: H = dh/dx, flat row-major m×n
! =============================================================================
subroutine fa_observe_jacobian(obs_id, n, x, m, t, obs_params, nop, H, info) &
    bind(C, name="fa_observe_jacobian")
    use iso_c_binding
    implicit none

    integer(c_int), value  :: obs_id, n, m, nop
    real(c_double), value  :: t
    real(c_double), intent(in)  :: x(n)
    real(c_double), intent(in)  :: obs_params(nop)
    real(c_double), intent(out) :: H(m*n)
    integer(c_int), intent(out) :: info

    info = 0
    H(1:m*n) = 0.0d0

    select case (obs_id)
    case (1)
        call jacobian_position(n, m, H, info)
    case (2)
        call jacobian_range(n, x, m, obs_params, nop, H, info)
    case (3)
        call jacobian_bearing(n, x, m, obs_params, nop, H, info)
    case (4)
        call jacobian_range_bearing(n, x, m, obs_params, nop, H, info)
    case (10)
        call jacobian_velocity(n, m, H, info)
    case (11)
        call jacobian_doppler(n, x, m, obs_params, nop, H, info)
    case (20)
        call jacobian_attitude_sensor(n, x, m, obs_params, nop, H, info)
    case (21)
        call jacobian_attitude_sensor(n, x, m, obs_params, nop, H, info)
    case (22)
        call jacobian_attitude_sensor(n, x, m, obs_params, nop, H, info)
    case (30)
        call jacobian_accelerometer(n, x, m, obs_params, nop, H, info)
    case (31)
        call jacobian_gyroscope(n, m, H, info)
    case (40)
        call jacobian_radar(n, x, m, obs_params, nop, H, info)
    case (50)
        call jacobian_scalar(n, m, obs_params, nop, H, info)
    case (60)
        call jacobian_pinhole(n, x, m, obs_params, nop, H, info)
    case (70)
        call jacobian_relative_position(n, m, H, info)
    case (71)
        call jacobian_relative_velocity(n, m, H, info)
    case default
        info = 3
        return
    end select

contains

    ! =========================================================================
    ! Jacobian: Position (obs_id = 1)
    ! H = [I_3 | 0], m=3, first 3 columns are identity
    ! Row-major: H((i-1)*n + j)
    ! =========================================================================
    subroutine jacobian_position(n, m, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info

        if (n < 3 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        ! H is already zeroed by caller
        ! Row 1: dz1/dx1 = 1
        H(1) = 1.0d0                       ! (0)*n + 1
        ! Row 2: dz2/dx2 = 1
        H(n + 2) = 1.0d0                   ! (1)*n + 2
        ! Row 3: dz3/dx3 = 1
        H(2*n + 3) = 1.0d0                 ! (2)*n + 3

    end subroutine jacobian_position

    ! =========================================================================
    ! Jacobian: Range (obs_id = 2)
    ! z = ||dx||, dz/dx_i = dx_i / ||dx||, m=1
    ! H is 1×n, only first 3 entries nonzero
    ! =========================================================================
    subroutine jacobian_range(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy, dz, r

        if (n < 3 .or. m /= 1 .or. nop < 3) then
            info = 3
            return
        end if

        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        dz = x(3) - obs_params(3)
        r = sqrt(dx*dx + dy*dy + dz*dz)

        if (r < 1.0d-15) then
            info = 3
            return
        end if

        info = 0
        H(1) = dx / r
        H(2) = dy / r
        H(3) = dz / r

    end subroutine jacobian_range

    ! =========================================================================
    ! Jacobian: Bearing (obs_id = 3)
    ! z = atan2(dy, dx)
    ! dz/dx1 = -dy / (dx^2+dy^2), dz/dx2 = dx / (dx^2+dy^2)
    ! H is 1×n, only first 2 entries nonzero
    ! =========================================================================
    subroutine jacobian_bearing(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy, r2

        if (n < 3 .or. m /= 1 .or. nop < 3) then
            info = 3
            return
        end if

        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        r2 = dx*dx + dy*dy

        if (r2 < 1.0d-30) then
            info = 3
            return
        end if

        info = 0
        H(1) = -dy / r2
        H(2) = dx / r2

    end subroutine jacobian_bearing

    ! =========================================================================
    ! Jacobian: Range+Bearing (obs_id = 4)
    ! Row 1: range Jacobian, Row 2: bearing Jacobian
    ! H is 2×n
    ! =========================================================================
    subroutine jacobian_range_bearing(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy, dz, r, r2

        if (n < 3 .or. m /= 2 .or. nop < 3) then
            info = 3
            return
        end if

        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        dz = x(3) - obs_params(3)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        r2 = dx*dx + dy*dy

        if (r < 1.0d-15 .or. r2 < 1.0d-30) then
            info = 3
            return
        end if

        info = 0
        ! Row 1: range partials
        H(1) = dx / r
        H(2) = dy / r
        H(3) = dz / r
        ! Row 2: bearing partials (offset by n)
        H(n + 1) = -dy / r2
        H(n + 2) = dx / r2
        ! H(n+3) = 0 (bearing doesn't depend on z)

    end subroutine jacobian_range_bearing

    ! =========================================================================
    ! Jacobian: Velocity (obs_id = 10)
    ! H = [0 | I_3 | 0], m=3, columns 4-6 are identity
    ! =========================================================================
    subroutine jacobian_velocity(n, m, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info

        if (n < 6 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        ! Row 1, col 4: dz1/dx4 = 1
        H(4) = 1.0d0                       ! (0)*n + 4
        ! Row 2, col 5: dz2/dx5 = 1
        H(n + 5) = 1.0d0                   ! (1)*n + 5
        ! Row 3, col 6: dz3/dx6 = 1
        H(2*n + 6) = 1.0d0                 ! (2)*n + 6

    end subroutine jacobian_velocity

    ! =========================================================================
    ! Jacobian: Doppler (obs_id = 11)
    ! z = dot(dx, dv) / r
    ! dz/dxi = dv_i/r - dot(dx,dv)*dx_i/r^3   (i=1,2,3 for position)
    ! dz/dxi = dx_(i-3)/r                       (i=4,5,6 for velocity)
    ! =========================================================================
    subroutine jacobian_doppler(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx(3), dv(3), r, r3, dot_dxdv
        integer :: i

        if (n < 6 .or. m /= 1 .or. nop < 6) then
            info = 3
            return
        end if

        dx(1) = x(1) - obs_params(1)
        dx(2) = x(2) - obs_params(2)
        dx(3) = x(3) - obs_params(3)
        dv(1) = x(4) - obs_params(4)
        dv(2) = x(5) - obs_params(5)
        dv(3) = x(6) - obs_params(6)

        r = sqrt(dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))
        if (r < 1.0d-15) then
            info = 3
            return
        end if

        info = 0
        r3 = r*r*r
        dot_dxdv = dx(1)*dv(1) + dx(2)*dv(2) + dx(3)*dv(3)

        ! Partials w.r.t. position states x(1:3)
        do i = 1, 3
            H(i) = dv(i)/r - dot_dxdv*dx(i)/r3
        end do

        ! Partials w.r.t. velocity states x(4:6)
        do i = 1, 3
            H(3 + i) = dx(i)/r
        end do

    end subroutine jacobian_doppler

    ! =========================================================================
    ! Jacobian: Attitude sensors (obs_id = 20, 21, 22)
    ! z = R(q) * v, where v is the reference vector (params)
    ! Jacobian w.r.t. quaternion states x(7:10)
    !
    ! For z_k = sum_j R_kj * v_j, we need dR_kj/dq_i for each quaternion
    ! component. The Jacobian H is 3×n with nonzero columns at 7,8,9,10.
    !
    ! dR/dq0, dR/dq1, dR/dq2, dR/dq3 computed analytically:
    !   dR/dq0 = 2*[q0, q3, -q2; -q3, q0, q1; q2, -q1, q0]
    !   dR/dq1 = 2*[q1, q2, q3; q2, -q1, q0; q3, -q0, -q1]
    !   dR/dq2 = 2*[-q2, q1, -q0; q1, q2, q3; q0, q3, -q2]
    !   dR/dq3 = 2*[-q3, q0, q1; -q0, -q3, q2; q1, q2, q3]
    ! =========================================================================
    subroutine jacobian_attitude_sensor(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: q0, q1, q2, q3, v(3)
        real(c_double) :: dR_dqi(3,3)
        integer :: row, qi

        if (n < 10 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        info = 0
        q0 = x(7)
        q1 = x(8)
        q2 = x(9)
        q3 = x(10)
        v(1) = obs_params(1)
        v(2) = obs_params(2)
        v(3) = obs_params(3)

        ! For each quaternion component qi (column 7,8,9,10 in H):
        ! H(row, 6+qi) = sum_j dR(row,j)/dq_{qi-1} * v(j)

        ! dR/dq0
        dR_dqi(1,1) = 2.0d0*q0;  dR_dqi(1,2) = 2.0d0*q3;  dR_dqi(1,3) = -2.0d0*q2
        dR_dqi(2,1) = -2.0d0*q3; dR_dqi(2,2) = 2.0d0*q0;  dR_dqi(2,3) = 2.0d0*q1
        dR_dqi(3,1) = 2.0d0*q2;  dR_dqi(3,2) = -2.0d0*q1; dR_dqi(3,3) = 2.0d0*q0
        do row = 1, 3
            H((row-1)*n + 7) = dR_dqi(row,1)*v(1) + dR_dqi(row,2)*v(2) + dR_dqi(row,3)*v(3)
        end do

        ! dR/dq1
        dR_dqi(1,1) = 2.0d0*q1;  dR_dqi(1,2) = 2.0d0*q2;  dR_dqi(1,3) = 2.0d0*q3
        dR_dqi(2,1) = 2.0d0*q2;  dR_dqi(2,2) = -2.0d0*q1; dR_dqi(2,3) = 2.0d0*q0
        dR_dqi(3,1) = 2.0d0*q3;  dR_dqi(3,2) = -2.0d0*q0; dR_dqi(3,3) = -2.0d0*q1
        do row = 1, 3
            H((row-1)*n + 8) = dR_dqi(row,1)*v(1) + dR_dqi(row,2)*v(2) + dR_dqi(row,3)*v(3)
        end do

        ! dR/dq2
        dR_dqi(1,1) = -2.0d0*q2; dR_dqi(1,2) = 2.0d0*q1;  dR_dqi(1,3) = -2.0d0*q0
        dR_dqi(2,1) = 2.0d0*q1;  dR_dqi(2,2) = 2.0d0*q2;  dR_dqi(2,3) = 2.0d0*q3
        dR_dqi(3,1) = 2.0d0*q0;  dR_dqi(3,2) = 2.0d0*q3;  dR_dqi(3,3) = -2.0d0*q2
        do row = 1, 3
            H((row-1)*n + 9) = dR_dqi(row,1)*v(1) + dR_dqi(row,2)*v(2) + dR_dqi(row,3)*v(3)
        end do

        ! dR/dq3
        dR_dqi(1,1) = -2.0d0*q3; dR_dqi(1,2) = 2.0d0*q0;  dR_dqi(1,3) = 2.0d0*q1
        dR_dqi(2,1) = -2.0d0*q0; dR_dqi(2,2) = -2.0d0*q3; dR_dqi(2,3) = 2.0d0*q2
        dR_dqi(3,1) = 2.0d0*q1;  dR_dqi(3,2) = 2.0d0*q2;  dR_dqi(3,3) = 2.0d0*q3
        do row = 1, 3
            H((row-1)*n + 10) = dR_dqi(row,1)*v(1) + dR_dqi(row,2)*v(2) + dR_dqi(row,3)*v(3)
        end do

    end subroutine jacobian_attitude_sensor

    ! =========================================================================
    ! Jacobian: Accelerometer (obs_id = 30)
    ! Same as attitude sensor but with negated gravity vector
    ! =========================================================================
    subroutine jacobian_accelerometer(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: neg_g(3)

        if (n < 10 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        ! Negate gravity and delegate to attitude sensor Jacobian
        neg_g(1) = -obs_params(1)
        neg_g(2) = -obs_params(2)
        neg_g(3) = -obs_params(3)

        call jacobian_attitude_sensor(n, x, m, neg_g, 3, H, info)

    end subroutine jacobian_accelerometer

    ! =========================================================================
    ! Jacobian: Gyroscope (obs_id = 31)
    ! H = [0...0 | I_3] at columns 11,12,13
    ! =========================================================================
    subroutine jacobian_gyroscope(n, m, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info

        if (n < 13 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        ! Row 1, col 11
        H(11) = 1.0d0                      ! (0)*n + 11
        ! Row 2, col 12
        H(n + 12) = 1.0d0                  ! (1)*n + 12
        ! Row 3, col 13
        H(2*n + 13) = 1.0d0                ! (2)*n + 13

    end subroutine jacobian_gyroscope

    ! =========================================================================
    ! Jacobian: Radar (obs_id = 40)
    ! z = [r, az, el] where r=||dx||, az=atan2(dy,dx), el=atan2(dz,rxy)
    !
    ! dr/dx1 = dx/r, dr/dx2 = dy/r, dr/dx3 = dz/r
    ! daz/dx1 = -dy/(dx^2+dy^2), daz/dx2 = dx/(dx^2+dy^2), daz/dx3 = 0
    ! del/dx1 = -dx*dz/(r^2*rxy), del/dx2 = -dy*dz/(r^2*rxy), del/dx3 = rxy/r^2
    ! =========================================================================
    subroutine jacobian_radar(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: dx, dy, dz, r, r2, rxy, rxy2

        if (n < 3 .or. m /= 3 .or. nop < 3) then
            info = 3
            return
        end if

        dx = x(1) - obs_params(1)
        dy = x(2) - obs_params(2)
        dz = x(3) - obs_params(3)

        r2 = dx*dx + dy*dy + dz*dz
        r = sqrt(r2)
        rxy2 = dx*dx + dy*dy
        rxy = sqrt(rxy2)

        if (r < 1.0d-15 .or. rxy < 1.0d-15) then
            info = 3
            return
        end if

        info = 0

        ! Row 1: range partials
        H(1) = dx / r
        H(2) = dy / r
        H(3) = dz / r

        ! Row 2: azimuth partials (offset by n)
        H(n + 1) = -dy / rxy2
        H(n + 2) = dx / rxy2
        ! H(n + 3) = 0 already

        ! Row 3: elevation partials (offset by 2*n)
        H(2*n + 1) = -dx*dz / (r2*rxy)
        H(2*n + 2) = -dy*dz / (r2*rxy)
        H(2*n + 3) = rxy / r2

    end subroutine jacobian_radar

    ! =========================================================================
    ! Jacobian: Scalar (obs_id = 50)
    ! H is 1×n with a single 1 at column = int(obs_params(1))
    ! =========================================================================
    subroutine jacobian_scalar(n, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        integer :: idx

        if (m /= 1 .or. nop < 1) then
            info = 3
            return
        end if

        idx = nint(obs_params(1))
        if (idx < 1 .or. idx > n) then
            info = 3
            return
        end if

        info = 0
        H(idx) = 1.0d0

    end subroutine jacobian_scalar

    ! =========================================================================
    ! Jacobian: Pinhole (obs_id = 60)
    ! z1 = fx*X/Z + cx, z2 = fy*Y/Z + cy
    ! dz1/dX = fx/Z, dz1/dY = 0, dz1/dZ = -fx*X/Z^2
    ! dz2/dX = 0, dz2/dY = fy/Z, dz2/dZ = -fy*Y/Z^2
    ! H is 2×n
    ! =========================================================================
    subroutine jacobian_pinhole(n, x, m, obs_params, nop, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m, nop
        real(c_double), intent(in)  :: x(n), obs_params(nop)
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info
        real(c_double) :: fx, fy, Z, inv_z, inv_z2

        if (n < 3 .or. m /= 2 .or. nop < 4) then
            info = 3
            return
        end if

        Z = x(3)
        if (abs(Z) < 1.0d-15) then
            info = 3
            return
        end if

        info = 0
        fx = obs_params(1)
        fy = obs_params(2)
        inv_z = 1.0d0 / Z
        inv_z2 = inv_z * inv_z

        ! Row 1: dz1/dX, dz1/dY, dz1/dZ
        H(1) = fx * inv_z                  ! (0)*n + 1
        ! H(2) = 0                         ! (0)*n + 2
        H(3) = -fx * x(1) * inv_z2         ! (0)*n + 3

        ! Row 2: dz2/dX, dz2/dY, dz2/dZ
        ! H(n + 1) = 0                     ! (1)*n + 1
        H(n + 2) = fy * inv_z              ! (1)*n + 2
        H(n + 3) = -fy * x(2) * inv_z2     ! (1)*n + 3

    end subroutine jacobian_pinhole

    ! =========================================================================
    ! Jacobian: Relative Position (obs_id = 70)
    ! H = [I_3 | 0], same as position Jacobian
    ! =========================================================================
    subroutine jacobian_relative_position(n, m, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info

        if (n < 3 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        H(1) = 1.0d0                       ! (0)*n + 1
        H(n + 2) = 1.0d0                   ! (1)*n + 2
        H(2*n + 3) = 1.0d0                 ! (2)*n + 3

    end subroutine jacobian_relative_position

    ! =========================================================================
    ! Jacobian: Relative Velocity (obs_id = 71)
    ! H = [0 | I_3 | 0], columns 4-6 are identity
    ! =========================================================================
    subroutine jacobian_relative_velocity(n, m, H, info)
        implicit none
        integer(c_int), intent(in)  :: n, m
        real(c_double), intent(out) :: H(m*n)
        integer(c_int), intent(out) :: info

        if (n < 6 .or. m /= 3) then
            info = 3
            return
        end if

        info = 0
        H(4) = 1.0d0                       ! (0)*n + 4
        H(n + 5) = 1.0d0                   ! (1)*n + 5
        H(2*n + 6) = 1.0d0                 ! (2)*n + 6

    end subroutine jacobian_relative_velocity

end subroutine fa_observe_jacobian
