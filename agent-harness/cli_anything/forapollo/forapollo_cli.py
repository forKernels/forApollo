"""Click CLI + REPL for forApollo -- guidance, navigation, control.

Groups: nav, guidance, orbital
"""

import json
import math
import click
import numpy as np

from .utils.forapollo_backend import (
    kf_predict, kf_update, ekf_predict, ukf_predict,
    eskf_predict, eskf_inject, if_predict,
    guidance_zem, guidance_zev, guidance_lambert, guidance_lqr,
    guidance_pure_pursuit, guidance_stanley, guidance_pn_pure, guidance_min_energy,
    coords_eci_to_ecef, coords_ecef_to_eci,
    coords_ecef_to_geodetic, coords_geodetic_to_ecef,
)


def _output(result, as_json):
    if as_json:
        if isinstance(result, np.ndarray):
            click.echo(json.dumps({"array": result.tolist(), "shape": list(result.shape)}))
        elif isinstance(result, dict):
            out = {}
            for k, v in result.items():
                out[k] = v.tolist() if isinstance(v, np.ndarray) else v
            click.echo(json.dumps(out, indent=2))
        elif isinstance(result, tuple):
            click.echo(json.dumps([
                r.tolist() if isinstance(r, np.ndarray) else r for r in result
            ]))
        else:
            click.echo(json.dumps(result))
    else:
        click.echo(result)


def _parse_vec(s):
    return np.array([float(v) for v in s.split(",")], dtype=np.float64)


# =========================================================================
# Main CLI group
# =========================================================================

@click.group(invoke_without_command=True)
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
@click.pass_context
def cli(ctx, as_json):
    """forApollo -- Guidance, Navigation & Control CLI."""
    ctx.ensure_object(dict)
    ctx.obj["json"] = as_json
    if ctx.invoked_subcommand is None:
        _repl(as_json)


# =========================================================================
# nav group -- Kalman filters
# =========================================================================

@cli.group()
@click.pass_context
def nav(ctx):
    """Navigation estimators: KF, EKF, UKF, ESKF, particle filter."""
    pass


@nav.command("kf-predict")
@click.option("-n", type=int, required=True, help="State dimension.")
@click.argument("state")
@click.argument("covariance")
@click.argument("transition")
@click.argument("process_noise")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def nav_kf_pred(ctx, n, state, covariance, transition, process_noise, as_json):
    """Kalman filter predict. Flat vectors: state(n), P(n*n), F(n*n), Q(n*n)."""
    as_json = as_json or ctx.obj.get("json", False)
    x, P, info = kf_predict(n, _parse_vec(state), _parse_vec(covariance),
                             _parse_vec(transition), _parse_vec(process_noise))
    _output({"state": x, "covariance_diag": np.diag(P.reshape(n, n)), "info": info}, as_json)


@nav.command("ekf-predict")
@click.option("-n", type=int, required=True, help="State dimension.")
@click.argument("state")
@click.argument("covariance")
@click.argument("process_noise")
@click.option("--dt", type=float, required=True, help="Time step.")
@click.option("--model", type=int, default=0, help="Built-in dynamics model ID.")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def nav_ekf_pred(ctx, n, state, covariance, process_noise, dt, model, as_json):
    """EKF predict with built-in dynamics model."""
    as_json = as_json or ctx.obj.get("json", False)
    x, P, info = ekf_predict(n, _parse_vec(state), _parse_vec(covariance),
                              _parse_vec(process_noise), dt, model)
    _output({"state": x, "covariance_diag": np.diag(P.reshape(n, n)), "info": info}, as_json)


@nav.command("ukf-predict")
@click.option("-n", type=int, required=True, help="State dimension.")
@click.argument("state")
@click.argument("covariance")
@click.argument("process_noise")
@click.option("--dt", type=float, required=True, help="Time step.")
@click.option("--model", type=int, default=0)
@click.option("--alpha", type=float, default=1e-3)
@click.option("--beta", type=float, default=2.0)
@click.option("--kappa", type=float, default=0.0)
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def nav_ukf_pred(ctx, n, state, covariance, process_noise, dt, model, alpha, beta, kappa, as_json):
    """UKF predict with sigma points."""
    as_json = as_json or ctx.obj.get("json", False)
    x, P, info = ukf_predict(n, _parse_vec(state), _parse_vec(covariance),
                              _parse_vec(process_noise), dt, model, alpha, beta, kappa)
    _output({"state": x, "covariance_diag": np.diag(P.reshape(n, n)), "info": info}, as_json)


@nav.command("eskf-predict")
@click.option("-n", type=int, required=True, help="State dimension.")
@click.argument("nominal")
@click.argument("error_state")
@click.argument("covariance")
@click.argument("process_noise")
@click.option("--dt", type=float, required=True)
@click.option("--model", type=int, default=0)
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def nav_eskf_pred(ctx, n, nominal, error_state, covariance, process_noise, dt, model, as_json):
    """Error-State KF predict (Apollo heritage ESKF)."""
    as_json = as_json or ctx.obj.get("json", False)
    x_nom, dx, P, info = eskf_predict(n, _parse_vec(nominal), _parse_vec(error_state),
                                       _parse_vec(covariance), _parse_vec(process_noise), dt, model)
    _output({"nominal": x_nom, "error_state": dx,
             "covariance_diag": np.diag(P.reshape(n, n)), "info": info}, as_json)


@nav.command("eskf-inject")
@click.option("-n", type=int, required=True)
@click.argument("nominal")
@click.argument("error_state")
@click.argument("covariance")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def nav_eskf_inj(ctx, n, nominal, error_state, covariance, as_json):
    """ESKF error-state injection (reset after update)."""
    as_json = as_json or ctx.obj.get("json", False)
    x_nom, P, info = eskf_inject(n, _parse_vec(nominal), _parse_vec(error_state),
                                  _parse_vec(covariance))
    _output({"nominal": x_nom, "covariance_diag": np.diag(P.reshape(n, n)), "info": info}, as_json)


# =========================================================================
# guidance group
# =========================================================================

@cli.group()
@click.pass_context
def guidance(ctx):
    """Guidance laws: ZEM/ZEV, Lambert, LQR, pursuit, proportional nav."""
    pass


@guidance.command("zem")
@click.option("-n", type=int, required=True, help="Spatial dimension (3 for 3D).")
@click.argument("state")
@click.option("--tgo", type=float, required=True, help="Time-to-go (s).")
@click.argument("gravity")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def guid_zem(ctx, n, state, tgo, gravity, as_json):
    """Zero-Effort Miss guidance (Apollo P63). State: pos+vel (2*n), gravity: (n)."""
    as_json = as_json or ctx.obj.get("json", False)
    a_cmd, info = guidance_zem(n, _parse_vec(state), tgo, _parse_vec(gravity))
    _output({"acceleration_cmd": a_cmd, "info": info, "law": "ZEM"}, as_json)


@guidance.command("zev")
@click.option("-n", type=int, required=True)
@click.argument("state")
@click.argument("v_target")
@click.option("--tgo", type=float, required=True)
@click.argument("gravity")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def guid_zev(ctx, n, state, v_target, tgo, gravity, as_json):
    """Zero-Effort Velocity guidance. Includes target velocity."""
    as_json = as_json or ctx.obj.get("json", False)
    a_cmd, info = guidance_zev(n, _parse_vec(state), _parse_vec(v_target),
                                tgo, _parse_vec(gravity))
    _output({"acceleration_cmd": a_cmd, "info": info, "law": "ZEV"}, as_json)


@guidance.command("lambert")
@click.argument("r1")
@click.argument("r2")
@click.option("--tof", type=float, required=True, help="Time of flight (s).")
@click.option("--mu", type=float, default=3.986004418e14, help="Gravitational parameter.")
@click.option("--direction", type=int, default=0, help="0=short, 1=long way.")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def guid_lambert(ctx, r1, r2, tof, mu, direction, as_json):
    """Lambert problem: compute transfer orbit. r1/r2 as x,y,z (m)."""
    as_json = as_json or ctx.obj.get("json", False)
    v1, v2, info = guidance_lambert(_parse_vec(r1), _parse_vec(r2), tof, mu, direction)
    dv1 = float(np.linalg.norm(v1))
    dv2 = float(np.linalg.norm(v2))
    _output({"v1": v1, "v2": v2, "delta_v1_mag": dv1, "delta_v2_mag": dv2,
             "total_delta_v": dv1 + dv2, "info": info, "law": "Lambert"}, as_json)


@guidance.command("lqr")
@click.option("-n", type=int, required=True, help="State dimension.")
@click.option("-m", type=int, required=True, help="Control dimension.")
@click.argument("A")
@click.argument("B")
@click.argument("Q")
@click.argument("R")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def guid_lqr(ctx, n, m, a, b, q, r, as_json):
    """LQR controller: compute gain K. Matrices as flat comma-separated values."""
    as_json = as_json or ctx.obj.get("json", False)
    K, S, info = guidance_lqr(n, m, _parse_vec(a), _parse_vec(b),
                               _parse_vec(q), _parse_vec(r))
    _output({"K": K, "S_diag": np.diag(S), "info": info, "law": "LQR"}, as_json)


@guidance.command("pure-pursuit")
@click.argument("position")
@click.argument("target")
@click.option("--lookahead", type=float, required=True)
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def guid_pp(ctx, position, target, lookahead, as_json):
    """Pure pursuit path tracking. Position and target as x,y."""
    as_json = as_json or ctx.obj.get("json", False)
    steer, info = guidance_pure_pursuit(_parse_vec(position), _parse_vec(target), lookahead)
    _output({"steering_angle": steer, "info": info, "law": "pure_pursuit"}, as_json)


@guidance.command("pn-pure")
@click.option("-n", type=int, required=True, help="Spatial dimension.")
@click.argument("pursuer")
@click.argument("target")
@click.option("--gain", type=float, default=3.0, help="Navigation gain N.")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def guid_pn(ctx, n, pursuer, target, gain, as_json):
    """Pure proportional navigation guidance."""
    as_json = as_json or ctx.obj.get("json", False)
    a_cmd, info = guidance_pn_pure(n, _parse_vec(pursuer), _parse_vec(target), gain)
    _output({"acceleration_cmd": a_cmd, "info": info, "law": "PN_pure"}, as_json)


# =========================================================================
# orbital group
# =========================================================================

@cli.group()
@click.pass_context
def orbital(ctx):
    """Orbital mechanics: coordinate transforms, geodetic conversions."""
    pass


@orbital.command("eci-to-ecef")
@click.argument("position")
@click.option("--gmst", type=float, required=True, help="Greenwich Mean Sidereal Time (rad).")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def orb_eci2ecef(ctx, position, gmst, as_json):
    """Convert ECI (J2000) to ECEF. Position as x,y,z (m)."""
    as_json = as_json or ctx.obj.get("json", False)
    result, info = coords_eci_to_ecef(_parse_vec(position), gmst)
    _output({"ecef": result, "info": info}, as_json)


@orbital.command("ecef-to-eci")
@click.argument("position")
@click.option("--gmst", type=float, required=True, help="GMST (rad).")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def orb_ecef2eci(ctx, position, gmst, as_json):
    """Convert ECEF to ECI (J2000). Position as x,y,z (m)."""
    as_json = as_json or ctx.obj.get("json", False)
    result, info = coords_ecef_to_eci(_parse_vec(position), gmst)
    _output({"eci": result, "info": info}, as_json)


@orbital.command("ecef-to-geodetic")
@click.argument("position")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def orb_ecef2geo(ctx, position, as_json):
    """Convert ECEF to geodetic (lat, lon, alt) using WGS84. Position as x,y,z (m)."""
    as_json = as_json or ctx.obj.get("json", False)
    result, info = coords_ecef_to_geodetic(_parse_vec(position))
    result["lat_deg"] = math.degrees(result["lat"])
    result["lon_deg"] = math.degrees(result["lon"])
    result["info"] = info
    _output(result, as_json)


@orbital.command("geodetic-to-ecef")
@click.option("--lat", type=float, required=True, help="Latitude (rad).")
@click.option("--lon", type=float, required=True, help="Longitude (rad).")
@click.option("--alt", type=float, required=True, help="Altitude (m).")
@click.option("--json", "as_json", is_flag=True)
@click.pass_context
def orb_geo2ecef(ctx, lat, lon, alt, as_json):
    """Convert geodetic (lat, lon, alt) to ECEF. WGS84 ellipsoid."""
    as_json = as_json or ctx.obj.get("json", False)
    result, info = coords_geodetic_to_ecef(lat, lon, alt)
    _output({"ecef": result, "info": info}, as_json)


# =========================================================================
# REPL
# =========================================================================

def _repl(as_json):
    """Interactive REPL for forApollo."""
    click.echo("forApollo REPL -- Guidance, Navigation & Control")
    click.echo("The same Kalman filter that guided Apollo to the Moon.")
    click.echo("Type 'help' for commands, 'quit' to exit.\n")

    while True:
        try:
            line = input("forapollo> ").strip()
        except (EOFError, KeyboardInterrupt):
            click.echo("\nBye.")
            break

        if not line:
            continue
        if line in ("quit", "exit", "q"):
            break
        if line == "help":
            click.echo("Groups: nav, guidance, orbital")
            click.echo("Run: cli-anything-forapollo <group> --help")
            click.echo("Examples:")
            click.echo("  forapollo guidance zem -n 3 0,0,1000,0,0,-50 --tgo 20 0,0,-9.81 --json")
            click.echo("  forapollo orbital ecef-to-geodetic 6378137,0,0 --json")
            click.echo("  forapollo guidance lambert 7000000,0,0 0,8000000,0 --tof 3600 --json")
            continue

        parts = line.split()
        try:
            cli(parts, standalone_mode=False, obj={"json": as_json})
        except click.exceptions.UsageError as e:
            click.echo(f"Error: {e}")
        except SystemExit:
            pass
        except Exception as e:
            click.echo(f"Error: {e}")
