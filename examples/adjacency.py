import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wick_contraction import (
    GAMMA_0,
    GAMMA_5,
    ColorIndex,
    SpinIndex,
    QuarkField,
    Quark,
    QuarkBilinear,
    Diquark,
    SpinProjector,
    Correlator,
)


def show_simplified(name, correlator, order=None, degenerate=False):
    flags = []
    if order is not None:
        flags.append(f"order={order}")
    if degenerate:
        flags.append("degenerate=True")
    label = ", ".join(flags) if flags else "default"
    print(f"{name} simplified with {label}:")
    print("  Before:")
    print("   ", correlator)
    correlator.simplify(degenerate=degenerate, order=order)
    print("  After:")
    print("   ", correlator)
    print("  Einsum:")
    for term in correlator.terms:
        factor, einsum, operands = term.to_einsum()
        print(f"    {factor} {einsum}  |  {', '.join(operands)}")
    print()


if __name__ == "__main__":
    u = QuarkField("u")
    d = QuarkField("d")
    u_w = QuarkField.link("u", "z")
    dbar_gamma5_u = QuarkBilinear(d, u, GAMMA_5)
    dbar_w_gamma5_u = QuarkBilinear(d, u_w, GAMMA_5)
    ubar_gamma5_u = QuarkBilinear(u, u, GAMMA_5)
    dbar_gamma5_d = QuarkBilinear(d, d, GAMMA_5)
    dbar_gamma0_d = QuarkBilinear(d, d, GAMMA_0)
    u_gamma5_d = Diquark(u, d, GAMMA_5)
    u_gamma0 = Quark(u, GAMMA_0)

    # Pion two-point: ū γ5 d at both ends.
    pion = Correlator(dbar_gamma5_u.at("x") * dbar_gamma5_u.at("y").adjoint())
    show_simplified("Pion two-point function", pion, order=["x", "y"], degenerate=True)

    # Non-local pion two-point: Wilson-line-displaced u at the sink.
    non_local_pion = Correlator(dbar_w_gamma5_u.at("x") * dbar_gamma5_u.at("y").adjoint())
    show_simplified("Non-local pion two-point function", non_local_pion, order=["x", "z", "y"], degenerate=True)

    # Pion three-point: source at y, sink at x, intermediate d̄ γ0 d at z.
    pion_3pt = Correlator(dbar_gamma5_u.at("x") * dbar_gamma0_d.at("z") * dbar_gamma5_u.at("y").adjoint())
    show_simplified("Pion three-point function", pion_3pt, order=["z", "y", "x"], degenerate=True)

    # Projected proton two-point: χ_p = ϵ_abc (u_a^T C γ5 d_b) u_c, with parity projector P+.
    c_snk, s_snk = ColorIndex.new(), SpinIndex.new()
    c_src, s_src = ColorIndex.new(), SpinIndex.new()
    proton_snk = u_gamma5_d.at("x", c_snk) * u_gamma0.at("x", s_snk, c_snk)
    proton_src = (u_gamma5_d.at("y", c_src) * u_gamma0.at("y", s_src, c_src)).adjoint()
    p_plus = SpinProjector.P_plus(s_src, s_snk)
    proton = Correlator(p_plus * proton_snk * proton_src)
    show_simplified("Projected proton two-point function", proton, order=["x", "y"])

    # Flavor-singlet pseudoscalar η: (ūu + d̄d) / √2 at both ends.
    eta_snk = 2**-0.5 * (ubar_gamma5_u.at("x") + dbar_gamma5_d.at("x"))
    eta_src = 2**-0.5 * (ubar_gamma5_u.at("y") + dbar_gamma5_d.at("y")).adjoint()
    eta = Correlator(eta_snk * eta_src)
    show_simplified("Flavor-singlet pseudoscalar eta two-point function", eta, order=["x", "y"], degenerate=True)
