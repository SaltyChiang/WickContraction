# WickContraction

Perform Wick contraction in lattice QCD. Input a correlation function in the
form of interpolating field operators, and the algorithm returns the contraction
in the form of adjacency matrices.

## Interpolating Field Operators

First, we define 5 blocks.

- Quark (column vector)

$$Q_{\alpha,a}(x,q,\Gamma)=\delta_{ab}\Gamma_{\alpha\beta}q_{\beta,b}(x)$$

- Anti-quark (row vector)

$$\bar{Q}_{\alpha,a}(x,q,\Gamma)=\delta_{ab}\bar{q}_{\beta,b}(x)\Gamma_{\beta\alpha}$$

- Quark bilinear

$$B(x,q,f,\Gamma)=\delta_{ab}\bar{q}_{\alpha,a}(x)\Gamma_{\alpha\beta}f_{\beta,b}(x)$$

- Diquark (column vector)

$$D_{c}(x,q,f,\Gamma)=\epsilon_{abc}q^T_{\alpha,a}(x)(C\Gamma)_{\alpha\beta}f_{\beta,b}$$

- Anti-diquark (row vector)

$$\bar{D}_{c}(x,q,f,\Gamma)=\epsilon_{abc}\bar{q}_{\alpha,a}(x)({\Gamma}C)_{\alpha\beta}\bar{f}^T_{\beta,b}$$

Here we have the convention that $q$ and $f$ are quark fields with different
flavors (we usually use $u,d,s,c,t,b$ as flavors), Latin letters $a,b,c,\dots$
represents for color indices, Greek letters $\alpha,\beta,\gamma,\dots$
represents spin indices. Here, the transpose operation $T$ only indicates the
shape (row or column vector) of the field, and does not change the order of the
indices.

And we need to define the corresponding Dirac conjugate (adjoint) operation:

- Quark

$$Q_{\beta,a}^\dagger(x,q,\Gamma)(\gamma_4)_{\beta\alpha}=\delta_{ab}\bar{q}_{\beta,b}(x)(\gamma_4\Gamma^\dagger\gamma_4)_{\beta\alpha}$$
$$Q_{\beta,a}^\dagger(x,q,\Gamma)(\gamma_4)_{\beta\alpha}=\bar{Q}_{\alpha,a}(x,q,\gamma_4\Gamma^\dagger\gamma_4)$$

- Anti-quark

$$(\gamma_4)_{\alpha\beta}\bar{Q}_{\beta,a}^\dagger(x,q,\Gamma)=\delta_{ab}(\gamma_4\Gamma^\dagger\gamma_4)_{\alpha\beta}q_{\beta,b}(x)$$
$$(\gamma_4)_{\alpha\beta}\bar{Q}_{\beta,a}^\dagger(x,q,\Gamma)=Q_{\alpha,a}(x,q,\gamma_4\Gamma^\dagger\gamma_4)$$

- Quark bilinear

$$B^\dagger(x,q,f,\Gamma)=\delta_{ab}\bar{f}_{\alpha,a}(x)(\gamma_4\Gamma^\dagger\gamma_4)_{\alpha\beta}q_{\beta,b}(x)$$
$$B^\dagger(x,q,f,\Gamma)=B(x,f,q,\gamma_4\Gamma^\dagger\gamma_4)$$

- Diquark

$$D_{c}^\dagger(x,q,f,\Gamma)=\epsilon_{abc}\bar{f}_{\alpha,a}(x)(\gamma_4\Gamma^\dagger\gamma_4C)_{\alpha\beta}\bar{q}^T_{\beta,b}(x)$$
$$D_{c}^\dagger(x,q,f,\Gamma)=\bar{D}_{c}(x,f,q,\gamma_4\Gamma^\dagger\gamma_4)$$

- Anti-diquark

$$\bar{D}_{c}^\dagger(x,q,f,\Gamma)=\epsilon_{abc}f^T_{\alpha,a}(x)(C\gamma_4\Gamma^{\dagger}\gamma_4)_{\alpha\beta}q_{\beta,b}(x)$$
$$\bar{D}_{c}^\dagger(x,q,f,\Gamma)=D_{c}(x,f,q,\gamma_4 \Gamma^\dagger\gamma_4)$$

Here we used the property $C^\dagger=\gamma_4C\gamma_4$.

Most local lattice QCD interpolating field operators can be constructed from
these elementary blocks. For example:

- A local proton interpolating field operator can be constructed by a diquark
  and a quark:

$$[\chi_p(x)]_{\gamma} = \epsilon_{abc}u^T_{\alpha,a}(x)C\gamma_5d_{\beta,b}(x)u_{\gamma,c}(x)$$
$$[\chi_p(x)]_{\gamma} = D_c(x,u,d,\gamma_5)Q_{\gamma,c}(x,u,\gamma_0)$$

- A local $T_{cc}^+$ interpolating field operator with $J^P=1^+$ can be
  constructed by a diquark and an anti-diquark:

$$\chi_{T_{cc}^+}(x) = \epsilon_{abc}\epsilon_{dec}c^T_{\alpha,a}(x)C{\gamma_\mu}c_{\beta,b}(x)\bar{u}_{\gamma,d}(x)(\gamma_5C)_{\gamma\delta}\bar{d}^T_{\delta,e}(x)$$
$$\chi_{T_{cc}^+}(x) = D_c(x,c,c,\gamma_\mu)\bar{D}_c(x,u,d,\gamma_5)$$

Sometimes, we also need to use non-local interpolating field operators, which
can be constructed by the above blocks after replacing local quark fields with
gauge covariant shifted quark fields. For example, the $\mu$ and $-\mu$ shifted
quark field is defined as

$$q^{\mu}_{\alpha,a}(x)=[U_\mu(x)]_{ab}q_{\alpha,b}(x+\hat{\mu}),\;q^{-\mu}_{\alpha,a}(x)=[U^\dagger_\mu(x-\hat{\mu})]_{ab}q_{\alpha,b}(x-\hat{\mu})$$

where $\hat{\mu}$ is the unit vector in $\mu$ direction. Now a non-local meson
interpolating field operator can be defined as

$$M^\mu(x,q,f,\Gamma)=\bar{q}_{\alpha,a}(x)\Gamma_{\alpha\beta}f^\mu_{\beta,a}(x)=M(x,q,f^\mu,\Gamma)$$

## Propagators

A two-point correlation function of pion is

$$M(x,q,f,\gamma_5)M^\dagger(y,q,f,\gamma_5) = \bar{q}_{\alpha,a}(x)(\gamma_5)_{\alpha\beta}f_{\beta,a}(x)\bar{f}_{\alpha',a'}(x)(\gamma_4\gamma_5^\dagger\gamma_4)_{\alpha'\beta'}q_{\beta',a'}(x)$$

Define a propagator as

$$S^{q}_{\alpha\beta,ab}(x,y) = q_{\alpha,a}(x)\bar{q}_{\beta,b}(y)$$

The order of quark and anti-quark field is important. Finally, the two-point
correlation function can be expressed in terms of propagators as

$$M(x,q,f,\gamma_5)M^\dagger(y,q,f,\gamma_5) = -S^q_{\beta'\alpha,a'a}(y,x)(\gamma_5)_{\alpha\beta}S^f_{\beta\alpha',aa'}(x,y)(\gamma_4\gamma_5^\dagger\gamma_4)_{\alpha'\beta'}$$

The minus sign comes from the anti-commuting nature of the quark fields.
Finally, we usually simplify this expression by applying $\gamma_5$-Hermiticity
and degenerating two light-flavor quarks to make the coordinate order consistent
in all propagators:

$$M(x,q,f,\gamma_5)M^\dagger(y,q,f,\gamma_5) = S^{l*}_{\alpha\beta',aa'}(x,y)(\gamma_0)_{\alpha\beta}S^l_{\beta\alpha',aa'}(x,y)(\gamma_0)_{\alpha'\beta'}$$

## Adjacency Matrices

Every term in the final contraction can be represented by an adjacency matrix.
Each row of the matrix corresponds to a quark field, and each column corresponds
to an anti-quark field in the original correlation function. Thus, an edge at
position $(i,j)$ in the adjacency matrix represents a propagator connecting the
sink quark field $i$ and the source anti-quark field $\bar{j}$. The propagator
has the spin and color indices $S_{\alpha_i\beta_j, a_i b_j}$. The residual
spin/color tensors are stored separately from the adjacency matrix.

## Python API

The active implementation lives in `wick_contraction/`.

In code, `QuarkField` stores only flavor and non-local dressing information.
The operator constructors such as `QuarkBilinear`, `Diquark`, and `Quark`
describe the flavor / gamma structure first, and are then placed at a concrete
coordinate by calling `.at(...)`. This avoids repeating the same block
definition when only the coordinate changes:

```python
from wick_contraction.gamma import GAMMA_0, GAMMA_5
from wick_contraction.index import ColorIndex, SpinIndex
from wick_contraction.operator import QuarkField, QuarkBilinear, Diquark, Quark, SpinProjector
from wick_contraction.correlator import Correlator

u = QuarkField("u")
d = QuarkField("d")
d_shift = QuarkField.shift("d", "y")
d_unlinked = QuarkField.unlinked_shift("d", "y")
d_smear = QuarkField.smearing("d", "ρ")
d_nabla = QuarkField.derivative("d", "μ")

dbar_gamma5_u = QuarkBilinear(d, u, GAMMA_5)
dbar_nabla_u = QuarkBilinear(d, d_nabla, GAMMA_5)
dbar_shift_u = QuarkBilinear(d, d_shift, GAMMA_5)
u_gamma5_d = Diquark(u, d, GAMMA_5)
u_gamma0 = Quark(u, GAMMA_0)

# the same bilinear at different coordinates
pion_sink = dbar_gamma5_u.at("x")
pion_src = dbar_gamma5_u.at("y").adjoint()

# linear combinations are supported after placing blocks at coordinates
ubar_gamma5_u = QuarkBilinear(u, u, GAMMA_5)
dbar_gamma5_d = QuarkBilinear(d, d, GAMMA_5)
eta_sink = 2**-0.5 * (ubar_gamma5_u.at("x") + dbar_gamma5_d.at("x"))
eta_src = 2**-0.5 * (ubar_gamma5_u.at("y") + dbar_gamma5_d.at("y")).adjoint()
```

The non-local field constructors currently mean:

- `QuarkField.shift("d", "y")` for the gauge-covariant shifted field `W(x,y)d(y)`
- `QuarkField.unlinked_shift("d", "y")` for the plain displaced field `d(y)`
- `QuarkField.smearing("d", "ρ")` for the smeared field `(φ_ρ)d(x)`
- `QuarkField.derivative("d", "μ")` for the derivative field `(∇_μ)d(x)`

An `Operator` is contracted by constructing a `Correlator`:

```python
op = dbar_gamma5_u.at("x") * dbar_gamma5_u.at("y").adjoint()
corr = Correlator(op)
print(corr)
```

This prints a sum of adjacency-matrix terms. Each term contains residual
spin/color tensors and an adjacency matrix `A`, whose entries are propagators
written as `flavor(sink,source)`.

You can then simplify the correlator by choosing a canonical coordinate order
for propagators, and optionally merging `u/d` into a degenerate light flavor
`l`:

```python
corr.simplify(order=["x", "y"], degenerate=True)
print(corr)
```

This is the main workflow used in `examples/adjacency_new.py`.

For baryon two-point functions, you will often want to insert a spin/parity
projector to close the source and sink spin indices before contracting:

```python
c_snk, s_snk = ColorIndex.new(), SpinIndex.new()
c_src, s_src = ColorIndex.new(), SpinIndex.new()

u_gamma5_d = Diquark(u, d, GAMMA_5)
u_gamma0 = Quark(u, GAMMA_0)

proton_snk = u_gamma5_d.at("x", c_snk) * u_gamma0.at("x", s_snk, c_snk)
proton_src = (u_gamma5_d.at("y", c_src) * u_gamma0.at("y", s_src, c_src)).adjoint()
p_plus = SpinProjector.P_plus(s_src, s_snk)

proton_corr = Correlator(proton_snk * proton_src * p_plus)
proton_corr.simplify(order=["x", "y"])
print(proton_corr)
```

The returned adjacency matrix uses the following convention: fermions are
numbered separately from right to left in the original correlator. Quarks label
the rows (sink side), anti-quarks label the columns (source side), and each
edge stores the propagator as `flavor(sink,source)`. In the residual tensors,
row-side color / spin indices are renamed to `a0, a1, ...` / `α0, α1, ...`,
while column-side indices become `b0, b1, ...` / `β0, β1, ...`.
