# RAINBOR 12-billion-address sparse meta-pixel experiment — 2026-11-12

**Hosted verdict: PASS (two independent instrumentation runs).**

- [Optimized hosted test](https://github.com/JesseBrown1980/RAINBOR/actions/runs/43304245335), status PASS, complete 144-screen HBP result and SHA-256 manifest in Actions artifact `rainbon-12b-exact-jumps-and-density-experiment`.
- [Original allocation-instrumented hosted test](https://github.com/JesseBrown1980/RAINBOR/actions/runs/43303628610), also PASS. All **186 result lines are byte-identical to the optimized run except for the one RESOURCE line**, which intentionally measures memory differently.
- The optimized HBP record SHA-256 is `60c5bc244f0a1b476106ea9fedd22d641c3b5a11e6dd616824b4f9137702e830`; the original record is `ebaa9c628e420c70b920348a9b9c522a1fb422d2af5d957faba55eeb893a8b65`. Both workflow-produced SHA256SUMS files passed, and deliberate one-byte tampering was refused.
- Fully reproducible script: [rainbon_12b_jump_density.py](rainbon_12b_jump_density.py), anchored to the original `RAINBOR-DOUBLE-SLIT.html` and `RAINBOR-COLOUR-SLITS-AND-SPHERES.html` SHA256 sidecars and exact `The-matrix-law-essentials` commit `e9eb98803d0474c31aca045723d1d896a03dc8f1`. The frozen [pre-run hypotheses](RAINBON-12B-JUMP-DENSITY-FREEZE-20261112.hbp) remain intact.

## What was computed and what was not

The detector **can address 12,000,000,000 virtual meta-pixel coordinates**, but it **sampled only 243** on each of **144 detector screens**. There was no 12-billion-item state allocation, no 12-billion-particle interaction and no quantum-computer simulation. Peak optimized *process resident memory* was **20,258,816 bytes**, with the original instrumented run reporting **1,436,718 peak Python traced allocation bytes** instead (different counters; do not compare them as identical metrics).

The dynamical state is the exact 27-site integer circulant `M=I+R−R² (mod 1,000,081)`. The source is copied without mathematical changes from `The-matrix-law-essentials/python/matrix_tricks/jump3_ntt.py`. It was independently checked against the original published hashes at the six requested power-of-three checkpoints.

| Checkpoint | Integer steps | Complete-state SHA-256 prefix | Result |
| --- | ---: | --- | --- |
| `3^3` | 27 | `7d05096afc24b08c` | direct = NTT = cube |
| `3^6` | 729 | `cc3581ab163de3a7` | direct = NTT = cube |
| `3^9` | 19,683 | `ee4b40fb753713ab` | NTT = cube |
| `3^10` | 59,049 | `81323c66a0e8a92f` | NTT = cube |
| `3^20` | 3,486,784,401 | `8ef1ccc86cd3ff0c` | NTT = cube |
| `3^40` | 12,157,665,459,056,928,801 | `177066f133a5c9ef` | NTT = cube |

Freeze/play at `3^40`: **forward +1, forward +3, backward −1 PASS**. Independently evaluated `2^40` binary time jump: SHA-256 `40895bdc21e2ee5a`. Mixed `3^40+2^20+3` step via frozen checkpoint versus a single NTT evaluation: **bit-for-bit equal**, SHA-256 `83137948f6537251`.

The model *also* directly reproduces the original RAINBOR single-, double- and triple-slit normalized Fraunhofer intensity, tested at **2,001 angles for each slit count**: maximum absolute error **7.34967642302×10⁻¹⁴**.

## Artificial density, temperature, solar-wind and emission controls

The experiments use **dimensionless, deliberately assumed** source rules:

```
source_strength = 0.06 * density * max(0,1 − T/T_cap)   [cutoff enabled]
source_strength = 0.06 * density                         [cutoff disabled]
```

This strength determines a softened, off-axis logarithmic **optical phase perturbation** between slits. It is *not* derived from general relativity, a stress-energy tensor, a neutron-star equation of state, or a black-hole solution. Its temperature gate is precisely the conjecture being tested **inside this software**; nature does not currently establish that gravity disappears above a temperature cap.

Six source treatments were tested at four frozen/jumped times and one, two or three slits, with monochromatic or **409-channel synthetic broad spectrum** light. The 409 channels sample model wavelengths 8–40; they are **not** 409 independently measured Asolaria stellar GGUF wavefunctions. Two additional control calculations per slit/time/spectrum separately varied a synthetic wind-phase direction and an incoherent emission-halo term.

The following are measured **RMS differences of sampled screen intensity from the matched vacuum screen** at the `3^40` checkpoint:

| Model source | 2 slits · 409 channels | 3 slits · 409 channels |
| --- | ---: | ---: |
| Ordinary low density | see full Actions HBP | see full Actions HBP |
| Neutron-like, density=12, T/Tcap=.25 | **0.026275710** | **0.017089760** |
| Black-hole-like, density=80, T/Tcap=.65 | **0.077457548** | **0.090221216** |
| Supercritical, density=250, T/Tcap=1.20, cutoff ON | **0.000000000** | **0.000000000** |
| The same supercritical source, cutoff OFF | **0.235227249** | **0.196220755** |

For the same checkpoint, model wind reversal changes the three-slit broad-spectrum footprint by RMS **0.016521908**; switching on the separately assumed radiation halo contributes RMS **0.001306317**. For a single slit, a pure phase perturbation does not change the intensity within rounding (RMS `0.000000000`)—an independent null control.

The strongest numeric result is **exact control consistency**, not a discovery about real compact objects: under the *programmed* temperature gate, high-density phase gradients cancel and the screen becomes exactly the matched vacuum screen. Disabling that gate causes a visible change, as engineered.

## The critical negative control: nonlinear/source-driven evolution

A separate 27-step run inserts a local `+7` modular kick **after every step**. Its state differs from `M^27 x₀`, so the experiment **refuses** to fast-forward that forced dynamics using the unmodified fixed matrix. A state-dependent density interaction needs its own valid reduction theorem, enlarged linear representation, exact composition law or interval-bounded approximation before a power-of-three shortcut is justified.

Similarly, SHA-256 hashes verify that the saved states changed; hashes do **not** generate forces, gravity or physical nonlocal interactions. The slit simulations do not establish physical rainbon particles, actual gluons, glueballs, protons, black-hole gravity, Hawking radiation or quantum entanglement.

## Next falsifier

Hold the 27-channel jump law fixed and try a **nonlinear, spatially localized source that changes the wave-state on each microscopic time step**, with paired direct-versus-jumped computations. If the model cannot prove a commuting/closed reduced update rule and the two diverge, the proposed cross-scale jump must stay restricted to unperturbed segments or piecewise externally forced segments with explicit checkpoint composition. The present result shows the *readout* can be sparse and changed at 12 billion virtual addresses, not that nonlinear many-body Yang–Mills evolution has been made matrix-free by hashing.

The protected Yang–Mills coupling `r=9` computation remains **HOLD**; no such run was attempted.
