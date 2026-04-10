# Jones index of prime endomorphisms in the Bost–Connes system: higher invariants, fusion rules, and logarithmic corrections to black hole entropy

**Author:** E. F. Perez-Eugenio  
**ORCID:** [0009-0006-3228-4847](https://orcid.org/0009-0006-3228-4847)  
**Date:** April 2026  
**DOI:** [10.5281/zenodo.19500605](https://doi.org/10.5281/zenodo.19500605)  
**Extends:** Paper IV — DOI: [10.5281/zenodo.19216302](https://doi.org/10.5281/zenodo.19216302)

---

## Abstract

We prove that the Jones index of the prime endomorphisms ρ_p in the Bost–Connes system equals p, provide complete proofs via the Takesaki crossed product, and compute the higher invariants: abelian fusion rules (ρ_n ∘ ρ_m = ρ_nm), one-dimensional intertwiner spaces, and group-type Ocneanu invariants for p=2 (ℤ/2ℤ) and p=3 (ℤ/3ℤ). The Riemann zeta function is expressed as ζ(β) = ∏_p (1 − [M:N]_p^{−β})^{−1}. We connect the Bost–Connes subfactors with minimal models of two-dimensional conformal field theory via the ADE classification, documenting both coincidences (shared Jones index [M:N]=2 with the Ising model) and differences (inequivalent fusion categories). Motivated by the recent quantitative verification that the critical Ising chain has a BTZ black hole dual [Wang and He, arXiv:2603.07639 (2026)], we identify an open question: the logarithmic correction coefficient k_Ising to the BTZ entropy for central charge c=1/2. We predict k_Ising ≠ 3/2.

---

## Key results

| Result | Status |
|--------|--------|
| [e_p R e_p : ρ_p(R)] = p | Proven (Theorem 1) |
| Fusion rules ρ_n ∘ ρ_m = ρ_nm (abelian) | Proven (§3.1) |
| dim Hom(ρ_n, ρ_m) = δ_nm, all ρ_n irreducible | Proven (Theorem 3) |
| Ocneanu invariant p=2: ℤ/2ℤ, graph A₃ | Computed (Appendix B) |
| Ocneanu invariant p=3: ℤ/3ℤ, graph A₅ | Computed (Appendix B) |
| C*-tensor category: Vec ⊗ ℕ× | Established (§3.4) |
| BC ↔ Ising: shared index, inequivalent categories | Established (§4) |
| k_Ising ≠ 3/2 | Predicted (§5.4) |

---

## Repository structure

```
jones-bc-higher-invariants/
│
├── paper_V_jones_bc_higher_invariants.tex
├── perez-eugenio_2026_jones-bc-higher-invariants.pdf
├── README.md
├── LICENSE
├── .gitignore
├── .zenodo.json
│
├── code/
│   └── TFI_critical_OTOC_v4.py
│
└── data/
    ├── TFI_v4_all_results.json
    └── TFI_v6_c_vs_hz.json
```

---

## Numerical data — OTOC simulations

### code/TFI_critical_OTOC_v4.py

OTOC scrambling simulation for the Transverse Field Ising (TFI) model.  
Produces `data/TFI_v4_all_results.json`.

**Convention (exact OTOC project standard):**
- `C(t) = |Tr(W(t) V W(t) V)| / dim` (infinite temperature)
- W = σ_z site 0, V = σ_z site 1 (r = 1, fixed)
- `D_MAX = 5N`, late average = second half of time series
- PBC, N = 4–12 (even and odd)
- `H = −J Σ σ^z_i σ^z_{i+1} − h_x Σ σ^x_i − h_z Σ σ^z_i`

**Dependencies:** `numpy`, `scipy`

```bash
pip install numpy scipy
python TFI_critical_OTOC_v4.py
```

### data/TFI_v4_all_results.json

Four models: TFI pure (h/J=1) and TFI+h_z ∈ {0.1, 0.3, 0.5}.

| Model | c | c_err | R² | Monotone |
|-------|---|-------|----|---------|
| TFI pure (h_z=0) | — | — | — | No (parity effect) |
| h_z=0.1 | 3.073 | 0.366 | 0.963 | Yes |
| h_z=0.3 | 2.118 | 0.359 | 0.907 | Yes |
| h_z=0.5 | — | — | — | No |

### data/TFI_v6_c_vs_hz.json

c(h_z) curve with t_max=40, N=4,6,8,10, eight values h_z ∈ {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5}.

| h_z | c | c_err | R² |
|-----|---|-------|----|
| 0.0 | 4.402 | 0.127 | 0.9996 |
| 0.05 | 4.658 | 0.022 | 0.99999 |
| 0.1 | 3.033 | 0.328 | 0.9907 |
| 0.15 | 3.232 | 0.265 | 0.9950 |
| **0.2** | **3.407** | **0.089** | **0.9995** |
| 0.3 | 2.512 | 0.070 | 0.9992 |
| 0.4 | 2.592 | 0.175 | 0.9950 |
| 0.5 | 2.682 | 0.439 | 0.9708 |

**Cleanest result:** h_z=0.2 → c = 3.407 ± 0.089, R² = 0.9995.  
Coincides with [M:N] = 2+√2 ≈ 3.414 at 0.2% (within one error bar).

---

## Citation

```bibtex
@misc{PerezEugenio2026V,
  author = {Perez-Eugenio, E. F.},
  title  = {Jones index of prime endomorphisms in the {Bost--Connes} system:
            higher invariants, fusion rules, and logarithmic corrections
            to black hole entropy},
  year   = {2026},
  doi    = {10.5281/zenodo.19500605}
}
```

---

## License

[CC BY 4.0](LICENSE)
