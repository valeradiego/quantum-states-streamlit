# ⚛️ Quantum States Explorer

Interactive Streamlit app built from the supplied merged Jupyter notebook. Explore angular spherical harmonics, hydrogen-atom radial probabilities, and energy eigenvalues of an infinite spherical potential well.

## Features

**Spherical harmonics**
- Choose `l = 0,...,4` and `m = -l,...,l`.
- Rotate an XYZ unit sphere, colored by **absolute magnitude** `|Y_l^m(θ, φ)|` (not by radial deformation).
- Uses the same color limits across selections; SciPy evaluates the harmonic.

**Radial probability & energy**
- **Hydrogen:** `n = 1,...,5`, `l = 0,...,n-1`; choose quantum numbers and view `a₀Pₙₗ` against `r/a₀` and `Eₙ = -13.6057/n² eV`. `Z=1`.
- **Infinite spherical well:** `n = 1,...,4`, `l = 0,...,4`, adjustable radius `a/a₀`; view `aPₙₗ` against `r/a`, the Bessel zero `xₙₗ`, and `Eₙₗ = ℏ²xₙₗ² / (2mₑa²)` in eV.
- Numerical normalization check and internal radial-node markers.
- Unit conventions are defined in the app; the hard-wall particle is an electron.

## Project structure

```text
quantum-states-streamlit/
├── streamlit_app.py              # Streamlit web app (deployment entrypoint)
├── quantum_physics.py            # Numerical functions from SciPy
├── requirements.txt              # Runtime dependencies
├── .streamlit/config.toml        # UI theme
├── .github/workflows/tests.yml   # CI numerical and syntax checks
├── tests/test_quantum_physics.py # Normalization, roots, energies
├── notebooks/
│   └── merged_quantum_notebook.ipynb  # Your original notebook
├── .gitignore
└── README.md
```

## Run locally

Recommended Python version: **3.11** (3.12 also supported by these dependencies).

```bash
python -m venv .venv
# On Windows: .venv\Scripts\activate
# On macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run streamlit_app.py
```

The local app opens at `http://localhost:8501`.

## Publish with GitHub + Streamlit Community Cloud

1. Create a new empty GitHub repository, e.g. **quantum-states-streamlit**, at <https://github.com/new>. Avoid generating another README in GitHub.
2. Unzip this project and open a terminal **inside** its `quantum-states-streamlit` folder.
3. Run (replace `YOUR_USERNAME`):

   ```bash
   git init
   git add .
   git commit -m "Build interactive quantum states explorer"
   git branch -M main
   git remote add origin https://github.com/YOUR_USERNAME/quantum-states-streamlit.git
   git push -u origin main
   ```

4. Open <https://share.streamlit.io/>, sign in with GitHub, and choose **Create app** → **I have an app**.
5. Select the GitHub repository, branch **main**, entrypoint **`streamlit_app.py`**, and click **Deploy**.
6. Community Cloud creates a public-facing `*.streamlit.app` URL. Future pushes to GitHub update the app.

The ZIP is a ready-to-push repository; it is **not already deployed**, and you must create the repository in your GitHub account.

## Equations

### Spherical harmonics

`Y_l^m(θ,φ) = scipy.special.sph_harm_y(l, m, θ, φ)` (θ polar, φ azimuth).

The plotted sphere is always `x = sin θ cos φ`, `y = sin θ sin φ`, `z = cos θ` and only the color is `|Y_l^m|`.

### Radial probability

`Pₙₗ(r) = r² |Rₙₗ(r)|²`, with `∫Pₙₗ(r)dr=1`.

For the infinite well, `Rₙₗ ∝ j_l(xₙₗr/a)` with **n-th strictly positive zero** `xₙₗ` and hard wall at `r=a`.

For hydrogen, `Rₙₗ` uses the confluent hypergeometric `₁F₁(l+1-n; 2l+2; 2r/(na₀))` as in the supplied notebook.

**Important:** In the hydrogen formula the denominator is `n²`, not `m²`: the energy depends on principal quantum number `n`. In the spherical-well formula the numerator is `ℏ²xₙₗ²`, not `h xₙₗ²`.

## Tests

```bash
python -m pip install pytest
python -m pytest -q
```

Covers all **20 well states**, all **15 hydrogen states**, energy scaling, nodes, and the angular harmonic conventions.

## Implementation notes

- `numpy`, `scipy`, `plotly`, `streamlit`; does not require Jupyter or `ipywidgets` for deployment.
- Uses an infinite-mass nucleus, nonrelativistic hydrogen approximation with an approximately 13.6057 eV ground-state binding magnitude, and an electron of mass `mₑ` in the infinite well.
- The plotted hydrogen range is finite (the wavefunction extends to infinity); the normalization displayed is therefore an integral over the **shown interval**, not exactly all space.
- Works directly from numerical library functions; does not hardcode a table of Bessel zeros.
