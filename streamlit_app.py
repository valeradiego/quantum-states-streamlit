"""Interactive Quantum States Explorer for Streamlit Community Cloud.

Run with: streamlit run streamlit_app.py
"""

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from scipy.integrate import trapezoid

from quantum_physics import (
    HYDROGEN_RYDBERG_EV,
    hydrogen_dimensionless_density,
    hydrogen_energy_ev,
    nth_spherical_bessel_zero,
    spherical_harmonic_magnitude,
    well_dimensionless_density,
    well_energy_ev,
)

st.set_page_config(
    page_title="Quantum States Explorer",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
      .block-container { max-width: 1220px; padding-top: 2.2rem; }
      h1, h2, h3 { letter-spacing: -0.02em; }
      [data-testid="stMetricValue"] { font-size: 1.65rem; }
      .eyebrow { font-size: 0.76rem; letter-spacing: .16em; font-weight: 700; 
                  text-transform: uppercase; color: #38bdf8; }
      .small-note { font-size: .86rem; opacity: .78; line-height: 1.5; }
    </style>
    """,
    unsafe_allow_html=True,
)

PLOT_BG = "rgba(0,0,0,0)"
GRID_COLOR = "rgba(148,163,184,0.20)"
ACCENT = "#38bdf8"


@st.cache_data(show_spinner=False)
def angular_grid():
    polar = np.linspace(0.0, np.pi, 110)
    azimuth = np.linspace(0.0, 2.0 * np.pi, 220)
    theta, phi = np.meshgrid(polar, azimuth, indexing="ij")
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return theta, phi, x, y, z


def create_harmonic_figure(l, m):
    theta, phi, x, y, z = angular_grid()
    magnitude = spherical_harmonic_magnitude(l, m, theta, phi)
    global_color_max = np.sqrt((2 * 4 + 1) / (4.0 * np.pi))

    fig = go.Figure(go.Surface(
        x=x, y=y, z=z,
        surfacecolor=magnitude,
        colorscale="Viridis",
        cmin=0, cmax=global_color_max,
        colorbar=dict(title=dict(text="|Y|"), thickness=14, len=0.75),
        hovertemplate=(
            "x=%{x:.3f}<br>y=%{y:.3f}<br>z=%{z:.3f}"
            "<br>|Y|=%{surfacecolor:.5f}<extra></extra>"
        ),
    ))
    fig.update_layout(
        margin=dict(l=0, r=0, t=12, b=0),
        height=570,
        paper_bgcolor=PLOT_BG,
        font=dict(color="#a4b3c7"),
        scene=dict(
            xaxis=dict(title="x", range=[-1.15, 1.15],
                       backgroundcolor=PLOT_BG, gridcolor=GRID_COLOR),
            yaxis=dict(title="y", range=[-1.15, 1.15],
                       backgroundcolor=PLOT_BG, gridcolor=GRID_COLOR),
            zaxis=dict(title="z", range=[-1.15, 1.15],
                       backgroundcolor=PLOT_BG, gridcolor=GRID_COLOR),
            aspectmode="cube",
            camera=dict(eye=dict(x=1.65, y=1.5, z=1.18)),
        ),
    )
    return fig


def create_radial_figure(radius, density, xlabel, ylabel, title, nodes):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=radius, y=density,
        mode="lines", name="Radial probability",
        line=dict(color=ACCENT, width=3),
        fill="tozeroy", fillcolor="rgba(56,189,248,0.10)",
        hovertemplate="Radius: %{x:.4f}<br>Density: %{y:.6f}<extra></extra>",
    ))
    if nodes:
        for root in nodes:
            fig.add_vline(
                x=float(root), line_width=1,
                line_dash="dot", line_color="rgba(251,191,36,0.65)",
            )
    fig.update_layout(
        title=dict(text=title, font=dict(size=16)),
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        paper_bgcolor=PLOT_BG,
        plot_bgcolor=PLOT_BG,
        font=dict(color="#a4b3c7"),
        hovermode="x unified",
        margin=dict(l=14, r=10, t=48, b=12),
        height=480,
        showlegend=False,
        xaxis=dict(range=[float(radius[0]), float(radius[-1])],
                   showgrid=True, gridcolor=GRID_COLOR, zeroline=False),
        yaxis=dict(rangemode="tozero", showgrid=True, gridcolor=GRID_COLOR, zeroline=False),
    )
    return fig


def hydrogen_nodes(n, l):
    """Analytic radial nodes from generalized Laguerre polynomial zeros."""
    from scipy.special import roots_genlaguerre

    count = n - l - 1
    if count == 0:
        return []
    rho_nodes, _ = roots_genlaguerre(count, 2*l + 1)
    return list(n * rho_nodes / 2.0)


def well_nodes(n, l):
    xnl = nth_spherical_bessel_zero(l, n)
    return [nth_spherical_bessel_zero(l, k) / xnl for k in range(1, n)]


with st.sidebar:
    st.markdown('<p class="eyebrow">QUANTUM PHYSICS LAB</p>', unsafe_allow_html=True)
    st.markdown("### Quantum States Explorer")
    st.caption("Interactive stationary-state solutions of the Schrödinger equation.")
    st.divider()
    st.markdown("**What you can explore**")
    st.markdown("• Spherical harmonics, on a true unit sphere")
    st.markdown("• Hydrogen-atom radial probabilities")
    st.markdown("• Infinite spherical-well radial probabilities")
    st.markdown("• Energy eigenvalues in electronvolts")
    st.divider()
    st.caption("Based on the supplied Jupyter notebook. Built with SciPy + Plotly + Streamlit.")

st.markdown('<p class="eyebrow">INTERACTIVE VISUALIZATION</p>', unsafe_allow_html=True)
st.title("Quantum States Explorer")
st.write("Explore angular wavefunctions, radial probability distributions, and their energy eigenvalues.")

tab_angular, tab_radial = st.tabs(["◉  Spherical harmonics", "◌  Radial probability & energy"])

with tab_angular:
    st.subheader("Spherical harmonics")
    st.markdown(
        "The surface stays a **unit sphere**. Only the color varies with the magnitude "
        r"$|Y_l^m(\theta,\phi)|$."
    )
    controls, plot = st.columns([1, 2.5], gap="large")
    with controls:
        st.markdown("#### Quantum numbers")
        l_ang = st.slider("Orbital quantum number  l", min_value=0, max_value=4, value=2, step=1)
        m_ang = st.select_slider(
            "Magnetic quantum number  m",
            options=list(range(-l_ang, l_ang + 1)),
            value=0,
            key=f"angular_m_{l_ang}",
        )
        st.metric("Selected harmonic", f"l = {l_ang}, m = {m_ang}")
        st.markdown("**Geometry**")
        st.latex(r"x=\sin\theta\cos\phi,\quad y=\sin\theta\sin\phi,\quad z=\cos\theta")
        st.markdown(
            '<p class="small-note">Color limits are fixed for all l ≤ 4 so that '
            'different harmonics can be compared directly. '
            'The absolute value is identical for +m and −m.</p>',
            unsafe_allow_html=True,
        )
    with plot:
        st.plotly_chart(create_harmonic_figure(l_ang, m_ang), use_container_width=True,
                        config={"displaylogo": False})
    st.caption("Spherical harmonics describe the angular part of the wavefunction; by themselves they do not determine an energy.")

with tab_radial:
    st.subheader("Radial probability density & energy")
    st.markdown(r"Plot the normalized $P_{nl}(r)=r^2|R_{nl}(r)|^2$ and its energy eigenvalue.")

    potential = st.radio("Choose a physical system", ["Hydrogen atom", "Infinite spherical well"],
                         horizontal=True, label_visibility="visible")

    settings, graph = st.columns([1, 2.5], gap="large")
    with settings:
        st.markdown("#### State parameters")
        if potential == "Hydrogen atom":
            n = st.slider("Principal quantum number  n", 1, 5, 3, key="hydrogen_n")
            l = st.selectbox("Orbital quantum number  l", options=list(range(n)),
                             index=min(1, n-1), key=f"hydrogen_l_n{n}")
            radius = np.linspace(0, 2.0*n*(n+4), 2400)
            density = hydrogen_dimensionless_density(radius, n, l)
            energy = hydrogen_energy_ev(n)
            internal_nodes = hydrogen_nodes(n, l)
            xmax = None
            label_x, label_y = "r / a₀", "a₀ Pₙₗ(r)"
            plot_title = f"Hydrogen · n = {n}, l = {l}"
        else:
            n = st.slider("Radial quantum number  n", 1, 4, 2, key="well_n")
            l = st.slider("Orbital quantum number  l", 0, 4, 1, key="well_l")
            a_bohr = st.number_input(
                "Well radius  a / a₀", min_value=0.2, max_value=25.0,
                value=1.0, step=0.2, format="%.2f",
                help="The radius a is given in Bohr radii. The plot uses r/a.",
            )
            radius = np.linspace(0, 1, 1800)
            density = well_dimensionless_density(radius, n, l)
            energy = well_energy_ev(n, l, a_bohr)
            internal_nodes = well_nodes(n, l)
            xmax = nth_spherical_bessel_zero(l, n)
            label_x, label_y = "r / a", "a Pₙₗ(r)"
            plot_title = f"Spherical well · n = {n}, l = {l}"

        st.metric("Energy eigenvalue", f"{energy:,.5f} eV")
        if xmax is not None:
            st.metric("Bessel zero  xₙₗ", f"{xmax:.6f}")
        st.metric("Internal radial nodes", f"{len(internal_nodes)}")

    with graph:
        fig = create_radial_figure(radius, density, label_x, label_y,
                                   plot_title, internal_nodes)
        st.plotly_chart(fig, use_container_width=True,
                        config={"displaylogo": False})
        integral = float(trapezoid(density, radius))
        st.caption(
            f"Numerical probability within the displayed interval: {integral:.7f}. "
            "Dotted lines mark internal radial nodes."
        )

    st.divider()
    with st.expander("Equations and conventions", expanded=True):
        if potential == "Hydrogen atom":
            st.latex(r"E_n=-\frac{13.6057\;\mathrm{eV}}{n^2}")
            st.latex(
                r"R_{nl}(r)=\frac{1}{(2l+1)!}"
                r"\left(\frac{2r}{na_0}\right)^l e^{-r/(na_0)}"
                r"\sqrt{\left(\frac{2}{na_0}\right)^3"
                r"\frac{(n+l)!}{2n(n-l-1)!}}\,"
                r"{}_1F_1\!\left(l+1-n;2l+2;\frac{2r}{na_0}\right)"
            )
            st.markdown(
                "The app uses the confluent hypergeometric function via "
                "`scipy.special.hyp1f1`, with **Z = 1** and the "
                "nonrelativistic fixed-nucleus approximation. The energy "
                "depends on **n**, not **l** or **m**. The radial axis is r/a₀."
            )
        else:
            st.latex(r"E_{nl}=\frac{\hbar^2x_{nl}^2}{2m_ea^2},\qquad j_l(x_{nl})=0")
            st.latex(
                r"R_{nl}(r)=\frac{\sqrt{2}}{a^{3/2}|j_{l+1}(x_{nl})|}"
                r"j_l\!\left(\frac{x_{nl}r}{a}\right)"
            )
            st.markdown(
                "An **electron** of mass mₑ in a hard-wall spherical cavity "
                "with radius a. Here n indexes strictly positive spherical-Bessel "
                "zeros (it is not the principal quantum number of hydrogen)."
            )
        st.latex(r"\int P_{nl}(r)\,dr=1")
        st.caption("Axis labels are dimensionless densities: a₀P versus r/a₀ for hydrogen, and aP versus r/a for the well.")

st.divider()
st.caption("Numerical functions: SciPy · Charts: Plotly · Interface: Streamlit · Source notebook included in the repository")
