import streamlit as st
import plotly.graph_objects as go
import pandas as pd
from rdkit import Chem

# local package import
from green_core import GreenAssess

st.set_page_config(page_title="Enviro-Mol Demo", page_icon="🌿", layout="wide")
st.title("Enviro-Mol · Çevre Zarar Tahmini (Demo)")

st.markdown(
    "Girilen molekül için biyobozunurluk, biyobirikim, akut toksisite vb. parametreleri hesaplar ve 0–100 arası **Enviro-Risk Index (ERI)** üretir. "
)

# -----------------------------------------------------------------------------
# Input section
# -----------------------------------------------------------------------------

smiles = st.text_input("**SMILES** girin (örn. CCO, DDT = 'Clc1ccc(c(c1)Cl)C(Cl)(Cl)Cl'):")

col1, col2 = st.columns(2)

if col1.button("Hesapla") and smiles:
    try:
        ga = GreenAssess(smiles)
    except ValueError as err:
        st.error(f"SMILES hatası: {err}")
        st.stop()

    # ------------------------------------------------------------------
    # Tablo
    # ------------------------------------------------------------------
    st.subheader("Parametre Sonuçları")
    df = ga.to_dataframe()
    st.dataframe(df, use_container_width=True)

    # ------------------------------------------------------------------
    # Radial (spider) chart
    # ------------------------------------------------------------------
    st.subheader("Radar Grafiği")
    radar_fig = go.Figure()
    categories = list(ga.scores.keys())
    radar_fig.add_trace(
        go.Scatterpolar(
            r=list(ga.scores.values()) + [list(ga.scores.values())[0]],
            theta=categories + [categories[0]],
            fill="toself",
        )
    )
    radar_fig.update_layout(polar=dict(radialaxis=dict(range=[0, 1])), showlegend=False)
    st.plotly_chart(radar_fig, use_container_width=True)

    # ------------------------------------------------------------------
    # Leaf Chart – simple prototype
    # ------------------------------------------------------------------
    st.subheader("Leaf Chart (hidrofobiklik × kalıcılık)")

    def leaf_chart(log_kow: float, t_half: float):
        # Prototype: draw three regions as background shapes
        fig = go.Figure()
        # Low-risk base leaf (polygon approx.)
        fig.add_shape(type="path", path="M0,0 Q1,2 2,0 Z", fillcolor="#d5e8d4", line_color="#d5e8d4", opacity=0.4)
        # Mid-risk midrib
        fig.add_shape(type="line", x0=0, y0=0, x1=2, y1=0, line=dict(color="#82b366", width=3, dash="dot"))
        # High-risk outline (placeholder rectangle)
        fig.add_shape(type="rect", x0=2, y0=0, x1=4, y1=4, fillcolor="#e6b8af", line=dict(color="#e6b8af"), opacity=0.3)
        # Data point
        fig.add_trace(go.Scatter(x=[log_kow], y=[t_half], mode="markers", marker=dict(size=12, color="#333333")))
        fig.update_layout(xaxis_title="log Kow", yaxis_title="t½ (gün)", showlegend=False, width=500, height=500)
        return fig

    leaf_fig = leaf_chart(ga.descriptors["logKow"], 15 + 165 * ga.scores["biodeg"])  # crude mapping
    st.plotly_chart(leaf_fig, use_container_width=True)

    # ------------------------------------------------------------------
    # Gauge
    # ------------------------------------------------------------------
    st.subheader("Enviro-Risk Index (ERI)")
    gauge = go.Figure(go.Indicator(
        mode="gauge+number",
        value=ga.eri,
        gauge={
            "axis": {"range": [0, 100]},
            "bar": {"color": "#2c5234"},
            "steps": [
                {"range": [0, 33], "color": "#c6e0b4"},
                {"range": [33, 66], "color": "#ffe699"},
                {"range": [66, 100], "color": "#f4b084"},
            ],
        },
        number={"suffix": "/100"},
    ))
    gauge.update_layout(height=300)
    st.plotly_chart(gauge, use_container_width=False)

else:
    col2.info("SMILES girip **Hesapla**'ya basın.")

