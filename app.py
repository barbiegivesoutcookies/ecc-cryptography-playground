"""
app.py
======

Interactive Streamlit playground for the library-free ECC engine.

Features
--------
1. Curve selection: a small visualisation prime field (user-tunable a, b, p)
   or the standard secp256k1 domain parameters.
2. Enumeration of every valid discrete coordinate on the selected curve and
   rendering as a native Streamlit scatter chart over the finite field.
3. A step-by-step interactive Elliptic Curve Diffie-Hellman (ECDH) handshake
   between Alice and Bob: private scalar selection, public point derivation,
   public key exchange, and verification of the identical shared secret.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from ecc_engine import (
    CURVE_SECP256K1,
    ECCPoint,
    ECDHExchange,
    ECDHResult,
    EllipticCurve,
    enumerate_curve_points,
)

# ---------------------------------------------------------------------------
# Page configuration and session-state initialisation
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="ECC Playground",
    page_icon=":material/key:",
    layout="wide",
)

st.session_state.setdefault("ecdh_step", 1)
st.session_state.setdefault("alice_private", 7)
st.session_state.setdefault("bob_private", 11)
st.session_state.setdefault("ecdh_result", None)


def reset_handshake() -> None:
    """Reset the ECDH workflow to its initial step."""
    st.session_state["ecdh_step"] = 1
    st.session_state["ecdh_result"] = None


# ---------------------------------------------------------------------------
# Cached curve enumeration
# ---------------------------------------------------------------------------


@st.cache_data(ttl="1h", max_entries=32, show_spinner=False)
def build_curve(a: int, b: int, p: int) -> EllipticCurve:
    """Construct an :class:`EllipticCurve` (cached by parameter tuple)."""
    return EllipticCurve(a=a, b=b, p=p, name=f"toy-p{p}")


@st.cache_data(ttl="1h", max_entries=32, show_spinner="Enumerating curve points...")
def enumerate_points_as_frame(a: int, b: int, p: int) -> pd.DataFrame:
    """Enumerate all affine points on y^2 = x^3 + ax + b (mod p).

    Returns:
        A DataFrame with columns ``x`` and ``y`` covering every valid
        discrete coordinate, sorted lexicographically.
    """
    curve: EllipticCurve = build_curve(a, b, p)
    points: list[tuple[int, int]] = enumerate_curve_points(curve)
    return pd.DataFrame(points, columns=["x", "y"])


# ---------------------------------------------------------------------------
# Sidebar: curve selection
# ---------------------------------------------------------------------------

st.sidebar.header("Curve selection", divider=True)

curve_mode: str = st.sidebar.segmented_control(
    "Domain parameters",
    ["Toy field (visualisation)", "secp256k1 (standard)"],
    default="Toy field (visualisation)",
    key="curve_mode",
)

toy_curve: EllipticCurve = build_curve(2, 3, 97)

if curve_mode == "secp256k1 (standard)":
    active_curve: EllipticCurve = CURVE_SECP256K1
    st.sidebar.caption(
        "secp256k1: y² = x³ + 7 over F_p with a 256-bit prime. "
        "Point enumeration is infeasible at this scale (≈2²⁵⁶ points); "
        "the scatter plot is only available for toy fields."
    )
else:
    col_a, col_b = st.sidebar.columns(2)
    coeff_a: int = col_a.number_input(
        "Coefficient a", min_value=0, max_value=500, value=2, key="coeff_a"
    )
    coeff_b: int = col_b.number_input(
        "Coefficient b", min_value=0, max_value=500, value=3, key="coeff_b"
    )
    prime_p: int = st.sidebar.number_input(
        "Prime p", min_value=11, max_value=997, value=97, key="prime_p"
    )
    try:
        active_curve = build_curve(int(coeff_a), int(coeff_b), int(prime_p))
    except ValueError as exc:
        st.sidebar.error(f"Invalid curve parameters: {exc}")
        st.stop()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

st.title("Elliptic curve cryptography playground")
st.caption(
    "Library-free ECC over prime fields: discrete point enumeration, "
    "the chord-and-tangent group law, and a step-by-step ECDH handshake."
)

# ---------------------------------------------------------------------------
# Section 1: discrete point scatter plot
# ---------------------------------------------------------------------------

st.subheader("Discrete coordinates over the finite field", divider=True)

if curve_mode == "secp256k1 (standard)":
    st.info(
        "secp256k1 contains roughly 2²⁵⁶ points — enumeration is "
        "computationally infeasible by design. Switch to a toy field to "
        "visualise the discrete point cloud.",
        icon=":material/info:",
    )
else:
    points_frame: pd.DataFrame = enumerate_points_as_frame(
        active_curve.a, active_curve.b, active_curve.p
    )
    point_count: int = len(points_frame)

    summary_left, summary_right = st.columns([3, 1])
    summary_left.markdown(
        f"Curve **y² = x³ + {active_curve.a}x + {active_curve.b} "
        f"(mod {active_curve.p})** contains **{point_count} affine points** "
        f"(plus the point at infinity)."
    )
    summary_right.metric("Affine points", point_count)

    if point_count > 0:
        st.scatter_chart(
            points_frame,
            x="x",
            y="y",
            color="#4b7bec",
            height=480,
        )
        st.caption(
            "Each marker is a solution of the curve equation in F_p. "
            "The cloud is symmetric about y = p/2 because (x, y) on the "
            "curve implies (x, −y mod p) is too."
        )
    else:
        st.warning(
            "The selected parameters yield no affine points; adjust a, b, "
            "or p.",
            icon=":material/warning:",
        )

# ---------------------------------------------------------------------------
# Section 2: step-by-step ECDH handshake
# ---------------------------------------------------------------------------

st.subheader("Interactive ECDH handshake", divider=True)

# Generator selection.
if curve_mode == "secp256k1 (standard)":
    generator_point: ECCPoint = active_curve.generator
    st.markdown(
        "Base point **G** is the standard secp256k1 generator "
        "(a 256-bit point of prime order n)."
    )
else:
    points_frame = enumerate_points_as_frame(
        active_curve.a, active_curve.b, active_curve.p
    )
    if points_frame.empty:
        st.warning(
            "No affine points exist, so no generator can be chosen.",
            icon=":material/warning:",
        )
        st.stop()
    point_labels: list[str] = [
        f"({row.x}, {row.y})"
        for row in points_frame.itertuples(index=False)
    ]
    selected_label: str = st.selectbox(
        "Base point G",
        point_labels,
        index=0,
        key="generator_select",
    )
    gx_str, gy_str = selected_label.strip("()").split(",")
    generator_point = ECCPoint(int(gx_str), int(gy_str), active_curve)

exchange: ECDHExchange = ECDHExchange(active_curve, generator=generator_point)
bound: int = exchange.order_bound

step: int = int(st.session_state["ecdh_step"])

# Progress indicator across the four protocol stages.
st.progress(min(step, 4) / 4, text=f"Protocol stage {min(step, 4)} of 4")

# -- Stage 1: private scalar selection --------------------------------------

stage_one, stage_two = st.columns(2)

with stage_one:
    with st.container(border=True):
        st.markdown("**Stage 1 — Alice selects a private scalar**")
        alice_input: int = st.number_input(
            "Alice's private scalar d_A",
            min_value=1,
            max_value=bound,
            value=int(st.session_state["alice_private"]),
            key="alice_scalar_input",
        )
        if st.button(
            "Randomise Alice's scalar",
            icon=":material/casino:",
            key="alice_random",
        ):
            alice_input = exchange.generate_private_scalar()
            st.session_state["alice_private"] = alice_input
            st.rerun()
        st.session_state["alice_private"] = int(alice_input)
        st.caption(f"d_A = {st.session_state['alice_private']} (secret)")

with stage_two:
    with st.container(border=True):
        st.markdown("**Stage 2 — Bob selects a private scalar**")
        bob_input: int = st.number_input(
            "Bob's private scalar d_B",
            min_value=1,
            max_value=bound,
            value=int(st.session_state["bob_private"]),
            key="bob_scalar_input",
        )
        if st.button(
            "Randomise Bob's scalar",
            icon=":material/casino:",
            key="bob_random",
        ):
            bob_input = exchange.generate_private_scalar()
            st.session_state["bob_private"] = bob_input
            st.rerun()
        st.session_state["bob_private"] = int(bob_input)
        st.caption(f"d_B = {st.session_state['bob_private']} (secret)")

if st.button(
    "Continue to public key derivation",
    icon=":material/arrow_forward:",
    type="primary",
    key="advance_stage_1",
):
    st.session_state["ecdh_step"] = 2
    st.rerun()

# -- Stage 2: public point derivation ----------------------------------------

if step >= 2:
    d_a: int = int(st.session_state["alice_private"])
    d_b: int = int(st.session_state["bob_private"])
    q_a: ECCPoint = exchange.compute_public_key(d_a)
    q_b: ECCPoint = exchange.compute_public_key(d_b)

    pub_left, pub_right = st.columns(2)
    with pub_left:
        with st.container(border=True):
            st.markdown("**Stage 3 — Alice publishes Q_A = d_A · G**")
            st.code(f"Q_A = {q_a}", language=None)
    with pub_right:
        with st.container(border=True):
            st.markdown("**Stage 4 — Bob publishes Q_B = d_B · G**")
            st.code(f"Q_B = {q_b}", language=None)

    st.caption(
        "Both public points are computed with the O(log k) double-and-add "
        "ladder. They may now be transmitted over an insecure channel."
    )

    if st.button(
        "Exchange public keys and derive the shared secret",
        icon=":material/sync_alt:",
        type="primary",
        key="advance_stage_2",
    ):
        result: ECDHResult = exchange.run_exchange(
            alice_scalar=d_a, bob_scalar=d_b
        )
        st.session_state["ecdh_result"] = result
        st.session_state["ecdh_step"] = 3
        st.rerun()

# -- Stage 3: shared secret verification --------------------------------------

result_state = st.session_state.get("ecdh_result")
if step >= 3 and result_state is not None:
    shared: ECCPoint = result_state.shared_point

    secret_left, secret_right = st.columns(2)
    with secret_left:
        with st.container(border=True):
            st.markdown("**Alice computes S = d_A · Q_B**")
            st.code(f"S = {shared}", language=None)
    with secret_right:
        with st.container(border=True):
            st.markdown("**Bob computes S = d_B · Q_A**")
            st.code(f"S = {shared}", language=None)

    st.success(
        "Both parties derived the identical shared point "
        "S = d_A · d_B · G — the ECDH invariant holds.",
        icon=":material/check_circle:",
    )

    st.markdown("**Derived 256-bit symmetric key (SHA-256 of S.x)**")
    st.code(result_state.shared_secret_hex, language=None)

    meta_left, meta_right, meta_third = st.columns(3)
    meta_left.metric("Alice scalar bits", result_state.alice_private_bits)
    meta_right.metric("Bob scalar bits", result_state.bob_private_bits)
    meta_third.metric("Shared point", "Infinity" if shared.is_infinity else "Finite")

    if shared.is_infinity:
        st.error(
            "The shared point is the point at infinity — this scalar pair "
            "is degenerate for the chosen generator. Pick different scalars.",
            icon=":material/error:",
        )

    if st.button(
        "Restart handshake",
        icon=":material/restart_alt:",
        key="restart_handshake",
    ):
        reset_handshake()
        st.rerun()

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------

st.divider()
st.caption(
    "Educational demonstration only. Toy-field scalars are trivially "
    "brute-forceable; production systems must use standardised curves "
    "with ≥256-bit primes."
)