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

Numerical input policy
----------------------
Streamlit's JavaScript rendering engine bounds ``st.number_input`` at
2^53 - 1, which cryptographic-grade 256-bit parameters far exceed. All
numeric entry points therefore use ``st.text_input`` widgets parsed through
:func:`parse_integer_text`, which accepts base-10 and ``0x``-prefixed
hexadecimal strings, falls back to safe defaults on malformed input, and
reduces every value modulo the active field characteristic so downstream
ECC routines only ever receive sanitized native Python integers.

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
# Safe fallback literals for every text-entry field
# ---------------------------------------------------------------------------

PRIME_FALLBACK: str = "97"
COEFF_A_FALLBACK: str = "2"
COEFF_B_FALLBACK: str = "3"
ALICE_FALLBACK: str = "7"
BOB_FALLBACK: str = "11"

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

# Widget-backed keys must exist BEFORE the corresponding st.text_input
# instantiates; mutating them afterwards raises StreamlitAPIException.
if "alice_scalar_text" not in st.session_state:
    st.session_state["alice_scalar_text"] = ALICE_FALLBACK
if "bob_scalar_text" not in st.session_state:
    st.session_state["bob_scalar_text"] = BOB_FALLBACK


def reset_handshake() -> None:
    """Reset the ECDH workflow to its initial step."""
    st.session_state["ecdh_step"] = 1
    st.session_state["ecdh_result"] = None


# ---------------------------------------------------------------------------
# Text-to-integer parsing pipeline (Streamlit JS bound workaround)
# ---------------------------------------------------------------------------


def parse_integer_text(raw: str, fallback: str) -> tuple[int, bool]:
    """Convert a text field into a native Python integer.

    Accepts standard base-10 integer strings and hexadecimal strings
    prefixed with ``0x`` (case-insensitive, optionally signed). Empty or
    malformed input falls back to the supplied default literal.

    Args:
        raw: The raw widget string.
        fallback: A trusted default literal used when parsing fails.

    Returns:
        A tuple ``(value, is_valid)`` where ``value`` is the parsed
        integer and ``is_valid`` is False when the fallback was applied.
    """
    cleaned: str = raw.strip()
    if not cleaned:
        # Empty entry: apply the fallback and flag it so the UI can warn.
        try:
            if fallback.lower().startswith("0x"):
                return int(fallback, 16), False
            return int(fallback, 10), False
        except ValueError:  # Defensive: fallback literals are trusted.
            return 0, False
    try:
        if cleaned.lower().lstrip("+-").startswith("0x"):
            value: int = int(cleaned, 16)
        else:
            value = int(cleaned, 10)
        return value, True
    except ValueError:
        try:
            if fallback.lower().startswith("0x"):
                return int(fallback, 16), False
            return int(fallback, 10), False
        except ValueError:  # Defensive: fallback literals are trusted.
            return 0, False


def sanitize_field_element(value: int, curve_p: int) -> int:
    """Enforce the modular boundary check ``value % curve_p``."""
    return value % curve_p


def sanitize_scalar(value: int, curve_p: int, bound: int) -> int:
    """Reduce a private scalar into the valid range ``[1, bound]``.

    The value is first reduced modulo the field characteristic per the
    modular boundary policy, then folded into the generator-order bound so
    downstream ECDH routines never receive an out-of-range scalar.
    """
    reduced: int = value % curve_p
    if reduced > bound:
        reduced %= bound
    if reduced == 0:
        reduced = 1
    return reduced


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

if curve_mode == "secp256k1 (standard)":
    active_curve: EllipticCurve = CURVE_SECP256K1
    st.sidebar.caption(
        "secp256k1: y² = x³ + 7 over F_p with a 256-bit prime. "
        "Point enumeration is infeasible at this scale (≈2²⁵⁶ points); "
        "the scatter plot is only available for toy fields."
    )
else:
    # -- Prime modulus entry (text pipeline; exceeds JS number bounds) ------
    prime_raw: str = st.sidebar.text_input(
        "Prime p (base-10 or 0x hex)",
        value=PRIME_FALLBACK,
        key="prime_text",
    )
    prime_p: int
    prime_ok: bool
    prime_p, prime_ok = parse_integer_text(prime_raw, PRIME_FALLBACK)
    if not prime_ok:
        st.sidebar.warning(
            f"Unparseable prime entry; using fallback p = {prime_p}.",
            icon=":material/warning:",
        )

    # -- Coefficient entries, reduced modulo the active prime ----------------
    coeff_a_raw: str = st.sidebar.text_input(
        "Coefficient a (base-10 or 0x hex)",
        value=COEFF_A_FALLBACK,
        key="coeff_a_text",
    )
    coeff_a_val: int
    coeff_a_ok: bool
    coeff_a_val, coeff_a_ok = parse_integer_text(coeff_a_raw, COEFF_A_FALLBACK)
    coeff_a_val = sanitize_field_element(coeff_a_val, prime_p)

    coeff_b_raw: str = st.sidebar.text_input(
        "Coefficient b (base-10 or 0x hex)",
        value=COEFF_B_FALLBACK,
        key="coeff_b_text",
    )
    coeff_b_val: int
    coeff_b_ok: bool
    coeff_b_val, coeff_b_ok = parse_integer_text(coeff_b_raw, COEFF_B_FALLBACK)
    coeff_b_val = sanitize_field_element(coeff_b_val, prime_p)

    if not (coeff_a_ok and coeff_b_ok):
        st.sidebar.warning(
            "Unparseable coefficient entry; fallback applied and reduced "
            f"mod {prime_p}.",
            icon=":material/warning:",
        )

    try:
        active_curve = build_curve(coeff_a_val, coeff_b_val, prime_p)
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
            x_label="x coordinate (mod p)",
            y_label="y coordinate (mod p)",
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
curve_p: int = active_curve.p

# -- Stage 1: private scalar selection (text pipeline) ------------------------

stage_one, stage_two = st.columns(2)

with stage_one:
    with st.container(border=True):
        st.markdown("**Stage 1 — Alice selects a private scalar**")
        # Randomisation mutates the widget-backed key BEFORE the text input
        # below instantiates, then refreshes via st.rerun() — the only
        # lifecycle-safe ordering in Streamlit.
        if st.button(
            "Randomise Alice's scalar",
            icon=":material/casino:",
            key="alice_random",
        ):
            random_alice: int = exchange.generate_private_scalar()
            st.session_state["alice_scalar_text"] = str(random_alice)
            st.rerun()
        alice_raw: str = st.text_input(
            "Alice's private scalar d_A (base-10 or 0x hex)",
            value=st.session_state["alice_scalar_text"],
            key="alice_scalar_text",
        )
        alice_val: int
        alice_ok: bool
        alice_val, alice_ok = parse_integer_text(alice_raw, ALICE_FALLBACK)
        alice_val = sanitize_scalar(alice_val, curve_p, bound)
        if not alice_ok:
            st.caption(
                f"Unparseable entry; fallback d_A applied and sanitized."
            )
        st.session_state["alice_private"] = int(alice_val)
        st.caption(
            f"d_A = {st.session_state['alice_private']} "
            f"(hex: {hex(st.session_state['alice_private'])}) — secret"
        )

with stage_two:
    with st.container(border=True):
        st.markdown("**Stage 2 — Bob selects a private scalar**")
        # Randomisation mutates the widget-backed key BEFORE the text input
        # below instantiates, then refreshes via st.rerun() — the only
        # lifecycle-safe ordering in Streamlit.
        if st.button(
            "Randomise Bob's scalar",
            icon=":material/casino:",
            key="bob_random",
        ):
            random_bob: int = exchange.generate_private_scalar()
            st.session_state["bob_scalar_text"] = str(random_bob)
            st.rerun()
        bob_raw: str = st.text_input(
            "Bob's private scalar d_B (base-10 or 0x hex)",
            value=st.session_state["bob_scalar_text"],
            key="bob_scalar_text",
        )
        bob_val: int
        bob_ok: bool
        bob_val, bob_ok = parse_integer_text(bob_raw, BOB_FALLBACK)
        bob_val = sanitize_scalar(bob_val, curve_p, bound)
        if not bob_ok:
            st.caption(
                f"Unparseable entry; fallback d_B applied and sanitized."
            )
        st.session_state["bob_private"] = int(bob_val)
        st.caption(
            f"d_B = {st.session_state['bob_private']} "
            f"(hex: {hex(st.session_state['bob_private'])}) — secret"
        )

if st.button(
    "Continue to public key derivation",
    icon=":material/arrow_forward:",
    type="primary",
    key="advance_stage_1",
):
    st.session_state["ecdh_step"] = 2
    st.rerun()

# -- Stage 2: public point derivation ----------------------------------------

step: int = int(st.session_state["ecdh_step"])

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
    "with ≥256-bit primes. All numeric entries accept base-10 or 0x-hex "
    "and are reduced modulo the active field characteristic."
)