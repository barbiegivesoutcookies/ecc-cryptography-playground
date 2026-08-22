# ECC Playground

A production-ready, **library-free** Elliptic Curve Cryptography (ECC) playground.
All cryptographic primitives — arbitrary-precision modular arithmetic, the
Extended Euclidean Algorithm, the chord-and-tangent group law, double-and-add
scalar multiplication, and the ECDH key exchange — are implemented from scratch
in pure Python. No third-party cryptography libraries are used.

## Project structure

| File | Purpose |
|------|---------|
| `ecc_engine.py` | Core ECC engine: modular arithmetic, `ECCPoint`, group law, O(log k) scalar multiplication, ECDH protocol |
| `app.py` | Interactive Streamlit interface: discrete point scatter plot and step-by-step ECDH handshake simulation |
| `test_ecc.py` | Native `unittest` suite: group axioms (identity, inverse, closure) plus 100 randomized ECDH agreement runs |
| `requirements.txt` | Runtime dependencies (`numpy`, `streamlit`) |
| `Chapter_3_Mathematical_Theory.md` | Dissertation chapter: finite fields, Abelian groups, slope equations, ECDLP as a trapdoor function |
| `Chapter_4_Performance_Analysis.md` | Dissertation chapter: Big-O analysis of double-and-add vs. brute force, BSGS, and Pollard's Rho |

## Prerequisites

- Python 3.10 or newer (the code uses modern type-hint syntax)
- Windows PowerShell / Command Prompt (commands below use `py`; on macOS/Linux substitute `python3`)

## Setup

### 1. Navigate to the project directory

```powershell
cd C:\Users\hp\Desktop\ecc_project
```

### 2. Create the virtual environment

```powershell
py -m venv .venv
```

### 3. Activate the virtual environment

PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Command Prompt:

```cmd
.venv\Scripts\activate.bat
```

> If PowerShell blocks script execution, run:
> `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`
> for the current session only.

### 4. Install dependencies

```powershell
pip install -r requirements.txt
```

## Usage

### Run the interactive playground

```powershell
streamlit run app.py
```

The app opens at `http://localhost:8501`. From there you can:

1. Select a small toy prime field (tunable `a`, `b`, `p`) or standard secp256k1.
2. View every valid discrete coordinate of the curve rendered as a scatter plot over F_p × F_p.
3. Walk through the four-stage ECDH handshake: private scalar selection → public point derivation → public key exchange → identical shared secret verification.

### Execute the test suite

```powershell
python -m unittest test_ecc -v
```

Expected outcome: all tests pass, including the 100-iteration randomized ECDH
agreement test on secp256k1.

### Run the engine standalone demo

```powershell
python ecc_engine.py
```

Prints toy-curve scalar multiples and a deterministic secp256k1 ECDH exchange
to the terminal.

## Security notice

This project is an educational demonstration. The toy field is trivially
brute-forceable by design, and the symmetric-key derivation uses plain SHA-256
rather than a full KDF such as HKDF. Do not deploy this code in production;
use audited libraries (e.g., `cryptography`) with standardized curves instead.

## License

MIT