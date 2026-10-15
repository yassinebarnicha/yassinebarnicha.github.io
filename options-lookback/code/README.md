# Unified Lookback Options Pricing Library (C++)

## Overview

This project provides a **unified pricing DLL** for **European floating-strike lookback options** under the **Black–Scholes** model.

The library integrates **two independent Monte Carlo approaches** developed separately and compiled into a **single shared library**:

- **Barnicha Engine**  
  One-step Monte Carlo with exact extremum sampling (**continuous monitoring**, **no discretization bias**)

- **Chivet Engine**  
  Time-stepping Monte Carlo with optional **BGK** bias correction

Both engines are accessible through a **single unified DLL interface** and can be used interchangeably for pricing and Greeks computation.

---

# Contract Specification

European floating-strike lookback options:

- **Call payoff**  
 $$ \( S_T - \min_{0\le t\le T} S_t \) $$

- **Put payoff**  
$$  \( \max_{0\le t\le T} S_t - S_T \) $$

Model assumptions:

- Black–Scholes dynamics
- Constant volatility $$ \( \sigma \) $$
- Constant risk-free rate  $$ \( r \) $$
- Continuous monitoring (exact in Barnicha engine)

---

# Engine A — Barnicha (Exact One-Step Monte Carlo)

## Core idea (no discretization bias)

Let $$ \( X_t = \log S_t \) $$. Under Black–Scholes:

$$
X_t = \log S_0 + \left(r-\frac{1}{2}\sigma^2\right)t + \sigma W_t.
$$

We simulate the terminal value exactly:

$$
X_T = \log S_0 + \left(r-\frac{1}{2}\sigma^2\right)T + \sigma \sqrt{T}\,Z,
\qquad Z\sim\mathcal N(0,1).
$$

Conditional on $$  \( (X_0, X_T) \) $$ , the process is a Brownian bridge; the conditional distribution of its running maximum is:

$$
\mathbb P\left(\max_{0\le t\le T}X_t \le m \mid X_0=x_0, X_T=x_T\right)
= 1-\exp\left(-\frac{2(m-x_0)(m-x_T)}{\sigma^2T}\right),\quad m\ge\max(x_0,x_T).
$$

Inversion gives a quadratic equation; we take the root $$\( m\ge \max(x_0,x_T) \)$$.  
The minimum is obtained by symmetry (the other root $$ \( \le \min(x_0,x_T) \)$$).

---

## Pricing estimator (per path)

Once $$\( S_T \), \( S_{\max} \), \( S_{\min} \) $$ are sampled for a path:

- **Call price estimator**
  $$
  \hat V_{\text{call}} = e^{-rT}\bigl(S_T - S_{\min}\bigr)
  $$

- **Put price estimator**
  $$
  \hat V_{\text{put}} = e^{-rT}\bigl(S_{\max} - S_T\bigr)
  $$

The final price is the Monte Carlo average of the discounted payoff estimators.

---

## Architecture (Barnicha engine)

- `Barnicha/include/lookback/` public headers
- `Barnicha/src/` implementation
- `Barnicha/app/` example executable

Main classes:

- **`OneStepBridgeSampler`**  
  Generates and stores random draws once (common random numbers), including:
  - terminal draw for $$\( X_T \)$$
  - uniforms (or equivalent) used to invert the bridge extremum distribution

- **`MonteCarloEngine`**  
  Prices using stored draws; supports variance reduction

- **`AnalyticLookback`**  
  Closed-form benchmark engine (special handling when `rate == 0` if applicable)

- **`PricingEngine`** (interface / base)  
  Common interface for pricing engines

- **`GreeksCalculator`** (abstract)  
  Interface for computing requested Greeks

- **`BumpAndRevalueGreeks`**  
  Finite-difference Greeks computed using CRN, controlled by boolean flags:
  `compute_price`, `compute_delta`, `compute_gamma`, `compute_theta`,
  `compute_rho`, `compute_vega`

---

## Variance reduction (Barnicha engine)

Antithetic variates are applied **only** to the normal draw used for the terminal value:

For each sample $$\( (Z, U_{\max}, U_{\min}) \)$$, we also evaluate  
$$\((-Z, U_{\max}, U_{\min}) \)$$and average the payoffs.

This reduces variance while preserving **exact** extremum sampling conditional on endpoints.

---

## Monte Carlo uncertainty (reported outputs)

Monte Carlo results can be reported with:

- standard error
- 95% confidence interval

Greeks also have uncertainty because they are computed **per path** using CRN.
# Engine B — Chivet (Time-Stepping Monte Carlo + BGK)

## Core idea (path simulation)

This approach simulates the asset price on a discrete grid:

$$
0=t_0 < t_1 < \dots < t_N = T,
\qquad \Delta t = \frac{T}{N}.
$$

The Black–Scholes time stepping scheme used is:

$$
S_{t_{i+1}}
=
S_{t_i}\exp\!\left(\left(r-\frac{1}{2}\sigma^2\right)\Delta t + \sigma\sqrt{\Delta t}\,Z_i\right),
\qquad Z_i\sim\mathcal N(0,1).
$$

For each simulated path, we compute:

- terminal value \( S_T \)
- discrete running extrema:
  $$
  S_{\max}^{\text{disc}}=\max_{0\le i\le N} S_{t_i},
  \qquad
  S_{\min}^{\text{disc}}=\min_{0\le i\le N} S_{t_i}.
  $$

---

## Pricing estimator (without BGK)

- **Call estimator**
  $$
  \hat V_{\text{call}} = e^{-rT}\bigl(S_T - S_{\min}^{\text{disc}}\bigr)
  $$

- **Put estimator**
  $$
  \hat V_{\text{put}} = e^{-rT}\bigl(S_{\max}^{\text{disc}} - S_T\bigr)
  $$

---

## BGK correction (optional, used for pricing)

Discrete monitoring misses extrema between grid points.  
With BGK enabled, the observed extrema are adjusted before computing the payoff:

Minimum (BGK correction):

$$
S_{\min}^{\mathrm{corr}} = S_{\min}^{\mathrm{disc}} \exp\left(-\beta \sigma \sqrt{\Delta t}\right)
$$

Maximum (BGK correction):

$$
S_{\max}^{\mathrm{corr}} = S_{\max}^{\mathrm{disc}} \exp\left(\beta \sigma \sqrt{\Delta t}\right)
$$

where the constant is:

$$
\beta \approx 0.5826.
$$

Then the payoff uses the corrected extrema:

- **Call (BGK)**
  $$
  \hat V_{\text{call}}^{\text{BGK}} = e^{-rT}\bigl(S_T - S_{\min}^{\text{corr}}\bigr)
  $$

- **Put (BGK)**
  $$
  \hat V_{\text{put}}^{\text{BGK}} = e^{-rT}\bigl(S_{\max}^{\text{corr}} - S_T\bigr)
  $$

---

# Greeks (common to both engines)

Greeks are computed with **bump-and-revalue**.

Example (central differences):

- **Delta**
  $$
  \Delta = \frac{V(S+\varepsilon)-V(S-\varepsilon)}{2\varepsilon}
  $$

- **Gamma**
  $$
  \Gamma = \frac{V(S+\varepsilon)-2V(S)+V(S-\varepsilon)}{\varepsilon^2}
  $$

- **Vega**
  $$
  \nu = \frac{V(\sigma+\varepsilon)-V(\sigma-\varepsilon)}{2\varepsilon}
  $$

Common Random Numbers (CRN) are used between bumped evaluations to reduce variance.

---

# Monte Carlo uncertainty

Reported metrics typically include:

- standard error
- 95% confidence interval

Greeks also report uncertainty since they are computed from Monte Carlo pathwise estimates.




# Build system

The project uses CMake to compile a single unified shared library.

Minimum required version:

CMake ≥ 3.16

C++ standard:

C++17


# Build instructions (Windows)

Step 1 — Create build directory

mkdir build
cd build


Step 2 — Configure project

cmake ..


Step 3 — Build shared library

cmake --build . --config Release


Output generated:

UnifiedLookback.dll


# Unified DLL architecture

The unified DLL combines both pricing engines into a single binary.

Logical structure:

Excel / External program
        ↓
unified/unified_dll.cpp
        ↓
Barnicha engine  OR  Chivet engine


The unified export layer isolates both implementations.


# Unified export layer

File:

unified/unified_dll.cpp

Role:

• provides exported functions  
• dispatches pricing calls to the correct engine  
• prevents symbol conflicts  
• ensures a stable external interface  


Typical exported functions:

PriceLookbackBarnicha(...)
PriceLookbackChivet(...)


These functions internally construct the appropriate engine objects.


# End of file
