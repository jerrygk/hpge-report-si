"""Shared analysis code for the PHAS0051 HPGe experiment.

Every session notebook imports this module, so there is one tested copy of each method.
The physics behind each function is explained in the notebooks; the docstrings here say
what each function does and which conventions it uses.

Conventions
-----------
* Channel i covers pulse heights from i - 1/2 to i + 1/2, so a fitted centroid mu is a
  channel-centre position and E(mu) = a0 + a1 * mu.
* InterSpec's polynomial gives the LOWER edge of channel i, so C0 = a0 - a1/2 and C1 = a1.
* Energies in keV, times in s, uncertainties are one standard deviation.

Sections
--------
1. Constants and paths          5. Calibration            9. Exports
2. Reading spectra              6. Resolution             10. Primordial nuclides
3. Reference data               7. Peak search
4. Peak model and fitting       8. Inventory and identification
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

# ---------------------------------------------------------------- 1. constants and paths
ROOT = Path(__file__).resolve().parents[1]          # repository root
REFERENCE = ROOT / "data" / "reference"
FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))       # 2.3548
M_E = 510.99895                                   # electron rest energy, keV (CODATA 2018)
EPS_GE = 2.96e-3                                  # mean energy per electron-hole pair in Ge at 77 K, keV
N_A = 6.02214076e23                               # Avogadro constant, mol^-1 (exact, SI 2019)
YEAR_S = 365.25 * 86400                           # Julian year, s
EARTH_AGE_Y = 4.54e9                              # age of the Earth (Dalrymple 2001), y
SQRT2 = np.sqrt(2.0)


# ---------------------------------------------------------------- 2. reading spectra
@dataclass
class Spectrum:
    counts: np.ndarray
    live: float                      # s
    real: float                      # s
    start: dt.datetime
    stored_cal: tuple                # (a0, a1, a2) written by Maestro, E = a0 + a1 i + a2 i^2
    path: Path
    sha256: str
    meta: dict = field(default_factory=dict)

    @property
    def channels(self) -> np.ndarray:
        return np.arange(len(self.counts), dtype=float)

    @property
    def end(self) -> dt.datetime:
        return self.start + dt.timedelta(seconds=self.real)

    @property
    def dead_fraction(self) -> float:
        return 1 - self.live / self.real

    def stored_energy(self, ch=None) -> np.ndarray:
        ch = self.channels if ch is None else np.asarray(ch, float)
        a0, a1, a2 = self.stored_cal
        return a0 + a1 * ch + a2 * ch**2


_MONTHS = {m: i for i, m in enumerate("JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC".split(), 1)}


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_chn(path) -> Spectrum:
    """Read an ORTEC integer .Chn file.

    Header (32 bytes, little endian): int16 -1 | uint16 MCA number | uint16 segment |
    char[2] start seconds | int32 real time | int32 live time (both in 20 ms ticks) |
    char[8] date 'DDMMMYY*' ('1' in the last place means 21st century) | char[4] time 'HHMM' |
    uint16 first channel | uint16 number of channels.
    Then one uint32 per channel, then a trailer whose first int16 is -101 or -102,
    followed by an int16 and the energy calibration coefficients as float32.
    """
    path = Path(path)
    b = path.read_bytes()
    ftype, mca, segment = struct.unpack("<hHH", b[0:6])
    if ftype != -1:
        raise ValueError(f"{path} is not an ORTEC integer .Chn file")
    seconds = b[6:8].decode("ascii")
    real_ticks, live_ticks = struct.unpack("<ii", b[8:16])
    date, hhmm = b[16:24].decode("ascii"), b[24:28].decode("ascii")
    first_channel, n_channels = struct.unpack("<HH", b[28:32])
    counts = np.frombuffer(b[32:32 + 4 * n_channels], dtype="<u4").astype(np.int64)
    year = (2000 if date[7] == "1" else 1900) + int(date[5:7])
    start = dt.datetime(year, _MONTHS[date[2:5].upper()], int(date[0:2]),
                        int(hhmm[0:2]), int(hhmm[2:4]), int(seconds))
    trailer = b[32 + 4 * n_channels:]
    trailer_type = struct.unpack("<h", trailer[0:2])[0]
    stored = tuple(float(v) for v in struct.unpack("<3f", trailer[4:16]))
    if trailer_type == -101:                      # old format: linear calibration only
        stored = (stored[0], stored[1], 0.0)
    return Spectrum(counts, live_ticks * 0.02, real_ticks * 0.02, start, stored, path, sha256(path),
                    dict(mca=mca, segment=segment, first_channel=first_channel, trailer_type=trailer_type))


def read_maestro_txt(path) -> np.ndarray:
    """Read the Maestro printout: lines of the form 'channel:  c c c c c c c'. Counts only."""
    counts = []
    with open(path, encoding="ascii", errors="replace") as f:
        for line in f:
            head, sep, rest = line.partition(":")
            if sep and head.strip().isdigit():
                counts += [int(v) for v in rest.split()]
    return np.array(counts, dtype=np.int64)


# ---------------------------------------------------------------- 3. reference data
def pretty(nuc: str) -> str:
    """'214Bi' -> 'Bi-214', '234mPa' -> 'Pa-234m'."""
    A, m, sym = re.fullmatch(r"(\d+)(m?)([A-Z][a-z]?)", nuc).groups()
    return f"{sym}-{A}{m}"


def load_ensdf(path=None) -> pd.DataFrame:
    df = pd.read_csv(path or REFERENCE / "gamma_lines_ENSDF.csv")
    df["source_name"] = df.nuclide.map(pretty)
    return df


def load_features(path=None) -> pd.DataFrame:
    return pd.read_csv(path or REFERENCE / "non_nuclear_features.csv")


def load_primordial(path=None) -> pd.DataFrame:
    return pd.read_csv(path or REFERENCE / "primordial_nuclides.csv")


def build_library(ensdf=None, features=None) -> list:
    """List of (energy, source, origin): every ENSDF line plus the non-nuclear features."""
    ensdf = load_ensdf() if ensdf is None else ensdf
    features = load_features() if features is None else features
    lib = [(r.E_keV, r.source_name, r.chain) for r in ensdf.itertuples()]
    lib += [(r.E_keV, r.label, r.origin) for r in features.itertuples()]
    return lib


# ---------------------------------------------------------------- 4. peak model and fitting
BASE = ["b0", "b1", "sigma", "step", "eta", "beta"]   # shared parameters, then (A_k, mu_k) per peak


def gauss_in_bins(x, mu, sigma):
    """Fraction of a unit-area Gaussian falling in channel x (channel x spans x-1/2 to x+1/2)."""
    return 0.5 * (special.erf((x + 0.5 - mu) / (SQRT2 * sigma)) - special.erf((x - 0.5 - mu) / (SQRT2 * sigma)))


def step_shape(x, mu, sigma):
    """Smoothed step: 1 well below the peak, 0 well above it."""
    return 0.5 * special.erfc((x - mu) / (SQRT2 * sigma))


def tail_shape(x, mu, sigma, beta):
    """Unit-area Gaussian convolved with an exponential tail on the low-energy side (decay length beta).
    exp(a) * erfc(z) is evaluated as erfcx(z) * exp(-u^2 / 2 sigma^2) when z > 0 to avoid overflow."""
    u = x - mu
    z = u / (SQRT2 * sigma) + sigma / (SQRT2 * beta)
    with np.errstate(over="ignore", invalid="ignore"):
        stable = special.erfcx(np.clip(z, 0, None)) * np.exp(-u**2 / (2 * sigma**2))
        direct = np.exp(np.minimum(u / beta + sigma**2 / (2 * beta**2), 700)) * special.erfc(z)
    return np.where(z > 0, stable, direct) / (2 * beta)


def roi_model(theta, x, n_peaks, x0):
    """Linear continuum + n_peaks Gaussians (shared sigma) + smoothed steps (+ optional tails)."""
    b0, b1, sigma, step, eta, beta = theta[:6]
    f = b0 + b1 * (x - x0)
    for k in range(n_peaks):
        A, mu = theta[6 + 2 * k], theta[7 + 2 * k]
        shape = (1 - eta) * gauss_in_bins(x, mu, sigma)
        if eta > 0:
            shape = shape + eta * tail_shape(x, mu, sigma, beta)
        f = f + A * shape + step * A / (np.sqrt(2 * np.pi) * sigma) * step_shape(x, mu, sigma)
    return f


def fit_roi(y_all, centres, sigma0, half=7.0, step=None, sigma_free=True, tail=False):
    """Poisson maximum-likelihood fit of one ROI.

    centres    : starting channel(s) of the peak(s) in the ROI
    sigma0     : starting (or fixed) Gaussian sigma in channels
    half       : ROI half-width in units of sigma0, measured from the outermost peaks
    step       : None to fit the step ratio, or a number to fix it
    sigma_free : fit sigma (True) or fix it at sigma0 (False)
    tail       : add a low-energy exponential tail (fraction eta, decay length beta >= sigma0)
    Returns a dict with parameters, uncertainties (inverse expected Fisher information),
    the Poisson deviance and its p-value.
    """
    centres = np.atleast_1d(np.asarray(centres, float))
    n = len(centres)
    lo = max(int(np.floor(centres.min() - half * sigma0)), 0)
    hi = min(int(np.ceil(centres.max() + half * sigma0)), len(y_all) - 1)
    x = np.arange(lo, hi + 1, dtype=float)
    y = np.asarray(y_all[lo:hi + 1], float)
    x0 = 0.5 * (lo + hi)

    m = max(3, int(sigma0))
    left, right = y[:m].mean(), y[-m:].mean()
    theta0 = [0.5 * (left + right), (right - left) / (x[-1] - x[0]), sigma0,
              0.0 if step is None else step, 0.1 if tail else 0.0, 2.0 * sigma0]
    for c in centres:
        near = np.abs(x - c) <= 2.5 * sigma0
        theta0 += [max((y[near] - 0.5 * (left + right)).sum(), 10.0), c]
    theta0 = np.array(theta0, float)

    free = np.ones(len(theta0), bool)
    free[2] = sigma_free
    free[3] = step is None
    free[4] = free[5] = tail
    lb = np.full(len(theta0), -np.inf)
    ub = np.full(len(theta0), np.inf)
    lb[2], ub[2] = 0.3 * sigma0, 3.0 * sigma0
    lb[3], ub[3] = -0.2, 0.5
    lb[4], ub[4] = 0.0, 0.6
    lb[5], ub[5] = 1.0 * sigma0, 15.0 * sigma0
    for k, c in enumerate(centres):
        lb[6 + 2 * k] = 0.0
        lb[7 + 2 * k], ub[7 + 2 * k] = c - 3 * sigma0, c + 3 * sigma0

    def full(p):
        t = theta0.copy()
        t[free] = p
        return t

    def deviance_residuals(p):
        f = np.clip(roi_model(full(p), x, n, x0), 1e-9, None)
        d = 2 * (f - y + np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / f), 0.0))
        return np.sign(y - f) * np.sqrt(np.clip(d, 0, None))

    p0 = np.clip(theta0[free], lb[free] + 1e-9, ub[free] - 1e-9)
    sol = optimize.least_squares(deviance_residuals, p0, bounds=(lb[free], ub[free]), x_scale="jac", max_nfev=20000)
    theta = full(sol.x)

    f = np.clip(roi_model(theta, x, n, x0), 1e-9, None)
    idx = np.flatnonzero(free)
    J = np.empty((len(x), len(idx)))
    for j, i in enumerate(idx):
        h = 1e-6 * max(abs(theta[i]), 1e-3)
        tp, tm = theta.copy(), theta.copy()
        tp[i] += h
        tm[i] -= h
        J[:, j] = (roi_model(tp, x, n, x0) - roi_model(tm, x, n, x0)) / (2 * h)
    cov = np.zeros((len(theta), len(theta)))
    cov[np.ix_(idx, idx)] = np.linalg.inv(J.T @ (J / f[:, None]))
    dev = float(np.sum(sol.fun**2))
    ndf = len(x) - len(idx)
    names = BASE + [f"{p}{k + 1}" for k in range(n) for p in ("A", "mu")]
    return dict(theta=theta, err=np.sqrt(np.diag(cov)), cov=cov, names=names, x=x, y=y, model=f, x0=x0,
                n_peaks=n, dev=dev, ndf=ndf, pval=float(stats.chi2.sf(dev, ndf)), ok=sol.success)


def peak(fit, k=0) -> dict:
    """Area, centroid and sigma (with uncertainties) of peak k in a fit."""
    t, e = fit["theta"], fit["err"]
    return dict(A=t[6 + 2 * k], sA=e[6 + 2 * k], mu=t[7 + 2 * k], smu=e[7 + 2 * k], sigma=t[2], ssigma=e[2])


def plot_fit(fit, ax, axr, title=""):
    """Data with Poisson error bars, fitted model and normalised residuals."""
    x, y, f = fit["x"], fit["y"], fit["model"]
    ax.errorbar(x, y, yerr=np.sqrt(np.maximum(y, 1)), fmt=".", ms=3, color="k", lw=0.6)
    ax.step(x, f, where="mid", color="tab:red", lw=1.2)
    ax.set_title(title, fontsize=9)
    axr.bar(x, (y - f) / np.sqrt(f), width=1, color="tab:gray")
    axr.axhline(0, color="k", lw=0.6)
    axr.set_ylim(-4, 4)
    axr.set_yticks([-3, 0, 3])
    axr.tick_params(labelsize=7)


def run_toys(truth, n_toys, rng, x_len=140, step=0.025) -> pd.DataFrame:
    """Simulate spectra from the model with known parameters, refit each one and return the pulls."""
    x = np.arange(x_len, dtype=float)
    theta = np.array([truth["b0"], 0.0, truth["sigma"], step, 0.0, 1.0, truth["A"], truth["mu"]])
    mean = roi_model(theta, x, 1, x.mean())
    rows = []
    for _ in range(n_toys):
        y = rng.poisson(mean)
        fit = fit_roi(y, [truth["mu"] + rng.uniform(-1, 1)], truth["sigma"] * rng.uniform(0.9, 1.1), step=step)
        p = peak(fit)
        rows.append(((p["mu"] - truth["mu"]) / p["smu"], (p["A"] - truth["A"]) / p["sA"],
                     (p["sigma"] - truth["sigma"]) / p["ssigma"], fit["dev"], fit["ndf"]))
    return pd.DataFrame(rows, columns=["pull_mu", "pull_A", "pull_sigma", "dev", "ndf"])


# ---------------------------------------------------------------- 5. calibration
def highest_channel(counts, guess, window) -> int:
    g = int(round(guess))
    return g - window + int(np.argmax(counts[g - window:g + window + 1]))


def bootstrap_map(spec: Spectrum, E1=1460.820, E2=2614.511, window=30):
    """Rough linear map from two strong isolated lines found on the stored energy scale.
    Returns (g0, g1, c1, c2) with E = g0 + g1 * channel."""
    E_st = spec.stored_energy()
    c1 = highest_channel(spec.counts, np.interp(E1, E_st, spec.channels), window)
    c2 = highest_channel(spec.counts, np.interp(E2, E_st, spec.channels), window)
    g1 = (E2 - E1) / (c2 - c1)
    return E1 - g1 * c1, g1, c1, c2


@dataclass
class Calibration:
    """Linear (or polynomial) energy calibration E(mu) = sum a_j mu^j with covariance C."""
    a: np.ndarray
    C: np.ndarray
    chi2: float
    ndf: int
    scale: float = 1.0               # sqrt(chi2/ndf) already applied to C when > 1

    def energy(self, mu):
        mu = np.asarray(mu, float)
        return sum(c * mu**j for j, c in enumerate(self.a))

    def sigma(self, mu):
        mu = np.asarray(mu, float)
        J = np.stack([mu**j for j in range(len(self.a))], axis=-1)
        return np.sqrt(np.einsum("...i,ij,...j->...", J, self.C, J))

    def channel(self, E):
        if len(self.a) != 2:
            raise NotImplementedError("inverse only implemented for a linear calibration")
        return (np.asarray(E, float) - self.a[0]) / self.a[1]

    @property
    def interspec(self):
        """(C0, C1) for InterSpec, whose polynomial gives the lower edge of each channel."""
        return self.a[0] - self.a[1] / 2, self.a[1]


def calibrate(mu, smu, E, sE, degree=1, n_iter=3, scale_by_chi2=True) -> Calibration:
    """Weighted polynomial fit E(mu) with effective variance s^2 = sE^2 + (dE/dmu * smu)^2.
    If scale_by_chi2, the covariance is multiplied by S^2 = chi2/ndf when S > 1 (PDG convention)."""
    mu, smu, E, sE = (np.asarray(v, float) for v in (mu, smu, E, sE))
    A = np.vander(mu, degree + 1, increasing=True)
    slope = np.polyfit(mu, E, 1)[0] * np.ones_like(mu)
    for _ in range(n_iter):
        s = np.sqrt(sE**2 + (slope * smu)**2)
        W = 1 / s**2
        C = np.linalg.inv(A.T @ (A * W[:, None]))
        a = C @ A.T @ (W * E)
        slope = sum(j * a[j] * mu**(j - 1) for j in range(1, degree + 1))
    chi2 = float(np.sum(((E - A @ a) / s)**2))
    ndf = len(E) - degree - 1
    S = max(1.0, np.sqrt(chi2 / ndf)) if (scale_by_chi2 and ndf > 0) else 1.0
    return Calibration(a=a, C=C * S**2, chi2=chi2, ndf=ndf, scale=S)


def f_test(lin: Calibration, quad: Calibration):
    """F statistic and p-value for adding one polynomial term."""
    F = (lin.chi2 - quad.chi2) / (quad.chi2 / quad.ndf)
    return F, float(stats.f.sf(F, 1, quad.ndf))


# ---------------------------------------------------------------- 6. resolution
def fwhm3(E, w0, w1, w2):
    """FWHM^2 = w0^2 + w1^2 E + w2^2 E^2 (noise, carrier statistics, charge collection and drift)."""
    E = np.clip(np.asarray(E, float), 0, None)
    return np.sqrt(w0**2 + w1**2 * E + w2**2 * E**2)


def fano_from_w1(w1, sw1=0.0):
    """Fano factor from w1^2 = (2 sqrt(2 ln 2))^2 F eps."""
    k = FWHM_PER_SIGMA**2 * EPS_GE
    return w1**2 / k, 2 * w1 * sw1 / k


# ---------------------------------------------------------------- 7. peak search
def peak_search(counts, sigma_ch, threshold=5.0, sigmas=np.arange(2.0, 7.01, 0.25)):
    """Matched second-derivative (Mexican-hat) filter search.
    sigma_ch : expected Gaussian sigma (channels) for every channel.
    Returns (candidate channels, significance array S)."""
    y = np.asarray(counts, float)
    sig = np.zeros((len(sigmas), len(y)))
    for i, s in enumerate(sigmas):
        j = np.arange(-int(4 * s), int(4 * s) + 1)
        k = (1 - j**2 / s**2) * np.exp(-j**2 / (2 * s**2))
        k -= k.mean()                                     # zero sum: flat and sloping continua give no response
        out = np.convolve(y, k[::-1], mode="same")
        var = np.convolve(y, (k**2)[::-1], mode="same")
        sig[i] = out / np.sqrt(np.clip(var, 1, None))
    pick = np.clip(np.round((sigma_ch - sigmas[0]) / (sigmas[1] - sigmas[0])).astype(int), 0, len(sigmas) - 1)
    S = sig[pick, np.arange(len(y))]
    first = np.flatnonzero(y)[0] + 10                    # skip the ADC threshold edge
    cands = [i for i in range(first, len(y) - 10)
             if S[i] > threshold and S[i] == S[max(i - int(1.5 * sigma_ch[i]), 0):i + int(1.5 * sigma_ch[i]) + 1].max()]
    return cands, S


# ---------------------------------------------------------------- 8. inventory and identification
def group_candidates(cands, sigma_ch, n_sigma=10.0) -> list:
    """Group candidates closer than n_sigma: a closer neighbour would sit inside the +-7 sigma ROI."""
    groups, current = [], [cands[0]]
    for i in cands[1:]:
        if i - current[-1] < n_sigma * sigma_ch[i]:
            current.append(i)
        else:
            groups.append(current)
            current = [i]
    groups.append(current)
    return groups


def fit_inventory(counts, cands, S, cal: Calibration, sigma_model_ch, step, library, live) -> pd.DataFrame:
    """Fit every candidate (sigma free only for strong isolated peaks) and match it to the library."""
    sig_all = sigma_model_ch(cal.energy(np.arange(len(counts))))
    lib_E = np.array([r[0] for r in library])
    rows = []
    for g in group_candidates(cands, sig_all):
        s_model = float(sigma_model_ch(cal.energy(g[0])))
        free_width = len(g) == 1 and S[g[0]] > 10
        fit = fit_roi(counts, g, s_model, step=step, sigma_free=free_width)
        for k, i in enumerate(g):
            p = peak(fit, k)
            E = float(cal.energy(p["mu"]))
            sE_stat = cal.a[1] * p["smu"]
            sE = float(np.hypot(sE_stat, cal.sigma(p["mu"])))
            near = np.flatnonzero(np.abs(lib_E - E) < max(0.5, 3 * sE))
            near = near[np.argsort(np.abs(lib_E[near] - E))]
            best = library[near[0]] if len(near) else None
            rows.append(dict(E_keV=E, sE_keV=sE, sE_stat=sE_stat, mu=p["mu"], smu=p["smu"], S=S[i],
                             area=p["A"], s_area=p["sA"],
                             width_ratio=p["sigma"] / s_model if free_width else np.nan,
                             dev_ndf=fit["dev"] / fit["ndf"], p_fit=fit["pval"],
                             source=best[1] if best else "unidentified", origin=best[2] if best else "",
                             E_lib=best[0] if best else np.nan, dE=E - best[0] if best else np.nan,
                             other_candidates=", ".join(f"{library[j][1]} {library[j][0]:.1f}" for j in near[1:])))
    df = pd.DataFrame(rows)
    df["rate_cpd"] = df.area / live * 86400
    df["s_rate_cpd"] = df.s_area / live * 86400
    return df


# ---------------------------------------------------------------- 9. exports
def write_calibrated_spectrum(path, spec: Spectrum, cal: Calibration, notes=()):
    """CSV with channel, energy at lower edge / centre / upper edge, calibration uncertainty,
    counts and count rate. Metadata goes in '#' comment lines at the top."""
    ch = spec.channels
    df = pd.DataFrame({"channel": ch.astype(int),
                       "E_low_keV": cal.energy(ch - 0.5), "E_centre_keV": cal.energy(ch),
                       "E_high_keV": cal.energy(ch + 0.5), "sigma_E_cal_keV": cal.sigma(ch),
                       "counts": spec.counts, "rate_per_s": spec.counts / spec.live})
    head = [f"# source file: {spec.path.name}  (sha256 {spec.sha256})",
            f"# start: {spec.start.isoformat()}  live time: {spec.live:.2f} s  real time: {spec.real:.2f} s",
            f"# calibration: E = {cal.a[0]:.6f} + {cal.a[1]:.10f} * channel  (channel centre)",
            f"# rate_per_s uses the live time"] + [f"# {n}" for n in notes]
    with open(path, "w") as f:
        f.write("\n".join(head) + "\n")
        df.to_csv(f, index=False, float_format="%.6f")
    return df


def calibration_record(spec: Spectrum, cal: Calibration, lines: pd.DataFrame, extra=None) -> dict:
    C0, C1 = cal.interspec
    rec = dict(
        created=dt.datetime.now().isoformat(timespec="seconds"),
        source_file=spec.path.name, source_sha256=spec.sha256,
        start=spec.start.isoformat(), live_time_s=spec.live, real_time_s=spec.real,
        model="E = a0 + a1 * mu, mu = channel-centre position (channel i spans i-1/2 to i+1/2)",
        a0_keV=float(cal.a[0]), a1_keV_per_channel=float(cal.a[1]),
        covariance=cal.C.tolist(), sigma_a0=float(np.sqrt(cal.C[0, 0])), sigma_a1=float(np.sqrt(cal.C[1, 1])),
        correlation=float(cal.C[0, 1] / np.sqrt(cal.C[0, 0] * cal.C[1, 1])),
        chi2=cal.chi2, ndf=cal.ndf, scale_factor_applied=cal.scale,
        interspec_lower_edge=dict(C0_keV=float(C0), C1_keV_per_channel=float(C1)),
        interspec_if_channel_centre=dict(C0_keV=float(cal.a[0]), C1_keV_per_channel=float(C1)),
        stored_calibration_in_file=list(spec.stored_cal),
        lines_used=lines[["nuclide", "E_ref", "sE_ref", "mu", "smu"]].to_dict(orient="records"))
    if extra:
        rec.update(extra)
    return rec


def write_calp(path, C0, C1, C2=0.0):
    """InterSpec / PeakEasy CALp file (polynomial, lower channel edge). Labels taken from the
    strings InterSpec itself writes ('#PeakEasy CALp File Ver:  4.00' ... '#END'). InterSpec writes
    6 significant figures, which would shift 2.6 MeV by about 0.006 keV, so 10 are written here."""
    lines = ["#PeakEasy CALp File Ver:  4.00",
             f"Offset (keV)           :  {C0:1.9e}",
             f"Gain (keV / Chan)      :  {C1:1.9e}",
             f"2nd Order Coef         :  {C2:1.9e}",
             "#END"]
    Path(path).write_text("\n".join(lines) + "\n")


def read_calp(path) -> dict:
    out = {}
    for line in Path(path).read_text().splitlines():
        if ":" in line and not line.startswith("#"):
            key, val = line.split(":", 1)
            out[key.strip()] = float(val)
    return out


# ---------------------------------------------------------------- 10. primordial nuclides
def primordial_table(df=None, age_y=EARTH_AGE_Y) -> pd.DataFrame:
    """Survival since the Earth formed, activity per kg of average upper crust and photon emission."""
    df = load_primordial() if df is None else df.copy()
    lam = np.log(2) / (df.half_life_y * YEAR_S)
    atoms_per_kg = df.ucc_ppm * 1e-6 * 1000 / df.molar_mass_g_mol * N_A * df.isotopic_abundance
    df["fraction_surviving"] = 2.0 ** (-age_y / df.half_life_y)
    df["crust_Bq_kg"] = lam * atoms_per_kg
    df["crust_photons_s_kg"] = df.crust_Bq_kg * df.gamma_yield
    df["Bq_per_g_element"] = lam * N_A / df.molar_mass_g_mol * df.isotopic_abundance
    return df


def decay_factor(t_years, half_life_y):
    return 2.0 ** (-np.asarray(t_years, float) / half_life_y)


def half_life_years(text) -> float:
    """'30.08 y', '53.22 d', '6.15 h', '3.053 m', '3.96 s' -> years (NaN if empty)."""
    if not isinstance(text, str) or not text.strip():
        return np.nan
    value, unit = text.split()
    scale = {"y": 1.0, "d": 1 / 365.25, "h": 1 / 8766.0, "m": 1 / 525960.0, "s": 1 / YEAR_S}[unit]
    return float(value) * scale


def load_calibration(path) -> Calibration:
    """Rebuild a Calibration from a JSON record written by calibration_record(), e.g. to check gain
    drift in a later session or to apply the background calibration to another spectrum."""
    rec = json.loads(Path(path).read_text())
    return Calibration(a=np.array([rec["a0_keV"], rec["a1_keV_per_channel"]]), C=np.array(rec["covariance"]),
                       chi2=rec["chi2"], ndf=rec["ndf"], scale=rec["scale_factor_applied"])


def read_interspec_csv(path) -> pd.DataFrame:
    """Read an InterSpec Peak Manager CSV (two header rows: column names, then units)."""
    import csv
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        rows = list(csv.reader(f))
    h = next(i for i, r in enumerate(rows) if r and r[0].strip() == "Centroid")
    names = [c.strip() for c in rows[h]]
    units = [c.strip() for c in rows[h + 1]] + [""] * len(names)
    cols = [f"{n}|{u}" if u else n for n, u in zip(names, units)]
    data = [r for r in rows[h + 2:] if r and r[0].strip()]
    df = pd.DataFrame([[c.strip() for c in r[:len(cols)]] for r in data], columns=cols[:len(data[0])])
    num = lambda c: pd.to_numeric(df[c], errors="coerce") if c in df else np.nan
    return pd.DataFrame({"E_interspec": num("Centroid|keV"), "area_interspec": num("Net_Area|Counts"),
                         "s_area_interspec": num("Net_Area|Uncertainty"), "fwhm_interspec": num("FWHM|keV"),
                         "nuclide": df.get("Nuclide", ""), "E_photopeak": num("Photopeak_Energy|keV")})
