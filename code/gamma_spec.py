"""gamma_spec: a general toolkit for HPGe (and other) gamma-ray spectra.

Written for PHAS0051 session 2 (27 Oct 2026). It replaces nothing in the session-1 module
`code/hpge_analysis.py`, which is left untouched so that session-1 results stay
reproducible. This module is meant to work for any similar experiment, so nothing in it assumes a
particular detector, MCA, file format, gain or source.

Why a new module was needed (each point is demonstrated in the session-2 notebook):
  1. The session-2 spectra come from a Mirion/Canberra Lynx MCA as CSV files. The old code reads
     only ORTEC .Chn and Maestro .Txt files.
  2. The Lynx files carry no energy calibration (offset 0, slope 1, so the "Energy (keV)" column is
     simply the channel number). The old code started its calibration from the calibration stored
     in the file, so it looked for K-40 at channel ~1461 and found nothing useful.
  3. A different detector has a different gain, so the session-1 calibration cannot be transferred.
     The new code finds the gain from the spectrum itself by pattern matching known lines.
  4. Photopeaks have a low-energy tail (incomplete charge collection). The old inventory fitted
     symmetric Gaussians. The new peak model is an exponentially modified Gaussian (EMG), the model
     that won the session-1 K-40 study, integrated exactly over each channel, and the peak position
     is its mode (the most probable pulse height).
  5. Sample work needs background subtraction, comparisons between spectra and exports that
     InterSpec can read with the calibration already applied (IAEA .Spe).

Conventions (the same as session 1, so results can be compared):
  * Channels are indexed from 0 in arrays. Index i spans the continuous pulse-height coordinate
    i - 1/2 to i + 1/2, so a fitted position is a "channel-centre" coordinate and
    E(x) = a0 + a1 x (+ a2 x^2).
  * Files that number their channels from 1 (the Lynx CSV) are shifted on reading: file channel n
    becomes index n - 1. The shift is recorded in Spectrum.meta["first_channel_label"].
  * InterSpec's polynomial gives the LOWER edge of channel i, so its coefficients are obtained by
    evaluating our polynomial at x = i - 1/2 (see Calibration.lower_edge_coefficients).
  * Energies in keV, times in s, rates per second of live time unless stated otherwise,
    uncertainties are one standard deviation.

Sections
  1. constants                     6. automatic first calibration (pattern matching)
  2. reading spectra               7. energy calibration fit
  3. reference data                8. resolution and tail models
  4. peak shapes                   9. peak search, inventory and identification
  5. fitting a region of interest  10. comparing spectra and background subtraction
                                   11. exports (CSV, JSON, CALp, IAEA .Spe) and InterSpec CSV
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

# ============================================================== 1. constants
FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))     # 2.3548, FWHM of a Gaussian in units of sigma
M_E = 510.99895                                       # electron rest energy, keV (CODATA 2018)
EPS_GE = 2.96e-3                                      # mean energy per electron-hole pair in Ge at 77 K, keV
SQRT2 = np.sqrt(2.0)
SQRT2PI = np.sqrt(2.0 * np.pi)

# The reference tables live in data/reference of this repository. A different experiment can
# pass its own file paths to the loaders instead.
ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "reference"


# ============================================================== 2. reading spectra
@dataclass
class Spectrum:
    """One measured spectrum and the metadata needed to turn counts into rates."""
    counts: np.ndarray               # counts per channel, index 0 = first channel in the file
    live: float                      # live time, s (time the electronics could accept a pulse)
    real: float                      # real (clock) time, s
    start: dt.datetime | None        # start of acquisition as written in the file (may be wrong)
    stored_cal: tuple                # (a0, a1, a2) found in the file, E = a0 + a1 n + a2 n^2 in FILE numbering
    path: Path
    sha256: str
    fmt: str                         # "lynx-csv", "ortec-chn", "iaea-spe", "maestro-txt", "columns"
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.counts)

    @property
    def channels(self) -> np.ndarray:
        return np.arange(len(self.counts), dtype=float)

    @property
    def end(self):
        return None if self.start is None else self.start + dt.timedelta(seconds=self.real)

    @property
    def dead_fraction(self) -> float:
        return 1.0 - self.live / self.real

    @property
    def label(self) -> str:
        return self.meta.get("label", self.path.stem)


def sha256(path) -> str:
    """Checksum of the raw file, stored with every result so that any change to the data shows up."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


_MONTHS = {m: i for i, m in enumerate("JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC".split(), 1)}


def parse_js_date(text):
    """Parse a JavaScript Date.toString() stamp such as 'Mon Oct 19 2026 11:33:35 GMT+0100 (BST)'.
    The Lynx web interface writes its start time in this form. Returns a timezone-aware datetime."""
    m = re.search(r"([A-Za-z]{3})\s+(\d{1,2})\s+(\d{4})\s+(\d{2}):(\d{2}):(\d{2})\s*GMT([+-])(\d{2})(\d{2})", text)
    if not m:
        return None
    mon, day, year, hh, mm, ss, sign, oh, om = m.groups()
    offset = (1 if sign == "+" else -1) * dt.timedelta(hours=int(oh), minutes=int(om))
    return dt.datetime(int(year), _MONTHS[mon.upper()], int(day), int(hh), int(mm), int(ss),
                       tzinfo=dt.timezone(offset))


def _lines(path):
    """All lines of a text file, whatever its line endings (\\n, \\r\\n or a bare \\r)."""
    text = Path(path).read_bytes().decode("utf-8", errors="replace")
    return text.replace("\r\n", "\n").replace("\r", "\n").split("\n")


def read_lynx_csv(path) -> Spectrum:
    """Read the CSV written by the Lynx web interface.

    Layout: a few 'key, value' header lines (Start Time, Energy calibration, Live Time (s),
    Real Time (s), Elapsed Computational), a line 'Spectrum', a column header
    'Channel, Energy (keV), Counts', then one row per channel. Channels are numbered from 1.
    Files re-saved by a spreadsheet may carry trailing commas, which are ignored."""
    path = Path(path)
    header, rows = {}, []
    for line in _lines(path):
        parts = [p.strip() for p in line.split(",")]
        while parts and parts[-1] == "":
            parts.pop()
        if not parts:
            continue
        if parts[0].isdigit() and len(parts) >= 3:
            rows.append((int(parts[0]), float(parts[1]), int(float(parts[2]))))
        elif len(parts) >= 2 and not rows:
            header[parts[0]] = ",".join(parts[1:])
    if not rows:
        raise ValueError(f"{path} has no channel rows")
    ch = np.array([r[0] for r in rows])
    e_col = np.array([r[1] for r in rows])
    counts = np.array([r[2] for r in rows], dtype=np.int64)
    if not np.array_equal(ch, np.arange(ch[0], ch[0] + len(ch))):
        raise ValueError(f"{path}: channel numbers are not consecutive")
    cal = header.get("Energy calibration", "")
    num = lambda key: float(re.search(key + r"\s*:\s*([-+0-9.eE]+)", cal).group(1)) if key in cal else 0.0
    stored = (num("Offset"), num("Slope"), num("Quadratic"))
    live = float(header["Live Time (s)"])
    real = float(header["Real Time (s)"])
    start = parse_js_date(header.get("Start Time", ""))
    meta = dict(first_channel_label=int(ch[0]), header=header,
                energy_column_is_channel_number=bool(np.allclose(e_col, ch)))
    return Spectrum(counts, live, real, start, stored, path, sha256(path), "lynx-csv", meta)


def read_chn(path) -> Spectrum:
    """Read an ORTEC integer .Chn file (header 32 bytes, uint32 counts, trailer with calibration).
    Real and live times are stored in 20 ms ticks; the trailer holds float32 calibration coefficients."""
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
    start = dt.datetime(year, _MONTHS[date[2:5].upper()], int(date[0:2]), int(hhmm[0:2]), int(hhmm[2:4]), int(seconds))
    trailer = b[32 + 4 * n_channels:]
    trailer_type = struct.unpack("<h", trailer[0:2])[0]
    stored = tuple(float(v) for v in struct.unpack("<3f", trailer[4:16]))
    if trailer_type == -101:
        stored = (stored[0], stored[1], 0.0)
    meta = dict(first_channel_label=int(first_channel), mca=mca, segment=segment)
    return Spectrum(counts, live_ticks * 0.02, real_ticks * 0.02, start, stored, path, sha256(path), "ortec-chn", meta)


def read_spe(path) -> Spectrum:
    """Read an IAEA/ORTEC ASCII .Spe file ($DATA, $MEAS_TIM, $DATE_MEA, $MCA_CAL or $ENER_FIT)."""
    path = Path(path)
    lines = [l.strip() for l in _lines(path)]
    blocks, key = {}, None
    for l in lines:
        if l.startswith("$"):
            key = l.rstrip(":")
            blocks[key] = []
        elif key:
            blocks[key].append(l)
    lo, hi = (int(v) for v in blocks["$DATA"][0].split()[:2])
    counts = np.array([int(float(v)) for row in blocks["$DATA"][1:] for v in row.split()][:hi - lo + 1], dtype=np.int64)
    live, real = (float(v) for v in blocks["$MEAS_TIM"][0].split()[:2])
    start = None
    if "$DATE_MEA" in blocks and blocks["$DATE_MEA"]:
        try:
            start = dt.datetime.strptime(blocks["$DATE_MEA"][0], "%m/%d/%Y %H:%M:%S")
        except ValueError:
            pass
    stored = (0.0, 1.0, 0.0)
    if "$MCA_CAL" in blocks:
        c = [float(v) for v in blocks["$MCA_CAL"][1].split()[:3]]
        stored = tuple(c + [0.0] * (3 - len(c)))
    elif "$ENER_FIT" in blocks:
        c = [float(v) for v in blocks["$ENER_FIT"][0].split()[:2]]
        stored = (c[0], c[1], 0.0)
    return Spectrum(counts, live, real, start, stored, path, sha256(path), "iaea-spe", dict(first_channel_label=lo))


def read_maestro_txt(path) -> np.ndarray:
    """Counts from a Maestro printout ('channel:  c c c c c c c'). It holds no times, so counts only."""
    counts = []
    for line in _lines(path):
        head, sep, rest = line.partition(":")
        if sep and head.strip().isdigit():
            counts += [int(v) for v in rest.split()]
    return np.array(counts, dtype=np.int64)


def read_spectrum(path, live=None, real=None, label=None) -> Spectrum:
    """Open any supported spectrum file, choosing the reader from the extension and the content.

    For a bare two- or three-column text file (channel[, energy], counts) the live and real times
    must be supplied, because without them counts cannot be turned into rates."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".chn":
        spec = read_chn(path)
    elif ext == ".spe":
        spec = read_spe(path)
    else:
        head = "\n".join(_lines(path)[:12])
        if "Live Time" in head and "Channel" in head:
            spec = read_lynx_csv(path)
        elif re.search(r"^\s*\d+:\s+\d+", head, re.M):
            if live is None or real is None:
                raise ValueError("a Maestro printout has no times; pass live= and real=")
            spec = Spectrum(read_maestro_txt(path), live, real, None, (0.0, 1.0, 0.0), path, sha256(path),
                            "maestro-txt", dict(first_channel_label=0))
        else:
            if live is None or real is None:
                raise ValueError("a plain column file has no times; pass live= and real=")
            data = np.loadtxt(path, delimiter="," if "," in head else None, comments="#")
            spec = Spectrum(data[:, -1].astype(np.int64), live, real, None, (0.0, 1.0, 0.0), path, sha256(path),
                            "columns", dict(first_channel_label=int(data[0, 0])))
    if label:
        spec.meta["label"] = label
    return spec


def describe(spectra) -> pd.DataFrame:
    """One row per spectrum: times, dead time, total counts and rate. Used to check the metadata."""
    rows = []
    for s in spectra:
        rows.append(dict(label=s.label, file=s.path.name, format=s.fmt, channels=s.n,
                         start=s.start.strftime("%a %d %b %Y %H:%M:%S") if s.start else "",
                         end_from_real_time=s.end.strftime("%a %d %b %H:%M") if s.start else "",
                         live_d=s.live / 86400, real_d=s.real / 86400, dead_pct=100 * s.dead_fraction,
                         total_counts=int(s.counts.sum()), rate_cps=s.counts.sum() / s.live,
                         stored_cal=s.stored_cal))
    return pd.DataFrame(rows)


def download_time(path):
    """When the file was downloaded, from the macOS 'where from' metadata the browser attaches
    (kMDItemDownloadedDate). An independent clock for checking the start time written in the file.
    Returns None on other systems or for files that were not downloaded."""
    import subprocess
    try:
        out = subprocess.run(["mdls", "-raw", "-name", "kMDItemDownloadedDate", str(path)],
                             capture_output=True, text=True, timeout=10).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r"(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) \+0000", out)
    if not m:
        return None
    return dt.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S").replace(tzinfo=dt.timezone.utc)


def timeline_check(spectra) -> list:
    """Look for impossible timing: two runs on one MCA cannot overlap or share a start time.
    Returns a list of warnings in plain words (empty if everything is consistent)."""
    warn = []
    ss = [s for s in spectra if s.start is not None]
    for i, a in enumerate(ss):
        for b in ss[i + 1:]:
            if a.start == b.start:
                warn.append(f"{a.label} and {b.label} have the same start time {a.start}; one header is wrong")
            elif a.start < b.start < a.end or b.start < a.start < b.end:
                warn.append(f"{a.label} ({a.start:%d %b %H:%M} to {a.end:%d %b %H:%M}) overlaps {b.label} "
                            f"({b.start:%d %b %H:%M} to {b.end:%d %b %H:%M})")
    return warn


# ============================================================== 3. reference data
def pretty(nuc: str) -> str:
    """'214Bi' -> 'Bi-214', '234mPa' -> 'Pa-234m'."""
    m = re.fullmatch(r"(\d+)(m?)([A-Z][a-z]?)", nuc)
    if not m:
        return nuc
    A, iso, sym = m.groups()
    return f"{sym}-{A}{iso}"


def load_ensdf(path=None) -> pd.DataFrame:
    df = pd.read_csv(path or REFERENCE / "gamma_lines_ENSDF.csv").drop_duplicates(["nuclide", "E_keV"])
    df["source_name"] = df.nuclide.map(pretty)
    return df.reset_index(drop=True)


def load_features(path=None) -> pd.DataFrame:
    return pd.read_csv(path or REFERENCE / "non_nuclear_features.csv")


def build_library(ensdf=None, features=None, min_intensity=0.0) -> pd.DataFrame:
    """Every ENSDF line above min_intensity (% per decay) plus the non-nuclear features
    (annihilation, escape peaks, X-rays, neutron-induced Ge lines), as one table."""
    ensdf = load_ensdf() if ensdf is None else ensdf
    features = load_features() if features is None else features
    a = ensdf[ensdf.I_pct >= min_intensity][["E_keV", "source_name", "chain", "I_pct"]]
    a = a.rename(columns={"source_name": "source", "chain": "origin"})
    b = features.rename(columns={"label": "source"}).assign(I_pct=np.nan)[["E_keV", "source", "origin", "I_pct"]]
    return pd.concat([a, b], ignore_index=True).sort_values("E_keV").reset_index(drop=True)


# Strong, single, well-known lines present in any natural-background, soil or sediment spectrum.
# Energies (keV) and uncertainties from ENSDF (Nuclear Data Sheets); 2021 A = 214 evaluation.
# These anchor the automatic calibration and are the default calibration lines. For a sealed
# calibration source pass its own list instead (e.g. Co-60 1173.228 and 1332.492 keV).
NATURAL_LINES = [
    ("Pb-212", 238.632, 0.002), ("Pb-214", 295.224, 0.002), ("Ac-228", 338.320, 0.003),
    ("Pb-214", 351.932, 0.021), ("Tl-208", 583.187, 0.002), ("Bi-214", 609.321, 0.007),
    ("Ac-228", 911.204, 0.004), ("Ac-228", 968.971, 0.017), ("Bi-214", 1120.294, 0.006),
    ("Bi-214", 1238.122, 0.010), ("K-40", 1460.820, 0.005), ("Bi-214", 1764.491, 0.014),
    ("Tl-208", 2614.511, 0.010)]


# ============================================================== 4. peak shapes
# A photopeak is the response of the detector to photons that deposit their full energy E.
# Ideal charge collection gives a Gaussian in pulse height (carrier statistics plus noise).
# Trapping removes a random, exponentially distributed amount of charge from a small fraction of
# events, which skews the peak towards LOW pulse height. Convolving the Gaussian with an
# exponential of mean beta on the low side gives the exponentially modified Gaussian (EMG):
#     g(x) = 1/(2 beta) exp(u/beta + sigma^2/(2 beta^2)) erfc(u/(sqrt2 sigma) + sigma/(sqrt2 beta)),
# with u = x - mu. It tends to the Gaussian as beta -> 0.

def gauss_pdf(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (SQRT2PI * sigma)


def emg_pdf(x, mu, sigma, beta):
    """Unit-area EMG with a LOW-energy tail of mean length beta (same units as x).

    exp(.) and erfc(.) can overflow and underflow separately, so for z > 0 the product is written as
    erfcx(z) exp(-u^2 / 2 sigma^2) with erfcx(z) = exp(z^2) erfc(z), which is the same number in a
    safe form. For z <= 0 the direct form is already safe (its exponent is then negative)."""
    x = np.asarray(x, float)
    if beta <= 1e-6 * sigma:
        return gauss_pdf(x, mu, sigma)
    u = x - mu
    z = u / (SQRT2 * sigma) + sigma / (SQRT2 * beta)
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        safe = special.erfcx(np.clip(z, 0, None)) * np.exp(-u**2 / (2 * sigma**2))
        direct = np.exp(np.minimum(u / beta + sigma**2 / (2 * beta**2), 700.0)) * special.erfc(z)
    return np.where(z > 0, safe, direct) / (2 * beta)


def emg_cdf(x, mu, sigma, beta):
    """Cumulative distribution of the EMG: F(x) = Phi((x - mu)/sigma) + beta g(x).

    Proof sketch: g(x) = int_0^inf (1/beta) e^(-t/beta) phi(x + t) dt. Differentiating and integrating
    by parts gives beta g'(x) = g(x) - phi(x), so d/dx [Phi + beta g] = phi + g - phi = g, and the
    limits are right (F -> 0 at -inf, F -> 1 at +inf). This closed form lets us integrate the peak
    exactly over each channel instead of sampling it at the channel centre."""
    x = np.asarray(x, float)
    base = special.ndtr((x - mu) / sigma)
    if beta <= 1e-6 * sigma:
        return base
    return base + beta * emg_pdf(x, mu, sigma, beta)


def peak_in_bins(x, mu, sigma, beta=0.0):
    """Fraction of a unit-area peak that falls in channel x, i.e. between x - 1/2 and x + 1/2."""
    return emg_cdf(x + 0.5, mu, sigma, beta) - emg_cdf(x - 0.5, mu, sigma, beta)


def step_shape(x, mu, sigma):
    """Smoothed step: 1 far below the peak, 0 far above it (small-angle scattering in front of the
    crystal and partial charge collection put counts just below the peak, never above it)."""
    return 0.5 * special.erfc((x - mu) / (SQRT2 * sigma))


def emg_mode(mu, sigma, beta):
    """Most probable pulse height of the EMG. From beta g' = g - phi, the slope is zero where the
    EMG density equals the underlying Gaussian density, which happens once, below mu."""
    if beta <= 1e-6 * sigma:
        return float(mu)
    f = lambda t: emg_pdf(t, mu, sigma, beta) - gauss_pdf(t, mu, sigma)
    lo = mu - sigma
    while f(lo) < 0:                     # step down until the tail density exceeds the Gaussian
        lo -= sigma
    return float(optimize.brentq(f, lo, mu, xtol=1e-12))


def emg_window_fraction(sigma, beta, half):
    """Fraction of a left-tailed EMG peak (Gaussian width sigma, tail length beta) that lies within mode +/- half."""
    d = emg_mode(0.0, sigma, beta)                     # mode relative to mu, negative for a left tail
    return float(emg_cdf(d + half, 0.0, sigma, beta) - emg_cdf(d - half, 0.0, sigma, beta))


def emg_fwhm(mu, sigma, beta):
    """Full width at half maximum of the EMG (equals 2.3548 sigma for beta = 0)."""
    if beta <= 1e-6 * sigma:
        return FWHM_PER_SIGMA * sigma
    m = emg_mode(mu, sigma, beta)
    half = 0.5 * emg_pdf(m, mu, sigma, beta)
    f = lambda t: emg_pdf(t, mu, sigma, beta) - half
    low = optimize.brentq(f, m - 20 * (sigma + beta), m)
    high = optimize.brentq(f, m, m + 10 * sigma)
    return float(high - low)


# ============================================================== 5. fitting a region of interest
SHARED = ["b0", "b1", "b2", "sigma", "beta", "step"]     # then (A_k, mu_k) for every peak


def roi_model(theta, x, n_peaks, x0):
    """Expected counts per channel: polynomial continuum + n_peaks EMG (or Gaussian) peaks with a
    common sigma and beta + a smoothed step under each peak proportional to its area."""
    b0, b1, b2, sigma, beta, step = theta[:6]
    f = b0 + b1 * (x - x0) + b2 * (x - x0) ** 2
    for k in range(n_peaks):
        A, mu = theta[6 + 2 * k], theta[7 + 2 * k]
        f = f + A * peak_in_bins(x, mu, sigma, beta)
        if step != 0:
            f = f + step * A / (SQRT2PI * sigma) * step_shape(x, mu, sigma)
    return f


def poisson_deviance(y, f):
    """D = 2 sum [f - y + y ln(y/f)]: -2 ln(likelihood ratio) against a perfect model. For a correct
    model it follows a chi-squared distribution with (channels - free parameters) degrees of freedom."""
    f = np.clip(f, 1e-12, None)
    return 2 * np.sum(f - y + np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / f), 0.0))


def fit_roi(counts, centres, sigma0, *, shape="emg", beta0=None, step="free", sigma_free=True,
            beta_free=True, background="linear", half=7.0, lo=None, hi=None, fixed_peaks=(), tail_ratio=None):
    """Poisson maximum-likelihood fit of one region of interest (ROI).

    counts      : whole spectrum (array of counts per channel)
    centres     : starting channel(s) of the peak(s) to fit together
    sigma0      : starting (or fixed) Gaussian sigma, channels
    shape       : "emg" (Gaussian with low-energy exponential tail) or "gauss"
    beta0       : starting (or fixed) tail length, channels (default 0.5 sigma0)
    step        : "free" to fit the step height, or a number to fix it (0 for no step)
    sigma_free, beta_free : fit these shape parameters, or hold them at sigma0 and beta0
    background  : "linear" or "quadratic" continuum
    half        : ROI half-width in units of sigma0 beyond the outermost peaks (ignored if lo, hi given)
    fixed_peaks : list of (area, channel) for known weak neighbours to include with fixed values
    tail_ratio  : if given, the tail length is tied to the width, beta = tail_ratio * sigma, instead of
                  being a free or fixed number (used when one shape is shared by many peaks)

    Why maximum likelihood with Poisson statistics: channels hold few counts in the wings and the
    continuum, where Neyman's or Pearson's chi-squared is biased (session-1 notebooks, Section 4).
    Minimising the Poisson deviance gives unbiased areas and centroids at any count level.
    Uncertainties come from the inverse of the expected Fisher information, J^T W J with W = 1/f.
    """
    centres = np.atleast_1d(np.asarray(centres, float))
    n = len(centres)
    beta0 = 0.5 * sigma0 if beta0 is None else beta0
    # The window depends only on sigma0 and half, never on the shape, so that fits with different
    # models of the same peak see exactly the same channels and their deviances can be compared.
    if lo is None:
        lo = max(int(np.floor(centres.min() - half * sigma0)), 0)
    if hi is None:
        hi = min(int(np.ceil(centres.max() + half * sigma0)), len(counts) - 1)
    x = np.arange(lo, hi + 1, dtype=float)
    y = np.asarray(counts[lo:hi + 1], float)
    x0 = 0.5 * (lo + hi)
    fixed_peaks = list(fixed_peaks)

    m = max(3, int(sigma0))
    left, right = y[:m].mean(), y[-m:].mean()
    theta0 = [0.5 * (left + right), (right - left) / (x[-1] - x[0]), 0.0, sigma0,
              beta0 if shape == "emg" else 0.0, 0.0 if step == "free" else float(step)]
    for c in centres:
        near = np.abs(x - c) <= 2.5 * sigma0
        theta0 += [max((y[near] - 0.5 * (left + right)).sum(), 5.0), c]
    theta0 = np.array(theta0, float)

    free = np.ones(len(theta0), bool)
    free[2] = background == "quadratic"
    free[3] = sigma_free
    free[4] = shape == "emg" and beta_free and tail_ratio is None
    free[5] = step == "free"
    if shape == "emg" and tail_ratio is not None:
        theta0[4] = tail_ratio * sigma0
    lb = np.full(len(theta0), -np.inf)
    ub = np.full(len(theta0), np.inf)
    lb[3], ub[3] = 0.3 * sigma0, 3.0 * sigma0
    lb[4], ub[4] = 0.02 * sigma0, 6.0 * sigma0
    lb[5], ub[5] = -0.1, 0.5
    for k, c in enumerate(centres):
        lb[6 + 2 * k] = 0.0
        lb[7 + 2 * k], ub[7 + 2 * k] = c - 3 * sigma0, c + 3 * sigma0

    def full(p):
        t = theta0.copy()
        t[free] = p
        if shape == "emg" and tail_ratio is not None:
            t[4] = tail_ratio * t[3]
        return t

    def model(t):
        f = roi_model(t, x, n, x0)
        for A_fix, c_fix in fixed_peaks:
            f = f + A_fix * peak_in_bins(x, c_fix, t[3], t[4])
        return f

    def deviance_residuals(p):
        f = np.clip(model(full(p)), 1e-12, None)
        d = 2 * (f - y + np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / f), 0.0))
        return np.sign(y - f) * np.sqrt(np.clip(d, 0, None))

    p0 = np.clip(theta0[free], lb[free] + 1e-9, ub[free] - 1e-9)
    sol = optimize.least_squares(deviance_residuals, p0, bounds=(lb[free], ub[free]), x_scale="jac", max_nfev=20000)
    theta = full(sol.x)
    f = np.clip(model(theta), 1e-12, None)

    def tie(t):
        if shape == "emg" and tail_ratio is not None:
            t = t.copy()
            t[4] = tail_ratio * t[3]
        return t

    idx = np.flatnonzero(free)
    J = np.empty((len(x), len(idx)))
    for j, i in enumerate(idx):
        h = 1e-6 * max(abs(theta[i]), 1e-3)
        tp, tm = theta.copy(), theta.copy()
        tp[i] += h
        tm[i] -= h
        J[:, j] = (model(tie(tp)) - model(tie(tm))) / (2 * h)
    cov = np.zeros((len(theta), len(theta)))
    try:
        cov[np.ix_(idx, idx)] = np.linalg.inv(J.T @ (J / f[:, None]))
    except np.linalg.LinAlgError:
        cov[np.ix_(idx, idx)] = np.linalg.pinv(J.T @ (J / f[:, None]))
    dev = float(poisson_deviance(y, f))
    ndf = len(x) - len(idx)
    n_nonpositive = int(np.sum(model(theta) <= 0))     # channels where the clip at 1e-12 would hide a negative model
    names = SHARED + [f"{p}{k + 1}" for k in range(n) for p in ("A", "mu")]
    at_bound = [names[i] for i in idx if np.isclose(theta[i], lb[i], atol=1e-6) or np.isclose(theta[i], ub[i], atol=1e-6)]
    return dict(theta=theta, err=np.sqrt(np.clip(np.diag(cov), 0, None)), cov=cov, names=names, free=free,
                x=x, y=y, model=f, x0=x0, n_peaks=n, shape=shape, tail_ratio=tail_ratio if shape == "emg" else None,
                dev=dev, ndf=ndf, k=len(idx), n_nonpositive=n_nonpositive,
                pval=float(stats.chi2.sf(dev, ndf)), ok=bool(sol.success), at_bound=at_bound,
                aic=dev + 2 * len(idx), bic=dev + len(idx) * np.log(len(x)))


def _gradient_error(func, theta, cov, idx):
    """First-order error propagation sigma_f^2 = g^T C g, with the gradient g taken numerically."""
    g = np.zeros(len(theta))
    for j in idx:
        h = 1e-5 * max(abs(theta[j]), 1.0)
        tp, tm = theta.copy(), theta.copy()
        tp[j] += h
        tm[j] -= h
        g[j] = (func(tp) - func(tm)) / (2 * h)
    return float(np.sqrt(max(g @ cov @ g, 0.0)))


def peak(fit, k=0) -> dict:
    """Numbers for peak k of a fit: area, the three position definitions, widths, with errors.

    position = mode (most probable pulse height), the definition used for calibration;
    mu       = centre of the underlying Gaussian (lies above the visible top for an EMG);
    mean     = mu - beta (the average pulse height, pulled down by the tail)."""
    t, C = fit["theta"], fit["cov"]
    iA, imu = 6 + 2 * k, 7 + 2 * k
    r = fit.get("tail_ratio")
    b_of = (lambda q: r * q[3]) if r is not None else (lambda q: q[4])   # beta, tied to sigma or not
    mode_of = lambda q: emg_mode(q[imu], q[3], b_of(q))
    fwhm_of = lambda q: emg_fwhm(q[imu], q[3], b_of(q))
    mean_of = lambda q: q[imu] - b_of(q)
    return dict(A=t[iA], sA=fit["err"][iA], mu=t[imu], smu=fit["err"][imu],
                position=mode_of(t), sposition=_gradient_error(mode_of, t, C, [imu, 3, 4]),
                mean=mean_of(t), smean=_gradient_error(mean_of, t, C, [imu, 4]),
                sigma=t[3], ssigma=fit["err"][3], beta=b_of(t),
                sbeta=fit["err"][3] * r if r is not None else fit["err"][4],
                fwhm=fwhm_of(t), sfwhm=_gradient_error(fwhm_of, t, C, [3, 4]), step=t[5], sstep=fit["err"][5])


def plot_fit(fit, ax, axr, title="", gauss_overlay=None):
    """Data with Poisson error bars, fitted model and deviance-free normalised residuals (n - f)/sqrt(f)."""
    x, y, f = fit["x"], fit["y"], fit["model"]
    ax.errorbar(x, y, yerr=np.sqrt(np.maximum(y, 1)), fmt=".", ms=3, color="k", lw=0.6)
    ax.step(x, f, where="mid", color="tab:red", lw=1.2, label=fit["shape"])
    if gauss_overlay is not None:
        ax.step(gauss_overlay["x"], gauss_overlay["model"], where="mid", color="tab:blue", lw=0.9, ls="--", label="gauss")
    ax.set_title(title, fontsize=8.5)
    axr.bar(x, (y - f) / np.sqrt(f), width=1, color="tab:gray")
    axr.axhline(0, color="k", lw=0.6)
    axr.set_ylim(-4, 4)
    axr.set_yticks([-3, 0, 3])
    axr.tick_params(labelsize=7)


def simulate_roi(theta, x_len, n_peaks, rng):
    """Poisson spectrum drawn from the model with known parameters (for testing the fitter)."""
    x = np.arange(x_len, dtype=float)
    return rng.poisson(roi_model(np.asarray(theta, float), x, n_peaks, x.mean()))


# ============================================================== 6. automatic first calibration
def matched_filter(counts, sigma):
    """Significance S of a peak of width sigma (channels) at every channel.

    The kernel is the negative second derivative of a Gaussian (a 'Mexican hat') with zero sum, so a
    flat or straight continuum gives no response while a peak gives a positive one. S is the filter
    output divided by its Poisson standard deviation, so S ~ N(0, 1) where there is no peak."""
    y = np.asarray(counts, float)
    j = np.arange(-int(4 * sigma), int(4 * sigma) + 1)
    k = (1 - j**2 / sigma**2) * np.exp(-j**2 / (2 * sigma**2))
    k -= k.mean()
    out = np.convolve(y, k[::-1], mode="same")
    var = np.convolve(y, (k**2)[::-1], mode="same")
    return out / np.sqrt(np.clip(var, 1, None))


def multiscale_search(counts, sigmas=np.arange(1.5, 12.01, 0.5), threshold=5.0, skip_low=30):
    """Peak search that needs no energy calibration: the matched filter is run for a range of widths,
    each channel keeps its best significance and the width that gave it, and local maxima above the
    threshold are returned as (channel, S, sigma) sorted by significance."""
    S_all = np.array([matched_filter(counts, s) for s in sigmas])
    best = S_all.argmax(axis=0)
    S = S_all[best, np.arange(len(counts))]
    width = sigmas[best]
    first = np.flatnonzero(counts)[0] + skip_low            # skip the ADC threshold edge
    out = []
    for i in range(first, len(counts) - 10):
        w = max(2, int(1.5 * width[i]))
        if S[i] > threshold and S[i] == S[max(i - w, 0):i + w + 1].max():
            out.append((i, float(S[i]), float(width[i])))
    return sorted(out, key=lambda r: -r[1]), S, width


def quick_centroid(counts, i, sigma):
    """Centre of a peak from the first moment of the counts above a straight floor (no fit).
    Good to a fraction of a channel, which is all that line matching needs."""
    w = int(np.ceil(2 * sigma))
    n = len(counts)
    xs = np.arange(max(i - w, 0), min(i + w, n - 1) + 1)          # clipped at the ends of the spectrum
    y = counts[xs].astype(float)
    side = np.r_[counts[max(i - 3 * w, 0):max(i - 2 * w, 0)], counts[min(i + 2 * w + 1, n):min(i + 3 * w + 1, n)]].astype(float)
    net = np.clip(y - (side.mean() if len(side) else 0.0), 0, None)
    return float((xs * net).sum() / net.sum()) if net.sum() > 0 else float(i)


def match_lines(peak_ch, peak_S, ref_E, gain_range=(0.01, 5.0), max_offset=30.0, tol_ch=3.0, tol_keV=0.8):
    """Find the linear map E = g0 + g1 * channel that puts the most reference lines on found peaks.

    This is the mentor's 'Delta E over Delta channel' idea made systematic. Any two found peaks and
    any two reference lines define a candidate gain g1 = (E_b - E_a)/(c_j - c_i) and offset g0. Only
    the right assignment makes many OTHER reference lines land on found peaks as well, so each
    candidate is scored by the number of lines it explains (ties broken by the summed significance
    of the matched peaks). Random coincidences are rare because peaks are sparse: with 30 peaks in
    16 384 channels a line has about a 1 % chance of landing within 3 channels of one by accident.
    Returns (g0, g1, matches) where matches lists (reference index, peak index)."""
    c = np.asarray(peak_ch, float)
    S = np.asarray(peak_S, float)
    E = np.asarray(ref_E, float)
    best = (-1, -np.inf, None)
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            dc = c[j] - c[i]
            if dc == 0:
                continue
            for a in range(len(E)):
                for b in range(a + 1, len(E)):
                    g1 = (E[b] - E[a]) / dc
                    if not gain_range[0] <= g1 <= gain_range[1]:
                        continue
                    g0 = E[a] - g1 * c[i]
                    if abs(g0) > max_offset:
                        continue
                    pred = (E - g0) / g1
                    d = np.abs(pred[:, None] - c[None, :])
                    hit = d.min(axis=1) < max(tol_ch, tol_keV / g1)
                    score = (int(hit.sum()), float(np.log(S[d.argmin(axis=1)][hit]).sum()))
                    if score > best[:2]:
                        best = (score[0], score[1], (g0, g1))
    g0, g1 = best[2]
    for _ in range(3):                                   # refine on all matches, then re-match
        pred = (E - g0) / g1
        d = np.abs(pred[:, None] - c[None, :])
        pairs = [(r, int(d[r].argmin())) for r in range(len(E)) if d[r].min() < max(tol_ch, tol_keV / g1)]
        cc = np.array([c[p] for _, p in pairs])
        ee = np.array([E[r] for r, _ in pairs])
        g1, g0 = np.polyfit(cc, ee, 1)
    return float(g0), float(g1), pairs


# Known neighbours that sit inside the fitting window of a calibration line and must be fitted with it.
NEIGHBOURS = {238.632: [241.995],      # Pb-214 242.0 (blended with Ra-224 241.0)
              295.224: [300.087],      # Pb-212 300.1
              968.971: [964.766]}      # Ac-228 964.8


def fit_lines(spec: Spectrum, rough: Calibration, lines=NATURAL_LINES, tail_ratio=0.0, step=0.0,
              neighbours=NEIGHBOURS, fwhm_guess=lambda E: 0.8 + 0.0005 * E, half=7.0):
    """Fit every reference line at the channel where the rough calibration puts it.

    Width free for each line; tail tied to the width by the shared tail_ratio (0 = Gaussian);
    step fixed (0 = none) or "free". Returns a table with the three position definitions, width,
    area and fit quality for each line, plus the fits themselves."""
    rows, fits = [], {}
    for name, E, sE in lines:
        x_e = float(rough.channel(E))
        s0 = fwhm_guess(E) / FWHM_PER_SIGMA / float(rough.gain(x_e))
        centres = [x_e] + [float(rough.channel(e)) for e in neighbours.get(E, [])]
        emg = tail_ratio > 0
        fit = fit_roi(spec.counts, centres, s0, shape="emg" if emg else "gauss", tail_ratio=tail_ratio if emg else None,
                      step=step, half=half)
        p = peak(fit)
        rows.append(dict(nuclide=name, E_ref=E, sE_ref=sE, position=p["position"], sposition=p["sposition"],
                         mu=p["mu"], smu=p["smu"], mean=p["mean"], sigma_ch=p["sigma"], ssigma_ch=p["ssigma"],
                         beta_ch=p["beta"], fwhm_ch=p["fwhm"], sfwhm_ch=p["sfwhm"], step=p["step"],
                         area=p["A"], sarea=p["sA"], signif=p["A"] / p["sA"] if p["sA"] > 0 else 0.0,
                         offset_ch=p["position"] - x_e, dev=fit["dev"], ndf=fit["ndf"], pval=fit["pval"],
                         at_bound=",".join(fit["at_bound"])))
        fits[E] = fit
    return pd.DataFrame(rows), fits


def select_lines(table, min_signif=5.0, max_offset_sigma=3.0, min_pval=1e-4):
    """Lines good enough to calibrate with: clearly significant, found close to where the rough map
    expected them (a large offset means the fit locked onto something else), an acceptable fit,
    and no shape parameter stuck at a bound."""
    return ((table.signif >= min_signif) & (np.abs(table.offset_ch) <= max_offset_sigma * table.sigma_ch)
            & (table.pval >= min_pval) & ~table.at_bound.str.contains("sigma"))


def auto_calibrate(counts, lines=NATURAL_LINES, n_peaks=30, threshold=5.0, **kw):
    """Rough linear calibration straight from the counts: multiscale search, quick centroids of the
    n_peaks most significant peaks, then pattern matching against the reference lines."""
    found, _, _ = multiscale_search(counts, threshold=threshold)
    found = found[:n_peaks]
    ch = [quick_centroid(counts, i, s) for i, _, s in found]
    S = [s for _, s, _ in found]
    E = [e for _, e, _ in lines]
    g0, g1, pairs = match_lines(ch, S, E, **kw)
    table = pd.DataFrame([dict(nuclide=lines[r][0], E_ref=lines[r][1], channel=ch[p], S=S[p],
                               E_rough=g0 + g1 * ch[p]) for r, p in pairs])
    return g0, g1, table, found


# ============================================================== 7. energy calibration fit
@dataclass
class Calibration:
    """Polynomial energy calibration E(x) = sum_j a_j x^j on the channel-centre coordinate x,
    with the covariance matrix C of the coefficients."""
    a: np.ndarray
    C: np.ndarray
    chi2: float
    ndf: int
    scale: float = 1.0               # sqrt(chi2/ndf) already applied to C when > 1 (PDG convention)
    sys: float = 0.0                 # per-line systematic (keV) added in quadrature, e.g. ADC non-linearity

    @property
    def degree(self):
        return len(self.a) - 1

    def energy(self, x):
        x = np.asarray(x, float)
        return sum(c * x**j for j, c in enumerate(self.a))

    def gain(self, x):
        """dE/dx, keV per channel, at position x."""
        x = np.asarray(x, float)
        return sum(j * c * x**(j - 1) for j, c in enumerate(self.a) if j > 0)

    def sigma(self, x):
        """Standard uncertainty of E(x) from the coefficient covariance."""
        x = np.asarray(x, float)
        Jm = np.stack([x**j for j in range(len(self.a))], axis=-1)
        return np.sqrt(np.einsum("...i,ij,...j->...", Jm, self.C, Jm))

    def channel(self, E):
        """Inverse: position x at which E(x) = E (Newton iterations from the linear estimate)."""
        E = np.asarray(E, float)
        x = (E - self.a[0]) / self.a[1]
        for _ in range(20):
            x = x - (self.energy(x) - E) / self.gain(x)
        return x

    def lower_edge_coefficients(self):
        """Coefficients of E at the LOWER edge of channel i, E(i - 1/2), as InterSpec and the IAEA
        .Spe/.CALp files expect. For E = a0 + a1 x + a2 x^2 this gives
        (a0 - a1/2 + a2/4, a1 - a2, a2)."""
        a = np.r_[self.a, np.zeros(3 - len(self.a))] if len(self.a) < 3 else self.a
        a0, a1, a2 = a[:3]
        return (a0 - a1 / 2 + a2 / 4, a1 - a2, a2)

    def to_dict(self):
        return dict(a=list(map(float, self.a)), C=np.asarray(self.C).tolist(), chi2=self.chi2, ndf=self.ndf,
                    scale=self.scale, sys=self.sys)

    @staticmethod
    def from_dict(d):
        return Calibration(np.array(d["a"]), np.array(d["C"]), d["chi2"], d["ndf"], d.get("scale", 1.0), d.get("sys", 0.0))


def linear_calibration(g0, g1):
    """A calibration object from a rough two-number map (no covariance)."""
    return Calibration(np.array([g0, g1]), np.zeros((2, 2)), np.nan, 0)


def fit_calibration(x, sx, E, sE, degree=1, n_iter=4, scale_by_chi2=True, sys_keV=0.0) -> Calibration:
    """Weighted least-squares polynomial E(x) with the 'effective variance' method.

    Each point has an uncertainty in E (the reference energy, small) and in x (the fitted position).
    The x error is turned into an energy error with the local gain, s^2 = sE^2 + (dE/dx sx)^2, and the
    fit is repeated with the updated gain until it settles.

    Extra scatter can be handled in one of two ways.
      * scale_by_chi2 (PDG): if chi2/ndf > 1 the covariance is multiplied by chi2/ndf.
      * sys_keV: a per-line systematic added in quadrature to every point, s^2 += sys^2. With
        sys_keV="auto" it is chosen so that chi2/ndf = 1, which is the natural description of a
        detector whose energy scale has small wiggles (ADC integral non-linearity) that no low-order
        polynomial can follow. The value is stored in Calibration.sys and must be added to the
        uncertainty of any energy measured with this calibration."""
    x, sx, E, sE = (np.asarray(v, float) for v in (x, sx, E, sE))
    A = np.vander(x, degree + 1, increasing=True)
    ndf = len(E) - degree - 1

    def solve(sys):
        slope = np.full_like(x, np.polyfit(x, E, 1)[0])
        for _ in range(n_iter):
            s = np.sqrt(sE**2 + (slope * sx)**2 + sys**2)
            W = 1 / s**2
            C = np.linalg.inv(A.T @ (A * W[:, None]))
            a = C @ A.T @ (W * E)
            slope = sum(j * a[j] * x**(j - 1) for j in range(1, degree + 1))
        return a, C, float(np.sum(((E - A @ a) / s)**2))

    sys = 0.0 if sys_keV == "auto" else float(sys_keV)
    a, C, chi2 = solve(sys)
    if sys_keV == "auto" and ndf > 0 and chi2 / ndf > 1:
        sys = optimize.brentq(lambda v: solve(v)[2] / ndf - 1.0, 0.0, 10.0, xtol=1e-6)
        a, C, chi2 = solve(sys)
        scale_by_chi2 = False
    S = max(1.0, np.sqrt(chi2 / ndf)) if (scale_by_chi2 and ndf > 0) else 1.0
    return Calibration(a=a, C=C * S**2, chi2=chi2, ndf=ndf, scale=S, sys=sys)


def fit_calibration_borrowed_shape(x, sx, E, sE, reference: Calibration, sys_keV="auto") -> Calibration:
    """Offset and gain only, with the non-linear terms (x^2 and above) taken from a reference
    calibration of the same electronics. A non-linearity of the ADC belongs to the hardware, not to
    the run, so a spectrum with too few strong lines to measure it can borrow it from one that has
    many. The uncertainty of the borrowed terms is neglected (it is small compared with the
    systematic floor)."""
    x, sx, E, sE = (np.asarray(v, float) for v in (x, sx, E, sE))
    higher = lambda xx: sum(c * np.asarray(xx, float)**j for j, c in enumerate(reference.a) if j >= 2)
    lin = fit_calibration(x, sx, E - higher(x), sE, degree=1, sys_keV=sys_keV)
    a = np.r_[lin.a, reference.a[2:]]
    C = np.zeros((len(a), len(a)))
    C[:2, :2] = lin.C
    return Calibration(a=a, C=C, chi2=lin.chi2, ndf=lin.ndf, scale=lin.scale, sys=lin.sys)


def f_test(simple: Calibration, rich: Calibration):
    """Does the extra polynomial term reduce chi2 by more than chance? F = (chi2_1 - chi2_2)/(chi2_2/ndf_2)."""
    F = (simple.chi2 - rich.chi2) / (rich.chi2 / rich.ndf)
    return float(F), float(stats.f.sf(F, simple.ndf - rich.ndf, rich.ndf))


# ============================================================== 8. resolution and tail models
def fwhm_model(E, w0, w1, w2):
    """FWHM^2 = w0^2 + w1^2 E + w2^2 E^2: electronic noise + carrier statistics + charge collection/drift."""
    E = np.clip(np.asarray(E, float), 0, None)
    return np.sqrt(w0**2 + w1**2 * E + w2**2 * E**2)


def fano_from_w1(w1, sw1=0.0):
    """Fano factor from the carrier-statistics term, w1^2 = 2.3548^2 F eps."""
    k = FWHM_PER_SIGMA**2 * EPS_GE
    return w1**2 / k, 2 * w1 * sw1 / k


def fit_resolution(E, fwhm, sfwhm):
    """Fit the three-term FWHM model (w terms forced positive by fitting their squares' roots)."""
    p, C = optimize.curve_fit(fwhm_model, E, fwhm, p0=[1.0, 0.04, 3e-4], sigma=sfwhm, absolute_sigma=True,
                              bounds=([0, 0, 0], [10, 1, 0.01]))
    chi2 = float(np.sum(((fwhm - fwhm_model(E, *p)) / sfwhm)**2))
    return p, C, chi2, len(E) - 3


def fit_tail(E, beta, sbeta):
    """Weighted straight line beta(E) = t0 + t1 E for the tail length in keV (an empirical model;
    trapping losses grow with the number of carriers, hence roughly with E)."""
    W = 1 / np.asarray(sbeta, float)**2
    A = np.vstack([np.ones_like(E), E]).T
    C = np.linalg.inv(A.T @ (A * W[:, None]))
    t = C @ A.T @ (W * beta)
    chi2 = float(np.sum(W * (beta - A @ t)**2))
    return t, C, chi2, len(E) - 2


# ============================================================== 9. peak search, inventory, identification
def group_candidates(cands, sigma_ch, n_sigma=10.0) -> list:
    """Candidates closer than n_sigma are fitted together; a closer neighbour would sit inside the
    other's ROI and distort its continuum if it were left out of the model."""
    if not cands:
        return []
    groups, current = [], [cands[0]]
    for i in cands[1:]:
        if i - current[-1] < n_sigma * sigma_ch[i]:
            current.append(i)
        else:
            groups.append(current)
            current = [i]
    groups.append(current)
    return groups


def calibrated_search(counts, cal: Calibration, sigma_keV, threshold=5.0, skip_low=30):
    """Matched-filter search with the width each channel should have from the resolution model."""
    E = cal.energy(np.arange(len(counts)))
    sig_ch = sigma_keV(E) / cal.gain(np.arange(len(counts)))
    sigmas = np.arange(1.0, 15.01, 0.25)
    S_all = np.array([matched_filter(counts, s) for s in sigmas])
    pick = np.clip(np.round((sig_ch - sigmas[0]) / 0.25).astype(int), 0, len(sigmas) - 1)
    S = S_all[pick, np.arange(len(counts))]
    first = np.flatnonzero(counts)[0] + skip_low
    cands = [i for i in range(first, len(counts) - 10)
             if S[i] > threshold and S[i] == S[max(i - int(1.5 * sig_ch[i]), 0):i + int(1.5 * sig_ch[i]) + 1].max()]
    return cands, S, sig_ch


def identify(E, sE, library: pd.DataFrame, window=None):
    """Library lines within max(0.5 keV, 3 sE) of E, nearest first."""
    w = max(0.5, 3 * sE) if window is None else window
    near = library[np.abs(library.E_keV - E) < w].copy()
    near["d"] = np.abs(near.E_keV - E)
    return near.sort_values(["d", "I_pct"], ascending=[True, False])


def fit_inventory(spec: Spectrum, cal: Calibration, sigma_keV, tail_ratio, library, threshold=5.0,
                  shape="emg", step=0.0, free_width_S=10.0):
    """Fit every candidate peak and match each fitted energy to the library.

    Weak peaks (and peaks in multiplets) are fitted with the width fixed from the resolution model and
    the tail tied to it (beta = tail_ratio * sigma), which keeps their fits stable. Strong isolated
    peaks (S > free_width_S) are fitted with the width free, and the ratio of the fitted width to the
    model is reported: a ratio well above 1 flags a peak that is not an ordinary photopeak, such as
    the Doppler-broadened 511 keV annihilation line, whose area a fixed width would underestimate."""
    counts = spec.counts
    cands, S, sig_ch = calibrated_search(counts, cal, sigma_keV, threshold)
    emg = shape == "emg" and tail_ratio > 0
    rows = []
    for g in group_candidates(cands, sig_ch):
        s0 = float(sig_ch[g[0]])
        free = len(g) == 1 and S[g[0]] > free_width_S
        fit = fit_roi(counts, g, s0, shape="emg" if emg else "gauss", tail_ratio=tail_ratio if emg else None,
                      step=step, sigma_free=free, half=7.0)
        for k, i in enumerate(g):
            p = peak(fit, k)
            E = float(cal.energy(p["position"]))
            sE_stat = float(cal.gain(p["position"]) * p["sposition"])
            sE = float(np.sqrt(sE_stat**2 + cal.sigma(p["position"])**2 + cal.sys**2))
            near = identify(E, sE, library)
            best = near.iloc[0] if len(near) else None
            rows.append(dict(E_keV=E, sE_keV=sE, sE_stat=sE_stat, position=p["position"], sposition=p["sposition"],
                             S=float(S[i]), area=p["A"], s_area=p["sA"],
                             width_ratio=p["sigma"] / s0 if free else np.nan,
                             dev_ndf=fit["dev"] / fit["ndf"], p_fit=fit["pval"],
                             source=best.source if best is not None else "unidentified",
                             origin=best.origin if best is not None else "",
                             E_lib=float(best.E_keV) if best is not None else np.nan,
                             dE=E - float(best.E_keV) if best is not None else np.nan,
                             other=", ".join(f"{r.source} {r.E_keV:.1f}" for r in near.iloc[1:4].itertuples())))
    df = pd.DataFrame(rows)
    if len(df):
        df["rate_cpd"] = df.area / spec.live * 86400
        df["s_rate_cpd"] = df.s_area / spec.live * 86400
    return df


def fit_line_at(spec: Spectrum, cal: Calibration, E, sigma_keV, tail_ratio, shape="emg", step=0.0,
                free_width=False, neighbours=(), half=7.0):
    """Fit one known line (plus optional known neighbours, given as energies) wherever it is
    expected, whether or not the search found it. Used for like-for-like comparisons between
    spectra: the same model at the same energy in every spectrum. The width comes from the
    resolution model (or is fitted if free_width) and the tail is tied to it."""
    x_e = float(cal.channel(E))
    s0 = float(sigma_keV(E) / cal.gain(x_e))
    emg = shape == "emg" and tail_ratio > 0
    centres = [x_e] + [float(cal.channel(e)) for e in neighbours]
    fit = fit_roi(spec.counts, centres, s0, shape="emg" if emg else "gauss", tail_ratio=tail_ratio if emg else None,
                  step=step, sigma_free=free_width, half=half)
    p = peak(fit, 0)
    p.update(E_fit=float(cal.energy(p["position"])), rate=p["A"] / spec.live, s_rate=p["sA"] / spec.live,
             dev=fit["dev"], ndf=fit["ndf"], pval=fit["pval"], n_nonpositive=fit["n_nonpositive"])
    return p, fit


# ============================================================== 10. comparing spectra
def net_rate(R_s, sR_s, R_b, sR_b):
    """Sample minus background rate (both already per second of their own live time) and its error.
    The two measurements are independent, so their variances add."""
    R = R_s - R_b
    s = np.hypot(sR_s, sR_b)
    return R, s, R / s


def ratio(a, sa, b, sb):
    """a/b with first-order error propagation (relative errors in quadrature)."""
    r = a / b
    return r, abs(r) * np.hypot(sa / a, sb / b)


def window_rate(spec: Spectrum, cal: Calibration, E_lo, E_hi):
    """Counts per second per keV in an energy window, with the Poisson error (for continua)."""
    i0 = int(np.ceil(cal.channel(E_lo)))
    i1 = int(np.floor(cal.channel(E_hi)))
    N = spec.counts[i0:i1 + 1].sum()
    width = float(cal.energy(i1 + 0.5) - cal.energy(i0 - 0.5))
    return N / spec.live / width, np.sqrt(N) / spec.live / width


def rebin_to_energy(spec: Spectrum, cal: Calibration, edges):
    """Counts per second per keV on a common energy grid, so spectra with different gains can be
    overlaid or divided. Each channel's counts are shared out by the overlap of its energy interval
    with the new bins (counts are conserved). Display only: fits always use raw channels."""
    lo = cal.energy(spec.channels - 0.5)
    hi = cal.energy(spec.channels + 0.5)
    cum = np.r_[0.0, np.cumsum(spec.counts)]
    E_edges = np.r_[lo, hi[-1]]
    c_at = np.interp(edges, E_edges, cum)
    N = np.diff(c_at)
    return N / spec.live / np.diff(edges), np.sqrt(np.clip(N, 0, None)) / spec.live / np.diff(edges)


# ============================================================== 11. exports
def write_calibrated_csv(path, spec: Spectrum, cal: Calibration, notes=()):
    """Every channel with its lower, centre and upper energy, the calibration uncertainty, counts and
    count rate. Metadata goes in '#' comment lines at the top."""
    ch = spec.channels
    df = pd.DataFrame({"channel_index": ch.astype(int),
                       "file_channel": (ch + spec.meta.get("first_channel_label", 0)).astype(int),
                       "E_low_keV": cal.energy(ch - 0.5), "E_centre_keV": cal.energy(ch),
                       "E_high_keV": cal.energy(ch + 0.5), "sigma_E_cal_keV": cal.sigma(ch),
                       "counts": spec.counts, "rate_per_s": spec.counts / spec.live})
    poly = " + ".join(f"{c:.10g} x^{j}" for j, c in enumerate(cal.a))
    head = [f"# source file: {spec.path.name}  (sha256 {spec.sha256})",
            f"# start (as written in file): {spec.start}  live time: {spec.live:.2f} s  real time: {spec.real:.2f} s",
            f"# calibration: E = {poly}  (x = channel index, channel-centre coordinate)",
            "# rate_per_s uses the live time"] + [f"# {n}" for n in notes]
    with open(path, "w") as f:
        f.write("\n".join(head) + "\n")
        df.to_csv(f, index=False, float_format="%.6f")
    return df


def write_calp(path, cal: Calibration):
    """InterSpec / PeakEasy CALp file (polynomial, lower channel edge, 10 significant figures)."""
    C0, C1, C2 = cal.lower_edge_coefficients()
    Path(path).write_text("\n".join(["#PeakEasy CALp File Ver:  4.00", f"Offset (keV)           :  {C0:1.9e}",
                                     f"Gain (keV / Chan)      :  {C1:1.9e}", f"2nd Order Coef         :  {C2:1.9e}",
                                     "#END"]) + "\n")


def write_spe(path, spec: Spectrum, cal: Calibration | None = None, title=None):
    """IAEA ASCII .Spe file that InterSpec (and most other codes) can open with the calibration
    already applied. The polynomial is written for the LOWER edge of each channel."""
    C0, C1, C2 = cal.lower_edge_coefficients() if cal is not None else (0.0, 1.0, 0.0)
    start = spec.start.strftime("%m/%d/%Y %H:%M:%S") if spec.start else "01/01/2000 00:00:00"
    out = ["$SPEC_ID:", title or spec.label, "$SPEC_REM:",
           f"converted from {spec.path.name} (sha256 {spec.sha256[:16]}...) by gamma_spec.write_spe",
           "$DATE_MEA:", start, "$MEAS_TIM:", f"{spec.live:.3f} {spec.real:.3f}",
           "$DATA:", f"0 {spec.n - 1}"]
    out += [f"{c:8d}" for c in spec.counts]
    out += ["$ENER_FIT:", f"{C0:.9e} {C1:.9e}", "$MCA_CAL:", "3", f"{C0:.9e} {C1:.9e} {C2:.9e} keV", "$ENDRECORD:"]
    Path(path).write_text("\n".join(out) + "\n")


def read_interspec_csv(path) -> pd.DataFrame:
    """Read an InterSpec Peak Manager CSV export (header row of names, then units, then peaks)."""
    import csv
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        rows = list(csv.reader(f))
    h = next(i for i, r in enumerate(rows) if r and r[0].strip() == "Centroid")
    names = [c.strip() for c in rows[h]]
    units = [c.strip() for c in rows[h + 1]] + [""] * len(names)
    cols = [f"{n}|{u}" if u else n for n, u in zip(names, units)]
    data = [r for r in rows[h + 2:] if r and r[0].strip()]
    df = pd.DataFrame([[c.strip() for c in r[:len(cols)]] for r in data], columns=cols[:len(data[0])])
    num = lambda c: pd.to_numeric(df[c], errors="coerce") if c in df else pd.Series(np.nan, index=df.index)
    return pd.DataFrame({"E_interspec": num("Centroid|keV"), "area_interspec": num("Net_Area|Counts"),
                         "s_area_interspec": num("Net_Area|Uncertainty"), "fwhm_interspec": num("FWHM|keV"),
                         "nuclide": df.get("Nuclide", ""), "E_photopeak": num("Photopeak_Energy|keV")})


def save_json(path, obj):
    """JSON writer that copes with numpy numbers and arrays."""
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (dt.datetime, dt.date, Path)):
            return str(o)
        raise TypeError(type(o))
    Path(path).write_text(json.dumps(obj, indent=2, default=conv))
