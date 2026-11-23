"""Quick checks of gamma_spec. Run from anywhere with:  python3 code/test_gamma_spec.py
Every test prints 'ok'. They use only synthetic data, except the reader test, which needs the session-2 files."""
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gamma_spec as G


def test_emg_cdf_is_integral_of_pdf():
    # F(x) = Phi + beta g must be the integral of g, and tend to the Gaussian for beta -> 0
    x = np.linspace(-40, 20, 60001)
    for beta in (0.3, 2.0, 6.0):
        g = G.emg_pdf(x, 0.0, 3.0, beta)
        F = G.emg_cdf(x, 0.0, 3.0, beta)
        num = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(x))])
        assert np.max(np.abs(num - (F - F[0]))) < 1e-6
        assert abs(F[-1] - 1) < 1e-8
    assert np.allclose(G.emg_pdf(x, 0.0, 3.0, 1e-9), G.gauss_pdf(x, 0.0, 3.0))


def test_emg_mode_mean_and_width():
    sig, beta = 4.0, 2.0
    x = np.linspace(-60, 30, 900001)
    g = G.emg_pdf(x, 0.0, sig, beta)
    assert abs(G.emg_mode(0.0, sig, beta) - x[np.argmax(g)]) < 1e-3
    assert abs(np.sum(x * g) / np.sum(g) - (-beta)) < 1e-3            # mean = mu - beta
    assert abs(G.emg_fwhm(0.0, sig, 1e-9) - G.FWHM_PER_SIGMA * sig) < 1e-9


def test_fit_recovers_emg_peak():
    rng = np.random.default_rng(3)
    theta = np.array([30.0, 0.0, 0.0, 4.0, 1.6, 0.0, 20000.0, 60.0])
    pulls_mu, pulls_A = [], []
    for _ in range(60):
        y = G.simulate_roi(theta, 121, 1, rng)
        fit = G.fit_roi(y, [60.0 + rng.uniform(-1, 1)], 4.0, shape="emg", tail_ratio=0.4, step=0.0, lo=0, hi=120)
        p = G.peak(fit)
        pulls_mu.append((p["position"] - G.emg_mode(60.0, 4.0, 1.6)) / p["sposition"])
        pulls_A.append((p["A"] - 20000.0) / p["sA"])
    assert abs(np.mean(pulls_mu)) < 0.5 and 0.6 < np.std(pulls_mu) < 1.5
    assert abs(np.mean(pulls_A)) < 0.5 and 0.6 < np.std(pulls_A) < 1.5


def test_auto_calibration_on_synthetic_spectrum():
    rng = np.random.default_rng(5)
    g0, g1, n = 0.3, 0.21, 16384
    x = np.arange(n, dtype=float)
    mean = 2000.0 * np.exp(-x / 3000.0) + 5.0
    for name, E, _ in G.NATURAL_LINES:
        c = (E - g0) / g1
        mean += 3000.0 * G.peak_in_bins(x, c, (0.9 + 0.0006 * E) / G.FWHM_PER_SIGMA / g1)
    counts = rng.poisson(mean)
    a0, a1, table, _ = G.auto_calibrate(counts)
    assert len(table) >= 11
    assert abs(a1 - g1) / g1 < 2e-4 and abs(a0 - g0) < 0.5


def test_calibration_and_lower_edge_conventions():
    x = np.array([100.0, 2000.0, 5000.0, 9000.0, 14000.0])
    a = np.array([0.2, 0.17, 6e-9])
    E = a[0] + a[1] * x + a[2] * x**2
    cal = G.fit_calibration(x, np.full(5, 0.01), E, np.full(5, 0.001), degree=2)
    assert np.allclose(cal.a, a, rtol=1e-6, atol=1e-9)
    C0, C1, C2 = cal.lower_edge_coefficients()
    i = np.arange(0, 16000, 1000.0)
    assert np.allclose(C0 + C1 * i + C2 * i**2, cal.energy(i - 0.5))  # InterSpec's polynomial = lower edge
    assert np.allclose(cal.channel(cal.energy(i)), i, atol=1e-6)


def test_systematic_floor_gives_unit_chi2():
    rng = np.random.default_rng(2)
    x = np.linspace(1000, 15000, 12)
    E = 0.1 + 0.168 * x + rng.normal(0, 0.08, x.size)               # 80 eV of extra scatter
    cal = G.fit_calibration(x, np.full(x.size, 0.05), E, np.full(x.size, 0.002), degree=1, sys_keV="auto")
    assert abs(cal.chi2 / cal.ndf - 1) < 1e-4 and 0.03 < cal.sys < 0.15


def test_spe_round_trip():
    s = G.Spectrum(np.arange(16, dtype=np.int64) * 3, 100.5, 101.25, None, (0.0, 1.0, 0.0), Path("x.csv"), "0" * 64,
                   "columns", dict(first_channel_label=1))
    cal = G.Calibration(np.array([0.5, 0.2, 1e-6]), np.zeros((3, 3)), 0.0, 1)
    with tempfile.TemporaryDirectory() as d:
        G.write_spe(Path(d) / "t.spe", s, cal)
        back = G.read_spe(Path(d) / "t.spe")
    assert np.array_equal(back.counts, s.counts) and abs(back.live - 100.5) < 1e-6
    assert np.allclose(back.stored_cal, cal.lower_edge_coefficients(), rtol=1e-8)


def test_lynx_reader_on_session2_files():
    data = HERE.parent / "data" / "raw"
    files = sorted(data.glob("*.csv"))
    if not files:
        print("   (session-2 data not found, reader test skipped)")
        return
    for f in files:
        s = G.read_spectrum(f)
        assert s.n == 16384 and s.meta["first_channel_label"] == 1 and s.meta["energy_column_is_channel_number"]
        assert s.live < s.real and s.start is not None



def test_emg_window_is_centred_on_the_mode():
    from scipy import stats
    sigma, beta, half = 1.0, 0.39751, 1.25 * 2.49712
    right = stats.exponnorm(beta / sigma)              # right-tailed law; its mirror image is the left-tailed peak
    d = G.emg_mode(0.0, sigma, beta)
    assert d < 0
    expected = right.sf(-(d + half)) - right.sf(-(d - half))
    assert abs(G.emg_window_fraction(sigma, beta, half) - expected) < 1e-6
    assert abs(G.emg_window_fraction(sigma, 1e-9, 2.0) - (stats.norm.cdf(2) - stats.norm.cdf(-2))) < 1e-6


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok ", name)
