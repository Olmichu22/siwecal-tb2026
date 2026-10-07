// Unit tests of GainAnchorFit.h on synthetic data: no ROOT, no Gaudi.
//   ctest --test-dir gaudi_source/build -R gain_anchor     (or run the binary)
#include "k4SiWEcalReco/GainAnchorFit.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>

namespace {
int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("%s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
  if (!ok) ++failures;
}

bool near(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

// Hits along y = k x + c with Gaussian scatter, x uniform in [lo, hi].
void fillLine(k4siwecal::GainAnchorHistogram& h, double k, double c, double sigma, int n, double lo, double hi,
              unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> ux(lo, hi);
  std::normal_distribution<double> g(0.0, sigma);
  for (int i = 0; i < n; ++i) {
    const double x = ux(rng);
    h.fill(x, k * x + c + g(rng));
  }
}
}  // namespace

int main() {
  using k4siwecal::GainAnchorHistogram;

  {  // a clean line is recovered
    GainAnchorHistogram h;
    fillLine(h, 0.0924, 2.35, 1.5, 400000, -50, 2000, 1);
    const auto r = k4siwecal::fitGainAnchor(h, 150, 1200, 50);
    check(r.ok, "clean line: fit succeeds");
    check(near(r.k, 0.0924, 0.0005), "clean line: k = 0.0924 (got " + std::to_string(r.k) + ")");
    check(near(r.c, 2.35, 0.3), "clean line: c = 2.35 (got " + std::to_string(r.c) + ")");
    check(r.sliceX.size() == 84, "clean line: 84 slices in [150, 1200] (got " + std::to_string(r.sliceX.size()) + ")");
  }

  {  // the median ignores a one-sided tail that would drag a least-squares fit
    GainAnchorHistogram h;
    fillLine(h, 0.0924, 2.35, 1.5, 400000, -50, 2000, 2);
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> ux(150, 1200), uy(0, 100);
    for (int i = 0; i < 60000; ++i) {  // 15 % junk far above the band
      const double x = ux(rng);
      h.fill(x, 0.0924 * x + 2.35 + 20 + uy(rng));
    }
    const auto r = k4siwecal::fitGainAnchor(h, 150, 1200, 50);
    check(r.ok && near(r.k, 0.0924, 0.002) && near(r.c, 2.35, 3.0),
          "15 % one-sided tail: line stays on the band (k " + std::to_string(r.k) + ", c " + std::to_string(r.c) + ")");
  }

  {  // only slices inside the fiducial band enter, and the band edges are inclusive on slice centres
    GainAnchorHistogram h;
    fillLine(h, 0.0924, 2.35, 1.5, 200000, 150, 1200, 4);
    fillLine(h, 0.2, 50.0, 1.5, 200000, 1300, 2000, 5);  // a different line outside the band
    const auto r = k4siwecal::fitGainAnchor(h, 150, 1200, 50);
    check(r.ok && near(r.k, 0.0924, 0.0008), "outside-band hits do not enter the fit");
    check(near(r.sliceX.front(), 156.25, 1e-9) && near(r.sliceX.back(), 1193.75, 1e-9),
          "first/last fitted slice centres are 156.25 / 1193.75");
  }

  {  // slices below MinEntriesPerSlice are skipped; too few slices fail cleanly
    GainAnchorHistogram h;
    for (int i = 0; i < 49; ++i) h.fill(500.0, 48.0);
    check(std::isnan(h.sliceMedian(static_cast<int>((500.0 + 50.0) / 12.5), 50)), "49 entries: no median");
    const auto r = k4siwecal::fitGainAnchor(h, 150, 1200, 50);
    check(!r.ok && !r.error.empty(), "no usable slice: fit fails with a message");
  }

  {  // median convention: centre of the first bin reaching half the slice
    GainAnchorHistogram h;
    for (int i = 0; i < 30; ++i) h.fill(500.0, 10.1);  // bin [10.0, 10.5)
    for (int i = 0; i < 30; ++i) h.fill(500.0, 20.1);  // bin [20.0, 20.5)
    const int ix = static_cast<int>((500.0 + 50.0) / 12.5);
    check(near(h.sliceMedian(ix, 50), 10.25, 1e-9), "median of 30+30 split is the lower bin's centre (10.25)");
  }

  {  // pedestal table: means kept only in (0, 500), missing cells are 0
    const std::string path = "test_gain_anchor_ped.txt";
    {
      std::ofstream f(path);
      f << "#header\n#layer chip channel ...\n";
      f << "3 5 7";
      for (int s = 0; s < 15; ++s) f << " " << (s == 2 ? 0.0 : (s == 4 ? 600.0 : 250.0 + s)) << " 0.1 1.4";
      f << "\n";
    }
    k4siwecal::PedestalGrid p;
    std::string err;
    check(p.read(path, err), "pedestal table reads");
    check(near(p.get(3, 5, 0, 7), 250.0, 1e-9) && near(p.get(3, 5, 14, 7), 264.0, 1e-9), "per-SCA means placed");
    check(p.get(3, 5, 2, 7) == 0.0 && p.get(3, 5, 4, 7) == 0.0, "mean 0 and mean >= 500 are dropped");
    check(p.get(0, 0, 0, 0) == 0.0, "absent channel is 0");
    std::remove(path.c_str());
  }

  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "OK", failures);
  return failures ? 1 : 0;
}
