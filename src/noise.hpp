#ifndef NOISE_HPP_
#define NOISE_HPP_
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "./rng.hpp"

namespace cubiomes {
namespace math {
inline constexpr double maintain_precision(
    double x) {  // This is a highly performance critical
                 // function that is used to correct
  // progressing errors from float-maths. However, since cubiomes uses
  // doubles anyway, this seems useless in practice.

  // return x - round(x / 33554432.0) * 33554432.0;

  // in C++ version:
  // it should be a static inline constexpr function and it does now.
  return x;
}
inline constexpr double indexed_lerp(int idx, double a, double b, double c) {
  switch (idx & 0xf) {
    case 0:
      return a + b;
    case 1:
      return -a + b;
    case 2:
      return a - b;
    case 3:
      return -a - b;
    case 4:
      return a + c;
    case 5:
      return -a + c;
    case 6:
      return a - c;
    case 7:
      return -a - c;
    case 8:
      return b + c;
    case 9:
      return -b + c;
    case 10:
      return b - c;
    case 11:
      return -b - c;
    case 12:
      return a + b;
    case 13:
      return -b + c;
    case 14:
      return -a + b;
    case 15:
      return -b - c;
  }
  std::unreachable();
  return 0;
}
}  // namespace math
struct perlin_noise {
  double sample_perlin(double d1, double d2, double d3, double yamp,
                       double ymin) const {
    d1 += a;
    d2 += b;
    d3 += c;
    const std::uint8_t* idx = d;
    int i1 = (int)floor(d1);
    int i2 = (int)floor(d2);
    int i3 = (int)floor(d3);
    d1 -= i1;
    d2 -= i2;
    d3 -= i3;
    double t1 = d1 * d1 * d1 * (d1 * (d1 * 6.0 - 15.0) + 10.0);
    double t2 = d2 * d2 * d2 * (d2 * (d2 * 6.0 - 15.0) + 10.0);
    double t3 = d3 * d3 * d3 * (d3 * (d3 * 6.0 - 15.0) + 10.0);

    if (yamp) {
      double yclamp = ymin < d2 ? ymin : d2;
      d2 -= floor(yclamp / yamp) * yamp;
    }

    i1 &= 0xff;
    i2 &= 0xff;
    i3 &= 0xff;

    int a1 = idx[i1] + i2;
    int b1 = idx[i1 + 1] + i2;

    int a2 = idx[a1] + i3;
    int a3 = idx[a1 + 1] + i3;
    int b2 = idx[b1] + i3;
    int b3 = idx[b1 + 1] + i3;

    double l1 = math::indexed_lerp(idx[a2], d1, d2, d3);
    double l2 = math::indexed_lerp(idx[b2], d1 - 1, d2, d3);
    double l3 = math::indexed_lerp(idx[a3], d1, d2 - 1, d3);
    double l4 = math::indexed_lerp(idx[b3], d1 - 1, d2 - 1, d3);
    double l5 = math::indexed_lerp(idx[a2 + 1], d1, d2, d3 - 1);
    double l6 = math::indexed_lerp(idx[b2 + 1], d1 - 1, d2, d3 - 1);
    double l7 = math::indexed_lerp(idx[a3 + 1], d1, d2 - 1, d3 - 1);
    double l8 = math::indexed_lerp(idx[b3 + 1], d1 - 1, d2 - 1, d3 - 1);

    l1 = math::lerp(t1, l1, l2);
    l3 = math::lerp(t1, l3, l4);
    l5 = math::lerp(t1, l5, l6);
    l7 = math::lerp(t1, l7, l8);

    l1 = math::lerp(t2, l1, l3);
    l5 = math::lerp(t2, l5, l7);

    return math::lerp(t3, l1, l5);
  }
  void sample_perlin_beta17_terrain(double* v, double d1, double d3,
                                    double y_lac_amp) const {
    // y_lac_amp stands for y lacunarity amplitude
    int genFlag = -1;
    double l1 = 0;
    double l3 = 0;
    double l5 = 0;
    double l7 = 0;

    d1 += a;
    d3 += c;
    const uint8_t* idx = d;
    int i1 = (int)floor(d1);
    int i3 = (int)floor(d3);
    d1 -= i1;
    d3 -= i3;
    double t1 = d1 * d1 * d1 * (d1 * (d1 * 6.0 - 15.0) + 10.0);
    double t3 = d3 * d3 * d3 * (d3 * (d3 * 6.0 - 15.0) + 10.0);

    i1 &= 0xff;
    i3 &= 0xff;

    double d2;
    int i2, yi, yic, gfCopy;
    for (yi = 0; yi <= 7; yi++) {
      d2 = yi * lacunarity * y_lac_amp + b;
      i2 = ((int)floor(d2)) & 0xff;
      if (yi == 0 || i2 != genFlag) {
        yic = yi;
        gfCopy = genFlag;
        genFlag = i2;
      }
    }
    genFlag = gfCopy;

    double t2;
    for (yi = yic; yi <= 8; yi++) {
      d2 = yi * lacunarity * y_lac_amp + b;
      i2 = (int)floor(d2);
      d2 -= i2;
      t2 = d2 * d2 * d2 * (d2 * (d2 * 6.0 - 15.0) + 10.0);

      i2 &= 0xff;

      if (yi == 0 || i2 != genFlag) {
        genFlag = i2;
        int a1 = idx[i1] + i2;
        int b1 = idx[i1 + 1] + i2;

        int a2 = idx[a1] + i3;
        int a3 = idx[a1 + 1] + i3;
        int b2 = idx[b1] + i3;
        int b3 = idx[b1 + 1] + i3;

        double m1 = math::indexed_lerp(idx[a2], d1, d2, d3);
        double l2 = math::indexed_lerp(idx[b2], d1 - 1, d2, d3);
        double m3 = math::indexed_lerp(idx[a3], d1, d2 - 1, d3);
        double l4 = math::indexed_lerp(idx[b3], d1 - 1, d2 - 1, d3);
        double m5 = math::indexed_lerp(idx[a2 + 1], d1, d2, d3 - 1);
        double l6 = math::indexed_lerp(idx[b2 + 1], d1 - 1, d2, d3 - 1);
        double m7 = math::indexed_lerp(idx[a3 + 1], d1, d2 - 1, d3 - 1);
        double l8 = math::indexed_lerp(idx[b3 + 1], d1 - 1, d2 - 1, d3 - 1);

        l1 = math::lerp(t1, m1, l2);
        l3 = math::lerp(t1, m3, l4);
        l5 = math::lerp(t1, m5, l6);
        l7 = math::lerp(t1, m7, l8);
      }

      if (yi >= 7) {
        double n1 = math::lerp(t2, l1, l3);
        double n5 = math::lerp(t2, l5, l7);

        v[yi - 7] += math::lerp(t3, n1, n5) * amplitude;
      }
    }
  }

  double sample_simple_x2d(double x, double y) const {
    const double SKEW = 0.5 * (sqrt(3) - 1.0);
    const double UNSKEW = (3.0 - sqrt(3)) / 6.0;

    double hf = (x + y) * SKEW;
    int hx = (int)floor(x + hf);
    int hz = (int)floor(y + hf);
    double mhxz = (hx + hz) * UNSKEW;
    double x0 = x - (hx - mhxz);
    double y0 = y - (hz - mhxz);
    int offx = (x0 > y0);
    int offz = !offx;
    double x1 = x0 - offx + UNSKEW;
    double y1 = y0 - offz + UNSKEW;
    double x2 = x0 - 1.0 + 2.0 * UNSKEW;
    double y2 = y0 - 1.0 + 2.0 * UNSKEW;
    int gi0 = d[0xff & (hz)];
    int gi1 = d[0xff & (hz + offz)];
    int gi2 = d[0xff & (hz + 1)];
    gi0 = d[0xff & (gi0 + hx)];
    gi1 = d[0xff & (gi1 + hx + offx)];
    gi2 = d[0xff & (gi2 + hx + 1)];
    double t = 0;
    t += simple_xgrad(gi0 % 12, x0, y0, 0.0, 0.5);
    t += simple_xgrad(gi1 % 12, x1, y1, 0.0, 0.5);
    t += simple_xgrad(gi2 % 12, x2, y2, 0.0, 0.5);
    return 70.0 * t;
  }
  perlin_noise(random& seed) {
    int i = 0;
    // memset(noise, 0, sizeof(*noise));
    a = seed.next_double() * 256.0;
    b = seed.next_double() * 256.0;
    c = seed.next_double() * 256.0;
    amplitude = 1.0;
    lacunarity = 1.0;

    uint8_t* idx = d;
    for (i = 0; i < 256; i++) {
      idx[i] = i;
    }
    for (i = 0; i < 256; i++) {
      int j = seed.next_int(256 - i) + i;
      uint8_t n = idx[i];
      idx[i] = idx[j];
      idx[j] = n;
      idx[i + 256] = idx[i];
    }
  }
  perlin_noise(xoroshiro& xr) {
    int i = 0;
    a = xr.next_double() * 256.0;
    b = xr.next_double() * 256.0;
    c = xr.next_double() * 256.0;
    amplitude = 1.0;
    lacunarity = 1.0;

    std::uint8_t* idx = d;
    for (i = 0; i < 256; i++) {
      idx[i] = i;
    }
    for (i = 0; i < 256; i++) {
      int j = xr.next_int(256 - i) + i;  // FIXME: should be std::size_t.
      std::uint8_t n = idx[i];
      idx[i] = idx[j];
      idx[j] = n;
      idx[i + 256] = idx[i];
    }
  }
  std::uint8_t d[512];
  double a, b, c;
  double amplitude;
  double lacunarity;

 private:
  static double simple_xgrad(int idx, double x, double y, double z, double d) {
    double con = d - x * x - y * y - z * z;
    if (con < 0) return 0;
    con *= con;
    return con * con * math::indexed_lerp(idx, x, y, z);
  }
};
struct octave_noise {
  double sample_octave_amp(double x, double y, double z, double yamp,
                           double ymin, int ydefault) const {
    double v = 0;
    for (std::size_t i = 0; i < octaves.size(); i++) {
      const perlin_noise* p = &octaves[i];
      double lf = p->lacunarity;
      double ax = math::maintain_precision(x * lf);
      double ay = ydefault ? -p->b : math::maintain_precision(y * lf);
      double az = math::maintain_precision(z * lf);
      double pv = p->sample_perlin(ax, ay, az, yamp * lf, ymin * lf);
      v += p->amplitude * pv;
    }
    return v;
  }

  double sample_octave(double x, double y, double z) const {
    double v = 0;
    int i;
    for (i = 0; i < octaves.size(); i++) {
      const perlin_noise* p = &octaves[i];
      double lf = p->lacunarity;
      double ax = math::maintain_precision(x * lf);
      double ay = math::maintain_precision(y * lf);
      double az = math::maintain_precision(z * lf);
      double pv = p->sample_perlin(ax, ay, az, 0, 0);
      v += p->amplitude * pv;
    }
    return v;
  }

  double sample_octave_beta17_biome(double x, double z) const {
    double v = 0;
    int i;
    for (i = 0; i < octaves.size(); i++) {
      const perlin_noise* p = &octaves[i];
      double lf = p->lacunarity;
      double ax = math::maintain_precision(x * lf) + p->a;
      double az = math::maintain_precision(z * lf) + p->b;
      double pv = p->sample_simple_x2d(ax, az);
      v += p->amplitude * pv;
    }
    return v;
  }

  void sample_octave_beta17_terrain(double* v, double x, double z,
                                    int y_lac_flag, double lacmin) const {
    v[0] = 0.0;
    v[1] = 0.0;
    int i;
    for (i = 0; i < octaves.size(); i++) {
      const perlin_noise* p = &octaves[i];
      double lf = p->lacunarity;
      if (lacmin && lf > lacmin) continue;
      double ax = math::maintain_precision(x * lf);
      double az = math::maintain_precision(z * lf);
      p->sample_perlin_beta17_terrain(v, ax, az, y_lac_flag ? 0.5 : 1.0);
    }
  }
  octave_noise() = default;
  octave_noise(random& seed, int omin, std::size_t len) {
    int i;
    int end = omin + len - 1;
    double persist = 1.0 / ((1LL << len) - 1.0);
    double lacuna = pow(2.0, end);
    if (len < 1 || end > 0) {
      throw std::invalid_argument(
          "octavePerlinInit(): unsupported octave range\n");
      return;
    }
    octaves.reserve(len);

    if (end == 0) {
      octaves.emplace_back(seed);
      octaves[0].amplitude = persist;
      octaves[0].lacunarity = lacuna;
      persist *= 2.0;
      lacuna *= 0.5;
      i = 1;
    } else {
      seed.skip(-end * 262);
      i = 0;
    }

    for (; i < len; i++) {
      octaves.emplace_back(seed);
      octaves[i].amplitude = persist;
      octaves[i].lacunarity = lacuna;
      persist *= 2.0;
      lacuna *= 0.5;
    }
  }
  octave_noise(random& seed, std::size_t octcnt, double lac, double lacMul,
               double persist, double persistMul) {
    int i;
    octaves.reserve(octcnt);
    for (i = 0; i < octcnt; i++) {
      octaves.emplace_back(seed);
      octaves[i].amplitude = persist;
      octaves[i].lacunarity = lac;
      persist *= persistMul;
      lac *= lacMul;
    }
  }
  octave_noise(xoroshiro& xr, const double* amplitudes, int omin,
               std::size_t len, int nmax) {
    constexpr std::uint64_t md5_octave_n[][2] = {
        {0xb198de63a8012672, 0x7b84cad43ef7b5a8},  // md5 "octave_-12"
        {0x0fd787bfbc403ec3, 0x74a4a31ca21b48b8},  // md5 "octave_-11"
        {0x36d326eed40efeb2, 0x5be9ce18223c636a},  // md5 "octave_-10"
        {0x082fe255f8be6631, 0x4e96119e22dedc81},  // md5 "octave_-9"
        {0x0ef68ec68504005e, 0x48b6bf93a2789640},  // md5 "octave_-8"
        {0xf11268128982754f, 0x257a1d670430b0aa},  // md5 "octave_-7"
        {0xe51c98ce7d1de664, 0x5f9478a733040c45},  // md5 "octave_-6"
        {0x6d7b49e7e429850a, 0x2e3063c622a24777},  // md5 "octave_-5"
        {0xbd90d5377ba1b762, 0xc07317d419a7548d},  // md5 "octave_-4"
        {0x53d39c6752dac858, 0xbcd1c5a80ab65b3e},  // md5 "octave_-3"
        {0xb4a24d7a84e7677b, 0x023ff9668e89b5c4},  // md5 "octave_-2"
        {0xdffa22b534c5f608, 0xb9b67517d3665ca9},  // md5 "octave_-1"
        {0xd50708086cef4d7c, 0x6e1651ecc7f43309},  // md5 "octave_0"
    };
    constexpr double lacuna_ini[] = {
        // -omin = 3..12
        1,        .5,       .25,      1. / 8,    1. / 16,   1. / 32,   1. / 64,
        1. / 128, 1. / 256, 1. / 512, 1. / 1024, 1. / 2048, 1. / 4096,
    };
    constexpr double persist_ini[] = {
        // len = 4..9
        0,        1,        2. / 3,    4. / 7,     8. / 15,
        16. / 31, 32. / 63, 64. / 127, 128. / 255, 256. / 511,
    };
    // extra range assert (pass by exception, so it could not be skipped in C++
    // version)
    if (-omin < 0 || -omin >= (int)(sizeof(lacuna_ini) / sizeof(double)) ||
        len < 0 || len >= (int)(sizeof(persist_ini) / sizeof(double))) {
      throw std::invalid_argument(
          "Fatal: octave initialization out of range\n");
    }
    double lacuna = lacuna_ini[-omin];
    double persist = persist_ini[len];
    std::uint64_t xlo = xr.next_long();
    std::uint64_t xhi = xr.next_long();
    int i = 0, n = 0;
    octaves.reserve(len);
    for (; i < len && n != nmax; i++, lacuna *= 2.0, persist *= 0.5) {
      if (amplitudes[i] == 0) continue;
      xoroshiro pxr;
      pxr.lo = xlo ^ md5_octave_n[12 + omin + i][0];
      pxr.hi = xhi ^ md5_octave_n[12 + omin + i][1];
      octaves.emplace_back(pxr);
      octaves[n].amplitude = amplitudes[i] * persist;
      octaves[n].lacunarity = lacuna;
      n++;
    }
  }
  std::vector<perlin_noise> octaves;
};
struct double_perlin_noise {
  double sample_double_perlin(double x, double y, double z) const {
    const double f = 337.0 / 331.0;
    double v = 0;

    v += octA.sample_octave(x, y, z);
    v += octB.sample_octave(x * f, y * f, z * f);

    return v * amplitude;
  }
  double_perlin_noise(random& seed, int omin, std::size_t len)
      : amplitude((10.0 / 6.0) * len / (len + 1)),
        octA(seed, omin, len),
        octB(seed, omin, len) {  // require: len >= 1 && omin+len <= 0
    // NOTE: C++ version feature
    if (len < 1 || omin + len > 0)
      throw std::invalid_argument("double_perlin_noise() out of range");
  }

  /**
   * @brief Sets up a DoublePerlinNoise generator (MC 1.18+).
   * @param noise      Object to be initialized
   * @param xr         Xoroshiro random object
   * @param amplitudes Octave amplitude, needs at least one non-zero
   * @param omin       First octave
   * @param len        Length of amplitudes array
   * @param nmax       Number of octaves available in buffer (can be <=0 to
   * ignore)
   */
  double_perlin_noise(xoroshiro& xr, const double* amplitudes, int omin,
                      std::size_t len, int nmax) {
    // NOTE: should return: octA.octaves.size() + octB.octaves.size() (from C
    // version)
    int i, n = 0, na = -1, nb = -1;
    if (nmax > 0) {
      na = (nmax + 1) >> 1;
      nb = nmax - na;
    }
    octA = octave_noise(xr, amplitudes, omin, len, na);
    octB = octave_noise(xr, amplitudes, omin, len, nb);

    // trim amplitudes of zero
    for (i = len - 1; i >= 0 && amplitudes[i] == 0.0; i--) len--;
    for (i = 0; amplitudes[i] == 0.0; i++) len--;
    static const double amp_ini[] = {
        // (5 ./ 3) * len / (len + 1), len = 2..9
        0,        5. / 6,   10. / 9,  15. / 12, 20. / 15,
        25. / 18, 30. / 21, 35. / 24, 40. / 27, 45. / 30,
    };
    amplitude = amp_ini[len];
  }
  double amplitude;
  octave_noise octA;
  octave_noise octB;
};
}  // namespace cubiomes
#endif