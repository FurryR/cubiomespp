#ifndef RNG_HPP_
#define RNG_HPP_

#include <cinttypes>
#include <cstdint>
#include <cstdlib>

namespace cubiomes {
///=============================================================================
///                      Compiler and Platform Features
///=============================================================================

constexpr inline std::uint32_t BSWAP32(std::uint32_t x) {
  return ((x & 0x000000ff) << 24) | ((x & 0x0000ff00) << 8) |
         ((x & 0x00ff0000) >> 8) | ((x & 0xff000000) >> 24);
}

///=============================================================================
///                    C++ implementation of Java Random
///=============================================================================

struct random {
  random(std::uint64_t value)
      : seed((value ^ 0x5deece66d) & ((1ULL << 48) - 1)) {}
  inline constexpr int next(int bits) {
    seed = (seed * 0x5deece66d + 0xb) & ((1ULL << 48) - 1);
    return (int)((std::int64_t)seed >> (48 - bits));
  }
  inline constexpr int next_int(int n) {
    int bits, val;
    const int m = n - 1;

    if ((m & n) == 0) {
      std::uint64_t x = n * (std::uint64_t)next(31);
      return (int)((std::int64_t)x >> 31);
    }

    do {
      bits = next(31);
      val = bits % n;
    } while (bits - val + m < 0);
    return val;
  }
  inline constexpr std::uint64_t next_long() {
    return ((std::uint64_t)next(32) << 32) + next(32);
  }
  inline constexpr float next_float() { return next(24) / (float)(1 << 24); }
  inline constexpr double next_double() {
    std::uint64_t x = (std::uint64_t)next(26);
    x <<= 27;
    x += next(27);
    return (std::int64_t)x / (double)(1ULL << 53);
  }
  inline constexpr int next_int24() {
    std::uint64_t a = (1ULL << 48) - 1;
    std::uint64_t c = 0x5deece66dULL * seed;
    c += 11;
    a &= c;
    seed = a;
    a = (std::uint64_t)((std::int64_t)a >> 17);
    c = 0xaaaaaaab * a;
    c = (std::uint64_t)((std::int64_t)c >> 36);
    return (int)a - (int)(c << 3) * 3;
  }
  /* Jumps forwards in the random number sequence by simulating 'n' calls to
   * next.
   */
  inline constexpr int skip(std::size_t n) {
    std::uint64_t m = 1;
    std::uint64_t a = 0;
    std::uint64_t im = 0x5deece66dULL;
    std::uint64_t ia = 0xb;
    std::uint64_t k;

    for (k = n; k; k >>= 1) {
      if (k & 1) {
        m *= im;
        a = im * a + ia;
      }
      ia = (im + 1) * ia;
      im *= im;
    }

    seed = seed * m + a;
    seed &= 0xffffffffffffULL;
  }

 private:
  std::uint64_t seed;
};

///=============================================================================
///                               Xoroshiro 128
///=============================================================================

struct xoroshiro {
  constexpr inline std::uint64_t next_long() {
    std::uint64_t l = lo;
    std::uint64_t h = hi;
    std::uint64_t n = rotl64(l + h, 17) + l;
    h ^= l;
    lo = rotl64(l, 49) ^ h ^ (h << 21);
    hi = rotl64(h, 28);
    return n;
  }
  constexpr inline int next_int(std::uint32_t n) {
    std::uint64_t r = (next_long() & 0xFFFFFFFF) * n;
    if ((std::uint32_t)r < n) {
      while ((std::uint32_t)r < (~n + 1) % n) {
        r = (next_long() & 0xFFFFFFFF) * n;
      }
    }
    return r >> 32;
  }
  constexpr inline double next_double() {
    return (next_long() >> (64 - 53)) * 1.1102230246251565E-16;
  }
  constexpr inline float next_float() {
    return (next_long() >> (64 - 24)) * 5.9604645E-8F;
  }
  constexpr inline void skip(std::size_t count) {
    while (count-- > 0) next_long();
  }
  constexpr inline std::uint64_t next_long_j() {
    std::uint32_t a = next_long() >> 32;  // FIXME: std::int32_t?
    std::uint32_t b = next_long() >> 32;
    return ((std::uint64_t)a << 32) + b;
  }
  constexpr inline int next_int_j(std::uint32_t n) {
    int bits, val;
    const int m = n - 1;

    if ((m & n) == 0) {
      std::uint64_t x = n * (next_long() >> 33);
      return (int)((std::int64_t)x >> 31);
    }

    do {
      bits = (next_long() >> 33);
      val = bits % n;
    } while (bits - val + m < 0);
    return val;
  }
  xoroshiro() = default;
  xoroshiro(std::uint64_t seed) {
    constexpr std::uint64_t XL = 0x9e3779b97f4a7c15ULL;
    constexpr std::uint64_t XH = 0x6a09e667f3bcc909ULL;
    constexpr std::uint64_t A = 0xbf58476d1ce4e5b9ULL;
    constexpr std::uint64_t B = 0x94d049bb133111ebULL;
    std::uint64_t l = seed ^ XH;
    std::uint64_t h = l + XL;
    l = (l ^ (l >> 30)) * A;
    h = (h ^ (h >> 30)) * A;
    l = (l ^ (l >> 27)) * B;
    h = (h ^ (h >> 27)) * B;
    l = l ^ (l >> 31);
    h = h ^ (h >> 31);
    lo = l;
    hi = h;
  }
  std::uint64_t lo, hi;

 private:
  /// imitate amd64/x64 rotate instructions
  static constexpr inline std::uint64_t rotl64(std::uint64_t x,
                                               std::uint8_t b) {
    return (x << b) | (x >> (64 - b));
  }

  static constexpr inline std::uint32_t rotr32(std::uint32_t a,
                                               std::uint8_t b) {
    return (a >> b) | (a << (32 - b));
  }
};

//==============================================================================
//                              MC Seed Helpers
//==============================================================================

/**
 * The seed pipeline:
 *
 * get_layer_salt(n)                -> layerSalt (layer_salt)
 * layerSalt (layer_salt), worldSeed (world_seed) -> startSalt (st), startSeed
 * (ss) startSeed (ss), coords (x,z)   -> chunkSeed (cs)
 *
 * The chunkSeed alone is enough to generate the first PRNG integer with:
 *   first_prng(cs, mod)
 * subsequent PRNG integers are generated by stepping the chunkSeed forwards,
 * salted with startSalt:
 *   cs_next = step_seed(cs, st)
 */
namespace seed_helper {
/**
 * @brief The basic function of calcuating seeds.
 *
 * @param seed seed.
 * @param salt seed salt.
 * @return constexpr std::uint64_t the final seed.
 */
constexpr inline std::uint64_t step_seed(std::uint64_t seed,
                                         std::uint64_t salt) {
  return seed * (seed * 6364136223846793005ULL + 1442695040888963407ULL) + salt;
}
/**
 * @brief get the first PRNG integer.
 *
 * @param seed chunk seed.
 * @param mod FIXME: add doc later
 * @return constexpr int the first PRNG integer.
 */
constexpr inline int first_prng(std::uint64_t seed, int mod) {
  int ret = (int)(((std::int64_t)seed >> 24) % mod);
  if (ret < 0) ret += mod;
  return ret;
}
/**
 * @brief Checks if the first PRNG integer is zero.
 *
 * @param seed chunk seed.
 * @param mod FIXME: add doc later
 * @return true The PRNG is zero.
 * @return false The PRNG is not zero.
 */
constexpr inline bool check_prng_zero(std::uint64_t seed, int mod) {
  return (int)(((std::int64_t)seed >> 24) % mod) == 0;
}
/**
 * @brief Get the seed of specify chunk.
 *
 * @param ss Start seed.
 * @param x The chunk's X position FIXME: fix this doc
 * @param z The chunk's Z position FIXME: fix this doc
 * @return constexpr std::uint64_t The seed of the chunk.
 */
constexpr inline std::uint64_t get_chunk_seed(std::uint64_t start_seed, int x,
                                              int z) {
  std::uint64_t cs = start_seed + x;  // start se
  cs = step_seed(cs, z);
  cs = step_seed(cs, x);
  cs = step_seed(cs, z);
  return cs;
}
/**
 * @brief get layer salt.
 *
 * @param salt Layer. FIXME: fix this doc later
 * @return constexpr std::uint64_t
 */
constexpr inline std::uint64_t get_layer_salt(std::uint64_t salt) {
  std::uint64_t layer_salt = step_seed(salt, salt);
  layer_salt = step_seed(layer_salt, salt);
  layer_salt = step_seed(layer_salt, salt);
  return layer_salt;
}
/**
 * @brief Get start salt.
 *
 * @param world_seed world seed.
 * @param layer_salt layer salt from get_layer_salt.
 * @return constexpr std::uint64_t
 */
constexpr inline std::uint64_t get_start_salt(std::uint64_t world_seed,
                                              std::uint64_t layer_salt) {
  std::uint64_t st = world_seed;
  st = step_seed(st, layer_salt);
  st = step_seed(st, layer_salt);
  st = step_seed(st, layer_salt);
  return st;
}
/**
 * @brief Get the start seed.
 *
 * @param world_seed The world seed.
 * @param layer_salt Layer salt from get_layer_salt.
 * @return constexpr std::uint64_t the start seed.
 */
constexpr inline std::uint64_t getStartSeed(std::uint64_t world_seed,
                                            std::uint64_t layer_salt) {
  return step_seed(get_start_salt(world_seed, layer_salt), 0);
}
}  // namespace seed_helper

///============================================================================
///                               Arithmatic
///============================================================================

namespace math {
// NOTE: no need to make doc.
/* Linear interpolations
 */
inline constexpr double lerp(double part, double from, double to) {
  return from + part * (to - from);
}

inline constexpr double lerp2(double dx, double dy, double v00, double v10,
                              double v01, double v11) {
  return lerp(dy, lerp(dx, v00, v10), lerp(dx, v01, v11));
}
// HACK: should pass 3d arguments using a class.
inline constexpr double lerp3(double dx, double dy, double dz, double v000,
                              double v100, double v010, double v110,
                              double v001, double v101, double v011,
                              double v111) {
  v000 = lerp2(dx, dy, v000, v100, v010, v110);
  v001 = lerp2(dx, dy, v001, v101, v011, v111);
  return lerp(dz, v000, v001);
}

inline constexpr double clamped_lerp(double part, double from, double to) {
  if (part <= 0) return from;
  if (part >= 1) return to;
  return lerp(part, from, to);
}

/* Find the modular inverse: (1/x) | mod m.
 * Assumes x and m are positive (less than 2^63), co-prime.
 */
inline constexpr std::uint64_t mul_inv(std::uint64_t x, std::uint64_t m) {
  std::uint64_t t, q, a, b, n;
  if ((std::int64_t)m <= 1) return 0;  // no solution

  n = m;
  a = 0;
  b = 1;

  while ((std::int64_t)x > 1) {
    if (m == 0) return 0;  // x and m are co-prime
    q = x / m;
    t = m;
    m = x % m;
    x = t;
    t = a;
    a = b - q * a;
    b = t;
  }

  if ((std::int64_t)b < 0) b += n;
  return b;
}
}  // namespace math
}  // namespace cubiomes

#endif /* RNG_HPP_ */