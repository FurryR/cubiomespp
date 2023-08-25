#include "./noise.hpp"
namespace cubiomes {
/**
 * @brief LAYER_INIT_SHA from original.
 */
inline constexpr auto LAYER_INIT_SHA = (~0ULL);
struct layer;
using mapfunc_t = int(const layer *, int *, int, int, int,
                      int);  // TODO: wtf?
                             /**
                              * @brief Defines internal ID of Minecraft versions.
                              * MC_1_X refers to the latest patch of the respective 1.X release.
                              * @note Development effort focuses on just the newest patch for each major
                              * release. Minor releases and major versions <= 1.0 are experimental.
                              */
enum version {
  MC_UNDEF,  // unspecified.
  MC_B1_7,
  MC_B1_8,
  MC_1_0_0,
  MC_1_0 = MC_1_0_0,
  MC_1_1_0,
  MC_1_1 = MC_1_1_0,
  MC_1_2_5,
  MC_1_2 = MC_1_2_5,
  MC_1_3_2,
  MC_1_3 = MC_1_3_2,
  MC_1_4_7,
  MC_1_4 = MC_1_4_7,
  MC_1_5_2,
  MC_1_5 = MC_1_5_2,
  MC_1_6_4,
  MC_1_6 = MC_1_6_4,
  MC_1_7_10,
  MC_1_7 = MC_1_7_10,
  MC_1_8_9,
  MC_1_8 = MC_1_8_9,
  MC_1_9_4,
  MC_1_9 = MC_1_9_4,
  MC_1_10_2,
  MC_1_10 = MC_1_10_2,
  MC_1_11_2,
  MC_1_11 = MC_1_11_2,
  MC_1_12_2,
  MC_1_12 = MC_1_12_2,
  MC_1_13_2,
  MC_1_13 = MC_1_13_2,
  MC_1_14_4,
  MC_1_14 = MC_1_14_4,
  MC_1_15_2,
  MC_1_15 = MC_1_15_2,
  MC_1_16_1,
  MC_1_16_5,
  MC_1_16 = MC_1_16_5,
  MC_1_17_1,
  MC_1_17 = MC_1_17_1,
  MC_1_18_2,
  MC_1_18 = MC_1_18_2,
  MC_1_19_2,
  MC_1_19,  // 1.19.3 - 1.19.4
  MC_1_20,
  MC_NEWEST = MC_1_20,
};
/**
 * @brief Defines internal ID of Minecraft dimensions.
 */
enum dimension {
  DIM_NETHER = -1,    // nether.
  DIM_OVERWORLD = 0,  // overworld.
  DIM_END = +1,       // the end.
  DIM_UNDEF = 1000,   // unspecified.
};
/**
 * @brief Defines internal ID of Minecrart biomes.
 * @note Extra unnecessary aliases are now removed. Please see the original
 * code for more information.
 * @see https://github.com/Cubitect/cubiomes/blob/master/layers.h
 */
enum biome_id {
  none = -1,
  // 0
  ocean = 0,
  plains,
  desert,
  mountains,
  forest,
  taiga,
  swamp,
  river,
  nether_wastes,
  the_end,
  // 10
  frozen_ocean,
  frozen_river,
  snowy_tundra,
  snowy_mountains,
  mushroom_fields,
  mushroom_field_shore,
  beach,
  desert_hills,
  wooded_hills,
  taiga_hills,
  // 20
  mountain_edge,
  jungle,
  jungle_hills,
  jungle_edge,
  deep_ocean,
  stone_shore,
  snowy_beach,
  birch_forest,
  birch_forest_hills,
  dark_forest,
  // 30
  snowy_taiga,
  snowy_taiga_hills,
  giant_tree_taiga,
  giant_tree_taiga_hills,
  wooded_mountains,
  savanna,
  savanna_plateau,
  badlands,
  wooded_badlands_plateau,
  badlands_plateau,
  // 40  --  1.13
  small_end_islands,
  end_midlands,
  end_highlands,
  end_barrens,
  warm_ocean,
  lukewarm_ocean,
  cold_ocean,
  deep_warm_ocean,
  deep_lukewarm_ocean,
  deep_cold_ocean,
  // 50
  deep_frozen_ocean,
  // Alpha 1.2 - Beta 1.7
  seasonal_forest,
  rainforest,
  shrubland,

  the_void = 127,

  // mutated variants
  sunflower_plains = plains + 128,
  desert_lakes = desert + 128,
  gravelly_mountains = mountains + 128,
  flower_forest = forest + 128,
  taiga_mountains = taiga + 128,
  swamp_hills = swamp + 128,
  ice_spikes = snowy_tundra + 128,
  modified_jungle = jungle + 128,
  modified_jungle_edge = jungle_edge + 128,
  tall_birch_forest = birch_forest + 128,
  tall_birch_hills = birch_forest_hills + 128,
  dark_forest_hills = dark_forest + 128,
  snowy_taiga_mountains = snowy_taiga + 128,
  giant_spruce_taiga = giant_tree_taiga + 128,
  giant_spruce_taiga_hills = giant_tree_taiga_hills + 128,
  modified_gravelly_mountains = wooded_mountains + 128,
  shattered_savanna = savanna + 128,
  shattered_savanna_plateau = savanna_plateau + 128,
  eroded_badlands = badlands + 128,
  modified_wooded_badlands_plateau = wooded_badlands_plateau + 128,
  modified_badlands_plateau = badlands_plateau + 128,
  // 1.14
  bamboo_jungle = 168,
  bamboo_jungle_hills = 169,
  // 1.16
  soul_sand_valley = 170,
  crimson_forest = 171,
  warped_forest = 172,
  basalt_deltas = 173,
  // 1.17
  dripstone_caves = 174,
  lush_caves = 175,
  // 1.18
  meadow = 177,
  grove = 178,
  snowy_slopes = 179,
  jagged_peaks = 180,
  frozen_peaks = 181,
  stony_peaks = 182,
  old_growth_birch_forest = tall_birch_forest,
  old_growth_pine_taiga = giant_tree_taiga,
  old_growth_spruce_taiga = giant_spruce_taiga,
  snowy_plains = snowy_tundra,
  sparse_jungle = jungle_edge,
  stony_shore = stone_shore,
  windswept_hills = mountains,
  windswept_forest = wooded_mountains,
  windswept_gravelly_hills = gravelly_mountains,
  windswept_savanna = shattered_savanna,
  wooded_badlands = wooded_badlands_plateau,
  // 1.19
  deep_dark = 183,
  mangrove_swamp = 184,
  // 1.20
  cherry_grove = 185,
};
/**
 * @brief Defines internal ID of the biomes.
 */
enum biome_category { oceanic, warm, lush, cold, freezing, special };

/**
 * @brief Enumeration of the layer indices in the layer stack.
 */
enum layer_id {
  L_CONTINENT_4096 = 0,
  L_ISLAND_4096 = L_CONTINENT_4096,
  L_ZOOM_4096,  // b1.8
  L_LAND_4096,  // b1.8
  L_ZOOM_2048,
  L_LAND_2048,
  L_ADD_ISLAND_2048 = L_LAND_2048,
  L_ZOOM_1024,
  L_LAND_1024_A,
  L_ADD_ISLAND_1024A = L_LAND_1024_A,
  L_LAND_1024_B,
  L_ADD_ISLAND_1024B = L_LAND_1024_B,  // 1.7+
  L_LAND_1024_C,
  L_ADD_ISLAND_1024C = L_LAND_1024_C,  // 1.7+
  L_ISLAND_1024,
  L_REMOVE_OCEAN_1024 = L_ISLAND_1024,  // 1.7+
  L_SNOW_1024,
  L_ADD_SNOW_1024 = L_SNOW_1024,
  L_LAND_1024_D,
  L_ADD_ISLAND_1024D = L_LAND_1024_D,  // 1.7+
  L_COOL_1024,
  L_COOL_WARM_1024 = L_COOL_1024,  // 1.7+
  L_HEAT_1024,
  L_HEAT_ICE_1024 = L_HEAT_1024,  // 1.7+
  L_SPECIAL_1024,                 // 1.7+
  L_ZOOM_512,
  L_LAND_512,  // 1.6-
  L_ZOOM_256,
  L_LAND_256,
  L_ADD_ISLAND_256 = L_LAND_256,
  L_MUSHROOM_256,
  L_ADD_MUSHROOM_256 = L_MUSHROOM_256,
  L_DEEP_OCEAN_256,  // 1.7+
  L_BIOME_256,
  L_BAMBOO_256,
  L14_BAMBOO_256 = L_BAMBOO_256,  // 1.14+
  L_ZOOM_128,
  L_ZOOM_64,
  L_BIOME_EDGE_64,
  L_NOISE_256,
  L_RIVER_INIT_256 = L_NOISE_256,
  L_ZOOM_128_HILLS,
  L_ZOOM_64_HILLS,
  L_HILLS_64,
  L_SUNFLOWER_64,
  L_RARE_BIOME_64 = L_SUNFLOWER_64,  // 1.7+
  L_ZOOM_32,
  L_LAND_32,
  L_ADD_ISLAND_32 = L_LAND_32,
  L_ZOOM_16,
  L_SHORE_16,        // NOTE: in 1.0 this slot is scale 1:32
  L_SWAMP_RIVER_16,  // 1.6-
  L_ZOOM_8,
  L_ZOOM_4,
  L_SMOOTH_4,
  L_ZOOM_128_RIVER,
  L_ZOOM_64_RIVER,
  L_ZOOM_32_RIVER,
  L_ZOOM_16_RIVER,
  L_ZOOM_8_RIVER,
  L_ZOOM_4_RIVER,
  L_RIVER_4,
  L_SMOOTH_4_RIVER,
  L_RIVER_MIX_4,
  L_OCEAN_TEMP_256,
  L13_OCEAN_TEMP_256 = L_OCEAN_TEMP_256,  // 1.13+
  L_ZOOM_128_OCEAN,
  L13_ZOOM_128 = L_ZOOM_128_OCEAN,  // 1.13+
  L_ZOOM_64_OCEAN,
  L13_ZOOM_64 = L_ZOOM_64_OCEAN,  // 1.13+
  L_ZOOM_32_OCEAN,
  L13_ZOOM_32 = L_ZOOM_32_OCEAN,  // 1.13+
  L_ZOOM_16_OCEAN,
  L13_ZOOM_16 = L_ZOOM_16_OCEAN,  // 1.13+
  L_ZOOM_8_OCEAN,
  L13_ZOOM_8 = L_ZOOM_8_OCEAN,  // 1.13+
  L_ZOOM_4_OCEAN,
  L13_ZOOM_4 = L_ZOOM_4_OCEAN,  // 1.13+
  L_OCEAN_MIX_4,
  L13_OCEAN_MIX_4 = L_OCEAN_MIX_4,  // 1.13+

  L_VORONOI_1,
  L_VORONOI_ZOOM_1 = L_VORONOI_1,

  // largeBiomes layers
  L_ZOOM_LARGE_A,
  L_ZOOM_LARGE_B,
  L_ZOOM_L_RIVER_A,
  L_ZOOM_L_RIVER_B,

  L_NUM
};
struct range {
  // Cuboidal range, given by a position, size and scaling in the horizontal
  // axes, used to define a generation range. The parameters for the vertical
  // control can be left at zero when dealing with versions without 3D volume
  // support. The vertical scaling is equal to 1:1 iff scale == 1, and 1:4
  // (default biome scale) in all other cases!
  //
  // @scale:  Horizontal scale factor, should be one of 1, 4, 16, 64, or 256
  //          additionally a value of zero bypasses scaling and expects a
  //          manual generation entry layer.
  // @x,z:    Horizontal position, i.e. coordinates of north-west corner.
  // @sx,sz:  Horizontal size (width and height for 2D), should be positive.
  // @y       Vertical position, 1:1 iff scale==1, 1:4 otherwise.
  // @sy      Vertical size. Values <= 0 are treated equivalent to 1.
  //
  // Volumes generated with a range are generally indexed as:
  //  out [ i_y*sx*sz + i_z*sx + i_x ]
  // where i_x, i_y, i_z are indecies in their respective directions.
  //
  // EXAMPLES
  // Area at normal biome scale (1:4):
  //  Range r_2d = {4, x,z, sx,sz};
  // (C99 syntax allows ommission of the trailing zero-initialization.)
  //
  // Area at block scale (1:1) at sea level:
  //  Range r_surf = {1, x,z, sx,sz, 63};
  // (Block level scale uses voronoi sampling with 1:1 vertical scaling.)
  //
  // Area at chunk scale (1:16) near sea level:
  //  Range r_surf16 = {16, x,z, sx,sz, 15};
  // (Note that the vertical scaling is always 1:4 for non-voronoi scales.)
  //
  // Volume at scale (1:4):
  //  Range r_vol = {4, x,z, sx,sz, y,sy};

  int scale;
  int x, z, sx, sz;
  int y, sy;
};
struct layer {
  layer(std::uint64_t worldSeed) {
    if (p2 != nullptr) *p2 = layer(worldSeed);

    if (p != nullptr) *p = layer(worldSeed);

    if (noise != nullptr) {
      random s(worldSeed);
      *((perlin_noise*)noise) = perlin_noise(s);
    }

    std::uint64_t ls = layerSalt;
    if (ls == 0) {  // Pre 1.13 the Hills branch stays zero-initialized
      startSalt = 0;
      startSeed = 0;
    } else if (ls == LAYER_INIT_SHA) {  // Post 1.14 Voronoi uses SHA256 for
                                        // initialization
      startSalt = getVoronoiSHA(worldSeed);
      startSeed = 0;
    } else {
      uint64_t st = worldSeed;
      st = seed_helper::step_seed(st, ls);
      st = seed_helper::step_seed(st, ls);
      st = seed_helper::step_seed(st, ls);

      startSalt = st;
      startSeed = seed_helper::step_seed(st, 0);
    }
  }
  mapfunc_t *getMap;

  int8_t mc;    // minecraft version
  int8_t zoom;  // zoom factor of layer
  int8_t edge;  // maximum border required from parent layer
  int scale;    // scale of this layer (cell = scale x scale blocks)

  uint64_t layerSalt;  // processed salt or initialization mode
  uint64_t startSalt;  // (depends on world seed) used to step PRNG forward
  uint64_t startSeed;  // (depends on world seed) start for chunk seeds

  void *noise;  // (depends on world seed) noise map data
  void *data;   // generic data for custom layers

  layer *p, *p2;  // parent layers
};

// Overworld biome generator up to 1.17
struct layer_stack {
  layer layers[L_NUM];
  layer *entry_1;  // entry scale (1:1) [L_VORONOI_1]
  layer *entry_4;  // entry scale (1:4) [L_RIVER_MIX_4|L_OCEAN_MIX_4]
  // unofficial entries for other scales (latest sensible layers):
  layer *entry_16;   // [L_SWAMP_RIVER_16|L_SHORE_16]
  layer *entry_64;   // [L_HILLS_64|L_SUNFLOWER_64]
  layer *entry_256;  // [L_BIOME_256|L_BAMBOO_256]
  perlin_noise oceanRnd;
};

// Nether biome generator 1.16+
struct nether_noise {  // altitude and wierdness don't affect nether biomes
  // and the weight is a 5th noise parameter which is constant
  void setNetherSeed(uint64_t seed) {
    random s(seed);
    temperature = double_perlin_noise(s, -7, 2);
    s = random(seed + 1);
    humidity = double_perlin_noise(s, -7, 2);
  }
  double_perlin_noise temperature;
  double_perlin_noise humidity;
};

// End biome generator 1.9+
struct end_noise {
  perlin_noise perlin;
  int mc;
};

struct surface_noise {
  double sample_surface_noise(int x, int y, int z) const {
    double xzScale = 684.412 * xzScale;
    double yScale = 684.412 * yScale;
    double xzStep = xzScale / xzFactor;
    double yStep = yScale / yFactor;

    double minNoise = 0;
    double maxNoise = 0;
    double mainNoise = 0;
    double persist = 1.0;
    double dx, dy, dz, sy, ty;

    for (std::size_t i = 0; i < 16; i++) {
      dx = math::maintain_precision(x * xzScale * persist);
      dy = math::maintain_precision(y * yScale * persist);
      dz = math::maintain_precision(z * xzScale * persist);
      sy = yScale * persist;
      ty = y * sy;

      octmin.octaves[i].sample_perlin(dx, dy, dz, sy, ty) / persist;
      octmax.octaves[i].sample_perlin(dx, dy, dz, sy, ty) / persist;

      if (i < 8) {
        dx = math::maintain_precision(x * xzStep * persist);
        dy = math::maintain_precision(y * yStep * persist);
        dz = math::maintain_precision(z * xzStep * persist);
        sy = yStep * persist;
        ty = y * sy;
        mainNoise +=
            octmain.octaves[i].sample_perlin(dx, dy, dz, sy, ty) / persist;
      }
      persist /= 2.0;
    }

    return math::clamped_lerp(0.5 + 0.05 * mainNoise, minNoise / 512.0,
                              maxNoise / 512.0);
  }
  surface_noise(int dim, uint64_t seed) {
    random s(seed);
    octmin = octave_noise(s, -15, 16);
    octmax = octave_noise(s, -15, 16);
    octmain = octave_noise(s, -7, 8);
    if (dim == DIM_END) {
      xzScale = 2.0;
      yScale = 1.0;
      xzFactor = 80;
      yFactor = 160;
    } else {
      // DIM_OVERWORLD
      octsurf = octave_noise(s, -3, 4);
      s.skip(262 * 10);
      octdepth = octave_noise(s, -15, 16);
      xzScale = 0.9999999814507745;
      yScale = 0.9999999814507745;
      xzFactor = 80;
      yFactor = 160;
    }
  }
  double xzScale, yScale;
  double xzFactor, yFactor;
  octave_noise octmin;
  octave_noise octmax;
  octave_noise octmain;
  octave_noise octsurf;
  octave_noise octdepth;
};

struct surface_noise_beta {
  surface_noise_beta(uint64_t seed) {
    random s(seed);
    octmin = octave_noise(s, 16, 684.412, 0.5, 1.0, 2.0);
    octmax = octave_noise(s, 16, 684.412, 0.5, 1.0, 2.0);
    octmain = octave_noise(s, 8, 684.412 / 80.0, 0.5, 1.0, 2.0);
    s.skip(262 * 8);
    octcontA = octave_noise(s, 10, 1.121, 0.5, 1.0, 2.0);
    octcontB = octave_noise(s, 16, 200.0, 0.5, 1.0, 2.0);
  }
  octave_noise octmin;
  octave_noise octmax;
  octave_noise octmain;
  octave_noise octcontA;
  octave_noise octcontB;
};

struct sea_level_column_noise_beta {
  double contASample;
  double contBSample;
  double minSample[2];
  double maxSample[2];
  double mainSample[2];
};

struct spline {
  int len, typ;
  float loc[12];
  float der[12];
  spline *val[12];
};

struct fix_spline {
  int len;
  float val;
};

struct spline_stack {  // the stack size here is just sufficient for overworld
                       // generation
  spline stack[42];
  fix_spline fstack[151];
  int len, flen;
};

enum {
  NP_TEMPERATURE = 0,
  NP_HUMIDITY = 1,
  NP_CONTINENTALNESS = 2,
  NP_EROSION = 3,
  NP_SHIFT = 4,
  NP_DEPTH = NP_SHIFT,  // not a real climate
  NP_WEIRDNESS = 5,
  NP_MAX
};
// Overworld biome generator for 1.18+
struct BiomeNoise {
  double_perlin_noise climate[NP_MAX];
  perlin_noise oct[2 * 23];  // buffer for octaves in double perlin noise
  spline *sp;
  spline_stack ss;
  int nptype;
  int mc;
};
// Overworld biome generator for pre-Beta 1.8
struct BiomeNoiseBeta {
  octave_noise climate[3];
  perlin_noise oct[10];
  int nptype;
  int mc;
};
class minecraft {
  /**
   * @brief Test if the biome exists in specify Minecraft version.
   *
   * @param id Biome ID.
   * @return true The biome exists in the specified version.
   * @return false The biome does not exist in the specified version.
   */
  constexpr bool biome_exists(biome_id id) const noexcept {
    if (ver >= MC_1_18) {
      if (id >= soul_sand_valley && id <= basalt_deltas) return true;
      if (id >= small_end_islands && id <= end_barrens) return true;

      if (id == cherry_grove) return ver >= MC_1_20;

      if (id == deep_dark || id == mangrove_swamp) return ver >= MC_1_19_2;

      switch (id) {
        case ocean:
        case plains:
        case desert:
        case mountains:  // windswept_hills
        case forest:
        case taiga:
        case swamp:
        case river:
        case nether_wastes:
        case the_end:
        case frozen_ocean:
        case frozen_river:
        case snowy_tundra:  // snowy_plains
        case mushroom_fields:
        case beach:
        case jungle:
        case jungle_edge:  // sparse_jungle
        case deep_ocean:
        case stone_shore:  // stony_shore
        case snowy_beach:
        case birch_forest:
        case dark_forest:
        case snowy_taiga:
        case giant_tree_taiga:  // old_growth_pine_taiga
        case wooded_mountains:  // windswept_forest
        case savanna:
        case savanna_plateau:
        case badlands:
        case wooded_badlands_plateau:  // wooded_badlands
        case warm_ocean:
        case lukewarm_ocean:
        case cold_ocean:
        case deep_warm_ocean:
        case deep_lukewarm_ocean:
        case deep_cold_ocean:
        case deep_frozen_ocean:
        case sunflower_plains:
        case gravelly_mountains:  // windswept_gravelly_hills
        case flower_forest:
        case ice_spikes:
        case tall_birch_forest:   // old_growth_birch_forest
        case giant_spruce_taiga:  // old_growth_spruce_taiga
        case shattered_savanna:   // windswept_savanna
        case eroded_badlands:
        case bamboo_jungle:
        case dripstone_caves:
        case lush_caves:
        case meadow:
        case grove:
        case snowy_slopes:
        case stony_peaks:
        case jagged_peaks:
        case frozen_peaks:
          return true;
        default:
          return false;
      }
    }

    if (ver <= MC_B1_7) {
      switch (id) {
        case plains:
        case desert:
        case forest:
        case taiga:
        case swamp:
        case snowy_tundra:
        case savanna:
        case seasonal_forest:
        case rainforest:
        case shrubland:
        // we treat areas below the sea level as oceans
        case ocean:
        case frozen_ocean:
          return true;
        default:
          return false;
      }
    }

    if (ver <= MC_B1_8) {
      switch (id) {
        case frozen_ocean:
        case frozen_river:
        case snowy_tundra:
        case mushroom_fields:
        case mushroom_field_shore:
        case the_end:
          return false;
      }
    }
    if (ver <= MC_1_0) {
      switch (id) {
        case snowy_mountains:
        case beach:
        case desert_hills:
        case wooded_hills:
        case taiga_hills:
        case mountain_edge:
          return false;
      }
    }

    if (id >= ocean && id <= mountain_edge) return true;
    if (id >= jungle && id <= jungle_hills) return ver >= MC_1_2;
    if (id >= jungle_edge && id <= badlands_plateau) return ver >= MC_1_7;
    if (id >= small_end_islands && id <= end_barrens) return ver >= MC_1_9;
    if (id >= warm_ocean && id <= deep_frozen_ocean) return ver >= MC_1_13;

    switch (id) {
      case the_void:
        return ver >= MC_1_9;
      case sunflower_plains:
      case desert_lakes:
      case gravelly_mountains:
      case flower_forest:
      case taiga_mountains:
      case swamp_hills:
      case ice_spikes:
      case modified_jungle:
      case modified_jungle_edge:
      case tall_birch_forest:
      case tall_birch_hills:
      case dark_forest_hills:
      case snowy_taiga_mountains:
      case giant_spruce_taiga:
      case giant_spruce_taiga_hills:
      case modified_gravelly_mountains:
      case shattered_savanna:
      case shattered_savanna_plateau:
      case eroded_badlands:
      case modified_wooded_badlands_plateau:
      case modified_badlands_plateau:
        return ver >= MC_1_7;
      case bamboo_jungle:
      case bamboo_jungle_hills:
        return ver >= MC_1_14;
      case soul_sand_valley:
      case crimson_forest:
      case warped_forest:
      case basalt_deltas:
        return ver >= MC_1_16_1;
      case dripstone_caves:
      case lush_caves:
        return ver >= MC_1_17;
      default:
        return false;
    }
  }
  /**
   * @brief Test if the biome is valid and belongs overworld in the specify
   * version.
   *
   * @param id Biome ID.
   * @return true The biome is valid and belongs overworld.
   * @return false The biome is invalid or does not belongs overworld.
   */
  bool is_overworld(biome_id id) const noexcept {
    if (!biome_exists(id)) return false;

    if (id >= small_end_islands && id <= end_barrens) return false;
    if (id >= soul_sand_valley && id <= basalt_deltas) return false;

    switch (id) {
      case nether_wastes:
      case the_end:
        return false;
      case frozen_ocean:
        return ver <= MC_1_6 || ver >= MC_1_13;
      case mountain_edge:
        return ver <= MC_1_6;
      case deep_warm_ocean:
      case the_void:
        return false;
      case tall_birch_forest:
        return ver <= MC_1_8 || ver >= MC_1_11;
      case dripstone_caves:
      case lush_caves:
        return ver >= MC_1_18;
    }
    return true;
  }
  dimension get_dimension(biome_id id) const noexcept {
    if (id >= small_end_islands && id <= end_barrens) return DIM_END;
    if (id >= soul_sand_valley && id <= basalt_deltas) return DIM_NETHER;
    if (id == the_end) return DIM_END;
    if (id == nether_wastes) return DIM_NETHER;
    return DIM_OVERWORLD;
  }

  biome_id get_mutated(int id) const noexcept {
    switch (id) {
      case plains:
        return sunflower_plains;
      case desert:
        return desert_lakes;
      case mountains:
        return gravelly_mountains;
      case forest:
        return flower_forest;
      case taiga:
        return taiga_mountains;
      case swamp:
        return swamp_hills;
      case snowy_tundra:
        return ice_spikes;
      case jungle:
        return modified_jungle;
      case jungle_edge:
        return modified_jungle_edge;
      // emulate MC-98995
      case birch_forest:
        return (ver >= MC_1_9 && ver <= MC_1_10) ? tall_birch_hills
                                                 : tall_birch_forest;
      case birch_forest_hills:
        return (ver >= MC_1_9 && ver <= MC_1_10) ? none : tall_birch_hills;
      case dark_forest:
        return dark_forest_hills;
      case snowy_taiga:
        return snowy_taiga_mountains;
      case giant_tree_taiga:
        return giant_spruce_taiga;
      case giant_tree_taiga_hills:
        return giant_spruce_taiga_hills;
      case wooded_mountains:
        return modified_gravelly_mountains;
      case savanna:
        return shattered_savanna;
      case savanna_plateau:
        return shattered_savanna_plateau;
      case badlands:
        return eroded_badlands;
      case wooded_badlands_plateau:
        return modified_wooded_badlands_plateau;
      case badlands_plateau:
        return modified_badlands_plateau;
      default:
        return none;
    }
  }

  biome_id get_category(biome_id id) const noexcept {
    // HACK: category type
    switch (id) {
      case beach:
      case snowy_beach:
        return beach;

      case desert:
      case desert_hills:
      case desert_lakes:
        return desert;

      case mountains:
      case mountain_edge:
      case wooded_mountains:
      case gravelly_mountains:
      case modified_gravelly_mountains:
        return mountains;

      case forest:
      case wooded_hills:
      case birch_forest:
      case birch_forest_hills:
      case dark_forest:
      case flower_forest:
      case tall_birch_forest:
      case tall_birch_hills:
      case dark_forest_hills:
        return forest;

      case snowy_tundra:
      case snowy_mountains:
      case ice_spikes:
        return snowy_tundra;

      case jungle:
      case jungle_hills:
      case jungle_edge:
      case modified_jungle:
      case modified_jungle_edge:
      case bamboo_jungle:
      case bamboo_jungle_hills:
        return jungle;

      case badlands:
      case eroded_badlands:
      case modified_wooded_badlands_plateau:
      case modified_badlands_plateau:
        return badlands;

      case wooded_badlands_plateau:
      case badlands_plateau:
        return ver <= MC_1_15 ? badlands : badlands_plateau;

      case mushroom_fields:
      case mushroom_field_shore:
        return mushroom_fields;

      case stone_shore:
        return stone_shore;

      case ocean:
      case frozen_ocean:
      case deep_ocean:
      case warm_ocean:
      case lukewarm_ocean:
      case cold_ocean:
      case deep_warm_ocean:
      case deep_lukewarm_ocean:
      case deep_cold_ocean:
      case deep_frozen_ocean:
        return ocean;

      case plains:
      case sunflower_plains:
        return plains;

      case river:
      case frozen_river:
        return river;

      case savanna:
      case savanna_plateau:
      case shattered_savanna:
      case shattered_savanna_plateau:
        return savanna;

      case swamp:
      case swamp_hills:
        return swamp;

      case taiga:
      case taiga_hills:
      case snowy_taiga:
      case snowy_taiga_hills:
      case giant_tree_taiga:
      case giant_tree_taiga_hills:
      case taiga_mountains:
      case snowy_taiga_mountains:
      case giant_spruce_taiga:
      case giant_spruce_taiga_hills:
        return taiga;

      case nether_wastes:
      case soul_sand_valley:
      case crimson_forest:
      case warped_forest:
      case basalt_deltas:
        return nether_wastes;

      default:
        return none;
    }
  }

  bool areSimilar(biome_id id1, biome_id id2) const noexcept {
    if (id1 == id2) return true;

    if (ver <= MC_1_15) {
      if (id1 == wooded_badlands_plateau || id1 == badlands_plateau)
        return true;
    }

    return get_category(id1) == get_category(id2);
  }

  bool isMesa(biome_id id) const noexcept {
    switch (id) {
      case badlands:
      case eroded_badlands:
      case modified_wooded_badlands_plateau:
      case modified_badlands_plateau:
      case wooded_badlands_plateau:
      case badlands_plateau:
        return true;
      default:
        return false;
    }
  }

  bool is_shallow_ocean(biome_id id) const noexcept {
    const std::uint64_t shallow_bits =
        (1ULL << ocean) | (1ULL << frozen_ocean) | (1ULL << warm_ocean) |
        (1ULL << lukewarm_ocean) | (1ULL << cold_ocean);
    return (std::uint32_t)id < 64 && ((1ULL << id) & shallow_bits);
  }

  bool is_deep_ocean(biome_id id) const noexcept {
    const uint64_t deep_bits =
        (1ULL << deep_ocean) | (1ULL << deep_warm_ocean) |
        (1ULL << deep_lukewarm_ocean) | (1ULL << deep_cold_ocean) |
        (1ULL << deep_frozen_ocean);
    return (uint32_t)id < 64 && ((1ULL << id) & deep_bits);
  }

  bool is_oceanic(biome_id id) const noexcept {
    const std::uint64_t ocean_bits =
        (1ULL << ocean) | (1ULL << frozen_ocean) | (1ULL << warm_ocean) |
        (1ULL << lukewarm_ocean) | (1ULL << cold_ocean) | (1ULL << deep_ocean) |
        (1ULL << deep_warm_ocean) | (1ULL << deep_lukewarm_ocean) |
        (1ULL << deep_cold_ocean) | (1ULL << deep_frozen_ocean);
    return (std::uint32_t)id < 64 && ((1ULL << id) & ocean_bits);
  }

  bool is_snowy(biome_id id) const noexcept {
    switch (id) {
      case frozen_ocean:
      case frozen_river:
      case snowy_tundra:
      case snowy_mountains:
      case snowy_beach:
      case snowy_taiga:
      case snowy_taiga_hills:
      case ice_spikes:
      case snowy_taiga_mountains:
        return true;
      default:
        return false;
    }
  }

 private:
  version ver;
};
}  // namespace cubiomes