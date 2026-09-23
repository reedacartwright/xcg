/*
# MIT License
#
# Copyright (c) 2026 Reed A. Cartwright <racartwright@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
*/

/*
## XCG: Expanded Congruential Generator

XCG is a family of random number generators inspired by and derived from PCG's
extended generation scheme <https://www.pcg-random.org/>. At its core, XCG is a
128-bit linear congruential generator that returns its high 64-bits. This
simple generator is enough to pass PractRand tests (other than TMFn which is
designed to detect LCGs.) For flexibility, the family provides generators that
vary in the sizes of their parameter spaces. For example, XCG-1280 has 1280
bits of parameter-space: 128 bits for the state of the LCG, 128 bits for the
increment of the LCG, and 1024 bits for a set of 64-bit salts that are mixed
into the output.

### Period

The period of XCG is 2^126 if using the MCG core generator or 2^128 if using
the LCG core generator. While the increment and salts expand the size of XCG's
parameter space, they do not extend XCG's period. The periods of the core
generators are big enough that no application will ever wrap around during
normal usage. This simplifies XCG's algorithm without sacrificing utility or
parameter-space flexibility.
*/

#ifndef XCG_RANDOM_H
#define XCG_RANDOM_H

#include <array>
#include <bit>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <ranges>
#include <span>
#include <sys/types.h>
#include <utility>

static_assert(__cplusplus >= 202002L, "Requires C++20");
static_assert(__SIZEOF_INT128__ == 16, "Requires __uint128_t");

namespace xcg {

using uint128_t = __uint128_t;

/*
This coefficient came from searching parameter space for the best 64-bit
multiplier that works for 128-bit LCG, 128-bit MCG, 64-bit LCG, and 64-bit MCG.
*/
constexpr uint64_t XCG_MULT = 0xf68a43306d8e0225U;
// 0b1111011010001010010000110011000001101101100011100000001000100101

template <bool USE_LCG_ = false, std::size_t SALT_N_ = 0> struct xcg;

/*
// Example simple implementations of XCG generators
struct xcg128_t {
  uint128_t state = 1U;
};

inline uint64_t get_u64(xcg128_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  return (u >> 64);
}

struct xcg256_t {
  uint128_t state = 0U;
  uint128_t inc = 1U;
};

inline uint64_t get_u64(xcg256_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  gen->state += gen->inc;
  return (u >> 64);
}

struct xcg320_t {
  uint128_t state = 0U;
  uint128_t inc = 1U;
  uint64_t salt = 0U;
};

inline uint64_t get_u64(xcg320_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  gen->state += gen->inc;
  return (u >> 64) ^ gen->salt;
}

struct xcg512_t {
  uint128_t state = 0U;
  uint128_t inc = 1U;
  uint64_t salts[4] = {};
};

inline uint64_t get_u64(xcg512_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  gen->state += gen->inc;
  uint64_t low = (uint64_t)u;
  return (u >> 64) ^ gen->salts[low >> 62];
}
*/

namespace detail {

// Magic number from the first 128-bits of the fractional part of Sqrt(3).
constexpr uint64_t MAGICB_64 = 0xbb67ae8584caa73b;
constexpr uint64_t MAGICB_64_LOW = 0x25742d7078b83b89;
constexpr uint128_t MAGICB_128 =
    (static_cast<uint128_t>(MAGICB_64) << 64) | MAGICB_64_LOW;

// Constant that can be changed to distinguish different applications. It should
// not be zero.
constexpr uint32_t HASH_SALT = 1U;

// Variant 4 of Stafford's mixing algorithms. This is the same mixing algorithm
// used in splitmix64's 32-bit algorithm.
// URL: http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html
static inline uint32_t finalmix(uint64_t u) {
  u = (u ^ (u >> 33)) * 0x62a9d9ed799705f5;
  u = (u ^ (u >> 28)) * 0xcb24d0a5c88c35b3;
  return u >> 32;
}

// Ironseed Algorithm B
//
// Ironseed hashing is used to initialize the parameter space of an XCG
// generator. See https://github.com/reedacartwright/ironseed for more
// information about ironseed. By using both hashing and mixing, this method
// generates random looking results with excellent avalanche properties.
uint32_t ironseed_hash_once(uint64_t &m, const auto &values) {
  m += MAGICB_64;
  uint64_t entropy = m * 1;
  for (uint64_t u : values) {
    m += MAGICB_64;
    entropy += m * static_cast<uint32_t>(u);
    m += MAGICB_64;
    entropy += m * static_cast<uint32_t>(u >> 32);
  }
  m += MAGICB_64;
  entropy += m * HASH_SALT;
  return finalmix(entropy);
}

consteval bool is_pow2(std::size_t n) { return (n & (n - 1)) == 0; }

template <bool USE_LCG_> struct base_rng {
  uint128_t state = 1U;
};

template <> struct base_rng<true> {
  uint128_t state = 0U;
  uint128_t inc = 1U;
};

template <std::size_t SALT_N_> struct salt_array {
  std::array<uint64_t, SALT_N_> salts{};
};

template <> struct salt_array<0> {};

// Helpers for XCG concept
template <typename> struct is_xcg : std::false_type {};

template <bool B, std::size_t N> struct is_xcg<xcg<B, N>> : std::true_type {};

} // namespace detail

template <typename T>
concept XCG = detail::is_xcg<std::remove_cvref_t<T>>::value;

void seed(auto &x, const auto &values);


// XCG template
//
//  USE_LCG_ = use an LCG (true) or MCG (false)
//  SALT_N_ = size of the salt array
//
template <bool USE_LCG_, std::size_t SALT_N_>
struct xcg : detail::base_rng<USE_LCG_>, detail::salt_array<SALT_N_> {
  static_assert(detail::is_pow2(SALT_N_), "SALT_N_ must be a power of two");

  using xcg_type = xcg<USE_LCG_, SALT_N_>;

  static constexpr uint64_t MULT = XCG_MULT;
  static constexpr bool USE_LCG = USE_LCG_;
  static constexpr std::size_t SALT_N = SALT_N_;
  static constexpr std::size_t PARAM_SIZE =
      sizeof(uint64_t) * (2 + 2 * USE_LCG_ + SALT_N_);

  constexpr xcg() = default;

  // Generate a uniformly random uint64_t between [0, 2^64)
  uint64_t operator()() {
    uint128_t u = advance_();
    uint64_t high = static_cast<uint64_t>(u >> 64);
    uint64_t low = static_cast<uint64_t>(u);

    return high ^ random_salt_(low);
  }

  // Generate a uniformly random uint64_t between [0, range) using a 128-bit
  // short product. Example algorithms:
  // https://github.com/apple/swift/pull/39143 and
  // https://github.com/KWillets/range_generator/blob/master/include/range_generator.hpp
  // Has been called Cannon's method after https://github.com/stephentyrone
  // <stephentyrone@gmail.com>.
  //
  // Result is floor(range * u / 2^128) where u uniform in [0, 2^128)
  //
  // 64-bit range: sample >> 2^162 values to detect a bias
  // 32-bit range: sample >> 2^210
  //  K-bit range: sample >> 2^(258 - 3K/2)
  //
  uint64_t operator()(uint64_t range) {
    uint128_t u = advance_();
    uint64_t high = static_cast<uint64_t>(u >> 64);
    uint64_t low = static_cast<uint64_t>(u);

    high ^= random_salt_(low);
    low ^= random_salt_(high);

    //  Let u = (high * 2^64 + low). Then
    // 
    //  (range * u / 2^128)
    //      = (range * high * 2^64) / 2^128 + (range * low) / 2^128
    //      = (a * 2^64 + b) / 2^64 + (c * 2^64 + d) / 2^128
    //      = (a + b / 2^64 + c / 2^64 + d / 2^128)
    //      = (a + (b + c) / 2^64 + d / 2^128)
    //
    //  Now
    //
    //  floor(range * u / 2^128)
    //      = a + floor( ((b + c) * 2^64 + d) / 2^128 )
    //
    //  Since ((b + c) * 2^64 + d) < 2^129 the remainder is either 0 or 1. And
    //  it is 0 if (b + c) < 2^64 and 1 otherwise. d has no impact on the floor
    //  operation.

    uint128_t x = high;
    x = x * range;
    uint64_t a = static_cast<uint64_t>(x >> 64);
    uint64_t b = static_cast<uint64_t>(x);

    uint128_t y = low;
    y = y * range;
    uint64_t c = static_cast<uint64_t>(y >> 64);

    // (b + c < b) compiles to carry flag.
    return a + ((b + c < b) ? 1 : 0);
  }

  inline uint64_t random_salt_(uint64_t value) {
    (void)value; // To silence any warnings that `value` is not used.
    if constexpr (SALT_N_ == 0) {
      return 0;
    } else if constexpr (SALT_N_ == 1) {
      return this->salts[0];
    } else {
      constexpr unsigned int shift = 64 - std::countr_zero(SALT_N_);
      return this->salts[value >> shift];
    }
  }

  inline uint128_t advance_() {
    uint128_t u = this->state;
    this->state *= MULT;
    if constexpr (USE_LCG_) {
      this->state += this->inc;
    }
    return u;
  }

};

// Generate a uniformly random uint64_t between [0, 2^64)
template <XCG xcg_t> uint64_t random_u64(xcg_t &gen) { return gen(); }

// Generate a uniformly random uint64_t between [0, range)
template <XCG xcg_t> uint64_t random_u64(xcg_t &gen, uint64_t range) {
  return gen(range);
}

// Generate a uniformly random uint64_t between [umin, umax)
template <XCG xcg_t>
uint64_t random_u64(xcg_t &gen, uint64_t umin, uint64_t umax) {
  return umin + gen(umax - umin);
}

// Concept for a range that is convertible to a uint64_t
template <typename R>
concept ULongRange =
    std::ranges::range<R> &&
    std::convertible_to<std::ranges::range_value_t<R>, uint64_t>;

// Seed using a range of values
template <XCG xcg_t, ULongRange Range>
void seed(xcg_t &gen, const Range &values) {
  // Fill a buffer with 64-bit data created by hashing `values`.
  constexpr std::size_t N_ = xcg_t::PARAM_SIZE / sizeof(uint32_t);
  std::array<uint32_t, N_> buffer;
  uint64_t magic = 0;
  for (auto &&a : buffer) {
    a = detail::ironseed_hash_once(magic, values);
  }
  // Copy buffer data into parameter space.
  std::memcpy(&gen, buffer.data(), xcg_t::PARAM_SIZE);
  // Fixup state.
  if constexpr (xcg_t::USE_LCG) {
    gen.inc |= 1;
  } else {
    gen.state |= 1;
  }
}

template <XCG xcg_t, std::convertible_to<uint64_t> T>
void seed(xcg_t &gen, std::initializer_list<T> il) {
  // Use a span here to avoid an infinite recursion.
  seed(gen, std::span<const T>{il});
}

template <XCG xcg_t, std::convertible_to<uint64_t>... Args>
void seed(xcg_t &gen, Args &&...args) {
  seed(gen, std::initializer_list<uint64_t>{
                static_cast<uint64_t>(std::forward<Args>(args))...});
}

// Jump XCG state by 2^64 steps. The method comes from Brown (1994) Random
// number generation with arbitrary stride. Transactions of the American Nuclear
// Society. 71. Code adapted from PCG.
template <XCG xcg_t> auto permute_state(xcg_t gen) {
  uint128_t cur_mult = xcg_t::MULT;
  uint128_t cur_plus = 0;
  if constexpr (xcg_t::USE_LCG) {
    cur_plus = gen.inc;
  }
  for (int i = 0; i < 64; ++i) {
    cur_plus *= (cur_mult + 1);
    cur_mult *= cur_mult;
  }
  gen.state = cur_mult * gen.state + cur_plus;
  return gen;
}

// Jump XCG state by 2^96 steps.
template <XCG xcg_t> auto permute_state_huge(xcg_t gen) {
  uint128_t cur_mult = xcg_t::MULT;
  uint128_t cur_plus = 0;
  if constexpr (xcg_t::USE_LCG) {
    cur_plus = gen.inc;
  }
  for (int i = 0; i < 96; ++i) {
    cur_plus *= (cur_mult + 1);
    cur_mult *= cur_mult;
  }
  gen.state = cur_mult * gen.state + cur_plus;
  return gen;
}

// Permute XCG increment using a Weyl sequence while keeping it odd.
template <XCG xcg_t> auto permute_increment(xcg_t gen) {
  static_assert(xcg_t::USE_LCG,
                "Updating increment requires an XCG with an increment.");
  // Make sure that our magic constant is 2*(an odd number).
  constexpr uint128_t w = 2 * detail::MAGICB_128;
  gen.inc += w;
  return gen;
}

// Permute XCG salts by treating them as an N-bit number. The number of possible
// permutations will scale with the length of the salt array. By using add-with-
// carry, the permutation forms a Weyl sequence. Permuting the salts does not
// change the underlying generator, and this can be detected by XORing two
// related streams together.
template <XCG xcg_t> auto permute_salts(xcg_t gen) {
  static_assert(xcg_t::SALT_N > 0,
                "Updating salt requires an XCG that uses salts.");

  unsigned long long carry = 0;
  for (std::size_t i = 0; i < gen.salts.size(); ++i) {
    gen.salts[i] =
        __builtin_addcll(gen.salts[i], detail::MAGICB_64, carry, &carry);
  }
  return gen;
}

namespace utility {

inline int64_t i64(uint64_t u) { return u; }
inline int64_t i63(uint64_t u) { return u >> 1; }
inline int64_t i54s(uint64_t u) { return ((int64_t)u) >> 10; }
inline int64_t i53(uint64_t u) { return u >> 11; }

// uniform in [0,1] with variable steps
inline double f64(uint64_t u) { return u * (1.0 / 18446744073709551616.0); }

// uniform in [-1,1] with variable steps
inline double f64s(uint64_t u) {
  return i64(u) * (1.0 / 9223372036854775808.0);
}

// uniform in [0,1] with variable steps
inline double f63(uint64_t u) { return i63(u) * (1.0 / 9223372036854775808.0); }

// uniform in [-1,1) with equal steps
inline double f54(uint64_t u) { return i54s(u) * (1.0 / 9007199254740992.0); }

// uniform in (-1,1] with equal steps
inline double f54a(uint64_t u) {
  return (i54s(u) + 1) * (1.0 / 9007199254740992.0);
}

// uniform in (-1,1) with equal steps
inline double f54b(uint64_t u) {
  return (i54s(u) | 1) * (1.0 / 9007199254740992.0);
}

// uniform in [0,1) with equal steps
inline double f53(uint64_t u) { return i53(u) * (1.0 / 9007199254740992.0); }

// uniform in (0,1] with equal steps
inline double f53a(uint64_t u) {
  return (i53(u) + 1) * (1.0 / 9007199254740992.0);
}

// uniform in (0,1) with equal steps
inline double f53b(uint64_t u) {
  return (i53(u) | 1) * (1.0 / 9007199254740992.0);
}

} // namespace utility

// Generate a value between [0, 1.0)
template <XCG xcg_t> double random_f53(xcg_t &gen) {
  return utility::f53(gen());
}

// Generate a value between (0, 1.0]
template <XCG xcg_t> double random_f53a(xcg_t &gen) {
  return utility::f53a(gen());
}

// Generate a value between (0, 1.0)
template <XCG xcg_t> double random_f53b(xcg_t &gen) {
  return utility::f53b(gen());
}

// Generate a value between [-1, 1.0)
template <XCG xcg_t> double random_f54(xcg_t &gen) {
  return utility::f54(gen());
}

// Generate a value between (-1, 1.0]
template <XCG xcg_t> double random_f54a(xcg_t &gen) {
  return utility::f54a(gen());
}

// Generate a value between (-1, 1.0)
template <XCG xcg_t> double random_f54b(xcg_t &gen) {
  return utility::f54b(gen());
}

// Generate a value between [0, range) with absolutely no bias.
template <XCG xcg_t>
uint64_t random_u64_bounded_exact(xcg_t &gen, uint64_t range) {
  uint128_t xx = gen();
  xx = xx * range;
  uint64_t y = static_cast<uint64_t>(xx >> 64);
  uint64_t f = static_cast<uint64_t>(xx);
  // Optimize for small ranges
  if (range + f >= f) [[likely]] {
    return y;
  }
  do {
    xx = gen();
    xx = xx * range;
    uint64_t z = static_cast<uint64_t>(xx >> 64);
    z = z + f;
    if (z < f) {
      // we have carried
      return y + 1;
    } else if (z != -1) {
      // we will never carry
      break;
    }
    f = static_cast<uint64_t>(xx);
  } while (range + f < f);
  return y;
}

template <XCG xcg_t>
uint64_t random_u64_bounded_exact(xcg_t &gen, uint64_t umin, uint64_t umax) {
  return umin + random_u64_bounded_exact(gen, umax - umin);
}

using xcg128_t = xcg<false, 0>;
using xcg256_t = xcg<true, 0>;
using xcg320_t = xcg<true, 1>;
using xcg384_t = xcg<true, 2>;
using xcg512_t = xcg<true, 4>;
using xcg768_t = xcg<true, 8>;
using xcg1280_t = xcg<true, 16>;

}; // namespace xcg

#endif
