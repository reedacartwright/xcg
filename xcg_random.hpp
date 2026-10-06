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
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>

static_assert(__cplusplus >= 202002L, "XCG requires C++20");

#if !defined(__SIZEOF_INT128__)
static_assert(false, "XCG requires __uint128_t support");
#endif

namespace xcg {

using uint128_t = __uint128_t;

/*
This coefficient came from searching parameter space for the best 64-bit
multiplier that works for 128-bit LCG, 128-bit MCG, 64-bit LCG, and 64-bit MCG.
*/
constexpr uint64_t XCG_MULT = 0xf68a43306d8e0225U;
// 0b1111011010001010010000110011000001101101100011100000001000100101

template <bool USE_LCG_ = false, std::size_t SALT_N_ = 0,
          bool USE_ROTR_ = false>
struct xcg;

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
    (static_cast<uint128_t>(MAGICB_64) << 64U) | MAGICB_64_LOW;

// Constant that can be changed to distinguish different applications. It should
// not be zero.
constexpr uint32_t HASH_SALT = 1U;

// Variant 4 of Stafford's mixing algorithms. This is the same mixing algorithm
// used in splitmix64's 32-bit algorithm.
// URL: http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html
constexpr uint32_t finalmix(uint64_t u) noexcept {
  u = (u ^ (u >> 33U)) * 0x62a9d9ed799705f5;
  u = (u ^ (u >> 28U)) * 0xcb24d0a5c88c35b3;
  return static_cast<uint32_t>(u >> 32U);
}

// Ironseed Algorithm B
//
// Ironseed hashing is used to initialize the parameter space of an XCG
// generator. See https://github.com/reedacartwright/ironseed for more
// information about ironseed. By using both hashing and mixing, this method
// generates random looking results with excellent avalanche properties.
constexpr uint32_t ironseed_hash_once(uint64_t &m, const auto &values) {
  m += MAGICB_64;
  uint64_t entropy = m * 1;
  for (uint64_t u : values) {
    m += MAGICB_64;
    entropy += m * static_cast<uint32_t>(u);
    m += MAGICB_64;
    entropy += m * static_cast<uint32_t>(u >> 32U);
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

template <bool B, std::size_t N, bool R>
struct is_xcg<xcg<B, N, R>> : std::true_type {};

} // namespace detail

template <typename T>
concept XCG = detail::is_xcg<std::remove_cvref_t<T>>::value;

// XCG template
//
//  USE_LCG_ = use an LCG (true) or MCG (false)
//  SALT_N_ = size of the salt array
//  USE_ROT_ = if true, permute high bits using a random rotation based on low
//
template <bool USE_LCG_, std::size_t SALT_N_, bool USE_ROT_>
struct xcg : detail::base_rng<USE_LCG_>, detail::salt_array<SALT_N_> {
  static_assert(detail::is_pow2(SALT_N_), "SALT_N_ must be a power of two");

  using xcg_type = xcg;

  static constexpr uint64_t MULT = XCG_MULT;
  static constexpr bool USE_LCG = USE_LCG_;
  static constexpr bool USE_ROT = USE_ROT_;
  static constexpr std::size_t SALT_N = SALT_N_;
  static constexpr std::size_t PARAM_SIZE =
      sizeof(uint64_t) * (2 + 2 * USE_LCG_ + SALT_N_);

  constexpr xcg() = default;

  // Generate a uniformly random uint64_t between [0, 2^64)
  constexpr uint64_t operator()() noexcept {
    uint128_t u = advance_and_salt_();
    return static_cast<uint64_t>(u >> 64U);
  }

  // compatibility with <random>
  using result_type = uint64_t;

  static constexpr result_type min() noexcept {
    return std::numeric_limits<result_type>::min();
  }

  static constexpr result_type max() noexcept {
    return std::numeric_limits<result_type>::max();
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
  constexpr uint64_t operator()(uint64_t range) noexcept {
    uint128_t u = advance_and_salt_();
    auto high = static_cast<uint64_t>(u >> 64U);
    auto low = static_cast<uint64_t>(u);

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
    auto a = static_cast<uint64_t>(x >> 64U);
    auto b = static_cast<uint64_t>(x);

    // Uncomment to optimize for "small" ranges
    // Note this can have worse performance if range > 2^58
    // if (range + b >= range) [[likely]] {
    //   return a;
    // }

    uint128_t y = low;
    y = y * range;
    auto c = static_cast<uint64_t>(y >> 64U);

    // (b + c < c) compiles to carry flag.
    return a + ((b + c < c) ? 1 : 0);
  }

  constexpr uint64_t random_salt_(uint64_t value) const noexcept {
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

  constexpr uint128_t advance_() noexcept {
    uint128_t u = this->state;
    this->state *= MULT;
    if constexpr (USE_LCG_) {
      this->state += this->inc;
    }
    return u;
  }

  constexpr uint128_t advance_and_salt_() noexcept {
    uint128_t u = advance_();
    auto high = static_cast<uint64_t>(u >> 64U);
    auto low = static_cast<uint64_t>(u);

    // Permute high bits. Don't permute the low bits because for MCGs it will
    // spread the bias of bits 0 and 1 to the rest. Intentionally use the
    // highest bits for both rotating and salting to avoid extra instructions.
    if constexpr (USE_ROT_) {
      high = std::rotr(high, static_cast<int>(low >> 58U));
    }
    // Salt `high` using `low`. Salt `low` using `0`. This allows for ILP. Avoid
    // cross-salting because that may not be bijective.
    high ^= random_salt_(low);
    low ^= random_salt_(0);

    return (static_cast<uint128_t>(high) << 64U) | low;
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
    std::ranges::forward_range<const R> &&
    std::convertible_to<std::ranges::range_value_t<R>, uint64_t>;

// Seed using a range of values
template <XCG xcg_t, ULongRange Range>
constexpr void seed(xcg_t &gen, const Range &values) {
  uint64_t magic = 0;

  // Helper functions
  auto next_u32 = [&]() { return detail::ironseed_hash_once(magic, values); };
  auto next_u64 = [&]() {
    uint64_t low = next_u32();
    uint64_t high = next_u32();

    return low | (high << 32U);
  };

  auto next_u128 = [&]() {
    uint128_t low = next_u64();
    uint128_t high = next_u64();

    return low | (high << 64U);
  };

  // Initialize state space.
  gen.state = next_u128();

  if constexpr (xcg_t::USE_LCG) {
    gen.inc = next_u128();
    gen.inc |= 1;
  } else {
    gen.state |= 1;
  }

  if constexpr (xcg_t::SALT_N > 0) {
    for (uint64_t &salt : gen.salts) {
      salt = next_u64();
    }
  }
}

template <XCG xcg_t, std::convertible_to<uint64_t> T>
constexpr void seed(xcg_t &gen, std::initializer_list<T> il) {
  // Use a span here to avoid an infinite recursion.
  seed(gen, std::span<const T>{il});
}

template <XCG xcg_t, std::convertible_to<uint64_t>... Args>
constexpr void seed(xcg_t &gen, Args &&...args) {
  seed(gen, std::initializer_list<uint64_t>{
                static_cast<uint64_t>(std::forward<Args>(args))...});
}

// Jump XCG state by 2^64 steps. The method comes from Brown (1994) Random
// number generation with arbitrary stride. Transactions of the American Nuclear
// Society. 71. Code adapted from PCG.
template <XCG xcg_t>
[[nodiscard]]
auto permute_state(xcg_t gen) {
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
template <XCG xcg_t>
[[nodiscard]]
auto permute_state_huge(xcg_t gen) {
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
template <XCG xcg_t>
[[nodiscard]]
auto permute_increment(xcg_t gen) {
  static_assert(xcg_t::USE_LCG,
                "Updating increment requires an XCG with an increment.");
  // Make sure that our magic constant is 2*(an odd number).
  constexpr uint128_t w = 2 * detail::MAGICB_128;
  gen.inc += w;
  return gen;
}

// Permute XCG salts by treating the array as a large number. The number of
// possible permutations will scale with the length of the salt array. By using
// add-with-carry, the permutation forms a Weyl sequence. Permuting the salts
// does not change the underlying generator, and this can be detected by XORing
// two related streams together.
template <XCG xcg_t>
[[nodiscard]]
auto permute_salts(xcg_t gen) {
  static_assert(xcg_t::SALT_N > 0,
                "Updating salt requires an XCG that uses salts.");
  uint64_t carry = 0;
  for (uint64_t &salt : gen.salts) {
    uint128_t sum = static_cast<uint128_t>(salt) + detail::MAGICB_64 + carry;
    salt = static_cast<uint64_t>(sum);
    carry = static_cast<uint64_t>(sum >> 64U);
  }

  return gen;
}

namespace utility {

constexpr int64_t i64(uint64_t u) noexcept { return static_cast<int64_t>(u); }
constexpr int64_t i63(uint64_t u) noexcept {
  return static_cast<int64_t>(u >> 1U);
}
constexpr int64_t i54s(uint64_t u) noexcept {
  return static_cast<int64_t>(u) >> 10U;
}
constexpr int64_t i53(uint64_t u) noexcept {
  return static_cast<int64_t>(u >> 11U);
}

// uniform in [0,1] with variable steps
constexpr double f64(uint64_t u) noexcept {
  return static_cast<double>(u) * (1.0 / 18446744073709551616.0);
}

// uniform in [-1,1] with variable steps
constexpr double f64s(uint64_t u) noexcept {
  return static_cast<double>(i64(u)) * (1.0 / 9223372036854775808.0);
}

// uniform in [0,1] with variable steps
constexpr double f63(uint64_t u) noexcept {
  return static_cast<double>(i63(u)) * (1.0 / 9223372036854775808.0);
}

// uniform in [-1,1) with equal steps
constexpr double f54(uint64_t u) noexcept {
  return static_cast<double>(i54s(u)) * (1.0 / 9007199254740992.0);
}

// uniform in (-1,1] with equal steps
constexpr double f54a(uint64_t u) noexcept {
  return static_cast<double>(i54s(u) + 1) * (1.0 / 9007199254740992.0);
}

// uniform in (-1,1) with equal steps
constexpr double f54b(uint64_t u) noexcept {
  return static_cast<double>(i54s(u) | 1) * (1.0 / 9007199254740992.0);
}

// uniform in [0,1) with equal steps
constexpr double f53(uint64_t u) noexcept {
  return static_cast<double>(i53(u)) * (1.0 / 9007199254740992.0);
}

// uniform in (0,1] with equal steps
constexpr double f53a(uint64_t u) noexcept {
  return static_cast<double>(i53(u) + 1) * (1.0 / 9007199254740992.0);
}

// uniform in (0,1) with equal steps
constexpr double f53b(uint64_t u) noexcept {
  return static_cast<double>(i53(u) | 1) * (1.0 / 9007199254740992.0);
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

namespace detail {
template <XCG xcg_t>
uint64_t random_u64_bounded_exact_tail(xcg_t &gen, uint64_t range, uint64_t h0,
                                       uint64_t f0);
}

/*
## Generate a value between [0, range) with absolutely no bias.

Let `U0` be uniformly distributed in [0, 1). Then `B = floor(R*U0)` is a
uniformly distributed integer bounded in [0, R).

Consider, `Xi` be a 64-bit integer uniformly distributed in [0, 2^64). Then

    U0 = X0 * 2^-64 + X1 * 2^-128 + X2 * 2^-192 + ...
       = 2^-64 * (sum Xi * 2^(-64*i)) where i is [0, infinity)

More generally, let `Uj = sum X{i+j} * 2^(-64*(i+1))` where i is [0, infinity).

And

    R*U0 = R * X0 * 2^-64 + R * X1 * 2^-128 + R * X2 * 2^-192 + ...

Next let `R * Xi = Hi * 2^64 + Fi` Then

    R*U0 = 2^-64 * (H0 * 2^64 + F0) + 2^-64 * R  * U1
         = H0 + F0 / 2^64 + R * U1 / 2^64

And

    floor(R*U0) = H0 + floor( (F0 + R * U1) / 2^64 ) = H0 + K

Since `F0 < 2^64` and `R * U1 < 2^64`, K is 0 or 1. Note that if
`F0 + R <= 2^64` then `F0 + R * U1 < 2^64` and `K = 0`.

Let

    F0 + R * U1 = F0 + 2^-64 * (H1 * 2^64 + F1) + 2^-64 * R  * U2
                = F0 + H1 + F1 / 2^64 + R * U2 / 2^64

Following similar logic to above `(F1 + R * U2) / 2^64 < 2`. Therefore,

    - If `F0 + H1 <= 2^64 - 2`, then `K = 0`
    - If `F0 + H1 >= 2^64` then `K = 1`.
    - If `F0 + H1 == 2^64 - 1` the `K = 0 or 1` and more data is needed.
*/
template <XCG xcg_t>
inline uint64_t random_u64_bounded_exact(xcg_t &gen, uint64_t range) {
  uint128_t x = gen();
  x *= range;
  auto h0 = static_cast<uint64_t>(x >> 64U);
  auto f0 = static_cast<uint64_t>(x);
  // Optimize for small ranges
  if (range + f0 >= range) [[likely]] {
    return h0;
  }
  // Use a tail call to allow loops to inline better.
  return detail::random_u64_bounded_exact_tail(gen, range, h0, f0);
}

namespace detail {

template <XCG xcg_t>
uint64_t random_u64_bounded_exact_tail(xcg_t &gen, uint64_t range, uint64_t h0,
                                       uint64_t f0) {
  do {
    uint128_t x = gen();
    x *= range;
    auto h1 = static_cast<uint64_t>(x >> 64U);
    f0 += h1;
    if (f0 < h1) { // if F0 + H1 >= 2^64
      // we have carried
      return h0 + 1;
    } else if (f0 < std::numeric_limits<uint64_t>::max()) {
      // we will never carry
      break;
    }
    f0 = static_cast<uint64_t>(x);
  } while (range + f0 < range);
  return h0;
}

} // namespace detail

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

} // namespace xcg

#endif
