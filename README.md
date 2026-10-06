# XCG: Expanded Congruential Generator

XCG is a family of random number generators inspired by and derived from PCG's
extended generation scheme <https://www.pcg-random.org/>. At its core, XCG is a
128-bit linear congruential generator that returns its high 64 bits. This
simple generator is enough to pass PractRand 0.96 tests at 32 TB except for
TMFn, which is designed to detect LCGs. For flexibility, the family provides
generators that vary in the sizes of their parameter spaces. For example,
XCG-1280 has a 1280-bit parameter space: 128 bits for the state of the LCG,
128 bits for the increment of the LCG, and 1024 bits for a set of 64-bit salts
that are mixed into the output.

XCG is likely not cryptographically secure and should not be used where output
must be unpredictable to an attacker.

## Period

The period of XCG is `2^126` if using the MCG core generator or `2^128` if using
the LCG core generator. While the increment and salts expand the size of XCG's
parameter space, they do not extend XCG's period. The periods of the core
generators are sufficiently large that exhausting a stream is impractical for
ordinary applications. This simplifies XCG's algorithm without sacrificing
utility or parameter-space flexibility.

## Rotation

The standard XCG algorithms fail the TMFn tests found in PractRand 0.96. These
tests are designed to detect LCGs. The rotating XCG algorithms add a random
permutation step (bit rotation) to the standard algorithms. This permutation
step is enough to pass the TMFn tests at 32 TB.

## Multiplier

XCG uses the 64-bit multiplier `0xf68a43306d8e0225`, which was identified
through searching for the best 64-bit multiplier that works for 128-bit LCGs,
128-bit MCGs, 64-bit LCGs, and 64-bit MCGs.

## Example Implementations

```cpp
uint64_t XCG_MULT = 0xf68a43306d8e0225U;

struct xcg128_t {
  uint128_t state = 1U;
};

inline uint64_t get_u64(xcg128_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  return (u >> 64);
}

inline uint64_t get_u64r(xcg128_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  uint64_t high = (uint64_t)(u >> 64);
  uint64_t low = (uint64_t)u;
  return std::rotr(high, low >> 58);
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

inline uint64_t get_u64r(xcg256_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  gen->state += gen->inc;
  uint64_t high = (uint64_t)(u >> 64);
  uint64_t low = (uint64_t)u;
  return std::rotr(high, low >> 58);
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

inline uint64_t get_u64r(xcg320_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  gen->state += gen->inc;
  uint64_t high = (uint64_t)(u >> 64);
  uint64_t low = (uint64_t)u;
  return std::rotr(high, low >> 58) ^ gen->salt;
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

inline uint64_t get_u64r(xcg512_t *gen) {
  uint128_t u = gen->state;
  gen->state *= XCG_MULT;
  gen->state += gen->inc;
  uint64_t high = (uint64_t)(u >> 64);
  uint64_t low = (uint64_t)u;
  return std::rotr(high, low >> 58) ^ gen->salts[low >> 62];
}
```
