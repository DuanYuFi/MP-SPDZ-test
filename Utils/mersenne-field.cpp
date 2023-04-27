#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "Tools/time-func.h"
#include "Tools/intrinsics.h"
#include "Math/gf2nlong.h"

using namespace std;

class Mersenne {
    public:
        typedef unsigned __int128 uint128_t;

        static const uint32_t PRIME_EXP = 61;
        static const uint64_t PR = 2305843009213693951;

        static uint64_t modp(uint64_t);
        static uint64_t modp_128(uint128_t);
        static uint64_t neg(uint64_t);
        static uint64_t add(uint64_t, uint64_t);
        static uint64_t sub(uint64_t, uint64_t);
        static uint64_t mul(uint64_t, uint64_t);
        static uint64_t inverse(uint64_t a);
        static uint64_t inner_product(uint64_t*, uint64_t*, uint64_t);
        static uint64_t batch_sum(uint64_t*, uint64_t);
        static uint64_t batch_add(uint64_t*, uint64_t*, uint64_t);
        static uint64_t get_rand();
};

inline uint64_t Mersenne::modp(uint64_t a) {
    uint64_t res = (a>>PRIME_EXP) + (a & PR);
    if (res >= PR) {
        res -= PR;
    }
    return res;
}

inline uint64_t Mersenne::modp_128(uint128_t a){
    uint64_t higher, middle, lower;
    higher = (a >> (2 * PRIME_EXP));
    middle = (a >> PRIME_EXP) & PR;
    lower = a & PR;
    return modp(higher + middle + lower);
}

inline uint64_t Mersenne::neg(uint64_t a) {
    if (a > 0) {
        return PR - a;
    } else {
        return 0;
    }
}

inline uint64_t Mersenne::add(uint64_t a, uint64_t b) {
    uint64_t res = a + b;
    if (res >= PR) {
        res -= PR;
    }
    return res;
}

inline uint64_t Mersenne::sub(uint64_t a, uint64_t b) {
    if (a >= b) {
        return a - b;
    } else {
        return PR - b + a;
    }
}

inline uint64_t Mersenne::mul(uint64_t a, uint64_t b) {
    int128 res = _mm_clmulepi64_si128((uint128_t) a, (uint128_t) b, 0);
    uint64_t higher = res.get_lower();
    uint64_t lower = res.get_upper();
    return add(higher, lower);
}

inline uint64_t Mersenne::inverse(uint64_t a) {
    uint64_t left = a;
    uint64_t right = PR;
    uint64_t x = 1, y = 0, u = 0, v = 1;
    // uint64_t gcd = a;
    uint64_t w, z;
    while(left != 0) {
        w = right / left;
        z = right % left;
        right = left;
        left = z;

        z = u - w * x;
        u = x;
        x = z;

        z = v - w * y;
        v = y;
        y = z;
    }
    if (u >= PR) {
        u += PR;
    }
    return u;
}

inline uint64_t Mersenne::inner_product(uint64_t* a, uint64_t* b, uint64_t size) {
    uint128_t result = 0;
    uint64_t bound = 63;
    uint64_t start, end;
    start = 0;
    while(true) {
        if (start + bound < size) {
            end = start + bound;
        }
        else {
            end = size;
        }
        for(uint64_t i = start; i < end; i++) {
            result += ((uint128_t)a[i]) * ((uint128_t)b[i]);
        }
        result = modp_128(result);
        start = end;
        if (start == size) break;
    }
    return result;
}

inline uint64_t Mersenne::batch_add(uint64_t* a, uint64_t* b, uint64_t size) {
    uint128_t result = 0;
    uint64_t bound = 63;
    uint64_t start, end;
    start = 0;
    while(true) {
        if (start + bound < size) {
            end = start + bound;
        }
        else {
            end = size;
        }
        for(uint64_t i = start; i < end; i++) {
            result += ((uint128_t)a[i]) + ((uint128_t)b[i]);
        }
        result = modp_128(result);
        start = end;
        if (start == size) break;
    }
    return result;
}

inline uint64_t Mersenne::batch_sum(uint64_t* a, uint64_t size) {
    uint128_t result = 0;
    uint64_t bound = 63;
    uint64_t start, end;
    start = 0;
    while(true) {
        if (start + bound < size) {
            end = start + bound;
        }
        else {
            end = size;
        }
        for(uint64_t i = start; i < end; i++) {
            result += (uint128_t)a[i];
        }
        result = modp_128(result);
        start = end;
        if (start == size) break;
    }
    return result;
}

inline static uint64_t get_rand() {
    uint64_t res = 0;
    res |= rand();
    res = res << 16;
    res |= rand();
    res = res << 16;
    res |= rand();
    res = res << 16;
    res |= rand();

    return res % Mersenne::PR;
}

inline uint64_t normal_inner_product(uint64_t* a, uint64_t* b, uint64_t size) {
    uint64_t result = 0;
    for(uint64_t i = 0; i < size; i++) {
        result += a[i] * b[i];
    }
    return result;
}

void inner_product_test(const int TEST_SIZE = 1e8) {
    uint64_t *a = new uint64_t[TEST_SIZE];
    uint64_t *b = new uint64_t[TEST_SIZE];

    for(int i = 0; i < TEST_SIZE; i++) {
        a[i] = get_rand();
        b[i] = get_rand();
    }

    Timer timer;
    timer.start();

    uint64_t res = Mersenne::inner_product(a, b, TEST_SIZE);

    timer.stop();

    cout << "Result: " << res << endl;
    cout << "Time: " << timer.elapsed() << endl;

    timer.reset();
    timer.start();

    uint64_t res2 = normal_inner_product(a, b, TEST_SIZE);

    timer.stop();

    cout << "Result: " << res2 << endl;
    cout << "Time: " << timer.elapsed() << endl;
}

int main() {

    uint64_t a = 0x1234567890abcdef;
    uint64_t b = 0x1234567890abcdef;

    uint64_t c = Mersenne::mul(a, b);

    cout << c << endl;

    return 0;


}