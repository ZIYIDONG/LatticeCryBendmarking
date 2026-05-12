#include "../include/frd.h"
#include <iostream>
#include "unified_params.h"

using namespace cryptolib;

int main() {
    // F_7 上 x^2 + 1 是不是不可约? 7 mod 4 == 3, 所以 -1 不是 QR, 是的
    auto __u_p = unified::default_mp12_params(); long q = __u_p.q;
    Vec f = {1, 0, 1};   // 1 + 0·x + x^2 = x^2 + 1
    std::cout << "f = x^2 + 1, deg = " << poly_deg(f) << "\n";
    std::cout << "is_irreducible(x^2+1, F_7) = " << is_irreducible(f, q) << "\n";

    // F_7 上 x^2 - 2 = x^2 + 5 (因为 -2 mod 7 = 5)
    // 7 上 2 是 QR? 2,4,1,2,4,1 -> 2 is QR (3^2=9=2). 所以 x^2-2 可约
    Vec f2 = {5, 0, 1};
    std::cout << "f2 = x^2 + 5, is_irreducible? " << is_irreducible(f2, q) << "\n";

    // x^2 + x + 3 ?
    Vec f3 = {3, 1, 1};
    std::cout << "f3 = x^2 + x + 3, is_irreducible? " << is_irreducible(f3, q) << "\n";

    // 试着搜索
    for (long c = 0; c < q; c++) {
        for (long b = 0; b < q; b++) {
            Vec f = {c, b, 1};
            if (is_irreducible(f, q)) {
                std::cout << "found: x^2 + " << b << "x + " << c << "\n";
            }
        }
    }
    return 0;
}
