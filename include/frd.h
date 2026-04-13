#pragma once
/**
 * FRD: Full-Rank Difference Encoding (Agrawal-Boneh-Boyen, EUROCRYPT 2010)
 * ────────────────────────────────────────────────────────────────────────
 *
 * 数学定义:
 *   FRD : Z_q^n  →  Z_q^{n×n}
 *
 *   性质:
 *     ① 高效可计算
 *     ② 对任意 id1 ≠ id2 ∈ Z_q^n,
 *        H = FRD(id1) - FRD(id2) ∈ Z_q^{n×n} 满秩(可逆)
 *
 * 标准构造(ABB10):
 *   选取 F_q 上 n 次不可约多项式 f(x),从而 F_q[x]/(f(x)) ≅ F_{q^n}
 *
 *   把 id = (a_0, a_1, ..., a_{n-1}) ∈ Z_q^n 看作多项式
 *     a(x) = a_0 + a_1·x + ... + a_{n-1}·x^{n-1}  ∈ F_{q^n}
 *
 *   FRD(id) = M_a := [a(x), x·a(x), x²·a(x), ..., x^{n-1}·a(x)]
 *   即 M_a 的第 j 列是 x^j · a(x) mod f(x) 的系数向量
 *
 *   M_a 是 F_{q^n} 上"乘以 a(x)" 这个线性映射在基 {1,x,...,x^{n-1}} 下的矩阵
 *
 *   全秩证明:
 *     M_{a-b} = M_a - M_b  (线性性)
 *     若 a ≠ b,则 a-b ≠ 0 ∈ F_{q^n}(域!)
 *     乘以非零元在域中是双射 → M_{a-b} 满秩 ✓
 *
 * 要求:
 *   q 必须是素数(或素数幂),通常我们用 q 为素数
 *   存在 F_q 上 n 次不可约多项式(对几乎所有 n 都存在)
 *
 * 用途:
 *   - ABB10 IBE 中把身份映射成 tag matrix
 *   - DelTrapGen 中为不同身份生成不同 tag
 *   - 任何需要"无碰撞 + 差分可逆"的格密码方案
 */

#include <vector>
#include <random>
#include <stdexcept>
#include <cstdint>
#include <numeric>

namespace cryptolib {

using Vec = std::vector<long>;
using Mat = std::vector<Vec>;

inline long mod_pos(long x, long q) {
    return ((x % q) + q) % q;
}

inline Mat make_zero_mat(int r, int c) {
    return Mat(r, Vec(c, 0));
}

/* ══════════════════════════════════════════════════
   §1  F_q 上的多项式运算
   ══════════════════════════════════════════════════ */
/**
 * 多项式表示: 用 vector<long>,索引 i 对应 x^i 的系数
 * 长度可变,但通常截断到次数 deg+1
 */

// 多项式次数(去掉前导零)
inline int poly_deg(const Vec& p) {
    for (int i = (int)p.size() - 1; i >= 0; i--)
        if (p[i] != 0) return i;
    return -1;  // 零多项式
}

// 把多项式截断到精确长度 n
inline Vec poly_trim(const Vec& p, int n) {
    Vec r(n, 0);
    int copy_len = std::min((int)p.size(), n);
    for (int i = 0; i < copy_len; i++) r[i] = p[i];
    return r;
}

// 多项式加法 mod q
inline Vec poly_add(const Vec& a, const Vec& b, long q) {
    int n = (int)std::max(a.size(), b.size());
    Vec r(n, 0);
    for (int i = 0; i < n; i++) {
        long va = (i < (int)a.size()) ? a[i] : 0;
        long vb = (i < (int)b.size()) ? b[i] : 0;
        r[i] = mod_pos(va + vb, q);
    }
    return r;
}

// 多项式减法 mod q
inline Vec poly_sub(const Vec& a, const Vec& b, long q) {
    int n = (int)std::max(a.size(), b.size());
    Vec r(n, 0);
    for (int i = 0; i < n; i++) {
        long va = (i < (int)a.size()) ? a[i] : 0;
        long vb = (i < (int)b.size()) ? b[i] : 0;
        r[i] = mod_pos(va - vb, q);
    }
    return r;
}

// 多项式乘法 mod q (普通卷积,不约化)
inline Vec poly_mul(const Vec& a, const Vec& b, long q) {
    int da = poly_deg(a), db = poly_deg(b);
    if (da < 0 || db < 0) return Vec{0};
    Vec r(da + db + 1, 0);
    for (int i = 0; i <= da; i++)
        if (a[i] != 0)
            for (int j = 0; j <= db; j++)
                r[i + j] = mod_pos(r[i + j] + a[i] * b[j], q);
    return r;
}

// 模逆元(扩展欧几里得)
inline long mod_inv_q(long a, long q) {
    a = mod_pos(a, q);
    long t = 0, newt = 1, r = q, newr = a;
    while (newr != 0) {
        long quot = r / newr;
        long tmp;
        tmp = newt; newt = t - quot * newt; t = tmp;
        tmp = newr; newr = r - quot * newr; r = tmp;
    }
    if (r > 1) throw std::runtime_error("not invertible");
    return mod_pos(t, q);
}

// 多项式取模: a(x) mod f(x), 都在 F_q[x]
// 经典的"长除法"
inline Vec poly_mod(Vec a, const Vec& f, long q) {
    int df = poly_deg(f);
    if (df < 0) throw std::runtime_error("mod by zero polynomial");
    long inv_lead = mod_inv_q(f[df], q);  // f 的首项系数的逆

    int da = poly_deg(a);
    while (da >= df) {
        long coef = mod_pos(a[da] * inv_lead, q);
        // 从 a 中减去 coef · x^{da-df} · f
        for (int i = 0; i <= df; i++) {
            int idx = da - df + i;
            a[idx] = mod_pos(a[idx] - coef * f[i], q);
        }
        da = poly_deg(a);
    }
    return a;
}

/* ══════════════════════════════════════════════════
   §2  寻找 F_q 上 n 次不可约多项式
   ══════════════════════════════════════════════════ */
/**
 * Rabin 不可约性测试:
 *   F_q 上次数为 n 的多项式 f(x) 不可约 当且仅当
 *     ① x^{q^n} ≡ x  (mod f(x))
 *     ② 对 n 的每个素因子 p_i,gcd(x^{q^{n/p_i}} - x, f(x)) = 1
 *
 * 这里我们用一个简化但可靠的策略:
 *   随机生成 monic 首一多项式,用 Rabin 测试,直到找到不可约的为止
 *   平均尝试次数 ≈ n(对于"几乎所有"q,n)
 */

// 计算 x^k mod f(x) in F_q[x] (快速幂)
inline Vec poly_xpow_mod(long k, const Vec& f, long q) {
    int df = poly_deg(f);
    Vec result(1, 1);          // 1
    Vec base(2, 0); base[1] = 1; // x
    while (k > 0) {
        if (k & 1) result = poly_mod(poly_mul(result, base, q), f, q);
        base = poly_mod(poly_mul(base, base, q), f, q);
        k >>= 1;
    }
    return poly_trim(result, df);
}

// 用大整数指数版本(支持 q^n,可能很大,这里用 long long 演示)
// 真实 q^n 可能溢出,需要按 mod 分解;此处用于 q,n 较小的演示
inline Vec poly_xpow_mod_bigexp(long long exp, const Vec& f, long q) {
    int df = poly_deg(f);
    Vec result(1, 1);
    Vec x_poly(2, 0); x_poly[1] = 1;
    Vec base = x_poly;
    while (exp > 0) {
        if (exp & 1) result = poly_mod(poly_mul(result, base, q), f, q);
        base = poly_mod(poly_mul(base, base, q), f, q);
        exp >>= 1;
    }
    return poly_trim(result, df);
}

// 两个多项式的 gcd(用于不可约性测试)
inline Vec poly_gcd(Vec a, Vec b, long q) {
    while (poly_deg(b) >= 0) {
        Vec r = poly_mod(a, b, q);
        a = b;
        b = r;
    }
    // 把 a 归一化(首一)
    int da = poly_deg(a);
    if (da >= 0 && a[da] != 1) {
        long inv = mod_inv_q(a[da], q);
        for (auto& c : a) c = mod_pos(c * inv, q);
    }
    return a;
}

// 整数素因子分解(小整数,试除法)
inline std::vector<long> small_prime_factors(long n) {
    std::vector<long> factors;
    for (long p = 2; p * p <= n; p++) {
        if (n % p == 0) {
            factors.push_back(p);
            while (n % p == 0) n /= p;
        }
    }
    if (n > 1) factors.push_back(n);
    return factors;
}

// 比较两个多项式是否相等(忽略尾部零)
inline bool poly_equal(const Vec& a, const Vec& b) {
    return poly_deg(a) == poly_deg(b) &&
           [&] {
               int d = poly_deg(a);
               for (int i = 0; i <= d; i++)
                   if (a[i] != b[i]) return false;
               return true;
           }();
}

// Rabin 不可约性测试: f ∈ F_q[x] 是否是 n 次不可约多项式
inline bool is_irreducible(const Vec& f, long q) {
    int n = poly_deg(f);
    if (n <= 0) return false;
    if (n == 1) return true;

    // 计算 base^q mod f (用快速幂)
    auto qpow = [&](const Vec& base) -> Vec {
        Vec result(1, 1);
        Vec b = base;
        long e = q;
        while (e > 0) {
            if (e & 1) result = poly_mod(poly_mul(result, b, q), f, q);
            b = poly_mod(poly_mul(b, b, q), f, q);
            e >>= 1;
        }
        return result;
    };

    Vec x_poly(2, 0); x_poly[1] = 1;       // x

    // 计算 x^{q^i} mod f 的序列, i = 0..n
    std::vector<Vec> xqi_seq;
    xqi_seq.push_back(x_poly);             // i = 0: x
    for (int i = 1; i <= n; i++) {
        xqi_seq.push_back(qpow(xqi_seq.back()));
    }

    // 条件 ①: x^{q^n} ≡ x (mod f)
    if (!poly_equal(xqi_seq[n], x_poly)) return false;

    // 条件 ②: 对 n 的每个素因子 p,gcd(x^{q^{n/p}} - x, f) = 1
    auto factors = small_prime_factors(n);
    for (long p : factors) {
        int idx = n / (int)p;
        Vec diff = poly_sub(xqi_seq[idx], x_poly, q);
        if (poly_deg(diff) < 0) return false;  // x^{q^{n/p}} = x → 有真因子
        Vec g = poly_gcd(f, diff, q);
        if (poly_deg(g) > 0) return false;     // 有公因子 → 可约
    }
    return true;
}

// 生成一个随机首一(monic)n 次多项式
inline Vec random_monic_poly(int n, long q, std::mt19937_64& rng) {
    std::uniform_int_distribution<long> dist(0, q - 1);
    Vec f(n + 1, 0);
    for (int i = 0; i < n; i++) f[i] = dist(rng);
    f[n] = 1;  // monic
    return f;
}

// 找一个 F_q 上次数为 n 的不可约多项式
inline Vec find_irreducible(int n, long q, uint64_t seed = 0) {
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    for (int attempt = 0; attempt < 10000; attempt++) {
        Vec f = random_monic_poly(n, q, rng);
        if (is_irreducible(f, q)) return f;
    }
    throw std::runtime_error("failed to find irreducible polynomial");
}

/* ══════════════════════════════════════════════════
   §3  FRD 编码主函数
   ══════════════════════════════════════════════════ */

/**
 * FRD context: 持有不可约多项式,可以缓存重复使用
 * 在一个方案的 Setup 中固定一次,后续 Extract/Encrypt 都用同一个 f
 */
struct FRDContext {
    long q;            // 模数
    int  n;            // 维度
    Vec  f;            // F_q 上的不可约多项式,长度 n+1

    static FRDContext setup(int n, long q, uint64_t seed = 0) {
        FRDContext ctx;
        ctx.q = q;
        ctx.n = n;
        ctx.f = find_irreducible(n, q, seed);
        return ctx;
    }
};

/**
 * FRD 主函数: id ∈ Z_q^n  →  H ∈ Z_q^{n×n}
 *
 * H 的第 j 列 = (x^j · a(x)) mod f(x),其中 a(x) 是 id 对应的多项式
 *
 * 实现技巧:
 *   不需要每列都做一次完整的多项式乘法
 *   只需要先算 a(x) (= id 本身的系数向量),
 *   然后第 j 列 = 第 j-1 列乘 x mod f(x)
 *   "乘 x mod f" 是 O(n) 的:左移一位,如果新的最高项非零,减去合适倍数的 f
 */
inline Mat frd_encode(const FRDContext& ctx, const Vec& id) {
    if ((int)id.size() != ctx.n)
        throw std::invalid_argument("id length must equal n");

    int n = ctx.n;
    long q = ctx.q;
    const Vec& f = ctx.f;
    // f 是 monic n 次,f[n] = 1
    // 因此 x^n ≡ -(f_0 + f_1 x + ... + f_{n-1} x^{n-1})  mod f

    Mat H = make_zero_mat(n, n);

    // 第 0 列就是 a(x) 的系数
    Vec col(n, 0);
    for (int i = 0; i < n; i++) col[i] = mod_pos(id[i], q);
    for (int i = 0; i < n; i++) H[i][0] = col[i];

    // 后续列: col_j = x · col_{j-1} mod f(x)
    for (int j = 1; j < n; j++) {
        // 左移一位:新的 col 是 (0, col[0], col[1], ..., col[n-2]),
        // 而 col[n-1] 变成了 x^n 的系数 c
        long c = col[n - 1];
        for (int i = n - 1; i > 0; i--) col[i] = col[i - 1];
        col[0] = 0;
        // 减去 c · (x^n - 实际值) → 等价于加 c · (-(f 除掉首一项))
        // 即:col -= c · (f[0], f[1], ..., f[n-1])
        if (c != 0) {
            for (int i = 0; i < n; i++)
                col[i] = mod_pos(col[i] - c * f[i], q);
        }
        for (int i = 0; i < n; i++) H[i][j] = col[i];
    }

    return H;
}

} // namespace cryptolib
