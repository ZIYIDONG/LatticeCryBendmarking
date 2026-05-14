#pragma once
/**
 * 带模数切换的 Powersof2 — IBE 私钥生成
 *
 * 公式:
 *   s̄ = (1, -e) ∈ Z_q^{m+1}
 *   sk_id = Powersof2(s̄)
 *         = ( Powersof2_p(1),  -(p/q) · Powersof2_q(e) )
 *
 * 两段独立处理:
 *   ① 第一段: 标量 1 在模 p 下展开,长度 k_p = ⌈log₂ p⌉
 *   ② 第二段: 向量 e 在模 q 下展开,然后按 p/q 缩放并取整
 *
 * 最终拼接为长度 k_p + m·k_q 的整数向量(注意:它是定义在 Z_p 上的!)
 */

#include "powersof_plain-LWE.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>

namespace cryptolib {

/* ══════════════════════════════════════════════════
   §1  辅助: 有理数 round(p·x/q) — 模数切换的核心
   ══════════════════════════════════════════════════ */
/**
 * 计算 round(p · x / q),用 128-bit 中间值防溢出
 *
 * 在模数切换中,这个操作把 Z_q 中的元素"按比例缩放"到 Z_p,
 * 误差由四舍五入引入,理论上 ≤ 1/2
 *
 * 注意: 输入 x 应当先归约到中心区间 (-q/2, q/2],
 *       这样负数也能被正确缩放(否则 x = q-1 会被当成大正数处理)
 */
inline long round_scale(long x, long p, long q) {
    // 中心化到 (-q/2, q/2]
    long xc = ((x % q) + q) % q;
    if (xc > q / 2) xc -= q;

    // 用 __int128 防止 p · x 溢出
    __int128 num = (__int128)p * xc;
    // 四舍五入除法: round(num / q)
    long result;
    if (num >= 0) result = (long)((num + q / 2) / q);
    else          result = (long)(-(((-num) + q / 2) / q));
    return result;
}

/* ══════════════════════════════════════════════════
   §2  分两段的 Powersof2 (modulus switching 版本)
   ══════════════════════════════════════════════════ */
/**
 * Powersof2_with_modswitch(e, p, q)
 *
 * 输入:
 *   e ∈ Z_q^m         —— LWE 噪声/密钥向量
 *   p                 —— 小模数(明文模数)
 *   q                 —— 大模数(密文模数)
 *
 * 输出:
 *   sk_id = ( Powersof2_p(1),  -(p/q) · Powersof2_q(e) )
 *   长度: k_p + m · k_q
 *   元素值在 Z_p 范围内(整数,可能为负,通常需要再 mod p)
 *
 * 算法:
 *   ① 第一段: 直接对标量 1 调用 Powersof2_p,得到 (1, 2, 4, ..., 2^{k_p-1}) mod p
 *   ② 第二段:
 *      a. 对每个 e_i,在模 q 下计算 Powersof2_q(e_i) = (e_i, 2e_i, ..., 2^{k_q-1} e_i) mod q
 *      b. 对每个分量做 round_scale(·, p, q),即 round(p·val/q)
 *      c. 取负号
 *   ③ 拼接两段
 */
inline Vec powers_of_2_with_modswitch(
    const Vec& e,         // 长度 m 的向量
    long p,               // 小模数
    long q,               // 大模数
    bool reduce_mod_p = false) // 是否最后再 mod p
{
    if (p <= 1 || q <= p)
        throw std::invalid_argument(std::string("Need 1 < p < q (p=") + std::to_string(p) + ", q=" + std::to_string(q) + ")");

    int m   = (int)e.size();
    int k_p = compute_k(p, 2);    // 第一段长度
    int k_q = compute_k(q, 2);    // 每个 e_i 展开后的长度

    Vec sk;
    sk.reserve(k_p + m * k_q);

    /* ───── 第一段: Powersof2_p(1) ───── */
    // 对标量 1 应用 Powersof2,模数为 p
    Vec part1 = powers_of_b_scalar(1, 2, p);   // (1, 2, 4, ..., 2^{k_p-1}) mod p
    for (long v : part1) sk.push_back(v);

    /* ───── 第二段: -(p/q) · Powersof2_q(e) ───── */
    for (int i = 0; i < m; i++) {
        // 在模 q 下展开 e_i
        long val = ((e[i] % q) + q) % q;
        for (int j = 0; j < k_q; j++) {
            // val 此时是 2^j · e_i mod q
            // 模数切换: round(p · val / q)
            long scaled = round_scale(val, p, q);
            // 取负号
            long entry = -scaled;
            // 可选: 归约到 Z_p
            if (reduce_mod_p) entry = ((entry % p) + p) % p;
            sk.push_back(entry);
            // 下一项 = val · 2 mod q
            val = (val * 2) % q;
        }
    }

    return sk;
}

/* ══════════════════════════════════════════════════
   §3  构造 s̄ = (1, -e) 并应用 Powersof2 (上层封装)
   ══════════════════════════════════════════════════ */
/**
 * 完整对应图片公式:
 *   s̄    = (1, -e)
 *   sk_id = Powersof2(s̄) 拆成两段
 *
 * 这里直接吃 e (不需要先构造 -e),因为第二段的负号
 * 已经在 powers_of_2_with_modswitch 内部处理了
 */
inline Vec ibe_extract_key(const Vec& e, long p, long q) {
    return powers_of_2_with_modswitch(e, p, q, /*reduce_mod_p=*/true);
}

/* ══════════════════════════════════════════════════
   §4  解密侧的"对偶恒等式"验证辅助
   ══════════════════════════════════════════════════ */
/**
 * 在模数切换的设定下,对偶恒等式变成近似的:
 *   <BitDecomp_2(c̄), sk_id> ≈ (p/q) · <c̄, s̄>  (mod p)
 *
 * 其中 c̄ = (c_0, c_1, ..., c_m) ∈ Z_q^{m+1} 是密文
 * 误差来自 round_scale 引入的舍入(每次约 1/2)
 *
 * 这个函数模拟解密侧的内积计算
 */
inline long ibe_decrypt_inner(const Vec& c_bar, const Vec& sk_id, long p) {
    // BitDecomp_2(c̄) 在模 q 下,但分解后的位数字 ∈ {0,1},
    // 与 sk_id (定义在 Z_p) 做内积时直接整数相乘即可
    // 注意: 调用者需要先做 BitDecomp,这里只算内积部分
    if (c_bar.size() != sk_id.size())
        throw std::invalid_argument("size mismatch");
    long s = 0;
    for (size_t i = 0; i < c_bar.size(); i++)
        s = ((s + c_bar[i] * sk_id[i]) % p + p) % p;
    return s;
}

} // namespace cryptolib
