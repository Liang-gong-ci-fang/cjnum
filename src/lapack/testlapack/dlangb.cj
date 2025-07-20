package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlangbTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (norm in [MaxAbs, MaxRowSum, MaxColumnSum, Frobenius]) {
        for (n in [0, 1, 2, 3, 4, 5, 10]) {
            for (m in [0, 1, 2, 3, 4, 5, 10]) {
                for (kl in [0, 1, 2, 3, 4, 5, 10]) {
                    for (ku in [0, 1, 2, 3, 4, 5, 10]) {
                        for (ldab in [kl + ku + 1, kl + ku + 1 + 7]) {
                            if (!dlangbTest(lFloat64, rnd, norm, m, n, kl, ku, ldab)) {
                                return false
                            }
                        }
                    }
                }
            }
        }
    }
    return true
}

func dlangbTest(lFloat64: LFloat64, rnd: Random, norm: MatrixNorm, m: Int64, n: Int64, kl: Int64, ku: Int64, ldab: Int64): Bool {
    const tol = 1e-14

    let name = "m=${m},n=${n},kl=${kl},ku=${ku},ldab=${ldab}"

    let ab = randomSlice(m * ldab, rnd)
    if (m > 0 && n > 0 && rnd.nextFloat64() < 0.5) {
        let i = rnd.nextInt64(m)
        ab[i * ldab + kl] = Float64.NaN
    }
    let abCopy = Array<Float64>(ab.size, repeat: 0.0)
    copy(abCopy, ab)

    if (m == 0 || n == 0) {
        var got = lFloat64.dlangb(norm, m, n, kl, ku, [], ldab)
        if (got != 0.0) {
            throw Exception("${name}: unexpected result for zero-sized matrix with nil input")
        }
        got = lFloat64.dlangb(norm, m, n, kl, ku, ab, ldab)
        if (!fSame(ab, abCopy)) {
            throw Exception("${name}: unexpected modification in dl")
        }
        if (got != 0.0) {
            throw Exception("${name}: unexpected result for zero-sized matrix with non-nil input")
        }
        return true
    }

    let got = lFloat64.dlangb(norm, m, n, kl, ku, ab, ldab)

    if (!fSame(ab, abCopy)) {
        throw Exception("${name}: unexpected modification in ab")
    }

    let a = zeros(m, n, n)
    for (i in 0..m) {
        for (j in max(0, i - kl)..min(i + ku + 1, n)) {
            a.data[i * a.stride + j] = ab[i * ldab + j - i + kl]
        }
    }
    let want = dlange(norm, a.rows, a.cols, a.data, a.stride)

    if (want.isNaN()) {
        if (!got.isNaN()) {
            throw Exception("${name}: unexpected result with NaN element; got ${got}, want ${want}")
        }
        return true
    }

    if (got.isNaN()) {
        throw Exception("${name}: unexpected NaN; want ${want}")
    }

    if (norm == MaxAbs) {
        if (got != want) {
            throw Exception("${name}: unexpected result; got ${got}, want ${want}")
        }
        return true
    }
    let diff = abs(got - want)
    if (diff > tol) {
        throw Exception("${name}: unexpected result; got ${got}, want ${want}, diff=${diff}")
    }
    return true
}
