package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlangtTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (norm in [MaxAbs, MaxRowSum, MaxColumnSum, Frobenius]) {
        for (n in [0, 1, 2, 3, 4, 5, 10]) {
            for (_ in 0..10) {
                if (!dlangtTest(lFloat64, rnd, norm, n)) {
                    return false
                }
            }
        }
    }
    return true
}

func dlangtTest(lFloat64: LFloat64, rnd: Random, norm: MatrixNorm, n: Int64): Bool {
    const tol = 1e-14
    const extra = 10

    let name = "n=${n}"

    let dl = randomSlice(n + extra, rnd)
    let dlCopy = Array<Float64>(dl.size, item: 0.0)
    copy(dlCopy, dl)

    let d = randomSlice(n + 1 + extra, rnd)
    if (n > 0 && rnd.nextFloat64() < 0.5) {
        d[rnd.nextInt64(n)] = Float64.NaN
    }
    let dCopy = Array<Float64>(d.size, item: 0.0)
    copy(dCopy, d)

    let du = randomSlice(n + extra, rnd)
    let duCopy = Array<Float64>(du.size, item: 0.0)
    copy(duCopy, du)

    if (n == 0) {
        var got = lFloat64.dlangt(norm, n, [], [], [])
        if (got != 0.0) {
            throw Exception("${name}: unexpected result for zero-sized matrix with nil input")
        }
        got = lFloat64.dlangt(norm, n, dl, d, du)
        if (!fSame(dl, dlCopy)) {
            throw Exception("${name}: unexpected modification in dl")
        }
        if (!fSame(d, dCopy)) {
            throw Exception("${name}: unexpected modification in d")
        }
        if (!fSame(du, duCopy)) {
            throw Exception("${name}: unexpected modification in du")
        }
        if (got != 0.0) {
            throw Exception("${name}: unexpected result for zero-sized matrix with non-nil input")
        }
        return true
    }

    let a = zeros(n, n, n)
    for (i in 0..n - 1) {
        a.data[i * a.stride + i] = d[i]
        a.data[i * a.stride + i + 1] = du[i]
        a.data[(i + 1) * a.stride + i] = dl[i]
    }
    a.data[(n - 1) * a.stride + n - 1] = d[n - 1]

    let got = lFloat64.dlangt(norm, n, dl, d, du)

    if (!fSame(dl, dlCopy)) {
        throw Exception("${name}: unexpected modification in dl")
    }
    if (!fSame(d, dCopy)) {
        throw Exception("${name}: unexpected modification in d")
    }
    if (!fSame(du, duCopy)) {
        throw Exception("${name}: unexpected modification in du")
    }

    let want = dlange(norm, n, n, a.data, a.stride)

    if (want.isNaN()) {
        if (!got.isNaN()) {
            throw Exception("${name}: unexpected result with NaN element; got ${got}, want ${want}")
        }
        return true
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
