package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.floats.*

public func dgeconTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 50]) {
        for (lda in [max(1, n), n + 3]) {
            _dgeconTest(lFloat64, rnd, n, lda)
        }
    }
    return true
}

func _dgeconTest(lFloat64: LFloat64, rnd: Random, n: Int64, lda: Int64): Unit {
    const ratioThresh = 10.0

    let a = Array<Float64>(max(0, (n - 1) * lda + n), { _ => 2.0 * rnd.nextFloat64() - 1.0 })
    let iwork = Array<Int64>(n, item: 0)
    let work = Array<Float64>(max(1, 4 * n), item: 0.0)

    let aFac = a.clone()
    let ipiv = Array<Int64>(n, item: 0)
    var ok = lFloat64.dgetrf(n, n, aFac, lda, ipiv)
    if (!ok) {
        throw Exception("n=${n},lda=${lda}: bad matrix, Dgetrf failed")
    }
    let aFacCopy = aFac.clone()

    let aInv = aFac.clone()
    ok = lFloat64.dgetri(n, aInv, lda, ipiv, work, work.size)
    if (!ok) {
        throw Exception("n=${n},lda=${lda}: bad matrix, Dgetri failed")
    }

    for (norm in [MaxColumnSum, MaxRowSum]) {
        let name = "norm=${norm.toString()},n=${n},lda=${lda}"

        let aNorm = lFloat64.dlange(norm, n, n, a, lda, work)
        let aInvNorm = lFloat64.dlange(norm, n, n, aInv, lda, work)

        let rcondWant = if (aNorm > 0.0 && aInvNorm > 0.0) {
            1.0 / aNorm / aInvNorm
        } else {
            1.0
        }

        let rcondGot = lFloat64.dgecon(norm, n, aFac, lda, aNorm, work, iwork)
        if (!fEqual(aFac, aFacCopy)) {
            throw Exception("${name}: unexpected modification of aFac")
        }

        let ratio = rCondTestRatio(rcondGot, rcondWant)
        if (ratio >= ratioThresh) {
            throw Exception("${name}: unexpected value of rcond; got=${rcondGot}, want=${rcondWant} (ratio=${ratio})")
        }

        for (anorm in [0.0, Float64.Inf, Float64.NaN]) {
            let rcondGot = lFloat64.dgecon(norm, n, aFac, lda, anorm, work, iwork)
            if (n == 0) {
                if (rcondGot != 1.0) {
                    throw Exception("${name}: unexpected rcond when anorm=${anorm}: got=${rcondGot}, want=1")
                }
                continue
            }
            if (anorm.isNaN()) {
                if (!rcondGot.isNaN()) {
                    throw Exception("${name}: NaN not propagated when anorm=NaN: got=${rcondGot}")
                }
                continue
            }
            if (rcondGot != 0.0) {
                throw Exception("${name}: unexpected rcond when anorm=${anorm}: got=%${rcondGot}, want=0")
            }
        }
    }
}