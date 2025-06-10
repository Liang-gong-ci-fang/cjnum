package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dpstf2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (uplo in [Upper, Lower]) {
        for (n in [0, 1, 2, 3, 4, 5, 10, 20, 50]) {
            for (lda in [max(1, n), n + 5]) {
                for (rank in [Int64(0.7 * Float64(n)), n]) {
                    if (!dpstf2Test(lFloat64, rnd, uplo, n, lda, rank)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func dpstf2Test(lFloat64: LFloat64, rnd: Random, uplo: Uplo, n: Int64, lda: Int64, rankWant: Int64): Bool {
    const tol = 1e-14

    let name = "n=${n},lda=${lda}"
    let nFloat64 = nFloat64Implementation()

    let a = Array<Float64>(n * lda, item: 0.0)
    for (_ in 0..rankWant) {
        let x = randomSlice(n, rnd)
        nFloat64.dsyr(uplo, n, 1.0, x, 1, a, lda)
    }

    let aFac = Array<Float64>(a.size, item: 0.0)
    copy(aFac, a)

    let piv = Array<Int64>(n, item: 0)
    for (i in 0..piv.size) {
        piv[i] = -1
    }

    let work = Array<Float64>(2 * n, item: 0.0)

    let (rank, ok) = lFloat64.dpstf2(uplo, n, aFac, lda, piv, -1.0, work)

    if (ok != (rank == n)) {
        throw Exception("${name}: unexpected ok; got ${ok}, want ${rank == n}")
    }
    if (rank != rankWant) {
        throw Exception("${name}: unexpected rank; got ${rank}, want ${rankWant}")
    }

    if (n == 0) {
        return true
    }

    let resid = residualDpstrf(uplo, n, a, aFac, lda, rank, piv)
    if (resid > tol || resid.isNaN()) {
        throw Exception("${name}: residual too large; got ${resid}, want<=${tol}")
    }
    return true
}
