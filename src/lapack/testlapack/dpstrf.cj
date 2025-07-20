package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dpstrfTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (uplo in [Upper, Lower]) {
        for (n in [0, 1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65, 127, 128, 129]) {
            for (lda in [max(1, n), n + 5]) {
                for (rank in [Int64(0.7 * Float64(n)), n]) {
                    if (!dpstrfTest(lFloat64, rnd, uplo, n, lda, rank)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func dpstrfTest(lFloat64: LFloat64, rnd: Random, uplo: Uplo, n: Int64, lda: Int64, rankWant: Int64): Bool {
    const tol = 1e-13

    let name = "n=${n},lda=${lda}"
    let nFloat64 = nFloat64Implementation()

    let a = Array<Float64>(n * lda, repeat: 0.0)
    for (_ in 0..rankWant) {
        let x = randomSlice(n, rnd)
        nFloat64.dsyr(uplo, n, 1.0, x, 1, a, lda)
    }

    let aFac = Array<Float64>(a.size, repeat: 0.0)
    copy(aFac, a)

    var piv = Array<Int64>(n, repeat: 0)
    for (i in 0..piv.size) {
        piv[i] = -1
    }

    let work = Array<Float64>(2 * n, repeat: 0.0)

    let (rank, ok) = lFloat64.dpstrf(uplo, n, aFac, lda, piv, -1.0, work)

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

func residualDpstrf(uplo: Uplo, n: Int64, a: Array<Float64>, aFac: Array<Float64>, lda: Int64, rank: Int64,
    piv: Array<Int64>): Float64 {
    let nFloat64 = nFloat64Implementation()
    let perm = zeros(n, n, n)
    if (uplo == Upper) {
        let u = aFac
        let ldu = lda
        for (i in rank..n) {
            for (j in i..n) {
                u[i * ldu + j] = 0.0
            }
        }
        let aRec = zeros(n, n, n)
        for (i in 0..n) {
            for (j in i..n) {
                aRec.data[i * aRec.stride + j] = u[i * ldu + j]
            }
        }
        nFloat64.dtrmm(Left, Upper, Trans, NonUnitDiag, n, n, 1.0, u, ldu, aRec.data, aRec.stride)
        for (i in 0..n) {
            for (j in 0..n) {
                if (piv[i] > piv[j]) {
                    continue
                }
                if (i <= j) {
                    perm.data[piv[i] * perm.stride + piv[j]] = aRec.data[i * aRec.stride + j]
                } else {
                    perm.data[piv[i] * perm.stride + piv[j]] = aRec.data[j * aRec.stride + i]
                }
            }
        }
        for (i in 0..n) {
            for (j in i..n) {
                perm.data[i * perm.stride + j] -= a[i * lda + j]
            }
        }
    } else {
        let l = aFac
        let ldl = lda
        for (i in rank..n) {
            for (j in rank..=i) {
                l[i * ldl + j] = 0.0
            }
        }
        let aRec = zeros(n, n, n)
        for (i in 0..n) {
            for (j in 0..=i) {
                aRec.data[i * aRec.stride + j] = l[i * ldl + j]
            }
        }
        nFloat64.dtrmm(Right, Lower, Trans, NonUnitDiag, n, n, 1.0, l, ldl, aRec.data, aRec.stride)
        for (i in 0..n) {
            for (j in 0..n) {
                if (piv[i] < piv[j]) {
                    continue
                }
                if (i >= j) {
                    perm.data[piv[i] * perm.stride + piv[j]] = aRec.data[i * aRec.stride + j]
                } else {
                    perm.data[piv[i] * perm.stride + piv[j]] = aRec.data[j * aRec.stride + i]
                }
            }
        }
        for (i in 0..n) {
            for (j in 0..=i) {
                perm.data[i * perm.stride + j] -= a[i * lda + j]
            }
        }
    }
    return dlansy(MaxColumnSum, uplo, n, perm.data, perm.stride) / Float64(n)
}
