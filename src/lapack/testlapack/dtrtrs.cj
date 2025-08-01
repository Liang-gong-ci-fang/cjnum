package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*
import cjnum.floats.*
import cjnum.util.*

public func dtrtrsTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)

    for (trans in [NoTrans, Trans, ConjTrans]) {
        for (uplo in [Upper, Lower]) {
            for (diag in [UnitDiag, NonUnitDiag]) {
                for (n in [0, 1, 2, 3, 4, 5, 10]) {
                    for (nrhs in [0, 1, 2, 3, 4, 5, 10, 15]) {
                        for (lda in [max(1, n), n + 3]) {
                            for (ldb in [max(1, nrhs), nrhs + 3]) {
                                if (diag == UnitDiag) {
                                    if (!dtrtrsTest(lFloat64, rnd, uplo, trans, diag, n, nrhs, lda, ldb, false)) {
                                        return false
                                    }
                                } else {
                                    if (!dtrtrsTest(lFloat64, rnd, uplo, trans, diag, n, nrhs, lda, ldb, true)) {
                                        return false
                                    }
                                    if (!dtrtrsTest(lFloat64, rnd, uplo, trans, diag, n, nrhs, lda, ldb, false)) {
                                        return false
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return true
}

func dtrtrsTest(lFloat64: LFloat64, rnd: Random, uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, nrhs: Int64, lda: Int64,
    ldb: Int64, singular: Bool): Bool {
    var singular_ = singular
    if (singular_ && diag == UnitDiag) {
        throw Exception("blas.Unit triangular matrix cannot be singular")
    }

    const tol = 1e-14

    if (n == 0) {
        singular_ = false
    }
    let name = "uplo=${uplo},diag=${diag},n=${n},nrhs=${nrhs},lda=${lda},ldb=${ldb},sing=${singular_}"

    let a = Array<Float64>(n * lda, repeat: 0.0)
    for (i in 0..a.size) {
        a[i] = rnd.nextGaussianFloat64()
    }
    if (singular_) {
        let i = rnd.nextInt64(n)
        a[i * lda + i] = 0.0
    }
    let aCopy = Array<Float64>(a.size, repeat: 0.0)
    copy(aCopy, a)

    let x = Array<Float64>(n * ldb, repeat: 0.0)
    for (i in 0..x.size) {
        x[i] = rnd.nextGaussianFloat64()
    }

    let b = Array<Float64>(x.size, repeat: 0.0)
    copy(b, x)
    let nFloat64 = nFloat64Implementation()
    nFloat64.dtrmm(Left, uplo, trans, diag, n, nrhs, 1.0, a, lda, b, ldb)

    let got = Array<Float64>(b.size, repeat: 0.0)
    copy(got, b)

    let ok = lFloat64.dtrtrs(uplo, trans, diag, n, nrhs, a, lda, got, ldb)

    if (!fEqual(a, aCopy)) {
        throw Exception("${name}: unexpected modification of A")
    }

    if (ok == singular_) {
        throw Exception("${name}: misdetected singular matrix, ok=${ok}")
    }

    if (!ok) {
        if (!fEqual(got, b)) {
            throw Exception("${name}: unexpected modification of B when singular")
        }
        return true
    }

    if (n == 0 || nrhs == 0) {
        return true
    }

    let work = Array<Float64>(n, repeat: 0.0)

    var aNorm: Float64
    if (trans == NoTrans) {
        aNorm = dlantr(MaxColumnSum, uplo, diag, n, n, a, lda, work)
    } else {
        aNorm = dlantr(MaxRowSum, uplo, diag, n, n, a, lda, work)
    }

    var resid: Float64 = 0.0
    for (j in 0..nrhs) {
        nFloat64.dcopy(n, got[j..], ldb, work, 1)
        nFloat64.dtrmv(uplo, trans, diag, n, a, lda, work, 1)
        nFloat64.daxpy(n, -1.0, b[j..], ldb, work, 1)
        let rjNorm = nFloat64.dasum(n, work, 1)
        let xNorm = nFloat64.dasum(n, got[j..], ldb)
        resid = nanMax(resid, rjNorm / aNorm / xNorm)
    }
    if (resid > tol) {
        throw Exception("${name}: unexpected result; resid=${resid},want<=${tol}")
    }
    return true
}
