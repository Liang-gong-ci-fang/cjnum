package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.util.*
import std.math.*
import std.random.*

public func dlatbsTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 50]) {
        for (kd in [0, (n + 1) / 4, (3 * n - 1) / 4, (5 * n + 1) / 4]) {
            for (uplo in [Upper, Lower]) {
                for (trans in [NoTrans, Trans, ConjTrans]) {
                    for (ldab in [kd + 1, kd + 1 + 7]) {
                        for (kind in [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18]) {
                            if (!dlatbsTest(lFloat64, rnd, kind, uplo, trans, n, kd, ldab)) {
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

func dlatbsTest(lFloat64: LFloat64, rnd: Random, kind: Int64, uplo: Uplo, trans: Transpose, n: Int64, kd: Int64,
    ldab: Int64): Bool {
    const eps = 1e-15

    var ab = Array<Float64>()
    if (n > 0) {
        ab = Array<Float64>((n - 1) * ldab + kd + 1, repeat: 0.0)
    }
    for (i in 0..ab.size) {
        ab[i] = rnd.nextGaussianFloat64()
    }

    let (diag, b) = dlattb(kind, uplo, trans, n, kd, ab, ldab, rnd)

    let abCopy = Array<Float64>(ab.size, repeat: 0.0)
    copy(abCopy, ab)

    let cnorm = Array<Float64>(n, repeat: 0.0)
    for (i in 0..cnorm.size) {
        cnorm[i] = -1.0
    }

    let x = Array<Float64>(n, repeat: 0.0)
    copy(x, b)
    var scale = lFloat64.dlatbs(uplo, trans, diag, false, n, kd, ab, ldab, x, cnorm)

    let name = "kind=${kind},uplo=${uplo},trans=${trans},diag=${diag},n=${n},kd=${kd},ldab=${ldab}"

    if (!fEqual(ab, abCopy)) {
        throw Exception("${name}: unexpected modification of ab")
    }
    if (fCount({v: Float64 => v == -1.0}, cnorm) > 0) {
        throw Exception("${name}: expected modification of cnorm")
    }

    var resid = dlatbsResidual(uplo, trans, diag, n, kd, ab, ldab, scale, cnorm, b, x)
    if (resid >= eps) {
        throw Exception("${name}: unexpected result when normin=false. residual=${resid}")
    }

    let cnormCopy = Array<Float64>(cnorm.size, repeat: 0.0)
    copy(cnormCopy, cnorm)
    copy(x, b)
    scale = lFloat64.dlatbs(uplo, trans, diag, true, n, kd, ab, ldab, x, cnorm)

    if (!fEqualApprox(cnorm, cnormCopy, 1e-15)) {
        throw Exception("${name}: unexpected modification of cnorm")
    }

    resid = dlatbsResidual(uplo, trans, diag, n, kd, ab, ldab, scale, cnorm, b, x)
    if (resid >= eps) {
        throw Exception("${name}: unexpected result when normin=true. residual=${resid}")
    }
    return true
}

func dlatbsResidual(uplo: Uplo, trans: Transpose, diag: Diag, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64,
    scale: Float64, cnorm: Array<Float64>, b: Array<Float64>, x: Array<Float64>): Float64 {
    if (n == 0) {
        return 0.0
    }

    var tnorm = 0.0
    if (diag == NonUnitDiag) {
        if (uplo == Upper) {
            for (j in 0..n) {
                tnorm = nanMax(tnorm, abs(ab[j * ldab]) + cnorm[j])
            }
        } else {
            for (j in 0..n) {
                tnorm = nanMax(tnorm, abs(ab[j * ldab + kd]) + cnorm[j])
            }
        }
    } else {
        for (j in 0..n) {
            tnorm = nanMax(tnorm, 1.0 + cnorm[j])
        }
    }

    const eps = dlamchE
    const tiny = safmin

    let nFloat64 = nFloat64Implementation()

    var ix = nFloat64.idamax(n, x, 1)
    let xNorm = nanMax(1.0, abs(x[ix]))
    let xScal = (1.0 / xNorm) / Float64(kd + 1)

    let resid = Array<Float64>(x.size, repeat: 0.0)
    copy(resid, x)
    nFloat64.dscal(n, xScal, resid, 1)
    nFloat64.dtbmv(uplo, trans, diag, n, kd, ab, ldab, resid, 1)
    nFloat64.daxpy(n, -scale * xScal, b, 1, resid, 1)

    ix = nFloat64.idamax(n, resid, 1)
    var residNorm = abs(resid[ix])
    if (residNorm * tiny <= xNorm) {
        if (xNorm > 0.0) {
            residNorm /= xNorm
        }
    } else if (residNorm > 0.0) {
        residNorm = 1.0 / eps
    }
    if (residNorm * tiny <= tnorm) {
        if (tnorm > 0.0) {
            residNorm /= tnorm
        }
    } else if (residNorm > 0.0) {
        residNorm = 1.0 / eps
    }

    return residNorm
}
