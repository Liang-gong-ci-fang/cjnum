package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.complex.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlag2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (lda in [2, 5]) {
        for (ldb in [2, 5]) {
            for (aKind in 0..=20) {
                for (bKind in 0..=20) {
                    if (!dlag2Test(lFloat64, rnd, lda, ldb, aKind, bKind)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func dlag2Test(lFloat64: LFloat64, rnd: Random, lda: Int64, ldb: Int64, aKind: Int64, bKind: Int64): Bool {
    const tol = 1e-14

    let a = makeDlag2TestMatrix(rnd, lda, aKind)
    let b = makeDlag2TestMatrix(rnd, ldb, bKind)

    let aCopy = cloneGeneral(a)
    let bCopy = cloneGeneral(b)

    let (scale1, scale2, wr1, wr2, wi) = lFloat64.dlag2(a.data, a.stride, b.data, b.stride)

    let name = "lda=${lda},ldb=${ldb},aKind=${aKind},bKind=${bKind}"
    let aStr = "A = [${a.data[0]},${a.data[1]}]\n    [${a.data[a.stride]},${a.data[a.stride+1]}"
    let bStr = "B = [${b.data[0]},${b.data[1]}]\n    [0.0,${b.data[b.stride+1]}"

    if (!fSame(a.data, aCopy.data)) {
        throw Exception("${name}: unexpected modification of a")
    }
    if (!fSame(b.data, bCopy.data)) {
        throw Exception("${name}: unexpected modification of b")
    }

    if (wi < 0.0) {
        throw Exception("${name}: wi is negative; wi=${wi},\n${aStr}\n${bStr}")
    }

    if (wi > 0.0) {
        if (wr1 != wr2) {
            throw Exception("${name}: complex eigenvalue but wr1 != wr2; wr1=${wr1}, wr2=${wr2},\n${aStr}\n${bStr}")
        }
        if (scale1 != scale2) {
            throw Exception(
                "${name}: complex eigenvalue but scale1 != scale2; scale1=${scale1}, scale2=${scale2},\n${aStr}\n${bStr}"
            )
        }
    }

    var (resid, err) = residualDlag2(a, b, scale1, Complex128(wr1, wi))
    match (err) {
        case Some(T) => throw Exception("${name}: invalid input data: ${err}\n${aStr}\n${bStr}")
        case None => ()
    }
    if (resid > tol || resid.isNaN()) {
        throw Exception(
            "${name}: unexpected first eigenvalue ${Complex128(wr1, wi)} with s=${scale1}; resid=${resid}, want<=${tol}\n${aStr}\n${bStr}"
        )
    }

    (resid, err) = residualDlag2(a, b, scale2, Complex128(wr2, -wi))
    match (err) {
        case Some(T) => throw Exception("${name}: invalid input data: ${err}\n${aStr}\n${bStr}")
        case None => ()
    }
    if (resid > tol || resid.isNaN()) {
        throw Exception(
            "${name}: unexpected second eigenvalue ${Complex128(wr2, -wi)} with s=${scale2}; resid=${resid}, want<=${tol}\n${aStr}\n${bStr}"
        )
    }
    return true
}

func makeDlag2TestMatrix(rnd: Random, ld: Int64, kind: Int64): General64 {
    var a = zeros(2, 2, ld)
    if (kind == 0) {
    } else if (kind == 1) {
        a.data[0] = 1.0
        a.data[a.stride + 1] = 1.0
    } else if (kind == 2) {
        a.data[0] = 2.0 * safmax
        a.data[a.stride + 1] = 2.0 * safmax
    } else if (kind == 3) {
        a.data[0] = safmin
        a.data[a.stride + 1] = safmin
    } else if (kind == 4) {
        a.data[0] = safmin
        a.data[a.stride + 1] = safmax
    } else if (kind == 5) {
        a.data[0] = safmax
        a.data[a.stride + 1] = safmin
    } else if (kind == 6) {
        a.data[0] = safmax
        a.data[1] = safmax
        a.data[a.stride] = -safmax
        a.data[a.stride + 1] = safmax
    } else if (kind == 7) {
        a.data[0] = safmin
        a.data[1] = safmin
        a.data[a.stride] = -safmin
        a.data[a.stride + 1] = safmin
    } else if (kind == 8) {
        a.data[0] = safmax * (2.0 * rnd.nextFloat64() - 1.0)
        a.data[1] = safmax * (2.0 * rnd.nextFloat64() - 1.0)
        a.data[a.stride] = safmax * (2.0 * rnd.nextFloat64() - 1.0)
        a.data[a.stride + 1] = safmax * (2.0 * rnd.nextFloat64() - 1.0)
    } else if (kind == 9) {
        a.data[0] = safmin * (2.0 * rnd.nextFloat64() - 1.0)
        a.data[1] = safmin * (2.0 * rnd.nextFloat64() - 1.0)
        a.data[a.stride] = safmin * (2.0 * rnd.nextFloat64() - 1.0)
        a.data[a.stride + 1] = safmin * (2.0 * rnd.nextFloat64() - 1.0)
    } else {
        a = randomGeneral(2, 2, ld, rnd)
    }
    return a
}

func residualDlag2(a: General64, b: General64, s: Float64, w: Complex128): (Float64, ?String) {
    var s_ = s
    var w_ = w
    const ulp = dlamchP

    let (a11, a12) = (a.data[0], a.data[1])
    let (a21, a22) = (a.data[a.stride], a.data[a.stride + 1])

    let (b11, b12) = (b.data[0], b.data[1])
    let b22 = b.data[b.stride + 1]

    var absw = zabs(w_)
    var anorm = max(abs(a11) + abs(a21), abs(a12) + abs(a22))
    anorm = max(anorm, safmin)
    var bnorm = max(abs(b11), abs(b12) + abs(b22))
    bnorm = max(bnorm, safmin)

    var temp = (safmin * anorm) * s_ + (safmin * bnorm) * absw
    if (temp >= 1.0) {
        s_ /= temp
        w_ = scale(1.0 / temp, w_)
        absw = zabs(w)
    }

    var s1 = max(ulp * max(s * anorm, absw * bnorm), safmin * max(s, absw))
    if (s1 < safmin) {
        if (s < safmin && absw < safmin) {
            return (1.0 / ulp, "ulp*max(s*|A|,|w|*|B|) < safmin and s and w could not be scaled; s=${s_}, |w|=${absw}")
        }
        temp = 1.0 / max(s * anorm + absw * bnorm, safmin)
        s_ *= temp
        w_ = scale(temp, w_)
        absw = zabs(w_)
        s1 = max(ulp * max(s * anorm, absw * bnorm), safmin * max(s, absw))
        if (s1 < safmin) {
            return (1.0 / ulp, "ulp*max(s*|A|,|w|*|B|) < safmin and s and w could not be scaled; s=${s_}, |w|=${absw}")
        }
    }

    let c11 = Complex128(s * a11, 0.0) - w * Complex128(b11, 0.0)
    let c12 = Complex128(s * a12, 0.0) - w * Complex128(b12, 0.0)
    let c21 = Complex128(s * a21, 0.0)
    let c22 = Complex128(s * a22, 0.0) - w * Complex128(b22, 0.0)
    let cnorm = max(zabs(c11) + zabs(c21), zabs(c12) + zabs(c22))
    let cs = 1.0 / sqrt(max(cnorm, safmin))
    let det = cmplxdet2x2(scale(cs, c11), scale(cs, c12), scale(cs, c21), scale(cs, c22))
    return (zabs(det) / s1 * ulp, None)
}

func zabs(z: Complex128): Float64 {
    return abs(z.real) + abs(z.imag)
}

func scale(f: Float64, c: Complex128): Complex128 {
    return Complex128(f * c.real, f * c.imag)
}

func cmplxdet2x2(a11: Complex128, a12: Complex128, a21: Complex128, a22: Complex128): Complex128 {
    return a11 * a22 - a12 * a21
}
