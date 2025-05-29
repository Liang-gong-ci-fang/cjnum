package cjnum.lapack.testlapack

import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*
import std.random.*

public func dlasy2Test(lFloat64: LFloat64) {
    let rnd = Random(1)
    for (tranl in [true, false]) {
        for (tranr in [true, false]) {
            for (isgn in [1, -1]) {
                for (n1 in [0, 1, 2]) {
                    for (n2 in [0, 1, 2]) {
                        for (extra in [0, 3]) {
                            for (cas in 0..100) {
                                var big: Bool = false
                                if (cas % 2 == 0) {
                                    big = true
                                }
                                if (!testDlasy2(lFloat64, tranl, tranr, isgn, n1, n2, extra, big, rnd)) {
                                    return false
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

func testDlasy2(lFloat64: LFloat64, tranl: Bool, tranr: Bool, isgn: Int64, n1: Int64, n2: Int64, extra: Int64, big: Bool,
    rnd: Random): Bool {
    const tol = 1e-14

    let name = "Case n1=${n1}, n2=${n2}, isgn=${isgn}, big=${big}"

    let tl = randomGeneral(n1, n1, n1 + extra, rnd)
    let tr = randomGeneral(n2, n2, n2 + extra, rnd)
    let x = randomGeneral(n1, n2, n2 + extra, rnd)
    let b = randomGeneral(n1, n2, n2 + extra, rnd)
    if (big) {
        for (i in 0..n1) {
            for (j in 0..n2) {
                b.data[i * b.stride + j] *= bignum
            }
        }
    }

    let tlCopy = cloneGeneral(tl)
    let trCopy = cloneGeneral(tr)
    let bCopy = cloneGeneral(b)

    let (scale, xnorm, ok) = lFloat64.dlasy2(tranl, tranr, isgn, n1, n2, tl.data, tl.stride, tr.data, tr.stride, b.data,
        b.stride, x.data, x.stride)

    if (!equalGeneral(tl, tlCopy)) {
        throw Exception("${name}: unexpected modification in TL")
    }
    if (!equalGeneral(tr, trCopy)) {
        throw Exception("${name}: unexpected modification in TR")
    }
    if (!equalGeneral(b, bCopy)) {
        throw Exception("${name}: unexpected modification in B")
    }

    if (!generalOutsideAllNaN(x)) {
        throw Exception("${name}: out-of-range write to x\n${x.data}")
    }

    if (n1 == 0 || n2 == 0) {
        return true
    }

    if (scale <= 0.0 || 1.0 < scale) {
        throw Exception("${name}: invalid value of scale, want in (0,1], got ${scale}")
    }

    let xnormWant = dlange(MaxRowSum, x.rows, x.cols, x.data, x.stride)
    if (xnormWant != xnorm) {
        throw Exception("${name}: unexpected xnorm: want ${xnormWant}, got ${xnorm}")
    }

    if (!ok) {
        // t.Logf("%v: Dlasy2 returned ok=false", name)  //TODO 待修改
        return true
    }

    let diff = zeros(n1, n2, n2)
    if (tranl) {
        gemm(Trans, NoTrans, 1.0, tl, x, 0.0, diff)
    } else {
        gemm(NoTrans, NoTrans, 1.0, tl, x, 0.0, diff)
    }
    if (tranr) {
        gemm(NoTrans, Trans, Float64(isgn), x, tr, 1.0, diff)
    } else {
        gemm(NoTrans, NoTrans, Float64(isgn), x, tr, 1.0, diff)
    }
    for (i in 0..n1) {
        for (j in 0..n2) {
            diff.data[i * diff.stride + j] -= scale * b.data[i * b.stride + j]
        }
    }

    let resid = dlange(MaxColumnSum, n1, n2, diff.data, diff.stride) / xnorm
    if (resid > tol) {
        throw Exception("${name}: unexpected result, resid=${resid}, want<=${tol}")
    }
    return true
}
