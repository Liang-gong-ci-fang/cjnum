package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*
import cjnum.util.*

class dlasclTestStruct {
    dlasclTestStruct(
        let m: Int64,
        let n: Int64
    ) {}
}

public func dlasclTest(lFloat64: LFloat64): Bool {
    const tol = 1e-15

    let rnd = Random(1)
    var ti = 0
    for (test in [
            dlasclTestStruct(0, 0),
            dlasclTestStruct(1, 1),
            dlasclTestStruct(1, 10),
            dlasclTestStruct(10, 1),
            dlasclTestStruct(2, 2),
            dlasclTestStruct(2, 11),
            dlasclTestStruct(11, 2),
            dlasclTestStruct(3, 3),
            dlasclTestStruct(3, 11),
            dlasclTestStruct(11, 3),
            dlasclTestStruct(11, 11),
            dlasclTestStruct(11, 100),
            dlasclTestStruct(100, 11)
        ]) {
        let m = test.m
        let n = test.n
        for (extra in [0, 11]) {
            for (kind in [General, UpperTri, LowerTri]) {
                let a = randomGeneral(m, n, n + extra, rnd)
                let aCopy = cloneGeneral(a)
                let cfrom = rnd.nextGaussianFloat64()
                let cto = rnd.nextGaussianFloat64()
                let scale = cto / cfrom

                lFloat64.dlascl(kind, -1, -1, cfrom, cto, m, n, a.data, a.stride)

                let prefix = "Case #${ti}: kind=${kind},m=${m},n=${n},extra=${extra}"
                if (!generalOutsideAllNaN(a)) {
                    throw Exception("${prefix}: out-of-range write to A")
                }
                if (kind == UpperTri) { // TODO 验证
                    for (i in 0..m) {
                        for (j in 0..min(i, n)) {
                            if (a.data[i * a.stride + j] != aCopy.data[i * aCopy.stride + j]) {
                                throw Exception("${prefix}: unexpected modification in lower triangle of A")
                            }
                        }
                    }
                } else if (kind == LowerTri) {
                    for (i in 0..m) {
                        for (j in i + 1..n) {
                            if (a.data[i * a.stride + j] != aCopy.data[i * aCopy.stride + j]) {
                                throw Exception("${prefix}: unexpected modification in upper triangle of A")
                            }
                        }
                    }
                }

                var resid: Float64 = 0.0
                if (kind == General) {
                    for (i in 0..m) {
                        for (j in 0..n) {
                            let want = scale * aCopy.data[i * aCopy.stride + j]
                            let got = a.data[i * a.stride + j]
                            resid = nanMax(resid, abs(want - got))
                        }
                    }
                } else if (kind == UpperTri) {
                    for (i in 0..m) {
                        for (j in i..n) {
                            let want = scale * aCopy.data[i * aCopy.stride + j]
                            let got = a.data[i * a.stride + j]
                            resid = nanMax(resid, abs(want - got))
                        }
                    }
                } else if (kind == LowerTri) {
                    for (i in 0..m) {
                        for (j in 0..=min(i, n - 1)) {
                            let want = scale * aCopy.data[i * aCopy.stride + j]
                            let got = a.data[i * a.stride + j]
                            resid = nanMax(resid, abs(want - got))
                        }
                    }
                }
                if (resid > tol * Float64(max(m, n))) {
                    throw Exception("${prefix}: unexpected result; residual=${resid}, want<=${tol*Float64(max(m, n))}")
                }
            }
        }
        ti++
    }
    return true
}
