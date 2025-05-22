package cjnum.lapack.testlapack

import std.random.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64

public func dlarfxTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (side in [Right, Left]) {
        for (m in 1..12) {
            for (n in 1..12) {
                for (extra in [0, 1, 11]) {
                    for (cas in 0..10) {
                        if (!testDlarfx(lFloat64, side, m, n, extra, rnd)) {
                            return false
                        }
                    }
                }
            }
        }
    }
    return true
}

func testDlarfx(lFloat64: LFloat64, side: Side, m: Int64, n: Int64, extra: Int64, rnd: Random): Bool {
    const tol = 1e-13

    var v: Array<Float64>
    if (side == Left) {
        v = randomSlice(m, rnd)
    } else {
        v = randomSlice(n, rnd)
    }
    let tau = rnd.nextGaussianFloat64()
    let ldc = n + extra
    let c = randomGeneral(m, n, ldc, rnd)

    var h: General64
    if (side == Left) {
        h = eye(m, m + extra)
    } else {
        h = eye(n, n + extra)
    }
    ger(-tau, Vector(inc: 1, data: v), Vector(inc: 1, data: v), h)

    let cWant = nanGeneral(m, n, ldc)
    if (side == Left) {
        gemm(NoTrans, NoTrans, 1.0, h, c, 0.0, cWant)
    } else {
        gemm(NoTrans, NoTrans, 1.0, c, h, 0.0, cWant)
    }

    var work: Array<Float64> = Array<Float64>()
    if (h.rows > 10) {
        if (side == Left) {
            work = Array<Float64>(n, item: 0.0)
        } else {
            work = Array<Float64>(m, item: 0.0)
        }
    }

    lFloat64.dlarfx(side, m, n, v, tau, c.data, c.stride, work)

    let prefix = "Case side=${side}, m=${m}, n=${n}, extra=${extra}"

    if (!generalOutsideAllNaN(c)) {
        throw Exception("${prefix}: out-of-range write to C\n${c.data}")
    }

    if (!equalApproxGeneral(c, cWant, tol)) {
        throw Exception("${prefix}: unexpected C\n${c.data}")
    }
    return true
}
