package cjnum.lapack.testlapack

import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dormlqTestStruct {
    dormlqTestStruct(
        let common: Int64,
        let adim: Int64,
        let cdim: Int64,
        let lda: Int64,
        let ldc: Int64
    ) {}
}

public func dormlqTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (side in [Left, Right]) {
        for (trans in [NoTrans, Trans]) {
            for (wl in [minimumWork, mediumWork, optimumWork]) {
                for (test in [
                        dormlqTestStruct(0, 0, 0, 0, 0),
                        dormlqTestStruct(6, 7, 8, 0, 0),
                        dormlqTestStruct(6, 8, 7, 0, 0),
                        dormlqTestStruct(7, 6, 8, 0, 0),
                        dormlqTestStruct(7, 8, 6, 0, 0),
                        dormlqTestStruct(8, 6, 7, 0, 0),
                        dormlqTestStruct(8, 7, 6, 0, 0),
                        dormlqTestStruct(100, 200, 300, 0, 0),
                        dormlqTestStruct(100, 300, 200, 0, 0),
                        dormlqTestStruct(200, 100, 300, 0, 0),
                        dormlqTestStruct(200, 300, 100, 0, 0),
                        dormlqTestStruct(300, 100, 200, 0, 0),
                        dormlqTestStruct(300, 200, 100, 0, 0),
                        dormlqTestStruct(100, 200, 300, 400, 500),
                        dormlqTestStruct(100, 300, 200, 400, 500),
                        dormlqTestStruct(200, 100, 300, 400, 500),
                        dormlqTestStruct(200, 300, 100, 400, 500),
                        dormlqTestStruct(300, 100, 200, 400, 500),
                        dormlqTestStruct(300, 200, 100, 400, 500),
                        dormlqTestStruct(100, 200, 300, 500, 400),
                        dormlqTestStruct(100, 300, 200, 500, 400),
                        dormlqTestStruct(200, 100, 300, 500, 400),
                        dormlqTestStruct(200, 300, 100, 500, 400),
                        dormlqTestStruct(300, 100, 200, 500, 400),
                        dormlqTestStruct(300, 200, 100, 500, 400)
                    ]) {
                    var ma: Int64
                    var na: Int64
                    var mc: Int64
                    var nc: Int64
                    if (side == Left) {
                        ma = test.adim
                        na = test.common
                        mc = test.common
                        nc = test.cdim
                    } else {
                        ma = test.adim
                        na = test.common
                        mc = test.cdim
                        nc = test.common
                    }
                    var lda = test.lda
                    if (lda == 0) {
                        lda = max(1, na)
                    }
                    let a = Array<Float64>(ma * lda, {_ => rnd.nextFloat64()})

                    var ldc = test.ldc
                    if (ldc == 0) {
                        ldc = nc
                    }
                    let c = Array<Float64>(mc * ldc, {_ => rnd.nextFloat64()})

                    let k = min(ma, na)
                    let tau = Array<Float64>(k, repeat: 0.0)
                    var work = Array<Float64>(1, repeat: 0.0)
                    lFloat64.dgelqf(ma, na, a, lda, tau, work, -1)
                    work = Array<Float64>(Int64(work[0]), repeat: 0.0)
                    lFloat64.dgelqf(ma, na, a, lda, tau, work, work.size)

                    let cCopy = Array<Float64>(c.size, repeat: 0.0)
                    copy(cCopy, c)
                    let ans = Array<Float64>(c.size, repeat: 0.0)
                    copy(ans, cCopy)

                    var nw: Int64
                    if (side == Left) {
                        nw = nc
                    } else {
                        nw = mc
                    }
                    work = Array<Float64>(max(1, nw), repeat: 0.0)
                    lFloat64.dorml2(side, trans, mc, nc, k, a, lda, tau, ans, ldc, work)

                    var lwork: Int64 = 0
                    if (wl == minimumWork) {
                        lwork = nw
                    } else if (wl == optimumWork) {
                        lFloat64.dormlq(side, trans, mc, nc, k, a, lda, tau, c, ldc, work, -1)
                        lwork = Int64(work[0])
                    } else if (wl == mediumWork) {
                        work = Array<Float64>(1, repeat: 0.0)
                        lFloat64.dormlq(side, trans, mc, nc, k, a, lda, tau, c, ldc, work, -1)
                        lwork = (Int64(work[0]) + nw) / 2
                    }
                    lwork = max(1, lwork)
                    work = Array<Float64>(lwork, repeat: 0.0)

                    lFloat64.dormlq(side, trans, mc, nc, k, a, lda, tau, c, ldc, work, lwork)
                    if (!fEqualApprox(c, ans, 1e-13)) {
                        throw Exception("Dormqr and Dorm2r results mismatch")
                    }
                }
            }
        }
    }
    return true
}
