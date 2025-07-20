package cjnum.lapack.testlapack

import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dormqrTestStruct {
    dormqrTestStruct(
        let common: Int64,
        let adim: Int64,
        let cdim: Int64,
        let lda: Int64,
        let ldc: Int64
    ) {}
}

public func dormqrTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (side in [Left, Right]) {
        for (trans in [NoTrans, Trans]) {
            for (test in [
                    dormqrTestStruct(6, 7, 8, 0, 0),
                    dormqrTestStruct(6, 8, 7, 0, 0),
                    dormqrTestStruct(7, 6, 8, 0, 0),
                    dormqrTestStruct(7, 8, 6, 0, 0),
                    dormqrTestStruct(8, 6, 7, 0, 0),
                    dormqrTestStruct(8, 7, 6, 0, 0),
                    dormqrTestStruct(100, 200, 300, 0, 0),
                    dormqrTestStruct(100, 300, 200, 0, 0),
                    dormqrTestStruct(200, 100, 300, 0, 0),
                    dormqrTestStruct(200, 300, 100, 0, 0),
                    dormqrTestStruct(300, 100, 200, 0, 0),
                    dormqrTestStruct(300, 200, 100, 0, 0),
                    dormqrTestStruct(100, 200, 300, 400, 500),
                    dormqrTestStruct(100, 300, 200, 400, 500),
                    dormqrTestStruct(200, 100, 300, 400, 500),
                    dormqrTestStruct(200, 300, 100, 400, 500),
                    dormqrTestStruct(300, 100, 200, 400, 500),
                    dormqrTestStruct(300, 200, 100, 400, 500),
                    dormqrTestStruct(100, 200, 300, 500, 400),
                    dormqrTestStruct(100, 300, 200, 500, 400),
                    dormqrTestStruct(200, 100, 300, 500, 400),
                    dormqrTestStruct(200, 300, 100, 500, 400),
                    dormqrTestStruct(300, 100, 200, 500, 400),
                    dormqrTestStruct(300, 200, 100, 500, 400)
                ]) {
                var ma: Int64
                var na: Int64
                var mc: Int64
                var nc: Int64
                if (side == Left) {
                    ma = test.common
                    na = test.adim
                    mc = test.common
                    nc = test.cdim
                } else {
                    ma = test.common
                    na = test.adim
                    mc = test.cdim
                    nc = test.common
                }

                var lda = test.lda
                if (lda == 0) {
                    lda = na
                }
                let a = Array<Float64>(ma * lda, repeat: 0.0)
                for (i in 0..a.size) {
                    a[i] = rnd.nextFloat64()
                }

                var ldc = test.ldc
                if (ldc == 0) {
                    ldc = nc
                }
                let c = Array<Float64>(mc * ldc, repeat: 0.0)
                for (i in 0..c.size) {
                    c[i] = rnd.nextFloat64()
                }

                let k = min(ma, na)
                let tau = Array<Float64>(k, repeat: 0.0)
                var work = Array<Float64>(1, repeat: 0.0)
                lFloat64.dgeqrf(ma, na, a, lda, tau, work, -1)
                work = Array<Float64>(Int64(work[0]), repeat: 0.0)
                lFloat64.dgeqrf(ma, na, a, lda, tau, work, work.size)

                let cCopy = Array<Float64>(c.size, repeat: 0.0)
                copy(cCopy, c)
                let ans = Array<Float64>(c.size, repeat: 0.0)
                copy(ans, cCopy)

                if (side == Left) {
                    work = Array<Float64>(nc, repeat: 0.0)
                } else {
                    work = Array<Float64>(mc, repeat: 0.0)
                }
                lFloat64.dorm2r(side, trans, mc, nc, k, a, lda, tau, ans, ldc, work)

                for (i in 0..work.size) {
                    work[i] = rnd.nextFloat64()
                }
                copy(c, cCopy)
                lFloat64.dormqr(side, trans, mc, nc, k, a, lda, tau, c, ldc, work, work.size)
                if (!fEqualApprox(c, ans, 1e-12)) {
                    throw Exception("Dormqr and Dorm2r mismatch for small work")
                }

                copy(c, cCopy)
                lFloat64.dormqr(side, trans, mc, nc, k, a, lda, tau, c, ldc, work, -1)
                work = Array<Float64>(Int64(work[0]), repeat: 0.0)
                for (i in 0..work.size) {
                    work[i] = rnd.nextFloat64()
                }
                lFloat64.dormqr(side, trans, mc, nc, k, a, lda, tau, c, ldc, work, work.size)
                if (!fEqualApprox(c, ans, 1e-12)) {
                    throw Exception("Dormqr and Dorm2r mismatch for full work")
                    // println("ccopy")
                    // for (i in 0..mc) {
                    //     println(cCopy[i * ldc..(i + 1) * ldc])
                    // }
                    // println("ans =")
                    // for (i in 0..mc) {
                    //     println(ans[i * ldc..(i + 1) * ldc])
                    // }
                    // println("c =")
                    // for (i in 0..mc) {
                    //     println(c[i * ldc..(i + 1) * ldc])
                    // }
                }

                copy(c, cCopy)
                if (side == Left) {
                    work = Array<Float64>(3 * nc, repeat: 0.0)
                } else {
                    work = Array<Float64>(3 * mc, repeat: 0.0)
                }
                lFloat64.dormqr(side, trans, mc, nc, k, a, lda, tau, c, ldc, work, work.size)
                if (!fEqualApprox(c, ans, 1e-12)) {
                    throw Exception("Dormqr and Dorm2r mismatch for medium work")
                }
            }
        }
    }
    return true
}
