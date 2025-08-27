package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dlarfbTestStruct {
    dlarfbTestStruct(
        let ma: Int64,
        let na: Int64,
        let cdim: Int64,
        let lda: Int64,
        let ldt: Int64,
        let ldc: Int64
    ) {}
}

public func dlarfbTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (store in [ColumnWise, RowWise]) {
        for (direct in [Forward, Backward]) {
            for (side in [Left, Right]) {
                for (trans in [Trans, NoTrans]) {
                    var cas = 0
                    for (test in [
                            dlarfbTestStruct(6, 6, 6, 0, 0, 0),
                            dlarfbTestStruct(6, 8, 10, 0, 0, 0),
                            dlarfbTestStruct(6, 10, 8, 0, 0, 0),
                            dlarfbTestStruct(8, 6, 10, 0, 0, 0),
                            dlarfbTestStruct(8, 10, 6, 0, 0, 0),
                            dlarfbTestStruct(10, 6, 8, 0, 0, 0),
                            dlarfbTestStruct(10, 8, 6, 0, 0, 0),
                            dlarfbTestStruct(6, 6, 6, 12, 15, 30),
                            dlarfbTestStruct(6, 8, 10, 12, 15, 30),
                            dlarfbTestStruct(6, 10, 8, 12, 15, 30),
                            dlarfbTestStruct(8, 6, 10, 12, 15, 30),
                            dlarfbTestStruct(8, 10, 6, 12, 15, 30),
                            dlarfbTestStruct(10, 6, 8, 12, 15, 30),
                            dlarfbTestStruct(10, 8, 6, 12, 15, 30),
                            dlarfbTestStruct(6, 6, 6, 15, 12, 30),
                            dlarfbTestStruct(6, 8, 10, 15, 12, 30),
                            dlarfbTestStruct(6, 10, 8, 15, 12, 30),
                            dlarfbTestStruct(8, 6, 10, 15, 12, 30),
                            dlarfbTestStruct(8, 10, 6, 15, 12, 30),
                            dlarfbTestStruct(10, 6, 8, 15, 12, 30),
                            dlarfbTestStruct(10, 8, 6, 15, 12, 30)
                        ]) {
                        let ma = test.ma
                        let na = test.na
                        var lda = test.lda
                        if (lda == 0) {
                            lda = na
                        }
                        let a = Array<Float64>(ma * lda, repeat: 0.0)
                        for (i in 0..ma) {
                            for (j in 0..lda) {
                                a[i * lda + j] = rnd.nextFloat64()
                            }
                        }
                        let k = min(ma, na)

                        var m: Int64
                        let n: Int64
                        let rowsWork: Int64
                        if (side == Left) {
                            m = test.ma
                            n = test.cdim
                            rowsWork = n
                        } else if (side == Right) {
                            m = test.cdim
                            n = test.ma
                            rowsWork = m
                        } else {
                            throw Exception("not implemented")
                        }

                        let tau = Array<Float64>(k, repeat: 0.0)
                        var work = Array<Float64>(na, repeat: 0.0)
                        lFloat64.dgeqr2(ma, k, a, lda, tau, work)

                        let vMatTmp = extractVMat(ma, na, a, lda, Forward, ColumnWise)
                        let vMat = constructVMat(vMatTmp, store, direct)
                        let v = vMat.data
                        let ldv = vMat.stride

                        var ldt = test.ldt
                        if (ldt == 0) {
                            ldt = k
                        }
                        let tm = Array<Float64>(k * ldt, repeat: 0.0)

                        lFloat64.dlarft(direct, store, ma, k, v, ldv, tau, tm, ldt)

                        var ldc = test.ldc
                        if (ldc == 0) {
                            ldc = n
                        }
                        let c = Array<Float64>(m * ldc, repeat: 0.0)
                        for (i in 0..m) {
                            for (j in 0..ldc) {
                                c[i * ldc + j] = rnd.nextFloat64()
                            }
                        }
                        let cCopy = Array<Float64>(c.size, repeat: 0.0)
                        copy(cCopy, c)

                        let ldwork = k
                        work = Array<Float64>(rowsWork * k, repeat: 0.0)

                        lFloat64.dlarfb(side, trans, direct, store, m, n, k, v, ldv, tm, ldt, c, ldc, work, ldwork)

                        let h = constructH(tau, vMat, store, direct)

                        let cMat = General64(
                            rows: m,
                            cols: n,
                            stride: ldc,
                            data: Array<Float64>(m * ldc, repeat: 0.0)
                        )
                        copy(cMat.data, cCopy)
                        let ans = General64(
                            rows: m,
                            cols: n,
                            stride: ldc,
                            data: Array<Float64>(m * ldc, repeat: 0.0)
                        )
                        copy(ans.data, cMat.data)
                        if (side == Left && trans == NoTrans) {
                            gemm(NoTrans, NoTrans, 1.0, h, cMat, 0.0, ans)
                        } else if (side == Left && trans == Trans) {
                            gemm(Trans, NoTrans, 1.0, h, cMat, 0.0, ans)
                        } else if (side == Right && trans == NoTrans) {
                            gemm(NoTrans, NoTrans, 1.0, cMat, h, 0.0, ans)
                        } else if (side == Right && trans == Trans) {
                            gemm(NoTrans, Trans, 1.0, cMat, h, 0.0, ans)
                        } else {
                            throw Exception("not implemented")
                        }
                        if (!fEqualApprox(ans.data, c, 1e-14)) {
                            throw Exception("Cas ${cas} mismatch. Want ${ans.data}, got ${c}.")
                        }
                        cas++
                    }
                }
            }
        }
    }
    return true
}
