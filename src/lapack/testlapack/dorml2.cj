package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dorml2TestStruct {
    dorml2TestStruct(
        let common: Int64,
        let adim: Int64,
        let cdim: Int64,
        let lda: Int64,
        let ldc: Int64
    ) {}
}

public func dorml2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (side in [Left, Right]) {
        for (trans in [NoTrans, Trans]) {
            for (test in [
                    dorml2TestStruct(3, 4, 5, 0, 0),
                    dorml2TestStruct(3, 5, 4, 0, 0),
                    dorml2TestStruct(4, 3, 5, 0, 0),
                    dorml2TestStruct(4, 5, 3, 0, 0),
                    dorml2TestStruct(5, 3, 4, 0, 0),
                    dorml2TestStruct(5, 4, 3, 0, 0),
                    dorml2TestStruct(3, 4, 5, 6, 20),
                    dorml2TestStruct(3, 5, 4, 6, 20),
                    dorml2TestStruct(4, 3, 5, 6, 20),
                    dorml2TestStruct(4, 5, 3, 6, 20),
                    dorml2TestStruct(5, 3, 4, 6, 20),
                    dorml2TestStruct(5, 4, 3, 6, 20),
                    dorml2TestStruct(3, 4, 5, 20, 6),
                    dorml2TestStruct(3, 5, 4, 20, 6),
                    dorml2TestStruct(4, 3, 5, 20, 6),
                    dorml2TestStruct(4, 5, 3, 20, 6),
                    dorml2TestStruct(5, 3, 4, 20, 6),
                    dorml2TestStruct(5, 4, 3, 20, 6)
                ]) {
                var ma = 0
                var na = 0
                var mc = 0
                var nc = 0
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
                    lda = na
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

                let q = constructQ("LQ", ma, na, a, lda, tau)

                let cMat = General64(
                    rows: mc,
                    cols: nc,
                    stride: ldc,
                    data: Array<Float64>(c.size, repeat: 0.0)
                )
                copy(cMat.data, c)
                let cMatCopy = General64(
                    rows: cMat.rows,
                    cols: cMat.cols,
                    stride: cMat.stride,
                    data: Array<Float64>(cMat.data.size, repeat: 0.0)
                )
                copy(cMatCopy.data, cMat.data)

                if (side == Left && trans == NoTrans) {
                    gemm(NoTrans, NoTrans, 1.0, q, cMatCopy, 0.0, cMat)
                } else if (side == Left && trans == Trans) {
                    gemm(Trans, NoTrans, 1.0, q, cMatCopy, 0.0, cMat)
                } else if (side == Right && trans == NoTrans) {
                    gemm(NoTrans, NoTrans, 1.0, cMatCopy, q, 0.0, cMat)
                } else if (side == Right && trans == Trans) {
                    gemm(NoTrans, Trans, 1.0, cMatCopy, q, 0.0, cMat)
                } else {
                    throw Exception("bad test")
                }

                if (side == Left) {
                    work = Array<Float64>(nc, repeat: 0.0)
                } else {
                    work = Array<Float64>(mc, repeat: 0.0)
                }
                let aCopy = Array<Float64>(a.size, repeat: 0.0)
                copy(aCopy, a)
                let tauCopy = Array<Float64>(tau.size, repeat: 0.0)
                copy(tauCopy, tau)
                lFloat64.dorml2(side, trans, mc, nc, k, a, lda, tau, c, ldc, work)

                if (!fEqual(a, aCopy)) {
                    throw Exception("a changed in call")
                }
                if (!fEqual(tau, tauCopy)) {
                    throw Exception("tau changed in call")
                }
                if (!fEqualApprox(cMat.data, c, 1e-14)) {
                    let isLeft = side == Left
                    let isTrans = trans == Trans
                    throw Exception("Multiplication mismatch. IsLeft = ${isLeft}. IsTrans = ${isTrans}")
                }
            }
        }
    }
    return true
}
