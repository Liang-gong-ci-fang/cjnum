package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dormbrTestStruct {
    dormbrTestStruct(
        let m: Int64,
        let n: Int64,
        let k: Int64,
        let lda: Int64,
        let ldc: Int64
    ) {}
}

public func dormbrTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    let nFloat64 = nFloat64Implementation()
    for (vect in [ApplyQ, ApplyP]) {
        for (side in [Left, Right]) {
            for (trans in [NoTrans, Trans]) {
                for (wl in [minimumWork, mediumWork, optimumWork]) {
                    for (test in [
                            dormbrTestStruct(3, 4, 5, 0, 0),
                            dormbrTestStruct(3, 5, 4, 0, 0),
                            dormbrTestStruct(4, 3, 5, 0, 0),
                            dormbrTestStruct(4, 5, 3, 0, 0),
                            dormbrTestStruct(5, 3, 4, 0, 0),
                            dormbrTestStruct(5, 4, 3, 0, 0),
                            dormbrTestStruct(3, 4, 5, 10, 12),
                            dormbrTestStruct(3, 5, 4, 10, 12),
                            dormbrTestStruct(4, 3, 5, 10, 12),
                            dormbrTestStruct(4, 5, 3, 10, 12),
                            dormbrTestStruct(5, 3, 4, 10, 12),
                            dormbrTestStruct(5, 4, 3, 10, 12),
                            dormbrTestStruct(150, 140, 130, 0, 0)
                        ]) {
                        let m = test.m
                        let n = test.n
                        let k = test.k
                        var ldc = test.ldc
                        if (ldc == 0) {
                            ldc = n
                        }
                        var nq = n
                        var nw = m
                        if (side == Left) {
                            nq = m
                            nw = n
                        }

                        var ma: Int64
                        var na: Int64
                        var a: Array<Float64>
                        if (vect == ApplyQ) {
                            ma = nq
                            na = k
                        } else {
                            ma = k
                            na = nq
                        }
                        var lda = test.lda
                        if (lda == 0) {
                            lda = na
                        }
                        a = Array<Float64>(ma * lda, {_ => rnd.nextGaussianFloat64()})

                        let nTau = min(nq, k)
                        let tauP = Array<Float64>(nTau, repeat: 0.0)
                        let tauQ = Array<Float64>(nTau, repeat: 0.0)
                        let d = Array<Float64>(nTau, repeat: 0.0)
                        let e = Array<Float64>(nTau, repeat: 0.0)

                        var work = Array<Float64>(1, repeat: 0.0)
                        lFloat64.dgebrd(ma, na, a, lda, d, e, tauQ, tauP, work, -1)
                        work = Array<Float64>(Int64(work[0]), repeat: 0.0)
                        lFloat64.dgebrd(ma, na, a, lda, d, e, tauQ, tauP, work, work.size)

                        let c = Array<Float64>(m * ldc, {_ => rnd.nextGaussianFloat64()})

                        let cCopy = Array<Float64>(c.size, repeat: 0.0)
                        copy(cCopy, c)

                        var lwork: Int64 = 0
                        if (wl == minimumWork) {
                            lwork = nw
                        } else if (wl == optimumWork) {
                            lFloat64.dormbr(vect, side, trans, m, n, k, a, lda, tauQ, c, ldc, work, -1)
                            lwork = Int64(work[0])
                        } else if (wl == mediumWork) {
                            work = Array<Float64>(1, repeat: 0.0)
                            lFloat64.dormbr(vect, side, trans, m, n, k, a, lda, tauQ, c, ldc, work, -1)
                            lwork = (Int64(work[0]) + nw) / 2
                        }
                        lwork = max(1, lwork)
                        work = Array<Float64>(lwork, repeat: 0.0)

                        if (vect == ApplyQ) {
                            lFloat64.dormbr(vect, side, trans, m, n, k, a, lda, tauQ, c, ldc, work, lwork)
                        } else {
                            lFloat64.dormbr(vect, side, trans, m, n, k, a, lda, tauP, c, ldc, work, lwork)
                        }

                        let cOrig = General64(
                            rows: m,
                            cols: n,
                            stride: ldc,
                            data: Array<Float64>(cCopy.size, repeat: 0.0)
                        )
                        copy(cOrig.data, cCopy)
                        let cAns = General64(
                            rows: m,
                            cols: n,
                            stride: ldc,
                            data: Array<Float64>(cCopy.size, repeat: 0.0)
                        )
                        copy(cAns.data, cCopy)
                        let nb = min(ma, na)
                        var mulMat: General64
                        if (vect == ApplyQ) {
                            mulMat = constructQPBidiagonal(ApplyQ, ma, na, nb, a, lda, tauQ)
                        } else {
                            mulMat = constructQPBidiagonal(ApplyP, ma, na, nb, a, lda, tauP)
                        }

                        let mulTrans = trans

                        if (side == Left) {
                            nFloat64.dgemm(mulTrans, NoTrans, m, n, m, 1.0, mulMat.data, mulMat.stride, cOrig.data,
                                cOrig.stride, 0.0, cAns.data, cAns.stride)
                        } else {
                            nFloat64.dgemm(NoTrans, mulTrans, m, n, n, 1.0, cOrig.data, cOrig.stride, mulMat.data,
                                mulMat.stride, 0.0, cAns.data, cAns.stride)
                        }

                        if (!fEqualApprox(cAns.data, c, 1e-13)) {
                            let isApplyQ = vect == ApplyQ
                            let isLeft = side == Left
                            let isTrans = trans == Trans

                            throw Exception(
                                "C mismatch. isApplyQ: ${isApplyQ}, isLeft: ${isLeft}, isTrans: ${isTrans}, m = ${m}, n = ${n}, k = ${k}, lda = ${lda}, ldc = ${ldc}"
                            )
                        }
                    }
                }
            }
        }
    }
    return true
}
