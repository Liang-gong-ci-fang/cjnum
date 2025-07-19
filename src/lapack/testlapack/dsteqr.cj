package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dsteqrTestStruct {
    dsteqrTestStruct(
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dsteqrTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (compz in [EVOrig, EVTridiag]) {
        for (test in [
                dsteqrTestStruct(1, 0),
                dsteqrTestStruct(4, 0),
                dsteqrTestStruct(8, 0),
                dsteqrTestStruct(10, 0),
                dsteqrTestStruct(2, 10),
                dsteqrTestStruct(8, 10),
                dsteqrTestStruct(10, 20)
            ]) {
            for (_ in 0..100) {
                let n = test.n
                let lda = if (test.lda == 0) {
                    n
                } else {
                    test.lda
                }
                let d = Array<Float64>(n, item: 0.0)
                for (i in 0..d.size) {
                    d[i] = rnd.nextFloat64()
                }
                let e = Array<Float64>(n - 1, item: 0.0)
                for (i in 0..e.size) {
                    e[i] = rnd.nextFloat64()
                }
                let a = Array<Float64>(n * lda, item: 0.0)
                for (i in 0..a.size) {
                    a[i] = rnd.nextFloat64()
                }
                let dCopy = Array<Float64>(d.size, item: 0.0)
                copy(dCopy, d)
                let eCopy = Array<Float64>(e.size, item: 0.0)
                copy(eCopy, e)
                let aCopy = Array<Float64>(a.size, item: 0.0)
                copy(aCopy, a)
                if (compz == EVOrig) {
                    let uplo = Upper
                    let tau = Array<Float64>(n, item: 0.0)
                    var work = Array<Float64>(1, item: 0.0)
                    lFloat64.dsytrd(Upper, n, a, lda, d, e, tau, work, -1)
                    work = Array<Float64>(Int64(work[0]), item: 0.0)
                    lFloat64.dsytrd(uplo, n, a, lda, d, e, tau, work, work.size)
                    lFloat64.dorgtr(uplo, n, a, lda, tau, work, work.size)
                } else {
                    for (i in 0..n) {
                        for (j in 0..n) {
                            a[i * lda + j] = 0.0
                            if (i == j) {
                                a[i * lda + j] = 1.0
                            }
                        }
                    }
                }
                let work = Array<Float64>(2 * n, item: 0.0)
                
                let aDecomp = Array<Float64>(a.size, item: 0.0)
                copy(aDecomp, a)
                let dDecomp = Array<Float64>(d.size, item: 0.0)
                copy(dDecomp, d)
                let eDecomp = Array<Float64>(e.size, item: 0.0)
                copy(eDecomp, e)
                lFloat64.dsteqr(compz, n, d, e, a, lda, work)
                let dAns = Array<Float64>(d.size, item: 0.0)
                copy(dAns, d)
                var truth: General64
                if (compz == EVOrig) {
                    truth = General64(
                        rows: n, 
                        cols: n, 
                        stride: n, 
                        data: Array<Float64>(n * n, item: 0.0)
                    )
                    for (i in 0..n) {
                        for (j in i..n) {
                            let v = aCopy[i * lda + j]
                            truth.data[i * truth.stride + j] = v
                            truth.data[j * truth.stride + i] = v
                        }
                    }
                } else {
                    truth = General64(
                        rows: n, 
                        cols: n, 
                        stride: n, 
                        data: Array<Float64>(n * n, item: 0.0)
                    )
                    for (i in 0..n) {
                        truth.data[i * truth.stride + i] = dCopy[i]
                        if (i != n - 1) {
                            truth.data[(i + 1) * truth.stride + i] = eCopy[i]
                            truth.data[i * truth.stride + i + 1] = eCopy[i]
                        }
                    }
                }
                let V = General64(rows: n, cols: n, stride: lda, data: a)
                if (!eigenDecompCorrect(d, truth, V)) {
                    throw Exception("Eigen reconstruction mismatch. fromFull = ${compz == EVOrig}, n = ${n}")
                }
                for (i in 0..work.size) {
                    work[i] = rnd.nextFloat64()
                }
                lFloat64.dsteqr(EVCompNone, n, dDecomp, eDecomp, aDecomp, lda, work)
                if (!fEqualApprox(d, dAns, 1e-8)) {
                    throw Exception("Eigenvalue mismatch when eigenvectors not computed")
                }
            }
        }
    }
    return true
}

func eigenDecompCorrect(values: Array<Float64>, A: General64, V: General64): Bool {
    let n = A.rows
    for (i in 0..n) {
        let lambda = values[i]
        let vector = Array<Float64>(n, item: 0.0)
        let ans2 = Array<Float64>(n, item: 0.0)
        for (j in 0..vector.size) {
            let v = V.data[j * V.stride + i]
            vector[j] = v
            ans2[j] = lambda * v
        }
        let v = Vector(inc: 1, data: vector)
        let ans1 = Vector(inc: 1, data: Array<Float64>(n, item: 0.0))
        gemv(NoTrans, 1.0, A, v, 0.0, ans1)
        if (!fEqualApprox(ans1.data, ans2, 1e-8)) {
            return false
        }
    }
    return true
}
