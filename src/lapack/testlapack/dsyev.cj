package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dsyevTestStruct {
    dsyevTestStruct(
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dsyevTest(lFloat64: LFloat64) {
    let rnd = Random(1)
    for (uplo in [Lower, Upper]) {
        for (test in [
                dsyevTestStruct(1, 0),
                dsyevTestStruct(2, 0),
                dsyevTestStruct(5, 0),
                dsyevTestStruct(10, 0),
                dsyevTestStruct(100, 0),
                dsyevTestStruct(1, 5),
                dsyevTestStruct(2, 5),
                dsyevTestStruct(5, 10),
                dsyevTestStruct(10, 20),
                dsyevTestStruct(100, 110)
            ]) {
            for (_ in 0..10) {
                let n = test.n
                var lda = test.lda
                if (lda == 0) {
                    lda = n
                }
                let a = Array<Float64>(n * lda, repeat: 0.0)
                for (i in 0..a.size) {
                    a[i] = rnd.nextGaussianFloat64()
                }
                let aCopy = Array<Float64>(a.size, repeat: 0.0)
                copy(aCopy, a)
                var w = Array<Float64>(n, repeat: 0.0)
                for (i in 0..w.size) {
                    w[i] = rnd.nextGaussianFloat64()
                }

                var work = Array<Float64>(1, repeat: 0.0)
                lFloat64.dsyev(EVCompute, uplo, n, a, lda, w, work, -1)
                work = Array<Float64>(Int64(work[0]), repeat: 0.0)
                lFloat64.dsyev(EVCompute, uplo, n, a, lda, w, work, work.size)

                let orig = General64(
                    rows: n,
                    cols: n,
                    stride: n,
                    data: Array<Float64>(n * n, repeat: 0.0)
                )
                if (uplo == Upper) {
                    for (i in 0..n) {
                        for (j in i..n) {
                            let v = aCopy[i * lda + j]
                            orig.data[i * orig.stride + j] = v
                            orig.data[j * orig.stride + i] = v
                        }
                    }
                } else {
                    for (i in 0..n) {
                        for (j in 0..=i) {
                            let v = aCopy[i * lda + j]
                            orig.data[i * orig.stride + j] = v
                            orig.data[j * orig.stride + i] = v
                        }
                    }
                }

                let V = General64(
                    rows: n,
                    cols: n,
                    stride: lda,
                    data: a,
                )

                if (!eigenDecompCorrect(w, orig, V)) {
                    throw Exception("Decomposition mismatch")
                }

                let wAns = Array<Float64>(w.size, repeat: 0.0)
                copy(wAns, w)
                copy(a, aCopy)
                for (i in 0..w.size) {
                    w[i] = rnd.nextFloat64()
                }
                for (i in 0..work.size) {
                    work[i] = rnd.nextFloat64()
                }
                lFloat64.dsyev(EVNone, uplo, n, a, lda, w, work, work.size)
                if (!fEqualApprox(w, wAns, 1e-8)) {
                    throw Exception("Eigenvalue mismatch when vectors not computed")
                }
            }
        }
    }
    return true
}
