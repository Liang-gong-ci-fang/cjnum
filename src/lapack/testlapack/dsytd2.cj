package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dsytd2TestStruct {
    dsytd2TestStruct(
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dsytd2Test(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    for (uplo in [Upper, Lower]) {
        for (test in [
                dsytd2TestStruct(3, 0),
                dsytd2TestStruct(4, 0),
                dsytd2TestStruct(5, 0),
                dsytd2TestStruct(3, 10),
                dsytd2TestStruct(4, 10),
                dsytd2TestStruct(5, 10)
            ]) {
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

            let d = Array<Float64>(n, repeat: 0.0)
            for (i in 0..d.size) {
                d[i] = Float64.NaN
            }
            let e = Array<Float64>(n - 1, repeat: 0.0)
            for (i in 0..e.size) {
                e[i] = Float64.NaN
            }
            let tau = Array<Float64>(n - 1, repeat: 0.0)
            for (i in 0..tau.size) {
                tau[i] = Float64.NaN
            }

            lFloat64.dsytd2(uplo, n, a, lda, d, e, tau)

            let qMat = General64(
                rows: n,
                cols: n,
                stride: n,
                data: Array<Float64>(n * n, repeat: 0.0)
            )
            let qCopy = General64(
                rows: n,
                cols: n,
                stride: n,
                data: Array<Float64>(qMat.data.size, repeat: 0.0)
            )
            for (i in 0..n) {
                qMat.data[i * qMat.stride + i] = 1.0
            }
            for (i in 0..n - 1) {
                let hMat = General64(
                    rows: n,
                    cols: n,
                    stride: n,
                    data: Array<Float64>(n * n, repeat: 0.0)
                )
                for (j in 0..n) {
                    hMat.data[j * hMat.stride + j] = 1.0
                }
                var vi: Vector
                if (uplo == Upper) {
                    vi = Vector(
                        inc: 1,
                        data: Array<Float64>(n, repeat: 0.0)
                    )
                    for (j in 0..i) {
                        vi.data[j] = a[j * lda + i + 1]
                    }
                    vi.data[i] = 1.0
                } else {
                    vi = Vector(
                        inc: 1,
                        data: Array<Float64>(n, repeat: 0.0)
                    )
                    vi.data[i + 1] = 1.0
                    for (j in i + 2..n) {
                        vi.data[j] = a[j * lda + i]
                    }
                }
                ger(-tau[i], vi, vi, hMat)
                copy(qCopy.data, qMat.data)

                if (uplo == Upper) {
                    gemm(NoTrans, NoTrans, 1.0, hMat, qCopy, 0.0, qMat)
                } else {
                    gemm(NoTrans, NoTrans, 1.0, qCopy, hMat, 0.0, qMat)
                }
            }

            let resid = residualOrthogonal(qMat, false)
            if (resid > tol) {
                throw Exception("Q is not orthogonal; resid=${resid}, want<=${tol}")
            }

            let aMat = General64(
                rows: n,
                cols: n,
                stride: n,
                data: Array<Float64>(a.size, repeat: 0.0),
            )

            for (i in 0..n) {
                for (j in i..n) {
                    var v = aCopy[i * lda + j]
                    if (uplo == Lower) {
                        v = aCopy[j * lda + i]
                    }
                    aMat.data[i * aMat.stride + j] = v
                    aMat.data[j * aMat.stride + i] = v
                }
            }

            let tmp = General64(
                rows: n,
                cols: n,
                stride: n,
                data: Array<Float64>(n * n, repeat: 0.0)
            )

            let ans = General64(
                rows: n,
                cols: n,
                stride: n,
                data: Array<Float64>(n * n, repeat: 0.0)
            )

            gemm(Trans, NoTrans, 1.0, qMat, aMat, 0.0, tmp)
            gemm(NoTrans, NoTrans, 1.0, tmp, qMat, 0.0, ans)

            let tMat = General64(
                rows: n,
                cols: n,
                stride: n,
                data: Array<Float64>(n * n, repeat: 0.0),
            )
            for (i in 0..n - 1) {
                tMat.data[i * tMat.stride + i] = d[i]
                tMat.data[i * tMat.stride + i + 1] = e[i]
                tMat.data[(i + 1) * tMat.stride + i] = e[i]
            }
            tMat.data[(n - 1) * tMat.stride + n - 1] = d[n - 1]

            var same = true
            for (i in 0..n) {
                for (j in 0..n) {
                    if (abs(ans.data[i * ans.stride + j] - tMat.data[i * tMat.stride + j]) > tol) {
                        same = false
                    }
                }
            }
            if (!same) {
                throw Exception("Matrix answer mismatch")
            }
        }
    }
    return true
}
