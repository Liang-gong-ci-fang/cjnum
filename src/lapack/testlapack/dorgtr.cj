package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dorgtrTestStruct {
    dorgtrTestStruct(
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dorgtrTest(lFloat64: LFloat64): Bool {
    const tol = 1e-13

    let rnd = Random(1)
    for (uplo in [Upper, Lower]) {
        for (wl in [minimumWork, mediumWork, optimumWork]) {
            for (test in [
                    dorgtrTestStruct(1, 0),
                    dorgtrTestStruct(2, 0),
                    dorgtrTestStruct(3, 0),
                    dorgtrTestStruct(6, 0),
                    dorgtrTestStruct(33, 0),
                    dorgtrTestStruct(100, 0),
                    dorgtrTestStruct(1, 3),
                    dorgtrTestStruct(2, 5),
                    dorgtrTestStruct(3, 7),
                    dorgtrTestStruct(6, 10),
                    dorgtrTestStruct(33, 50),
                    dorgtrTestStruct(100, 120)
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
                let e = Array<Float64>(n - 1, repeat: 0.0)
                let tau = Array<Float64>(n - 1, repeat: 0.0)

                var work = Array<Float64>(1, repeat: 0.0)
                lFloat64.dsytrd(uplo, n, a, lda, d, e, tau, work, -1)
                work = Array<Float64>(Int64(work[0]), repeat: 0.0)

                lFloat64.dsytrd(uplo, n, a, lda, d, e, tau, work, work.size)

                var lwork: Int64 = 0
                if (wl == minimumWork) {
                    lwork = max(1, n - 1)
                } else if (wl == mediumWork) {
                    work = Array<Float64>(1, repeat: 0.0)
                    lFloat64.dorgtr(uplo, n, a, lda, tau, work, -1)
                    lwork = (Int64(work[0]) + n - 1) / 2
                    lwork = max(1, lwork)
                } else if (wl == optimumWork) {
                    work = Array<Float64>(1, repeat: 0.0)
                    lFloat64.dorgtr(uplo, n, a, lda, tau, work, -1)
                    lwork = Int64(work[0])
                }
                work = nanSlice(lwork)

                lFloat64.dorgtr(uplo, n, a, lda, tau, work, work.size)
                let q = General64(
                    rows: n,
                    cols: n,
                    stride: lda,
                    data: a
                )

                let name = "uplo=${uplo},n=${n},lda=${lda},work=${wl}"

                let resid = residualOrthogonal(q, false)
                if (resid > tol * Float64(n)) {
                    throw Exception("Case ${name}: Q is not orthogonal; resid=${resid}, want<=${tol*Float64(n)}")
                }

                let tri = General64(
                    rows: n,
                    cols: n,
                    stride: n,
                    data: Array<Float64>(n * n, repeat: 0.0)
                )
                for (i in 0..n) {
                    tri.data[i * tri.stride + i] = d[i]
                    if (i != n - 1) {
                        tri.data[i * tri.stride + i + 1] = e[i]
                        tri.data[(i + 1) * tri.stride + i] = e[i]
                    }
                }

                let aMat = General64(
                    rows: n,
                    cols: n,
                    stride: n,
                    data: Array<Float64>(n * n, repeat: 0.0)
                )
                if (uplo == Upper) {
                    for (i in 0..n) {
                        for (j in i..n) {
                            let v = aCopy[i * lda + j]
                            aMat.data[i * aMat.stride + j] = v
                            aMat.data[j * aMat.stride + i] = v
                        }
                    }
                } else {
                    for (i in 0..n) {
                        for (j in 0..=i) {
                            let v = aCopy[i * lda + j]
                            aMat.data[i * aMat.stride + j] = v
                            aMat.data[j * aMat.stride + i] = v
                        }
                    }
                }

                let tmp = General64(rows: n, cols: n, stride: n, data: Array<Float64>(n * n, repeat: 0.0))
                gemm(NoTrans, NoTrans, 1.0, aMat, q, 0.0, tmp)
                let ans = General64(rows: n, cols: n, stride: n, data: Array<Float64>(n * n, repeat: 0.0))
                gemm(Trans, NoTrans, 1.0, q, tmp, 0.0, ans)

                if (!fEqualApprox(ans.data, tri.data, tol)) {
                    // "uplo=85,n=1,lda=1,work=0"
                    throw Exception("Case ${name}: Recombination mismatch")
                }
            }
        }
    }
    return true
}
