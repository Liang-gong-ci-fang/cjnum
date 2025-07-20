package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.random.*
import std.math.*

class dsytrdTestStruct {
    dsytrdTestStruct(
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dsytrdTest(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    var tc = 0
    for (test in [
            dsytrdTestStruct(1, 0),
            dsytrdTestStruct(2, 0),
            dsytrdTestStruct(3, 0),
            dsytrdTestStruct(4, 0),
            dsytrdTestStruct(10, 0),
            dsytrdTestStruct(50, 0),
            dsytrdTestStruct(100, 0),
            dsytrdTestStruct(150, 0),
            dsytrdTestStruct(300, 0),
            dsytrdTestStruct(1, 3),
            dsytrdTestStruct(2, 3),
            dsytrdTestStruct(3, 7),
            dsytrdTestStruct(4, 9),
            dsytrdTestStruct(10, 20),
            dsytrdTestStruct(50, 70),
            dsytrdTestStruct(100, 120),
            dsytrdTestStruct(150, 170),
            dsytrdTestStruct(300, 320)
        ]) {
        for (uplo in [Upper, Lower]) {
            for (wl in [minimumWork, mediumWork, optimumWork]) {
                let n = test.n
                var lda = test.lda
                if (lda == 0) {
                    lda = n
                }
                let a = randomGeneral(n, n, lda, rnd)
                for (i in 1..n) {
                    for (j in 0..i) {
                        a.data[i * a.stride + j] = a.data[j * a.stride + i]
                    }
                }
                let aCopy = cloneGeneral(a)

                let d = nanSlice(n)
                let e = nanSlice(n - 1)
                let tau = nanSlice(n - 1)

                var lwork: Int64 = 0
                if (wl == minimumWork) {
                    lwork = 1
                } else if (wl == mediumWork) {
                    let work = Array<Float64>(1, item: 0.0)
                    lFloat64.dsytrd(uplo, n, a.data, a.stride, d, e, tau, work, -1)
                    lwork = (Int64(work[0]) + 1) / 2
                    lwork = max(1, lwork)
                } else if (wl == optimumWork) {
                    let work = Array<Float64>(1, item: 0.0)
                    lFloat64.dsytrd(uplo, n, a.data, a.stride, d, e, tau, work, -1)
                    lwork = Int64(work[0])
                }
                var work = Array<Float64>(lwork, item: 0.0)

                lFloat64.dsytrd(uplo, n, a.data, a.stride, d, e, tau, work, lwork)

                let prefix = "Case #${tc}: uplo=${uplo},n=${n},lda=${lda},work=${wl}"

                if (!generalOutsideAllNaN(a)) {
                    throw Exception("${prefix}: out-of-range write to A")
                }

                let q = cloneGeneral(a)
                if (uplo == Upper) {
                    for (j in 0..n - 1) {
                        for (i in 0..j) {
                            q.data[i * q.stride + j] = q.data[i * q.stride + j + 1]
                        }
                        q.data[(n - 1) * q.stride + j] = 0.0
                    }
                    for (i in 0..n - 1) {
                        q.data[i * q.stride + n - 1] = 0.0
                    }
                    q.data[(n - 1) * q.stride + n - 1] = 1.0
                    if (n > 1) {
                        work = Array<Float64>(n - 1, item: 0.0)
                        lFloat64.dorgql(n - 1, n - 1, n - 1, q.data, q.stride, tau, work, work.size)
                    }
                } else {
                    for (j in n - 1..0 : -1) {
                        q.data[j] = 0.0
                        for (i in j + 1..n) {
                            q.data[i * q.stride + j] = q.data[i * q.stride + j - 1]
                        }
                    }
                    q.data[0] = 1.0
                    for (i in 1..n) {
                        q.data[i * q.stride] = 0.0
                    }
                    if (n > 1) {
                        work = Array<Float64>(n - 1, item: 0.0)
                        lFloat64.dorgqr(n - 1, n - 1, n - 1, q.data[q.stride + 1..], q.stride, tau, work, work.size)
                    }
                }
                var resid = residualOrthogonal(q, false)
                if (resid > tol * Float64(n)) {
                    throw Exception("${prefix}: Q is not orthogonal; resid=${resid}, want<=${tol*Float64(n)}")
                }

                let tMat = zeros(n, n, n)
                for (i in 0..n) {
                    tMat.data[i * tMat.stride + i] = d[i]
                }
                if (uplo == Upper) {
                    for (j in 1..n) {
                        tMat.data[(j - 1) * tMat.stride + j] = e[j - 1]
                        tMat.data[j * tMat.stride + j - 1] = e[j - 1]
                    }
                } else {
                    for (j in 0..n - 1) {
                        tMat.data[(j + 1) * tMat.stride + j] = e[j]
                        tMat.data[j * tMat.stride + j + 1] = e[j]
                    }
                }
                let qa = zeros(n, n, n)
                gemm(Trans, NoTrans, 1.0, q, aCopy, 0.0, qa)
                gemm(NoTrans, NoTrans, 1.0, qa, q, -1.0, tMat)
                resid = dlange(MaxColumnSum, n, n, tMat.data, tMat.stride)
                if (resid > tol * Float64(n)) {
                    throw Exception("${prefix}: |Qᵀ*A*Q - T|=${resid}, want<=${tol*Float64(n)}")
                }
            }
        }
        tc++
    }
    return true
}
