package cjnum.lapack.testlapack

import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.floats.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64

public func dgeql2Test(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    for (test in [
		[5, 5, 0],
		[5, 3, 0],
		[5, 4, 0]
    ]) {
        let m = test[0]
        let n = test[1]
        let lda = if (test[2] == 0) {
            n
        } else {
            test[2]
        }
        let a = Array<Float64>(m * lda, { _ => rnd.nextGaussianFloat64() })
        let tau = nanSlice(min(m, n))
        let work = nanSlice(n)

        let aCopy = a.clone()
        lFloat64.dgeql2(m, n, a, lda, tau, work)

        let k = min(m, n)
        let q = General64(
            rows: m,
            cols: m,
            stride: m,
            data: Array<Float64>(m * m, repeat: 0.0)
        )
        for (i in 0..m) {
            q.data[i * q.stride + i] = 1.0
        }
        for (i in 0..k) {
            let h = General64(
                rows: m,
                cols: m,
                stride: m,
                data: Array<Float64>(m * m, repeat: 0.0)
            )
            for (j in 0..m) {
                h.data[j * h.stride + j] = 1.0
            }
            let v = Vector(
                inc: 1,
                data: Array<Float64>(m, repeat: 0.0)
            )
            v.data[m - k + i] = 1.0
            for (j in 0..m - k + i) {
                v.data[j] = a[j * lda + n - k + i]
            }
            ger(-tau[i], v, v, h)
            let qTmp = General64(
                rows: q.rows,
                cols: q.cols,
                stride: q.stride,
                data: Array<Float64>(q.data.size, repeat: 0.0)
            )
            copy(qTmp.data, q.data)
            gemm(NoTrans, NoTrans, 1.0, h, qTmp, 0.0, q)
        }

        let resid = residualOrthogonal(q, false)
        if (resid > tol) {
            throw Exception("Q is not orthogonal; resid=${resid}, want<=${tol}")
        }

        let l = General64(
            rows: m,
            cols: n,
            stride: n,
            data: Array<Float64>(m * n, repeat: 0.0)
        )
        if (m >= n) {
            for (i in m - n..m) {
                for (j in 0..=min(i - (m - n), n - 1)) {
                    l.data[i * l.stride + j] = a[i * lda + j]
                }
            }
        } else {
            throw Exception("untested")
        }

        let ans = General64(
            rows: m,
            cols: n,
            stride: lda,
            data: Array<Float64>(a.size, repeat: 0.0)
        )
        gemm(NoTrans, NoTrans, 1.0, q, l, 0.0, ans)
        if (!fEqualApprox(ans.data, aCopy, tol)) {
            throw Exception("Reconstruction mismatch: m = ${m}, n = ${n}")
        }
    }
    return true
}