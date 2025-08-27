package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.random.*

public func dorgqlTest(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    for (m in [0, 1, 2, 3, 4, 5, 7, 10, 15, 30, 50, 150]) {
        for (extra in [0, 11]) {
            for (wl in [minimumWork, mediumWork, optimumWork]) {
                var m_ = m
                var k = 0
                if (m_ >= 129) {
                    k = 129 + rnd.nextInt64(m_ - 129 + 1)
                } else {
                    k = rnd.nextInt64(m + 1)
                }
                var n = k + rnd.nextInt64(m - k + 1)
                if (m_ == 0 || n == 0) {
                    m_ = 0
                    n = 0
                    k = 0
                }

                let a = nanGeneral(m_, n, n + extra)
                let tau = Array<Float64>(k, repeat: 0.0)
                for (l in 0..k) {
                    let jj = m_ - k + l
                    let v = randomSlice(jj, rnd)
                    (_, tau[l]) = lFloat64.dlarfg(v.size + 1, rnd.nextGaussianFloat64(), v, 1)
                    let j = n - k + l
                    for (i in 0..jj) {
                        a.data[i * a.stride + j] = v[i]
                    }
                }
                let aCopy = cloneGeneral(a)

                let q = eye(m_, m_)
                let qCopy = eye(m_, m_)
                for (l in 0..k) {
                    let h = eye(m_, m_)
                    let jj = m_ - k + l
                    let j = n - k + l
                    let v = Vector(data: Array<Float64>(m_, repeat: 0.0), inc: 1)
                    for (i in 0..jj) {
                        v.data[i] = a.data[i * a.stride + j]
                    }
                    v.data[jj] = 1.0
                    ger(-tau[l], v, v, h)
                    copy(qCopy.data, q.data)
                    gemm(NoTrans, NoTrans, 1.0, h, qCopy, 0.0, q)
                }

                let want = General64(
                    rows: m_,
                    cols: n,
                    stride: q.stride,
                    data: q.data[m_ - n..]
                )

                var lwork = 0
                if (wl == minimumWork) {
                    lwork = max(1, n)
                } else if (wl == mediumWork) {
                    let work = Array<Float64>(1, repeat: 0.0)
                    lFloat64.dorgql(m_, n, k, a.data, a.stride, tau, work, -1)
                    lwork = (Int64(work[0]) + n) / 2
                    lwork = max(1, lwork)
                } else if (wl == optimumWork) {
                    let work = Array<Float64>(1, repeat: 0.0)
                    lFloat64.dorgql(m_, n, k, a.data, a.stride, tau, work, -1)
                    lwork = Int64(work[0])
                }
                let work = Array<Float64>(lwork, repeat: 0.0)

                lFloat64.dorgql(m_, n, k, a.data, a.stride, tau, work, work.size)

                let prefix = "m=${m_},n=${n},k=${k},wl=${wl}"
                if (!generalOutsideAllNaN(a)) {
                    throw Exception("${prefix}: out-of-range write to A")
                }
                if (!equalApproxGeneral(want, a, tol)) {
                    throw Exception("${prefix}: unexpected Q")
                }

                // if (!hasDorg2l) {
                //     continue
                // }
                lFloat64.dorg2l(m_, n, k, aCopy.data, aCopy.stride, tau, work)
                if (!equalApproxGeneral(aCopy, a, tol)) {
                    throw Exception("${prefix}: mismatch between Dorgql and Dorg2l")
                }
            }
        }
    }
    return true
}
