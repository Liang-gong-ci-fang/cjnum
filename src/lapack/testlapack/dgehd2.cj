package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dgehd2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [1, 2, 3, 4, 5, 7, 10, 30]) {
        for (extra in [0, 1, 13]) {
            for (_ in 0..100) {
                if (!testDgehd2(lFloat64, n, extra, rnd)) {
                    return false
                }
            }
        }
    }
    return true
}

func testDgehd2(lFloat64: LFloat64, n: Int64, extra: Int64, rnd: Random): Bool {
    let tol = 1e-14
    var ilo = rnd.nextInt64(n)
    var ihi = rnd.nextInt64(n)
    if (ilo > ihi) {
        (ilo, ihi) = (ihi, ilo)
    }
    let tau = nanSlice(n - 1)
    let work = nanSlice(n)
    let a = randomGeneral(n, n, n + extra, rnd)
    for (i in 1..=ihi) {
        for (j in 0..min(ilo, i)) {
            a.data[i * a.stride + j] = Float64.NaN
        }
    }
    for (i in ihi + 1..n) {
        for (j in 0..i) {
            a.data[i * a.stride + j] = Float64.NaN
        }
    }
    let aCopy = General64(rows: a.rows, cols: a.cols, stride: a.stride, data: Array<Float64>(a.data.size, repeat: 0.0))
    copy(aCopy.data, a.data)
    lFloat64.dgehd2(n, ilo, ihi, a.data, a.stride, tau, work)
    let prefix = "Case n=${n}, ilo=${ilo}, ihi=${ihi}, extra=${extra}"
    if (!generalOutsideAllNaN(a)) {
        throw Exception("${prefix}: out-of-range write to A\n${a.data}")
    }
    for (i in ilo..=ihi) {
        for (j in 0..min(ilo, i)) {
            if (!a.data[i * a.stride + j].isNaN()) {
                throw Exception("${prefix}: expected NaN at A[${i},${j}]")
            }
        }
    }
    for (i in ihi + 1..n) {
        for (j in 0..i) {
            if (!a.data[i * a.stride + j].isNaN()) {
                throw Exception("${prefix}: expected NaN at A[${i},${j}]")
            }
        }
    }
    for (i in 0..=ilo) {
        for (j in i..ilo + 1) {
            if (a.data[i * a.stride + j] != aCopy.data[i * aCopy.stride + j]) {
                throw Exception("${prefix}: unexpected modification at A[${i},${j}]")
            }
        }
        for (j in ihi + 1..n) {
            if (a.data[i * a.stride + j] != aCopy.data[i * aCopy.stride + j]) {
                throw Exception("${prefix}: unexpected modification at A[${i},${j}]")
            }
        }
    }
    for (i in ihi + 1..n) {
        for (j in i..n) {
            if (a.data[i * a.stride + j] != aCopy.data[i * aCopy.stride + j]) {
                throw Exception("${prefix}: unexpected modification at A[${i},${j}]")
            }
        }
    }
    for (i in 0..tau.size) {
        let v = tau[i]
        if (i < ilo || i >= ihi) {
            if (!v.isNaN()) {
                throw Exception("${prefix}: expected NaN at tau[${i}]")
            }
        } else {
            if (v.isNaN()) {
                throw Exception("${prefix}: unexpected NaN at tau[${i}]")
            }
        }
    }

    let q = General64(
        rows: n,
        cols: n,
        stride: n,
        data: Array<Float64>(n * n, repeat: 0.0)
    )
    for (i in 0..q.rows) {
        q.data[i * q.stride + i] = 1.0
    }
    let qCopy = General64(rows: q.rows, cols: q.cols, stride: q.stride, data: Array<Float64>(a.data.size, repeat: 0.0))
    for (j in ilo..ihi) {
        let h = General64(
            rows: n,
            cols: n,
            stride: n,
            data: Array<Float64>(n * n, repeat: 0.0)
        )
        for (i in 0..h.rows) {
            h.data[i * h.stride + i] = 1.0
        }
        let v = Vector(
            inc: 1,
            data: Array<Float64>(n, repeat: 0.0)
        )
        v.data[j + 1] = 1.0
        for (i in j + 2..=ihi) {
            v.data[i] = a.data[i * a.stride + j]
        }
        ger(-tau[j], v, v, h)
        copy(qCopy.data, q.data)
        gemm(NoTrans, NoTrans, 1.0, qCopy, h, 0.0, q)
    }

    let resid = residualOrthogonal(q, false)
    if (resid > tol) {
        throw Exception("${prefix}: Q is not orthogonal; resid=${resid}, want<=${tol}")
    }

    for (i in 1..=ihi) {
        for (j in 0..min(ilo, i)) {
            aCopy.data[i * aCopy.stride + j] = 0.0
        }
    }
    for (i in ihi + 1..n) {
        for (j in 0..i) {
            aCopy.data[i * aCopy.stride + j] = 0.0
        }
    }

    let aq = General64(
        rows: n,
        cols: n,
        stride: n,
        data: Array<Float64>(n * n, repeat: 0.0)
    )
    gemm(NoTrans, NoTrans, 1.0, aCopy, q, 0.0, aq)

    let qaq = General64(
        rows: n,
        cols: n,
        stride: n,
        data: Array<Float64>(n * n, repeat: 0.0)
    )
    gemm(Trans, NoTrans, 1.0, q, aq, 0.0, qaq)

    for (i in ilo..=ihi) {
        for (j in ilo..=ihi) {
            let qaqij = qaq.data[i * qaq.stride + j]
            if (j < i - 1) {
                if (abs(qaqij) > tol) {
                    throw Exception("${prefix}: Qᵀ*A*Q is not upper Hessenberg, [${i},${j}]=${qaqij}")
                }
                continue
            }
            let diff = qaqij - a.data[i * a.stride + j]
            if (abs(diff) > tol) {
                throw Exception("${prefix}: Qᵀ*AOrig*Q and A are not equal, diff at [${i},${j}]=${diff}")
            }
        }
    }

    return true
}
