package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dlahr2TestStruct {
    dlahr2TestStruct(
        let n: Int64,
        let k: Int64,
        let nb: Int64
    ) {}
}

public func dlahr2Test(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    for (test in [
            dlahr2TestStruct(3, 0, 3),
            dlahr2TestStruct(3, 1, 2),
            dlahr2TestStruct(3, 1, 1),
            dlahr2TestStruct(5, 0, 5),
            dlahr2TestStruct(5, 1, 4),
            dlahr2TestStruct(5, 1, 3),
            dlahr2TestStruct(5, 1, 2),
            dlahr2TestStruct(5, 1, 1),
            dlahr2TestStruct(5, 2, 3),
            dlahr2TestStruct(5, 2, 2),
            dlahr2TestStruct(5, 2, 1),
            dlahr2TestStruct(5, 3, 2),
            dlahr2TestStruct(5, 3, 1),
            dlahr2TestStruct(7, 3, 4),
            dlahr2TestStruct(7, 3, 3),
            dlahr2TestStruct(7, 3, 2),
            dlahr2TestStruct(7, 3, 1),
            dlahr2TestStruct(10, 0, 10),
            dlahr2TestStruct(10, 1, 9),
            dlahr2TestStruct(10, 1, 5),
            dlahr2TestStruct(10, 1, 1),
            dlahr2TestStruct(10, 5, 5),
            dlahr2TestStruct(10, 5, 3),
            dlahr2TestStruct(10, 5, 1)
        ]) {
        for (_ in 0..100) {
            for (extraStride in [0, 1, 10]) {
                let n = test.n
                let k = test.k
                let nb = test.nb

                let a = randomGeneral(n, n - k + 1, n - k + 1 + extraStride, rnd)
                let aCopy = General64(rows: a.rows, cols: a.cols, stride: a.stride,
                    data: Array<Float64>(a.data.size, item: 0.0))
                copy(aCopy.data, a.data)
                let tmat = nanTriangular(Upper, nb, nb + extraStride)
                let y = nanGeneral(n, nb, nb + extraStride)
                let tau = nanSlice(nb)

                lFloat64.dlahr2(n, k, nb, a.data, a.stride, tau, tmat.data, tmat.stride, y.data, y.stride)

                let prefix = "n=${n},k=${k},nb=${nb},ldex=${extraStride}"

                if (!generalOutsideAllNaN(a)) {
                    throw Exception("${prefix}: out-of-range write to A\n${a.data}")
                }
                if (!triangularOutsideAllNaN(tmat)) {
                    throw Exception("${prefix}: out-of-range write to T\n${tmat.data}")
                }
                if (!generalOutsideAllNaN(y)) {
                    throw Exception("${prefix}: out-of-range write to Y\n${y.data}")
                }

                for (i in 0..n) {
                    for (j in 0..n - k + 1) {
                        if (i >= k && j < nb) {
                            continue
                        }
                        if (a.data[i * a.stride + j] != aCopy.data[i * aCopy.stride + j]) {
                            throw Exception("${prefix}: unexpected write to A[${i},${j}]")
                        }
                    }
                }

                for (i in 0..tau.size) {
                    let v = tau[i]
                    if (v.isNaN()) {
                        throw Exception("${prefix}: tau[${i}] not assigned")
                    }
                }

                let v = General64(
                    rows: n - k + 1,
                    cols: nb,
                    stride: nb,
                    data: Array<Float64>((n - k + 1) * nb, item: 0.0)
                )
                for (j in 0..v.cols) {
                    v.data[(j + 1) * v.stride + j] = 1.0
                    for (i in j + 2..v.rows) {
                        v.data[i * v.stride + j] = a.data[(i + k - 1) * a.stride + j]
                    }
                }

                let vt = General64(rows: v.rows, cols: v.cols, stride: v.stride,
                    data: Array<Float64>(v.data.size, item: 0.0))
                copy(vt.data, v.data)
                trmm(Right, NoTrans, 1.0, tmat, vt)
                let ywant = General64(
                    rows: n,
                    cols: nb,
                    stride: nb,
                    data: Array<Float64>(n * nb, item: 0.0)
                )
                gemm(NoTrans, NoTrans, 1.0, aCopy, vt, 0.0, ywant)

                for (i in 0..n) {
                    for (j in 0..nb) {
                        let diff = abs(ywant.data[i * ywant.stride + j] - y.data[i * y.stride + j])
                        if (diff > tol) {
                            throw Exception("${prefix}: unexpected Y[${i},${j}], diff=${diff}")
                        }
                    }
                }

                let q = constructQ("QR", n - k, nb, a.data[k * a.stride..], a.stride, tau)
                var resid = residualOrthogonal(q, false)
                if (resid > tol * Float64(n)) {
                    throw Exception("Case ${prefix}: Q is not orthogonal; resid=${resid}, want<=${tol*Float64(n)}")
                }
                let qwant = General64(
                    rows: n - k + 1,
                    cols: n - k + 1,
                    stride: n - k + 1,
                    data: Array<Float64>((n - k + 1) * (n - k + 1), item: 0.0)
                )
                for (i in 0..qwant.rows) {
                    qwant.data[i * qwant.stride + i] = 1.0
                }
                gemm(NoTrans, Trans, -1.0, vt, v, 1.0, qwant)
                resid = residualOrthogonal(qwant, false)
                if (resid > tol * Float64(n)) {
                    throw Exception(
                        "Case ${prefix}: Q = I - V*T*Vᵀ is not orthogonal; resid=${resid}, want<=${tol*Float64(n)}")
                }

                for (i in 0..n - k) {
                    for (j in 0..n - k) {
                        let diff = abs(q.data[i * q.stride + j] - qwant.data[(i + 1) * qwant.stride + j + 1])
                        if (diff > tol) {
                            throw Exception("${prefix}: unexpected Q[${i},${j}], diff=${diff}")
                        }
                    }
                }
            }
        }
    }
    return true
}
