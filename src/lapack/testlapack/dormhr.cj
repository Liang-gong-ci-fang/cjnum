package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dormhrTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)

    for (side in [Left, Right]) {
        for (trans in [NoTrans, Trans]) {
            for (m in [1, 2, 3, 4, 5, 8, 9, 10, 23]) {
                for (n in [1, 2, 3, 4, 5, 8, 9, 10, 23]) {
                    for (extra in [0, 1, 13]) {
                        for (_ in 0..10) {
                            var nq = m
                            if (side == Right) {
                                nq = n
                            }
                            var ilo = rnd.nextInt64(nq)
                            var ihi = rnd.nextInt64(nq)
                            if (ilo > ihi) {
                                (ilo, ihi) = (ihi, ilo)
                            }
                            if (!testDormhr(lFloat64, side, trans, m, n, ilo, ihi, extra, true, rnd)) {
                                return false
                            }
                            if (!testDormhr(lFloat64, side, trans, m, n, ilo, ihi, extra, false, rnd)) {
                                return false
                            }
                        }
                    }
                }
            }
        }
    }
    for (side in [Left, Right]) {
        for (trans in [NoTrans, Trans]) {
            if (!testDormhr(lFloat64, side, trans, 0, 0, 0, -1, 0, true, rnd)) {
                return false
            }
            if (!testDormhr(lFloat64, side, trans, 0, 0, 0, -1, 0, false, rnd)) {
                return false
            }
        }
    }
    return true
}

func testDormhr(lFloat64: LFloat64, side: Side, trans: Transpose, m: Int64, n: Int64, ilo: Int64, ihi: Int64,
    extra: Int64, optwork: Bool, rnd: Random): Bool {
    const tol = 1e-14

    var nq: Int64 = 0
    var nw: Int64 = 0
    if (side == Left) {
        nq = m
        nw = n
    } else if (side == Right) {
        nq = n
        nw = m
    }

    let a = randomGeneral(nq, nq, nq + extra, rnd)
    var tau: Array<Float64> = Array<Float64>(0, repeat: 0.0)
    if (nq > 1) {
        tau = nanSlice(nq - 1)
    }
    var work = nanSlice(max(1, nq))
    lFloat64.dgehrd(nq, ilo, ihi, a.data, a.stride, tau, work, work.size)

    let q = eye(nq, nq)
    let qCopy = eye(nq, nq)
    for (j in ilo..ihi) {
        let h = eye(nq, nq)
        let v = Vector(
            inc: 1,
            data: Array<Float64>(nq, repeat: 0.0)
        )
        v.data[j + 1] = 1.0
        for (i in j + 2..ihi + 1) {
            v.data[i] = a.data[i * a.stride + j]
        }
        ger(-tau[j], v, v, h)
        copy(qCopy.data, q.data)
        gemm(NoTrans, NoTrans, 1.0, qCopy, h, 0.0, q)
    }

    let c = randomGeneral(m, n, n + extra, rnd)

    let qc = randomGeneral(m, n, n + extra, rnd)
    if (side == Left) {
        gemm(trans, NoTrans, 1.0, q, c, 0.0, qc)
    } else {
        gemm(NoTrans, trans, 1.0, c, q, 0.0, qc)
    }

    if (optwork) {
        work = nanSlice(1)
        lFloat64.dormhr(side, trans, m, n, ilo, ihi, [], a.stride, [], [], c.stride, work, -1)
        work = nanSlice(Int64(work[0]))
    } else {
        work = nanSlice(max(1, nw))
    }
    lFloat64.dormhr(side, trans, m, n, ilo, ihi, a.data, a.stride, tau, c.data, c.stride, work, work.size)

    let prefix = "n=${m},lda=${n},ilo=${ilo},ihi=${ihi},extra=${extra},optwork=${optwork}"
    if (!generalOutsideAllNaN(c)) {
        throw Exception("${prefix}: out-of-range write to C\n${c.data}")
    }
    for (i in 0..m) {
        for (j in 0..n) {
            let cij = c.data[i * c.stride + j]
            let qcij = qc.data[i * qc.stride + j]
            if (abs(cij - qcij) > tol) {
                throw Exception(
                    "${prefix}: unexpected value of the QC product at [${i},${j}]: want ${qcij}, got ${cij}")
            }
        }
    }
    return true
}
