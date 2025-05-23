package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import std.collection.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64
import cjnum.lapack.*
import cjnum.complex.*

const dlamchE = 0x1p-53
const dlamchB = 2.0
const dlamchP = dlamchB * dlamchE
const dlamchS = 0x1p-1022
const safmin = dlamchS
const safmax = 1.0 / safmin
const ulp = dlamchP
const smlnum = safmin / ulp
const bignum = safmax * ulp

type worklen = Int64

const minimumWork: worklen = 0
const mediumWork: worklen = 1
const optimumWork: worklen = 2

func worklenTOString(work: worklen): String {
    return if (work == minimumWork) {
        "minimum"
    } else if (work == mediumWork) {
        "medium"
    } else if (work == optimumWork) {
        "optimum"
    } else {
        ""
    }
}

func diagToString(diag: Diag): String {
    if (diag == NonUnitDiag) {
        return "NonUnitDiag"
    } else if (diag == UnitDiag) {
        return "UnitDiag"
    } else {
        throw Exception("invalid diag")
    }
}

func normToString(norm: MatrixNorm): String {
    return if (norm == MaxAbs) {
        "MaxAbs"
    } else if (norm == MaxRowSum) {
        "MaxRowSum"
    } else if (norm == MaxColumnSum) {
        "MaxColSum"
    } else if (norm == Frobenius) {
        "Frobenius"
    } else {
        throw Exception("invalid norm")
    }
}

func uploToString(uplo: Uplo): String {
    return if (uplo == Lower) {
        "Lower"
    } else if (uplo == Upper) {
        "Upper"
    } else {
        throw Exception("invalid uplo")
    }
}

func sideToString(side: Side): String {
    return if (side == Left) {
        "Left"
    } else if (side == Right) {
        "Right"
    } else {
        throw Exception("invalid side")
    }
}

func nanSlice(n: Int64): Array<Float64> {
    return Array<Float64>(
        n,
        {
            _ => Float64.NaN
        }
    )
}

func randomSlice(n: Int64, rnd: Random): Array<Float64> {
    return Array<Float64>(
        n,
        {
            _ => rnd.nextFloat64()
        }
    )
}

func nanGeneral(r: Int64, c: Int64, stride: Int64): General64 {
    if (r < 0 || c < 0) {
        throw Exception("bad matrix size")
    }
    if (r == 0 || c == 0) {
        return General64(stride: max(1, stride))
    }
    if (stride < 0) {
        throw Exception("bad stride")
    }
    return General64(
        rows: r,
        cols: c,
        stride: stride,
        data: nanSlice((r - 1) * stride + c)
    )
}

func randomGeneral(r: Int64, c: Int64, stride: Int64, rnd: Random): General64 {
    let ans = nanGeneral(r, c, stride)
    for (i in 0..r) {
        for (j in 0..c) {
            ans.data[i * ans.stride + j] = rnd.nextGaussianFloat64()
        }
    }
    return ans
}

func randomHessenberg(n: Int64, stride: Int64, rnd: Random): General64 {
    let ans = nanGeneral(n, n, stride)
    for (i in 0..n) {
        for (j in 0..i - 1) {
            ans.data[i * ans.stride + j] = 0.0
        }
        for (j in max(0, n - 1)..n) {
            ans.data[i * ans.stride + j] = rnd.nextGaussianFloat64()
        }
    }
    return ans
}

func randomSchurCanonical(n: Int64, stride: Int64, bad: Bool, rnd: Random): (General64, Array<Float64>, Array<Float64>) {
    let t = randomGeneral(n, n, stride, rnd)
    let wr = ArrayList<Float64>()
    let wi = ArrayList<Float64>()

    for (i in 0..t.rows) {
        for (j in 0..=i) {
            t.data[i * t.stride + j] = 0.0
        }
    }
    var i = 0
    while (i < t.rows) {
        var a = rnd.nextGaussianFloat64()
        if (bad && rnd.nextFloat64() < 0.5) {
            if (rnd.nextFloat64() < 0.5) {
                a = dlamchS
            } else {
                a = 0.0
            }
        }

        if (rnd.nextFloat64() < 0.5 || i == t.rows - 1) {
            t.data[i * t.stride + i] = a
            wr.append(a)
            wi.append(0.0)
            i++
            continue
        }

        let d = a
        var c = rnd.nextGaussianFloat64()
        let b = if (bad && rnd.nextFloat64() < 0.5) {
            dlamchS
        } else {
            rnd.nextGaussianFloat64()
        }
        if (signbit(b) == signbit(c)) {
            c *= -1.0
        }

        t.data[i * t.stride + i] = a
        t.data[i * t.stride + i + 1] = b
        t.data[(i + 1) * t.stride + i] = c
        t.data[(i + 1) * t.stride + i + 1] = d

        wr.appendAll([a, a])
        let im = sqrt(abs(b)) * sqrt(abs(c))
        wi.appendAll([im, -im])
        i += 2
    }
    return (t, wr.toArray(), wi.toArray())
}

func blockedUpperTriGeneral(r: Int64, c: Int64, k: Int64, l: Int64, stride: Int64, kblock: Bool, rnd: Random): General64 {
    var t = l
    if (kblock) {
        t += k
    }
    let ans = zeros(r, c, stride)
    for (i in 0..min(r, t)) {
        var v = 0.0
        while (v == 0.0) {
            v = rnd.nextGaussianFloat64()
        }
        ans.data[i * ans.stride + i + (c - t)] = v
    }
    for (i in 0..min(r, t)) {
        for (j in i + (c - t) + 1..c) {
            ans.data[i * ans.stride + j] = rnd.nextGaussianFloat64()
        }
    }
    return ans
}

func nanTriangular(uplo: Uplo, n: Int64, stride: Int64): Triangular {
    if (n < 0) {
        throw Exception("bad matrix size")
    }
    if (n == 0) {
        return Triangular(
            stride: max(1, stride),
            uplo: uplo,
            diag: NonUnitDiag
        )
    }
    if (stride < 0) {
        throw Exception("bad stride")
    }
    return Triangular(
        n: n,
        stride: stride,
        data: nanSlice((n - 1) * stride + n),
        uplo: uplo,
        diag: NonUnitDiag
    )
}

func generalOutsideAllNaN(a: General64): Bool {
    for (i in 0..a.rows - 1) {
        for (v in a.data[i * a.stride + a.cols..i * a.stride + a.stride]) {
            if (!v.isNaN()) {
                return false
            }
        }
    }

    var last = (a.rows - 1) * a.stride + a.cols
    if (a.rows == 0 || a.cols == 0) {
        last = 0
    }
    for (v in a.data[last..]) {
        if (!v.isNaN()) {
            return false
        }
    }
    return true
}

func triangularOutsideAllNaN(a: Triangular): Bool {
    if (a.uplo == Upper) {
        for (i in 0..a.n) {
            for (v in a.data[i * a.stride..i * a.stride + i]) {
                if (!v.isNaN()) {
                    return false
                }
            }
        }

        for (i in 0..a.n - 1) {
            for (v in a.data[i * a.stride + a.n..i * a.stride + a.stride]) {
                if (!v.isNaN()) {
                    return false
                }
            }
        }
    } else {
        for (i in 0..a.n - 1) {
            for (v in a.data[i * a.stride + a.n..i * a.stride + a.stride]) {
                if (!v.isNaN()) {
                    return false
                }
            }
        }
    }

    for (v in a.data[max(0, a.n - 1) * a.stride + a.n..]) {
        if (!v.isNaN()) {
            return false
        }
    }
    return true
}

func transposeGeneral(a: General64): General64 {
    let ans = General64(
        rows: a.cols,
        cols: a.rows,
        stride: a.rows,
        data: Array<Float64>(a.cols * a.rows, item: 0.0)
    )
    for (i in 0..a.rows) {
        for (j in 0..a.cols) {
            ans.data[j * ans.stride + i] = a.data[i * a.stride + j]
        }
    }
    return ans
}

func columnNorms(m: Int64, n: Int64, a: Array<Float64>, lda: Int64): Array<Float64> {
    let nFloat64 = nFloat64Implementation()
    return Array<Float64>(n, {j => nFloat64.dnrm2(m, a[j..], lda)})
}

func extractVMat(m: Int64, n: Int64, a: Array<Float64>, lda: Int64, direct: Direct, store: StoreV): General64 {
    let k = min(m, n)
    if (direct == Forward && store == ColumnWise) {
        let v = General64(
            rows: m,
            cols: k,
            stride: k,
            data: Array<Float64>(m * k, item: 0.0)
        )
        for (i in 0..k) {
            for (j in 0..i) {
                v.data[j * v.stride + i] = 0.0
            }
            v.data[i * v.stride + i] = 1.0
            for (j in i + 1..m) {
                v.data[j * v.stride + i] = a[j * lda + i]
            }
        }
        return v
    } else if (direct == Forward && store == RowWise) {
        let v = General64(
            rows: k,
            cols: n,
            stride: n,
            data: Array<Float64>(k * n, item: 0.0)
        )
        for (i in 0..k) {
            for (j in 0..i) {
                v.data[i * v.stride + j] = 0.0
            }
            v.data[i * v.stride + i] = 1.0
            for (j in i + 1..n) {
                v.data[i * v.stride + j] = a[i * lda + j]
            }
        }
        return v
    } else {
        // TODO
        throw Exception("not implemented")
    }
}

func constructBidiagonal(uplo: Uplo, n: Int64, d: Array<Float64>, e: Array<Float64>): General64 {
    let bMat = General64(
        rows: n,
        cols: n,
        stride: n,
        data: Array<Float64>(n * n, item: 0.0)
    )

    for (i in 0..n - 1) {
        bMat.data[i * bMat.stride + i] = d[i]
        if (uplo == Upper) {
            bMat.data[i * bMat.stride + i + 1] = e[i]
        } else {
            bMat.data[(i + 1) * bMat.stride + i] = e[i]
        }
    }
    bMat.data[(n - 1) * bMat.stride + n - 1] = d[n - 1]
    return bMat
}

func constructVMat(vMat: General64, store: StoreV, direct: Direct): General64 {
    let m = vMat.rows
    let k = vMat.cols
    if (store == ColumnWise && direct == Forward) {
        let ldv = k
        let v = Array<Float64>(m * k, item: 0.0)
        for (i in 0..m) {
            for (j in 0..k) {
                if (j > i) {
                    v[i * ldv + j] = 0.0
                } else if (j == i) {
                    v[i * ldv + i] = 1.0
                } else {
                    v[i * ldv + j] = vMat.data[i * vMat.stride + j]
                }
            }
        }
        return General64(
            rows: m,
            cols: k,
            stride: k,
            data: v
        )
    } else if (store == RowWise && direct == Forward) {
        let ldv = m
        let v = Array<Float64>(m * k, item: 0.0)
        for (i in 0..m) {
            for (j in 0..k) {
                if (j > i) {
                    v[j * ldv + i] = 0.0
                } else if (j == i) {
                    v[j * ldv + i] = 1.0
                } else {
                    v[j * ldv + i] = vMat.data[i * vMat.stride + j]
                }
            }
        }
        return General64(
            rows: k,
            cols: m,
            stride: m,
            data: v
        )
    } else if (store == ColumnWise && direct == Backward) {
        let rowsv = m
        let ldv = k
        let v = Array<Float64>(m * k, item: 0.0)
        for (i in 0..m) {
            for (j in 0..k) {
                let vrow = rowsv - i - 1
                let vcol = k - j - 1
                if (j > i) {
                    v[vrow * ldv + vcol] = 0.0
                } else if (j == i) {
                    v[vrow * ldv + vcol] = 1.0
                } else {
                    v[vrow * ldv + vcol] = vMat.data[i * vMat.stride + j]
                }
            }
        }
        return General64(
            rows: rowsv,
            cols: ldv,
            stride: ldv,
            data: v
        )
    } else if (store == RowWise && direct == Backward) {
        let rowsv = k
        let ldv = m
        let v = Array<Float64>(m * k, item: 0.0)
        for (i in 0..m) {
            for (j in 0..k) {
                let vcol = ldv - i - 1
                let vrow = k - j - 1
                if (j > i) {
                    v[vrow * ldv + vcol] = 0.0
                } else if (j == i) {
                    v[vrow * ldv + vcol] = 1.0
                } else {
                    v[vrow * ldv + vcol] = vMat.data[i * vMat.stride + j]
                }
            }
        }
        return General64(
            rows: rowsv,
            cols: ldv,
            stride: ldv,
            data: v
        )
    } else {
        // TODO
        throw Exception("not implemented")
    }
}

func constructH(tau: Array<Float64>, v: General64, store: StoreV, direct: Direct): General64 {
    let (m, k) = if (store == RowWise) {
        (v.cols, v.rows)
    } else {
        (v.rows, v.cols)
    }
    let h = General64(
        rows: m,
        cols: m,
        stride: m,
        data: Array<Float64>(m * m, item: 0.0)
    )
    for (i in 0..m) {
        h.data[i * m + i] = 1.0
    }
    for (i in 0..k) {
        let vecData = Array<Float64>(m, item: 0.0)
        if (store == ColumnWise) {
            for (j in 0..m) {
                vecData[j] = v.data[j * v.cols + i]
            }
        } else {
            for (j in 0..m) {
                vecData[j] = v.data[i * v.cols + j]
            }
        }
        let vec = Vector(
            inc: 1,
            data: vecData
        )

        let hi = General64(
            rows: m,
            cols: m,
            stride: m,
            data: Array<Float64>(m * m, item: 0.0)
        )
        for (i in 0..m) {
            hi.data[i * m + i] = 1.0
        }
        ger(-tau[i], vec, vec, hi)

        let hcopy = General64(
            rows: m,
            cols: m,
            stride: m,
            data: Array<Float64>(m * m, item: 0.0)
        )
        copy(hcopy.data, h.data)
        if (direct == Forward) {
            gemm(NoTrans, NoTrans, 1.0, hcopy, hi, 0.0, h)
        } else {
            gemm(NoTrans, NoTrans, 1.0, hi, hcopy, 0.0, h)
        }
    }
    return h
}

func constructQ(kind: String, m: Int64, n: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>): General64 {
    let k = min(m, n)
    return constructQK(kind, m, n, k, a, lda, tau)
}

func constructQK(kind: String, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, tau: Array<Float64>): General64 {
    let sz = if (kind == "QR") {
        m
    } else if (kind == "LQ" || kind == "RQ") {
        n
    } else {
        0
    }

    let q = General64(
        rows: sz,
        cols: sz,
        stride: max(1, sz),
        data: Array<Float64>(sz * sz, item: 0.0)
    )
    for (i in 0..k) {
        q.data[i * sz + i] = 1.0
    }
    let qCopy = General64(
        rows: q.rows,
        cols: q.cols,
        stride: q.stride,
        data: Array<Float64>(q.data.size, item: 0.0)
    )
    for (i in 0..k) {
        let h = General64(
            rows: sz,
            cols: sz,
            stride: max(1, sz),
            data: Array<Float64>(sz * sz, item: 0.0)
        )
        for (j in 0..sz) {
            h.data[j * sz + j] = 1.0
        }
        let vVec = Vector(
            inc: 1,
            data: Array<Float64>(sz, item: 0.0)
        )
        if (kind == "QR") {
            vVec.data[i] = 1.0
            for (j in i + 1..sz) {
                vVec.data[j] = a[j * lda + i]
            }
        } else if (kind == "LQ") {
            vVec.data[i] = 1.0
            for (j in i + 1..sz) {
                vVec.data[j] = a[i * lda + j]
            }
        } else if (kind == "RQ") {
            for (j in 0..n - k + i) {
                vVec.data[j] = a[(m - k + i) * lda + j]
            }
            vVec.data[n - k + i] = 1.0
        }

        ger(-tau[i], vVec, vVec, h)
        copy(qCopy.data, q.data)
        if (kind == "QR" || kind == "RQ") {
            gemm(NoTrans, NoTrans, 1.0, qCopy, h, 0.0, q)
        } else if (kind == "LQ") {
            gemm(NoTrans, NoTrans, 1.0, h, qCopy, 0.0, q)
        }
    }
    return q
}

func checkBidiagonal(m: Int, n: Int, nb: Int, a: Array<Float64>, lda: Int, d: Array<Float64>,
    e: Array<Float64>, tauP: Array<Float64>, tauQ: Array<Float64>, aCopy: Array<Float64>): Unit {
    let qMat = constructQPBidiagonal(ApplyQ, m, n, nb, a, lda, tauQ)
    let pMat = constructQPBidiagonal(ApplyP, m, n, nb, a, lda, tauP)

    var aMat = General64(
        rows: m,
        cols: n,
        stride: lda,
        data: Array<Float64>(aCopy.size, item: 0.0)
    )
    copy(aMat.data, aCopy)

    var tmp1 = General64(
        rows: m,
        cols: n,
        stride: n,
        data: Array<Float64>(m * n, item: 0.0)
    )
    gemm(Trans, NoTrans, 1.0, qMat, aMat, 0.0, tmp1)
    var tmp2 = General64(
        rows: m,
        cols: n,
        stride: n,
        data: Array<Float64>(m * n, item: 0.0)
    )
    gemm(NoTrans, NoTrans, 1.0, tmp1, pMat, 0.0, tmp2)

    var correctDiag = true
    var matchD = true
    var matchE = true
    for (i in 0..m) {
        for (j in 0..n) {
            if (i >= nb && j >= nb) {
                continue
            }
            let v = tmp2.data[i * tmp2.stride + j]
            if (i == j) {
                if (abs(d[i] - v) > 1e-12) {
                    matchD = false
                }
                continue
            }
            if (m >= n && i == j - 1) {
                if (abs(e[j - 1] - v) > 1e-12) {
                    matchE = false
                }
                continue
            }
            if (m < n && i - 1 == j) {
                if (abs(e[i - 1] - v) > 1e-12) {
                    matchE = false
                }
                continue
            }
            if (abs(v) > 1e-12) {
                correctDiag = false
            }
        }
    }

    if (!correctDiag) {
        throw Exception("Updated A not bi-diagonal")
    }
    if (!matchD) {
        throw Exception("D Mismatch")
    }
    if (!matchE) {
        throw Exception("E mismatch")
    }
}

public func constructQPBidiagonal(vect: ApplyOrtho, m: Int, n: Int, nb: Int, a: Array<Float64>, lda: Int,
    tau: Array<Float64>): General64 {
    var sz = n
    if (vect == ApplyQ) {
        sz = m
    }

    var ldv: Int
    var v: General64
    if (vect == ApplyQ) {
        ldv = nb
        v = General64(
            rows: m,
            cols: nb,
            stride: ldv,
            data: Array<Float64>(m * ldv, item: 0.0)
        )
    } else {
        ldv = n
        v = General64(
            rows: nb,
            cols: n,
            stride: ldv,
            data: Array<Float64>(m * ldv, item: 0.0)
        )
    }

    if (vect == ApplyQ) {
        if (m >= n) {
            for (i in 0..m) {
                for (j in 0..min(nb - 1, i)) {
                    if (i == j) {
                        v.data[i * ldv + j] = 1.0
                        continue
                    }
                    v.data[i * ldv + j] = a[i * lda + j]
                }
            }
        } else {
            for (i in 1..m) {
                for (j in 0..min(nb - 1, i - 1)) {
                    if (i - 1 == j) {
                        v.data[i * ldv + j] = 1.0
                        continue
                    }
                    v.data[i * ldv + j] = a[i * lda + j]
                }
            }
        }
    } else {
        if (m < n) {
            for (i in 0..nb) {
                for (j in i..n) {
                    if (i == j) {
                        v.data[i * ldv + j] = 1.0
                        continue
                    }
                    v.data[i * ldv + j] = a[i * lda + j]
                }
            }
        } else {
            for (i in 0..nb) {
                for (j in i + 1..n) {
                    if (j - 1 == i) {
                        v.data[i * ldv + j] = 1.0
                        continue
                    }
                    v.data[i * ldv + j] = a[i * lda + j]
                }
            }
        }
    }

    var qMat = General64(
        rows: sz,
        cols: sz,
        stride: sz,
        data: Array<Float64>(sz * sz, item: 0.0)
    )
    var hMat = General64(
        rows: sz,
        cols: sz,
        stride: sz,
        data: Array<Float64>(sz * sz, item: 0.0)
    )

    for (i in 0..sz) {
        qMat.data[i * qMat.stride + i] = 1.0
    }

    for (i in 0..nb) {
        var qCopy = General64(rows: qMat.rows, cols: qMat.cols, stride: qMat.stride,
            data: Array<Float64>(qMat.data.clone()))
        copy(qCopy.data, qMat.data)

        for (i in 0..sz) {
            for (j in 0..sz) {
                if (i == j) {
                    hMat.data[i * sz + j] = 1.0
                } else {
                    hMat.data[i * sz + j] = 0.0
                }
            }
        }

        var vi: Vector

        if (vect == ApplyQ) {
            vi = Vector(
                inc: v.stride,
                data: v.data[i..]
            )
        } else {
            vi = Vector(
                inc: 1,
                data: v.data[i * v.stride..]
            )
        }

        ger(-tau[i], vi, vi, hMat)

        gemm(NoTrans, NoTrans, 1.0, qCopy, hMat, 0.0, qMat)
    }

    return qMat
}

//TODO func printRowise(a []float64, m, n, lda int, beyond bool) {

//TODO func copyGeneral(dst: General64, src: General64): Unit {  //TODO

//TODO func cloneGeneral(a: General64): General64 {  //TODO

//TODO func equalGeneral(a, b blas64.General) bool {

func equalApproxGeneral(a: General64, b: General64, tol: Float64): Bool {
    if (a.rows != b.rows || a.cols != b.cols) {
        throw Exception("bad input")
    }
    for (i in 0..a.rows) {
        for (j in 0..a.cols) {
            let diff = a.data[i * a.stride + j] - b.data[i * b.stride + j]
            if (diff.isNaN() || abs(diff) > tol) {
                return false
            }
        }
    }
    return true
}

//TODO func intsEqual(a, b []int) bool {

//TODO func randSymBand(uplo blas.Uplo, n, kd, ldab int, rnd *rand.Rand) []float64 {

//TODO func distSymBand(uplo blas.Uplo, n, kd int, a []float64, lda int, b []float64, ldb int) float64 {

func eye(n: Int, stride: Int): General64 {
    let ans = nanGeneral(n, n, stride)
    for (i in 0..n) {
        for (j in 0..n) {
            ans.data[i * ans.stride + j] = 0.0
        }
        ans.data[i * ans.stride + i] = 1.0
    }
    return ans
}

func zeros(m: Int64, n: Int64, stride: Int64): General64 {
    let a = nanGeneral(m, n, stride)
    for (i in 0..m) {
        for (j in 0..n) {
            a.data[i * a.stride + j] = 0.0
        }
    }
    return a
}

//TODO func extract2x2Block(t []float64, ldt int) (a, b, c, d float64) {

//TODO func isSchurCanonical(a, b, c, d float64) bool {

//TODO func isSchurCanonicalGeneral(t blas64.General) bool {

//TODO func schurBlockEigenvalues(a, b, c, d float64) (ev1, ev2 complex128) {

//TODO func schurBlockSize(t blas64.General, i int) (size int, first bool) {

//TODO func containsComplex(v: Array<Complex128>, z: Complex128, tol: Float64): (Bool, Int) {  //TODO

//TODO func isAllNaN(x: Array<Float64>): Bool {  //TODO

//TODO func isUpperHessenberg(h: General64): Bool {  //TODO

//TODO func isUpperTriangular(a blas64.General) bool {

//TODO func unbalancedSparseGeneral(m, n, stride int, nonzeros int, rnd *rand.Rand) blas64.General {

//TODO func rootsOfUnity(n int) []complex128 {

//TODO func constructGSVDresults(n, p, m, k, l int, a, b blas64.General, alpha, beta []float64) (zeroR, d1, d2 blas64.General) {

//TODO func constructGSVPresults(n, p, m, k, l int, a, b blas64.General) (zeroA, zeroB blas64.General) {

func distFromIdentity(n: Int64, a: Array<Float64>, lda: Int64): Float64 {
    var dist = 0.0
    for (i in 0..n) {
        for (j in 0..n) {
            let aij = a[i * lda + j]
            if (aij.isNaN()) {
                return Float64.Inf
            }
            if (i == j) {
                dist = max(dist, abs(aij - 1.0))
            } else {
                dist = max(dist, abs(aij))
            }
        }
    }
    return dist
}

//TODO func sameFloat64(a, b float64) bool {

//TODO func sameLowerTri(n int, a []float64, lda int, b []float64, ldb int) bool {

//TODO func sameUpperTri(n int, a []float64, lda int, b []float64, ldb int) bool {

//TODO func svdJobString(job lapack.SVDJob) string {

func rCondTestRatio(rcond: Float64, rcondc: Float64): Float64 {
    const eps = dlamchE
    return if (rcond > 0.0 && rcondc > 0.0) {
        max(rcond, rcondc) / min(rcond, rcondc) - (1.0 - eps)
    } else if (rcond > 0.0) {
        rcond / eps
    } else if (rcondc > 0.0) {
        rcondc / eps
    } else {
        0.0
    }
}

func residualOrthogonal(q: General64, rowwise: Bool): Float64 {
    let (m, n) = (q.rows, q.cols)
    if (m == 0 || n == 0) {
        return 0.0
    }
    let transq = if (m < n || (m == n && rowwise)) {
        NoTrans
    } else {
        Trans
    }
    let minmn = min(m, n)

    let work = Symmetric(
        uplo: Upper,
        n: minmn,
        data: Array<Float64>(minmn * minmn, item: 0.0),
        stride: minmn,
    )
    for (i in 0..minmn) {
        work.data[i * work.stride + i] = 1.0
    }

    syrk(transq, -1.0, q, 1.0, work)
    return dlansy(MaxColumnSum, Upper, work.n, work.data, work.stride)
}

func perm(n: Int64, rnd: Random): Array<Int64> {
    let a = Array<Int64>(n, {i => i})

    var i = n - 1
    while (i > 0) {
        let j = rnd.nextInt64(i + 1)

        let tmp = a[i]
        a[i] = a[j]
        a[j] = tmp
        i--
    }

    return a
}

func copysign(f: Float64, sign: Float64): Float64 {
    const signBit: UInt64 = 1 << 63
    return Float64.fromBits((f.toBits() & (!signBit)) | sign.toBits() & signBit)
}

func signbit(x: Float64): Bool {
    return (x.toBits() & (1 << 63)) != 0
}

func real(x: Complex128): Float64 {
    return x.real
}

func hypot(p: Float64, q: Float64): Float64 {
    var p_ = abs(p)
    var q_ = abs(q)

    if (p_ == Float64.Inf || q_ == Float64.Inf) {
        return Float64.Inf
    }

    if (p_.isNaN() || q_.isNaN()) {
        return Float64.NaN
    }

    if (p_ < q_) {
        let tmp = p_
        p_ = q_
        q_ = tmp
    }

    if (p_ == 0.0) {
        return 0.0
    }

    q_ = q_ / p_

    return p_ * sqrt(1.0 + pow(q_, 2))
}
