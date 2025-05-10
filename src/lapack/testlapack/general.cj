package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64
import cjnum.lapack.*

const dlamchE = 0x1p-53
const dlamchB = 2.0
const dlamchP = dlamchB * dlamchE

const dlamchS = 0x1p-1022

const safmin = dlamchS
const safmax = 1.0 / safmin
const ulp    = dlamchP
const smlnum = safmin / ulp
const bignum = safmax * ulp

type worklen = Int64

const minimumWork: worklen = 0
const mediumWork:  worklen = 1
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

func nanSlice(n : Int64): Array<Float64> {
    return Array<Float64>(n, { _ => 
        Float64.NaN
    })
}

func randomSlice(n: Int64, rnd: Random): Array<Float64> {
    return Array<Float64>(n, { _ =>
        rnd.nextFloat64()
    })
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

func zeros(m: Int64, n: Int64, stride: Int64): General64 {
    let a = nanGeneral(m, n, stride)
    for (i in 0..m) {
        for (j in 0..n) {
            a.data[i * a.stride + j] = 0.0
        }
    }
    return a
}

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

func perm(n: Int64, rnd: Random): Array<Int64> {
    let a = Array<Int64>(n, { i => i })

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
