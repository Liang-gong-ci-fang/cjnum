package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

// TODO dgeev.cj

func residualRightEV(a: General64, e: General64, wr: Array<Float64>, wi: Array<Float64>): Float64 {
    let n = a.rows
    if (n == 0) {
        return 0.0
    }

    let nFloat64 = nFloat64Implementation()
    let ldr = n
    let r = Array<Float64>(n * ldr, repeat: 0.0)
    let wmat = Array<Float64>(4, repeat: 0.0)
    var ipair = 0
    for (j in 0..n) {
        if (ipair == 0 && wi[j] != 0.0) {
            ipair = 1
        }
        if (ipair == 0) {
            nFloat64.daxpy(n, wr[j], e.data[j..], e.stride, r[j..], ldr)
        } else if (ipair == 1) {
            wmat[0] = wr[j]
            wmat[1] = wi[j]
            wmat[2] = -wi[j]
            wmat[3] = wr[j]
            nFloat64.dgemm(NoTrans, NoTrans, n, 2, 2, 1.0, e.data[j..], e.stride, wmat[0..], 2, 0.0, r[j..], ldr)
            ipair = 2
        } else if (ipair == 2) {
            ipair = 0
        }
    }
    nFloat64.dgemm(NoTrans, NoTrans, n, n, n, 1.0, a.data, a.stride, e.data, e.stride, -1.0, r, ldr)

    const eps = dlamchE
    let anorm = max(dlange(MaxColumnSum, n, n, a.data, a.stride), safmin)
    let enorm = max(dlange(MaxColumnSum, n, n, e.data, e.stride), eps)
    let errnorm = dlange(MaxColumnSum, n, n, r, ldr) / enorm
    if (anorm > errnorm) {
        return errnorm / anorm
    }
    if (anorm < 1.0) {
        return min(errnorm, anorm) / anorm
    }
    return min(errnorm / anorm, 1.0)
}

func residualLeftEV(a: General64, e: General64, wr: Array<Float64>, wi: Array<Float64>): Float64 {
    let n = a.rows
    if (n == 0) {
        return 0.0
    }

    let nFloat64 = nFloat64Implementation()
    let ldr = n
    let r = Array<Float64>(n * ldr, repeat: 0.0)
    let wmat = Array<Float64>(4, repeat: 0.0)
    var ipair = 0
    for (j in 0..n) {
        if (ipair == 0 && wi[j] != 0.0) {
            ipair = 1
        }
        if (ipair == 0) {
            nFloat64.daxpy(n, wr[j], e.data[j..], e.stride, r[j..], ldr)
        } else if (ipair == 1) {
            wmat[0] = wr[j]
            wmat[1] = wi[j]
            wmat[2] = -wi[j]
            wmat[3] = wr[j]
            nFloat64.dgemm(NoTrans, Trans, n, 2, 2, 1.0, e.data[j..], e.stride, wmat[0..], 2, 0.0, r[j..], ldr)
            ipair = 2
        } else if (ipair == 2) {
            ipair = 0
        }
    }
    nFloat64.dgemm(Trans, NoTrans, n, n, n, 1.0, a.data, a.stride, e.data, e.stride, -1.0, r, ldr)

    const eps = dlamchE
    let anorm = max(dlange(MaxRowSum, n, n, a.data, a.stride), safmin)
    let enorm = max(dlange(MaxColumnSum, n, n, e.data, e.stride), eps)
    let errnorm = dlange(MaxColumnSum, n, n, r, ldr) / enorm
    if (anorm > errnorm) {
        return errnorm / anorm
    }
    if (anorm < 1.0) {
        return min(errnorm, anorm) / anorm
    }
    return min(errnorm / anorm, 1.0)
}
