package cjnum.lapack.testlapack

import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlaexcTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [1, 2, 3, 4, 5, 6, 10, 18, 31, 53]) {
        for (extra in [0, 3]) {
            for (_ in 0..100) {
                if (!testDlaexc(lFloat64, rnd, n, extra)) {
                    return false
                }
            }
        }
    }
    return true
}

func testDlaexc(lFloat64: LFloat64, rnd: Random, n: Int, extra: Int): Bool {
    const tol = 1e-14

    let (tmat, _, _) = randomSchurCanonical(n, n + extra, true, rnd)
    let tmatCopy = cloneGeneral(tmat)

    var j1 = rnd.nextInt64(n)
    if (j1 > 0 && tmat.data[j1 * tmat.stride + j1 - 1] != 0.0) {
        j1--
    }

    var n1: Int
    var n2: Int
    match {
        case j1 == n - 1 => (n1, n2) = (1, 0)
        case j1 == n - 2 => if (tmat.data[(j1 + 1) * tmat.stride + j1] == 0.0) {
            (n1, n2) = (1, 1)
        } else {
            (n1, n2) = (2, 0)
        }
        case j1 == n - 3 => if (tmat.data[(j1 + 1) * tmat.stride + j1] == 0.0) {
            (n1, n2) = (1, 2)
        } else {
            (n1, n2) = (2, 1)
        }
        case _ => if (tmat.data[(j1 + 1) * tmat.stride + j1] == 0.0) {
            n1 = 1
            if (tmat.data[(j1 + 2) * tmat.stride + j1 + 1] == 0.0) {
                n2 = 1
            } else {
                n2 = 2
            }
        } else {
            n1 = 2
            if (tmat.data[(j1 + 3) * tmat.stride + j1 + 2] == 0.0) {
                n2 = 1
            } else {
                n2 = 2
            }
        }
    }

    let name = "Case n=${n},j1=${j1},n1=${n1},n2=${n2},extra=${extra}"

    var wantq = false

    var work = nanSlice(n)

    let ok = lFloat64.dlaexc(wantq, n, tmat.data, tmat.stride, [], 1, j1, n1, n2, work)

    wantq = true

    let tmat2 = cloneGeneral(tmatCopy)
    let q = eye(n, n + extra)
    let qCopy = cloneGeneral(q)
    work = nanSlice(n)

    let ok2 = lFloat64.dlaexc(wantq, n, tmat2.data, tmat2.stride, q.data, q.stride, j1, n1, n2, work)

    if (!generalOutsideAllNaN(tmat)) {
        throw Exception("${name}: out-of-range write to T")
    }
    if (!generalOutsideAllNaN(tmat2)) {
        throw Exception("${name}v: out-of-range write to T2")
    }
    if (!generalOutsideAllNaN(q)) {
        throw Exception("${name}: out-of-range write to Q")
    }

    if (ok != ok2) {
        throw Exception("${name}: ok != ok2")
    }
    if (!equalGeneral(tmat, tmat2)) {
        throw Exception("${name}: T != T2")
    }

    if (!ok) {
        if (n1 == 1 && n2 == 1) {
            throw Exception("${name}: unexpected failure")
        } else {
            // t.Logf("%v: Dlaexc returned false", name)
            println("${name}: Dlaexc returned false")
        }
    }

    if (!ok || n1 == 0 || n2 == 0 || j1 + n1 >= n) {
        if (!equalGeneral(tmat, tmatCopy)) {
            throw Exception("${name}: unexpected modification of T")
        }

        if (!equalGeneral(q, qCopy)) {
            throw Exception("${name}: unexpected modification of Q")
        }
        return true
    }

    for (i in 0..n) {
        if (j1 <= i && i < j1 + n1 + n2) {
            continue
        }
        for (j in 0..n) {
            if (j1 <= j && j < j1 + n1 + n2) {
                continue
            }
            let diff = tmat.data[i * tmat.stride + j] - tmatCopy.data[i * tmatCopy.stride + j]
            if (diff != 0.0) {
                throw Exception("${name}: unexpected modification of T[${i},${j}]")
            }
        }
    }

    if (!isSchurCanonicalGeneral(tmat)) {
        throw Exception("${name}: T is not in Schur canonical form")
    }

    var resid = residualOrthogonal(q, false)
    if (resid > tol) {
        throw Exception("${name}: Q is not orthogonal; resid=${resid}, want<=${tol}")
    }

    for (i in 0..n) {
        for (j in 0..n) {
            if (j1 <= j && j < j1 + n1 + n2) {
                continue
            }
            let diff = q.data[i * q.stride + j] - qCopy.data[i * qCopy.stride + j]
            if (diff != 0.0) {
                throw Exception("${name}: unexpected modification of Q[${i},${j}]")
            }
        }
    }

    let qt = zeros(n, n, n)
    gemm(Trans, NoTrans, 1.0, q, tmatCopy, 0.0, qt)
    let qtq = cloneGeneral(tmat)
    gemm(NoTrans, NoTrans, -1.0, qt, q, 1.0, qtq)
    resid = dlange(MaxColumnSum, n, n, qtq.data, qtq.stride)
    if (resid > Float64(n) * tol) {
        throw Exception(
            "${name}: mismatch between Qᵀ*(initial T)*Q and (final T); resid=${resid}, want<=${Float64(n)*tol}")
    }
    return true
}
