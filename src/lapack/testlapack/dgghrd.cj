package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dgghrdTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    let comps = [OrthoExplicit, OrthoNone, OrthoPostmul]
    for (compq in 0..comps.size) {
        for (compz in 0..comps.size) {
            for (n in [0, 1, 2, 3, 4, 15]) {
                for (ld in [max(1, n), n + 5]) {
                    if (!testDgghrd(lFloat64, rnd, comps[compq], comps[compz], n, 0, n - 1, ld, ld, ld, ld)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func testDgghrd(lFloat64: LFloat64, rnd: Random, compq: OrthoComp, compz: OrthoComp, n: Int, ilo: Int, ihi: Int,
    lda: Int, ldb: Int, ldq: Int, ldz: Int): Bool {
    let tol = 1e-13

    let a = randomGeneral(n, n, lda, rnd)
    let b = randomGeneral(n, n, ldb, rnd)

    var q: General64 = General64()
    var q1: General64 = General64()
    if (compq == OrthoExplicit) {
        q = randomGeneral(n, n, ldq, rnd)
    } else if (compq == OrthoPostmul) {
        q = randomOrthogonal(n, rnd)
        q1 = cloneGeneral(q)
    }

    var z: General64 = General64()
    var z1: General64 = General64()
    if (compz == OrthoExplicit) {
        z = randomGeneral(n, n, ldz, rnd)
    } else if (compz == OrthoPostmul) {
        z = randomOrthogonal(n, rnd)
        z1 = cloneGeneral(z)
    }

    let hGot = cloneGeneral(a)
    let tGot = cloneGeneral(b)
    lFloat64.dgghrd(compq, compz, n, ilo, ihi, hGot.data, hGot.stride, tGot.data, tGot.stride, q.data, max(1, q.stride),
        z.data, max(1, z.stride))

    if (n == 0) {
        return true
    }

    let name = "Case compq=${compq},compz=${compz},n=${n},ilo=${ilo},ihi=${ihi}"

    if (!isUpperHessenberg(hGot)) {
        throw Exception("${name}: H is not upper Hessenberg")
    }
    if (!isUpperTriangular(tGot)) {
        throw Exception("${name}: T is not upper triangular")
    }
    if (compq != OrthoNone) {
        let resid = residualOrthogonal(q, true)
        if (resid > tol) {
            throw Exception("${name}: Q is not orthogonal, resid=${resid}")
        }
    }
    if (compz != OrthoNone) {
        let resid = residualOrthogonal(z, true)
        if (resid > tol) {
            throw Exception("${name}: Z is not orthogonal, resid=${resid}")
        }
    }

    if (compq != compz) {
        return true
    }

    for (i in 1..n) {
        for (j in 0..i) {
            b.data[i * b.stride + j] = 0.0
        }
    }

    let aux = zeros(n, n, n)
    if (compq == OrthoExplicit) {
        let hCalc = zeros(n, n, n)
        gemm(Trans, NoTrans, 1.0, q, a, 0.0, aux)
        gemm(NoTrans, NoTrans, 1.0, aux, z, 1.0, hCalc)
        if (!equalApproxGeneral(hGot, hCalc, tol)) {
            throw Exception("${name}: Qᵀ*A*Z != H")
        }

        let tCalc = zeros(n, n, n)
        gemm(Trans, NoTrans, 1.0, q, b, 0.0, aux)
        gemm(NoTrans, NoTrans, 1.0, aux, z, 1.0, tCalc)
        if (!equalApproxGeneral(tGot, tCalc, tol)) {
            throw Exception("${name}: Qᵀ*B*Z != T")
        }
    } else if (compq == OrthoPostmul) {
        let lhs = zeros(n, n, n)
        gemm(NoTrans, NoTrans, 1.0, q1, a, 0.0, aux)
        gemm(NoTrans, Trans, 1.0, aux, z1, 0.0, lhs)

        let rhs = zeros(n, n, n)
        gemm(NoTrans, NoTrans, 1.0, q, hGot, 0.0, aux)
        gemm(NoTrans, Trans, 1.0, aux, z, 0.0, rhs)
        if (!equalApproxGeneral(lhs, rhs, tol)) {
            throw Exception("${name}: Q1 * A * Z1ᵀ != (Q1*Q) * H * (Z1*Z)ᵀ")
        }

        gemm(NoTrans, NoTrans, 1.0, q1, b, 0.0, aux)
        gemm(NoTrans, Trans, 1.0, aux, z1, 0.0, lhs)

        gemm(NoTrans, NoTrans, 1.0, q, tGot, 0.0, aux)
        gemm(NoTrans, Trans, 1.0, aux, z, 0.0, rhs)
        if (!equalApproxGeneral(lhs, rhs, tol)) {
            throw Exception("${name}: Q1 * B * Z1ᵀ != (Q1*Q) * T * (Z1*Z)ᵀ")
        }
    }
    return true
}
