package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.blas.blas64.*

public func dgetriTest(lFloat64: LFloat64): Bool {
    const tol = 1e-13
    let rnd = Random(1)
    let nFloat64 = nFloat64Implementation()
    for (test in [
		[5, 0],
		[5, 8],
		[45, 0],
		[45, 50],
		[63, 70],
		[64, 70],
		[65, 0],
		[65, 70],
		[66, 70],
		[150, 0],
		[150, 250]
    ]) {
        let n = test[0]
        let lda = if (test[1] == 0) {
            n
        } else {
            test[1]
        }
        
        let perm = perm(n, rnd)
        let a = Array<Float64>(n * lda, repeat: 0.0)
        for (i in 0..n) {
            a[i * lda + perm[i]] = 1.0
        }
        for (i in 0..a.size) {
            a[i] += 0.01 * rnd.nextFloat64()
        }
        let aCopy = a.clone()
        let ipiv = Array<Int64>(n, repeat: 0)
        lFloat64.dgetrf(n, n, a, lda, ipiv)

        for (wl in [minimumWork, mediumWork, optimumWork]) {
            let ainv = a.clone()

            var lwork = 0
            if (wl == minimumWork) {
                lwork = max(1, n)
            } else if (wl == mediumWork) {
                let work = Array<Float64>(1, repeat: 0.0)
                lFloat64.dgetri(n, ainv, lda, ipiv, work, -1)
                lwork = max(Int64(work[0]) - 2 * n, n)
            } else if (wl == optimumWork) {
                let work = Array<Float64>(1, repeat: 0.0)
                lFloat64.dgetri(n, ainv, lda, ipiv, work, -1)
                lwork = Int64(work[0])
            }

            let work = Array<Float64>(lwork, repeat: 0.0)

            let ok = lFloat64.dgetri(n, ainv, lda, ipiv, work, lwork)
            if (!ok) {
                throw Exception("Unexpected singular matrix.")
            }

            let ans = Array<Float64>(ainv.size, repeat: 0.0)
            nFloat64.dgemm(NoTrans, NoTrans, n, n, n, 1.0, aCopy, lda, ainv, lda, 0.0, ans, lda)
            let dist = distFromIdentity(n, ans, lda)
            if (dist > tol) {
                throw Exception("|Inv(A) * A - I|_inf = ${dist} is too large. n = ${n}, lda = ${lda}")
            }
        }
    }
    return true
}