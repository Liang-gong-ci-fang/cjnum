package cjnum.lapack.testlapack

import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.blas.blas64.*

public func dtrtriTest(lFloat64: LFloat64): Bool {
    const tol = 1e-10
    let rnd = Random(1)
    let nFloat64 = nFloat64Implementation()
    for (uplo in [Upper, Lower]) {
        for (diag in [NonUnitDiag, UnitDiag]) {
            for (test in [
                [3, 0], [70, 0], [200, 0], [3, 5], [70, 92], [200, 205]
            ]) {
                let n = test[0]
                let lda = if (test[1] == 0) {
                    n
                } else {
                    test[1]
                }
                let a = Array<Float64>(n * lda, { _ => rnd.nextFloat64() })
                for (i in 0..n) {
                    a[i * lda + i] += Float64(n)
                }

                let aCopy = a.clone()
                lFloat64.dtrtri(uplo, diag, n, a, lda)
                if (uplo == Upper) {
                    for (i in 1..n) {
                        for (j in 0..i) {
                            aCopy[i * lda + j] = 0.0
                            a[i * lda + j] = 0.0
                        }
                    }
                } else {
                    for (i in 0..n) {
                        for (j in i + 1..n) {
                            aCopy[i * lda + j] = 0.0
                            a[i * lda + j] = 0.0
                        }
                    }
                }
                if (diag == UnitDiag) {
                    for (i in 0..n) {
                        a[i * lda + i] = 1.0
                        aCopy[i * lda + i] = 1.0
                    }
                }
                let ans = Array<Float64>(a.size, repeat: 0.0)
                nFloat64.dgemm(NoTrans, NoTrans, n, n, n, 1.0, a, lda, aCopy, lda, 0.0, ans, lda)
                let dist = distFromIdentity(n, ans, lda)
                if (dist > tol) {
                    throw Exception("|inv(A) * A - I| = ${dist} is too large. Upper = ${uplo == Upper}, unit = ${diag == UnitDiag}, n = ${n}, lda = ${lda}")
                }
            }
        }
    }
    return true
}

