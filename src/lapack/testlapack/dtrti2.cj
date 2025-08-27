package cjnum.lapack.testlapack

import std.random.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.blas.blas64.*

class dtrti2TestStruct {
    dtrti2TestStruct(
        let a!: Array<Float64>,
        let n!: Int64,
        let uplo!: Uplo,
        let diag!: Diag,
        let ans!: Array<Float64>
    ) { }
}

public func dtrti2Test(lFloat64: LFloat64): Bool {
    const tol = 1e-14
    let tests = [
        dtrti2TestStruct(
            a: [2.0, 3.0, 4.0, 0.0, 5.0, 6.0, 8.0, 0.0, 8.0],
            n: 3,
            uplo: Upper,
            diag: NonUnitDiag,
            ans: [0.5, -0.3, -0.025, 0.0, 0.2, -0.15, 8.0, 0.0, 0.125]
        ),
        dtrti2TestStruct(
            a: [5.0, 3.0, 4.0, 0.0, 7.0, 6.0, 10.0, 0.0, 8.0],
            n: 3,
            uplo: Upper,
            diag: UnitDiag,
            ans: [5.0, -3.0, 14.0, 0.0, 7.0, -6.0, 10.0, 0.0, 8.0]
        ),
        dtrti2TestStruct(
            a: [2.0, 0.0, 0.0, 3.0, 5.0, 0.0, 4.0, 6.0, 8.0],
            n: 3,
            uplo: Lower,
            diag: NonUnitDiag,
            ans: [0.5, 0.0, 0.0, -0.3, 0.2, 0.0, -0.025, -0.15, 0.125]
        ),
        dtrti2TestStruct(
            a: [1.0, 0.0, 0.0, 3.0, 1.0, 0.0, 4.0, 6.0, 1.0],
            n: 3,
            uplo: Lower,
            diag: UnitDiag,
            ans: [1.0, 0.0, 0.0, -3.0, 1.0, 0.0, 14.0, -6.0, 1.0]
        )
    ]
    for (test in tests) {
        lFloat64.dtrti2(test.uplo, test.diag, test.n, test.a, test.n)
        if (!fEqualApprox(test.ans, test.a, tol)) {
            throw Exception("Matrix inverse mismatch. Want ${test.ans}, got ${test.a}.")
        }
    }
    let rnd = Random(1)
    let nFloat64 = nFloat64Implementation()
    for (uplo in [Upper, Lower]) {
        for (diag in [NonUnitDiag, UnitDiag]) {
            for (test in [
                [1, 0], [2, 0], [3, 0], [1, 5], [2, 5], [3, 5]
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
                lFloat64.dtrti2(uplo, diag, n, a, lda)
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
                    throw Exception("|inv(A) * A - I| = ${dist}. Upper = ${uplo == Upper}, unit = ${diag == UnitDiag}, ans = ${ans}")
                }
            }
        }
    }
    return true
}