package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*

public func dlacpyTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (uplo in [Upper, Lower, All]) {
        for (test in [
			[3, 5, 0, 0],
			[5, 5, 0, 0],
			[7, 5, 0, 0],

			[3, 5, 10, 12],
			[5, 5, 10, 12],
			[7, 5, 10, 12]
        ]) {
            let m = test[0]
            let n = test[1]
            let lda = if (test[2] == 0) {
                n
            } else {
                test[2]
            }
            let ldb = if (test[3] == 0) {
                n
            } else {
                test[3]
            }

            let a = Array<Float64>(m * lda, { _ => rnd.nextFloat64() })
            let b = Array<Float64>(m * ldb, { _ => rnd.nextFloat64() })

            lFloat64.dlacpy(uplo, m, n, a, lda, b, ldb)
            var equal = true
            if (uplo == Upper) {
                for (i in 0..m) {
                    for (j in i..n) {
                        if (b[i * ldb + j] != a[i * lda + j]) {
                            equal = false
                            doneCheck(equal, uplo, m, n)
                        }
                    }
                }
            } else if (uplo == Lower) {
                for (i in 0..m) {
                    for (j in 0..min(i, n)) {
                        if (b[i * ldb + j] != a[i * lda + j]) {
                            equal = false
                            doneCheck(equal, uplo, m, n)
                        }
                    }
                }
            } else if (uplo == All) {
                for (i in 0..m) {
                    for (j in 0..n) {
                        if (b[i * ldb + j] != a[i * lda + j]) {
                            equal = false
                            doneCheck(equal, uplo, m, n)
                        }
                    }
                }
            }
        }
    }
    return true
}

func doneCheck(equal: Bool, uplo: Uplo, m: Int64, n: Int64): Unit {
    if (!equal) {
        throw Exception("Matrices not equal after copy. Uplo = ${uplo}, m = ${m}, n = ${n}")
    }
}