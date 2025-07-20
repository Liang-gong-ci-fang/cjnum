package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*

public func dlansyTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (norm in [MaxAbs, MaxColumnSum, MaxRowSum, Frobenius]) {
        for (uplo in [Lower, Upper]) {
            for (test in [
				[1, 0],
				[3, 0],

				[1, 10],
				[3, 10]
            ]) {
                for (trial in 0..100) {
                    let n = test[0]
                    let lda = if (test[1] == 0) {
                        n
                    } else {
                        test[1]
                    }

                    let a = Array<Float64>(lda * n, { i => 
                        if (trial == 0) {
                            Float64(i)
                        } else {
                            rnd.nextGaussianFloat64()
                        }
                    })

                    let aDense = Array<Float64>(n * n, repeat: 0.0)
                    if (uplo == Upper) {
                        for (i in 0..n) {
                            for (j in i..n) {
                                let v = a[i * lda + j]
								aDense[i * n + j] = v
								aDense[j * n + i] = v
                            }
                        }
                    } else {
                        for (i in 0..n) {
                            for (j in 0..=i) {
                                let v = a[i * lda + j]
								aDense[i * n + j] = v
								aDense[j * n + i] = v
                            }
                        }
                    }

                    let work = Array<Float64>(n, repeat: 0.0)
                    let got = lFloat64.dlansy(norm, uplo, n, a, lda, work)
                    let want = lFloat64.dlange(norm, n, n, aDense, n, work)
                    if (abs(want - got) > 1e-14) {
                        throw Exception("Norm mismatch. norm = ${norm}, upper = ${uplo == Upper}, n = ${n}, lda = ${lda}, want ${want}, got ${got}.")
                    }
                }
            }
        }
    }
    return true
}