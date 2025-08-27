package cjnum.lapack.testlapack

import std.random.*
import cjnum.lapack.*
import cjnum.blas.*

public func dlasetTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    let tests = [
		[0, 0],
		[1, 1],
		[1, 10],
		[10, 1],
		[2, 2],
		[2, 10],
		[10, 2],
		[11, 11],
		[11, 100],
		[100, 11]
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let m = test[0]
        let n = test[1]

        for (uplo in [Upper, Lower, All]) {
            for (extra in [0, 10]) {
                let a = randomGeneral(m, n, n + extra, rnd)
                let alpha = 1.0
                let beta = 2.0

                lFloat64.dlaset(uplo, m, n, alpha, beta, a.data, a.stride)

                let prefix = "Case: ${i}: m=${m},n=${n},uplo=${uplo},extra=${extra}"

                if (!generalOutsideAllNaN(a)) {
                    throw Exception("${prefix}: out-of-range write to A")
                }
                for (i in 0..min(m, n)) {
                    if (a.data[i * a.stride + i] != beta) {
                        throw Exception("${prefix}: unexpected diagonal of A")
                    }
                }
                if (uplo == Upper || uplo == All) {
                    for (i in 0..m) {
                        for (j in i + 1..n) {
                            if (a.data[i * a.stride + j] != alpha) {
                                throw Exception("${prefix}: unexpected upper triangle of A")
                            }
                        }
                    }
                }
                if (uplo == Lower || uplo == All) {
                    for (i in 1..m) {
                        for (j in 0..min(i, n)) {
                            if (a.data[i * a.stride + j] != alpha) {
                                throw Exception("${prefix}: unexpected lower triangle of A")
                            }
                        }
                    }
                }
            }
        }
    }
    return true
}