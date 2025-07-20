package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.random.*
import std.math.*

public func dlauu2Test(lFloat64: LFloat64): Bool {
    for (uplo in [Upper, Lower]) {
        let ns = [0, 1, 2, 3, 4, 5, 10, 25]
        if (!dlauuTest(lFloat64.dlauu2, uplo, ns)) {
            return false
        }
    }
    return true
}

func dlauuTest(f: (uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64) -> Unit, uplo: Uplo, ns: Array<Int64>): Bool {
    const tol = 1e-13

    let nFloat64 = nFloat64Implementation()
    let rnd = Random(1)

    for (n in ns) {
        for (lda in [max(1, n), n + 11]) {
            let prefix = "n=${n},lda=${lda}"

            let a = Array<Float64>(n * lda, repeat: 0.0)
            for (i in 0..a.size) {
                a[i] = rnd.nextGaussianFloat64()
            }

            let aCopy = Array<Float64>(a.size, repeat: 0.0)
            copy(aCopy, a)

            f(uplo, n, a, lda)

            if (n == 0) {
                continue
            }

            if (uplo == Upper) {
                if (!sameLowerTri(n, aCopy, lda, a, lda)) {
                    throw Exception("${prefix}: unexpected modification in lower triangle")
                }
                for (i in 1..n) {
                    for (j in 0..i) {
                        a[i * lda + j] = a[j * lda + i]
                        aCopy[i * lda + j] = 0.0
                    }
                }
            } else {
                if (!sameUpperTri(n, aCopy, lda, a, lda)) {
                    throw Exception("${prefix}: unexpected modification in upper triangle")
                }
                for (i in 0..n - 1) {
                    for (j in i + 1..n) {
                        a[i * lda + j] = a[j * lda + i]
                        aCopy[i * lda + j] = 0.0
                    }
                }
            }

            let r = cloneGeneral(General64(rows: n, cols: n, data: a, stride: lda))
            if (uplo == Upper) {
                let u = aCopy
                let ldu = lda
                nFloat64.dgemm(NoTrans, Trans, n, n, n, 1.0, u, ldu, u, ldu, -1.0, r.data, r.stride)
            } else {
                let l = aCopy
                let ldl = lda
                nFloat64.dgemm(Trans, NoTrans, n, n, n, 1.0, l, ldl, l, ldl, -1.0, r.data, r.stride)
            }
            let resid = dlange(MaxColumnSum, r.rows, r.cols, r.data, r.stride)
            if (resid > tol * Float64(n)) {
                throw Exception("${prefix}: unexpected result")
            }
        }
    }
    return true
}
