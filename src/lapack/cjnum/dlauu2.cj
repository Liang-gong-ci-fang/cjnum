// Copyright ©2018 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.blas.blas64.*
import std.math.*

extend LFloat64Implementation {
    public func dlauu2(uplo: Uplo, n: Int, a: Array<Float64>, lda: Int): Unit {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }

        let nFloat64 = nFloat64Implementation()

        if (uplo == Upper) {
            for (i in 0..n) {
                let aii = a[i * lda + i]
                if (i < n - 1) {
                    a[i * lda + i] = nFloat64.ddot(n - i, a[i * lda + i..], 1, a[i * lda + i..], 1)
                    nFloat64.dgemv(NoTrans, i, n - i - 1, 1.0, a[i + 1..], lda, a[i * lda + i + 1..], 1, aii, a[i..],
                        lda)
                } else {
                    nFloat64.dscal(i + 1, aii, a[i..], lda)
                }
            }
        } else {
            for (i in 0..n) {
                let aii = a[i * lda + i]
                if (i < n - 1) {
                    a[i * lda + i] = nFloat64.ddot(n - i, a[i * lda + i..], lda, a[i * lda + i..], lda)
                    nFloat64.dgemv(Trans, n - i - 1, i, 1.0, a[(i + 1) * lda..], lda, a[(i + 1) * lda + i..], lda, aii,
                        a[i * lda..], 1)
                } else {
                    nFloat64.dscal(i + 1, aii, a[i * lda..], 1)
                }
            }
        }
    }
}
