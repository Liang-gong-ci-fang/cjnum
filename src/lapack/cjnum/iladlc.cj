package cjnum.lapack.cjnum

extend LFloat64Implementation {
    public func iladlc(m: Int64, n: Int64, a: Array<Float64>, lda: Int64): Int64 {
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0 || m == 0) {
            return -1
        }

        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }

        if (a[n - 1] != 0.0 || a[(m - 1) * lda + (n - 1)] != 0.0) {
            return n - 1
        }

        var highest = -1
        for (i in 0..m) {
            var j = n - 1
            while (j >= 0) {
                if (a[i * lda + j] != 0.0) {
                    highest = max(highest, j)
                    break
                }
                j--
            }
        }
        return highest
    }
}