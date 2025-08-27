package cjnum.lapack.cjnum

import cjnum.blas.*

extend LFloat64Implementation {
    public func dlagtm(trans: Transpose, m: Int64, n: Int64, alpha: Float64, dl: Array<Float64>, d: Array<Float64>,
        du: Array<Float64>, b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit {
        var dl_: Array<Float64> = dl
        var du_: Array<Float64> = du

        if (trans != NoTrans && trans != Trans && trans != ConjTrans) {
            throw Exception(badTrans)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (ldb < max(1, n)) {
            throw Exception(badLdB)
        }
        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (m == 0 || n == 0) {
            return
        }

        if (dl_.size < m - 1) {
            throw Exception(shortDL)
        }
        if (d.size < m) {
            throw Exception(shortD)
        }
        if (du_.size < m - 1) {
            throw Exception(shortDU)
        }
        if (b.size < (m - 1) * ldb + n) {
            throw Exception(shortB)
        }
        if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }

        if (beta != 1.0) {
            if (beta == 0.0) {
                for (i in 0..m) {
                    let ci = c[i * ldc..i * ldc + n]
                    for (j in 0..ci.size) {
                        ci[j] = 0.0
                    }
                }
            } else {
                for (i in 0..m) {
                    let ci = c[i * ldc..i * ldc + n]
                    for (j in 0..ci.size) {
                        ci[j] *= beta
                    }
                }
            }
        }

        if (alpha == 0.0) {
            return
        }

        if (m == 1) {
            if (alpha == 1.0) {
                for (j in 0..n) {
                    c[j] += d[0] * b[j]
                }
            } else {
                for (j in 0..n) {
                    c[j] += alpha * d[0] * b[j]
                }
            }
            return
        }

        if (trans != NoTrans) {
            (dl_, du_) = (du_, dl_)
        }

        if (alpha == 1.0) {
            for (j in 0..n) {
                c[j] += d[0] * b[j] + du_[0] * b[ldb + j]
            }
            for (i in 1..m - 1) {
                for (j in 0..n) {
                    c[i * ldc + j] += dl_[i - 1] * b[(i - 1) * ldb + j] + d[i] * b[i * ldb + j] + du_[i] * b[(i + 1) *
                            ldb + j]
                }
            }
            for (j in 0..n) {
                c[(m - 1) * ldc + j] += dl_[m - 2] * b[(m - 2) * ldb + j] + d[m - 1] * b[(m - 1) * ldb + j]
            }
        } else {
            for (j in 0..n) {
                c[j] += alpha * (d[0] * b[j] + du_[0] * b[ldb + j])
            }
            for (i in 1..m - 1) {
                for (j in 0..n) {
                    c[i * ldc + j] += alpha * (dl_[i - 1] * b[(i - 1) * ldb + j] + d[i] * b[i * ldb + j] + du_[i] * b[(i +
                            1) * ldb + j])
                }
            }
            for (j in 0..n) {
                c[(m - 1) * ldc + j] += alpha * (dl_[m - 2] * b[(m - 2) * ldb + j] + d[m - 1] * b[(m - 1) * ldb + j])
            }
        }
    }
}
