package cjnum.lapack.cjnum

import std.math.*

extend LFloat64Implementation {
    public func dgtsv(n: Int64, nrhs: Int64, dl: Array<Float64>, d: Array<Float64>, du: Array<Float64>, b: Array<Float64>,
        ldb: Int64): Bool {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (nrhs < 0) {
            throw Exception(nrhsLT0)
        } else if (ldb < max(1, nrhs)) {
            throw Exception(badLdB)
        }

        if (n == 0 || nrhs == 0) {
            return true
        }

        if (dl.size < n - 1) {
            throw Exception(shortDL)
        } else if (d.size < n) {
            throw Exception(shortD)
        } else if (du.size < n - 1) {
            throw Exception(shortDU)
        } else if (b.size < (n - 1) * ldb + nrhs) {
            throw Exception(shortB)
        }

        let dl_ = dl[0..n - 1]
        let d_ = d[0..n]
        let du_ = du[0..n - 1]

        for (i in 0..n - 1) {
            if (abs(d_[i]) >= abs(dl_[i])) {
                if (d_[i] == 0.0) {
                    return false
                }
                let fact = dl_[i] / d_[i]
                d_[i + 1] = d_[i + 1] - fact * du_[i]
                for (j in 0..nrhs) {
                    b[(i + 1) * ldb + j] = b[(i + 1) * ldb + j] - fact * b[i * ldb + j]
                }
                dl_[i] = 0.0
            } else {
                let fact = d_[i] / dl_[i]
                d_[i] = dl_[i]
                let tmp = d_[i + 1]
                d_[i + 1] = du_[i] - fact * tmp
                du_[i] = tmp
                if (i + 1 < n - 1) {
                    dl_[i] = du_[i + 1]
                    du_[i + 1] = -fact * dl_[i]
                }
                for (j in 0..nrhs) {
                    let tmp = b[i * ldb + j]
                    b[i * ldb + j] = b[(i + 1) * ldb + j]
                    b[(i + 1) * ldb + j] = tmp - fact * b[(i + 1) * ldb + j]
                }
            }
        }

        if (d_[n - 1] == 0.0) {
            return false
        }

        for (j in 0..nrhs) {
            b[(n - 1) * ldb + j] = b[(n - 1) * ldb + j] / d_[n - 1]
            if (n > 1) {
                b[(n - 2) * ldb + j] = (b[(n - 2) * ldb + j] - du_[n - 2] * b[(n - 1) * ldb + j]) / d_[n - 2]
            }
            for (i in n - 3..=0 : -1) {
                b[i * ldb + j] = (b[i * ldb + j] - du_[i] * b[(i + 1) * ldb + j] - dl_[i] * b[(i + 2) * ldb + j]) / d_[i]
            }
        }

        return true
    }
}
