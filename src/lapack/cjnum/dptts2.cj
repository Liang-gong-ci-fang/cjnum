package cjnum.lapack.cjnum

import cjnum.blas.blas64.*

extend LFloat64Implementation {
    public func dptts2(n: Int64, nrhs: Int64, d: Array<Float64>, e: Array<Float64>, b: Array<Float64>, ldb: Int64): Unit {
        if (n <= 1) {
            if (n == 1) {
                let nFloat64 = nFloat64Implementation()
                nFloat64.dscal(nrhs, 1.0 / d[0], b, 1)
            }
            return 
        }

        for (j in 0..nrhs) {
            for (i in 1..n) {
                b[i * ldb + j] -= b[(i - 1) * ldb + j] * e[i - 1]
                // let y = b[(i - 1) * ldb + j] * e[i - 1] - com
                // let t = b[i * ldb + j] - y
                // com = (t - b[i * ldb + j]) + y
                // b[i * ldb + j] = t
            }
            b[(n - 1) * ldb + j] /= d[n - 1]
            var i = n - 2
            while (i >= 0) {
                b[i * ldb + j] = b[i * ldb + j] / d[i] - b[(i + 1) * ldb + j] * e[i]
                // let div = b[i * ldb + j] / d[i]
                // let sub = b[(i + 1) * ldb + j] * e[i]
                // // 使用 FMA（Fused Multiply-Add）减少误差
                // b[i * ldb + j] = div - sub
                i--
            }
        }
    }
}