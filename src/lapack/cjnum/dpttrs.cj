package cjnum.lapack.cjnum

extend LFloat64Implementation {
    // public func dpttrs(n: Int64, nrhs: Int64, d: Array<Float64>, e: Array<Float64>, b: Array<Float64>, ldb: Int64): Unit {
    public func dpttrs(_: Int64, _: Int64, _: Array<Float64>, _: Array<Float64>, _: Array<Float64>, _: Int64): Unit {
        // TODO
        throw Exception("not implemented")
        // if (n < 0) {
        //     throw Exception(nLT0)
        // }
        // if (nrhs < 0) {
        //     throw Exception(nrhsLT0)
        // }
        // if (ldb < max(1, nrhs)) {
        //     throw Exception(badLdB)
        // }

        // if (n == 0 || nrhs == 0) {
        //     return 
        // }

        // if (d.size < n) {
        //     throw Exception(shortD)
        // }
        // if (e.size < n - 1) {
        //     throw Exception(shortE)
        // }
        // if (b.size < (n - 1) * ldb + nrhs) {
        //     throw Exception(shortB)
        // }

        // let nb = if (nrhs > 1) {
        //     max(1, ilaenv(1, "DPTTRS", " ", n, nrhs, -1, -1))
        // } else {
        //     1
        // }

        // if (nb >= nrhs) {
        //     dptts2(n, nrhs, d, e, b, ldb)
        // } else {
        //     var j = 0
        //     while (j < nrhs) {
        //         let jb = min(nrhs - j, nb)
        //         dptts2(n, jb, d, e, b[j..], ldb)
        //         j += nb
        //     }
        // }
    }
}
