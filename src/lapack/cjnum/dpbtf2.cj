package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import std.math.*

extend LFloat64Implementation {
    public func dpbtf2(uplo: Uplo, n: Int, kd: Int, ab: Array<Float64>, ldab: Int): Bool {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (kd < 0) {
            throw Exception(kdLT0)
        } else if (ldab < kd + 1) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return true
        }

        if (ab.size < (n - 1) * ldab + kd + 1) {
            throw Exception(shortAB)
        }

        let nFloat64 = nFloat64Implementation()
        let kld = max(1, ldab - 1)

        if (uplo == Upper) {
            for (j in 0..n) {
                var ajj = ab[j * ldab]
                if (ajj <= 0.0) {
                    return false
                }
                ajj = sqrt(ajj)
                ab[j * ldab] = ajj

                let kn = min(kd, n - j - 1)
                if (kn > 0) {
                    nFloat64.dscal(kn, 1.0 / ajj, ab[j * ldab + 1..], 1)
                    nFloat64.dsyr(Upper, kn, -1.0, ab[j * ldab + 1..], 1, ab[(j + 1) * ldab..], kld)
                }
            }
            return true
        }

        for (j in 0..n) {
            var ajj = ab[j * ldab + kd]
            if (ajj <= 0.0) {
                return false
            }
            ajj = sqrt(ajj)
            ab[j * ldab + kd] = ajj
            let kn = min(kd, n - j - 1)
            if (kn > 0) {
                nFloat64.dscal(kn, 1.0 / ajj, ab[(j + 1) * ldab + kd - 1..], kld)
                nFloat64.dsyr(Lower, kn, -1.0, ab[(j + 1) * ldab + kd - 1..], kld, ab[(j + 1) * ldab + kd..], kld)
            }
        }
        return true
    }
}
