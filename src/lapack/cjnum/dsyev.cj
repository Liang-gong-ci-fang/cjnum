package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dsyev(jobz: EVJob, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64, w: Array<Float64>,
        work: Array<Float64>, lwork: Int64): Bool {
        var ok = false
        if (jobz != EVNone && jobz != EVCompute) {
            throw Exception(badEVJob)
        } else if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (lda < max(1, n)) {
            throw Exception(badLdA)
        } else if (lwork < max(1, 3 * n - 1) && lwork != -1) {
            throw Exception(badLWork)
        } else if (work.size < max(1, lwork)) {
            throw Exception(shortWork)
        }

        if (n == 0) {
            return true
        }

        var opts: String
        if (uplo == Upper) {
            opts = "U"
        } else {
            opts = "L"
        }
        let nb = ilaenv(1, "DSYTRD", opts, n, -1, -1, -1)
        let lworkopt = max(1, (nb + 2) * n)
        if (lwork == -1) {
            work[0] = Float64(lworkopt)
            return false
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        } else if (w.size < n) {
            throw Exception(shortW)
        }

        if (n == 1) {
            w[0] = a[0]
            work[0] = 2.0
            if (jobz == EVCompute) {
                a[0] = 1.0
            }
            return true
        }

        let safmin = dlamchS
        let eps = dlamchP
        let smlnum = safmin / eps
        let bignum = 1.0 / smlnum
        let rmin = sqrt(smlnum)
        let rmax = sqrt(bignum)

        let anrm = dlansy(MaxAbs, uplo, n, a, lda, work)
        var scaled = false
        var sigma = 0.0
        if (anrm > 0.0 && anrm < rmin) {
            scaled = true
            sigma = rmin / anrm
        } else if (anrm > rmax) {
            scaled = true
            sigma = rmax / anrm
        }
        if (scaled) {
            var kind = LowerTri
            if (uplo == Upper) {
                kind = UpperTri
            }
            dlascl(kind, 0, 0, 1.0, sigma, n, n, a, lda)
        }
        let inde = 0
        let indtau = inde + n
        let indwork = indtau + n
        let llwork = lwork - indwork
        dsytrd(uplo, n, a, lda, w, work[inde..], work[indtau..], work[indwork..], llwork)

        if (jobz == EVNone) {
            ok = dsterf(n, w, work[inde..])
        } else {
            dorgtr(uplo, n, a, lda, work[indtau..], work[indwork..], llwork)
            ok = dsteqr(UInt8(jobz), n, w, work[inde..], a, lda, work[indtau..])
        }
        if (!ok) {
            return false
        }

        if (scaled) {
            let nFloat64 = nFloat64Implementation()
            nFloat64.dscal(n, 1.0 / sigma, w, 1)
        }
        work[0] = Float64(lworkopt)
        return true
    }
}
