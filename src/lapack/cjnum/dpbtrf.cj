package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dpbtrf(uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64): Bool {
        const nbmax = 32

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

        var opts = Upper.toString()
        if (uplo == Lower) {
            opts = Lower.toString()
        }
        var nb = ilaenv(1, "DPBTRF", opts, n, kd, -1, -1)
        nb = min(nb, nbmax)

        if (nb <= 1 || kd < nb) {
            return dpbtf2(uplo, n, kd, ab, ldab)
        }

        let ldwork = nb
        let work = Array<Float64>(nb * ldwork, repeat: 0.0)
        let nFloat64 = nFloat64Implementation()

        if (uplo == Upper) {
            for (i in 0..n : nb) {
                let ib = min(nb, n - i)
                let ok = dpotf2(uplo, ib, ab[i * ldab..], ldab - 1)
                if (!ok) {
                    return false
                }
                if (i + ib >= n) {
                    continue
                }
                let i2 = min(kd - ib, n - i - ib)
                if (i2 > 0) {
                    nFloat64.dtrsm(Left, Upper, Trans, NonUnitDiag, ib, i2, 1.0, ab[i * ldab..], ldab - 1,
                        ab[i * ldab + ib..], ldab - 1)
                    nFloat64.dsyrk(Upper, Trans, i2, ib, -1.0, ab[i * ldab + ib..], ldab - 1, 1.0, ab[(i + ib) * ldab..],
                        ldab - 1)
                }
                let i3 = min(ib, n - i - kd)
                if (i3 > 0) {
                    for (ii in 0..ib) {
                        for (jj in 0..=min(ii, i3 - 1)) {
                            work[ii * ldwork + jj] = ab[(i + ii) * ldab + kd - ii + jj]
                        }
                    }
                    nFloat64.dtrsm(Left, Upper, Trans, NonUnitDiag, ib, i3, 1.0, ab[i * ldab..], ldab - 1, work, ldwork)
                    if (i2 > 0) {
                        nFloat64.dgemm(Trans, NoTrans, i2, i3, ib, -1.0, ab[i * ldab + ib..], ldab - 1, work, ldwork,
                            1.0, ab[(i + ib) * ldab + kd - ib..], ldab - 1)
                    }
                    nFloat64.dsyrk(Upper, Trans, i3, ib, -1.0, work, ldwork, 1.0, ab[(i + kd) * ldab..], ldab - 1)
                    for (ii in 0..ib) {
                        for (jj in 0..=min(ii, i3 - 1)) {
                            ab[(i + ii) * ldab + kd - ii + jj] = work[ii * ldwork + jj]
                        }
                    }
                }
            }
        } else {
            for (i in 0..n : nb) {
                let ib = min(nb, n - i)
                let ok = dpotf2(uplo, ib, ab[i * ldab + kd..], ldab - 1)
                if (!ok) {
                    return false
                }
                if (i + ib >= n) {
                    continue
                }
                let i2 = min(kd - ib, n - i - ib)
                if (i2 > 0) {
                    nFloat64.dtrsm(Right, Lower, Trans, NonUnitDiag, i2, ib, 1.0, ab[i * ldab + kd..], ldab - 1,
                        ab[(i + ib) * ldab + kd - ib..], ldab - 1)
                    nFloat64.dsyrk(Lower, NoTrans, i2, ib, -1.0, ab[(i + ib) * ldab + kd - ib..], ldab - 1, 1.0,
                        ab[(i + ib) * ldab + kd..], ldab - 1)
                }
                let i3 = min(ib, n - i - kd)
                if (i3 > 0) {
                    for (ii in 0..i3) {
                        for (jj in ii..ib) {
                            work[ii * ldwork + jj] = ab[(ii + i + kd) * ldab + jj - ii]
                        }
                    }
                    nFloat64.dtrsm(Right, Lower, Trans, NonUnitDiag, i3, ib, 1.0, ab[i * ldab + kd..], ldab - 1, work,
                        ldwork)
                    if (i2 > 0) {
                        nFloat64.dgemm(NoTrans, Trans, i3, i2, ib, -1.0, work, ldwork, ab[(i + ib) * ldab + kd - ib..],
                            ldab - 1, 1.0, ab[(i + kd) * ldab + ib..], ldab - 1)
                    }
                    nFloat64.dsyrk(Lower, NoTrans, i3, ib, -1.0, work, ldwork, 1.0, ab[(i + kd) * ldab + kd..], ldab - 1
                    )
                    for (ii in 0..i3) {
                        for (jj in ii..ib) {
                            ab[(ii + i + kd) * ldab + jj - ii] = work[ii * ldwork + jj]
                        }
                    }
                }
            }
        }
        return true
    }
}
