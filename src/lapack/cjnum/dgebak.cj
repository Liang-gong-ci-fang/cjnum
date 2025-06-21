package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dgebak(job: BalanceJob, side: EVSide, n: Int64, ilo: Int64, ihi: Int64, scale: Array<Float64>, m: Int64,
        v: Array<Float64>, ldv: Int64): Unit {
        if (job != BalanceNone && job != Permute && job != Scale && job != PermuteScale) {
            throw Exception(badBalanceJob)
        } else if (side != EVLeft && side != EVRight) {
            throw Exception(badEVSide)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (ilo < 0 || max(0, n - 1) < ilo) {
            throw Exception(badIlo)
        } else if (ihi < min(ilo, n - 1) || n <= ihi) {
            throw Exception(badIhi)
        } else if (m < 0) {
            throw Exception(mLT0)
        } else if (ldv < max(1, m)) {
            throw Exception(badLdV)
        }

        if (n == 0 || m == 0) {
            return
        }

        if (scale.size < n) {
            throw Exception(shortScale)
        }
        if (v.size < (n - 1) * ldv + m) {
            throw Exception(shortV)
        }

        if (job == BalanceNone) {
            return
        }

        let nFloat64 = nFloat64Implementation()
        if (ilo != ihi && job != Permute) {
            if (side == EVRight) {
                for (i in ilo..=ihi) {
                    nFloat64.dscal(m, scale[i], v[i * ldv..], 1)
                }
            } else {
                for (i in ilo..=ihi) {
                    nFloat64.dscal(m, 1.0 / scale[i], v[i * ldv..], 1)
                }
            }
        }
        if (job == Scale) {
            return
        }
        for (i in ilo - 1..=0: -1) {
            let k = Int64(scale[i])
            if (k == i) {
                continue
            }
            nFloat64.dswap(m, v[i * ldv..], 1, v[k * ldv..], 1)
        }
        for (i in ihi + 1..n) {
            let k = Int64(scale[i])
            if (k == i) {
                continue
            }
            nFloat64.dswap(m, v[i * ldv..], 1, v[k * ldv..], 1)
        }
    }
}
