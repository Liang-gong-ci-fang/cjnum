package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlarft(direct: Direct, store: StoreV, n: Int64, k: Int64, v: Array<Float64>, ldv: Int64,
        tau: Array<Float64>, t: Array<Float64>, ldt: Int64): Unit {
        var (mv, nv) = (n, k)
        if (store == RowWise) {
            (mv, nv) = (k, n)
        }
        if (direct != Forward && direct != Backward) {
            throw Exception(badDirect)
        } else if (store != RowWise && store != ColumnWise) {
            throw Exception(badStoreV)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (k < 1) {
            throw Exception(kLT1)
        } else if (ldv < max(1, nv)) {
            throw Exception(badLdV)
        } else if (tau.size < k) {
            throw Exception(shortTau)
        } else if (ldt < max(1, k)) {
            throw Exception(shortT)
        }

        if (n == 0) {
            return
        }

        if (v.size < (mv - 1) * ldv + nv) {
            throw Exception(shortV)
        } else if (t.size < (k - 1) * ldt + k) {
            throw Exception(shortT)
        }

        let nFloat64 = nFloat64Implementation()

        if (direct == Forward) {
            var prevlastv = n - 1
            for (i in 0..k) {
                prevlastv = max(i, prevlastv)
                if (tau[i] == 0.0) {
                    for (j in 0..=i) {
                        t[j * ldt + i] = 0.0
                    }
                    continue
                }
                var lastv: Int64 = 0
                if (store == ColumnWise) {
                    lastv = n - 1
                    while (lastv >= i + 1) {
                        if (v[lastv * ldv + i] != 0.0) {
                            break
                        }
                        lastv--
                    }
                    for (j in 0..i) {
                        t[j * ldt + i] = -tau[i] * v[i * ldv + j]
                    }
                    let j = min(lastv, prevlastv)
                    nFloat64.dgemv(Trans, j - i, i, -tau[i], v[(i + 1) * ldv..], ldv, v[(i + 1) * ldv + i..], ldv, 1.0,
                        t[i..], ldt)
                } else {
                    lastv = n - 1
                    while (lastv >= i + 1) {
                        if (v[i * ldv + lastv] != 0.0) {
                            break
                        }
                        lastv--
                    }
                    for (j in 0..i) {
                        t[j * ldt + i] = -tau[i] * v[j * ldv + i]
                    }
                    let j = min(lastv, prevlastv)
                    nFloat64.dgemv(NoTrans, i, j - i, -tau[i], v[i + 1..], ldv, v[i * ldv + i + 1..], 1, 1.0, t[i..],
                        ldt)
                }
                nFloat64.dtrmv(Upper, NoTrans, NonUnitDiag, i, t, ldt, t[i..], ldt)
                t[i * ldt + i] = tau[i]
                if (i > 1) {
                    prevlastv = max(prevlastv, lastv)
                } else {
                    prevlastv = lastv
                }
            }
            return
        }
        var prevlastv = 0
        for (i in k - 1..=0 : -1) {
            if (tau[i] == 0.0) {
                for (j in i..k) {
                    t[j * ldt + i] = 0.0
                }
                continue
            }
            var lastv: Int64 = 0
            if (i < k - 1) {
                if (store == ColumnWise) {
                    lastv = 0
                    while (lastv < i) {
                        if (v[lastv * ldv + i] != 0.0) {
                            break
                        }
                        lastv++
                    }
                    for (j in i + 1..k) {
                        t[j * ldt + i] = -tau[i] * v[(n - k + i) * ldv + j]
                    }
                    let j = max(lastv, prevlastv)
                    nFloat64.dgemv(Trans, n - k + i - j, k - i - 1, -tau[i], v[j * ldv + i + 1..], ldv, v[j * ldv + i..],
                        ldv, 1.0, t[(i + 1) * ldt + i..], ldt)
                } else {
                    lastv = 0
                    while (lastv < i) {
                        if (v[i * ldv + lastv] != 0.0) {
                            break
                        }
                        lastv++
                    }
                    for (j in i + 1..k) {
                        t[j * ldt + i] = -tau[i] * v[j * ldv + n - k + i]
                    }
                    let j = max(lastv, prevlastv)
                    nFloat64.dgemv(NoTrans, k - i - 1, n - k + i - j, -tau[i], v[(i + 1) * ldv + j..], ldv,
                        v[i * ldv + j..], 1, 1.0, t[(i + 1) * ldt + i..], ldt)
                }
                nFloat64.dtrmv(Lower, NoTrans, NonUnitDiag, k - i - 1, t[(i + 1) * ldt + i + 1..], ldt,
                    t[(i + 1) * ldt + i..], ldt)
                if (i > 0) {
                    prevlastv = min(prevlastv, lastv)
                } else {
                    prevlastv = lastv
                }
            }
            t[i * ldt + i] = tau[i]
        }
    }
}
