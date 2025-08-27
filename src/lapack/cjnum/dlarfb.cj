package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlarfb(side: Side, trans: Transpose, direct: Direct, store: StoreV, m: Int64, n: Int64, k: Int64,
        v: Array<Float64>, ldv: Int64, t: Array<Float64>, ldt: Int64, c: Array<Float64>, ldc: Int64,
        work: Array<Float64>, ldwork: Int64): Unit {
        var nv = m
        if (side == Right) {
            nv = n
        }
        if (side != Left && side != Right) {
            throw Exception(badSide)
        } else if (trans != Trans && trans != NoTrans) {
            throw Exception(badTrans)
        } else if (direct != Forward && direct != Backward) {
            throw Exception(badDirect)
        } else if (store != ColumnWise && store != RowWise) {
            throw Exception(badStoreV)
        } else if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (k < 0) {
            throw Exception(kLT0)
        } else if (store == ColumnWise && ldv < max(1, k)) {
            throw Exception(badLdV)
        } else if (store == RowWise && ldv < max(1, nv)) {
            throw Exception(badLdV)
        } else if (ldt < max(1, k)) {
            throw Exception(badLdT)
        } else if (ldc < max(1, n)) {
            throw Exception(badLdC)
        } else if (ldwork < max(1, k)) {
            throw Exception(badLdWork)
        }

        if (m == 0 || n == 0) {
            return
        }

        var nw = n
        if (side == Right) {
            nw = m
        }
        if (store == ColumnWise && v.size < (nv - 1) * ldv + k) {
            throw Exception(shortV)
        } else if (store == RowWise && v.size < (k - 1) * ldv + nv) {
            throw Exception(shortV)
        } else if (t.size < (k - 1) * ldt + k) {
            throw Exception(shortT)
        } else if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        } else if (work.size < (nw - 1) * ldwork + k) {
            throw Exception(shortWork)
        }

        let nFloat64 = nFloat64Implementation()

        var transt = Trans
        if (trans == Trans) {
            transt = NoTrans
        }
        if (store == ColumnWise) {
            if (direct == Forward) {
                if (side == Left) {
                    for (j in 0..k) {
                        nFloat64.dcopy(n, c[j * ldc..], 1, work[j..], ldwork)
                    }
                    nFloat64.dtrmm(Right, Lower, NoTrans, UnitDiag, n, k, 1.0, v, ldv, work, ldwork)
                    if (m > k) {
                        nFloat64.dgemm(Trans, NoTrans, n, k, m - k, 1.0, c[k * ldc..], ldc, v[k * ldv..], ldv, 1.0, work,
                            ldwork)
                    }
                    nFloat64.dtrmm(Right, Upper, transt, NonUnitDiag, n, k, 1.0, t, ldt, work, ldwork)
                    if (m > k) {
                        nFloat64.dgemm(NoTrans, Trans, m - k, n, k, -1.0, v[k * ldv..], ldv, work, ldwork, 1.0,
                            c[k * ldc..], ldc)
                    }
                    nFloat64.dtrmm(Right, Lower, Trans, UnitDiag, n, k, 1.0, v, ldv, work, ldwork)
                    for (i in 0..n) {
                        for (j in 0..k) {
                            c[j * ldc + i] -= work[i * ldwork + j]
                        }
                    }
                    return
                }
                for (i in 0..k) {
                    nFloat64.dcopy(m, c[i..], ldc, work[i..], ldwork)
                }
                nFloat64.dtrmm(Right, Lower, NoTrans, UnitDiag, m, k, 1.0, v, ldv, work, ldwork)
                if (n > k) {
                    nFloat64.dgemm(NoTrans, NoTrans, m, k, n - k, 1.0, c[k..], ldc, v[k * ldv..], ldv, 1.0, work, ldwork
                    )
                }
                nFloat64.dtrmm(Right, Upper, trans, NonUnitDiag, m, k, 1.0, t, ldt, work, ldwork)
                if (n > k) {
                    nFloat64.dgemm(NoTrans, Trans, m, n - k, k, -1.0, work, ldwork, v[k * ldv..], ldv, 1.0, c[k..], ldc)
                }
                nFloat64.dtrmm(Right, Lower, Trans, UnitDiag, m, k, 1.0, v, ldv, work, ldwork)
                for (i in 0..m) {
                    for (j in 0..k) {
                        c[i * ldc + j] -= work[i * ldwork + j]
                    }
                }
                return
            }
            if (side == Left) {
                for (j in 0..k) {
                    nFloat64.dcopy(n, c[(m - k + j) * ldc..], 1, work[j..], ldwork)
                }
                nFloat64.dtrmm(Right, Upper, NoTrans, UnitDiag, n, k, 1.0, v[(m - k) * ldv..], ldv, work, ldwork)
                if (m > k) {
                    nFloat64.dgemm(Trans, NoTrans, n, k, m - k, 1.0, c, ldc, v, ldv, 1.0, work, ldwork)
                }
                nFloat64.dtrmm(Right, Lower, transt, NonUnitDiag, n, k, 1.0, t, ldt, work, ldwork)
                if (m > k) {
                    nFloat64.dgemm(NoTrans, Trans, m - k, n, k, -1.0, v, ldv, work, ldwork, 1.0, c, ldc)
                }
                nFloat64.dtrmm(Right, Upper, Trans, UnitDiag, n, k, 1.0, v[(m - k) * ldv..], ldv, work, ldwork)
                for (i in 0..n) {
                    for (j in 0..k) {
                        c[(m - k + j) * ldc + i] -= work[i * ldwork + j]
                    }
                }
                return
            }
            for (j in 0..k) {
                nFloat64.dcopy(m, c[n - k + j..], ldc, work[j..], ldwork)
            }

            nFloat64.dtrmm(Right, Upper, NoTrans, UnitDiag, m, k, 1.0, v[(n - k) * ldv..], ldv, work, ldwork)
            if (n > k) {
                nFloat64.dgemm(NoTrans, NoTrans, m, k, n - k, 1.0, c, ldc, v, ldv, 1.0, work, ldwork)
            }
            nFloat64.dtrmm(Right, Lower, trans, NonUnitDiag, m, k, 1.0, t, ldt, work, ldwork)
            if (n > k) {
                nFloat64.dgemm(NoTrans, Trans, m, n - k, k, -1.0, work, ldwork, v, ldv, 1.0, c, ldc)
            }
            nFloat64.dtrmm(Right, Upper, Trans, UnitDiag, m, k, 1.0, v[(n - k) * ldv..], ldv, work, ldwork)
            for (i in 0..m) {
                for (j in 0..k) {
                    c[i * ldc + n - k + j] -= work[i * ldwork + j]
                }
            }
            return
        }
        if (direct == Forward) {
            if (side == Left) {
                for (j in 0..k) {
                    nFloat64.dcopy(n, c[j * ldc..], 1, work[j..], ldwork)
                }
                nFloat64.dtrmm(Right, Upper, Trans, UnitDiag, n, k, 1.0, v, ldv, work, ldwork)
                if (m > k) {
                    nFloat64.dgemm(Trans, Trans, n, k, m - k, 1.0, c[k * ldc..], ldc, v[k..], ldv, 1.0, work, ldwork)
                }
                nFloat64.dtrmm(Right, Upper, transt, NonUnitDiag, n, k, 1.0, t, ldt, work, ldwork)
                if (m > k) {
                    nFloat64.dgemm(Trans, Trans, m - k, n, k, -1.0, v[k..], ldv, work, ldwork, 1.0, c[k * ldc..], ldc)
                }
                nFloat64.dtrmm(Right, Upper, NoTrans, UnitDiag, n, k, 1.0, v, ldv, work, ldwork)
                for (i in 0..n) {
                    for (j in 0..k) {
                        c[j * ldc + i] -= work[i * ldwork + j]
                    }
                }
                return
            }
            for (j in 0..k) {
                nFloat64.dcopy(m, c[j..], ldc, work[j..], ldwork)
            }
            nFloat64.dtrmm(Right, Upper, Trans, UnitDiag, m, k, 1.0, v, ldv, work, ldwork)
            if (n > k) {
                nFloat64.dgemm(NoTrans, Trans, m, k, n - k, 1.0, c[k..], ldc, v[k..], ldv, 1.0, work, ldwork)
            }
            nFloat64.dtrmm(Right, Upper, trans, NonUnitDiag, m, k, 1.0, t, ldt, work, ldwork)
            if (n > k) {
                nFloat64.dgemm(NoTrans, NoTrans, m, n - k, k, -1.0, work, ldwork, v[k..], ldv, 1.0, c[k..], ldc)
            }
            nFloat64.dtrmm(Right, Upper, NoTrans, UnitDiag, m, k, 1.0, v, ldv, work, ldwork)
            for (i in 0..m) {
                for (j in 0..k) {
                    c[i * ldc + j] -= work[i * ldwork + j]
                }
            }
            return
        }
        if (side == Left) {
            for (j in 0..k) {
                nFloat64.dcopy(n, c[(m - k + j) * ldc..], 1, work[j..], ldwork)
            }
            nFloat64.dtrmm(Right, Lower, Trans, UnitDiag, n, k, 1.0, v[m - k..], ldv, work, ldwork)
            if (m > k) {
                nFloat64.dgemm(Trans, Trans, n, k, m - k, 1.0, c, ldc, v, ldv, 1.0, work, ldwork)
            }
            nFloat64.dtrmm(Right, Lower, transt, NonUnitDiag, n, k, 1.0, t, ldt, work, ldwork)
            if (m > k) {
                nFloat64.dgemm(Trans, Trans, m - k, n, k, -1.0, v, ldv, work, ldwork, 1.0, c, ldc)
            }
            nFloat64.dtrmm(Right, Lower, NoTrans, UnitDiag, n, k, 1.0, v[m - k..], ldv, work, ldwork)
            for (i in 0..n) {
                for (j in 0..k) {
                    c[(m - k + j) * ldc + i] -= work[i * ldwork + j]
                }
            }
            return
        }
        for (j in 0..k) {
            nFloat64.dcopy(m, c[n - k + j..], ldc, work[j..], ldwork)
        }
        nFloat64.dtrmm(Right, Lower, Trans, UnitDiag, m, k, 1.0, v[n - k..], ldv, work, ldwork)
        if (n > k) {
            nFloat64.dgemm(NoTrans, Trans, m, k, n - k, 1.0, c, ldc, v, ldv, 1.0, work, ldwork)
        }
        nFloat64.dtrmm(Right, Lower, trans, NonUnitDiag, m, k, 1.0, t, ldt, work, ldwork)
        if (n > k) {
            nFloat64.dgemm(NoTrans, NoTrans, m, n - k, k, -1.0, work, ldwork, v, ldv, 1.0, c, ldc)
        }
        nFloat64.dtrmm(Right, Lower, NoTrans, UnitDiag, m, k, 1.0, v[n - k..], ldv, work, ldwork)
        for (i in 0..m) {
            for (j in 0..k) {
                c[i * ldc + n - k + j] -= work[i * ldwork + j]
            }
        }
    }
}
