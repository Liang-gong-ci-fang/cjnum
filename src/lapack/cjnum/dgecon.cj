package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.blas.blas64.*

extend LFloat64Implementation <: LFloat64{
    public func dgecon(norm: MatrixNorm, n: Int64, a: Array<Float64>, lda: Int64, anorm: Float64, work: Array<Float64>, iwork: Array<Int64>): Float64 {
        if (norm != MaxColumnSum && norm != MaxRowSum) {
            throw Exception(badNorm)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }
        if (anorm < 0.0) {
            throw Exception(negANorm)
        }

        if (n == 0) {
            return 1.0
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (work.size < 4 * n) {
            throw Exception(shortWork)
        }
        if (iwork.size < n) {
            throw Exception(shortIWork)
        }

        if (anorm == 0.0) {
            return 0.0
        }
        if (anorm.isNaN()) {
            return anorm
        }
        if (anorm == Float64.Inf) {
            return 0.0
        }

        let nFloat64 = nFloat64Implementation()
        var rcond = 0.0
        var ainvnm = 0.0
        var kase = 0
        var normin = false
        let isave = Array<Int64>(3, repeat: 0)
        let onenrm = norm == MaxColumnSum
        let smlnum = dlamchS
        let kase1 = if (onenrm) {
            1
        } else {
            2
        }
        while (true) {
            (ainvnm, kase) = dlacn2(n, work[n..], work, iwork, ainvnm, kase, isave)
            if (kase == 0) {
                if (ainvnm != 0.0) {
                    rcond = (1.0 / ainvnm) / anorm
                }
                return rcond
            }
            var sl = 0.0
            var su = 0.0
            if (kase == kase1) {
                sl = dlatrs(Lower, NoTrans, UnitDiag, normin, n, a, lda, work, work[2 * n..])
                su = dlatrs(Upper, NoTrans, NonUnitDiag, normin, n, a, lda, work, work[3 * n..])
            } else {
                sl = dlatrs(Upper, Trans, NonUnitDiag, normin, n, a, lda, work, work[3 * n..])
                su = dlatrs(Lower, Trans, UnitDiag, normin, n, a, lda, work, work[2 * n..])
            }
            let scale = sl * su
            normin = true
            if (scale != 1.0) {
                let ix = nFloat64.idamax(n, work, 1)
                if (scale == 0.0 || scale < abs(work[ix]) * smlnum) {
                    return rcond
                }
                drscl(n, scale, work, 1)
            }
        }
        throw Exception("RuntimeError")
    }
}