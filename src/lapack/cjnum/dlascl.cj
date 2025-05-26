package cjnum.lapack.cjnum

import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlascl(kind: MatrixType, _: Int, _: Int, cfrom: Float64, cto: Float64, m: Int, n: Int,
        a: Array<Float64>, lda: Int): Unit {
        match {
            case kind == b'H' || kind == b'B' || kind == b'Q' || kind == b'Z' => //UInt8 var a = b'x' // a is 120 with type UInt8
                throw Exception("not implemented")
            case kind == General || kind == UpperTri || kind == LowerTri => if (lda < max(1, n)) {
                throw Exception(badLdA)
            }
            case _ => throw Exception(badMatrixType)
        }
        if (cfrom == 0.0) {
            throw Exception(zeroCFrom)
        } else if (cfrom.isNaN()) {
            throw Exception(nanCFrom)
        } else if (cto.isNaN()) {
            throw Exception(nanCTo)
        } else if (m < 0) {
            throw Exception(mLT0)
        } else if (n < 0) {
            throw Exception(nLT0)
        }

        if (n == 0 || m == 0) {
            return
        }

        if (kind == General || kind == UpperTri || kind == LowerTri) {
            if (a.size < (m - 1) * lda + n) {
                throw Exception(shortA)
            }
        }

        let smlnum = dlamchS
        let bignum = 1.0 / smlnum
        var cfromc = cfrom
        var ctoc = cto
        let cfrom1 = cfromc * smlnum
        while (true) {
            var done: Bool
            var mul: Float64
            var ctol: Float64
            if (cfrom1 == cfromc) {
                mul = ctoc / cfromc
                done = true
                ctol = ctoc
            } else {
                ctol = ctoc / bignum
                if (ctol == ctoc) {
                    mul = ctoc
                    done = true
                    cfromc = 1.0
                } else if (abs(cfrom1) > abs(ctoc) && ctoc != 0.0) {
                    mul = smlnum
                    done = false
                    cfromc = cfrom1
                } else if (abs(ctol) > abs(cfromc)) {
                    mul = bignum
                    done = false
                    ctoc = ctol
                } else {
                    mul = ctoc / cfromc
                    done = true
                }
            }
            if (kind == General) {
                for (i in 0..m) {
                    for (j in 0..n) {
                        a[i * lda + j] *= mul
                    }
                }
            } else if (kind == UpperTri) {
                for (i in 0..m) {
                    for (j in i..n) {
                        a[i * lda + j] *= mul
                    }
                }
            } else if (kind == LowerTri) {
                for (i in 0..m) {
                    for (j in 0..=min(i, n - 1)) {
                        a[i * lda + j] *= mul
                    }
                }
            }
            if (done) {
                break
            }
        }
    }
}
