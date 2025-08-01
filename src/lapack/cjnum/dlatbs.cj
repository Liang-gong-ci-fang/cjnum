package cjnum.lapack.cjnum

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.util.*
import std.math.*

extend LFloat64Implementation {
    public func dlatbs(uplo: Uplo, trans: Transpose, diag: Diag, normin: Bool, n: Int64, kd: Int64, ab: Array<Float64>,
        ldab: Int64, x: Array<Float64>, cnorm: Array<Float64>): Float64 {
        var scale = 0.0

        let noTran = trans == NoTrans
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        } else if (!noTran && trans != Trans && trans != ConjTrans) {
            throw Exception(badTrans)
        } else if (diag != NonUnitDiag && diag != UnitDiag) {
            throw Exception(badDiag)
        } else if (n < 0) {
            throw Exception(nLT0)
        } else if (kd < 0) {
            throw Exception(kdLT0)
        } else if (ldab < kd + 1) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return 1.0
        }

        if (ab.size < (n - 1) * ldab + kd + 1) {
            throw Exception(shortAB)
        } else if (x.size < n) {
            throw Exception(shortX)
        } else if (cnorm.size < n) {
            throw Exception(shortCNorm)
        }

        let smlnum = dlamchS / dlamchP
        let bignum = 1.0 / smlnum

        let nFloat64 = nFloat64Implementation()
        let kld = max(1, ldab - 1)
        if (!normin) {
            if (uplo == Upper) {
                for (j in 0..n) {
                    let jlen = min(j, kd)
                    if (jlen > 0) {
                        cnorm[j] = nFloat64.dasum(jlen, ab[(j - jlen) * ldab + jlen..], kld)
                    } else {
                        cnorm[j] = 0.0
                    }
                }
            } else {
                for (j in 0..n) {
                    let jlen = min(n - j - 1, kd)
                    if (jlen > 0) {
                        cnorm[j] = nFloat64.dasum(jlen, ab[(j + 1) * ldab + kd - 1..], kld)
                    } else {
                        cnorm[j] = 0.0
                    }
                }
            }
        }

        var jFirst: Int64
        var jLast: Int64
        var jInc: Int64
        var maind: Int64

        if (noTran) {
            if (uplo == Upper) {
                jFirst = n - 1
                jLast = -1
                jInc = -1
                maind = 0
            } else {
                jFirst = 0
                jLast = n
                jInc = 1
                maind = kd
            }
        } else {
            if (uplo == Upper) {
                jFirst = 0
                jLast = n
                jInc = 1
                maind = 0
            } else {
                jFirst = n - 1
                jLast = -1
                jInc = -1
                maind = kd
            }
        }

        let tmax = cnorm[nFloat64.idamax(n, cnorm, 1)]
        var tscal = 1.0
        if (tmax > bignum) {
            tscal = 1.0 / (smlnum * tmax)
            nFloat64.dscal(n, tscal, cnorm, 1)
        }

        var xMax = abs(x[nFloat64.idamax(n, x, 1)])
        var xBnd = xMax
        var grow = 0.0
        var skipComputeGrow = true
        if (tscal != 1.0) {
            skipComputeGrow = false
        }
        if (skipComputeGrow) {
            if (noTran) {
                if (diag == NonUnitDiag) {
                    grow = 1.0 / nanMax(xBnd, smlnum)
                    xBnd = grow
                    var j = jFirst
                    while (j != jLast) {
                        if (grow <= smlnum) {
                            skipComputeGrow = false
                        }
                        if (!skipComputeGrow) {
                            break
                        }
                        let tjj = abs(ab[j * ldab + maind])
                        xBnd = min(xBnd, min(1.0, tjj) * grow)
                        if (tjj + cnorm[j] >= smlnum) {
                            grow *= tjj / (tjj + cnorm[j])
                        } else {
                            grow = 0.0
                        }
                        j += jInc
                    }
                    grow = xBnd
                } else {
                    var grow = min(1.0, 1.0 / nanMax(xBnd, smlnum))
                    var j = jFirst
                    while (j != jLast) {
                        if (grow <= smlnum) {
                            skipComputeGrow = false
                        }
                        if (!skipComputeGrow) {
                            break
                        }
                        grow /= 1.0 + cnorm[j]
                        j += jInc
                    }
                }
            } else {
                if (diag == NonUnitDiag) {
                    grow = 1.0 / nanMax(xBnd, smlnum)
                    xBnd = grow
                    var j = jFirst
                    while (j != jLast) {
                        if (grow <= smlnum) {
                            skipComputeGrow = false
                        }
                        if (!skipComputeGrow) {
                            break
                        }
                        let xj = 1.0 + cnorm[j]
                        grow = min(grow, xBnd / xj)
                        let tjj = abs(ab[j * ldab + maind])
                        if (xj > tjj) {
                            xBnd *= tjj / xj
                        }
                        j += jInc
                    }
                    grow = min(grow, xBnd)
                } else {
                    grow = min(1.0, 1.0 / nanMax(xBnd, smlnum))
                    var j = jFirst
                    while (j != jLast) {
                        if (grow <= smlnum) {
                            skipComputeGrow = false
                        }
                        if (!skipComputeGrow) {
                            break
                        }
                        grow /= 1.0 + cnorm[j]
                        j += jInc
                    }
                }
            }
        }

        if (grow * tscal > smlnum) {
            nFloat64.dtbsv(uplo, trans, diag, n, kd, ab, ldab, x, 1)
            if (tscal != 1.0) {
                nFloat64.dscal(n, 1.0 / tscal, cnorm, 1)
            }
            return 1.0
        }

        scale = 1.0
        if (xMax > bignum) {
            scale = bignum / xMax
            nFloat64.dscal(n, scale, x, 1)
            xMax = bignum
        }

        if (noTran) {
            var j = jFirst
            while (j != jLast) {
                var xj = abs(x[j])
                var tjjs = tscal
                if (diag == NonUnitDiag) {
                    tjjs *= ab[j * ldab + maind]
                }
                let tjj = abs(tjjs)
                if (tjj > smlnum) {
                    if (tjj < 1.0 && xj > tjj * bignum) {
                        let rec = 1.0 / xj
                        nFloat64.dscal(n, rec, x, 1)
                        scale *= rec
                        xMax *= rec
                    }
                    x[j] /= tjjs
                    xj = abs(x[j])
                } else if (tjj > 0.0) {
                    if (xj > tjj * bignum) {
                        var rec = tjj * bignum / xj
                        if (cnorm[j] > 1.0) {
                            rec /= cnorm[j]
                        }
                        nFloat64.dscal(n, rec, x, 1)
                        scale *= rec
                        xMax *= rec
                    }
                    x[j] /= tjjs
                    xj = abs(x[j])
                } else {
                    for (i in 0..n) {
                        x[i] = 0.0
                    }
                    x[j] = 1.0
                    xj = 1.0
                    scale = 0.0
                    xMax = 0.0
                }

                if (xj > 1.0) {
                    var rec = 1.0 / xj
                    if (cnorm[j] > (bignum - xMax) * rec) {
                        rec *= 0.5
                        nFloat64.dscal(n, rec, x, 1)
                        scale *= rec
                    }
                } else if (xj * cnorm[j] > bignum - xMax) {
                    nFloat64.dscal(n, 0.5, x, 1)
                    scale *= 0.5
                }

                if (uplo == Upper) {
                    if (j > 0) {
                        let jlen = min(j, kd)
                        if (jlen > 0) {
                            nFloat64.daxpy(jlen, -x[j] * tscal, ab[(j - jlen) * ldab + jlen..], kld, x[j - jlen..], 1)
                        }
                        let i = nFloat64.idamax(j, x, 1)
                        xMax = abs(x[i])
                    }
                } else if (j < n - 1) {
                    let jlen = min(kd, n - j - 1)
                    if (jlen > 0) {
                        nFloat64.daxpy(jlen, -x[j] * tscal, ab[(j + 1) * ldab + kd - 1..], kld, x[j + 1..], 1)
                    }
                    let i = j + 1 + nFloat64.idamax(n - j - 1, x[j + 1..], 1)
                    xMax = abs(x[i])
                }
                j += jInc
            }
        } else {
            var j = jFirst
            while (j != jLast) {
                var xj = abs(x[j])
                var tjjs = tscal
                if (diag == NonUnitDiag) {
                    tjjs *= ab[j * ldab + maind]
                }
                let tjj = abs(tjjs)
                var rec = 1.0 / nanMax(1.0, xMax)
                var uscal = tscal
                if (cnorm[j] > (bignum - xj) * rec) {
                    rec *= 0.5
                    if (tjj > 1.0) {
                        rec = min(1.0, rec * tjj)
                        uscal /= tjjs
                    }
                    if (rec < 1.0) {
                        nFloat64.dscal(n, rec, x, 1)
                        scale *= rec
                        xMax *= rec
                    }
                }

                var sumj = 0.0
                if (uscal == 1.0) {
                    if (uplo == Upper) {
                        let jlen = min(j, kd)
                        if (jlen > 0) {
                            sumj = nFloat64.ddot(jlen, ab[(j - jlen) * ldab + jlen..], kld, x[j - jlen..], 1)
                        }
                    } else {
                        let jlen = min(n - j - 1, kd)
                        if (jlen > 0) {
                            sumj = nFloat64.ddot(jlen, ab[(j + 1) * ldab + kd - 1..], kld, x[j + 1..], 1)
                        }
                    }
                } else {
                    if (uplo == Upper) {
                        let jlen = min(j, kd)
                        for (i in 0..jlen) {
                            sumj += (ab[(j - jlen + i) * ldab + jlen - i] * uscal) * x[j - jlen + i]
                        }
                    } else {
                        let jlen = min(n - j - 1, kd)
                        for (i in 0..jlen) {
                            sumj += (ab[(j + 1 + i) * ldab + kd - 1 - i] * uscal) * x[j + i + 1]
                        }
                    }
                }

                if (uscal == tscal) {
                    x[j] -= sumj
                    xj = abs(x[j])
                    if (tjj > smlnum) {
                        if (tjj < 1.0 && xj > tjj * bignum) {
                            rec = 1.0 / xj
                            nFloat64.dscal(n, rec, x, 1)
                            scale *= rec
                            xMax *= rec
                        }
                        x[j] /= tjjs
                    } else if (tjj > 0.0) {
                        if (xj > tjj * bignum) {
                            rec = (tjj * bignum) / xj
                            nFloat64.dscal(n, rec, x, 1)
                            scale *= rec
                            xMax *= rec
                        }
                        x[j] /= tjjs
                    } else {
                        for (i in 0..n) {
                            x[i] = 0.0
                        }
                        x[j] = 1.0
                        scale = 0.0
                        xMax = 0.0
                    }
                } else {
                    x[j] = x[j] / tjjs - sumj
                }
                xMax = nanMax(xMax, abs(x[j]))
                j += jInc
            }
            scale /= tscal
        }

        if (tscal != 1.0) {
            nFloat64.dscal(n, 1.0 / tscal, cnorm, 1)
        }
        return scale
    }
}
