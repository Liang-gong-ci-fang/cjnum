package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlatrs(uplo: Uplo, trans: Transpose, diag: Diag, normin: Bool, n: Int64, a: Array<Float64>, lda: Int64,
        x: Array<Float64>, cnorm: Array<Float64>): Float64 {
        if (uplo != Upper && uplo != Lower) {
            throw Exception(badUplo)
        }
        if (trans != NoTrans && trans != Trans && trans != ConjTrans) {
            throw Exception(badTrans)
        }
        if (diag != NonUnitDiag && diag != UnitDiag) {
            throw Exception(badDiag)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (n == 0) {
            return 1.0
        }

        if (a.size < (n - 1) * lda + n) {
            throw Exception(shortA)
        }
        if (x.size < n) {
            throw Exception(shortX)
        }
        if (cnorm.size < n) {
            throw Exception(shortCNorm)
        }

        let upper = uplo == Upper
        let nonUnitDiag = diag == NonUnitDiag
        let smlnum = dlamchS / dlamchP
        let bignum = 1.0 / smlnum
        let scale = 1.0

        let nFloat64 = nFloat64Implementation()

        if (!normin) {
            if (upper) {
                cnorm[0] = 0.0
                for (j in 1..n) {
                    cnorm[j] = nFloat64.dasum(j, a[j..], lda)
                }
            } else {
                for (j in 0..n - 1) {
                    cnorm[j] = nFloat64.dasum(n - j - 1, a[(j + 1) * lda + j..], lda)
                }
                cnorm[n - 1] = 0.0
            }
        }

        let imax = nFloat64.idamax(n, cnorm, 1)
        var tscal = 0.0
        if (cnorm[imax] <= bignum) {
            tscal = 1.0
        } else {
            var tmax = cnorm[imax]
            if (tmax <= Float64.Max) {
                tscal = 1.0 / (smlnum * tmax)
                nFloat64.dscal(n, tscal, cnorm, 1)
            } else {
                tmax = 0.0
                if (upper) {
                    for (j in 1..n) {
                        tmax = max(dlange(MaxAbs, j, 1, a[j..], lda, Array<Float64>()), tmax)
                    }
                } else {
                    for (j in 0..n - 1) {
                        tmax = max(dlange(MaxAbs, n - j - 1, 1, a[(j + 1) * lda + j..], lda, Array<Float64>()), tmax)
                    }
                }

                if (tmax <= Float64.Max) {
                    tscal = 1.0 / (smlnum * tmax)
                    for (j in 0..n) {
                        if (cnorm[j] <= Float64.Max) {
                            cnorm[j] *= tscal
                        } else {
                            cnorm[j] = 0.0
                            if (upper) {
                                for (i in 0..j) {
                                    cnorm[j] += tscal * abs(a[i * lda + j])
                                }
                            } else {
                                for (i in j + 1..n) {
                                    cnorm[j] += tscal * abs(a[i * lda + j])
                                }
                            }
                        }
                    }
                } else {
                    nFloat64.dtrsv(uplo, trans, diag, n, a, lda, x, 1)
                    return scale
                }
            }
        }

        let j = nFloat64.idamax(n, x, 1)
        let xmax = abs(x[j])
        var xband = xmax
        var grow = 0.0
        var jfirst = 0
        var jlast = 0
        var jinc = 0
        if (trans == NoTrans) {
            if (upper) {
                jfirst = n - 1
                jlast = -1
                jinc = -1
            } else {
                jfirst = 0
                jlast = n
                jinc = 1
            }
            if (tscal != 1.0) {
                grow = 0.0
                return __solve__(xmax, scale, nonUnitDiag, upper, grow, tscal, smlnum, jfirst, jlast, jinc, uplo, trans, diag,
                    n, a, lda, x, cnorm)
            }
            if (nonUnitDiag) {
                grow = 1.0 / max(xband, smlnum)
                xband = grow
                var jj = jfirst
                while (jj != jlast) {
                    if (grow <= smlnum) {
                        return __solve__(xmax, scale, nonUnitDiag, upper, grow, tscal, smlnum, jfirst, jlast, jinc, uplo, trans,
                            diag, n, a, lda, x, cnorm)
                    }
                    let tjj = abs(a[jj * lda + jj])
                    xband = min(xband, min(1.0, tjj) * grow)
                    if (tjj + cnorm[jj] >= smlnum) {
                        grow *= tjj / (tjj + cnorm[jj])
                    } else {
                        grow = 0.0
                    }
                    jj += jinc
                }
                grow = xband
            } else {
                grow = min(1.0, 1.0 / max(xband, smlnum))
                var jj = jfirst
                while (jj != jlast) {
                    if (grow <= smlnum) {
                        return __solve__(xmax, scale, nonUnitDiag, upper, grow, tscal, smlnum, jfirst, jlast, jinc, uplo, trans,
                            diag, n, a, lda, x, cnorm)
                    }
                    grow *= 1.0 / (1.0 + cnorm[jj])
                    jj += jinc
                }
            }
        } else {
            if (upper) {
                jfirst = 0
                jlast = n
                jinc = 1
            } else {
                jfirst = n - 1
                jlast = -1
                jinc = -1
            }
            if (tscal != 1.0) {
                grow = 0.0
                return __solve__(xmax, scale, nonUnitDiag, upper, grow, tscal, smlnum, jfirst, jlast, jinc, uplo, trans, diag,
                    n, a, lda, x, cnorm)
            }
            if (nonUnitDiag) {
                grow = 1.0 / max(xband, smlnum)
                xband = grow
                var jj = jfirst
                while (jj != jlast) {
                    if (grow <= smlnum) {
                        return __solve__(xmax, scale, nonUnitDiag, upper, grow, tscal, smlnum, jfirst, jlast, jinc, uplo, trans,
                            diag, n, a, lda, x, cnorm)
                    }
                    let xj = 1.0 + cnorm[jj]
                    grow = min(grow, xband / xj)
                    let tjj = abs(a[jj * lda + jj])
                    if (xj > tjj) {
                        xband *= tjj / xj
                    }
                    jj += jinc
                }
                grow = min(grow, xband)
            } else {
                grow = min(1.0, 1.0 / max(xband, smlnum))
                var jj = jfirst
                while (jj != jlast) {
                    if (grow <= smlnum) {
                        return __solve__(xmax, scale, nonUnitDiag, upper, grow, tscal, smlnum, jfirst, jlast, jinc, uplo, trans,
                            diag, n, a, lda, x, cnorm)
                    }
                    let xj = 1.0 + cnorm[jj]
                    grow /= xj
                    jj += jinc
                }
            }
        }

        return __solve__(xmax, scale, nonUnitDiag, upper, grow, tscal, smlnum, jfirst, jlast, jinc, uplo, trans, diag, n, a,
            lda, x, cnorm)
    }
}

private func __solve__(xmax_: Float64, scale_: Float64, nonUnitDiag: Bool, upper: Bool, gorw: Float64, tscal: Float64,
    smlnum: Float64, jfirst: Int64, jlast: Int64, jinc: Int64, uplo: Uplo, trans: Transpose, diag: Diag, n: Int64,
    a: Array<Float64>, lda: Int64, x: Array<Float64>, cnorm: Array<Float64>): Float64 {
    var xmax = xmax_
    var scale = scale_
    let bignum = 1.0 / smlnum
    let nFloat64 = nFloat64Implementation()

    if (gorw * tscal > smlnum) {
        nFloat64.dtrsv(uplo, trans, diag, n, a, lda, x, 1)
        if (tscal != 1.0) {
            nFloat64.dscal(n, 1.0 / tscal, cnorm, 1)
        }
        return scale
    }

    if (xmax > bignum) {
        scale = bignum / xmax
        nFloat64.dscal(n, scale, x, 1)
        xmax = bignum
    }
    if (trans == NoTrans) {
        var jj = jfirst
        while (jj != jlast) {
            var xj = abs(x[jj])
            var tjjs = 0.0
            if (nonUnitDiag) {
                tjjs = a[jj * lda + jj] * tscal
            } else {
                tjjs = tscal
                if (tscal == 1.0) {
                    if (xj > 1.0) {
                        var rec = 1.0 / xj
                        if (cnorm[jj] > (bignum - xmax) * rec) {
                            rec *= 0.5
                            nFloat64.dscal(n, 0.5, x, 1)
                            scale *= rec
                        }
                    } else if (xj * cnorm[jj] > bignum - xmax) {
                        nFloat64.dscal(n, 0.5, x, 1)
                        scale *= 0.5
                    }
                    if (upper) {
                        if (jj > 0) {
                            nFloat64.daxpy(jj, -x[jj] * tscal, a[jj..], lda, x, 1)
                            let i = nFloat64.idamax(jj, x, 1)
                            xmax = abs(x[i])
                        }
                    } else {
                        if (jj < n - 1) {
                            nFloat64.daxpy(n - jj - 1, -x[jj] * tscal, a[(jj + 1) * lda + jj..], lda, x[jj + 1..], 1)
                            let i = jj + nFloat64.idamax(n - jj - 1, x[jj + 1..], 1)
                            xmax = abs(x[i])
                        }
                    }
                    jj += jinc
                    continue
                }
            }
            let tjj = abs(tjjs)
            if (tjj > smlnum) {
                if (tjj < 1.0) {
                    if (xj > tjj * bignum) {
                        let rec = 1.0 / xj
                        nFloat64.dscal(n, rec, x, 1)
                        scale *= rec
                        xmax *= rec
                    }
                }
                x[jj] /= tjjs
                xj = abs(x[jj])
            } else if (tjj > 0.0) {
                if (xj > tjj * bignum) {
                    var rec = (tjj * bignum) / xj
                    if (cnorm[jj] > 1.0) {
                        rec /= cnorm[jj]
                    }
                    nFloat64.dscal(n, rec, x, 1)
                    scale *= rec
                    xmax *= rec
                }
                x[jj] /= tjjs
                xj = abs(x[jj])
            } else {
                for (i in 0..n) {
                    x[i] = 0.0
                }
                x[jj] = 1.0
                xj = 1.0
                scale = 0.0
                xmax = 0.0
            }
            if (xj > 1.0) {
                var rec = 1.0 / xj
                if (cnorm[jj] > (bignum - xmax) * rec) {
                    rec *= 0.5
                    nFloat64.dscal(n, 0.5, x, 1)
                    scale *= rec
                }
            } else if (xj * cnorm[jj] > bignum - xmax) {
                nFloat64.dscal(n, 0.5, x, 1)
                scale *= 0.5
            }
            if (upper) {
                if (jj > 0) {
                    nFloat64.daxpy(jj, -x[jj] * tscal, a[jj..], lda, x, 1)
                    let i = nFloat64.idamax(jj, x, 1)
                    xmax = abs(x[i])
                }
            } else {
                if (jj < n - 1) {
                    nFloat64.daxpy(n - jj - 1, -x[jj] * tscal, a[(jj + 1) * lda + jj..], lda, x[jj + 1..], 1)
                    let i = jj + nFloat64.idamax(n - jj - 1, x[jj + 1..], 1)
                    xmax = abs(x[i])
                }
            }
            jj += jinc
        }
    } else {
        var jj = jfirst
        while (jj != jlast) {
            let xj = abs(x[jj])
            var uscal = tscal
            var rec = 1.0 / max(xmax, 1.0)
            var tjjs = 0.0
            if (cnorm[jj] > (bignum - xj) * rec) {
                rec *= 0.5
                if (nonUnitDiag) {
                    tjjs = a[jj * lda + jj] * tscal
                } else {
                    tjjs = tscal
                }
                let tjj = abs(tjjs)
                if (tjj > 1.0) {
                    rec = min(1.0, rec * tjj)
                    uscal /= tjjs
                } 
                if (rec < 1.0) {
                    nFloat64.dscal(n, rec, x, 1)
                    scale *= rec
                    xmax *= rec
                }
            }

            var sumj = 0.0
            if (uscal == 1.0) {
                if (upper) {
                    sumj = nFloat64.ddot(jj, a[jj..], lda, x, 1)
                } else if (jj < n - 1) {
                    sumj = nFloat64.ddot(n - jj - 1, a[(jj + 1) * lda + jj..], lda, x[jj + 1..], 1)
                }
            } else {
                if (upper) {
                    for (i in 0..jj) {
                        sumj += (a[i * lda + jj] * uscal) * x[i]
                    }
                } else if (jj < n) {
                    for (i in jj + 1..n) {
                        sumj += (a[i * lda + jj] * uscal) * x[i]
                    }
                }
            }
            if (uscal == tscal) {
                x[jj] -= sumj
                let xj = abs(x[jj])
                var tjjs = 0.0
                if (nonUnitDiag) {
                    tjjs = a[jj * lda + jj] * tscal
                } else {
                    tjjs = tscal
                    if (tscal == 1.0) {
                        xmax = max(xmax, abs(x[jj]))
                        jj += jinc
                        continue
                    }
                }
                let tjj = abs(tjjs)
                if (tjj > smlnum) {
                    if (tjj < 1.0) {
                        if (xj > tjj * bignum) {
                            rec = 1.0 / xj
                            nFloat64.dscal(n ,rec, x, 1)
                            scale *= rec
                            xmax *= rec
                        }
                    }
                    x[jj] /= tjjs
                } else if (tjj > 0.0) {
                    if (xj > tjj * bignum) {
                        let rec = (tjj * bignum) / xj
                        nFloat64.dscal(n, rec, x, 1)
                        scale *= rec
                        xmax *= rec
                    }
                    x[jj] /= tjjs
                } else {
                    for (i in 0..n) {
                        x[i] = 0.0
                    }
                    x[jj] = 1.0
                    scale = 0.0
                    xmax = 0.0
                }
            } else {
                x[jj] = x[jj] / tjjs - sumj
            }
            xmax = max(xmax, abs(x[jj]))
            jj += jinc
        }
    }

    scale /= tscal
    if (tscal != 1.0) {
        nFloat64.dscal(n, 1.0 / tscal, cnorm, 1)
    }
    return scale
}

// private func __skip1__(xj: Float64, xmax_: Float64, scale_: Float64, jj: Int64, upper: Bool, tscal: Float64,
//     smlnum: Float64, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, cnorm: Array<Float64>): Float64 {
//     var xmax = xmax_
//     var scale = scale_
//     let bignum = 1.0 / smlnum
//     let nFloat64 = nFloat64Implementation()
//     if (xj > 1.0) {
//         var rec = 1.0 / xj
//         if (cnorm[jj] > (bignum - xmax) * rec) {
//             rec *= 0.5
//             nFloat64.dscal(n, 0.5, x, 1)
//             scale *= rec
//         }
//     } else if (xj * cnorm[jj] > bignum - xmax) {
//         nFloat64.dscal(n, 0.5, x, 1)
//         scale *= 0.5
//     }
//     if (upper) {
//         if (jj > 0) {
//             nFloat64.daxpy(jj, -x[jj] * tscal, a[jj..], lda, x, 1)
//             let i = nFloat64.idamax(jj, x, 1)
//             xmax = abs(x[i])
//         }
//     } else {
//         if (jj < n - 1) {
//             nFloat64.daxpy(n - jj - 1, -x[jj] * tscal, a[(jj + 1) * lda + jj..], lda, x[jj + 1..], 1)
//             let i = jj + nFloat64.idamax(n - jj - 1, x[jj + 1..], 1)
//             xmax = abs(x[i])
//         }
//     }

//     return scale
// }
