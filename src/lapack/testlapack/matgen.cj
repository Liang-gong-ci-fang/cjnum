package cjnum.lapack.testlapack

import std.random.*
import std.math.*
import cjnum.blas.*
import cjnum.blas.blas64.*

// TODO

func dlarnv(dst: Array<Float64>, dist: Int64, rnd: Random): Unit {
    match (dist) {
        case 1 => 
            for (i in 0..dst.size) {
                dst[i] = rnd.nextFloat64()
            }
        case 2 =>
            for (i in 0..dst.size) {
                dst[i] = 2.0 * rnd.nextFloat64() - 1.0
            }
        case 3 =>
            for (i in 0..dst.size) {
                dst[i] = rnd.nextGaussianFloat64()
            }
        case _ => throw Exception("testlapack: invalid dist")
    }
}

func dlattr(imat: Int64, uplo: Uplo, trans: Transpose, n: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>,
    work: Array<Float64>, rnd: Random): Diag {
    checkMatrix(n, n, a, lda)
    if (b.size != n) {
        throw Exception("testlapack: bad length of b")
    }
    if (work.size < 3 * n) {
        throw Exception("testlapack: insufficient length of work")
    }
    if (uplo != Upper && uplo != Lower) {
        throw Exception("testlapack: bad uplo")
    }
    if (trans != NoTrans && trans != Trans) {
        throw Exception("testlapack: bad trans")
    }

    if (n == 0) {
        return NonUnitDiag
    }

    const tiny = safmin
    const huge = (1.0 - ulp) / tiny

    let nFloat64 = nFloat64Implementation()

    var diag: Diag = NonUnitDiag

    match (imat) {
        case 7 => 
            diag = UnitDiag
            if (uplo == Upper) {
                for (i in 0..n) {
                    a[i * lda + i] = Float64.NaN
                    for (j in i + 1..n) {
                        a[i * lda + j] = 0.0
                    }
                }
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    for (j in 0..i) {
                        a[i * lda + j] = 0.0
                    }
                    a[i * lda + i] = Float64.NaN
                }
            }
        case 11 =>
            diag = NonUnitDiag
            if (uplo == Upper) {
                for (i in 0..n - 1) {
                    dlarnv(a[i * lda + i..i * lda + n], 2, rnd)
                }
            } else if (uplo == Lower) {
                for (i in 1..n) {
                    dlarnv(a[i * lda..i * lda + i + 1], 2, rnd)
                }
            }
            for (i in 0..n) {
                a[i * lda + i] = copysign(2.0, a[i * lda + i])
            }
            dlarnv(b, 2, rnd)
            let imax = nFloat64.idamax(n, b, 1)
            let bscal = huge / max(1.0, b[imax])
            nFloat64.dscal(n, bscal, b, 1)
        case 12 =>
            diag = NonUnitDiag
            let tscal = 1.0 / max(1.0, Float64(n - 1))
            if (uplo == Upper) {
                for (i in 0..n) {
                    dlarnv(a[i * lda + i..i * lda + n], 2, rnd)
                    nFloat64.dscal(n - i - 1, tscal, a[i * lda + i + 1..], 1)
                    a[i * lda + i] = copysign(1.0, a[i * lda + i])
                }
                a[(n - 1) * lda + n - 1] *= tiny
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    dlarnv(a[i * lda..i * lda + i + 1], 2, rnd)
                    nFloat64.dscal(i, tscal, a[i * lda..], 1)
                    a[i * lda + i] = copysign(1.0, a[i * lda + i])
                }
                a[0] *= tiny
            }
            dlarnv(b, 2, rnd)
        case 13 =>
            diag = NonUnitDiag
            if (uplo == Upper) {
                for (i in 0..n) {
                    dlarnv(a[i * lda + i..i * lda + n], 3, rnd)
                    a[i * lda + i] = copysign(1.0, a[i * lda + i])
                }
                a[(n - 1) * lda + n - 1] *= tiny
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    dlarnv(a[i * lda..i * lda + i + 1], 3, rnd)
                    a[i * lda + i] = copysign(1.0, a[i * lda + i])
                }
                a[0] *= tiny
            }
            dlarnv(b, 2, rnd)
        case 14 =>
            diag == NonUnitDiag
            if (uplo == Upper) {
                for (i in 0..n) {
                    for (j in i + 1..n) {
                        a[i * lda + j] = 0.0
                    }
                    if (((n - 1 - i) & 0x2) == 0) {
                        a[i * lda + i] = tiny
                    } else {
                        a[i * lda + i] = 1.0
                    }
                }
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    for (i in 0..n) {
                        for (j in 0..i) {
                            a[i * lda + j] = 0.0
                        }
                        if ((i & 0x2) == 0) {
                            a[i * lda + i] = tiny
                        } else {
                            a[i * lda + i] = 1.0
                        }
                    }
                }
            }

            if (uplo == Upper) {
                b[0] = 0.0
                var i = n - 1
                while (i > 0) {
                    b[i] = 0.0
                    b[i - 1] = tiny
                    i -= 2
                }
            } else if (uplo == Lower) {
                var i = 0
                while (i < n - 1) {
                    b[i] = 0.0
                    b[i + 1] = tiny
                    i += 2
                }
                b[n - 1] = 0.0
            }
        case 15 =>
            diag = NonUnitDiag
            let texp = 1.0 / max(1.0, Float64(n - 1))
            let tscal = pow(tiny, texp)
            if (uplo == Upper) {
                for (i in 0..n) {
                    a[i * lda + i] = tscal
                    if (i < n - 1) {
                        a[i * lda + i + 1] = -1.0
                    }
                    for (j in i + 2..n) {
                        a[i * lda + j] = 0.0
                    }
                }
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    for (j in 0..i - 1) {
                        a[i * lda + j] = 0.0
                    }
                    if (i > 0) {
                        a[i * lda + i - 1] = -1.0
                    }
                    a[i * lda + i] = tscal
                }
            }
            dlarnv(b, 2, rnd)
        case 16 =>
            diag = NonUnitDiag
            if (uplo == Upper) {
                for (i in 0..n) {
                    dlarnv(a[i * lda + i..i * lda + n], 2, rnd)
                    a[i * lda + i] = copysign(2.0, a[i * lda + i])
                }
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    dlarnv(a[i * lda..i * lda + i + 1], 2, rnd)
                    a[i * lda + i] = copysign(2.0, a[i * lda + i])
                }
            }
            let iy = n / 2
            a[iy * lda + iy] = 0.0
            dlarnv(b, 2, rnd)
            nFloat64.dscal(n, 2.0, b, 1)
        case 17 =>
            diag = NonUnitDiag
            let tscal = (1.0 - ulp) / tiny
            var texp = 1.0
            if (uplo == Upper) {
                for (i in 0..n) {
                    for (j in i..n) {
                        a[i * lda + j] = 0.0
                    }
                }
                var j = n - 1
                while (j >= 1) {
                    a[j] = -tscal / Float64(n + 1)
                    a[j * lda + j] = 1.0
                    b[j] = texp * (1.0 - ulp)
                    a[j - 1] = -tscal / Float64(n + 1) / Float64(n + 2)
                    a[(j - 1) * lda + j - 1] = 1.0
                    b[j - 1] = texp * Float64(n * n + n - 1)
                    texp *= 2.0
                    j -= 2
                }
                b[0] = Float64(n + 1) / Float64(n + 2) * tscal
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    for (j in 0..=i) {
                        a[i * lda + j] = 0.0
                    }
                }
                var j = 0
                while (j < n - 1) {
                    a[(n - 1) * lda + j] = -tscal / Float64(n + 1)
                    a[j * lda + j] = 1.0
                    b[j] = texp * (1.0 - ulp)
                    a[(n - 1) * lda + j + 1] = -tscal / Float64(n + 1) / Float64(n + 2)
                    a[(j + 1) * lda + j + 1] = 1.0
                    b[j + 1] = texp * Float64(n * n + n - 1)
                    texp *= 2.0
                    j += 2
                }
                b[n - 1] = Float64(n + 1) / Float64(n + 2) * tscal
            }
        case 18 =>
            diag = UnitDiag
            if (uplo == Upper) {
                for (i in 0..n) {
                    a[i * lda + i] = Float64.NaN
                    dlarnv(a[i * lda + i + 1..i * lda + n], 2, rnd)
                }
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    dlarnv(a[i * lda..i * lda + i], 2, rnd)
                    a[i * lda + i] = Float64.NaN
                }
            }
            dlarnv(b, 2, rnd)
            let iy = nFloat64.idamax(n, b, 1)
            let bnorm = abs(b[iy])
            let bscal = huge / max(1.0, bnorm)
            nFloat64.dscal(n, bscal, b, 1)
        case 19 =>
            diag = NonUnitDiag
            let tleft = huge / max(1.0, Float64(n - 1))
            let tscal = huge * (Float64(n - 1) / max(1.0, Float64(n)))
            if (uplo == Upper) {
                for (i in 0..n) {
                    dlarnv(a[i * lda + i..i * lda + n], 2, rnd)
                    for (j in i..n) {
                        let aij = a[i * lda + j]
                        a[i * lda + j] = copysign(tleft, aij) + tscal * aij
                    }
                }
            } else if (uplo == Lower) {
                for (i in 0..n) {
                    dlarnv(a[i * lda..i * lda + i + 1], 2, rnd)
                    for (j in 0..=i) {
                        let aij = a[i * lda + j]
                        a[i * lda + j] = copysign(tleft, aij) + tscal * aij
                    }
                }
            }
            dlarnv(b, 2, rnd)
            nFloat64.dscal(n, 2.0, b, 1)
        case _ => throw Exception("testlapack: invalid or unimplemented imat")
    }

    if (trans == Trans) {
        if (uplo == Upper) {
            for (j in 0.. n / 2) {
                nFloat64.dswap(n - 2 * j - 1, a[j * lda+ j ..], 1, a[(j + 1) * lda + n - j - 1..], -lda)
            }
        } else if (uplo == Lower) {
            for (j in 0..n / 2) {
                nFloat64.dswap(n - 2 * j - 1, a[j * lda + j..], lda, a[(n - j - 1) * lda + j + 1..], -1)
            }
        }
    }

    return diag
}

func checkMatrix(m: Int64, n: Int64, a: Array<Float64>, lda: Int64) {
    if (m < 0) {
        throw Exception("testlapack: m < 0")
    }
    if (n < 0) {
        throw Exception("testlapack: n < 0")
    }
    if (lda < max(1, n)) {
        throw Exception("testlapack: lda < max(1, n)")
    }
    if (a.size < (m - 1) * lda + n) {
        throw Exception("testlapack: insufficient matrix slice length")
    }
}

func newRandomSymTridiag(n: Int64, rnd: Random): (Array<Float64>, Array<Float64>) {
    if (n == 0) {
        return ([], [])
    }
    if (n == 1) {
        let d = Array<Float64>(1, { _ => rnd.nextFloat64() })
        return (d, [])
    }

    let d = Array<Float64>(n, item: 0.0)
    dlarnv(d, 1, rnd)
    let e = Array<Float64>(n - 1, item: 0.0)
    dlarnv(e, 2, rnd)

    d[0] += abs(e[0])
    for (i in 1..n - 1) {
        d[i] += abs(e[i]) + abs(e[i - 1])
    }
    d[n - 1] += abs(e[n - 2])

    return (d, e)
}

// TODO