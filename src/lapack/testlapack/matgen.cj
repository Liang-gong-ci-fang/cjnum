package cjnum.lapack.testlapack

import std.random.*
import std.math.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64
import cjnum.floats.*
import cjnum.util.*

// TODO

// func dlagge(m: Int64, n: Int64, kl: Int64, ku: Int64, d: Array<Float64>, a: Array<Float64>, lda: Int64, rnd: Random,
//     work: Array<Float64>): Unit {
//     checkMatrix(m, n, a, lda)
//     if (kl < 0 || max(0, m - 1) < kl) {
//         throw Exception("testlapack: invalid value of kl")
//     }
//     if (ku < 0 || max(0, n - 1) < ku) {
//         throw Exception("testlapack: invalid value of ku")
//     }
//     if (d.size != min(m, n)) {
//         throw Exception("testlapack: bad length of d")
//     }
//     if (work.size < m + n) {
//         throw Exception("testlapack: insufficient work length")
//     }

//     for (i in 0..m) {
//         for (j in 0..n) {
//             a[i * lda + j] = 0.0
//         }
//     }
//     for (i in 0..min(m, n)) {
//         a[i * lda + i] = d[i]
//     }

//     let nFloat64 = nFloat64Implementation()

//     for (i in min(m, n) - 1..=0 : -1) {
//         if (i < m - 1) {
//             for (j in 0..m - i) {
//                 work[j] = rnd.nextGaussianFloat64()
//             }
//             let wn = nFloat64.dnrm2(m - i, work[0..m - i], 1)
//             let wa = copysign(wn, work[0])
//             var tau: Float64 = 0.0
//             if (wn != 0.0) {
//                 let wb = work[0] + wa
//                 nFloat64.dscal(m - i - 1, 1.0 / wb, work[1..m - i], 1)
//                 work[0] = 1.0
//                 tau = wb / wa
//             }

//             nFloat64.dgemv(Trans, m - i, n - i, 1.0, a[i * lda + i..], lda, work[0..m - i], 1, 0.0, work[m..m + n - i],
//                 1)
//             nFloat64.dger(m - i, n - i, -tau, work[0..m - i], 1, work[m..m + n - i], 1, a[i * lda + i..], lda)
//         }
//         if (i < n - 1) {
//             for (j in 0..n - i) {
//                 work[j] = rnd.nextGaussianFloat64()
//             }
//             let wn = nFloat64.dnrm2(n - i, work[0..n - i], 1)
//             let wa = copysign(wn, work[0])
//             var tau: Float64 = 0.0
//             if (wn != 0.0) {
//                 let wb = work[0] + wa
//                 nFloat64.dscal(n - i - 1, 1.0 / wb, work[1..n - i], 1)
//                 work[0] = 1.0
//                 tau = wb / wa
//             }

//             nFloat64.dgemv(NoTrans, m - i, n - i, 1.0, a[i * lda + i..], lda, work[0..n - i], 1, 0.0, work[n..n + m - i],
//                 1)
//             nFloat64.dger(m - i, n - i, -tau, work[n..n + m - i], 1, work[0..n - i], 1, a[i * lda + i..], lda)
//         }
//     }
// }

func dlarnv(dst: Array<Float64>, dist: Int64, rnd: Random): Unit {
    match (dist) {
        case 1 => for (i in 0..dst.size) {
            dst[i] = rnd.nextFloat64()
        }
        case 2 => for (i in 0..dst.size) {
            dst[i] = 2.0 * rnd.nextFloat64() - 1.0
        }
        case 3 => for (i in 0..dst.size) {
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
            diag = NonUnitDiag
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
            for (j in 0..n / 2) {
                nFloat64.dswap(n - 2 * j - 1, a[j * lda + j..], 1, a[(j + 1) * lda + n - j - 1..], -lda)
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

func randomOrthogonal(n: Int64, rnd: Random): General64 {
    let q = eye(n, n)
    let x = Array<Float64>(n, item: 0.0)
    let v = Array<Float64>(n, item: 0.0)
    for (j in 0..n - 1) {
        for (i in 0..j) {
            x[i] = 0.0
        }
        for (i in j..n) {
            x[i] = rnd.nextGaussianFloat64()
        }
        reflector(v, x, j)
        applyReflector(q, q, v)
    }
    return q
}

func reflector(v: Array<Float64>, col: Array<Float64>, j: Int64): Unit {
    let n = col.size
    if (v.size != n) {
        throw Exception("slice length mismatch")
    }
    if (j < 0 || n <= j) {
        throw Exception("invalid column index")
    }

    for (i in 0..v.size) {
        v[i] = 0.0
    }
    if (j == n - 1) {
        return
    }
    var s = fNorm(col[j..], 2.0)
    if (s == 0.0) {
        return
    }
    v[j] = col[j] + copysign(s, col[j])
    copy(v[j + 1..], col[j + 1..])
    s = fNorm(v[j..], 2.0)
    fScale(1.0 / s, v[j..])
}

func applyReflector(qh: General64, q: General64, v: Array<Float64>): Unit {
    let n = v.size
    if (qh.rows != n || qh.cols != n) {
        throw Exception("bad size of qh")
    }
    if (q.rows != n || q.cols != n) {
        throw Exception("bad size of q")
    }
    let qv = Array<Float64>(n, item: 0.0)
    gemv(NoTrans, 1.0, q, Vector(data: v, inc: 1), 0.0, Vector(data: qv, inc: 1))
    for (i in 0..n) {
        for (j in 0..n) {
            qh.data[i * qh.stride + j] = q.data[i * q.stride + j]
        }
    }
    for (i in 0..n) {
        for (j in 0..n) {
            qh.data[i * qh.stride + j] -= 2.0 * qv[i] * v[j]
        }
    }
    var norm2: Float64 = 0.0
    for (vi in v) {
        norm2 += vi * vi
    }
    let norm2inv = 1.0 / norm2
    for (i in 0..n) {
        for (j in 0..n) {
            qh.data[i * qh.stride + j] *= norm2inv
        }
    }
}

func newRandomSymTridiag(n: Int64, rnd: Random): (Array<Float64>, Array<Float64>) {
    if (n == 0) {
        return ([], [])
    }
    if (n == 1) {
        let d = Array<Float64>(1, {_ => rnd.nextFloat64()})
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

func dlattb(kind: Int64, uplo: Uplo, trans: Transpose, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64, rnd: Random): (Diag, 
    Array<Float64>) {
    var diag: Diag = 0
    var b: Array<Float64> = []
    if (kind < 1 || 18 < kind) {
        throw Exception("bad matrix kind")
    } else if ((6 <= kind && kind <= 9) || kind == 17) {
        diag = UnitDiag
    } else {
        diag = NonUnitDiag
    }

    if (n == 0) {
        return (diag, b)
    }

    const tiny = safmin
    const huge = (1.0 - ulp) / tiny

    const small = 0.25 * (safmin / ulp)
    const large = 1.0 / small
    const badc2 = 0.1 / ulp

    let badc1 = sqrt(badc2)

    var cndnum: Float64
    if (kind == 2 || kind == 8) {
        cndnum = badc1
    } else if (kind == 3 || kind == 9) {
        cndnum = badc2
    } else {
        cndnum = 2.0
    }

    let uniformM11 = {=> 2.0 * rnd.nextFloat64() - 1.0}

    b = Array<Float64>(n, item: 0.0)
    for (i in 0..b.size) {
        b[i] = uniformM11()
    }

    let nFloat64 = nFloat64Implementation()
    if (kind == 1 || kind == 2 || kind == 3 || kind == 4 || kind == 5) {
        var kl: Int64 = 0
        var ku: Int64 = 0
        if (uplo == Upper) {
            ku = kd
            kl = 0
        } else if (uplo == Lower) {
            ku = 0
            kl = kd
        }
        var anorm = 1.0
        if (kind == 4) {
            anorm = small
        } else if (kind == 5) {
            anorm = large
        }
        (_, _, _) = (kl, ku, anorm)
        throw Exception("test matrix type not implemented")
    } else if (kind == 6) {
        if (uplo == Upper) {
            for (i in 0..n) {
                ab[i * ldab] = Float64(i + 2)
                for (j in 1..min(n - i, kd + 1)) {
                    ab[i * ldab + j] = 0.0
                }
            }
        } else {
            for (i in 0..n) {
                for (j in max(0, kd - i)..kd) {
                    ab[i * ldab + j] = 0.0
                }
                ab[i * ldab + kd] = Float64(i + 2)
            }
        }
    } else if (kind == 7 || kind == 8 || kind == 9) {
        let tnorm = sqrt(cndnum)

        if (uplo == Upper) {
            for (i in 0..n) {
                ab[i * ldab] = Float64(i + 2)
                for (j in 1..min(n - i, kd + 1)) {
                    ab[i * ldab + j] = 0.0
                }
            }
        } else {
            for (i in 0..n) {
                for (j in max(0, kd - i)..kd) {
                    ab[i * ldab + j] = 0.0
                }
                ab[i * ldab + kd] = Float64(i + 2)
            }
        }

        if (kd == 0) {
        // Do nothing
        } else if (kd == 1) {
            if (n > 1) {
                if (uplo == Upper) {
                    ab[1] = copysign(tnorm, uniformM11())
                    for (i in 2..n - 1 : 2) {
                        ab[i * ldab + 1] = tnorm * uniformM11()
                    }
                } else {
                    ab[ldab] = copysign(tnorm, uniformM11())
                    for (i in 3..n : 2) {
                        ab[i * ldab] = tnorm * uniformM11()
                    }
                }
            }
        } else {
            let work1 = Array<Float64>(n, item: 0.0)
            let work2 = Array<Float64>(n, item: 0.0)
            var star1 = copysign(tnorm, uniformM11())
            let sfac = sqrt(tnorm)
            var plus1 = copysign(sfac, uniformM11())
            for (i in 0..n : 2) {
                work1[i] = plus1
                work2[i] = star1
                if (i + 1 == n) {
                    continue
                }
                let plus2 = star1 / plus1
                work1[i + 1] = plus2
                plus1 = star1 / plus2
                let rexp = uniformM11()
                if (rexp < 0.0) {
                    star1 = -pow(sfac, 1.0 - rexp)
                } else {
                    star1 = pow(sfac, 1.0 + rexp)
                }
            }
            if (uplo == Upper) {
                nFloat64.dcopy(n - 1, work1, 1, ab[1..], ldab)
                if (n > 2) {
                    nFloat64.dcopy(n - 2, work2, 1, ab[2..], ldab)
                }
            } else {
                nFloat64.dcopy(n - 1, work1, 1, ab[ldab + kd - 1..], ldab)
                if (n > 2) {
                    nFloat64.dcopy(n - 2, work2, 1, ab[2 * ldab + kd - 2..], ldab)
                }
            }
        }
    } else if (kind == 10) {
        if (uplo == Upper) {
            for (i in 0..n) {
                var j = 0
                while (j < min(n - j, kd + 1)) {
                    ab[i * ldab + j] = uniformM11()
                    j++
                }
                ab[i * ldab] = copysign(2.0, ab[i * ldab])
            }
        } else {
            for (i in 0..n) {
                for (j in max(0, kd - i)..kd + 1) {
                    ab[i * ldab + j] = uniformM11()
                }
                ab[i * ldab + kd] = copysign(2.0, ab[i * ldab + kd])
            }
        }
        let bnorm = abs(b[nFloat64.idamax(n, b, 1)])
        let bscal = huge / max(1.0, bnorm)
        nFloat64.dscal(n, bscal, b, 1)
    } else if (kind == 11) {
        let tscal = 1.0 / Float64(kd + 1)
        if (uplo == Upper) {
            for (i in 0..n) {
                let jlen = min(n - i, kd + 1)
                let arow = ab[i * ldab..i * ldab + jlen]
                dlarnv(arow, 2, rnd)
                if (jlen > 1) {
                    nFloat64.dscal(jlen - 1, tscal, arow[1..], 1)
                }
                ab[i * ldab] = copysign(1.0, ab[i * ldab])
            }
            ab[(n - 1) * ldab] *= tiny
        } else {
            for (i in 0..n) {
                let jlen = min(i + 1, kd + 1)
                let arow = ab[i * ldab + kd + 1 - jlen..i * ldab + kd + 1]
                dlarnv(arow, 2, rnd)
                if (jlen > 1) {
                    nFloat64.dscal(jlen - 1, tscal, arow[..jlen - 1], 1)
                }
                ab[i * ldab + kd] = copysign(1.0, ab[i * ldab + kd])
            }
            ab[kd] *= tiny
        }
    } else if (kind == 12) {
        if (uplo == Upper) {
            for (i in 0..n) {
                let jlen = min(n - i, kd + 1)
                let arow = ab[i * ldab..i * ldab + jlen]
                dlarnv(arow, 2, rnd)
                ab[i * ldab] = copysign(1.0, ab[i * ldab])
            }
            ab[(n - 1) * ldab] *= tiny
        } else {
            for (i in 0..n) {
                let jlen = min(i + 1, kd + 1)
                let arow = ab[i * ldab + kd + 1 - jlen..i * ldab + kd + 1]
                dlarnv(arow, 2, rnd)
                ab[i * ldab + kd] = copysign(1.0, ab[i * ldab + kd])
            }
            ab[kd] *= tiny
        }
    } else if (kind == 13) {
        if (uplo == Upper) {
            var icount = 1
            for (i in n - 1..=0) {
                if (icount <= 2) {
                    ab[i * ldab] = tiny
                } else {
                    ab[i * ldab] = 1.0
                }
                for (j in 1..min(n - i, kd + 1)) {
                    ab[i * ldab + j] = 0.0
                }
                icount++
                if (icount > 4) {
                    icount = 1
                }
            }
        } else {
            var icount = 1
            for (i in 0..n) {
                for (j in max(0, kd - i)..kd) {
                    ab[i * ldab + j] = 0.0
                }
                if (icount <= 2) {
                    ab[i * ldab + kd] = tiny
                } else {
                    ab[i * ldab + kd] = 1.0
                }
                icount++
                if (icount > 4) {
                    icount = 1
                }
            }
        }
        if (uplo == Upper) {
            b[0] = 0.0
            for (i in n - 1..1 : -2) {
                b[i] = 0.0
                b[i - 1] = tiny
            }
        } else {
            b[n - 1] = 0.0
            for (i in 0..n - 1 : 2) {
                b[i] = 0.0
                b[i + 1] = tiny
            }
        }
    } else if (kind == 14) {
        let tscal = pow(tiny, 1.0 / Float64(kd + 1))
        if (uplo == Upper) {
            for (i in 0..n) {
                ab[i * ldab] = tscal
                if (i < n - 1 && kd > 0) {
                    ab[i * ldab + 1] = -1.0
                }
                for (j in 2..min(n - i, kd + 1)) {
                    ab[i * ldab + j] = 0.0
                }
            }
            b[n - 1] = 1.0
        } else {
            for (i in 0..n) {
                for (j in max(0, kd - i)..kd - 1) {
                    ab[i * ldab + j] = 0.0
                }
                if (i > 0 && kd > 0) {
                    ab[i * ldab + kd - 1] = -1.0
                }
                ab[i * ldab + kd] = tscal
            }
            b[0] = 1.0
        }
    } else if (kind == 15) {
        let iy = n / 2
        if (uplo == Upper) {
            for (i in 0..n) {
                let jlen = min(n - i, kd + 1)
                dlarnv(ab[i * ldab..i * ldab + jlen], 2, rnd)
                if (i != iy) {
                    ab[i * ldab] = copysign(2.0, ab[i * ldab])
                } else {
                    ab[i * ldab] = 0.0
                }
            }
        } else {
            for (i in 0..n) {
                let jlen = min(i + 1, kd + 1)
                dlarnv(ab[i * ldab + kd + 1 - jlen..i * ldab + kd + 1], 2, rnd)
                if (i != iy) {
                    ab[i * ldab + kd] = copysign(2.0, ab[i * ldab + kd])
                } else {
                    ab[i * ldab + kd] = 0.0
                }
            }
        }
        nFloat64.dscal(n, 2.0, b, 1)
    } else if (kind == 17) {
        if (uplo == Upper) {
            for (i in 0..n) {
                ab[i * ldab] = Float64(i + 2)
                let jlen = min(n - i - 1, kd)
                if (jlen > 0) {
                    dlarnv(ab[i * ldab + 1..i * ldab + 1 + jlen], 2, rnd)
                }
            }
        } else {
            for (i in 0..n) {
                let jlen = min(i, kd)
                if (jlen > 0) {
                    dlarnv(ab[i * ldab + kd - jlen..i * ldab + kd], 2, rnd)
                }
                ab[i * ldab + kd] = Float64(i + 2)
            }
        }
        let bnorm = abs(b[nFloat64.idamax(n, b, 1)])
        let bscal = huge / max(1.0, bnorm)
        nFloat64.dscal(n, bscal, b, 1)
    } else if (kind == 18) {
        let tleft = huge / max(1.0, Float64(kd))
        let tscal = huge / max(1.0, Float64(kd))
        if (uplo == Upper) {
            for (i in 0..n) {
                for (j in 0..min(n - i, kd + 1)) {
                    let r = uniformM11()
                    ab[i * ldab + j] = copysign(tleft, r) + tscal * r
                }
            }
        } else {
            for (i in 0..n) {
                for (j in max(0, kd - i)..kd + 1) {
                    let r = uniformM11()
                    ab[i * ldab + j] = copysign(tleft, r) + tscal * r
                }
            }
        }
        nFloat64.dscal(n, 2.0, b, 1)
    }

    if (trans != NoTrans) {
        if (uplo == Upper) {
            for (j in 0..n / 2) {
                let jlen = min(n - 2 * j - 1, kd + 1)
                nFloat64.dswap(jlen, ab[j * ldab..], 1, ab[(n - j - jlen) * ldab + jlen - 1..], min(-ldab + 1, -1))
            }
        } else {
            for (j in 0..n / 2) {
                let jlen = min(n - 2 * j - 1, kd + 1)
                nFloat64.dswap(jlen, ab[j * ldab + kd..], max(ldab - 1, 1), ab[(n - j - 1) * ldab + kd + 1 - jlen..], -1
                )
            }
        }
    }

    return (diag, b)
}
