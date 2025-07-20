package cjnum.lapack.testlapack

import std.math.*
import cjnum.blas.*
import cjnum.internal.asm.f64.*
import cjnum.lapack.*

func dlagtm(trans: Transpose, m: Int64, n: Int64, alpha: Float64, dl: Array<Float64>, d: Array<Float64>, du: Array<Float64>,
    b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64) {
    if (m == 0 || n == 0) {
        return
    }

    if (beta != 1.0) {
        if (beta == 0.0) {
            for (i in 0..m) {
                let ci = c[i * ldc..i * ldc + n]
                for (j in 0..ci.size) {
                    ci[j] = 0.0
                }
            }
        } else {
            for (i in 0..m) {
                let ci = c[i * ldc..i * ldc + n]
                for (j in 0..ci.size) {
                    ci[j] *= beta
                }
            }
        }
    }

    if (alpha == 0.0) {
        return
    }

    if (m == 1) {
        if (alpha == 1.0) {
            for (j in 0..n) {
                c[j] += d[0] * b[j]
            }
        } else {
            for (j in 0..n) {
                c[j] += alpha * d[0] * b[j]
            }
        }
        return
    }

    let (dl_, du_) = if (trans != NoTrans) {
        (du, dl)
    } else {
        (dl, du)
    }

    if (alpha == 1.0) {
        for (j in 0..n) {
            c[j] += d[0] * b[j] + du_[0] * b[ldb + j]
        }

        for (i in 1..m - 1) {
            for (j in 0..n) {
                c[i * ldc + j] += dl_[i - 1] * b[(i - 1) * ldb + j] + d[i] * b[i * ldb + j] + du_[i] * b[(i + 1) * ldb +
                        j]
            }
        }

        for (j in 0..n) {
            c[(m - 1) * ldc + j] += dl_[m - 2] * b[(m - 2) * ldb + j] + d[m - 1] * b[(m - 1) * ldb + j]
        }
    } else {
        for (j in 0..n) {
            c[j] += alpha * (d[0] * b[j] + du_[0] * b[ldb + j])
        }

        for (i in 1..m - 1) {
            for (j in 0..n) {
                c[i * ldc + j] += alpha * (dl_[i - 1] * b[(i - 1) * ldb + j] + d[i] * b[i * ldb + j] + du_[i] * b[(i + 1) *
                        ldb + j])
            }
        }

        for (j in 0..n) {
            c[(m - 1) * ldc + j] += alpha * (dl_[m - 2] * b[(m - 2) * ldb + j] + d[m - 1] * b[(m - 1) * ldb + j])
        }
    }
}

func dlangt(norm: MatrixNorm, n: Int64, dl: Array<Float64>, d: Array<Float64>, du: Array<Float64>): Float64 {
    if (n == 0) {
        return 0.0
    }

    let dl_ = dl[0..n - 1]
    let d_ = d[0..n]
    let du_ = du[0..n - 1]

    var anorm = 0.0

    if (norm == MaxAbs) {
        for (diag in [dl_, d_, du_]) {
            for (i in 0..diag.size) {
                var di = diag[i]
                if (di.isNaN()) {
                    return di
                }
                di = abs(di)
                if (di > anorm) {
                    anorm = di
                }
            }
        }
    } else if (norm == MaxColumnSum) {
        if (n == 1) {
            return abs(d_[0])
        }

        anorm = abs(d_[0]) + abs(dl_[0])
        if (anorm.isNaN()) {
            return anorm
        }

        var tmp = abs(du_[n - 2]) + abs(d_[n - 1])
        if (tmp.isNaN()) {
            return tmp
        }
        if (tmp > anorm) {
            anorm = tmp
        }

        for (i in 1..n - 1) {
            tmp = abs(du_[i - 1]) + abs(d_[i]) + abs(dl_[i])
            if (tmp.isNaN()) {
                return tmp
            }
            if (tmp > anorm) {
                anorm = tmp
            }
        }
    } else if (norm == MaxRowSum) {
        if (n == 1) {
            return abs(d_[0])
        }

        anorm = abs(d_[0]) + abs(du_[0])
        if (anorm.isNaN()) {
            return anorm
        }

        var tmp = abs(dl_[n - 2]) + abs(d_[n - 1])
        if (tmp.isNaN()) {
            return tmp
        }
        if (tmp > anorm) {
            anorm = tmp
        }

        for (i in 1..n - 1) {
            tmp = abs(dl_[i - 1]) + abs(d_[i]) + abs(du_[i])
            if (tmp.isNaN()) {
                return tmp
            }
            if (tmp > anorm) {
                anorm = tmp
            }
        }
    } else if (norm == Frobenius) {
        throw Exception("not implemented")
    } else {
        throw Exception("invalid norm")
    }

    return anorm
}

func dlansy(norm: MatrixNorm, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Float64 {
    if (n == 0) {
        return 0.0
    }
    let work = Array<Float64>(n, repeat: 0.0)
    match {
        case norm == MaxAbs =>
            if (uplo == Upper) {
                var max = 0.0
                for (i in 0..n) {
                    for (j in i..n) {
                        let v = abs(a[i * lda + j])
                        if (v.isNaN()) {
                            return Float64.NaN
                        }
                        if (v > max) {
                            max = v
                        }
                    }
                }
                return max
            }
            var max = 0.0
            for (i in 0..n) {
                for (j in 0..=i) {
                    let v = abs(a[i * lda + j])
                    if (v.isNaN()) {
                        return Float64.NaN
                    }
                    if (v > max) {
                        max = v
                    }
                }
            }
            return max
        case norm == MaxRowSum || norm == MaxColumnSum =>
            for (i in 0..n) {
                work[i] = 0.0
            }
            if (uplo == Upper) {
                for (i in 0..n) {
                    work[i] += abs(a[i * lda + i])
                    for (j in i + 1..n) {
                        let v = abs(a[i * lda + j])
                        work[i] += v
                        work[j] += v
                    }
                }
            } else {
                for (i in 0..n) {
                    for (j in 0..i) {
                        let v = abs(a[i * lda + j])
                        work[i] += v
                        work[j] += v
                    }
                    work[i] += abs(a[i * lda + i])
                }
            }
            var max = 0.0
            for (i in 0..n) {
                let v = work[i]
                if (v.isNaN()) {
                    return Float64.NaN
                }
                if (v > max) {
                    max = v
                }
            }
            return max
        case norm == Frobenius => throw Exception("not implemented") // TODO
        case _ => throw Exception("invalid norm")
    }
}

func dlange(norm: MatrixNorm, m: Int64, n: Int64, a: Array<Float64>, lda: Int64): Float64 {
    if (m == 0 || n == 0) {
        return 0.0
    }
    var value = 0.0
    if (norm == MaxAbs) {
        for (i in 0..m) {
            for (j in 0..n) {
                if (value.isNaN()) {  // TODO cj1.0.0
                    continue
                }
                value = max(value, abs(a[i * lda + j]))
            }
        }
    } else if (norm == MaxColumnSum) {
        let work = Array<Float64>(n, repeat: 0.0)
        for (i in 0..m) {
            for (j in 0..n) {
                work[j] += abs(a[i * lda + j])
            }
        }
        for (i in 0..n) {
            value = max(value, work[i])
        }
    } else if (norm == MaxRowSum) {
        for (i in 0..m) {
            var sum = 0.0
            for (j in 0..n) {
                sum += abs(a[i * lda + j])
            }
            value = max(value, sum)
        }
    } else if (norm == Frobenius) {
        for (i in 0..m) {
            let row = l2NormUnitary(a[i * lda..i * lda + n])
            value = hypot(value, row)
        }
    } else {
        throw Exception("invalid norm")
    }
    return value
}

//TODO func dlansb(norm lapack.MatrixNorm, uplo blas.Uplo, n, kd int, ab []float64, ldab int, work []float64) float64

func dlantr(norm: MatrixNorm, uplo: Uplo, diag: Diag, m: Int64, n: Int64, a: Array<Float64>, lda: Int64,
    work: Array<Float64>): Float64 {
    let minmn = min(m, n)
    if (minmn == 0) {
        return 0.0
    }
    match {
        case norm == MaxAbs =>
            if (diag == UnitDiag) {
                var value: Float64 = 1.0
                if (uplo == Upper) {
                    for (i in 0..m) {
                        for (j in i + 1..n) {
                            let tmp = abs(a[i * lda + j])
                            if (tmp.isNaN()) {
                                return tmp
                            }
                            if (tmp > value) {
                                value = tmp
                            }
                        }
                    }
                    return value
                }
                for (i in 1..m) {
                    for (j in 0..min(i, n)) {
                        let tmp = abs(a[i * lda + j])
                        if (tmp.isNaN()) {
                            return tmp
                        }
                        if (tmp > value) {
                            value = tmp
                        }
                    }
                }
                return value
            }
            var value: Float64 = 0.0
            if (uplo == Upper) {
                for (i in 0..m) {
                    for (j in i..n) {
                        let tmp = abs(a[i * lda + j])
                        if (tmp.isNaN()) {
                            return tmp
                        }
                        if (tmp > value) {
                            value = tmp
                        }
                    }
                }
                return value
            }
            for (i in 0..m) {
                for (j in 0..=min(i, n - 1)) {
                    let tmp = abs(a[i * lda + j])
                    if (tmp.isNaN()) {
                        return tmp
                    }
                    if (tmp > value) {
                        value = tmp
                    }
                }
            }
            return value

        case norm == MaxColumnSum =>
            if (diag == UnitDiag) {
                for (i in 0..minmn) {
                    work[i] = 1.0
                }
                for (i in minmn..n) {
                    work[i] = 0.0
                }
                if (uplo == Upper) {
                    for (i in 0..m) {
                        for (j in i + 1..n) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                } else {
                    for (i in 1..m) {
                        for (j in 0..min(i, n)) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                }
            } else {
                for (i in 0..n) {
                    work[i] = 0.0
                }
                if (uplo == Upper) {
                    for (i in 0..m) {
                        for (j in i..n) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                } else {
                    for (i in 0..m) {
                        for (j in 0..=min(i, n - 1)) {
                            work[j] += abs(a[i * lda + j])
                        }
                    }
                }
            }
            var max: Float64 = 0.0
            for (v in work[..n]) {
                if (v.isNaN()) {
                    return Float64.NaN
                }
                if (v > max) {
                    max = v
                }
            }
            return max

        case norm == MaxRowSum =>
            var maxsum: Float64 = 0.0
            if (diag == UnitDiag) {
                if (uplo == Upper) {
                    for (i in 0..m) {
                        var sum: Float64 = 0.0
                        if (i < minmn) {
                            sum = 1.0
                        }
                        for (j in i + 1..n) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return Float64.NaN
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                } else {
                    for (i in 0..m) {
                        var sum: Float64 = 0.0
                        if (i < minmn) {
                            sum = 1.0
                        }
                        for (j in 0..min(i, n)) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return Float64.NaN
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                }
            } else {
                if (uplo == Upper) {
                    for (i in 0..m) {
                        var sum: Float64 = 0.0
                        for (j in i..n) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return sum
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                } else {
                    for (i in 0..m) {
                        var sum: Float64 = 0.0
                        for (j in 0..=min(i, n - 1)) {
                            sum += abs(a[i * lda + j])
                        }
                        if (sum.isNaN()) {
                            return sum
                        }
                        if (sum > maxsum) {
                            maxsum = sum
                        }
                    }
                    return maxsum
                }
            }

        case norm == Frobenius => throw Exception("not implemented") // TODO
        case _ => throw Exception("invalid norm")
    }
}

//TODO func dlantb(norm: MatrixNorm, uplo: Uplo, diag: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, work: Array<Float64>): Float64

func dlanst(norm: MatrixNorm, n: Int64, d: Array<Float64>, e: Array<Float64>): Float64 {
    if (n == 0) {
        return 0.0
    }
    var value = 0.0
    match {
        case norm == MaxAbs => if (n == 1) {
            value = abs(d[0])
        } else {
            for (di in d[..n]) {
                value = max(value, abs(di))
            }
            for (ei in e[..n - 1]) {
                value = max(value, abs(ei))
            }
        }
        case norm == MaxColumnSum || norm == MaxRowSum => if (n == 1) {
            value = abs(d[0])
        } else {
            value = abs(d[0]) + abs(e[0])
            value = max(value, abs(d[n - 1]) + abs(e[n - 2]))
            for (i in 1..n - 1) {
                let sum = abs(d[i]) + abs(e[i]) + abs(e[i - 1])
                value = max(value, sum)
            }
        }
        case norm == Frobenius => throw Exception("not implemented") // TODO
        case _ => throw Exception("invalid norm")
    }
    return value
}
