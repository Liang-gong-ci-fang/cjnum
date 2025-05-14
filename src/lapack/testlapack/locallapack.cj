package cjnum.lapack.testlapack

import std.math.*
import cjnum.blas.*
import cjnum.internal.asm.f64.*
import cjnum.lapack.*

// TODO

func dlansy(norm: MatrixNorm, uplo: Uplo, n: Int64, a: Array<Float64>, lda: Int64): Float64 {
    if (n == 0) {
        return 0.0
    }
    let work = Array<Float64>(n, item: 0.0)
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
	case norm == Frobenius =>
		throw Exception("not implemented")
	case _ =>
		throw Exception("invalid norm")
	}
}

func dlange(norm: MatrixNorm, m: Int64, n: Int64, a: Array<Float64>, lda: Int64): Float64 {
    if (m == 0 || n == 0) {
        return 0.0
    }
    var value = 0.0
    match {
        case norm == MaxAbs =>
            for (i in 0..m) {
                for (j in 0..n) {
                    value = max(value, abs(a[i * lda + j]))
                }
            }
        case norm == MaxColumnSum =>
            let work = Array<Float64>(n, item: 0.0)
            for (i in 0..m) {
                for (j in 0..n) {
                    work[j] += abs(a[i * lda + j])
                }
            }
            for (i in 0..n) {
                value = max(value, work[i])
            }
        case norm == MaxRowSum =>
            for (i in 0..m) {
                var sum = 0.0
                for (j in 0..n) {
                    sum += abs(a[i * lda + j])
                }
                value = max(value, sum)
            }
        case norm == Frobenius =>
            for (i in 0..m) {
                let row = l2NormUnitary(a[i * lda..i * lda + n])
                value = hypot(value, row)
            }
        case _ => throw Exception("invalid norm")
    }
    return value
}

func dlanst(norm: MatrixNorm, n: Int64, d: Array<Float64>, e: Array<Float64>): Float64 {
    if (n == 0) {
        return 0.0
    }
    var value = 0.0
    match {
        case norm == MaxAbs =>
            if (n == 1) {
                value = abs(d[0])
            } else {
                for (di in d[..n]) {
                    value = max(value, abs(di))
                }
                for (ei in e[..n - 1]) {
                    value = max(value, abs(ei))
                }
            }
        case norm == MaxColumnSum || norm == MaxRowSum =>
            if (n == 1) {
                value = abs(d[0])
            } else {
                value = abs(d[0]) + abs(e[0])
                value = max(value, abs(d[n - 1]) + abs(e[n - 2]))
                for (i in 1..n - 1) {
                    let sum = abs(d[i]) + abs(e[i]) + abs(e[i - 1])
                    value = max(value, sum)
                }
            }
        case norm == Frobenius =>
            throw Exception("not implemented")
        case _ => throw Exception("invalid norm")
    }
    return value
}

// TODO
