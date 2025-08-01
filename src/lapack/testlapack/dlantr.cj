package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.floats.*
import cjnum.util.*

public func dlantrTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)
	for (m in [0, 1, 2, 3, 4, 5, 10]) {
		for (n in [0, 1, 2, 3, 4, 5, 10]) {
			for (uplo in [Lower, Upper]) {
				if ((uplo == Upper) && m > n) {
					continue
				}
				if (uplo == Lower && n > m) {
					continue
				}
				for (diag in [NonUnitDiag, UnitDiag]) {
					for (lda in [max(1, n), n + 3]) {
						_dlantrTest(lFloat64, rnd, uplo, diag, m, n, lda)
					}
				}
			}
		}
	}
    return true
}

func _dlantrTest( lFloat64: LFloat64, rnd: Random, uplo: Uplo, diag: Diag, m: Int64, n: Int64, lda: Int64): Unit {
	const tol = 1e-14

	let a = Array<Float64>(max(0, (m - 1) * lda + n), { _ => rnd.nextGaussianFloat64() })
	let rowsum = Array<Float64>(m, repeat: 0.0)
	let colsum = Array<Float64>(n, repeat: 0.0)
	var (frobWant, maxabsWant) = (0.0, 0.0)
	if (diag == UnitDiag) {
		for (i in 0..min(m, n)) {
			rowsum[i] = 1.0
			colsum[i] = 1.0
		}
		frobWant = Float64(min(m, n))
		if (min(m, n) > 0) {
			maxabsWant = 1.0
		}
	}
	if (uplo == Upper) {
		for (i in 0..min(m, n)) {
			var start = i
			if (diag == UnitDiag) {
				start = i + 1
			}
			for (j in start..n) {
				let aij = 2.0 * rnd.nextFloat64() - 1.0
				a[i*lda+j] = aij
				rowsum[i] += abs(aij)
				colsum[j] += abs(aij)
				maxabsWant = nanMax(maxabsWant, abs(aij))
				frobWant += aij * aij
			}
		}
	} else {
		for (i in 0..m) {
			var end = i
			if (diag == UnitDiag) {
				end = i - 1
			}
			for (j in 0..=min(end, n - 1)) {
				let aij = 2.0 * rnd.nextFloat64() - 1.0
				a[i*lda+j] = aij
				rowsum[i] += abs(aij)
				colsum[j] += abs(aij)
				maxabsWant = nanMax(maxabsWant, abs(aij))
				frobWant += aij * aij
			}
		}
	}
	frobWant = sqrt(frobWant)
	var (maxcolsumWant, maxrowsumWant) = (0.0, 0.0)
	if (n > 0) {
		maxcolsumWant = fMax(colsum)
	}
	if (m > 0) {
		maxrowsumWant = fMax(rowsum)
	}

	let aCopy = a.clone()

	for (norm in [MaxAbs, MaxColumnSum, MaxRowSum, Frobenius]) {
		let name = "norm=${norm},uplo=${uplo},diag=${diag},m=${m},n=${n},lda=${lda}"

		let work: Array<Float64> = if (norm == MaxColumnSum) {
			Array<Float64>(n, repeat: 0.0)
		} else {
            []
        }

		let normGot = lFloat64.dlantr(norm, uplo, diag, m, n, a, lda, work)

		if (!fEqual(a, aCopy)) {
			throw Exception("${name}: unexpected modification of a")
		}

		if (norm == MaxAbs) {
			if (normGot != maxabsWant) {
				throw Exception("${name}: unexpected result; got ${normGot}, want ${maxabsWant}")
			}
			continue
		}

		var normWant = 0.0
		if (norm == MaxColumnSum) {
            normWant = maxcolsumWant
        } else if (norm == MaxRowSum) {
            normWant = maxrowsumWant
        } else if (norm == Frobenius) {
            normWant = frobWant
        }

		if (abs(normGot - normWant) >= tol) {
			throw Exception("${name}: unexpected result; got ${normGot}, want ${normWant}")
		}
	}
}
