package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.floats.*

public func dlantbTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)
	for (norm in [MaxAbs, MaxRowSum, MaxColumnSum, Frobenius]) {
		for (uplo in [Lower, Upper]) {
			for (diag in [NonUnitDiag, UnitDiag]) {
				for (n in [0, 1, 2, 3, 4, 5, 10]) {
					for (k in [0, 1, 2, 3, n, n + 2]) {
						for (lda in [k + 1, k + 3]) {
							for (_ in 0..10) {
								_dlantbTest(lFloat64, rnd, norm, uplo, diag, n, k, lda)
							}
						}
					}
				}
			}
		}
	}
    return true
}

func _dlantbTest(lFloat64: LFloat64, rnd: Random, norm: MatrixNorm, uplo: Uplo, diag: Diag, n: Int64, k: Int64, lda: Int64): Unit {
	const tol = 1e-14

	let name = "n=${n},${k},lda=${lda}"

	if (n == 0) {
		let got = lFloat64.dlantb(norm, uplo, diag, n, k, [], lda, [])
		if (got != 0.0) {
			throw Exception("${name}: unexpected result for zero-sized matrix")
		}
		return
	}

	let a = Array<Float64>(max(0, (n - 1) * lda + k + 1), repeat: 0.0)
	if (rnd.nextFloat64() < 0.5) {
		for (i in 0..a.size) {
			a[i] = rnd.nextFloat64() - 0.5
		}
	} else {
		for (i in 0..a.size) {
			a[i] = 4.0 * rnd.nextFloat64() - 2.0
		}
	}

	if (rnd.nextFloat64() < 0.5) {
		a[Int64(rnd.nextUInt64(UInt64(a.size)))] = Float64.NaN
	}

	let aCopy = a.clone()
	var work: Array<Float64> = []
	if (norm == MaxColumnSum) {
		work = Array<Float64>(n, repeat: 0.0)
	}

	for (i in 0..work.size) {
		work[i] = rnd.nextGaussianFloat64()
	}

	let got = lFloat64.dlantb(norm, uplo, diag, n, k, a, lda, work)

	if (!fSame(a, aCopy)) {
		throw Exception("${name}: unexpected modification of a")
	}

	let ldaGen = n
	let aGen = Array<Float64>(n * ldaGen, repeat: 0.0)
	if (uplo == Upper) {
		for (i in 0..n) {
			for (j in 0..min(n - i, k + 1)) {
				aGen[i * ldaGen + i + j] = a[i * lda + j]
			}
		}
	} else {
		for (i in 0..n) {
			for (j in max(0, k - i)..k + 1) {
				aGen[i * ldaGen + i - (k - j)] = a[i * lda + j]
			}
		}
	}
	if (diag == UnitDiag) {
		for (i in 0..n) {
			aGen[i * ldaGen + i] = 1.0
		}
	}
	let want = dlange(norm, n, n, aGen, ldaGen)

	if (want.isNaN()) {
		if (!got.isNaN()) {
			throw Exception("${name}: unexpected result with NaN element; got ${got}, want ${want}")
		}
		return
	}

	if (norm == MaxAbs) {
		if (got != want) {
			throw Exception("${name}: unexpected result; got ${got}, want ${want}")
		}
		return
	}
	let diff = abs(got - want)
	if (diff > tol) {
		throw Exception("${name}: unexpected result; got ${got}, want ${want}, diff=${diff}")
	}
}
