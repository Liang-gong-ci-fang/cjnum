package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.floats.*
import cjnum.util.*

public func dlansbTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
	for (n in [0, 1, 2, 3, 4, 5, 10]) {
		for (kd in [0, (n + 1) / 4, (3*n - 1) / 4, (5*n + 1) / 4]) {
			for (uplo in [Upper, Lower]) {
				for (ldab in [kd + 1, kd + 1 + 7]) {
					_dlansbTest(lFloat64, rnd, uplo, n, kd, ldab)
				}
			}
		}
	}
    return true
}

func _dlansbTest(lFloat64: LFloat64, rnd: Random, uplo: Uplo, n: Int64, kd: Int64, ldab: Int64): Unit {
    const tol = 1e-15

	let ab = Array<Float64>(max(0, (n - 1) * ldab + kd + 1), repeat: 0.0)
	let rowsum = Array<Float64>(n, repeat: 0.0)
	let colsum = Array<Float64>( n, repeat: 0.0)
	var (frobWant, maxabsWant) = (0.0, 0.0)
	if (uplo == Upper) {
		for (i in 0..n) {
			for (jb in 0..min(n - i, kd + 1)) {
				let aij = 2.0 * rnd.nextFloat64() - 1.0
				ab[i * ldab + jb] = aij

				let j = jb + i
				colsum[j] += abs(aij)
				rowsum[i] += abs(aij)
				maxabsWant = nanMax(maxabsWant, abs(aij))
				frobWant += aij * aij
				if (i != j) {
					colsum[i] += abs(aij)
					rowsum[j] += abs(aij)
					frobWant += aij * aij
				}
			}
		}
	} else {
		for (i in 0..n) {
			for (jb in max(0, kd-i)..kd + 1) {
				let aij = 2.0 * rnd.nextFloat64() - 1.0
				ab[i * ldab + jb] = aij

				let j = jb - kd + i
				colsum[j] += abs(aij)
				rowsum[i] += abs(aij)
				maxabsWant = nanMax(maxabsWant, abs(aij))
				frobWant += aij * aij
				if (i != j) {
					colsum[i] += abs(aij)
					rowsum[j] += abs(aij)
					frobWant += aij * aij
				}
			}
		}
	}
	frobWant = sqrt(frobWant)
	var (maxcolsumWant, maxrowsumWant) = (0.0, 0.0)
	if (n > 0) {
		maxcolsumWant = fMax(colsum)
		maxrowsumWant = fMax(rowsum)
	}

	let abCopy = Array<Float64>(ab.size, repeat: 0.0)
	copy(abCopy, ab)

	let work = Array<Float64>(n, repeat: 0.0)
	var (maxcolsumGot, maxrowsumGot) = (0.0, 0.0)
	for (norm in [MaxAbs, MaxColumnSum, MaxRowSum, Frobenius]) {
		let name = "norm=${norm},uplo=${uplo},n=${n},kd=${kd},ldab=${ldab}"

		let normGot = lFloat64.dlansb(norm, uplo, n, kd, ab, ldab, work)

		if (!fEqual(ab, abCopy)) {
			throw Exception("${name}: unexpected modification of ab")
		}

		if (norm == MaxAbs) {
			if (normGot != maxabsWant) {
				throw Exception("${name}: unexpected result; got ${normGot}, want ${maxabsWant}")
			}
			continue
		}

		var normWant  = 0.0
        if (norm == MaxColumnSum) {
			normWant = maxcolsumWant
			maxcolsumGot = normGot
        } else if (norm == MaxRowSum) {
			normWant = maxrowsumWant
			maxrowsumGot = normGot
        } else if (norm == Frobenius) {
			normWant = frobWant
        }
		if (abs(normGot-normWant) > tol * Float64(n)) {
			throw Exception("${name}: unexpected result; got ${normGot}, want ${normWant}")
		}
	}

	if (maxcolsumGot != maxrowsumGot) {
		let name = "uplo=${uplo},n=${n},kd=${kd},ldab=${ldab}"
		throw Exception("${name}: unexpected mismatch between MaxColSum and MaxRowSum norms of A; MaxColSum ${maxcolsumGot}, MaxRowSum ${maxrowsumWant}")
	}
}