package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*

public func dlanstTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)
	for (norm in [MaxAbs, MaxColumnSum, MaxRowSum, Frobenius]) {
		for (n in [1, 3, 10, 100]) {
			for (cas in 0..100) {
				let d = Array<Float64>(n, { _ => rnd.nextGaussianFloat64()})
				let e = Array<Float64>(n - 1, { _ => rnd.nextGaussianFloat64()})

				let m = n
				let lda = n
				let a = Array<Float64>(m * lda, repeat: 0.0)
				for (i in 0..n) {
					a[i * lda + i] = d[i]
				}
				for (i in 0..n - 1) {
					a[i * lda + i + 1] = e[i]
					a[(i + 1) * lda + i] = e[i]
				}

				let work = Array<Float64>(n, repeat: 0.0)
				let syNorm = lFloat64.dlanst(norm, n, d, e)
				let geNorm = lFloat64.dlange(norm, m, n, a, lda, work)
				if (abs(syNorm - geNorm) > 1e-12) {
					throw Exception("Norm mismatch: norm = ${norm}, cas = ${cas}, n = ${n}. Want ${geNorm}, got ${syNorm}.")
				}
			}
		}
	}
    return true
}