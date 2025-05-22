package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*
import cjnum.floats.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64

public func dlasrTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)
	const pi = 3.14159265358979323846264338327950288419716939937510582097494459
	for (side in [Left, Right]) {
		for (pivot in [Variable, Top, Bottom]) {
			for (direct in [Forward, Backward]) {
				for (test in [
					[5, 5, 0],
					[5, 10, 0],
					[10, 5, 0],

					[5, 5, 20],
					[5, 10, 20],
					[10, 5, 20]
				]) {
					let m = test[0]
					let n = test[1]
					var lda = test[2]
					if (lda == 0) {
						lda = n
					}
					let a = Array<Float64>(m * lda, { _ => rnd.nextFloat64() })

					let (s, c) = if (side == Left) {
						(Array<Float64>(m - 1, item: 0.0), Array<Float64>(m - 1, item: 0.0))
					} else {
						(Array<Float64>(n - 1, item: 0.0), Array<Float64>(n - 1, item: 0.0))
					}
					for (k in 0..s.size) {
						let theta = rnd.nextFloat64() * 2.0 * pi
						s[k] = sin(theta)
						c[k] = cos(theta)
					}
					let aCopy = Array<Float64>(a.size, item: 0.0)
					copy(a, aCopy)

					lFloat64.dlasr(side, pivot, direct, m, n, c, s, a, lda)

					var pSize = m
					if (side == Right) {
						pSize = n
					}
					let p = General64(
						rows:   pSize,
						cols:   pSize,
						stride: pSize,
						data:   Array<Float64>(pSize * pSize, item: 0.0)
					)
					let pk = General64(
						rows:   pSize,
						cols:   pSize,
						stride: pSize,
						data:   Array<Float64>(pSize * pSize, item: 0.0)
					)
					let ptmp = General64(
						rows:   pSize,
						cols:   pSize,
						stride: pSize,
						data:   Array<Float64>(pSize * pSize, item: 0.0)
					)

					for (i in 0..pSize) {
						p.data[i * p.stride + i] = 1.0
						ptmp.data[i * p.stride + i] = 1.0
					}

					for (k in 0..s.size) {
						for (i in 0..p.data.size) {
							pk.data[i] = .0
						}
						for (i in 0..pSize) {
							pk.data[i * p.stride + i] = 1.0
						}
						if (pivot == Variable) {
							pk.data[k * p.stride + k] = c[k]
							pk.data[k * p.stride + k + 1] = s[k]
							pk.data[(k + 1) * p.stride + k] = -s[k]
							pk.data[(k + 1) * p.stride + k + 1] = c[k]
						} else if (pivot == Top) {
							pk.data[0] = c[k]
							pk.data[k + 1] = s[k]
							pk.data[(k + 1) * p.stride] = -s[k]
							pk.data[(k + 1) * p.stride + k + 1] = c[k]
						} else if (pivot == Bottom) {
							pk.data[(pSize - 1 - k) * p.stride + pSize - k - 1] = c[k]
							pk.data[(pSize - 1 - k) * p.stride + pSize - 1] = s[k]
							pk.data[(pSize - 1) * p.stride + pSize - 1 - k] = -s[k]
							pk.data[(pSize - 1) * p.stride + pSize - 1] = c[k]
						}
						if (direct == Forward) {
							gemm(NoTrans, NoTrans, 1.0, pk, ptmp, 0.0, p)
						} else {
							gemm(NoTrans, NoTrans, 1.0, ptmp, pk, 0.0, p)
						}
						copy(ptmp.data, p.data)
					}

					let aMat = General64(
						rows:   m,
						cols:   n,
						stride: lda,
						data:   Array<Float64>(m * lda, item: 0.0)
					)
					copy(a, aCopy)
					let newA = General64(
						rows:   m,
						cols:   n,
						stride: lda,
						data:   Array<Float64>(m * lda, item: 0.0)
					)
					if (side == Left) {
						gemm(NoTrans, NoTrans, 1.0, p, aMat, 0.0, newA)
					} else {
						gemm(NoTrans, NoTrans, 1.0, aMat, p, 0.0, newA)
					}
					if (!fEqualApprox(newA.data, a, 1e-12)) {
						throw Exception("A update mismatch")
					}
				}
			}
		}
	}
    return true
}
