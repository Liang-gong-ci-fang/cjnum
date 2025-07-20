package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*

public func dlarfTest(lFloat64: LFloat64): Bool {
    for (side in [Left, Right]) {
        runDlarfTest(lFloat64, side)
    }
    return true
}

func runDlarfTest(lFloat64: LFloat64, side: Side): Unit {
	let rnd = Random(1)
	for (m in [0, 1, 2, 3, 4, 5, 10]) {
		for (n in [0, 1, 2, 3, 4, 5, 10]) {
			for (incv in [1, 4]) {
				for (ldc in [max(1, n), n + 3]) {
					for (nnzv in [0, 1, 2]) {
						for (nnzc in [0, 1, 2]) {
							for (tau in [0.0, rnd.nextGaussianFloat64()]) {
								_dlarfTest(lFloat64, rnd, side, m, n, incv, ldc, nnzv, nnzc, tau)
							}
						}
					}
				}
			}
		}
	}
}

func _dlarfTest(lFloat64: LFloat64, rnd: Random, side: Side, m: Int64, n: Int64, incv: Int64, ldc: Int64, nnzv: Int64, nnzc: Int64, tau: Float64): Unit {
	const tol = 1e-14

	let c = Array<Float64>(m * ldc, { _ => rnd.nextGaussianFloat64() })
	if (nnzc == 0) {
		for (i in 0..m) {
			for (j in 0..n) {
				c[i * ldc + j] = 0.0
			}
		}
    } else if (nnzc == 1) {
		if (side == Left) {
			for (i in 0..m) {
				for (j in n / 2..n) {
					c[i * ldc + j] = 0.0
				}
			}
		} else {
			for (i in m / 2..m) {
				for (j in 0..n) {
					c[i * ldc + j] = 0.0
				}
			}
		}
    }

	let cCopy = c.clone()

	let work = if (side == Left) {
		Array<Float64>(n, repeat: 0.0)
	} else {
		Array<Float64>(m, repeat: 0.0)
	}

	var vlen = n
	if (side == Left) {
		vlen = m
	}
	vlen = max(1, vlen)
	let v = Array<Float64>(1 + (vlen - 1) * incv, { _ => rnd.nextGaussianFloat64() })

	if (nnzv == 0) {
        for (i in 0..vlen) {
            v[i * incv] = 0.0
        }
    } else if (nnzv == 1) {
        for (i in vlen / 2..vlen) {
            v[i * incv] = 0.0
        }
    }

	let vCopy = v.clone()

	lFloat64.dlarf(side, m, n, v, incv, tau, c, ldc, work)
	let got = c

	let name = "m=${m},n=${n},incv=${incv},tau=${tau},ldc=${ldc}"

	if (!fEqual(v, vCopy)) {
		throw Exception("${name}: unexpected modification of v")
	}
	if (tau == 0.0 && !fEqual(got, cCopy)) {
		throw Exception("${name}: unexpected modification of C")
	}

	if (m == 0 || n == 0 || tau == 0.0) {
		return
	}

	let nFloat64 = nFloat64Implementation()

	let want = Array<Float64>(cCopy.size, repeat: 0.0)
	if (side == Left) {
		let vtc = Array<Float64>(n, repeat: 0.0)
		nFloat64.dgemv(Trans, m, n, -tau, cCopy, ldc, v, incv, 0.0, vtc, 1)

		for (i in 0..m) {
			for (j in 0..n) {
				want[i * ldc + j] = cCopy[i * ldc + j] + v[i * incv] * vtc[j]
			}
		}
	} else {
		let cv = Array<Float64>(m, repeat: 0.0)
		nFloat64.dgemv(NoTrans, m, n, -tau, cCopy, ldc, v, incv, 0.0, cv, 1)

		for (i in 0..m) {
			for (j in 0..n) {
				want[i * ldc + j] = cCopy[i * ldc + j] + cv[i] * v[j * incv]
			}
		}
	}
	let diff = Array<Float64>(m * n, repeat: 0.0)
    for (i in 0..m) {
        for (j in 0..n) {
			diff[i * n + j] = got[i * ldc + j] - want[i * ldc + j]
		}
	}
	let resid = dlange(MaxColumnSum, m, n, diff, n)
	if (resid > tol * Float64(max(m, n))) {
		throw Exception("${name}: unexpected result; resid=${resid}, want<=${tol * Float64(max(m, n))}")
	}
}