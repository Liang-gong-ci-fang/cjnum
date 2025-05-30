package cjnum.lapack.testlapack

import cjnum.blas.blas64.nFloat64Implementation
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dgtsvTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)
	for (n in [0, 1, 2, 3, 4, 5, 10, 25, 50]) {
		for (nrhs in [0, 1, 2, 3, 4, 10]) {
			for (ldb in [max(1, nrhs), nrhs + 3]) {
				if (!dgtsvTest(lFloat64, rnd, n, nrhs, ldb)) {
                    return false
                }
			}
		}
	}
    return true
}

func dgtsvTest(lFloat64: LFloat64, rnd: Random, n: Int64, nrhs: Int64, ldb: Int64): Bool {
	let tol = 1e-14
	let extra = 10

	let name = "Case n=${n},nrhs=${nrhs},ldb=${ldb}"

	if (n == 0) {
		let ok = lFloat64.dgtsv(n, nrhs, [], [], [], [], ldb)
		if (!ok) {
			throw Exception("${name}: unexpected failure for zero size matrix")
		}
		return true
	}

	let d = randomSlice(n + 1 + extra, rnd)
	let dCopy = Array<Float64>(d.size, item: 0.0)
	copy(dCopy, d)
	var (dl, dlCopy, du, duCopy) = (Array<Float64>(), Array<Float64>(), Array<Float64>(), Array<Float64>())
	if (n > 1) {
		dl = randomSlice(n + extra, rnd)
		dlCopy = Array<Float64>(dl.size, item: 0.0)
		copy(dlCopy, dl)

		du = randomSlice(n + extra, rnd)
		duCopy = Array<Float64>(du.size, item: 0.0)
		copy(duCopy, du)
	}

	let b = randomGeneral(n, nrhs, ldb, rnd)
	let got = cloneGeneral(b)

	let ok = lFloat64.dgtsv(n, nrhs, dl, d, du, got.data, got.stride)
	if (!ok) {
		throw Exception("${name}: unexpected failure in Dgtsv")
	}

	dlagtm(NoTrans, n, nrhs, 1.0, dlCopy, dCopy, duCopy, got.data, got.stride, -1.0, b.data, b.stride)

	let anorm = dlangt(MaxColumnSum, n, dlCopy, dCopy, duCopy)
	let nFloat64 = nFloat64Implementation()
	var resid = 0.0
	for (j in 0..nrhs) {
		let bnorm = nFloat64.dasum(n, b.data[j..], b.stride)
		let xnorm = nFloat64.dasum(n, got.data[j..], got.stride)
		resid = max(resid, bnorm / anorm / xnorm)
	}
	if (resid > tol) {
		throw Exception("${name}: unexpected result; resid=${resid},want<=${tol}")
	}
    return true
}
