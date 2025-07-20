package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dgehrdTestStruct {
    dgehrdTestStruct(
        let n: Int64, 
        let ilo: Int64, 
        let ihi: Int64
    ) {}
}

public func dgehrdTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)

	for (n in [1, 2, 3, 4, 5, 10, 34]) {
		for (extra in [0, 13]) {
			for (optwork in [true, false]) {
				for (_ in 0..10) {
					var ilo = rnd.nextInt64(n)
					var ihi = rnd.nextInt64(n)
					if (ilo > ihi) {
						(ilo, ihi) = (ihi, ilo)
					}
					if (!testDgehrd(lFloat64, n, ilo, ihi, extra, optwork, rnd)) {
                        return false
                    }
				}
			}
		}
	}

	for (test in [
		dgehrdTestStruct(0, 0, -1),

		dgehrdTestStruct(68, 0, 63),
		dgehrdTestStruct(68, 0, 64),
		dgehrdTestStruct(68, 0, 65),
		dgehrdTestStruct(68, 0, 66),
		dgehrdTestStruct(68, 0, 67),

		dgehrdTestStruct(132, 2, 129),
		dgehrdTestStruct(132, 1, 129),
		dgehrdTestStruct(132, 0, 129),
		dgehrdTestStruct(132, 1, 130),
		dgehrdTestStruct(132, 0, 130),
		dgehrdTestStruct(132, 1, 131),
		dgehrdTestStruct(132, 0, 131),

		dgehrdTestStruct(260, 2, 257),
		dgehrdTestStruct(260, 1, 257),
		dgehrdTestStruct(260, 0, 257),
		dgehrdTestStruct(260, 0, 258),
		dgehrdTestStruct(260, 0, 259)
	]) {
		for (extra in [0, 13]) {
			for (optwork in [true, false]) {
				if (!testDgehrd(lFloat64, test.n, test.ilo, test.ihi, extra, optwork, rnd)) {
                    return false
                }
			}
		}
	}
    return true
}

func testDgehrd(lFloat64: LFloat64, n: Int64, ilo: Int64, ihi: Int64, extra: Int64, optwork: Bool, rnd: Random): Bool {
	const tol = 1e-13

	let a = randomGeneral(n, n, n+extra, rnd)
	let aCopy = General64(rows: a.rows, cols: a.cols, stride: a.stride, data: Array<Float64>(a.data.size, item: 0.0))
	copy(aCopy.data, a.data)

	var tau: Array<Float64> = Array<Float64>(0, item: 0.0)
	if (n > 1) {
		tau = nanSlice(n - 1)
	}

	var work: Array<Float64>
	if (optwork) {
		work = nanSlice(1)
		lFloat64.dgehrd(n, ilo, ihi, a.data, a.stride, tau, work, -1)
		work = nanSlice(Int64(work[0]))
	} else {
		work = nanSlice(max(1, n))
	}

	lFloat64.dgehrd(n, ilo, ihi, a.data, a.stride, tau, work, work.size)

	if (n == 0) {
		return true
	}

	let prefix = "n=${n},ilo=${ilo},ihi=${ihi},extra=${extra}"

	if (!generalOutsideAllNaN(a)) {
		throw Exception("${prefix}: out-of-range write to A\n${a.data}")
	}
	for (i in ilo..=ihi) {
		for (j in 0..min(ilo, i)) {
			if (a.data[i*a.stride+j] != aCopy.data[i*aCopy.stride+j]) {
				throw Exception("${prefix}: unexpected modification of A[${i},${j}]")
			}
		}
	}
	for (i in ihi+1..n) {
		for (j in 0..i) {
			if (a.data[i*a.stride+j] != aCopy.data[i*aCopy.stride+j]) {
				throw Exception("${prefix}: unexpected modification of A[${i},${j}]")
			}
		}
	}
	for (i in 0..=ilo) {
		for (j in i..ilo+1) {
			if (a.data[i*a.stride+j] != aCopy.data[i*aCopy.stride+j]) {
				throw Exception("${prefix}: unexpected modification at A[${i},${j}]")
			}
		}
		for (j in ihi+1..n) {
			if (a.data[i*a.stride+j] != aCopy.data[i*aCopy.stride+j]) {
				throw Exception("${prefix}: unexpected modification at A[${i},${j}]")
			}
		}
	}
	for (i in ihi+1..n) {
		for (j in i..n) {
			if (a.data[i*a.stride+j] != aCopy.data[i*aCopy.stride+j]) {
				throw Exception("${prefix}: unexpected modification at A[${i},${j}]")
			}
		}
	}

	for (i in 0..tau.size) {
        let v = tau[i]
		if (v.isNaN()) {
			throw Exception("${prefix}: unexpected NaN at tau[${i}]")
		}
	}

	let q = eye(n, n)
	if (ilo != ihi) {
		for (i in ilo+2..=ihi) {
			for (j in ilo+1..ihi) {
				q.data[i*q.stride+j] = a.data[i*a.stride+j-1]
			}
		}
		let nh = ihi - ilo
		lFloat64.dorgqr(nh, nh, nh, q.data[(ilo+1)*q.stride+ilo+1..], q.stride, tau[ilo..ihi], work, work.size)
	}
    let resid = residualOrthogonal(q, false)
	if (resid > tol) {
		throw Exception("${prefix}: Q is not orthogonal; resid=${resid}, want<=${tol}")
	}

	let aq = General64(
		rows:   n,
		cols:   n,
		stride: n,
		data:   Array<Float64>(n*n, item: 0.0)
    )
	gemm(NoTrans, NoTrans, 1.0, aCopy, q, 0.0, aq)
	let qaq = General64(
		rows:   n,
		cols:   n,
		stride: n,
		data:   Array<Float64>(n*n, item: 0.0)
    )
	gemm(Trans, NoTrans, 1.0, q, aq, 0.0, qaq)
	for (i in 0..=ilo) {
		for (j in ilo+1..=ihi) {
			let qaqij = qaq.data[i*qaq.stride+j]
			let diff = qaqij - a.data[i*a.stride+j]
			if (abs(diff) > tol) {
				throw Exception("${prefix}: Qᵀ*AOrig*Q and A are not equal, diff at [${i},${j}]=${diff}")
			}
		}
	}
	for (i in ilo+1..=ihi) {
		for (j in ilo..n) {
			let qaqij = qaq.data[i*qaq.stride+j]
			if (j < i-1) {
				if (abs(qaqij) > tol) {
					throw Exception("${prefix}: Qᵀ*AOrig*Q is not upper Hessenberg, [${i},${j}]=${qaqij}")
				}
				continue
			}
			let diff = qaqij - a.data[i*a.stride+j]
			if (abs(diff) > tol) {
				throw Exception("${prefix}: Qᵀ*AOrig*Q and A are not equal, diff at [${i},${j}]=${diff}")
			}
		}
	}
    return true
}