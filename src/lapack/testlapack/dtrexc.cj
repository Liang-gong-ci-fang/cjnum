package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dtrexcTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)

	for (n in [0, 1, 2, 3, 4, 5, 6, 10, 18, 31, 53]) {
		for (extra in [0, 3]) {
			for (_ in 0..100) {
				var ifst: Int64 = 0
				var ilst: Int64 = 0
				if (n > 0) {
					ifst = rnd.nextInt64(n)
					ilst = rnd.nextInt64(n)
				}
				if (!dtrexcTest(lFloat64, rnd, n, ifst, ilst, extra)) {
					return false
				}
			}
		}
	}
	return true
}

func dtrexcTest(lFloat64: LFloat64, rnd: Random, n: Int64, ifst: Int64, ilst: Int64, extra: Int64): Bool {
    const tol = 1e-13

    let (tmat, _, _) = randomSchurCanonical(n, n+extra, true, rnd)
    let tmatCopy = cloneGeneral(tmat)

    let (fstSize, fstFirst) = schurBlockSize(tmat, ifst)
    let (lstSize, lstFirst) = schurBlockSize(tmat, ilst)

    let name = "Case n=${n},ifst=${ifst},nbfst=${fstSize},ilst=${ilst},nblst=${lstSize},extra=${extra}"

    var compq = UpdateSchurNone

    var work = nanSlice(n)

    let (ifstGot, ilstGot, ok) = lFloat64.dtrexc(compq, n, tmat.data, tmat.stride, [], 1, ifst, ilst, work)

    if (!generalOutsideAllNaN(tmat)) {
        throw Exception("${name}: out-of-range write to T")
    }

    compq = UpdateSchur

    let tmat2 = cloneGeneral(tmatCopy)
    let q = eye(n, n+extra)
    let qCopy = cloneGeneral(q)
    work = nanSlice(n)

    let (ifstGot2, ilstGot2, ok2) = lFloat64.dtrexc(compq, n, tmat2.data, tmat2.stride, q.data, q.stride, ifst, ilst, work)

    if (!generalOutsideAllNaN(tmat2)) {
        throw Exception("${name}: out-of-range write to T2")
    }
    if (!generalOutsideAllNaN(q)) {
        throw Exception("${name}: out-of-range write to Q")
    }

    if (ifstGot != ifstGot2) {
        throw Exception("${name}: ifstGot != ifstGot2")
    }
    if (ilstGot != ilstGot2) {
        throw Exception("${name}: ilstGot != ilstGot2")
    }
    if (ok != ok2) {
        throw Exception("${name}: ok != ok2")
    }
    if (!equalGeneral(tmat, tmat2)) {
        throw Exception("${name}: T != T2")
    }

    var ifstWant = ifst
    if (!fstFirst) {
        ifstWant = ifst - 1
    }
    if (ifstWant != ifstGot) {
        throw Exception("${name}: unexpected ifst=${ifstGot}, want ${ifstWant}")
    }

    var ilstWant = ilst
    if (!lstFirst) {
        ilstWant--
    }
    if (ok) {
        if (ifstWant < ilstWant) {
            if (fstSize == 2 && lstSize == 1) {
                ilstWant--
            } else if (fstSize == 1 && lstSize == 2) {
                ilstWant++
            }
        }
        if (ilstWant != ilstGot) {
            throw Exception("${name}: unexpected ilst=${ilstGot}, want ${ilstWant}")
        }
    }

    if (n <= 1 || ifstGot == ilstGot) {
        if (!equalGeneral(tmat, tmatCopy)) {
            throw Exception("${name}: unexpected modification of T when no swapping")
        }
        if (!equalGeneral(q, qCopy)) {
            throw Exception("${name}: unexpected modification of Q when no swapping")
        }
        return true
    }

    if (!isSchurCanonicalGeneral(tmat)) {
        throw Exception("${name}: T is not in Schur canonical form")
    }

    let modMin = min(ifstGot, ilstGot)
    let modMax = max(ifstGot, ilstGot) + fstSize
    for (i in 0..n) {
        for (j in 0..n) {
            if (modMin <= i && i < modMax && j+1 >= i) {
                continue
            }
            if (modMin <= j && j < modMax && j+1 >= i) {
                continue
            }
            let diff = tmat.data[i*tmat.stride+j] - tmatCopy.data[i*tmatCopy.stride+j]
            if (diff != 0.0) {
                throw Exception("${name}: unexpected modification at T[${i},${j}]")
            }
        }
    }

    var resid = residualOrthogonal(q, false)
    if (resid > tol) {
        throw Exception("${name}: Q is not orthogonal; resid=${resid}, want<=${tol}")
    }

    for (i in 0..n) {
        for (j in 0..n) {
            if (modMin <= j && j < modMax) {
                continue
            }
            if (q.data[i*q.stride+j] != qCopy.data[i*qCopy.stride+j]) {
                throw Exception("${name}: unexpected modification of Q[${i},${j}]")
            }
        }
    }

    let qt = zeros(n, n, n)
    gemm(Trans, NoTrans, 1.0, q, tmatCopy, 0.0, qt)
    let qtq = cloneGeneral(tmat)
    gemm(NoTrans, NoTrans, -1.0, qt, q, 1.0, qtq)
    resid = dlange(MaxColumnSum, n, n, qtq.data, qtq.stride)
    if (resid > tol) {
        throw Exception("${name}: mismatch between Qᵀ*(initial T)*Q and (final T); resid=${resid}, want<=${tol}")
    }
	return true
}
