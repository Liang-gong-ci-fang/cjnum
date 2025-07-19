package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.complex.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dlaqr04TestStruct {
    dlaqr04TestStruct(
        let h!: General64 = General64(),
        let ilo!: Int64 = 0, 
        let ihi!: Int64 = 0,
        let iloz!: Int64 = 0, 
        let ihiz!: Int64 = 0,
        var wantt!: Bool = false, 
        var wantz!: Bool = false,

        let evWant!: Array<Complex128> = []
    ) {}
	

	
}

public func dlaqr04Test(lFloat64: LFloat64): Bool {
	let rnd = Random(1)

    for (wantt in [true, false]) {
        for (wantz in [true, false]) {
            for (n in [1, 2, 3, 4, 5, 6, 10, 11, 12, 18, 29]) {
                for (extra in [0, 11]) {
                    for (recur in 0..=2) {
                        for (_ in 0..n) {
                            var ilo = rnd.nextInt64(n)
                            var ihi = rnd.nextInt64(n)
                            if (ilo > ihi) {
                                (ilo, ihi) = (ihi, ilo)
                            }
                            let iloz = rnd.nextInt64(ilo + 1)
                            let ihiz = ihi + rnd.nextInt64(n - ihi)
                            let h = randomHessenberg(n, n + extra, rnd)
                            if (ilo - 1 >= 0) {
                                h.data[ilo * h.stride + ilo - 1] = 0.0
                            }
                            if (ihi + 1 < n) {
                                h.data[(ihi + 1) * h.stride + ihi] = 0.0
                            }
                            let test = dlaqr04TestStruct(
                                h: h,
                                ilo: ilo,
                                ihi: ihi,
                                iloz: iloz,
                                ihiz: ihiz,
                                wantt: wantt,
                                wantz: wantz,
                            )
                            if (!testDlaqr04(lFloat64, test, false, recur)) {
								return false
							}
                            if (!testDlaqr04(lFloat64, test, true, recur)) {
								return false
							}
                        }
                    }
                }
            }
        }
    }

	for (n in [100, 500]) {
		for (_ in 0..5) {
			let h = randomHessenberg(n, n, rnd)
			let test = dlaqr04TestStruct(
				h:     h,
				ilo:   0,
				ihi:   n - 1,
				iloz:  0,
				ihiz:  n - 1,
				wantt: true,
				wantz: true,
            )
			if (!testDlaqr04(lFloat64, test, true, 1)) {
				return false
			}
		}
	}

	for (wantt in [true, false]) {
		for (wantz in [true, false]) {
			for (extra in [0, 1, 11]) {
				for (test in [
					dlaqr04TestStruct(
						h:    randomHessenberg(0, extra, rnd),
						ilo:  0,
						ihi:  -1,
						iloz: 0,
						ihiz: -1
                    ),
					dlaqr04TestStruct(
						h:    randomHessenberg(1, 1+extra, rnd),
						ilo:  0,
						ihi:  0,
						iloz: 0,
						ihiz: 0
					),
					dlaqr04TestStruct(
						h:    randomHessenberg(2, 2+extra, rnd),
						ilo:  1,
						ihi:  1,
						iloz: 1,
						ihiz: 1
					),
					dlaqr04TestStruct(
						h:    randomHessenberg(2, 2+extra, rnd),
						ilo:  0,
						ihi:  1,
						iloz: 0,
						ihiz: 1
					),
					dlaqr04TestStruct(
						h:    randomHessenberg(10, 10+extra, rnd),
						ilo:  0,
						ihi:  0,
						iloz: 0,
						ihiz: 0
					),
					dlaqr04TestStruct(
						h:    randomHessenberg(10, 10+extra, rnd),
						ilo:  0,
						ihi:  9,
						iloz: 0,
						ihiz: 9
					),
					dlaqr04TestStruct(
						h:    randomHessenberg(10, 10+extra, rnd),
						ilo:  0,
						ihi:  1,
						iloz: 0,
						ihiz: 1
					),
					dlaqr04TestStruct(
						h:    randomHessenberg(10, 10+extra, rnd),
						ilo:  0,
						ihi:  1,
						iloz: 0,
						ihiz: 9
                    ),
					dlaqr04TestStruct(
						h:    randomHessenberg(10, 10+extra, rnd),
						ilo:  9,
						ihi:  9,
						iloz: 0,
						ihiz: 9
                    )
                 ]) {
					if (test.ilo-1 >= 0) {
						test.h.data[test.ilo*test.h.stride+test.ilo-1] = 0.0
					}
					if (test.ihi+1 < test.h.rows) {
						test.h.data[(test.ihi+1)*test.h.stride+test.ihi] = 0.0
					}
					test.wantt = wantt
					test.wantz = wantz
					if (!testDlaqr04(lFloat64, test, false, 1)) {
						return false
					}
					if (!testDlaqr04(lFloat64, test, true, 1)) {
						return false
					}
				}
			}
		}
	}

	for (test in [
		dlaqr04TestStruct(
			h: General64(
				rows:   1,
				cols:   1,
				stride: 1,
				data:   Array<Float64>([7.09965484086874e-1])
			),
			ilo:    0,
			ihi:    0,
			iloz:   0,
			ihiz:   0,
			evWant: [Complex128(7.09965484086874e-1)]
		),
		dlaqr04TestStruct(
			h: General64(
				rows:   2,
				cols:   2,
				stride: 2,
				data: Array<Float64>([0.0, -1.0, 1.0, 0.0])
			),
			ilo:    0,
			ihi:    1,
			iloz:   0,
			ihiz:   1,
			evWant: [Complex128(0.0, 1.0), Complex128(0.0, -1.0)]
		),
		dlaqr04TestStruct(
			h: General64(
				rows:   2,
				cols:   2,
				stride: 2,
				data: Array<Float64>([
					6.25219991450918e-1, 8.17510791994361e-1,
					3.31218891622294e-1, 1.24103744878131e-1])
			),
			ilo:    0,
			ihi:    1,
			iloz:   0,
			ihiz:   1,
			evWant: [Complex128(9.52203547663447e-1), Complex128(-2.02879811334398e-1)]
		),
		dlaqr04TestStruct(
			h: General64(
				rows:   4,
				cols:   4,
				stride: 4,
				data: Array<Float64>([
					1.0, 0.0, 0.0, 0.0,
					0.0, 6.25219991450918e-1, 8.17510791994361e-1, 0.0,
					0.0, 3.31218891622294e-1, 1.24103744878131e-1, 0.0,
					0.0, 0.0, 0.0, 1.0])
			),
			ilo:    1,
			ihi:    2,
			iloz:   0,
			ihiz:   3,
			evWant: [Complex128(9.52203547663447e-1), Complex128(-2.02879811334398e-1)]
		),
		dlaqr04TestStruct(
			h: General64(
				rows:   2,
				cols:   2,
				stride: 2,
				data: Array<Float64>([
					-1.1219562276608, 6.85473513349362e-1,
					-8.19951061145131e-1, 1.93728523178888e-1])
			),
			ilo:  0,
			ihi:  1,
			iloz: 0,
			ihiz: 1,
			evWant: [
				Complex128(-4.64113852240958e-1, 3.59580510817350e-1),
				Complex128(-4.64113852240958e-1, -3.59580510817350e-1)]
        ),
		dlaqr04TestStruct(
			h: General64(
				rows:   5,
				cols:   5,
				stride: 5,
				data: Array<Float64>([
					9.57590178533658e-1, -5.10651295522708e-1, 9.24974510015869e-1, -1.30016306879522e-1, 2.92601986926954e-2,
					-1.08084756637964, 1.77529701001213, -1.36480197632509, 2.23196371219601e-1, 1.12912853063308e-1,
					0.0, -8.44075612174676e-1, 1.067867614486, -2.55782915176399e-1, -2.00598563137468e-1,
					0.0, 0.0, -5.67097237165410e-1, 2.07205057427341e-1, 6.54998340743380e-1,
					0.0, 0.0, 0.0, -1.89441413886041e-1, -4.18125416021786e-1])
            ),
			ilo:  0,
			ihi:  4,
			iloz: 0,
			ihiz: 4,
			evWant: [
				Complex128(2.94393309555622),
				Complex128(4.97029793606701e-1, 3.63041654992384e-1),
				Complex128(4.97029793606701e-1, -3.63041654992384e-1),
				Complex128(-1.74079119166145e-1, 2.01570009462092e-1),
				Complex128(-1.74079119166145e-1, -2.01570009462092e-1)]
        )
     ]) {
		test.wantt = true
		test.wantz = true
		if (!testDlaqr04(lFloat64, test, false, 1)) {
			return false
		}
		if (!testDlaqr04(lFloat64, test, true, 1)) {
			return false
		}
	}
	return true
}

func testDlaqr04(lFloat64: LFloat64, test: dlaqr04TestStruct, optwork: Bool, recur: Int64): Bool {
    const tol = 1e-14

    let h = cloneGeneral(test.h)
    let n = h.cols
    let extra = h.stride - h.cols
    let wantt = test.wantt
    let wantz = test.wantz
    let ilo = test.ilo
    let ihi = test.ihi
    let iloz = test.iloz
    let ihiz = test.ihiz

    var z = General64()
    var zCopy = General64()
    if (wantz) {
        z = eye(n, n + extra)
        zCopy = cloneGeneral(z)
    }

    let wr = nanSlice(ihi + 1)
    let wi = nanSlice(ihi + 1)

    var work = Array<Float64>([])
    if (optwork) {
        work = nanSlice(1)
        lFloat64.dlaqr04(wantt, wantz, n, ilo, ihi, h.data, h.stride, wr, wi, iloz, ihiz, z.data, max(1, z.stride), work, -1, recur)
        work = nanSlice(Int64(work[0]))
    } else {
        work = nanSlice(max(1, n))
    }

    let unconverged = lFloat64.dlaqr04(wantt, wantz, n, ilo, ihi, h.data, h.stride, wr, wi, iloz, ihiz, z.data, max(1, z.stride), work, work.size, recur)

    let prefix = "Case wantt=${wantt}, wantz=${wantz}, n=${n}, ilo=${ilo}, ihi=${ihi}, iloz=${iloz}, ihiz=${ihiz}, extra=${extra}, opt=${optwork}"

    if (!generalOutsideAllNaN(h)) {
        throw Exception("${prefix}: out-of-range write to H\n${h.data}")
    }
    if (!generalOutsideAllNaN(z)) {
        throw Exception("${prefix}: out-of-range write to Z\n${z.data}")
    }

    var start = ilo
    if (unconverged != 0) {
        start = unconverged
        if (start == ihi + 1) {
            throw Exception("${prefix}: no eigenvalue has converged")
        }
    }

    if (!isAllNaN(wr[0..start])) {
        throw Exception("${prefix}: unexpected modification of wr")
    }
    if (!isAllNaN(wi[0..start])) {
        throw Exception("${prefix}: unexpected modification of wi")
    }

    var hasReal = false
	var i = start
    while (i <= ihi) {
        if (wi[i] == 0.0) {
            hasReal = true

            if (wantt && wr[i] != h.data[i * h.stride + i]) {
                throw Exception("${prefix}: wr[${i}] != H[${i},${i}]")
            }
			i++
            continue
        }

        if (wr[i] != wr[i + 1]) {
            throw Exception("${prefix}: real part of conjugate pair not equal, i=${i}")
        }

        if (wi[i] < 0.0) {
            throw Exception("${prefix}: wi[${i}] not positive")
        }

        if (wi[i] != -wi[i + 1]) {
            throw Exception("${prefix}: wi[${i}] != wi[${i + 1}]")
        }
        if (wantt) {
            if (wr[i] != h.data[i * h.stride + i]) {
                throw Exception("${prefix}: wr[${i}] != H[${i},${i}]")
            }
            if (wr[i] != h.data[(i + 1) * h.stride + i + 1]) {
                throw Exception("${prefix}: wr[${i}] != H[${i + 1},${i + 1}]")
            }
            let im = sqrt(abs(h.data[(i + 1) * h.stride + i])) * sqrt(abs(h.data[i * h.stride + i + 1]))
            if (abs(im - wi[i]) > tol) {
                throw Exception("${prefix}: unexpected value of wi[${i}]: want ${im}, got ${wi[i]}")
            }
        }
		i += 2
    }

    if ((ihi + 1 - start) % 2 != 0 && !hasReal) {
        throw Exception("${prefix}: expected at least one real eigenvalue")
    }

    if (!test.evWant.isEmpty()) {
        for (i in start..=ihi) {
            let ev = Complex128(wr[i], wi[i])
            let (found, _) = containsComplex(test.evWant, ev, tol)
            if (!found) {
                throw Exception("${prefix}: unexpected eigenvalue ${ev}")
            }
        }
    }

    if (!wantz) {
        return true
    }

	let resid = residualOrthogonal(z, false)
    if (resid > tol * Float64(n)) {
        throw Exception("Case ${prefix}: Z is not orthogonal; resid=${resid}, want<=${tol * Float64(n)}")
    }

    for (i in 0..n) {
        for (j in 0..n) {
            if (iloz <= i && i <= ihiz && ilo <= j && j <= ihi) {
                continue
            }
            if (z.data[i * z.stride + j] != zCopy.data[i * zCopy.stride + j]) {
                throw Exception("${prefix}: Z modified outside of [iloz:ihiz+1,ilo:ihi+1] block")
            }
        }
    }
    if (wantt) {
        for (i in 2..n) {
            for (j in 0..i - 1) {
                h.data[i * h.stride + j] = 0.0
            }
        }
        let hz = eye(n, n)
        gemm(NoTrans, NoTrans, 1.0, test.h, z, 0.0, hz)
        let zhz = eye(n, n)
        gemm(Trans, NoTrans, 1.0, z, hz, 0.0, zhz)
        if (!equalApproxGeneral(zhz, h, 10.0 * tol)) {
            throw Exception("${prefix}: Zᵀ*(initial H)*Z and (final H) are not equal")
        }
    }
	return true
}
