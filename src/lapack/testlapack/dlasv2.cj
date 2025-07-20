package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

public func dlasv2Test(lFloat64: LFloat64): Bool {
	let rnd = Random(1)
	for (_ in 0..100) {
		let f = rnd.nextGaussianFloat64()
		let g = rnd.nextGaussianFloat64()
		let h = rnd.nextGaussianFloat64()

		let (ssmin, ssmax, snr, csr, snl, csl) = lFloat64.dlasv2(f, g, h)

		let tmp11 = csl * f
		let tmp12 = csl*g + snl*h
		let tmp21 = -snl * f
		let tmp22 = -snl*g + csl*h

		let ans11 = tmp11*csr + tmp12*snr
		let ans12 = tmp11*-snr + tmp12*csr
		let ans21 = tmp21*csr + tmp22*snr
		let ans22 = tmp21*-snr + tmp22*csr

		let lhs: Array<Float64> = [ans11, ans12, ans21, ans22]
		let rhs: Array<Float64> = [ssmax, 0.0, 0.0, ssmin]
		if (!fEqualApprox(rhs, lhs, 1e-12)) {
			throw Exception("SVD mismatch. f = ${f}, g = ${g}, h = ${h}.\nLHS: ${lhs}\nRHS: ${rhs}")
		}
	}
    return true
}
