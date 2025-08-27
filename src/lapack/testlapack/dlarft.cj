package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.random.*

class dlarftTestStruct {
    dlarftTestStruct(
        let m: Int64, 
        let n: Int64, 
        let ldv: Int64, 
        let ldt: Int64
    ) {}
}

public func dlarftTest(lFloat64: LFloat64): Bool {
	let rnd = Random(1)
	for (store in [ColumnWise, RowWise]) {
		for (direct in [Forward, Backward]) {
			for (test in [
				dlarftTestStruct(6, 6, 0, 0),
				dlarftTestStruct(8, 6, 0, 0),
				dlarftTestStruct(6, 8, 0, 0),
				dlarftTestStruct(6, 6, 10, 15),
				dlarftTestStruct(8, 6, 10, 15),
				dlarftTestStruct(6, 8, 10, 15),
				dlarftTestStruct(6, 6, 15, 10),
				dlarftTestStruct(8, 6, 15, 10),
				dlarftTestStruct(6, 8, 15, 10)
			]) {

				let m = test.m
				let n = test.n
				var lda = n
				if (lda == 0) {
					lda = n
				}

				let a = Array<Float64>(m * lda, repeat: 0.0)
				for (i in 0..m) {
					for (j in 0..lda) {
						a[i*lda+j] = rnd.nextFloat64()
					}
				}

				let k = min(m, n)
				let tau = Array<Float64>(k, repeat: 0.0)
				let work = Array<Float64>(n, repeat: 0.0)
				lFloat64.dgeqr2(m, n, a, lda, tau, work)

				let vMatTmp = extractVMat(m, n, a, lda, Forward, ColumnWise)
				let vMat = constructVMat(vMatTmp, store, direct)
				let v = vMat.data
				let ldv = vMat.stride

				let h = constructH(tau, vMat, store, direct)

				var ldt = test.ldt
				if (ldt == 0) {
					ldt = k
				}

				let tm = Array<Float64>(k*ldt, repeat: 0.0)
				for (i in 0..tm.size) {
					tm[i] = 100.0 + rnd.nextFloat64()
				}

				lFloat64.dlarft(direct, store, m, k, v, ldv, tau, tm, ldt)

				let tData = Array<Float64>(tm.size, repeat: 0.0)
				copy(tData, tm)
				if (direct == Forward) {

					for (i in 0..k) {
						for (j in 0..i) {
							tData[i*ldt+j] = 0.0
						}
					}
				} else {

					for (i in 0..k) {
						for (j in i + 1..k) {
							tData[i*ldt+j] = 0.0
						}
					}
				}

				let T = General64(
					rows:   k,
					cols:   k,
					stride: ldt,
					data:   tData
                )

				let vMatT = General64(
					rows:   vMat.cols,
					cols:   vMat.rows,
					stride: vMat.rows,
					data:   Array<Float64>(vMat.cols*vMat.rows, repeat: 0.0)
                )
				for (i in 0..vMat.rows) {
					for (j in 0..vMat.cols) {
						vMatT.data[j*vMatT.stride+i] = vMat.data[i*vMat.stride+j]
					}
				}
				var comp: General64 = General64()
				if (store == ColumnWise) {

					let tmp = General64(
						rows:   T.rows,
						cols:   vMatT.cols,
						stride: vMatT.cols,
						data:   Array<Float64>(T.rows*vMatT.cols, repeat: 0.0)
                    )

					gemm(NoTrans, NoTrans, 1.0, T, vMatT, 0.0, tmp)
					comp = General64(
						rows:   vMat.rows,
						cols:   tmp.cols,
						stride: tmp.cols,
						data:   Array<Float64>(vMat.rows*tmp.cols, repeat: 0.0)
                    )

					gemm(NoTrans, NoTrans, 1.0, vMat, tmp, 0.0, comp)
				} else {

					let tmp = General64(
						rows:   T.rows,
						cols:   vMat.cols,
						stride: vMat.cols,
						data:   Array<Float64>(T.rows*vMat.cols, repeat: 0.0),
                    )

					gemm(NoTrans, NoTrans, 1.0, T, vMat, 0.0, tmp)
					comp = General64(
						rows:   vMatT.rows,
						cols:   tmp.cols,
						stride: tmp.cols,
						data:   Array<Float64>(vMatT.rows*tmp.cols, repeat: 0.0)
                    )

					gemm(NoTrans, NoTrans, 1.0, vMatT, tmp, 0.0, comp)
				}

				for (i in 0..comp.rows) {
					for (j in 0..comp.cols) {
						comp.data[i*m+j] *= -1.0
						if (i == j) {
							comp.data[i*m+j] += 1.0
						}
					}
				}
				if (!fEqualApprox(comp.data, h.data, 1e-14)) {
					throw Exception("T does not construct proper H. Store = ${store}, Direct = ${direct}.\nWant ${h.data}\ngot ${comp.data}.")
				}
			}
		}
	}
	return true
}
