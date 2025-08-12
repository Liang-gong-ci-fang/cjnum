package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.blas.blas64.*

public func dgetrfTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (test in [
		[10, 5, 0],
		[5, 10, 0],
		[10, 10, 0],
		[300, 5, 0],
		[3, 500, 0],
		[4, 5, 0],
		[300, 200, 0],
		[204, 300, 0],
		[1, 3000, 0],
		[3000, 1, 0],
		[10, 5, 20],
		[5, 10, 20],
		[10, 10, 20],
		[300, 5, 400],
		[3, 500, 600],
		[200, 200, 300],
		[300, 200, 300],
		[204, 300, 400],
		[1, 3000, 4000],
		[3000, 1, 4000]
    ]) {
        let m = test[0]
        let n = test[1]
        let lda = if (test[2] == 0) {
            n
        } else {
            test[2]
        }
        let a = Array<Float64>(m * lda, { _ => rnd.nextFloat64() })
        let mn = min(m, n)
        let ipiv = Array<Int64>(mn, { _ => rnd.nextInt64(Int64.Max) })

        let aCopy = a.clone()
        let ok = lFloat64.dgetrf(m, n, a, lda, ipiv)
        checkPLU(ok, m, n, lda, ipiv, a, aCopy, 1e-10, false)
    }
    return true
}
