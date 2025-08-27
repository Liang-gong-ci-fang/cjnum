package cjnum.lapack.testlapack

import std.random.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.floats.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64

class Dgetf2TestStruct {
    Dgetf2TestStruct(
        let m!: Int64,
        let n!: Int64,
        let lda!: Int64,
        let a!: Array<Float64>
    ) {}
}

public func dgetf2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (test in [
        [10, 10, 0],
		[10, 5, 0],
		[10, 5, 0],
		[10, 10, 20],
		[5, 10, 20],
		[10, 5, 20]
    ]) {
        let m = test[0]
        let n = test[1]
        let lda = if (test[2] == 0) {
            n
        } else {
            test[2]
        }
        let a = Array<Float64>(m * lda, { _ => rnd.nextFloat64() })
        let aCopy = a.clone()
        let mn = min(m, n)
        let ipiv = Array<Int64>(mn, { _ => rnd.nextInt64() })

        let ok = lFloat64.dgetf2(m, n, a, lda, ipiv)
        checkPLU(ok, m, n, lda, ipiv, a, aCopy, 1e-14, true)
    }

    let tests = [
        Dgetf2TestStruct(
            m: 2,
            n: 2,
            lda: 2,
            a: [1.0, 0.0, 0.0, 0.0]
        ),
        Dgetf2TestStruct(
            m: 2,
            n: 2,
            lda: 2,
            a: [1.0, 5.0, 2.0, 10.0]
        ),
        Dgetf2TestStruct(
            m: 3,
            n: 3,
            lda: 3,
            a: [1.0, 5.0, 7.0, 2.0, 10.0, -3.0, 5.0, 25.0, 1.0]
        ),
        Dgetf2TestStruct(
            m: 3,
            n: 4,
            lda: 4,
            a: [1.0, 5.0, 7.0, 9.0, 2.0, 10.0, -3.0, 11.0, 5.0, 25.0, 1.0, 31.0]
        )
    ]
    for (n in 0..tests.size) {
        let test = tests[n]
        if (lFloat64.dgetf2(test.m, test.n, test.a, test.lda, Array<Int64>(min(test.m, test.n), repeat: 0))) {
            throw Exception("${n}: Returned ok with singular matrix.")
        }
    }
    return true
}

func checkPLU(ok: Bool, m: Int64, n: Int64, lda: Int64, ipiv: Array<Int64>, factorized: Array<Float64>, original: Array<Float64>, tol: Float64, print: Bool): Unit {
    var hasZeroDiagonal = false
    for (i in 0..min(m, n)) {
        if (factorized[i * lda + i] == 0.0) {
            hasZeroDiagonal = true
            break
        }
    }
    if (hasZeroDiagonal && ok) {
        throw Exception("Has a zero diagonal but returned ok")
    }
    if (!hasZeroDiagonal && !ok) {
        throw Exception("Non-zero diagonal but returned !ok")
    }

    let mn = min(m, n)
    let l = Array<Float64>(m * mn, repeat: 0.0)
    let ldl = mn
    let u = Array<Float64>(n * mn, repeat: 0.0)
    let ldu = n
    for (i in 0..m) {
        for (j in 0..n) {
            let v = factorized[i * lda + j]
            if (i == j) {
                l[i * ldl + i] = 1.0
                u[i * ldu + i] = v
            } else if (i > j) {
                l[i * ldl + j] = v
            } else if (i < j) {
                u[i * ldu + j] = v
            }
        }
    }

    let LU = General64(
        rows: m,
        cols: n,
        stride: n,
        data: Array<Float64>(m * n, repeat: 0.0)
    )
    let U = General64(
        rows: mn,
        cols: n,
        stride: ldu,
        data: u
    )
    let L = General64(
        rows: m,
        cols: mn,
        stride: ldl,
        data: l
    )
    gemm(NoTrans, NoTrans, 1.0, L, U, 0.0, LU)

    let p = Array<Float64>(m * m, repeat: 0.0) // a
    let ldp = m
    for (i in 0..m) {
        p[i * ldp + i] = 1.0
    }
    var i = ipiv.size - 1
    while (i >= 0) {
        let v = ipiv[i]
        swap(Vector(n: m, inc: 1, data: p[i * ldp..]), Vector(n: m, inc: 1, data: p[v * ldp..]))
        i--
    }
    let P = General64(
        rows: m,
        cols: m,
        stride: m,
        data: p
    )
    let aComp = General64(
        rows: m,
        cols: n,
        stride: lda,
        data: Array<Float64>(m * lda, repeat: 0.0)
    )
    copy(aComp.data, factorized)
    gemm(NoTrans, NoTrans, 1.0, P, LU, 0.0, aComp)
    if (!fEqualApprox(aComp.data, original, tol)) {
        if (print) {
            throw Exception("PLU multiplication does not match original matrix.\nWant: ${original}\nGot: ${aComp.data}")
        }
        throw Exception("PLU multiplication does not match original matrix.")
    }
}