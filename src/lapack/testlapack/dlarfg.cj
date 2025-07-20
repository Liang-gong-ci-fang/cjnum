package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64

class DlarfgTestStruct {
    DlarfgTestStruct(
        let alpha!: Float64,
        let n!: Int64,
        let x!: Array<Float64> = []
    ) {}
}

public func dlarfgTest(lFloat64: LFloat64): Bool {
    const tol = 1e-14
    let rnd = Random(1)
    let tests = [
        DlarfgTestStruct(
            alpha: 4.0,
            n: 3
        ),
        DlarfgTestStruct(
            alpha: -2.0,
            n: 3
        ),
        DlarfgTestStruct(
            alpha: 0.0,
            n: 3
        ),
        DlarfgTestStruct(
            alpha: 1.0,
            n: 1
        ),
        DlarfgTestStruct(
            alpha: 1.0,
            n: 4,
            x: [4.0, 5.0, 6.0]
        ),
        DlarfgTestStruct(
            alpha: 1.0,
            n: 4,
            x: [0.0, 0.0, 0.0]
        ),
        DlarfgTestStruct(
            alpha: dlamchS,
            n: 4,
            x: [dlamchS, dlamchS, dlamchS]
        )
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let n = test.n
        let incX = 1
        var x: Array<Float64> = []
        if (test.x.size == 0) {
            x = Array<Float64>(n - 1, { _ => rnd.nextFloat64() })
        } else {
            if (test.x.size != n - 1) {
                throw Exception("bad test")
            }
            x = Array<Float64>(n - 1, repeat: 0.0)
            copy(x, test.x)
        }
        let xcopy = x.clone()
        let alpha = test.alpha
        let (beta, tau) = lFloat64.dlarfg(n, alpha, x, incX)

        let h = Array<Float64>(n * n, repeat: 0.0)
        for (i in 0..n) {
            h[i * n + i] = 1.0
        }
        let hmat = General64(
            rows: n,
            cols: n,
            stride: n,
            data: h
        )
        let v = Array<Float64>(n, repeat: 0.0)
        copy(v[1..], x)
        v[0] = 1.0
        let vVec = Vector(
            inc: 1,
            data: v
        )
        ger(-tau, vVec, vVec, hmat)
        let eye = General64(
            rows: n,
            cols: n,
            stride: n,
            data: Array<Float64>(n * n, repeat: 0.0)
        )
        gemm(Trans, NoTrans, 1.0, hmat, hmat, 0.0, eye)
        let dist = distFromIdentity(n, eye.data, n)
        if (dist > tol) {
            throw Exception("Hᵀ * H is not close to I, dist=${dist}")
        } 

        let xVec = Vector(
            inc: 1,
            data: Array<Float64>(n, repeat: 0.0)
        )
        xVec.data[0] = test.alpha
        copy(xVec.data[1..], xcopy)

        let ans = Array<Float64>(n, repeat: 0.0)
        let ansVec = Vector(
            inc: 1,
            data: ans
        )
        gemv(NoTrans, 1.0, hmat, xVec, 0.0, ansVec)
        if (abs(ans[0] - beta) > tol) {
            throw Exception("Case ${i}, beta mismatch. Want ${ans[0]}, got ${beta}")
        }
        if (fNorm(ans[1..n], Float64.Inf) > tol) {
            throw Exception("Case ${i}, nonzero answer ${ans[1..n]}")
        }
    }
    return true
}