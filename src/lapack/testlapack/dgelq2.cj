package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dgelq2TestStruct {
    dgelq2TestStruct(
        let m: Int64,
        let n: Int64,
        let lda: Int64
    ) {}
}

public func dgelq2Test(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    var c = 0
    for (test in [
            dgelq2TestStruct(1, 1, 0),
            dgelq2TestStruct(2, 2, 0),
            dgelq2TestStruct(3, 2, 0),
            dgelq2TestStruct(2, 3, 0),
            dgelq2TestStruct(1, 12, 0),
            dgelq2TestStruct(2, 6, 0),
            dgelq2TestStruct(3, 4, 0),
            dgelq2TestStruct(4, 3, 0),
            dgelq2TestStruct(6, 2, 0),
            dgelq2TestStruct(1, 12, 0),
            dgelq2TestStruct(1, 1, 20),
            dgelq2TestStruct(2, 2, 20),
            dgelq2TestStruct(3, 2, 20),
            dgelq2TestStruct(2, 3, 20),
            dgelq2TestStruct(1, 12, 20),
            dgelq2TestStruct(2, 6, 20),
            dgelq2TestStruct(3, 4, 20),
            dgelq2TestStruct(4, 3, 20),
            dgelq2TestStruct(6, 2, 20),
            dgelq2TestStruct(1, 12, 20)
        ]) {
        let n = test.n
        let m = test.m
        var lda = test.lda
        if (lda == 0) {
            lda = test.n
        }
        let k = min(m, n)
        let tau = Array<Float64>(k, {_ => rnd.nextFloat64()})

        let work = Array<Float64>(m, {_ => rnd.nextFloat64()})

        let a = Array<Float64>(m * lda, {_ => rnd.nextFloat64()})

        let aCopy = Array<Float64>(a.size, repeat: 0.0)
        copy(aCopy, a)
        lFloat64.dgelq2(m, n, a, lda, tau, work)

        let Q = constructQ("LQ", m, n, a, lda, tau)

        let resid = residualOrthogonal(Q, false)
        if (resid > tol) {
            throw Exception("Case ${c}: Q not orthogonal; resid=${resid}, want<=${tol}")
        }

        let L = General64(
            rows: m,
            cols: n,
            stride: n,
            data: Array<Float64>(m * n, repeat: 0.0)
        )
        for (i in 0..m) {
            for (j in 0..=min(i, n - 1)) {
                L.data[i * L.stride + j] = a[i * lda + j]
            }
        }

        let ans = General64(
            rows: m,
            cols: n,
            stride: lda,
            data: Array<Float64>(m * lda, repeat: 0.0)
        )
        copy(ans.data, aCopy)
        gemm(NoTrans, NoTrans, 1.0, L, Q, 0.0, ans)
        if (!fEqualApprox(aCopy, ans.data, tol)) {
            throw Exception("Case ${c}, LQ mismatch. Want ${aCopy}, got ${ans.data}.")
        }
        c++
    }
    return true
}
