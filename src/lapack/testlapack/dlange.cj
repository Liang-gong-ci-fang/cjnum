package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.lapack.*
import cjnum.blas.blas64.*

class dlangeTestStruct {
    dlangeTestStruct(
        let m: Int64,
        let n: Int64,
        let lda: Int64
    ) { }
}

public func dlangeTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    let tests = [
        dlangeTestStruct(4, 3, 0),
        dlangeTestStruct(3, 4, 0),
        dlangeTestStruct(4, 3, 100),
        dlangeTestStruct(3, 4, 100)
    ]
    for (i in 0..tests.size) {
        let test = tests[i]
        let m = test.m
        let n = test.n
        let lda = if (test.lda == 0) {
            n
        } else {
            test.lda
        }
        let a = Array<Float64>(m * lda, { _ => 
            rnd.nextFloat64() - 0.5
        })
        let aCopy = a.clone()
        let work = Array<Float64>(n, { _ => 
            rnd.nextFloat64()
        })

        var norm = lFloat64.dlange(MaxAbs, m, n, a, lda, work)
        var ans = 0.0
        for (i in 0..m) {
            let idx = iamax(Vector(n: n, inc: 1, data: aCopy[i * lda..]))
            ans = max(ans, abs(a[i * lda + idx]))
        }

        if (ans != norm) {
            throw Exception("MaxAbs mismatch. Want ${ans}, got ${norm}.")
        }

        norm = lFloat64.dlange(MaxColumnSum, m, n, a, lda, work)
        ans = 0.0
        for (i in 0..n) {
            let sum = asum(Vector(n: m, inc: lda, data: aCopy[i..]))
            ans = max(ans, sum)
        }
        if (abs(norm - ans) > 1e-14) {
            throw Exception("MaxColumnSum mismatch. Want ${ans}, got ${norm}.")
        }

        norm = lFloat64.dlange(MaxRowSum, m, n, a, lda, work)
        ans = 0.0
        for (i in 0..m) {
            let sum = asum(Vector(n: n, inc: 1, data: aCopy[i * lda..]))
            ans = max(ans, sum)
        }
        if (abs(norm - ans) > 1e-14) {
            throw Exception("MaxRowSum mismatch. Want ${ans}, got ${norm}.")
        }

        norm = lFloat64.dlange(Frobenius, m, n, a, lda, work)
        ans = 0.0
        for (i in 0..m) {
            let sum = nrm2(Vector(n: n, inc: 1, data: aCopy[i * lda..]))
            ans += sum * sum
        }
        ans = sqrt(ans)
        if (abs(norm - ans) > 1e-14) {
            throw Exception("Frobenius norm mismatch. Want ${ans}, got ${norm}.")
        }
    }
    return true
}