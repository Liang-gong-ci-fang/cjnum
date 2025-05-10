package cjnum.lapack.testlapack

import std.math.*
import std.random.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.blas.*
import cjnum.blas.blas64.*

public func dlacn2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [1, 2, 3, 4, 5, 7, 10, 15, 20, 100]) {
        for (cas in 0..10) {
            let a = randomGeneral(n, n, n, rnd)

            var norm1 = 0.0
            for ( j in 0..n) {
                var sum = 0.0
                for (i in 0..n) {
                    sum += abs(a.data[i * a.stride + j])
                }
                if (sum > norm1) {
                    norm1 = sum
                }
            }

            let x = Array<Float64>(n, item: 0.0)
            let work = Array<Float64>(n, item: 0.0)
            let v = Array<Float64>(n, item: 0.0)
            let isgn = Array<Int64>(n, item: 0)
            var kase = 0
            var isave = Array<Int64>(3, item: 0)
            var got = 0.0

            while (true) {
                (got, kase) = lFloat64.dlacn2(n, v, x, isgn, got, kase, isave)
                if (kase == 0) {
                    break
                } else if (kase == 1) {
                    gemv(NoTrans, 1.0, a, Vector(data: x, inc: 1), 0.0, Vector(data: work, inc: 1))
                    copy(x, work)
                } else if (kase == 2) {
                    gemv(Trans, 1.0, a, Vector(data: x, inc: 1), 0.0, Vector(data: work, inc: 1))
                    copy(x, work)
                } else {
                    throw Exception("Dlacn2 returned invalid value of kase")
                }
            }

            if (abs(got - norm1) > 1e-8 && got > norm1) {
                throw Exception("Case n=${n}: not lower estimate. 1-norm ${norm1}, estimate ${got}")
            }
        }
    }
    return true
}