package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.floats.scalar.*
import std.math.*
import std.random.*

public func dlags2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (upper in [true, false]) {
        for (_ in 0..100) {
            let a1 = rnd.nextFloat64()
            let a2 = rnd.nextFloat64()
            let a3 = rnd.nextFloat64()
            let b1 = rnd.nextFloat64()
            let b2 = rnd.nextFloat64()
            let b3 = rnd.nextFloat64()

            let (csu, snu, csv, snv, csq, snq) = lFloat64.dlags2(upper, a1, a2, a3, b1, b2, b3)

            let detU = det2x2(csu, snu, -snu, csu)
            if (!equalWithinAbsOrRel(abs(detU), 1.0, 1e-14, 1e-14)) {
                throw Exception("U not orthogonal: det(U)=${detU}")
            }
            let detV = det2x2(csv, snv, -snv, csv)
            if (!equalWithinAbsOrRel(abs(detV), 1.0, 1e-14, 1e-14)) {
                throw Exception("V not orthogonal: det(V)=${detV}")
            }
            let detQ = det2x2(csq, snq, -snq, csq)
            if (!equalWithinAbsOrRel(abs(detQ), 1.0, 1e-14, 1e-14)) {
                throw Exception("Q not orthogonal: det(Q)=${detQ}")
            }

            let u = General64(
                rows:   2,
                cols:   2,
                stride: 2,
                data:   [csu, snu, -snu, csu],
            )
            let v = General64(
                rows:   2,
                cols:   2,
                stride: 2,
                data:   [csv, snv, -snv, csv],
            )
            let q = General64(
                rows:   2,
                cols:   2,
                stride: 2,
                data:   [csq, snq, -snq, csq],
            )

            let a = General64(rows: 2, cols: 2, stride: 2)
            let b = General64(rows: 2, cols: 2, stride: 2)
            if (upper) {
                a.data = [a1, a2, 0.0, a3]
                b.data = [b1, b2, 0.0, b3]
            } else {
                a.data = [a1, 0.0, a2, a3]
                b.data = [b1, 0.0, b2, b3]
            }

            let tmp = General64(rows: 2, cols: 2, stride: 2, data: Array<Float64>(4, item: 0.0))
            gemm(Trans, NoTrans, 1.0, u, a, 0.0, tmp)
            gemm(NoTrans, NoTrans, 1.0, tmp, q, 0.0, a)
            gemm(Trans, NoTrans, 1.0, v, b, 0.0, tmp)
            gemm(NoTrans, NoTrans, 1.0, tmp, q, 0.0, b)

            var gotA = 0.0
            var gotB = 0.0
            if (upper) {
                gotA = a.data[1]
                gotB = b.data[1]
            } else {
                gotA = a.data[2]
                gotB = b.data[2]
            }
            if (!equalWithinAbsOrRel(gotA, 0.0, 1e-14, 1e-14)) {
                throw Exception("unexpected non-zero value for zero triangle of Uᵀ*A*Q: ${gotA}")
            }
            if (!equalWithinAbsOrRel(gotB, 0.0, 1e-14, 1e-14)) {
                throw Exception("unexpected non-zero value for zero triangle of Vᵀ*B*Q: ${gotB}")
            }
        }
    }
    return true
}

func det2x2(a: Float64, b: Float64, c: Float64, d: Float64): Float64 {
    return a*d - b*c
}