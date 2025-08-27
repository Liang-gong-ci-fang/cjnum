package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.random.*

public func dgebakTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)

    for (job in [BalanceNone, Permute, Scale, PermuteScale]) {
        for (side in [EVLeft, EVRight]) {
            for (n in [0, 1, 2, 3, 4, 5, 6, 10, 18, 31, 53]) {
                for (extra in [0, 11]) {
                    for (_ in 0..100) {
                        let m = rnd.nextInt64(n + 1)
                        let v = randomGeneral(n, m, m + extra, rnd)
                        var ilo: Int64 = 0
                        var ihi: Int64 = 0
                        if (v.rows > 0) {
                            ihi = rnd.nextInt64(n)
                            ilo = rnd.nextInt64(ihi + 1)
                        } else {
                            ihi = -1
                        }
                        if (!testDgebak(lFloat64, job, side, ilo, ihi, v, rnd)) {
                            return false
                        }
                    }
                }
            }
        }
    }
    return true
}

func testDgebak(lFloat64: LFloat64, job: BalanceJob, side: EVSide, ilo: Int64, ihi: Int64, v: General64, rnd: Random): Bool {
    const tol = 1e-15
    let n = v.rows
    let m = v.cols
    let extra = v.stride - v.cols

    let d = eye(n, n)
    let dinv = eye(n, n)
    var scale = nanSlice(n)
    if (job == Scale || job == PermuteScale) {
        if (ilo == ihi) {
            scale[ilo] = 1.0
        } else {
            for (i in ilo..=ihi) {
                scale[i] = 2.0 * rnd.nextFloat64()
                d.data[i * d.stride + i] = scale[i]
                dinv.data[i * dinv.stride + i] = 1.0 / scale[i]
            }
        }
    }

    let p = eye(n, n)
    if (job == Permute || job == PermuteScale) {
        for (i in n - 1..ihi: -1) {
            scale[i] = Float64(rnd.nextInt64(i + 1))
            swap(Vector(n: n, data: p.data[i..], inc: p.stride),
                Vector(n: n, data: p.data[Int64(scale[i])..], inc: p.stride))
        }
        for (i in 0..ilo) {
            scale[i] = Float64(i + rnd.nextInt64(ihi - i + 1))
            swap(Vector(n: n, data: p.data[i..], inc: p.stride),
                Vector(n: n, data: p.data[Int64(scale[i])..], inc: p.stride))
        }
    }

    let got = cloneGeneral(v)
    lFloat64.dgebak(job, side, n, ilo, ihi, scale, m, got.data, got.stride)

    let prefix = "job=${job},side=${side},n=${n},ilo=${ilo},ihi=${ihi},m=${m},extra=${extra}"

    if (!generalOutsideAllNaN(got)) {
        throw Exception("${prefix}: out-of-range write to V\n${got.data}")
    }

    let dv = zeros(n, m, m)
    if (side == EVRight) {
        gemm(NoTrans, NoTrans, 1.0, d, v, 0.0, dv)
    } else {
        gemm(NoTrans, NoTrans, 1.0, dinv, v, 0.0, dv)
    }
    let want = zeros(n, m, m)
    gemm(NoTrans, NoTrans, 1.0, p, dv, 0.0, want)

    if (!equalApproxGeneral(want, got, tol)) {
        throw Exception("${prefix}: unexpected value of V")
    }
    return true
}
