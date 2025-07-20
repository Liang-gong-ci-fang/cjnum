package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlaqr5Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [1, 2, 3, 4, 5, 6, 10, 30]) {
        for (extra in [0, 1, 20]) {
            for (kacc22 in [0, 1, 2]) {
                for (_ in 0..100) {
                    if (!testDlaqr5(lFloat64, n, extra, kacc22, rnd)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func testDlaqr5(lFloat64: LFloat64, n: Int64, extra: Int64, kacc22: Int64, rnd: Random): Bool {
    let tol = 1e-14

    let wantt = true
    let wantz = true
    let nshfts = 2 * n
    let sr = Array<Float64>(nshfts, item: 0.0)
    let si = Array<Float64>(nshfts, item: 0.0)
    var i = 0
    while (i < nshfts) {
        if (i == nshfts - 1 || rnd.nextFloat64() < 0.5) {
            let re = rnd.nextGaussianFloat64()
            (sr[i], si[i]) = (re, 0.0)
            i++
            continue
        }
        let re = rnd.nextGaussianFloat64()
        let im = rnd.nextGaussianFloat64()
        (sr[i], sr[i + 1]) = (re, re)
        (si[i], si[i + 1]) = (im, -im)
        i += 2
    }
    var ktop = rnd.nextInt64(n)
    var kbot = rnd.nextInt64(n)
    if (kbot < ktop) {
        (ktop, kbot) = (kbot, ktop)
    }

    let v = randomGeneral(nshfts / 2, 3, 3 + extra, rnd)
    let u = randomGeneral(2 * nshfts, 2 * nshfts, 2 * nshfts + extra, rnd)
    let nh = n
    let wh = randomGeneral(2 * nshfts, n, n + extra, rnd)
    let nv = n
    let wv = randomGeneral(n, 2 * nshfts, 2 * nshfts + extra, rnd)

    let h = randomHessenberg(n, n + extra, rnd)
    if (ktop > 0) {
        h.data[ktop * h.stride + ktop - 1] = 0.0
    }
    if (kbot < n - 1) {
        h.data[(kbot + 1) * h.stride + kbot] = 0.0
    }
    let hCopy = General64(rows: h.rows, cols: h.cols, stride: h.stride, data: Array<Float64>(h.data.size, item: 0.0))
    copy(hCopy.data, h.data)

    let z = eye(n, n + extra)

    lFloat64.dlaqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h.data, h.stride, 0, n - 1, z.data, z.stride,
        v.data, v.stride, u.data, u.stride, nv, wv.data, wv.stride, nh, wh.data, wh.stride)

    let prefix = "n=${n},lda=${extra},big=${kacc22}"

    if (!generalOutsideAllNaN(h)) {
        throw Exception("${prefix}: out-of-range write to H\n${h.data}")
    }
    if (!generalOutsideAllNaN(z)) {
        throw Exception("${prefix}: out-of-range write to Z\n${z.data}")
    }
    if (!generalOutsideAllNaN(u)) {
        throw Exception("${prefix}: out-of-range write to U\n${u.data}")
    }
    if (!generalOutsideAllNaN(v)) {
        throw Exception("${prefix}: out-of-range write to V\n${v.data}")
    }
    if (!generalOutsideAllNaN(wh)) {
        throw Exception("${prefix}: out-of-range write to WH\n${wh.data}")
    }
    if (!generalOutsideAllNaN(wv)) {
        throw Exception("${prefix}: out-of-range write to WV\n${wv.data}")
    }

    for (i in 0..n) {
        for (j in 0..i - 1) {
            if (h.data[i * h.stride + j] != 0.0) {
                throw Exception("${prefix}: H is not Hessenberg, H[${i},${j}]!=0")
            }
        }
    }

    var resid = residualOrthogonal(z, false)
    if (resid > tol * Float64(n)) {
        throw Exception(
            "Case ${prefix}: Z is not orthogonal; resid=${residualOrthogonal(z, false)}, want<=${tol * Float64(n)}")
    }

    let hz = zeros(n, n, n)
    gemm(NoTrans, NoTrans, 1.0, hCopy, z, 0.0, hz)
    let zhz = cloneGeneral(h)
    gemm(Trans, NoTrans, 1.0, z, hz, -1.0, zhz)
    resid = dlange(MaxColumnSum, n, n, zhz.data, zhz.stride)
    if (resid > tol * Float64(n)) {
        throw Exception("${prefix}: |Zᵀ*HOrig*Z - H|=${resid}, want<=${tol * Float64(n)}")
    }
    return true
}
