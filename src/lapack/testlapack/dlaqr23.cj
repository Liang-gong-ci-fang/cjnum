package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.complex.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dlaqr23TestStruct {
    dlaqr23TestStruct(
        var wantt!: Bool = false,
        var wantz!: Bool = false,
        let ktop!: Int64 = 0,
        let kbot!: Int64 = 0,
        var nw!: Int64 = 0,
        let h!: General64 = General64(),
        let iloz!: Int64 = 0,
        let ihiz!: Int64 = 0,
        let evWant!: Array<Complex128> = Array<Complex128>([])
    ) {}
}

func newDlaqr23TestCase(wantt: Bool, wantz: Bool, n: Int64, ldh: Int64, rnd: Random): dlaqr23TestStruct {
    var nw: Int64
    if (n <= 75) {
        nw = rnd.nextInt64(n) + 1
    } else {
        nw = 76 + rnd.nextInt64(n - 75)
    }

    let h = randomHessenberg(n, ldh, rnd)

    let ktop = rnd.nextInt64(n - nw + 1)
    var kbot = ktop + nw - 1
    kbot += rnd.nextInt64(n - kbot)

    if (ktop - 1 >= 0) {
        h.data[ktop * h.stride + ktop - 1] = 0.0
    }
    if (kbot + 1 < n) {
        h.data[(kbot + 1) * h.stride + kbot] = 0.0
    }

    let iloz = rnd.nextInt64(ktop + 1)
    let ihiz = kbot + rnd.nextInt64(n - kbot)
    return dlaqr23TestStruct(
        wantt: wantt,
        wantz: wantz,
        ktop: ktop,
        kbot: kbot,
        nw: nw,
        h: h,
        iloz: iloz,
        ihiz: ihiz,
    )
}

public func dlaqr23Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)

    for (wantt in [true, false]) {
        for (wantz in [true, false]) {
            for (n in [1, 2, 3, 4, 5, 6, 10, 18, 31, 100]) {
                for (extra in [0, 11]) {
                    for (_ in 0..10) {
                        let test = newDlaqr23TestCase(wantt, wantz, n, n + extra, rnd)
                        if (!testDlaqr23(lFloat64, test, false, 1, rnd)) {
                            return false
                        }
                        if (!testDlaqr23(lFloat64, test, true, 1, rnd)) {
                            return false
                        }
                        if (!testDlaqr23(lFloat64, test, false, 0, rnd)) {
                            return false
                        }
                        if (!testDlaqr23(lFloat64, test, true, 0, rnd)) {
                            return false
                        }
                    }
                }
            }
        }
    }

    for (wantt in [true, false]) {
        for (wantz in [true, false]) {
            for (extra in [0, 1, 11]) {
                let test = dlaqr23TestStruct(
                    wantt: wantt,
                    wantz: wantz,
                    h: randomHessenberg(0, extra, rnd),
                    ktop: 0,
                    kbot: -1,
                    iloz: 0,
                    ihiz: -1,
                    nw: 0,
                )
                if (!testDlaqr23(lFloat64, test, true, 1, rnd)) {
                    return false
                }
                if (!testDlaqr23(lFloat64, test, false, 1, rnd)) {
                    return false
                }
                if (!testDlaqr23(lFloat64, test, true, 0, rnd)) {
                    return false
                }
                if (!testDlaqr23(lFloat64, test, false, 0, rnd)) {
                    return false
                }
            }
        }
    }

    for (test in [
            dlaqr23TestStruct(
                h: General64(
                    rows: 1,
                    cols: 1,
                    stride: 1,
                    data: Array<Float64>([7.09965484086874e-1])
                ),
                ktop: 0,
                kbot: 0,
                iloz: 0,
                ihiz: 0,
                evWant: [Complex128(7.09965484086874e-1)]
            ),
            dlaqr23TestStruct(
                h: General64(
                    rows: 2,
                    cols: 2,
                    stride: 2,
                    data: Array<Float64>([0.0, -1.0, 1.0, 0.0])
                ),
                ktop: 0,
                kbot: 1,
                iloz: 0,
                ihiz: 1,
                evWant: [Complex128(0.0, 1.0), Complex128(0.0, -1.0)]
            ),
            dlaqr23TestStruct(
                h: General64(
                    rows: 2,
                    cols: 2,
                    stride: 2,
                    data: Array<Float64>(
                        [6.25219991450918e-1, 8.17510791994361e-1, 3.31218891622294e-1, 1.24103744878131e-1])
                ),
                ktop: 0,
                kbot: 1,
                iloz: 0,
                ihiz: 1,
                evWant: [Complex128(9.52203547663447e-1), Complex128(-2.02879811334398e-1)]
            ),
            dlaqr23TestStruct(
                h: General64(
                    rows: 4,
                    cols: 4,
                    stride: 4,
                    data: Array<Float64>(
                        [1.0, 0.0, 0.0, 0.0, 0.0, 6.25219991450918e-1, 8.17510791994361e-1, 0.0, 0.0,
                            3.31218891622294e-1, 1.24103744878131e-1, 0.0, 0.0, 0.0, 0.0, 1.0])
                ),
                ktop: 1,
                kbot: 2,
                iloz: 0,
                ihiz: 3,
                evWant: [Complex128(9.52203547663447e-1), Complex128(-2.02879811334398e-1)]
            ),
            dlaqr23TestStruct(
                h: General64(
                    rows: 2,
                    cols: 2,
                    stride: 2,
                    data: Array<Float64>(
                        [-1.1219562276608, 6.85473513349362e-1, -8.19951061145131e-1, 1.93728523178888e-1])
                ),
                ktop: 0,
                kbot: 1,
                iloz: 0,
                ihiz: 1,
                evWant: [Complex128(-4.64113852240958e-1, 3.59580510817350e-1),
                    Complex128(-4.64113852240958e-1, -3.59580510817350e-1)]
            ),
            dlaqr23TestStruct(
                h: General64(
                    rows: 5,
                    cols: 5,
                    stride: 5,
                    data: Array<Float64>(
                        [9.57590178533658e-1, -5.10651295522708e-1, 9.24974510015869e-1, -1.30016306879522e-1,
                            2.92601986926954e-2, -1.08084756637964, 1.77529701001213, -1.36480197632509,
                            2.23196371219601e-1, 1.12912853063308e-1, 0.0, -8.44075612174676e-1, 1.067867614486,
                            -2.55782915176399e-1, -2.00598563137468e-1, 0.0, 0.0, -5.67097237165410e-1,
                            2.07205057427341e-1, 6.54998340743380e-1, 0.0, 0.0, 0.0, -1.89441413886041e-1,
                            -4.18125416021786e-1])
                ),
                ktop: 0,
                kbot: 4,
                iloz: 0,
                ihiz: 4,
                evWant: [
                    Complex128(2.94393309555622),
                    Complex128(4.97029793606701e-1, 3.63041654992384e-1),
                    Complex128(4.97029793606701e-1, -3.63041654992384e-1),
                    Complex128(-1.74079119166145e-1, 2.01570009462092e-1),
                    Complex128(-1.74079119166145e-1, -2.01570009462092e-1)
                ]
            )
        ]) {
        test.wantt = true
        test.wantz = true
        test.nw = test.kbot - test.ktop + 1
        if (!testDlaqr23(lFloat64, test, true, 1, rnd)) {
            return false
        }
        if (!testDlaqr23(lFloat64, test, false, 1, rnd)) {
            return false
        }
        if (!testDlaqr23(lFloat64, test, true, 0, rnd)) {
            return false
        }
        if (!testDlaqr23(lFloat64, test, false, 0, rnd)) {
            return false
        }
    }
    return true
}

func testDlaqr23(lFloat64: LFloat64, test: dlaqr23TestStruct, opt: Bool, recur: Int64, rnd: Random): Bool {
    const tol = 1e-14

    let h = cloneGeneral(test.h)

    let n = h.cols
    let extra = h.stride - h.cols
    let wantt = test.wantt
    let wantz = test.wantz
    let ktop = test.ktop
    let kbot = test.kbot
    let nw = test.nw
    let iloz = test.iloz
    let ihiz = test.ihiz

    var z: General64 = General64()
    var zCopy: General64 = General64()
    if (wantz) {
        z = eye(n, n + extra)
        zCopy = cloneGeneral(z)
    }

    let sr = nanSlice(kbot + 1)
    let si = nanSlice(kbot + 1)

    let v = randomGeneral(nw, nw, nw + extra, rnd)
    var nh: Int = 0
    if (nw > 0) {
        nh = nw + rnd.nextInt64(nw)
    }
    let tmat = randomGeneral(nw, nh, nh + extra, rnd)
    var nv: Int = 0
    if (nw > 0) {
        nv = rnd.nextInt64(nw) + 1
    }
    let wv = randomGeneral(nv, nw, nw + extra, rnd)

    var work: Array<Float64> = Array<Float64>()
    if (opt) {
        work = nanSlice(1)
        lFloat64.dlaqr23(wantt, wantz, n, ktop, kbot, nw, h.data, h.stride, iloz, ihiz, z.data, max(1, z.stride), sr, si,
            v.data, v.stride, tmat.cols, tmat.data, tmat.stride, wv.rows, wv.data, wv.stride, work, -1, recur)
        work = nanSlice(Int64(work[0]))
    } else {
        work = nanSlice(max(1, 2 * nw))
    }

    let (ns, nd) = lFloat64.dlaqr23(wantt, wantz, n, ktop, kbot, nw, h.data, h.stride, iloz, ihiz, z.data,
        max(1, z.stride), sr, si, v.data, v.stride, tmat.cols, tmat.data, tmat.stride, wv.rows, wv.data, wv.stride, work,
        work.size, recur)

    let prefix = "Case wantt=${wantt}, wantz=${wantz}, n=${n}, ktop=${ktop}, kbot=${kbot}, nw=${nw}, iloz=${iloz}, ihiz=${ihiz}, extra=${extra}"

    if (!generalOutsideAllNaN(h)) {
        throw Exception("${prefix}: out-of-range write to H\n${h.data}")
    }
    if (!generalOutsideAllNaN(z)) {
        throw Exception("${prefix}: out-of-range write to Z\n${z.data}")
    }
    if (!generalOutsideAllNaN(v)) {
        throw Exception("${prefix}: out-of-range write to V\n${v.data}")
    }
    if (!generalOutsideAllNaN(tmat)) {
        throw Exception("${prefix}: out-of-range write to T\n${tmat.data}")
    }
    if (!generalOutsideAllNaN(wv)) {
        throw Exception("${prefix}: out-of-range write to WV\n${wv.data}")
    }
    if (!isAllNaN(sr[0..kbot - nd - ns + 1]) || !isAllNaN(sr[kbot + 1..])) {
        throw Exception("${prefix}: out-of-range write to sr")
    }
    if (!isAllNaN(si[0..kbot - nd - ns + 1]) || !isAllNaN(si[kbot + 1..])) {
        throw Exception("${prefix}: out-of-range write to si")
    }

    if (!isUpperHessenberg(h)) {
        throw Exception("${prefix}: H is not upper Hessenberg")
    }

    if (!test.evWant.isEmpty()) {
        for (i in kbot - nd + 1..=kbot) {
            let ev = Complex128(sr[i], si[i])
            let (found, _) = containsComplex(test.evWant, ev, tol)
            if (!found) {
                throw Exception("${prefix}: unexpected eigenvalue ${ev}")
            }
        }
    }

    if (!wantz) {
        return true
    }

    var zmod = false
    for (i in 0..n) {
        for (j in 0..n) {
            if (z.data[i * z.stride + j] == zCopy.data[i * zCopy.stride + j]) {
                continue
            }
            if (i < iloz || ihiz < i || j < kbot - nw + 1 || kbot < j) {
                zmod = true
            }
        }
    }
    if (zmod) {
        throw Exception("${prefix}: unexpected modification of Z")
    }
    let resid = residualOrthogonal(z, false)
    if (resid > tol * Float64(n)) {
        throw Exception("Case ${prefix}: Z is not orthogonal; resid=${resid}, want<=${tol*Float64(n)}")
    }
    if (wantt) {
        let hz = zeros(n, n, n)
        gemm(NoTrans, NoTrans, 1.0, test.h, z, 0.0, hz)
        let r = cloneGeneral(h)
        gemm(Trans, NoTrans, 1.0, z, hz, -1.0, r)
        let resid = dlange(MaxColumnSum, r.rows, r.cols, r.data, r.stride)
        if (resid > tol * Float64(n)) {
            throw Exception("${prefix}: |Zᵀ*(initial H)*Z - (final H)|=${resid}, want<=${tol*Float64(n)}")
        }
    }
    return true
}
