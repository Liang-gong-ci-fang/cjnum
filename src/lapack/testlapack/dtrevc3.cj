package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import cjnum.util.*
import std.math.*
import std.random.*

public func dtrevc3Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)

    for (side in [EVRight, EVLeft, EVBoth]) {
        if (!runDtrevc3Test(lFloat64, rnd, side)) {
            return false
        }
    }
    return true
}

func runDtrevc3Test(lFloat64: LFloat64, rnd: Random, side: EVSide) {
    for (n in [0, 1, 2, 3, 4, 5, 6, 7, 10, 34]) {
        for (extra in [0, 11]) {
            for (optwork in [true, false]) {
                for (_ in 0..10) {
                    if (!dtrevc3Test(lFloat64, side, n, extra, optwork, rnd)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func dtrevc3Test(lFloat64: LFloat64, side: EVSide, n: Int64, extra: Int64, optwork: Bool, rnd: Random): Bool {
    const tol = 1e-15

    let name = "n=${n},extra=${extra},optwk=${optwork}"

    let right = side != EVLeft
    let left = side != EVRight

    let (tmat, wr, wi) = randomSchurCanonical(n, n + extra, true, rnd)
    let tmatCopy = cloneGeneral(tmat)

    var howmny = EVAll

    var vr: General64 = General64()
    var vl: General64 = General64()
    if (right) {
        vr = nanGeneral(n, n, n + extra)
    }
    if (left) {
        vl = nanGeneral(n, n, n + extra)
    }

    var work: Array<Float64>
    if (optwork) {
        work = [0.0]
        lFloat64.dtrevc3(side, howmny, [], n, tmat.data, tmat.stride, vl.data, max(1, vl.stride), vr.data,
            max(1, vr.stride), n, work, -1)
        work = Array<Float64>(Int64(work[0]), repeat: 0.0)
    } else {
        work = Array<Float64>(max(1, 3 * n), repeat: 0.0)
    }

    var mGot = lFloat64.dtrevc3(side, howmny, [], n, tmat.data, tmat.stride, vl.data, max(1, vl.stride), vr.data,
        max(1, vr.stride), n, work, work.size)

    if (!generalOutsideAllNaN(tmat)) {
        throw Exception("${name}: out-of-range write to T")
    }
    if (!equalGeneral(tmat, tmatCopy)) {
        throw Exception("${name}: unexpected modification of T")
    }
    if (!generalOutsideAllNaN(vr)) {
        throw Exception("${name}: out-of-range write to VR")
    }
    if (!generalOutsideAllNaN(vl)) {
        throw Exception("${name}: out-of-range write to VL")
    }

    var mWant = n
    if (mGot != mWant) {
        throw Exception("${name}: unexpected value of m=${mGot}, want ${mWant}")
    }

    if (right) {
        var resid = residualRightEV(tmat, vr, wr, wi)
        if (resid > tol) {
            throw Exception("${name}: unexpected right eigenvectors; residual=${resid}, want<=${tol}")
        }
        resid = residualEVNormalization(vr, wi)
        if (resid > tol) {
            throw Exception("${name}: unexpected normalization of right eigenvectors; residual=${resid}, want<=${tol}")
        }
    }
    if (left) {
        var resid = residualLeftEV(tmat, vl, wr, wi)
        if (resid > tol) {
            throw Exception("${name}: unexpected left eigenvectors; residual=${resid}, want<=${tol}")
        }
        resid = residualEVNormalization(vl, wi)
        if (resid > tol) {
            throw Exception("${name}: unexpected normalization of left eigenvectors; residual=${resid}, want<=${tol}")
        }
    }

    howmny = EVSelected

    let selected = Array<Bool>(n, repeat: true)
    let selectedWant = Array<Bool>(n, repeat: true)
    var nselr: Int64 = 0
    var nselc: Int64 = 0
    var j = n - 1
    while (j > 0) {
        if (wi[j] == 0.0) {
            if (nselr < max(1, n / 4)) {
                nselr++
                selected[j] = true
                selectedWant[j] = true
            }
            j--
        } else {
            if (nselc < max(1, n / 4)) {
                nselc++
                selected[j] = true
                selected[j - 1] = true
                selectedWant[j] = false
                selectedWant[j - 1] = true
            }
            j -= 2
        }
    }
    mWant = nselr + 2 * nselc

    var vrSel: General64 = General64()
    var vlSel: General64 = General64()
    if (right) {
        vrSel = nanGeneral(n, mWant, n + extra)
    }
    if (left) {
        vlSel = nanGeneral(n, mWant, n + extra)
    }

    if (optwork) {
        work = [0.0]
        lFloat64.dtrevc3(side, howmny, selected, n, tmat.data, tmat.stride, vlSel.data, max(1, vlSel.stride), vrSel.data,
            max(1, vrSel.stride), mWant, work, -1)
        work = Array<Float64>(Int64(work[0]), repeat: 0.0)
    }

    mGot = lFloat64.dtrevc3(side, howmny, selected, n, tmat.data, tmat.stride, vlSel.data, max(1, vlSel.stride),
        vrSel.data, max(1, vrSel.stride), mWant, work, work.size)

    if (!generalOutsideAllNaN(tmat)) {
        throw Exception("${name}: out-of-range write to T")
    }
    if (!equalGeneral(tmat, tmatCopy)) {
        throw Exception("${name}: unexpected modification of T")
    }
    if (!generalOutsideAllNaN(vrSel)) {
        throw Exception("${name}: out-of-range write to selected VR")
    }
    if (!generalOutsideAllNaN(vlSel)) {
        throw Exception("${name}: out-of-range write to selected VL")
    }

    if (mGot != mWant) {
        throw Exception("${name}: unexpected value of selected m=${mGot}, want ${mWant}")
    }

    for (i in 0..selected.size) {
        if (selected[i] != selectedWant[i]) {
            throw Exception("${name}: unexpected selected[${i}]")
        }
    }

    var k = 0
    if (right) {
        for (j in 0..n) {
            if (selected[j] && wi[j] == 0.0) {
                for (i in 0..n) {
                    if (vrSel.data[i * vrSel.stride + k] != vr.data[i * vr.stride + j]) {
                        throw Exception("${name}: unexpected selected VR")
                    }
                }
                k++
            } else if (selected[j] && wi[j] != 0.0) {
                for (i in 0..n) {
                    if (vrSel.data[i * vrSel.stride + k] != vr.data[i * vr.stride + j] ||
                        vrSel.data[i * vrSel.stride + k + 1] != vr.data[i * vr.stride + j + 1]) {
                        throw Exception("${name}: unexpected selected VR")
                    }
                }

                k += 2
            }
        }
    }

    k = 0
    if (left) {
        for (j in 0..n) {
            if (selected[j] && wi[j] == 0.0) {
                for (i in 0..n) {
                    if (vlSel.data[i * vlSel.stride + k] != vl.data[i * vl.stride + j]) {
                        throw Exception("${name}: unexpected selected VL")
                    }
                }
                k++
            } else if (selected[j] && wi[j] != 0.0) {
                for (i in 0..n) {
                    if (vlSel.data[i * vlSel.stride + k] != vl.data[i * vl.stride + j] ||
                        vlSel.data[i * vlSel.stride + k + 1] != vl.data[i * vl.stride + j + 1]) {
                        throw Exception("${name}: unexpected selected VL")
                    }
                }
                k += 2
            }
        }
    }

    howmny = EVAllMulQ

    var vrMul: General64 = General64()
    var qr: General64 = General64()
    var vlMul: General64 = General64()
    var ql: General64 = General64()
    if (right) {
        vrMul = randomGeneral(n, n, n + extra, rnd)
        qr = cloneGeneral(vrMul)
    }
    if (left) {
        vlMul = randomGeneral(n, n, n + extra, rnd)
        ql = cloneGeneral(vlMul)
    }

    if (optwork) {
        work = [0.0]
        lFloat64.dtrevc3(side, howmny, [], n, tmat.data, tmat.stride, vlMul.data, max(1, vlMul.stride), vrMul.data,
            max(1, vrMul.stride), n, work, -1)
        work = Array<Float64>(Int64(work[0]), repeat: 0.0)
    }

    mGot = lFloat64.dtrevc3(side, howmny, selected, n, tmat.data, tmat.stride, vlMul.data, max(1, vlMul.stride),
        vrMul.data, max(1, vrMul.stride), n, work, work.size)

    if (!generalOutsideAllNaN(tmat)) {
        throw Exception("${name}: out-of-range write to T")
    }
    if (!equalGeneral(tmat, tmatCopy)) {
        throw Exception("${name}: unexpected modification of T")
    }
    if (!generalOutsideAllNaN(vrMul)) {
        throw Exception("${name}: out-of-range write to VRMul")
    }
    if (!generalOutsideAllNaN(vlMul)) {
        throw Exception("${name}: out-of-range write to VLMul")
    }

    mWant = n
    if (mGot != mWant) {
        throw Exception("${name}: unexpected value of m=${mGot}, want ${mWant}")
    }

    if (right) {
        let qvWant = zeros(n, n, n)
        gemm(NoTrans, NoTrans, 1.0, qr, vr, 0.0, qvWant)
        normalizeEV(qvWant, wi)

        let r = zeros(n, n, n)
        for (i in 0..n) {
            for (j in 0..n) {
                r.data[i * r.stride + j] = vrMul.data[i * vrMul.stride + j] - qvWant.data[i * qvWant.stride + j]
            }
        }
        let qvNorm = dlange(MaxColumnSum, n, n, qvWant.data, qvWant.stride)
        let resid = dlange(MaxColumnSum, n, n, r.data, r.stride) / qvNorm / Float64(n)
        if (resid > tol) {
            throw Exception("${name}: unexpected VRMul; resid=${resid}, want <=${tol}")
        }
    }
    if (left) {
        let qvWant = zeros(n, n, n)
        gemm(NoTrans, NoTrans, 1.0, ql, vl, 0.0, qvWant)
        normalizeEV(qvWant, wi)

        let r = zeros(n, n, n)
        for (i in 0..n) {
            for (j in 0..n) {
                r.data[i * r.stride + j] = vlMul.data[i * vlMul.stride + j] - qvWant.data[i * qvWant.stride + j]
            }
        }
        let qvNorm = dlange(MaxColumnSum, n, n, qvWant.data, qvWant.stride)
        let resid = dlange(MaxColumnSum, n, n, r.data, r.stride) / qvNorm / Float64(n)
        if (resid > tol) {
            throw Exception("${name}: unexpected VLMul; resid=${resid}, want <=${tol}")
        }
    }
    return true
}

func residualEVNormalization(emat: General64, wi: Array<Float64>): Float64 {
    let n = emat.rows
    if (n == 0) {
        return 0.0
    }
    let e = emat.data
    let lde = emat.stride
    var enrmin = Float64.Inf
    var enrmax = 0.0
    var ipair = 0
    for (j in 0..n) {
        if (ipair == 0 && j < n - 1 && wi[j] != 0.0) {
            ipair = 1
        }
        var nrm = 0.0
        if (ipair == 0) {
            for (i in 0..n) {
                nrm = nanMax(nrm, abs(e[i * lde + j]))
            }
            enrmin = nanMin(enrmin, nrm)
            enrmax = nanMax(enrmax, nrm)
        } else if (ipair == 1) {
            for (i in 0..n) {
                nrm = nanMax(nrm, abs(e[i * lde + j]) + abs(e[i * lde + j + 1]))
            }
            enrmin = nanMin(enrmin, nrm)
            enrmax = nanMax(enrmax, nrm)
            ipair = 2
        } else if (ipair == 2) {
            ipair = 0
        }
    }
    return nanMax(abs(enrmin - 1.0), abs(enrmax - 1.0))
}

func normalizeEV(emat: General64, wi: Array<Float64>): Unit {
    let n = emat.rows
    if (n == 0) {
        return
    }
    let nFloat64 = nFloat64Implementation()
    let e = emat.data
    let lde = emat.stride
    var ipair = 0
    for (j in 0..n) {
        if (ipair == 0 && j < n - 1 && wi[j] != 0.0) {
            ipair = 1
        }
        if (ipair == 0) {
            let ii = nFloat64.idamax(n, e[j..], lde)
            let remax = 1.0 / abs(e[ii * lde + j])
            nFloat64.dscal(n, remax, e[j..], lde)
        } else if (ipair == 1) {
            var emax = 0.0
            for (i in 0..n) {
                emax = nanMax(emax, abs(e[i * lde + j]) + abs(e[i * lde + j + 1]))
            }
            nFloat64.dscal(n, 1.0 / emax, e[j..], lde)
            nFloat64.dscal(n, 1.0 / emax, e[j + 1..], lde)
            ipair = 2
        } else if (ipair == 2) {
            ipair = 0
        }
    }
}
