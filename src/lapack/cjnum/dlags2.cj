package cjnum.lapack.cjnum

import std.math.*

extend LFloat64Implementation {
    public func dlags2(upper: Bool, a1: Float64, a2: Float64, a3: Float64, b1: Float64, b2: Float64, b3: Float64): (Float64, 
        Float64, Float64, Float64, Float64, Float64) {
        var csu: Float64
        var snu: Float64
        var csv: Float64
        var snv: Float64
        var csq: Float64
        var snq: Float64

        if (upper) {
            let a = a1 * b3
            let d = a3 * b1
            let b = a2 * b1 - a1 * b2
            let (_, _, snr, csr, snl, csl) = dlasv2(a, b, d)
            if (abs(csl) >= abs(snl) || abs(csr) >= abs(snr)) {
                let ua11r = csl * a1
                let ua12 = csl * a2 + snl * a3
                let vb11r = csr * b1
                let vb12 = csr * b2 + snr * b3
                let aua12 = abs(csl) * abs(a2) + abs(snl) * abs(a3)
                let avb12 = abs(csr) * abs(b2) + abs(snr) * abs(b3)
                if (abs(ua11r) + abs(ua12) != 0.0) {
                    if (aua12 / (abs(ua11r) + abs(ua12)) <= avb12 / (abs(vb11r) + abs(vb12))) {
                        (csq, snq, _) = dlartg(-ua11r, ua12)
                    } else {
                        (csq, snq, _) = dlartg(-vb11r, vb12)
                    }
                } else {
                    (csq, snq, _) = dlartg(-vb11r, vb12)
                }
                csu = csl
                snu = -snl
                csv = csr
                snv = -snr
            } else {
                let ua21 = -snl * a1
                let ua22 = -snl * a2 + csl * a3
                let vb21 = -snr * b1
                let vb22 = -snr * b2 + csr * b3
                let aua22 = abs(snl) * abs(a2) + abs(csl) * abs(a3)
                let avb22 = abs(snr) * abs(b2) + abs(csr) * abs(b3)
                if (abs(ua21) + abs(ua22) != 0.0) {
                    if (aua22 / (abs(ua21) + abs(ua22)) <= avb22 / (abs(vb21) + abs(vb22))) {
                        (csq, snq, _) = dlartg(-ua21, ua22)
                    } else {
                        (csq, snq, _) = dlartg(-vb21, vb22)
                    }
                } else {
                    (csq, snq, _) = dlartg(-vb21, vb22)
                }
                csu = snl
                snu = csl
                csv = snr
                snv = csr
            }
        } else {
            let a = a1 * b3
            let d = a3 * b1
            let c = a2 * b3 - a3 * b2
            let (_, _, snr, csr, snl, csl) = dlasv2(a, c, d)
            if (abs(csr) >= abs(snr) || abs(csl) >= abs(snl)) {
                let ua21 = -snr * a1 + csr * a2
                let ua22r = csr * a3
                let vb21 = -snl * b1 + csl * b2
                let vb22r = csl * b3
                let aua21 = abs(snr) * abs(a1) + abs(csr) * abs(a2)
                let avb21 = abs(snl) * abs(b1) + abs(csl) * abs(b2)
                if ((abs(ua21) + abs(ua22r)) != 0.0) {
                    if (aua21 / (abs(ua21) + abs(ua22r)) <= avb21 / (abs(vb21) + abs(vb22r))) {
                        (csq, snq, _) = dlartg(ua22r, ua21)
                    } else {
                        (csq, snq, _) = dlartg(vb22r, vb21)
                    }
                } else {
                    (csq, snq, _) = dlartg(vb22r, vb21)
                }
                csu = csr
                snu = -snr
                csv = csl
                snv = -snl
            } else {
                let ua11 = csr * a1 + snr * a2
                let ua12 = snr * a3
                let vb11 = csl * b1 + snl * b2
                let vb12 = snl * b3
                let aua11 = abs(csr) * abs(a1) + abs(snr) * abs(a2)
                let avb11 = abs(csl) * abs(b1) + abs(snl) * abs(b2)
                if ((abs(ua11) + abs(ua12)) != 0.0) {
                    if (aua11 / (abs(ua11) + abs(ua12)) <= avb11 / (abs(vb11) + abs(vb12))) {
                        (csq, snq, _) = dlartg(ua12, ua11)
                    } else {
                        (csq, snq, _) = dlartg(vb12, vb11)
                    }
                } else {
                    (csq, snq, _) = dlartg(vb12, vb11)
                }
                csu = snr
                snu = csr
                csv = snl
                snv = csl
            }
        }
        return (csu, snu, csv, snv, csq, snq)
    }
}
