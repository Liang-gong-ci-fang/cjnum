package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

public func dlapmrTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (fwd in [true, false]) {
        for (m in [0, 1, 2, 3, 4, 5, 10]) {
            for (n in [0, 1, 4]) {
                for (ldx in [max(1, n), n + 3]) {
                    if (!dlapmrTest(lFloat64, rnd, fwd, m, n, ldx)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func dlapmrTest(lFloat64: LFloat64, rnd: Random, fwd: Bool, m: Int64, n: Int64, ldx: Int64): Bool {
    let name = "forwrd=${fwd},m=${m},n=${n},ldx=${ldx}"

    let nFloat64 = nFloat64Implementation()

    let k = Array<Int64>(m, item: 0)
    for (i in 0..k.size) {
        k[i] = i
    }
    let p = eye(m, m)
    for (i in 0..m - 1) {
        let j = i + rnd.nextInt64(m - i)
        (k[i], k[j]) = (k[j], k[i])
        nFloat64.dswap(m, p.data[i * p.stride..], 1, p.data[j * p.stride..], 1)
    }
    let kCopy = Array<Int64>(k.size, item: 0)
    copy(kCopy, k)

    let x = randomGeneral(m, n, ldx, rnd)

    var trans = NoTrans
    if (!fwd) {
        trans = Trans
    }
    let want = zeros(m, n, n)
    nFloat64.dgemm(trans, NoTrans, m, n, m, 1.0, p.data, p.stride, x.data, x.stride, 0.0, want.data, want.stride)

    lFloat64.dlapmr(fwd, m, n, x.data, x.stride, k)
    let got = x

    var kmod = false
    for (i in 0..k.size) {
        let ki = k[i]
        if (ki != kCopy[i]) {
            kmod = true
            break
        }
    }
    if (kmod) {
        throw Exception("${name}: unexpected modification of k")
    }

    if (!equalGeneral(got, want)) {
        throw Exception("${name}: unexpected result")
    }
    return true
}
