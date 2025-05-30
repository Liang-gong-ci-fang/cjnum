package cjnum.lapack.testlapack

import cjnum.blas.*
import cjnum.lapack.*
import std.random.*

public func dpbtf2Test(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 20]) {
        for (kd in [0, (n + 1) / 4, (3 * n - 1) / 4, (5 * n + 1) / 4]) {
            for (uplo in [Upper, Lower]) {
                for (ldab in [kd + 1, kd + 1 + 7]) {
                    if (!dpbtf2Test(lFloat64, rnd, uplo, n, kd, ldab)) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func dpbtf2Test(lFloat64: LFloat64, rnd: Random, uplo: Uplo, n: Int64, kd: Int64, ldab: Int64): Bool {
    const tol = 1e-12

    let name = "uplo=${uplo},n=${n},kd=${kd},ldab=${ldab}"

    let ab = randSymBand(uplo, n, kd, ldab, rnd)

    let abFac = Array<Float64>(ab.size, item: 0.0)
    copy(abFac, ab)
    let ok = lFloat64.dpbtf2(uplo, n, kd, abFac, ldab)
    if (!ok) {
        throw Exception("${name}: bad test matrix, Dpbtf2 failed")
    }

    dsbmm(uplo, n, kd, abFac, ldab)

    let dist = distSymBand(uplo, n, kd, abFac, ldab, ab, ldab)
    if (dist > tol) {
        throw Exception("${name}: unexpected result, diff=${dist}")
    }
    return true
}
