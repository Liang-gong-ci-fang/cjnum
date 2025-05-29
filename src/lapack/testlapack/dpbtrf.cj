package cjnum.lapack.testlapack

import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

// public func dpbtrfTest(lFloat64: LFloat64): Bool {
//     let rnd = Random(1)
//     for (n in [0, 1, 2, 3, 4, 5, 64, 65, 66, 91, 96, 97, 101, 128, 130]) {
//         for (kd in [(n + 1) / 4, (3 * n - 1) / 4, (5 * n + 1) / 4]) {
//             for (uplo in [Upper, Lower]) {
//                 for (ldab in [kd + 1, kd + 1 + 7]) {
//                     if (!dpbtrfTest(lFloat64, uplo, n, kd, ldab, rnd)) {
//                         return false
//                     }
//                 }
//             }
//         }
//     }
//     return true
// }

// func dpbtrfTest(lFloat64: LFloat64, uplo: Uplo, n: Int64, kd: Int64, ldab: Int64, rnd: Random): Bool {
//     const tol = 1e-12

//     let name = "uplo=${uplo.toString()},n=${n},kd=${kd},ldab=${ldab}"

//     let ab = randSymBand(uplo, n, kd, ldab, rnd)

//     let abFac = Array<Float64>(ab.size, item: 0.0)
//     copy(abFac, ab)
//     let ok = lFloat64.dpbtrf(uplo, n, kd, abFac, ldab)
//     if (!ok) {
//         throw Exception("${name}: bad test matrix, Dpbtrf failed")
//     }

//     dsbmm(uplo, n, kd, abFac, ldab)

//     let dist = distSymBand(uplo, n, kd, abFac, ldab, ab, ldab)
//     if (dist > tol * Float64(n)) {
//         throw Exception("${name}: unexpected result, diff=${dist}")
//     }
//     return true
// }

func dsbmm(uplo: Uplo, n: Int64, kd: Int64, ab: Array<Float64>, ldab: Int64) {
    let nFloat64 = nFloat64Implementation()
    if (uplo == Upper) {
        for (k in n - 1..=0 : -1) {
            let klen = min(kd, n - k - 1)
            if (klen > 0) {
                nFloat64.dsyr(Upper, klen, 1.0, ab[k * ldab + 1..], 1, ab[(k + 1) * ldab..], ldab - 1)
            }
            nFloat64.dscal(klen + 1, ab[k * ldab], ab[k * ldab..], 1)
        }
    } else if (uplo == Lower) {
        for (k in n - 1..=0 : -1) {
            let kc = max(0, kd - k)
            let klen = kd - kc
            ab[k * ldab + kd] = nFloat64.ddot(klen + 1, ab[k * ldab + kc..], 1, ab[k * ldab + kc..], 1)
            if (klen > 0) {
                nFloat64.dtrmv(Lower, NoTrans, NonUnitDiag, klen, ab[(k - klen) * ldab + kd..], ldab - 1,
                    ab[k * ldab + kc..], 1)
            }
        }
    }
}
