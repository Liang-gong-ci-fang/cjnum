package cjnum.lapack.testlapack

import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*

class dpotf2TestStruct {
    dpotf2TestStruct(
        let a!: Array<Array<Float64>>,
        let pos!: Bool,
        let U!: Array<Array<Float64>>
    ) {}
}

public func dpotf2Test(lFloat64: LFloat64): Bool {
    for (test in [
            dpotf2TestStruct(
                a: [
                    [23.0, 37.0, 34.0, 32.0],
                    [108.0, 71.0, 48.0, 48.0],
                    [109.0, 109.0, 67.0, 58.0],
                    [106.0, 107.0, 106.0, 63.0]
                ],
                pos: true,
                U: [
                    [4.795831523312719, 7.715033320111766, 7.089490077940543, 6.672461249826393],
                    [0.0, 3.387958215439679, -1.976308959006481, -1.026654004678691],
                    [0.0, 0.0, 3.582364210034111, 2.419258947036024],
                    [0.0, 0.0, 0.0, 3.401680257083044]
                ],
            ),
            dpotf2TestStruct(
                a: [[8.0, 2.0], [2.0, 4.0]],
                pos: true,
                U: [[2.82842712474619, 0.707106781186547], [0.0, 1.870828693386971]],
            )
        ]) {
        if (!testDpotf2(lFloat64, test.pos, test.a, test.U, test.a[0].size, Upper)) {
            return false
        }
        if (!testDpotf2(lFloat64, test.pos, test.a, test.U, test.a[0].size + 5, Upper)) {
            return false
        }
        let aT = transpose(test.a)
        let L = transpose(test.U)
        if (!testDpotf2(lFloat64, test.pos, aT, L, test.a[0].size, Lower)) {
            return false
        }
        if (!testDpotf2(lFloat64, test.pos, aT, L, test.a[0].size + 5, Lower)) {
            return false
        }
    }
    return true
}

func testDpotf2(lFloat64: LFloat64, testPos: Bool, a: Array<Array<Float64>>, ans: Array<Array<Float64>>, stride: Int64,
    ul: Uplo): Bool {
    let aFlat = flattenTri(a, stride, ul)
    let ansFlat = flattenTri(ans, stride, ul)
    let pos = lFloat64.dpotf2(ul, a[0].size, aFlat, stride)
    if (pos != testPos) {
        throw Exception("Positive definite mismatch: Want ${testPos}, Got ${pos}")
        // return
    }
    if (testPos && !fEqualApprox(ansFlat, aFlat, 1e-14)) {
        throw Exception("Result mismatch: Want ${ansFlat}, Got  ${aFlat}")
    }
    return true
}

func flattenTri(a: Array<Array<Float64>>, stride: Int64, ul: Uplo): Array<Float64> {
    let m = a.size
    let n = a[0].size
    if (stride < n) {
        throw Exception("bad stride")
    }
    let upper = ul == Upper
    let v = Array<Float64>(m * stride, repeat: 0.0)
    var count = 1000.0
    for (i in 0..m) {
        for (j in 0..stride) {
            if (j >= n || (upper && j < i) || (!upper && j > i)) {
                v[i * stride + j] = count
                count += 1.0
                continue
            }
            v[i * stride + j] = a[i][j]
        }
    }
    return v
}

func transpose(a: Array<Array<Float64>>): Array<Array<Float64>> {
    let m = a.size
    let n = a[0].size
    if (m != n) {
        throw Exception("not square")
    }
    let aNew = Array<Array<Float64>>(m, repeat: Array<Float64>(n, repeat: 0.0))
    for (i in 0..m) {
        aNew[i] = Array<Float64>(n, repeat: 0.0)
    }
    for (i in 0..m) {
        if (a[i].size != n) {
            throw Exception("bad n size")
        }
        for (j in 0..n) {
            aNew[j][i] = a[i][j]
        }
    }
    return aNew
}
