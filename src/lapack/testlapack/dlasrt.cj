package cjnum.lapack.testlapack

import cjnum.floats.*
import cjnum.lapack.*

class dlasrtTestStruct {
    dlasrtTestStruct(
        let data!: Array<Float64>,
        let wantInc!: Array<Float64>,
        let wantDec!: Array<Float64>
    ) {}
}

public func dlasrtTest(lFloat64: LFloat64): Bool {
    var ti = 0
    for (test in [
            dlasrtTestStruct(
                data: [],
                wantInc: [],
                wantDec: []
            ),
            dlasrtTestStruct(
                data: Array<Float64>([]),
                wantInc: Array<Float64>([]),
                wantDec: Array<Float64>([])
            ),
            dlasrtTestStruct(
                data: [1.0],
                wantInc: [1.0],
                wantDec: [1.0]
            ),
            dlasrtTestStruct(
                data: [1.0, 2.0],
                wantInc: [1.0, 2.0],
                wantDec: [2.0, 1.0]
            ),
            dlasrtTestStruct(
                data: [1.0, 2.0, -3.0],
                wantInc: [-3.0, 1.0, 2.0],
                wantDec: [2.0, 1.0, -3.0]
            ),
            dlasrtTestStruct(
                data: [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                wantInc: [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                wantDec: [5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0, -5.0]
            ),
            dlasrtTestStruct(
                data: [5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0, -5.0],
                wantInc: [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                wantDec: [5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0, -5.0]
            ),
            dlasrtTestStruct(
                data: [-2.0, 4.0, -1.0, 2.0, -4.0, 0.0, 3.0, 5.0, -5.0, 1.0, -3.0],
                wantInc: [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                wantDec: [5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0, -5.0]
            )
        ]) {
        let n = test.data.size
        let ds = Array<Float64>(n, item: 0.0)

        copy(ds, test.data)
        lFloat64.dlasrt(SortIncreasing, n, ds)
        if (!fEqual(ds, test.wantInc)) {
            throw Exception("Case #${ti}: unexpected result of SortIncreasing")
        }

        copy(ds, test.data)
        lFloat64.dlasrt(SortDecreasing, n, ds)
        if (!fEqual(ds, test.wantDec)) {
            throw Exception("Case #${ti}: unexpected result of SortIncreasing")
        }
        ti++
    }
    return true
}
