package cjnum.internal.asm.f32

import std.unittest.*
import std.unittest.testmacro.*

class l2NormUnitaryTestStruct {
    l2NormUnitaryTestStruct(
        let want: Float32,
        var x: Array<Float32>
    ) {}
}

@Test
public func testL2NormUnitary() {
    const tol: Float32 = 1e-7

    let src_gd: Float32 = 1.0
    let l2NormUnitaryTests = [
        l2NormUnitaryTestStruct(0.0, []),
        l2NormUnitaryTestStruct(2.0, [2.0]),
        l2NormUnitaryTestStruct(3.7416573867739413, [1.0, 2.0, 3.0]),
        l2NormUnitaryTestStruct(3.7416573867739413, [-1.0, -2.0, -3.0]),
        l2NormUnitaryTestStruct(nan, [nan]),
        l2NormUnitaryTestStruct(17.88854381999832, [8.0, -8.0, 8.0, -8.0, 8.0]),
        l2NormUnitaryTestStruct(2.23606797749979, [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0])
    ]
    for (j in 0..l2NormUnitaryTests.size) {
        let v = l2NormUnitaryTests[j]
        let g_ln = 4 + j % 2
        v.x = guardVector(v.x, src_gd, g_ln)
        let src = v.x[g_ln..v.x.size - g_ln]
        let ret = l2NormUnitary(src)
        @Expect(sameApprox(ret, v.want, tol))
        @Expect(isValidGuard(v.x, src_gd, g_ln))
    }
}

class l2NormIncTestStruct {
    l2NormIncTestStruct(
        let inc: Int64,
        let want: Float32,
        var x: Array<Float32>
    ) {}
}

@Test
public func testL2NormInc() {
    const tol: Float32 = 1e-7

    let src_gd: Float32 = 1.0
    let l2NormIncTests = [
        l2NormIncTestStruct(2, 0.0, []),
        l2NormIncTestStruct(3, 2.0, [2.0]),
        l2NormIncTestStruct(10, 3.7416573867739413, [1.0, 2.0, 3.0]),
        l2NormIncTestStruct(5, 3.7416573867739413, [-1.0, -2.0, -3.0]),
        l2NormIncTestStruct(3, nan, [nan]),
        l2NormIncTestStruct(15, 17.88854381999832, [8.0, -8.0, 8.0, -8.0, 8.0]),
        l2NormIncTestStruct(1, 2.23606797749979, [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0])
    ]
    for (j in 0..l2NormIncTests.size) {
        let v = l2NormIncTests[j]
        let g_ln = 4 + j % 2
        let ln = v.x.size
        v.x = guardIncVector(v.x, src_gd, v.inc, g_ln)
        let src = v.x[g_ln..v.x.size - g_ln]
        let ret = l2NormInc(src, ln, v.inc)
        @Expect(sameApprox(ret, v.want, tol))
        @Expect(checkValidIncGuard(v.x, src_gd, v.inc, g_ln))
    }
}

class l2DistanceUnitaryTestStruct {
    l2DistanceUnitaryTestStruct(
        let want: Float32,
        var x: Array<Float32>,
        var y: Array<Float32>
    ) {}
}

@Test
public func testL2DistanceUnitary() {
    const tol: Float32 = 1e-7

    let src_gd: Float32 = 1.0
    let l2DistanceUnitaryTests = [
        l2DistanceUnitaryTestStruct(0.0, [], []),
        l2DistanceUnitaryTestStruct(2.0, [3.0], [1.0]),
        l2DistanceUnitaryTestStruct(3.7416573867739413, [2.0, 4.0, 6.0], [1.0, 2.0, 3.0]),
        l2DistanceUnitaryTestStruct(3.7416573867739413, [1.0, 2.0, 3.0], [2.0, 4.0, 6.0]),
        l2DistanceUnitaryTestStruct(nan, [nan], [0.0]),
        l2DistanceUnitaryTestStruct(17.88854381999832, [9.0, -9.0, 9.0, -9.0, 9.0], [1.0, -1.0, 1.0, -1.0, 1.0]),
        l2DistanceUnitaryTestStruct(2.23606797749979, [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0],
            [0.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 2.0])
    ]
    for (j in 0..l2DistanceUnitaryTests.size) {
        let v = l2DistanceUnitaryTests[j]
        let g_ln = 4 + j % 2
        v.x = guardVector(v.x, src_gd, g_ln)
        v.y = guardVector(v.y, src_gd, g_ln)
        let srcX = v.x[g_ln..v.x.size - g_ln]
        let srcY = v.y[g_ln..v.y.size - g_ln]
        let ret = l2DistanceUnitary(srcX, srcY)
        @Expect(sameApprox(ret, v.want, tol))
        @Expect(isValidGuard(v.x, src_gd, g_ln))
    }
}
