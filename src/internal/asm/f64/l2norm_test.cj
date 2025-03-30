package cjnum.internal.asm.f64

import cjnum.floats.scalar.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

const tol = 1e-15

class l2NormUnitaryTestStruct {
    l2NormUnitaryTestStruct(
        let want: Float64,
        let x: Array<Float64>
    ) {}
}

@Test
func TestL2NormUnitary(){  
	var src_gd: Float64 = 1.0
    let l2NormUnitaryTests=[
        l2NormUnitaryTestStruct(0.0, Array<Float64>([])),
        l2NormUnitaryTestStruct(2.0, Array<Float64>([2.0])),
        l2NormUnitaryTestStruct(3.7416573867739413, Array<Float64>([1.0, 2.0, 3.0])),
        l2NormUnitaryTestStruct(3.7416573867739413, Array<Float64>([-1.0, -2.0, -3.0])),
        l2NormUnitaryTestStruct(nan, Array<Float64>([nan])),
        l2NormUnitaryTestStruct(nan, Array<Float64>([1.0, inf, 3.0, naNWith(25), 5.0])),
        l2NormUnitaryTestStruct(17.88854381999832, Array<Float64>([8.0, -8.0, 8.0, -8.0, 8.0])),
        l2NormUnitaryTestStruct(2.23606797749979, Array<Float64>([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0]))
    ]
    for (j in 0..l2NormUnitaryTests.size){
        let v= l2NormUnitaryTests[j]
        let g_ln: Int64 = 4 + j % 2
        let guardedX = guardVector(v.x, src_gd, g_ln)
        let src = guardedX[g_ln..guardedX.size - g_ln]
        let ret = L2NormUnitary(src)
        @Expect(sameApprox(ret, v.want, tol))
        @Expect(isValidGuard(guardedX, src_gd, g_ln))
    }
}

class l2NormIncTestStruct {
    l2NormIncTestStruct(
        let inc: Int64,
        let want: Float64,
        let x: Array<Float64>
    ) {}
}

@Test
public func TestL2NormInc() {   
    var src_gd: Float64 = 1.0
    let l2NormIncTests = [
        l2NormIncTestStruct(2, 0.0, Array<Float64>([])),
        l2NormIncTestStruct(3, 2.0, Array<Float64>([2.0])),
        l2NormIncTestStruct(10, 3.7416573867739413, Array<Float64>([1.0, 2.0, 3.0])),
        l2NormIncTestStruct(5, 3.7416573867739413, Array<Float64>([-1.0, -2.0, -3.0])),
        l2NormIncTestStruct(3, nan, Array<Float64>([nan])),
        l2NormIncTestStruct(15, 17.88854381999832, Array<Float64>([8.0, -8.0, 8.0, -8.0, 8.0])),
        l2NormIncTestStruct(1, 2.23606797749979, Array<Float64>([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0]))
    ]

    for (j in 0..l2NormIncTests.size) {
        let v = l2NormIncTests[j]
        let g_ln: Int64 = 4 + j % 2
        let ln = v.x.size
        let guardedX = guardIncVector(v.x, src_gd, v.inc, g_ln)
        let src = guardedX[g_ln..guardedX.size - g_ln]
        let ret = L2NormInc(src, ln, v.inc)
        @Expect(sameApprox(ret, v.want, tol))
        @Expect(checkValidIncGuard(guardedX, src_gd, v.inc, g_ln))
    }
}

class l2DistanceUnitaryTestStruct {
    l2DistanceUnitaryTestStruct(
        let want: Float64,
        let x: Array<Float64>,
        let y: Array<Float64>
    ) {}
}

@Test
func TestL2DistanceUnitary(){
	var src_gd: Float64 = 1.0
    let l2DistanceUnitaryTests=[
        l2DistanceUnitaryTestStruct(0.0, Array<Float64>([]),Array<Float64>([])),
        l2DistanceUnitaryTestStruct(2.0, Array<Float64>([3.0]),Array<Float64>([1.0])),
        l2DistanceUnitaryTestStruct(3.7416573867739413, Array<Float64>([2.0, 4.0, 6.0]),Array<Float64>([1.0, 2.0, 3.0])),
        l2DistanceUnitaryTestStruct(3.7416573867739413, Array<Float64>([1.0, 2.0, 3.0]),Array<Float64>([2.0, 4.0, 6.0])),
        l2DistanceUnitaryTestStruct(nan, Array<Float64>([nan]),Array<Float64>([0.0])),
        l2DistanceUnitaryTestStruct(17.88854381999832, Array<Float64>([9.0, -9.0, 9.0, -9.0, 9.0]),Array<Float64>([1.0, -1.0, 1.0, -1.0, 1.0])),
        l2DistanceUnitaryTestStruct(2.23606797749979, Array<Float64>([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0]), Array<Float64>([0.0, 2.0, 0.0, -2.0, 0.0, 2.0, 0.0, -2.0, 0.0, 2.0]))
    ]
    for (j in 0..l2DistanceUnitaryTests.size){
        let v= l2DistanceUnitaryTests[j]
        let g_ln: Int64 = 4 + j % 2
        let guardedX = guardVector(v.x, src_gd, g_ln)      
        let guardedY = guardVector(v.y, src_gd, g_ln)
        let srcX = guardedX[g_ln..guardedX.size - g_ln]
        let srcY = guardedY[g_ln..guardedY.size - g_ln]
        let ret = L2DistanceUnitary(srcX,srcY)
        @Expect(sameApprox(ret, v.want, tol))
        @Expect(isValidGuard(guardedX, src_gd, g_ln))
    }
}

func netlib(x: Array<Float64>) :Float64 {
        var scale: Float64 = 0.0
        var sumSquares: Float64 = 1.0
        for (v in x) {
            if (v == 0.0) {
                continue
            }
            var absxi = abs(v)
            if (absxi.isNaN()) {
                return Float64.NaN
            }
            if (scale < absxi) {
                var s = scale / absxi
                sumSquares = 1.0 + sumSquares * s * s
                scale = absxi
            } else {
                var s = scale / absxi
                sumSquares += s * s
            }
        }
        if (scale.isInf()) {
            return Float64.Inf
        }
        return scale * sqrt(sumSquares)
    }
   
@Test 
@Configure[baseline: "testNetlib"]
class L2NormNetlibBenchmarks{
    let num = randomSlice(1, 1)[0]
    let lengths: Array<Int64> = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000] 
     
    @Bench[ln in lengths]
    func testNetlib(ln:Int64): Unit 
    {
        var x = Array(ln) { i => Float64(i) * 1.0 }
        x[0] = num
        netlib(x)
    }

    @Bench[ln in lengths]
    func testL2NormUnitary(ln:Int64): Unit 
    {
        var x = Array(ln) { i => Float64(i) * 1.0 }
        x[0] = num
        L2NormUnitary(x)
    }
}