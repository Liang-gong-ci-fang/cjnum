package cjnum.blas.internal.f64

import std.unittest.*
import std.unittest.testmacro.*

internal class axpyTestStruct {
    internal axpyTestStruct(
        alpha: Float64,
        x: Array<Float64>,
        y: Array<Float64>,
        want: Array<Float64>,
        wantRev: Array<Float64>
    ) {}
}

let axpyTests = [
    axpyTestStruct(0.0, Array<Float64>(), Array<Float64>(), Array<Float64>(), Array<Float64>()),
    axpyTestStruct(0.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([-3.0]), Array<Float64>([-3.0])),
    axpyTestStruct(1.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([-1.0]), Array<Float64>([-1.0])),
    axpyTestStruct(3.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([3.0]), Array<Float64>([3.0])),
    axpyTestStruct(-3.0, Array<Float64>([2.0]), Array<Float64>([-3.0]), Array<Float64>([-9.0]), Array<Float64>([-9.0])),
    axpyTestStruct(1.0, Array<Float64>([1.0, 5.0]), Array<Float64>([2.0, -3.0]), Array<Float64>([3.0, 2.0]), Array<Float64>([7.0, -2.0])),
    axpyTestStruct(1.0, Array<Float64>([2.0, 3.0, 4.0]), Array<Float64>([-3.0, -2.0, -1.0]), Array<Float64>([-1.0, 1.0, 3.0]), 
        Array<Float64>([1.0, 1.0, 1.0])),
    axpyTestStruct(0.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]), Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), 
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0])),
    axpyTestStruct(1.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]), Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), 
        Array<Float64>([0.0, 1.0, 1.0, 4.0, -2.0, 2.0, -10.0]), Array<Float64>([-4.0, -2.0, 2.0, 4.0, -3.0, 5.0, -6.0])),
    axpyTestStruct(3.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]), Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), 
        Array<Float64>([0.0, 1.0, 3.0, 6.0, 2.0, -4.0, -18.0]), Array<Float64>([-12.0, -8.0, 6.0, 6.0, -1.0, 5.0, -6.0])),
    axpyTestStruct(-3.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0, 0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0]), 
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0, 0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0]), 
        Array<Float64>([0.0, 1.0, -3.0, 0.0, -10.0, 14.0, 6.0, 0.0, 1.0, -3.0, 0.0, -10.0, 14.0, 6.0]), 
        Array<Float64>([12.0, 10.0, -6.0, 0.0, -7.0, 5.0, -6.0, 12.0, 10.0, -6.0, 0.0, -7.0, 5.0, -6.0])),
    axpyTestStruct(-5.0, Array<Float64>([0.0, 0.0, 1.0, 1.0, 2.0, -3.0, -4.0, 5.0, 1.0, 2.0, -3.0, -4.0, 5.0]), 
        Array<Float64>([0.0, 1.0, 0.0, 3.0, -4.0, 5.0, -6.0, 7.0, 3.0, -4.0, 5.0, -6.0, 7.0]), 
        Array<Float64>([0.0, 1.0, -5.0, -2.0, -14.0, 20.0, 14.0, -18.0, -2.0, -14.0, 20.0, 14.0, -18.0]), 
        Array<Float64>([-25.0, 21.0, 15.0, -7.0, -9.0, -20.0, 14.0, 22.0, -7.0, -9.0, 0.0, -6.0, 7.0]))
]

@Test
public func testAxpyUnitary() {
    // TODO 
}