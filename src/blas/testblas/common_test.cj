package cjnum.blas.testblas

import std.math.*
import cjnum.blas.*
import cjnum.floats.*
import std.unittest.*
import std.unittest.testmacro.*

class flattenBandedStruct {
    flattenBandedStruct(
        let dense!: Array<Array<Float64>>,
        let ku!: Int64,
        let kl!: Int64,
        let condensed!: Array<Array<Float64>>
    ) {}
}

@Test
public func testFlattenBanded() {
    let flattenBandedTests = [
        flattenBandedStruct(
            dense: [[3.0]],
            ku: 0,
            kl: 0,
            condensed: [[3.0]]
        ),
        flattenBandedStruct(
            dense: [[3.0, 4.0, 0.0]],
            ku: 1,
            kl: 0,
            condensed: [[3.0, 4.0]]
        ),
        flattenBandedStruct(
            dense: [[3.0, 4.0, 0.0, 0.0, 0.0, 0.0]],
            ku: 1,
            kl: 0,
            condensed: [[3.0, 4.0]]
        ),
        flattenBandedStruct(
            dense: [
                [3.0, 4.0, 0.0],
                [0.0, 5.0, 8.0],
                [0.0, 0.0, 2.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]
            ],
            ku: 1,
            kl: 0,
            condensed: [
                [3.0, 4.0],
                [5.0, 8.0],
                [2.0, Float64.NaN],
                [Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN]
            ]
        ),
        flattenBandedStruct(
            dense: [
                [3.0, 4.0, 6.0],
                [0.0, 5.0, 8.0],
                [0.0, 0.0, 2.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]
            ],
            ku: 2,
            kl: 0,
            condensed: [
                [3.0, 4.0, 6.0],
                [5.0, 8.0, Float64.NaN],
                [2.0, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN]
            ]
        ),
        flattenBandedStruct(
            dense: [
                [3.0, 4.0, 6.0],
                [1.0, 5.0, 8.0],
                [0.0, 6.0, 2.0],
                [0.0, 0.0, 7.0],
                [0.0, 0.0, 0.0]
            ],
            ku: 2,
            kl: 1,
            condensed: [
                [Float64.NaN, 3.0, 4.0, 6.0],
                [1.0, 5.0, 8.0, Float64.NaN],
                [6.0, 2.0, Float64.NaN, Float64.NaN],
                [7.0, Float64.NaN, Float64.NaN, Float64.NaN],
                [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]
            ]
        ),
        flattenBandedStruct(
            dense: [
                [1.0, 2.0, 0.0],
                [3.0, 4.0, 5.0],
                [6.0, 7.0, 8.0],
                [0.0, 9.0, 10.0],
                [0.0, 0.0, 11.0]
            ],
            ku: 1,
            kl: 2,
            condensed: [
                [Float64.NaN, Float64.NaN, 1.0, 2.0],
                [Float64.NaN, 3.0, 4.0, 5.0],
                [6.0, 7.0, 8.0, Float64.NaN],
                [9.0, 10.0, Float64.NaN, Float64.NaN],
                [11.0, Float64.NaN, Float64.NaN, Float64.NaN]
            ]
        ),
        flattenBandedStruct(
            dense: [
                [1.0, 0.0, 0.0],
                [3.0, 4.0, 0.0],
                [6.0, 7.0, 8.0],
                [0.0, 9.0, 10.0],
                [0.0, 0.0, 11.0]
            ],
            ku: 0,
            kl: 2,
            condensed: [
                [Float64.NaN, Float64.NaN, 1.0],
                [Float64.NaN, 3.0, 4.0],
                [6.0, 7.0, 8.0],
                [9.0, 10.0, Float64.NaN],
                [11.0, Float64.NaN, Float64.NaN]
            ]
        ),
        flattenBandedStruct(
            dense: [
                [1.0, 0.0, 0.0, 0.0, 0.0],
                [3.0, 4.0, 0.0, 0.0, 0.0],
                [1.0, 3.0, 5.0, 0.0, 0.0]
            ],
            ku: 0,
            kl: 2,
            condensed: [
                [Float64.NaN, Float64.NaN, 1.0],
                [Float64.NaN, 3.0, 4.0],
                [1.0, 3.0, 5.0]
            ]
        )
    ]
    var count = 0
    for (test in flattenBandedTests) {
        let condensed = flattenBanded(test.dense, test.ku, test.kl)
        let correct = flatten(test.condensed)
        @Expect(fSame(condensed, correct))
        count++
    }
}

class flattenTriangularStruct {
    flattenTriangularStruct(
        let a!: Array<Array<Float64>>,
        let ans!: Array<Float64>,
        let ul!: Uplo
    ) {}
}

@Test
public func TestFlattenTriangular() {
    let flattenTriangularTests = [
        flattenTriangularStruct(
            a: [
                [1.0, 2.0, 3.0],
                [0.0, 4.0, 5.0],
                [0.0, 0.0, 6.0]
            ],
            ans: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            ul: Upper),
        flattenTriangularStruct(
            a: [
                [1.0, 0.0, 0.0],
                [2.0, 3.0, 0.0],
                [4.0, 5.0, 6.0]
            ],
            ans: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            ul: Lower)
    ]
    for (test in flattenTriangularTests) {
        let a = flattenTriangular(test.a, test.ul)
        @Expect(fEqual(a, test.ans))
    }
}

// TODO complext test