package cjnum.mat

import cjnum.blas.blas64.*
import std.unittest.*
import std.unittest.testmacro.*

private class NewDiagDenseTest {
    let data: Array<Float64>
    let n: Int64
    let mat: DiagDense
    let dense: Dense
    public init(data: Array<Float64>, n: Int64, mat: DiagDense, dense: Dense) {
        this.data = data
        this.n = n
        this.mat = mat
        this.dense = dense
    }
}

private let newDiagDenseTests: Array<NewDiagDenseTest> = [
    NewDiagDenseTest(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        6,
        DiagDense(Vector(n: 6, inc: 1, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    )
]

@Test
public func testNewDiagDense(): Unit {
    for (i in 0..newDiagDenseTests.size) {
        let test = newDiagDenseTests[i]
        let band = newDiagDense(test.n, test.data)
        let (rows, cols) = band.dims()
        @Expect(rows,test.n)
        @Expect(cols,test.n)
        @Expect(band==test.mat)
        @Expect(mEqual(band,test.mat))
        @Expect(mEqual(band,test.dense))
    }
}

private let diagDenseZeroTests: Array<DiagDense> = [
    DiagDense(
        Vector(
            n: 5,
            inc: 2,
            data: [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
        )
    )
]

@Test
public func testDiagDenseZero(): Unit {
    // Elements that equal 1 should be set to zero, elements that equal -1
    // should remain unchanged.
    for (test in diagDenseZeroTests) {
        var dataCopy = Array<Float64>(test.mat.data.size, repeat: 0.0)
        dataCopy = test.mat.data[0..]
        test.zero()
        for (i in 0..test.mat.data.size) {
            let v = test.mat.data[i]
            @Expect(!(dataCopy[i] != -1.0 && v != 0.0))
            @Expect(!(dataCopy[i] == -1.0 && v != -1.0))
        }
    }
}

private class DiagonalStrideTest {
    let diag: DiagDense
    let dense: Dense
    public init(diag: DiagDense, dense: Dense) {
        this.diag = diag
        this.dense = dense
    }
}

private let diagonalStrideTests: Array<DiagonalStrideTest> = [
    DiagonalStrideTest(
        DiagDense(Vector(n: 6, inc: 1, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagonalStrideTest(
        DiagDense(Vector(n: 6, inc: 2, data: [1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0, 0.0, 6.0]),),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagonalStrideTest(
        DiagDense(
            Vector(n: 6, inc: 5,
                data: [1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0,
                    0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 6.0]),),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    )
]

@Test
public func TestDiagonalStride(): Unit {
    for (test in diagonalStrideTests) {
        @Expect(mEqual(test.diag, test.dense))
    }
}

private class DiagFromTest {
    let mat: Matrix
    let want: Dense
    public init(mat: Matrix, want: Dense) {
        this.mat = mat
        this.want = want
    }
}

private let diagFromTests: Array<DiagFromTest> = [
    DiagFromTest(
        newDiagDense(6, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newBandDense(6, 6, 1, 1,
            [Float64.NaN, 1.0, Float64.NaN, Float64.NaN, 2.0, Float64.NaN, Float64.NaN, 3.0, Float64.NaN, Float64.NaN,
                4.0, Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, 6.0, Float64.NaN]),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newDense(6, 6,
            [1.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 2.0, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 3.0, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 4.0, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, 6.0]),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newDense(6, 4,
            [1.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 2.0, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, 3.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 4.0, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]),
        newDense(4, 4, [1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 4.0]),
    ),
    DiagFromTest(
        newDense(4, 6,
            [1.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 2.0, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 3.0, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 4.0, Float64.NaN, Float64.NaN]),
        newDense(4, 4, [1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 4.0]),
    ),
    DiagFromTest(
        newSymBandDense(6, 1,
            [1.0, Float64.NaN, 2.0, Float64.NaN, 3.0, Float64.NaN, 4.0, Float64.NaN, 5.0, Float64.NaN, 6.0, Float64.NaN]
        ),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newSymDense(
            6,
            [1.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 2.0, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 3.0, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 4.0, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, 6.0]
        ),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newTriBandDense(6, 2, mUpper,
            [1.0, Float64.NaN, Float64.NaN, 2.0, Float64.NaN, Float64.NaN, 3.0, Float64.NaN, Float64.NaN, 4.0,
                Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, 6.0, Float64.NaN, Float64.NaN]),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newTriBandDense(6, 2, mLower,
            [Float64.NaN, Float64.NaN, 1.0, Float64.NaN, Float64.NaN, 2.0, Float64.NaN, Float64.NaN, 3.0, Float64.NaN,
                Float64.NaN, 4.0, Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, 6.0]),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newTriDense(6, mUpper,
            [1.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 2.0, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 3.0, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 4.0, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, 6.0]),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newTriDense(6, mLower,
            [1.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 2.0, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 3.0, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 4.0, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN,
                Float64.NaN, Float64.NaN, 6.0]),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    DiagFromTest(
        newVecDense(6, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        newDense(1, 1, [1.0]),
    ),
    DiagFromTest(
        BasicMatrix(
            General(
                rows: 6,
                cols: 6,
                stride: 6,
                data: [1.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 2.0,
                    Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 3.0, Float64.NaN,
                    Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 4.0, Float64.NaN, Float64.NaN,
                    Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 5.0, Float64.NaN, Float64.NaN, Float64.NaN,
                    Float64.NaN, Float64.NaN, Float64.NaN, 6.0],
            ),
            6,
            6,
        ),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    )
]

@Test
public func testDiagFrom(): Unit {
    for (i in 0..diagFromTests.size) {
        let test = diagFromTests[i]
        var got = DiagDense(Vector())
        got.diagFrom(test.mat)
        @Expect(mEqual(got, test.want))
    }
}

public func testDiagView(_: Int64, test: DiagDenseViewer): Unit {
    // Check the DiagView matches the Diagonal.
    let (r, c) = test.dims()
    let diagView = test.diagView()
    for (i in 0..min(r, c)) {
        @Expect(diagView.at(i, i),test.at(i, i))
    }

    // Check that changes to the diagonal are reflected.
    let offset = 10.0
    let diag = diagView as DiagDense
    for (i in 0..min(r, c)) {
        let v = test.at(i, i)
        if (!diag.isNone()) {
            diag.getOrThrow().setDiag(i, v + offset)
            @Expect(test.at(i, i),v+offset)
        }
    }

    // Check that DiagView and DiagFrom match.
    let diag2: DiagDense = DiagDense(Vector())
    diag2.diagFrom(test)
    @Expect(mEqual(diag.getOrThrow(), diag2))
}

@Test
public func testDiagonalAtSet(): Unit {
    for (n in [1, 3, 8]) {
        for (nilstart in [true, false]) {
            var diag = DiagDense(Vector())
            if (nilstart) {
                diag = newDiagDense(n, Array<Float64>(0, repeat: 0.0))
            } else {
                var data = Array<Float64>(n, repeat: 0.0)
                diag = newDiagDense(n, data)
                for (i in 0..data.size) {
                    data[i] = -Float64(i) - 1.0
                    let v = diag.at(i, i)
                    @Expect(v == data[i])
                }
            }
            for (i in 0..n) {
                for (j in 0..n) {
                    if (i != j) {
                        @Expect(diag.at(i, j) == 0.0)
                    }
                    let v = Float64(i) + 1.0
                    diag.setDiag(i, v)
                    let v2 = diag.at(i, i)
                    @Expect(v2 == v)
                }
            }
        }
    }
}

/*
因为缺乏对randDiagDense方法的测试，所以这里没有实现该方法
 */
// TODO func randDiagDense(size int, rnd *rand.Rand) *DiagDense

/*
因为缺乏类对象的比较函数，所以在这里让DiagDense类实现了Equatable<DiagDense>接口，
完成比较运算符的实现
 */
extend DiagDense {
    public operator func ==(rhs: DiagDense): Bool {
        var judge = true
        for (i in 0..max(this.mat.data.size, rhs.mat.data.size)) {
            if (this.mat.data[i] != rhs.mat.data[i]) {
                judge = false
            }
        }
        return (this.mat.n == rhs.mat.n && this.mat.inc == rhs.mat.inc && judge)
    }

    public operator func !=(rhs: DiagDense): Bool {
        var judge = false
        for (i in 0..max(this.mat.data.size, rhs.mat.data.size)) {
            if (this.mat.data[i] != rhs.mat.data[i]) {
                judge = true
            }
        }
        return (this.mat.n == rhs.mat.n || this.mat.inc == rhs.mat.inc || judge)
    }
}
