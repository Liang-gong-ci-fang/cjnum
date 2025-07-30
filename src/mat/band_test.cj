package cjnum.mat

import cjnum.blas.blas64.*
import cjnum.floats.scalar.*
import std.reflect.*
import std.math.*
import std.core.*
import std.unittest.*
import std.unittest.testmacro.*

private class NewBandTest {
    let data: Array<Float64>
    let r: Int64
    let c: Int64
    let kl: Int64
    let ku: Int64
    let mat: BandDense
    let dense: Dense
    public init(data: Array<Float64>, r: Int64, c: Int64, kl: Int64, ku: Int64, mat: BandDense, dense: Dense) {
        this.data = data
        this.r = r
        this.c = c
        this.kl = kl
        this.ku = ku
        this.mat = mat
        this.dense = dense
    }
}

private let newBandTests = [
    NewBandTest(
        [-1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, -1.0,
            19.0, 20.0, -1.0, -1.0],
        6,
        6,
        1,
        2,
        BandDense(
            Band(
                rows: 6,
                cols: 6,
                kl: 1,
                ku: 2,
                stride: 4,
                data: [-1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
                    18.0, -1.0, 19.0, 20.0, -1.0, -1.0]
            )
        ),
        newDense(6, 6,
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 0.0, 0.0, 0.0, 8.0, 9.0, 10.0, 11.0, 0.0, 0.0, 0.0, 12.0,
                13.0, 14.0, 15.0, 0.0, 0.0, 0.0, 16.0, 17.0, 18.0, 0.0, 0.0, 0.0, 0.0, 19.0, 20.0]),
    ),
    NewBandTest(
        [-1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, -1.0,
            19.0, 20.0, -1.0, -1.0, 21.0, -1.0, -1.0, -1.0],
        10,
        6,
        1,
        2,
        BandDense(
            Band(
                rows: 10,
                cols: 6,
                kl: 1,
                ku: 2,
                stride: 4,
                data: [-1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
                    18.0, -1.0, 19.0, 20.0, -1.0, -1.0, 21.0, -1.0, -1.0, -1.0]
            )
        ),
        newDense(10, 6,
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 0.0, 0.0, 0.0, 8.0, 9.0, 10.0, 11.0, 0.0, 0.0, 0.0, 12.0,
                13.0, 14.0, 15.0, 0.0, 0.0, 0.0, 16.0, 17.0, 18.0, 0.0, 0.0, 0.0, 0.0, 19.0, 20.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    ),
    NewBandTest(
        [-1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
            20.0, 21.0, 22.0, 23.0],
        6,
        10,
        1,
        2,
        BandDense(
            Band(
                rows: 6,
                cols: 10,
                kl: 1,
                ku: 2,
                stride: 4,
                data: [-1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
                    18.0, 19.0, 20.0, 21.0, 22.0, 23.0]
            )
        ),
        newDense(6, 10,
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                8.0, 9.0, 10.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 13.0, 14.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 16.0, 17.0, 18.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 21.0, 22.0, 23.0, 0.0, 0.0]),
    )
]

@Test
public func TestNewBand(): Unit {
    for (test in newBandTests) {
        let band = newBandDense(test.r, test.c, test.kl, test.ku, test.data)
        let (rows, cols) = band.dims()
        @Expect(rows ,test.r)
        @Expect(cols ,test.c)
        @Expect(band==test.mat)
        @Expect(mEqual(band, test.mat))
        @Expect(mEqual(band, test.dense))
    }
}

private class NewDiagonalRectTest {
    let data: Array<Float64>
    let r: Int64
    let c: Int64
    let mat: BandDense
    let dense: Dense
    public init(data: Array<Float64>, r: Int64, c: Int64, mat: BandDense, dense: Dense) {
        this.data = data
        this.r = r
        this.c = c
        this.mat = mat
        this.dense = dense
    }
}

private let newDiagonalRectTests: Array<NewDiagonalRectTest> = [
    NewDiagonalRectTest(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        6,
        6,
        BandDense(
            Band(
                rows: 6,
                cols: 6,
                stride: 1,
                data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            ),
        ),
        newDense(6, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]),
    ),
    NewDiagonalRectTest(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        7,
        6,
        BandDense(
            Band(
                rows: 7,
                cols: 6,
                stride: 1,
                data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            ),
        ),
        newDense(7, 6,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    ),
    NewDiagonalRectTest(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        6,
        7,
        BandDense(
            Band(
                rows: 6,
                cols: 7,
                stride: 1,
                data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            ),
        ),
        newDense(6, 7,
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0]
        ),
    )
]

@Test
public func TestNewDiagonalRect(): Unit {
    for (test in newDiagonalRectTests) {
        let band = newDiagonalRect(test.r, test.c, test.data)
        let (rows, cols) = band.dims()

        @Expect(rows,test.r)
        @Expect(cols ,test.c)
        @Expect(band==test.mat)
        @Expect(mEqual(band,test.mat))
        @Expect(mEqual(band,test.dense))
    }
}

private let bandDenseZeroTests: Array<BandDense> = [
    BandDense(
        Band(
            rows: 6,
            cols: 7,
            stride: 8,
            kl: 1,
            ku: 2,
            data: [-1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0,
                1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0,
                -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        )
    ),
    BandDense(
        Band(
            rows: 6,
            cols: 7,
            stride: 8,
            kl: 2,
            ku: 1,
            data: [-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0,
                1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0,
                -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
        )
    )
]

public func testBandDenseZero(): Unit {
    for (test in bandDenseZeroTests) {
        let dataCopy = Array<Float64>(test.mat.data.size, repeat: 0.0)
        test.mat.data.copyTo(dataCopy, 0, 0, test.mat.data.size)
        test.zero()
        for (i in 0..test.mat.data.size) {
            let v = test.mat.data[i]
            @Expect(!(dataCopy[i]!=-1.0&&v!=0.0))
            @Expect(!(dataCopy[i] == -1.0 && v != -1.0))
        }
    }
}

private let bandDiagViewTests: Array<BandDense> = [
    newBandDense(1, 1, 0, 0, [1.0]),
    newBandDense(6, 6, 1, 2,
        [-1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, -1.0,
            21.0, 22.0, -1.0, -1.0]),
    newBandDense(6, 6, 2, 1,
        [-1.0, -1.0, 1.0, 2.0, -1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
            18.0, 19.0, 20.0, -1.0])
]

@Test
public func testBandDiagView(): Unit {
    for (cas in 0..bandDiagViewTests.size) {
        let test = bandDiagViewTests[cas]
        testDiagView(cas, test)
    }
}

// bandImplicit is an implicit band matrix returning val(i, j)
// for the value at (i, j).
private class BandImplicit <: Matrix {
    let r: Int64
    let c: Int64
    let kl: Int64
    let ku: Int64
    let val: (i: Int64, j: Int64) -> Float64
    public init(r: Int64, c: Int64, kl: Int64, ku: Int64, val: (i: Int64, j: Int64) -> Float64) {
        this.r = r
        this.c = c
        this.kl = kl
        this.ku = ku
        this.val = val
    }
    public func dims(): (Int64, Int64) {
        return (this.r, this.c)
    }

    public func t(): Matrix {
        return Transpose(this)
    }

    public func at(i: Int64, j: Int64): Float64 {
        if (i < 0 || this.r <= i) {
            Exception("row")
        }
        if (j < 0 || this.c <= j) {
            Exception("col")
        }
        if (j < i - this.kl || i + this.ku < j) {
            return 0.0
        }
        return this.val(i, j)
    }
}

@Test
public func testBandAtSet(): Unit {
    // 2  3  4  0  0  0
    // 5  6  7  8  0  0
    // 0  9 10 11 12  0
    // 0  0 13 14 15 16
    // 0  0  0 17 18 19
    // 0  0  0  0 21 22
    let band = newBandDense(6, 6, 1, 2,
        [-1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, -1.0,
            21.0, 22.0, -1.0, -1.0])

    let (rows, cols) = band.dims()
    let (kl, ku) = band.bandwidth()

    // Explicitly test all indexes.
    let want = BandImplicit(rows, cols, kl, ku, {
        i: Int64, j: Int64 => return Float64(i * (kl + ku) + j + kl + 1)
    })
    for (i in 0..6) {
        for (j in 0..6) {
            @Expect(band.at(i, j) , want.at(i, j))
        }
    }
    // Do that same thing via a call to Equal.
    @Expect(mEqual(band, want))

    // Check At out of bounds
    for (row in [-1, rows, rows + 1]) {
        let (panicked, message) = mPanics({=> band.at(row, 0)})
        if (!panicked || message != errRowAccess) {
            throw Exception("expected panic for invalid row access N=${rows} r=${row}")
        }
    }

    for (col in [-1, cols, cols + 1]) {
        let (panicked, message) = mPanics({=> band.at(0, col)})
        if (!panicked || message != errColAccess) {
            throw Exception("expected panic for invalid column access N=${cols} c=${col}")
        }
    }

    // Check Set out of bounds
    for (row in [-1, rows, rows + 1]) {
        let (panicked, message) = mPanics({=> band.setBand(row, 0, 1.2)})
        if (!panicked || message != errRowAccess) {
            throw Exception("expected panic for invalid row access N=${rows} r=${row}")
        }
    }
    for (col in [-1, cols, cols + 1]) {
        let (panicked, message) = mPanics({=> band.setBand(0, col, 1.2)})
        if (!panicked || message != errColAccess) {
            throw Exception("expected panic for invalid column access N=${cols} c=${col}")
        }
    }

    let sts: Array<(Int64, Int64)> = [
        (0, 3),
        (0, 4),
        (0, 5),
        (1, 4),
        (1, 5),
        (2, 5),
        (2, 0),
        (3, 1),
        (4, 2),
        (5, 3)
    ]

    for (st in sts) {
        let (panicked, message) = mPanics({=> band.setBand(st[0], st[1], 1.2)})
        @Expect(!(!panicked || message != errBandSet))
    }

    let sts2: Array<(Int64, Int64, Float64, Float64)> = [(1, 2, 7.0, 15.0), (2, 3, 11.0, 15.0)]

    for (st in sts2) {
        var e = band.at(st[0], st[1])
        @Expect(e == st[2])
        band.setBand(st[0], st[1], st[3])
        e = band.at(st[0], st[1])
        @Expect(e == st[3])
    }
}

extend BandDense {
    public operator func ==(rhs: BandDense): Bool {
        var judge = true
        for (i in 0..max(this.mat.data.size, rhs.mat.data.size)) {
            if (this.mat.data[i] != rhs.mat.data[i]) {
                judge = false
            }
        }
        return (this.mat.rows == rhs.mat.rows && this.mat.cols == rhs.mat.cols && this.mat.kl == rhs.mat.kl &&
            this.mat.ku == rhs.mat.ku && this.mat.stride == rhs.mat.stride && judge)
    }

    public operator func !=(rhs: BandDense): Bool {
        var judge = false
        for (i in 0..max(this.mat.data.size, rhs.mat.data.size)) {
            if (this.mat.data[i] != rhs.mat.data[i]) {
                judge = true
            }
        }
        return (this.mat.rows == rhs.mat.rows || this.mat.cols == rhs.mat.cols || this.mat.kl == rhs.mat.kl ||
            this.mat.ku == rhs.mat.ku || this.mat.stride == rhs.mat.stride || judge)
    }
}

@Test
@Parallel
public func testBendResetAndIsEmpty(): Unit {
    for (test in [
            newBandDense(1, 1, 0, 0, [1.0]),
            newBandDense(6, 6, 1, 2,
                [-1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
                    19.0, -1.0, 21.0, 22.0, -1.0, -1.0]),
            newBandDense(6, 6, 2, 1,
                [-1.0, -1.0, 1.0, 2.0, -1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
                    17.0, 18.0, 19.0, 20.0, -1.0])
        ]) {
        @Expect(!test.isEmpty())
        test.reset()
        @Expect(test.isEmpty())
    }
}

@Test
@Parallel
public func testBandDenseTrace(): Unit {
    let data = [
        -1.0, 1.0, 2.0, 3.0,
        4.0, 5.0, 6.0, 7.0,
        8.0, 9.0, 10.0, 11.0,
        12.0, 13.0, 14.0, 15.0,
        16.0, 17.0, 18.0, -1.0,
        19.0, 20.0, -1.0, -1.0
    ]
    let b = newBandDense(6, 6, 1, 2, data)
    let expect = 1.0 + 5.0 + 9.0 + 13.0 + 17.0 + 20.0
    let tr = b.trace()
    if (!equalWithinAbsOrRel(tr, expect, 1e-14, 1e-14)) {
        throw Exception("Trace error: got ${tr}, want ${expect}")
    }
}

@Test
@Parallel
public func testBandDenseNorm(): Unit {
    let data = [
        -1.0, 1.0, 2.0, 3.0,
        4.0, 5.0, 6.0, 7.0,
        8.0, 9.0, 10.0, 11.0,
        12.0, 13.0, 14.0, 15.0,
        16.0, 17.0, 18.0, -1.0,
        19.0, 20.0, -1.0, -1.0
    ]

    let b = newBandDense(6, 6, 1, 2, data)

    let dense = newDense(6, 6, [
        1.0, 2.0, 3.0, 0.0, 0.0, 0.0,
        4.0, 5.0, 6.0, 7.0, 0.0, 0.0,
        0.0, 8.0, 9.0, 10.0, 11.0, 0.0,
        0.0, 0.0, 12.0, 13.0, 14.0, 15.0,
        0.0, 0.0, 0.0, 16.0, 17.0, 18.0,
        0.0, 0.0, 0.0, 0.0, 19.0, 20.0
    ])
    let want1 = norm(dense, 1.0)
    let got1 = b.norm(1.0)
    if (!equalWithinAbsOrRel(got1, want1, 1e-14, 1e-14)) {
        throw Exception("Norm(1) error: got ${got1}, want ${want1}")
    }

    let wantInf = norm(dense, Float64.Inf)
    let gotInf = b.norm(Float64.Inf)
    if (!equalWithinAbsOrRel(gotInf, wantInf, 1e-14, 1e-14)) {
        throw Exception("Norm(Inf) error: got ${gotInf}, want ${wantInf}")
    }

    let wantF = norm(dense, 2.0)
    let gotF = b.norm(2.0)
    if (!equalWithinAbsOrRel(gotF, wantF, 1e-14, 1e-14)) {
        throw Exception("Norm(2) error: got ${gotF}, want ${wantF}")
    }
}