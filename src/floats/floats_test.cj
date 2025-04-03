package cjnum.floats

import std.math.*
import std.random.*
import std.unittest.*
import std.unittest.testmacro.*
import std.collection.*
import cjnum.floats.scalar.equalWithinRel
import cjnum.floats.scalar.equalWithinAbsOrRel
import cjnum.floats.scalar.same as scalar_fSame

const EqTolerance = 1e-14
const Small       = 10
const Medium      = 1000
const Large       = 100000
const Huge        = 10000000

func areSlicesEqual(truth: Array<Float64>, comp: Array<Float64>): Bool {
    return fEqualApprox(truth, comp, EqTolerance)
}

func areSlicesSame(truth: Array<Float64>, comp: Array<Float64>): Bool {
    if (truth.size == comp.size) {
        for (i in 0..truth.size) {
            if (!equalWithinAbsOrRel(truth[i], comp[i], EqTolerance, EqTolerance) && !scalar_fSame(truth[i], comp[i])) {
                return false
            }
        }
    }
    return true
}

func exception(f: () -> Unit): Bool {
    try {
        f()
    } catch (e: Exception) {
        return true
    }
    return false
}

@Test
@Parallel
public func testAdd() {
    let a = [1.0, 2.0, 3.0]
    let b = [4.0, 5.0, 6.0]
    let c = [7.0, 8.0, 9.0]
    let truth = [12.0, 15.0, 18.0]
    let n = Array<Float64>(a.size, item: 0.0)

    fAdd(n, a)
    fAdd(n, b)
    fAdd(n, c)
    @Expect(areSlicesEqual(truth, n))
    fAdd(a, b)
    fAdd(a, c)
    @Expect(areSlicesEqual(truth, n))
    @Expect(exception({ => fAdd(Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0)) }))
}

@Test
@Parallel
public func testAddTo() {
    let a = [1.0, 2.0, 3.0]
    let b = [4.0, 5.0, 6.0]
    let truth = [5.0, 7.0, 9.0]
    let n1 = Array<Float64>(a.size, item: 0.0)

    let n2 = fAddTo(n1, a, b)
    @Expect(areSlicesEqual(truth, n1))
    @Expect(areSlicesEqual(truth, n2))

    @Expect(exception( { => fAddTo(Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0)) }))
    @Expect(exception( { => fAddTo(Array<Float64>(3, item: 0.0), Array<Float64>(3, item: 0.0)) }))
}

@Test
@Parallel
public func testAddConst() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let c = 6.0
    let truth = [9.0, 10.0, 7.0, 13.0, 11.0]
    fAddConst(c, s)
    @Expect(areSlicesEqual(truth, s))
}

@Test
@Parallel
public func testAddScaled() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let alpha = 6.0
    let dst = [1.0, 2.0, 3.0, 4.0, 5.0]
    let ans = [19.0, 26.0, 9.0, 46.0, 35.0]
    fAddScaled(dst, alpha, s)
    @Expect(fEqualApprox(dst, ans, EqTolerance))
    let short = [1.0]
    @Expect(exception({ => fAddScaled(dst, alpha, short) }))
    @Expect(exception({ => fAddScaled(short, alpha, s) }))
}

@Test
@Parallel
public func testAddScaledTo() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let alpha = 6.0
    let y = [1.0, 2.0, 3.0, 4.0, 5.0]
    let dst1 = Array<Float64>(5, item: 0.0)
    let ans = [19.0, 26.0, 9.0, 46.0, 35.0]
    let dst2 = fAddScaledTo(dst1, y, alpha, s)
    @Expect(fEqualApprox(dst1, ans, EqTolerance))
    @Expect(fEqualApprox(dst2, ans, EqTolerance))
    fAddScaledTo(dst1, y, alpha, s)
    @Expect(fEqualApprox(dst1, ans, EqTolerance))
    let short = [1.0]
    @Expect(exception({ => fAddScaledTo(dst1, y, alpha, short) }))
    @Expect(exception({ => fAddScaledTo(short, y, alpha, s) }))
    @Expect(exception({ => fAddScaledTo(dst1, short, alpha, s) }))
}

@Test
@Parallel
public func testArgsort() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let inds = Array<Int64>(s.size, item: 0)

    fArgsort(s, inds)

    let sortedS = [1.0, 3.0, 4.0, 5.0, 7.0]
    let trueInds = [2, 0, 1, 4, 3]

    @Expect(fEqual(s, sortedS))
    for (i in 0..trueInds.size) {
        @Expect(trueInds[i], inds[i])
    }

    let indx = [1, 2]
    @Expect(exception({ => fArgsort(s, indx) }))
}

@Test
@Parallel
public func testArgsortStable() {
    for (i in 0..=100) {
        let data = Array<Float64>(i + 2, item: 0.0)
        data[0] = 2.0
        data[data.size - 1] = 2.0
        for (j in 1..=i) {
            data[j] = 1.0
        }
        let idx = Array<Int64>(data.size, item: 0)
        fArgsortStable(data, idx)
        @Expect(isSort(data))
        @Expect(isSort(idx[..i]))
    }

    @Expect(exception({ => fArgsortStable([1.0, 2.0, 3.0], [0, 0]) }))
}

@Test
@Parallel
public func testCount() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let f = { v: Float64 => v > 3.5}
    @Expect(fCount(f, s), 3)
}

@Test
@Parallel
public func testCumProd() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let receiver = Array<Float64>(s.size, item: 0.0)
    let result = fCumProd(receiver, s)
    var truth = [3.0, 12.0, 12.0, 84.0, 420.0]
    @Expect(areSlicesEqual(truth, receiver))
    @Expect(areSlicesEqual(truth, result))
    fCumProd(receiver, s)
    @Expect(areSlicesEqual(truth, receiver))

    @Expect(exception({ => fCumProd(Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0)) }))

    let emptyReceiver = Array<Float64>()
    truth = []
    fCumProd(emptyReceiver, emptyReceiver)
    @Expect(areSlicesEqual(truth, emptyReceiver))
}

@Test
@Parallel
public func testCumSum() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let receiver = Array<Float64>(s.size, item: 0.0)
    let result = fCumSum(receiver, s)
    var truth = [3.0, 7.0, 8.0, 15.0, 20.0]
    @Expect(areSlicesEqual(truth, receiver))
    @Expect(areSlicesEqual(truth, result))
    fCumSum(receiver, s)
    @Expect(areSlicesEqual(truth, receiver))
    @Expect(exception({ => fCumSum(Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0)) }))
    let emptyReceiver = Array<Float64>()
    truth = []
    fCumSum(emptyReceiver, emptyReceiver)
    @Expect(areSlicesEqual(truth, emptyReceiver))
}

class DistanceTestStruct {
    DistanceTestStruct(
        let s!: Array<Float64> = [],
        let t!: Array<Float64> = []
    ) {}
}

@Test
@Parallel
public func testDistance() {
    let fNorms = [1.0, 2.0, 4.0, Float64.Inf]
    let tests = [
        DistanceTestStruct(),
        DistanceTestStruct(s: [8.0, 9.0, 10.0, -12.0], t: [8.0, 9.0, 10.0, -12.0]),
        DistanceTestStruct(s: [1.0, 2.0, 3.0, -4.0, -5.0, 8.0], t: [-9.2, -6.8, 9.0, -3.0, -2.0, 1.0])
    ]
    for (test in tests) {
        let tmp = Array<Float64>(test.s.size, item: 0.0)
        for (L in fNorms) {
            let dist = fDistance(L, s: test.s, t: test.t)
            copy(tmp, test.s)
            fSub(tmp, test.t)
            @Expect(dist, fNorm(tmp, L))
        }
    }
}

func copy<T>(dst: Array<T>, src: Array<T>): Int64 {
    let minSize = min(src.size, dst.size)
    dst[0..minSize] = src[0..minSize]
    return minSize
}

@Test
@Parallel
public func testDiv() {
    let s1 = [5.0, 12.0, 27.0]
    let s2 = [1.0, 2.0, 3.0]
    let ans = [5.0, 6.0, 9.0]
    fDiv(s1, s2)
    @Expect(fEqualApprox(s1, ans, EqTolerance))
    let s1Short = [1.0]
    @Expect(exception({ => fDiv(s1Short, s2) }))
    let s2Short = [1.0]
    @Expect(exception({ => fDiv(s1, s2Short) }))
}

@Test
@Parallel
public func testDivTo() {
    let s1 = [5.0, 12.0, 27.0]
    let s1orig = [5.0, 12.0, 27.0]
    let s2 = [1.0, 2.0, 3.0]
    let s2orig = [1.0, 2.0, 3.0]
    let dst1 = Array<Float64>(3, item: 0.0)
    let ans = [5.0, 6.0, 9.0]
    let dst2 = fDivTo(dst1, s1, s2)
    @Expect(fEqualApprox(dst1, ans, EqTolerance))
    @Expect(fEqualApprox(dst2, ans, EqTolerance))
    @Expect(fEqualApprox(s1, s1orig, EqTolerance))
    @Expect(fEqualApprox(s2, s2orig, EqTolerance))
    fDivTo(dst1, s1, s2)
    @Expect(fEqualApprox(dst1, ans, EqTolerance))
    let dstShort = [1.0]
    @Expect(exception{ => fDivTo(dstShort, s1, s2) })
    let s1Short = [1.0]
    @Expect(exception{ => fDivTo(dst1, s1Short, s2) })
    let s2Short = [1.0]
    @Expect(exception{ => fDivTo(dst1, s1, s2Short) })
}

@Test
@Parallel
public func testDot() {
    let s1 = [1.0, 2.0, 3.0, 4.0]
    let s2 = [-3.0, 4.0, 5.0, -6.0]
    let truth = -4.0
    let ans = fDot(s1, s2)
    @Expect(ans, truth)
    @Expect(exception({ => fDot(Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0)) }))
}

@Test
@Parallel
public func testEquals() {
    var s1 = [1.0, 2.0, 3.0, 4.0]
    var s2 = [1.0, 2.0, 3.0, 4.0]
    @Expect(fEqual(s1, s2))
    s2 = [1.0, 2.0, 3.0, 4.0 + 1e-14]
    @Expect(!fEqual(s1, s2))
    @Expect(!fEqual(s1, []))
}

@Test
@Parallel
public func testEqualApprox() {
    var s1 = [1.0, 2.0, 3.0, 4.0]
    var s2 = [1.0, 2.0, 3.0, 4.0 + 1e-10]
    @Expect(!fEqualApprox(s1, s2, 1e-13))
    @Expect(fEqualApprox(s1, s2, 1e-5))
    s1 = [1.0, 2.0, 3.0, 1000.0]
    s2 = [1.0, 2.0, 3.0, 1000.0 * (1.0 + 1e-7)]
    @Expect(!fEqualApprox(s1, s2, 1e-8))
    @Expect(fEqualApprox(s1, s2, 1e-5))
    @Expect(!fEqualApprox(s1, [], 1e-5))
}

@Test
@Parallel
public func testEqualFunc() {
    let s1 = [1.0, 2.0, 3.0, 4.0]
    var s2 = [1.0, 2.0, 3.0, 4.0]
    let eq = { x: Float64, y: Float64 => x == y }
    @Expect(fEqualFunc(s1, s2, eq))
    s2 = [1.0, 2.0, 3.0, 4.0 + 1e-14]
    @Expect(!fEqualFunc(s1, s2, eq))
    @Expect(!fEqualFunc(s1, [], eq))
}

const smallestNonzeroFloat64 = 4.9406564584124654e-324

class fEqualsRelativeTestStruct {
    public fEqualsRelativeTestStruct(
        let a: Float64,
        let b: Float64,
        var tol: Float64,
        let fEqual: Bool
    ) {}
}

@Test
@Parallel
public func testEqualsRelative() {
    let fEqualityTests = [
        fEqualsRelativeTestStruct(1000000.0, 1000001.0, 0.0, true),
        fEqualsRelativeTestStruct(1000001.0, 1000000.0, 0.0, true),
        fEqualsRelativeTestStruct(10000.0, 10001.0, 0.0, false),
        fEqualsRelativeTestStruct(10001.0, 10000.0, 0.0, false),
        fEqualsRelativeTestStruct(-1000000.0, -1000001.0, 0.0, true),
        fEqualsRelativeTestStruct(-1000001.0, -1000000.0, 0.0, true),
        fEqualsRelativeTestStruct(-10000.0, -10001.0, 0.0, false),
        fEqualsRelativeTestStruct(-10001.0, -10000.0, 0.0, false),
        fEqualsRelativeTestStruct(1.0000001, 1.0000002, 0.0, true),
        fEqualsRelativeTestStruct(1.0000002, 1.0000001, 0.0, true),
        fEqualsRelativeTestStruct(1.0002, 1.0001, 0.0, false),
        fEqualsRelativeTestStruct(1.0001, 1.0002, 0.0, false),
        fEqualsRelativeTestStruct(-1.000001, -1.000002, 0.0, true),
        fEqualsRelativeTestStruct(-1.000002, -1.000001, 0.0, true),
        fEqualsRelativeTestStruct(-1.0001, -1.0002, 0.0, false),
        fEqualsRelativeTestStruct(-1.0002, -1.0001, 0.0, false),
        fEqualsRelativeTestStruct(0.000000001000001, 0.000000001000002, 0.0, true),
        fEqualsRelativeTestStruct(0.000000001000002, 0.000000001000001, 0.0, true),
        fEqualsRelativeTestStruct(0.000000000001002, 0.000000000001001, 0.0, false),
        fEqualsRelativeTestStruct(0.000000000001001, 0.000000000001002, 0.0, false),
        fEqualsRelativeTestStruct(-0.000000001000001, -0.000000001000002, 0.0, true),
        fEqualsRelativeTestStruct(-0.000000001000002, -0.000000001000001, 0.0, true),
        fEqualsRelativeTestStruct(-0.000000000001002, -0.000000000001001, 0.0, false),
        fEqualsRelativeTestStruct(-0.000000000001001, -0.000000000001002, 0.0, false),
        fEqualsRelativeTestStruct(0.0, 0.0, 0.0, true),
        fEqualsRelativeTestStruct(0.0, -0.0, 0.0, true),
        fEqualsRelativeTestStruct(-0.0, -0.0, 0.0, true),
        fEqualsRelativeTestStruct(0.00000001, 0.0, 0.0, false),
        fEqualsRelativeTestStruct(0.0, 0.00000001, 0.0, false),
        fEqualsRelativeTestStruct(-0.00000001, 0.0, 0.0, false),
        fEqualsRelativeTestStruct(0.0, -0.00000001, 0.0, false),
        fEqualsRelativeTestStruct(0.0, 1e-310, 0.01, true),
        fEqualsRelativeTestStruct(1e-310, 0.0, 0.01, true),
        fEqualsRelativeTestStruct(1e-310, 0.0, 0.000001, false),
        fEqualsRelativeTestStruct(0.0, 1e-310, 0.000001, false),
        fEqualsRelativeTestStruct(0.0, -1e-310, 0.1, true),
        fEqualsRelativeTestStruct(-1e-310, 0.0, 0.1, true),
        fEqualsRelativeTestStruct(-1e-310, 0.0, 0.00000001, false),
        fEqualsRelativeTestStruct(0.0, -1e-310, 0.00000001, false),
        fEqualsRelativeTestStruct(Float64.Inf, Float64.Inf, 0.0, true),
        fEqualsRelativeTestStruct(-Float64.Inf, -Float64.Inf, 0.0, true),
        fEqualsRelativeTestStruct(-Float64.Inf, Float64.Inf, 0.0, false),
        fEqualsRelativeTestStruct(Float64.Inf, Float64.Max, 0.0, false),
        fEqualsRelativeTestStruct(-Float64.Inf, -Float64.Max, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, 0.0, 0.0, false),
        fEqualsRelativeTestStruct(-0.0, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, -0.0, 0.0, false),
        fEqualsRelativeTestStruct(0.0, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, Float64.Inf, 0.0, false),
        fEqualsRelativeTestStruct(Float64.Inf, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, -Float64.Inf, 0.0, false),
        fEqualsRelativeTestStruct(-Float64.Inf, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, Float64.Max, 0.0, false),
        fEqualsRelativeTestStruct(Float64.Max, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, -Float64.Max, 0.0, false),
        fEqualsRelativeTestStruct(-Float64.Max, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, smallestNonzeroFloat64, 0.0, false),
        fEqualsRelativeTestStruct(smallestNonzeroFloat64, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(Float64.NaN, -smallestNonzeroFloat64, 0.0, false),
        fEqualsRelativeTestStruct(-smallestNonzeroFloat64, Float64.NaN, 0.0, false),
        fEqualsRelativeTestStruct(1.000000001, -1.0, 0.0, false),
        fEqualsRelativeTestStruct(-1.0, 1.000000001, 0.0, false),
        fEqualsRelativeTestStruct(-1.000000001, 1.0, 0.0, false),
        fEqualsRelativeTestStruct(1.0, -1.000000001, 0.0, false),
        fEqualsRelativeTestStruct(10.0 * smallestNonzeroFloat64, 10.0 * (-smallestNonzeroFloat64), 0.0, true),
        fEqualsRelativeTestStruct(1e11 * smallestNonzeroFloat64, 1e11 * (-smallestNonzeroFloat64), 0.0, false),
        fEqualsRelativeTestStruct(smallestNonzeroFloat64, -smallestNonzeroFloat64, 0.0, true),
        fEqualsRelativeTestStruct(-smallestNonzeroFloat64, smallestNonzeroFloat64, 0.0, true),
        fEqualsRelativeTestStruct(smallestNonzeroFloat64, 0.0, 0.0, true),
        fEqualsRelativeTestStruct(0.0, smallestNonzeroFloat64, 0.0, true),
        fEqualsRelativeTestStruct(-smallestNonzeroFloat64, 0.0, 0.0, true),
        fEqualsRelativeTestStruct(0.0, -smallestNonzeroFloat64, 0.0, true),
        fEqualsRelativeTestStruct(0.000000001, -smallestNonzeroFloat64, 0.0, false),
        fEqualsRelativeTestStruct(0.000000001, smallestNonzeroFloat64, 0.0, false),
        fEqualsRelativeTestStruct(smallestNonzeroFloat64, 0.000000001, 0.0, false),
        fEqualsRelativeTestStruct(-smallestNonzeroFloat64, 0.000000001, 0.0, false)
    ]

    for (ts in fEqualityTests) {
        if (ts.tol == 0.0) {
            ts.tol = 1e-5
        }
        @Expect(equalWithinRel(ts.a, ts.b, ts.tol) == ts.fEqual)
    }
}

@Test
@Parallel
public func testEqualLengths() {
    let s1 = [1.0, 2.0, 3.0, 4.0]
    let s2 = [1.0, 2.0, 3.0, 4.0]
    let s3 = [1.0, 2.0, 3.0]
    @Expect(fEqualLengths(s1, s2))
    @Expect(!fEqualLengths(s1, s3))
    @Expect(fEqualLengths(s1))
    @Expect(fEqualLengths())
}

func eqIntSlice(one: Array<Int64>, two: Array<Int64>): String {
    if (one.size != two.size) {
        return "Length mismatch"
    }
    for (i in 0..one.size) {
        if (one[i] != two[i]) {
            return "Index " + i.toString() + " mismatch"
        }
    }
    return ""
}

@Test
@Parallel
public func testFind() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let f = { v: Float64 => v > 3.5}
    let allTrueInds = [1, 3, 4]
    
    var inds = fFind(f, s, 2)
    let trueInds = allTrueInds[0..2]
    @Expect(eqIntSlice(inds.toArray(), trueInds), "")

    inds = ArrayList<Int64>([1, 2, 3, 4, 5, 6])
    inds = fFind(f, s, 0, inds: inds)
    @Expect(eqIntSlice(inds.toArray(), []), "")

    inds = ArrayList<Int64>([1, 2, 3, 4, 5, 6])
    inds = fFind(f, s, 2, inds: inds)
    @Expect(eqIntSlice(inds.toArray(), trueInds), "")

    try {
        inds = fFind(f, s, 4, inds: inds)
        @Expect(false)
    } catch (e: Exception) {
        // fall
    }
    @Expect(eqIntSlice(inds.toArray(), allTrueInds), "")
    
    inds = fFind(f, s, -1)
    @Expect(eqIntSlice(inds.toArray(), allTrueInds), "")
}

class HasNaNTestStruct {
    HasNaNTestStruct(
        let s!: Array<Float64> = [],
        let ans!: Bool = false
    ) {}
}

@Test
@Parallel
public func testHasNaN() {
    let tests = [
        HasNaNTestStruct(),
        HasNaNTestStruct(s: [1.0, 2.0, 3.0, 4.0]),
        HasNaNTestStruct(s: [1.0, Float64.NaN, 3.0, 4.0], ans: true),
        HasNaNTestStruct(s: [1.0, 2.0, 3.0, Float64.NaN], ans: true)
    ]
    for (test in tests) {
        @Expect(fHasNaN(test.s), test.ans)
    }
}

@Test
@Parallel
public func testLogSpan() {
    let receiver1 = Array<Float64>(6, item: 0.0)
    let truth = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
    let receiver2 = fLogSpan(0.001, 100.0, dst: receiver1)
    let tst = Array<Float64>(6, item: 0.0)
    for (i in 0..truth.size) {
        tst[i] = receiver1[i] / truth[i]
    }
    let comp = Array<Float64>(6, item: 1.0)
    @Expect(areSlicesEqual(comp, tst))
    for (i in 0..truth.size) {
        tst[i] = receiver2[i] / truth[i]
    }
    @Expect(areSlicesEqual(comp, tst))
    @Expect(exception({ => fLogSpan(1.0, 5.0)}))
    @Expect(exception({ => fLogSpan(1.0, 5.0, dst: Array<Float64>(1, item: 0.0)) }))
}

@Test
@Parallel
public func testLogSumExp() {
    var s = [1.0, 2.0, 3.0, 4.0, 5.0]
    var truth = 5.4519143959375933331957225109748087179338972737576824
    var val = fLogSumExp(s)
    @Expect(abs(val - truth) <= EqTolerance)

    s = [1.0, 2.0]
    truth = 2.3132616875182228340489954949678556419152800856703483
    val = fLogSumExp(s)
    @Expect(abs(val - truth) <= EqTolerance)

    s = [-1001.0, -1002.0, -1003.0, -1004.0, -1005.0]
    truth = -1000.54808560406240666680427748902519128206610272624
    val = fLogSumExp(s)
    @Expect(abs(val - truth) <= EqTolerance)

    s = [1.0, 2.0, 3.0, 4.0, 5.0, Float64.Inf]
    truth = Float64.Inf
    val = fLogSumExp(s)
    @Expect(val, truth)

    s = [1.0, 2.0, 3.0, 4.0, 5.0, -Float64.Inf]
    truth = 5.4519143959375933331957225109748087179338972737576824
    val = fLogSumExp(s)
    @Expect(abs(val - truth) <= EqTolerance)
}

class IdxTestStruct {
    IdxTestStruct(
        let in_!: Array<Float64>,
        let wantIdx!: Int64,
        let wantVal!: Float64
    ) {}
}

@Test
@Parallel
public func testMaxAndIdx() {
    let tests = [
        IdxTestStruct(in_: [3.0, 4.0, 1.0, 7.0, 5.0], wantIdx: 3, wantVal: 7.0),
        IdxTestStruct(in_: [Float64.NaN, 4.0, 1.0, 7.0, 5.0], wantIdx: 3, wantVal: 7.0),
        IdxTestStruct(in_: [Float64.NaN, Float64.NaN, Float64.NaN], wantIdx: 0, wantVal: Float64.NaN),
        IdxTestStruct(in_: [Float64.NaN, -Float64.Inf], wantIdx: 1, wantVal: -Float64.Inf),
        IdxTestStruct(in_: [Float64.NaN, Float64.Inf], wantIdx: 1, wantVal: Float64.Inf)
    ]

    for (test in tests) {
        @Expect(fMaxIdx(test.in_), test.wantIdx)
        @Expect(scalar_fSame(fMax(test.in_), test.wantVal))
    }
    @Expect(exception({ => fMaxIdx([]) }))
}

@Test
@Parallel
public func testMinAndIdx() {
    let tests = [
        IdxTestStruct(in_: [3.0, 4.0, 1.0, 7.0, 5.0], wantIdx: 2, wantVal: 1.0),
        IdxTestStruct(in_: [Float64.NaN, 4.0, 1.0, 7.0, 5.0], wantIdx: 2, wantVal: 1.0),
        IdxTestStruct(in_: [Float64.NaN, Float64.NaN, Float64.NaN], wantIdx: 0, wantVal: Float64.NaN),
        IdxTestStruct(in_: [Float64.NaN, -Float64.Inf], wantIdx: 1, wantVal: -Float64.Inf),
        IdxTestStruct(in_: [Float64.NaN, Float64.Inf], wantIdx: 1, wantVal: Float64.Inf)
    ]

    for (test in tests) {
        @Expect(fMinIdx(test.in_), test.wantIdx)
        @Expect(scalar_fSame(fMin(test.in_), test.wantVal))
    }
    @Expect(exception({ => fMinIdx([]) }))
}

@Test
@Parallel
public func testMul() {
    let s1 = [1.0, 2.0, 3.0]
    let s2 = [1.0, 2.0, 3.0]
    let ans = [1.0, 4.0, 9.0]
    fMul(s1, s2)
    @Expect(fEqualApprox(s1, ans, EqTolerance))
    let s1Short = [1.0]
    @Expect(exception({ => fMul(s1Short, s2) }))
    let s2Short = [1.0]
    @Expect(exception({ => fMul(s1, s2Short) }))
}

@Test
@Parallel
public func testMulTo() {
    let s1 = [1.0, 2.0, 3.0]
    let s1orig = [1.0, 2.0, 3.0]
    let s2 = [1.0, 2.0, 3.0]
    let s2orig = [1.0, 2.0, 3.0]
    let dst1 = Array<Float64>(3, item: 0.0)
    let ans = [1.0, 4.0, 9.0]
    let dst2 = fMulTo(dst1, s1, s2)
    @Expect(fEqualApprox(dst1, ans, EqTolerance))
    @Expect(fEqualApprox(dst2, ans, EqTolerance))
    @Expect(fEqualApprox(s1, s1orig, EqTolerance))
    @Expect(fEqualApprox(s2, s2orig, EqTolerance))
    fMulTo(dst1, s1, s2)
    @Expect(fEqualApprox(dst1, ans, EqTolerance))
    let dstShort = [1.0]
    @Expect(exception({ => fMulTo(dstShort, s1, s2) }))
    let s1Short = [1.0]
    @Expect(exception({ => fMulTo(dst1, s1Short, s2) }))
    let s2Short = [1.0]
    @Expect(exception({ => fMulTo(dst1, s1, s2Short) }))
}

class NearestIdxTestStruct {
    NearestIdxTestStruct(
        let in_!: Array<Float64>,
        let query!: Float64,
        let want!: Int64
    ) {}
}

@Test
@Parallel
public func testNearestIdx() {
    let tests = [
        NearestIdxTestStruct(in_: [6.2, 3.0, 5.0, 6.2, 8.0], query: 2.0, want: 1),
        NearestIdxTestStruct(in_: [6.2, 3.0, 5.0, 6.2, 8.0], query: 9.0, want: 4),
        NearestIdxTestStruct(in_: [6.2, 3.0, 5.0, 6.2, 8.0], query: 3.1, want: 1),
        NearestIdxTestStruct(in_: [6.2, 3.0, 5.0, 6.2, 8.0], query: 2.9, want: 1),
        NearestIdxTestStruct(in_: [6.2, 3.0, 5.0, 6.2, 8.0], query: 3.0, want: 1),
        NearestIdxTestStruct(in_: [6.2, 3.0, 5.0, 6.2, 8.0], query: 6.2, want: 0),
        NearestIdxTestStruct(in_: [6.2, 3.0, 5.0, 6.2, 8.0], query: 4.0, want: 1),
        NearestIdxTestStruct(in_: [Float64.NaN, 3.0, 2.0, -1.0], query: 2.0, want: 2),
        NearestIdxTestStruct(in_: [0.0, Float64.NaN, -1.0, 2.0], query: Float64.NaN, want: 0),
        NearestIdxTestStruct(in_: [0.0, Float64.NaN, -1.0, 2.0], query: Float64.Inf, want: 3),
        NearestIdxTestStruct(in_: [0.0, Float64.NaN, -1.0, 2.0], query: -Float64.Inf, want: 2),
        NearestIdxTestStruct(in_: [Float64.NaN, Float64.NaN, Float64.NaN], query: 1.0, want: 0),
        NearestIdxTestStruct(in_: [Float64.NaN, Float64.Inf], query: 1.0, want: 1)
    ]

    for (test in tests) {
        @Expect(fNearestIdx(test.in_, test.query), test.want)
    }
    @Expect(exception({ => fNearestIdx([], 0.0)}))
}

class NearestIdxForSpanTestStruct {
    NearestIdxForSpanTestStruct(
        let length!: Int64,
        let lower!: Float64,
        let upper!: Float64,
        let value!: Float64,
        let idx!: Int64
    ) {}
}

@Test
@Parallel
public func testNearestIdxForSpan() {
    let tests = [
        NearestIdxForSpanTestStruct(length: 13, lower: 7.0, upper: 8.2, value: 6.0, idx: 0),
        NearestIdxForSpanTestStruct(length: 13, lower: 7.0, upper: 8.2, value: 10.0, idx: 12),
        NearestIdxForSpanTestStruct(length: 13, lower: 7.0, upper: 8.2, value: 7.19, idx: 2),
        NearestIdxForSpanTestStruct(length: 13, lower: 7.0, upper: 8.2, value: 7.21, idx: 2),
        NearestIdxForSpanTestStruct(length: 13, lower: 7.0, upper: 8.2, value: 7.2, idx: 2),
        NearestIdxForSpanTestStruct(length: 13, lower: 7.0, upper: 8.2, value: 7.151, idx: 2),
        NearestIdxForSpanTestStruct(length: 13, lower: 7.0, upper: 8.2, value: 7.249, idx: 2),
        NearestIdxForSpanTestStruct(length: 4, lower: -Float64.Inf, upper: Float64.Inf, value: copysign(0.0, -1.0), idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: -Float64.Inf, upper: Float64.Inf, value: 0.0, idx: 2),
        NearestIdxForSpanTestStruct(length: 5, lower: -Float64.Inf, upper: Float64.Inf, value: -Float64.Inf, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: -Float64.Inf, upper: Float64.Inf, value: Float64.Inf, idx: 3),
        NearestIdxForSpanTestStruct(length: 4, lower: -Float64.Inf, upper: Float64.Inf, value: 0.0, idx: 2),
        NearestIdxForSpanTestStruct(length: 4, lower: -Float64.Inf, upper: Float64.Inf, value: Float64.Inf, idx: 2),
        NearestIdxForSpanTestStruct(length: 4, lower: -Float64.Inf, upper: Float64.Inf, value: -Float64.Inf, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: Float64.Inf, upper: Float64.Inf, value: 1.0, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: Float64.NaN, upper: Float64.NaN, value: 1.0, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: 0.0, upper: 1.0, value: Float64.NaN, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: Float64.NaN, upper: 1.0, value: 0.0, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: Float64.NaN, upper: 1.0, value: 0.0, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: -Float64.Inf, upper: 1.0, value: -Float64.Inf, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: -Float64.Inf, upper: 1.0, value: 0.0, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: Float64.Inf, upper: 1.0, value: Float64.Inf, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: Float64.Inf, upper: 1.0, value: 0.0, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: 100.0, upper: -Float64.Inf, value: -Float64.Inf, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: 100.0, upper: -Float64.Inf, value: 200.0, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: 100.0, upper: Float64.Inf, value: Float64.Inf, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: 100.0, upper: Float64.Inf, value: 200.0, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: -1.0, upper: 2.0, value: -Float64.Inf, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: -1.0, upper: 2.0, value: Float64.Inf, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: 1.0, upper: -2.0, value: -Float64.Inf, idx: 4),
        NearestIdxForSpanTestStruct(length: 5, lower: 1.0, upper: -2.0, value: Float64.Inf, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: 2.0, upper: 0.0, value: 3.0, idx: 0),
        NearestIdxForSpanTestStruct(length: 5, lower: 2.0, upper: 0.0, value: -1.0, idx: 4)
    ]
    for (test in tests) {
        @Expect(fNearestIdxForSpan(test.length, test.lower, test.upper, test.value) == test.idx)
    }
    @Expect(exception({ => fNearestIdxForSpan(1, 0.0, 1.0, 0.5) }))
}

@Test
@Parallel
public func testNorm() {
    let s = [-1.0, -3.4, 5.0, -6.0]
    var val = fNorm(s, Float64.Inf)
    var truth = 6.0
    @Expect(abs(val - truth) <= EqTolerance)

    val = fNorm(s, 2.0)
    truth = 8.5767126569566267590651614132751986658027271236078592
    @Expect(abs(val - truth) <= EqTolerance)

    val = fNorm(s, 3.0)
    truth = 7.2514321388020228478109121239004816430071237369356233
    @Expect(abs(val - truth) <= EqTolerance)

    val = fNorm(s, 1.0)
    truth = 15.4
    @Expect(abs(val - truth) <= EqTolerance)
}

@Test
@Parallel
public func testProd() {
    var s: Array<Float64> = []
    @Expect(fProd(s), 1.0)

    s = [3.0, 4.0, 1.0, 7.0, 5.0]
    @Expect(fProd(s), 420.0)
}

@Test
@Parallel
public func testReverse() {
    let ss = [
        [0.0],
        [1.0, 0.0],
        [2.0, 1.0, 0.0],
        [3.0, 2.0, 1.0, 0.0],
        [9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]
    ]
    for (s in ss) {
        fReverse(s)
        for (i in 0..s.size) {
            @Expect(s[i], Float64(i))
        }
    }
}

@Test
@Parallel
public func testScale() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let c = 5.0
    let truth = [15.0, 20.0, 5.0, 35.0, 25.0]
    fScale(c, s)
    @Expect(areSlicesEqual(truth, s))
}

@Test
@Parallel
public func testScaleTo() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let sCopy = s.clone()
    let c = 5.0
    let truth = [15.0, 20.0, 5.0, 35.0, 25.0]
    let dst = Array<Float64>(s.size, item: 0.0)
    fScaleTo(dst, c, s)
    @Expect(fSame(dst, truth))
    @Expect(fSame(s, sCopy))
    @Expect(exception({ => fScaleTo(dst, 0.0, [1.0]) }))
}

@Test
@Parallel
public func testSame() {
    var s1 = [1.0, 2.0, 3.0, 4.0]
    var s2 = [1.0, 2.0, 3.0, 4.0]
    @Expect(fSame(s1, s2))
    s2 = [1.0, 2.0, 3.0, 4.0 + 1e-14]
    @Expect(!fSame(s1, s2))
    @Expect(!fSame(s1, []))
    s1 = [1.0, 2.0, Float64.NaN, 4.0]
    s2 = [1.0, 2.0, Float64.NaN, 4.0]
    @Expect(fSame(s1, s2))
    s1 = [1.0, 2.0, Float64.NaN, 4.0]
    s2 = [1.0, Float64.NaN, 3.0, 4.0]
    @Expect(!fSame(s1, s2))
} 

class fSpanTestStruct {
    fSpanTestStruct(
        let n!: Int64,
        let l!: Float64,
        let u!: Float64,
        let want!: Array<Float64>
    ) {}
}

@Test
@Parallel
public func testSpan() {
    var receiver1 = Array<Float64>(5, item: 0.0)
    var truth = [1.0, 2.0, 3.0, 4.0, 5.0]
    let receiver2 = fSpan(receiver1, 1.0, 5.0)
    @Expect(areSlicesEqual(truth, receiver1))
    @Expect(areSlicesEqual(truth, receiver2))
    receiver1 = Array<Float64>(6, item: 0.0)
    truth = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    fSpan(receiver1, 0.0, 1.0)
    @Expect(areSlicesEqual(truth, receiver1))
    @Expect(exception({ => fSpan(Array<Float64>(1, item: 0.0), 1.0, 5.0) }))

    let tests = [
        fSpanTestStruct(n: 4, l: -Float64.Inf, u: Float64.Inf, want: [-Float64.Inf, -Float64.Inf, Float64.Inf, Float64.Inf]),
        fSpanTestStruct(n: 4, l: Float64.Inf, u: -Float64.Inf, want: [Float64.Inf, Float64.Inf, -Float64.Inf, -Float64.Inf]),
        fSpanTestStruct(n: 5, l: -Float64.Inf, u: Float64.Inf, want: [-Float64.Inf, -Float64.Inf, 0.0, Float64.Inf, Float64.Inf]),
        fSpanTestStruct(n: 5, l: Float64.Inf, u: -Float64.Inf, want: [Float64.Inf, Float64.Inf, 0.0, -Float64.Inf, -Float64.Inf]),
        fSpanTestStruct(n: 5, l: Float64.Inf, u: Float64.Inf, want: [Float64.Inf, Float64.Inf, Float64.Inf, Float64.Inf, Float64.Inf]),
        fSpanTestStruct(n: 5, l: -Float64.Inf, u: -Float64.Inf, want: [-Float64.Inf, -Float64.Inf, -Float64.Inf, -Float64.Inf, -Float64.Inf]),
        fSpanTestStruct(n: 5, l: -Float64.Inf, u: Float64.NaN, want: [-Float64.Inf, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]),
        fSpanTestStruct(n: 5, l: Float64.Inf, u: Float64.NaN, want: [Float64.Inf, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]),
        fSpanTestStruct(n: 5, l: Float64.NaN, u: -Float64.Inf, want: [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, -Float64.Inf]),
        fSpanTestStruct(n: 5, l: Float64.NaN, u: Float64.Inf, want: [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, Float64.Inf]),
        fSpanTestStruct(n: 5, l: 42.0, u: -Float64.Inf, want: [42.0, -Float64.Inf, -Float64.Inf, -Float64.Inf, -Float64.Inf]),
        fSpanTestStruct(n: 5, l: 42.0, u: Float64.Inf, want: [42.0, Float64.Inf, Float64.Inf, Float64.Inf, Float64.Inf]),
        fSpanTestStruct(n: 5, l: 42.0, u: Float64.NaN, want: [42.0, Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN]),
        fSpanTestStruct(n: 5, l: -Float64.Inf, u: 42.0, want: [-Float64.Inf, -Float64.Inf, -Float64.Inf, -Float64.Inf, 42.0]),
        fSpanTestStruct(n: 5, l: Float64.Inf, u: 42.0, want: [Float64.Inf, Float64.Inf, Float64.Inf, Float64.Inf, 42.0]),
        fSpanTestStruct(n: 5, l: Float64.NaN, u: 42.0, want: [Float64.NaN, Float64.NaN, Float64.NaN, Float64.NaN, 42.0])
    ]
    for (test in tests) {
        let got = fSpan(Array<Float64>(test.n, item: 0.0), test.l, test.u)
        @Expect(areSlicesSame(test.want, got))
    }
}

@Test
@Parallel
public func testSub() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let v = [1.0, 2.0, 3.0, 4.0, 5.0]
    var truth = [2.0, 2.0, -2.0, 3.0, 0.0]
    fSub(s, v)
    @Expect(areSlicesEqual(truth, s))
    @Expect(exception({ => fSub(Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0)) }))
}

@Test
@Parallel
public func testSubTo() {
    let s = [3.0, 4.0, 1.0, 7.0, 5.0]
    let v = [1.0, 2.0, 3.0, 4.0, 5.0]
    var truth = [2.0, 2.0, -2.0, 3.0, 0.0]
    var dst1 = Array<Float64>(s.size, item: 0.0)
    var dst2 = fSubTo(dst1, s, v)
    @Expect(areSlicesEqual(truth, dst1))
    @Expect(areSlicesEqual(truth, dst2))
    @Expect(exception({ => fSubTo(Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0), Array<Float64>(3, item: 0.0)) }))
    @Expect(exception({ => fSubTo(Array<Float64>(3, item: 0.0), Array<Float64>(2, item: 0.0), Array<Float64>(3, item: 0.0)) }))
    @Expect(exception({ => fSubTo(Array<Float64>(3, item: 0.0), Array<Float64>(3, item: 0.0), Array<Float64>(2, item: 0.0)) }))
}

@Test
@Parallel
public func testSum() {
    var s: Array<Float64> = []
    @Expect(fSum(s), 0.0)

    s = [3.0, 4.0, 1.0, 7.0, 5.0]
    @Expect(fSum(s), 20.0)
}

class WithInTestStruct {
    WithInTestStruct(
        let s!: Array<Float64> = [],
        let v!: Float64 = 0.0,
        let idx!: Int64 = 0,
        let exception!: Bool = false
    ) {}
}

@Test
@Parallel
public func testWithIn() {
    let tests = [
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: 1.0, idx: 0),
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: 9.0, idx: -1),
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: 1.5, idx: 0),
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: 2.0, idx: 1),
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: 2.5, idx: 1),
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: -3.0, idx: -1),
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: 15.0, idx: -1),
        WithInTestStruct(s: [1.0, 2.0, 5.0, 9.0], v: Float64.NaN, idx: -1),
        WithInTestStruct(s: [5.0, 2.0, 6.0], exception: true),
        WithInTestStruct(exception: true),
        WithInTestStruct(s: [1.0], exception: true)
    ]
    for (test in tests) {
        var idx: Int64 = 0
        try {
            idx = fWithin(test.s, test.v)
        } catch (e: Exception) {
            @Expect(test.exception) 
            continue
        }
        if (test.exception) {
            @Expect(false)
            continue
        }
        @Expect(idx, test.idx)
    }
}

class SumCompensatedTestStruct {
    SumCompensatedTestStruct(
        let s!: Array<Float64>,
        let want!: Float64
    ) {}
}

@Test
@Parallel
public func testSumCompensated() {
    let k = 100000
    let s1 = Array<Float64>(2 * k + 1, item: 0.0)
    for (i in -k..=k) {
        s1[i + k] = 0.2 * Float64(i)
    }
    let s2 = Array<Float64>(k + 1, item: 0.0)
    for (i in 0..k) {
        s2[i] = 10.0 / Float64(k)
    }
    s2[k] = -10.0

    let tests = [
        SumCompensatedTestStruct(s: s1, want: 0.0),
        SumCompensatedTestStruct(s: s2, want: 0.0),
        SumCompensatedTestStruct(s: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], want: 55.0),
        SumCompensatedTestStruct(s: [1.2e20, 0.1, -2.4e20, -0.1, 1.2e20, 0.2, 0.2], want: 0.4),
        SumCompensatedTestStruct(s: [1.0, 1e100, 1.0, -1e100], want: 2.0)
    ]
    for (test in tests) {
        let got = fSumCompensated(test.s)
        @Expect(abs(got - test.want) <= EqTolerance, true)
    }
}
