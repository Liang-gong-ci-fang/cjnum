package cjnum.floats

import std.math.*
import std.sort.*
import std.collection.*
import cjnum.floats.scalar.*
import cjnum.internal.asm.f64.*
import cjnum.internal.asm.f64.cumProd as f64_cumProd
import cjnum.internal.asm.f64.cumSum as f64_cumSum
import cjnum.internal.asm.f64.div as f64_div
import cjnum.internal.asm.f64.divTo as f64_divTo
import cjnum.internal.asm.f64.sum as f64_sum
import cjnum.internal.asm.f64.addConst as f64_addConst

const zeroLength   = "floats: zero length slice"
const shortSpan    = "floats: slice length less than 2"
const badLength    = "floats: slice lengths do not match"
const badDstLength = "floats: destination slice length does not match input"

public func add(dst: Array<Float64>, s: Array<Float64>): Unit {
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    axpyUnitaryTo(dst, 1.0, s, dst)
}

public func addTo(dst: Array<Float64>, s: Array<Float64>, t: Array<Float64>): Array<Float64> {
    if (s.size != t.size) {
        throw Exception(badLength)
    }
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    axpyUnitaryTo(dst, 1.0, s, t)
    return dst
}

public func addConst(c: Float64, dst: Array<Float64>): Unit {
    f64_addConst(c, dst)
}

public func addScaled(dst: Array<Float64>, alpha: Float64, s: Array<Float64>): Unit{
    if (dst.size != s.size) {
        throw Exception(badLength)
    }
    axpyUnitaryTo(dst, alpha, s, dst)
}

public func addScaledTo(dst: Array<Float64>, y: Array<Float64>, alpha: Float64, s: Array<Float64>): Array<Float64> {
    if (s.size != y.size) {
        throw Exception(badLength)
    }
    if (dst.size != y.size) {
        throw Exception(badDstLength)
    }
    axpyUnitaryTo(dst, alpha, s, y)
    return dst
}

class Argsort{ 
    let s_inds: Array<(Float64, Int64)>

    Argsort(let s!: Array<Float64>, let inds!: Array<Int64>) {
        if (s.size != inds.size) {
            throw Exception(badDstLength)
        }
        this.s_inds = Array<(Float64, Int64)>(s.size, item: (0.0, 0))
        for (i in 0..s.size) {
            this.s_inds[i] = (s[i], inds[i])
        }
    }

    public func len(): Int64 {
        return this.s_inds.size
    }

    public func sort(stable!: Bool = false): Unit {
        this.s_inds.sortBy(stable: stable, comparator: { ts, ta => 
            if (ts[0] == ta[0]) {
                Ordering.EQ
            } else if (ts[0] > ta[0]) {
                Ordering.GT
            } else {
                Ordering.LT
            }
        })
        for (i in 0..s.size) {
            this.s[i] = s_inds[i][0]
            this.inds[i] = s_inds[i][1]
        }
    }
}

public func argsort(dst: Array<Float64>, inds: Array<Int64>): Unit { 
    if (dst.size != inds.size) {
        throw Exception(badDstLength)
    }
    for (i in 0..dst.size) {
        inds[i] = i
    }

    let a = Argsort(s: dst, inds: inds)
    a.sort()
}

public func argsortStable(dst: Array<Float64>, inds: Array<Int64>): Unit { 
    if (dst.size != inds.size) {
        throw Exception(badDstLength)
    }
    for (i in 0..dst.size) {
        inds[i] = i
    }

    let a = Argsort(s: dst, inds: inds)
    a.sort(stable: true)
}

public func count(f: (Float64) -> Bool, s: Array<Float64>): Int64 {
    var n: Int64 = 0
    for (val in s) {
        if (f(val)) {
            n++
        }
    }
    return n
}

public func cumProd(dst: Array<Float64>, s: Array<Float64>): Array<Float64> {
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    if (dst.size == 0) {
        return dst
    }
    return f64_cumProd(dst, s)
}

public func cumSum(dst: Array<Float64>, s: Array<Float64>): Array<Float64> {
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    if (dst.size == 0) {
        return dst
    }
    return f64_cumSum(dst, s)
}

public func distance(L: Float64, s!: Array<Float64> = [], t!: Array<Float64> = []): Float64 {
    if (s.size != t.size) {
        throw Exception(badLength)
    }
    if (s.size == 0) {
        return 0.0
    }
    if (L == 2.0) {
        return l2DistanceUnitary(s, t)
    }
    var norm: Float64 = 0.0
    if (L == 1.0) {
        for (i in 0..s.size) {
            norm += abs(t[i] - s[i])
        }
        return norm
    }
    if (L == Float64.Inf) {
        for (i in 0..s.size) {
            let absDiff = abs(t[i]- s[i])
            if (absDiff > norm) {
                norm = absDiff
            }
        }
        return norm
    }
    for (i in 0..s.size) {
        norm += pow(abs(t[i]- s[i]), L)
    }
    return pow(norm, 1.0 / L)
}

public func div(dst: Array<Float64>, s: Array<Float64>): Unit {
    if (dst.size != s.size) {
        throw Exception(badLength)
    }
    f64_div(dst, s)
}

public func divTo(dst: Array<Float64>, s: Array<Float64>, t: Array<Float64>): Array<Float64> {
    if (s.size != t.size) {
        throw Exception(badLength)
    }
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    return f64_divTo(dst, s, t)
}

public func dot(s1: Array<Float64>, s2: Array<Float64>): Float64 {
    if (s1.size != s2.size) {
        throw Exception(badLength)
    }
    return dotUnitary(s1, s2)
}

public func dotUnitary(x: Array<Float64>, y: Array<Float64>): Float64 {
    var sum : Float64 = 0.0
    for (i in 0..x.size){
        sum += y[i] * x[i]
    }
    return sum
}

public func equal(s1: Array<Float64>, s2: Array<Float64>): Bool {
    if (s1.size != s2.size) {
        return false
    }
    for (i in 0..s1.size) {
        if (s2[i] != s1[i]) {
            return false
        }
    }
    return true
}

public func equalApprox(s1: Array<Float64>, s2: Array<Float64>, tol: Float64): Bool {
    if (s1.size != s2.size) {
        return false
    }
    for (i in 0..s1.size) {
        if (!equalWithinAbsOrRel(s1[i], s2[i], tol, tol)) {
            return false
        }
    }
    return true
}

public func equalFunc(s1: Array<Float64>, s2: Array<Float64>, f: (Float64, Float64) -> Bool): Bool {
    if (s1.size != s2.size) {
        return false
    }
    for (i in 0..s1.size) {
        if (!f(s1[i], s2[i])) {
            return false
        }
    }
    return true
}

public func equalLengths(slices: Array<Array<Float64>>): Bool {
    if (slices.size == 0) {
        return true
    }
    let l = slices[0].size
    for (i in 1..slices.size) {
        if (slices[i].size != l) {
            return false
        }
    }
    return true
}

public func find(f: (Float64) -> Bool, s: Array<Float64>, k: Int64, inds!: ArrayList<Int64> = ArrayList<Int64>()): ArrayList<Int64> {
    inds.clear()

    if (k == 0) {
        return inds
    }

    if (k < 0) {
        for (i in 0..s.size) {
            if (f(s[i])) {
                inds.append(i)
            }
        }
        return inds
    }

    var nFound = 0
    for (i in 0..s.size) {
        if (f(s[i])) {
            inds.append(i)
            nFound++
            if (nFound == k) {
                return inds
            }
        }
    }

    throw Exception("floats: insufficient elements found")
}

public func hasNaN(s: Array<Float64>): Bool {
    for (v in s) {
        if (v.isNaN()) {
            return true
        }
    }
    return false
}

public func logSpan(l: Float64, u: Float64, dst!: Array<Float64> = []): Array<Float64> {
    span(dst, log(l), log(u))
    for (i in 0..dst.size) {
        dst[i] = exp(dst[i])
    }
    return dst
}

public func logSumExp(s: Array<Float64>): Float64 {
    let fMaxVal = fMax(s)
    if (fMaxVal.isInf()) {
        return fMaxVal
    }
    var lse: Float64 = 0.0
    for (val in s) {
        lse += exp(val - fMaxVal)
    }
    return log(lse) + fMaxVal
}

public func fMax(s: Array<Float64>): Float64 {
    return s[fMaxIdx(s)]
}

public func fMaxIdx(s: Array<Float64>): Int64 {
    if (s.isEmpty()) {
        throw Exception(zeroLength)
    }
    var max_ = Float64.NaN
    var ind: Int64 = 0
    for (i in 0..s.size) {
        if (s[i].isNaN()) {
            continue
        }
        if (s[i] > max_ || max_.isNaN()) {
            max_ = s[i]
            ind = i
        }
    }
    return ind
}

public func fMin(s: Array<Float64>): Float64 {
    return s[fMinIdx(s)]
}

public func fMinIdx(s: Array<Float64>): Int64 {
    if (s.isEmpty()) {
        throw Exception(zeroLength)
    }
    var min_ = Float64.NaN
    var ind: Int64 = 0
    for (i in 0..s.size) {
        if (s[i].isNaN()) {
            continue
        }
        if (s[i] < min_ || min_.isNaN()) {
            min_ = s[i]
            ind = i
        }
    }
    return ind
}

public func mul(dst: Array<Float64>, s: Array<Float64>): Unit {
    if (dst.size != s.size) {
        throw Exception(badLength)
    }
    for (i in 0..s.size) {
        dst[i] *= s[i]
    }
}

public func mulTo(dst: Array<Float64>, s: Array<Float64>, t: Array<Float64>): Array<Float64> {
    if (s.size  != t.size) {
        throw Exception(badLength)
    }
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    for (i in 0..t.size) {
        dst[i] = t[i] * s[i]
    }
    return dst
}

public func nearestIdx(s: Array<Float64>, v:Float64): Int64 {
    if (s.isEmpty()) {
        throw Exception(zeroLength)
    }
    if (v.isNaN()) {
        return 0
    } else if (v == Float64.Inf) {
        return fMaxIdx(s)
    } else if (v == -Float64.Inf) {
        return fMinIdx(s)
    }
    var ind: Int64 = 0
    var dist = Float64.NaN
    for (i in 0..s.size) {
        let newDist = abs(v - s[i])
        if (newDist.isNaN()) {
            continue
        }
        if (newDist < dist || dist.isNaN()) {
            dist = newDist
            ind = i
        }
    }
    return ind
}

public func nearestIdxForSpan(n: Int64, l: Float64, u: Float64, v: Float64): Int64 {
    if (n < 2) {
        throw Exception(shortSpan)
    }
    if (v.isNaN()) {
        return 0
    }

    if (l.isNaN() && !u.isNaN()) {
        return n - 1
    } else if (u.isNaN()) {
        return 0
    } else if (l.isInf() && u.isInf()) {
        if (l == u) {
            return 0
        }
        if (n % 2 == 1) {
            if (!v.isInf()) {
                return n / 2
            }
            if (copysign(1.0, v) == copysign(1.0, l)) {
                return 0
            }
            return n / 2 + 1
        }
        if (copysign(1.0, v) == copysign(1.0, l)) {
            return 0
        }
        return n / 2
    } else if (l.isInf()) {
        if (v == l) {
            return 0
        }
        return n- 1
    } else if (u.isInf()) {
        if (v == u) {
            return n - 1
        }
        return 0
    } else if (v == -Float64.Inf) {
        if (l <= u) {
            return 0
        }
        return n - 1
    } else if (v == Float64.Inf) {
        if (u <= l) {
            return 0
        }
        return n - 1
    }

    if (l < u) {
        if (v <= l) {
            return 0
        }
        if (v >= u) {
            return n - 1
        }
    } else if (l > u) {
        if (v >= l) {
            return 0
        }
        if (v <= u) {
            return n - 1
        }
    } else {
        return 0
    }

    return Int64((Float64(n) - 1.0) / (u - l) * (v - l) + 0.5)
} 

func copysign(f: Float64, sign: Float64): Float64 {
    const signBit: UInt64 = 1 << 63
    return Float64.fromBits((f.toBits() & (!signBit)) | sign.toBits() & signBit)
}

public func norm(s: Array<Float64>, L: Float64): Float64 {
    if (s.size == 0) {
        return 0.0
    }
    if (L == 2.0) {
        return l2NormUnitary(s)
    }
    var norm_: Float64 = 0.0
    if (L == 1.0) {
        for (val in s) {
            norm_ += abs(val)
        }
        return norm_
    }
    if (L == Float64.Inf) {
        for (val in s) {
            norm_ = fMax(norm_, abs(val))
        }
        return norm_
    }
    for (val in s) {
        norm_ += pow(abs(val), L)
    }
    return pow(norm_, 1.0 / L)
}

public func prod(s: Array<Float64>): Float64 {
    var prod_ = 1.0
    for (val in s) {
        prod_ *= val
    }
    return prod_
}

public func reverse(s: Array<Float64>): Unit {
    s.reverse()
}

public func fSame(s: Array<Float64>, t: Array<Float64>): Bool {
    if (s.size != t.size) {
        return false
    }
    for (i in 0..s.size) {
        let w = t[i]
        if (s[i] != w && !(s[i].isNaN() && w.isNaN())) {
            return false
        }
    }
    return true
}

public func scale(c: Float64, dst: Array<Float64>): Unit {
    if (dst.size > 0) {
        scalUnitary(c, dst)
    }
}

public func scaleTo(dst: Array<Float64>, c: Float64, s: Array<Float64>): Array<Float64> {
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    if (dst.size > 0) {
        scalUnitaryTo(dst, c, s)
    }
    return dst
}

public func span(dst: Array<Float64>, l: Float64, u: Float64): Array<Float64> {
    let n = dst.size
    if (n < 2) {
        throw Exception(shortSpan)
    }

    if (l.isNaN()) {
        for (i in 0..dst[0..dst.size - 1].size) {
            dst[i] = Float64.NaN
        }
        dst[dst.size - 1] = u
        return dst
    }
    else if (u.isNaN()) {
        for (i in 0..dst[1..].size) {
            dst[i + 1] = Float64.NaN
        }
        dst[0] = l
        return dst 
    }
    else if (l.isInf() && u.isInf()) {
        for (i in 0..dst[0..dst.size / 2].size) {
            dst[i] = l
            dst[dst.size - i - 1] = u
        }
        if (dst.size % 2 == 1) {
            if (l != u) {
                dst[dst.size / 2] = 0.0
            } else {
                dst[dst.size / 2] = l
            }
        }
        return dst
    }
    else if (l.isInf()) {
        for (i in 0..dst[0..dst.size - 1].size) {
            dst[i] = l
        }
        dst[dst.size - 1] = u
        return dst
    } else if (u.isInf()) {
        for (i in 0..dst[1..].size) {
            dst[i + 1] = u
        }
        dst[0] = l
        return dst
    }

    let step = (u - l) / Float64(n - 1)
    for (i in 0..dst.size) {
        dst[i] = l + step * Float64(i)
    }
    return dst
}

public func sub(dst: Array<Float64>, s: Array<Float64>): Unit {
    if (dst.size != s.size) {
        throw Exception(badLength)
    }
    axpyUnitaryTo(dst, -1.0, s, dst)
}

public func subTo(dst: Array<Float64>, s: Array<Float64>, t: Array<Float64>): Array<Float64> {
    if (s.size != t.size) {
        throw Exception(badLength)
    }
    if (dst.size != s.size) {
        throw Exception(badDstLength)
    }
    axpyUnitaryTo(dst, -1.0, t, s)
    return dst
}

public func sum(s: Array<Float64>): Float64 {
    return f64_sum(s)
}

public func within(s: Array<Float64>, v: Float64): Int64 {
    if (s.size < 2) {
        throw Exception(shortSpan)
    }
    if (!isSort(s)) {
        throw Exception("floats: input slice not sorted")
    }
    if (v < s[0] || v >= s[s.size - 1] || v.isNaN()) {
        return -1
    }
    for (i in 0..s[1..].size) {
        if (v < s[i + 1]) {
            return i
        }
    }
    return -1
}

func isSort<T>(arr: Array<T>): Bool where T <: Comparable<T> {
    if (arr.size <= 1) {
        return true
    }
    for (i in 1..arr.size) {
        if (arr[i - 1] > arr[i]) {
            return false
        }
    }
    return true
}

public func sumCompensated(s: Array<Float64>): Float64 {
    var sum: Float64 = 0.0
    var c: Float64 = 0.0
    for (x in s) {
        let t = Float64(sum + x)
        if (abs(sum) >= abs(x)) {
            c += (sum - t) + x
        } else {
            c += (x - t) + sum
        }
        sum = t
    }
    return sum + c
}

