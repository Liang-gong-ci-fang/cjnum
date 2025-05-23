package cjnum.internal.asm.f64

import std.math.*
import std.random.*
import cjnum.floats.scalar.*

const msgVal = "%v: unexpected value at %v Got: %v Expected: %v"
const msgGuard = "%v: Guard violated in %s vector %v %v"
const msgReadOnly = "%v: modified read-only %v argument"
let nan = Float64.NaN
let inf = Float64.Inf

func newGuardedVector(data: Array<Float64>, inc: Int64): (Array<Float64>, Array<Float64>, Array<Float64>) {
    let inc_ = abs(inc)
    let guard = 2 * inc_
    let size = (data.size - 1) * inc_ + 1
    let whole = Array<Float64>(size + 2 * guard, item: 0.0)
    let v = whole[guard..(whole.size - guard)]
    for (i in 0..whole.size) {
        whole[i] = Float64.NaN
    }
    for (i in 0..data.size) {
        v[i * inc_] = data[i]
    }
    return (v, whole[..guard], whole[(whole.size - guard)..])
}

func allNaN(x: Array<Float64>): Bool {
    for (v in x) {
        if (!v.isNaN()) {
            return false
        }
    }
    return true
}

func equalStrided(ref: Array<Float64>, x: Array<Float64>, inc: Int64): Bool {
    let inc_ = abs(inc)
    for (i in 0..ref.size) {
        if (!same(x[i * inc_], ref[i])) {
            return false
        }
    }
    return true
}

func nonStridedWrite(x: Array<Float64>, inc: Int64): Bool {
    let inc_ = abs(inc)
    for (i in 0..x.size) {
        if (i % inc_ != 0 && !x[i].isNaN()) {
            return true
        }
    }
    return false
}

func guardVector(vec: Array<Float64>, gdVal: Float64, gdLn: Int64): Array<Float64> {
    let guarded = Array<Float64>(
        vec.size + gdLn * 2,
        {
            i => if (i >= gdLn && i < vec.size + gdLn) {
                vec[i - gdLn]
            } else {
                0.0
            }
        }
    )
    for (i in 0..gdLn) {
        guarded[i] = gdVal
        guarded[guarded.size - 1 - i] = gdVal
    }
    return guarded
}

func isValidGuard(vec: Array<Float64>, gdVal: Float64, gdLn: Int64): Bool {
    for (i in 0..gdLn) {
        if (!same(vec[i], gdVal) || !same(vec[vec.size - 1 - i], gdVal)) {
            return false
        }
    }
    return true
}

func guardIncVector(vec: Array<Float64>, gdVal: Float64, inc: Int64, gdLen: Int64): Array<Float64> {
    let inc_ = abs(inc)
    let inrLen = vec.size * inc_
    let guarded = Array<Float64>(inrLen + gdLen * 2, item: gdVal)
    for (i in 0..vec.size) {
        guarded[gdLen + i * inc_] = vec[i]
    }
    return guarded
}

func checkValidIncGuard(vec: Array<Float64>, gdVal: Float64, inc: Int64, gdLen: Int64): Bool {
    let srcLn = vec.size - 2 * gdLen
    for (i in 0..vec.size) {
        if (same(vec[i], gdVal)) {
            return true
        } else if ((i - gdLen) % inc == 0 && (i - gdLen) / inc < vec.size) {
        // ignore
        } else if (i < gdLen) {
            return false
        } else if (i > gdLen + srcLn) {
            return false
        } else {
            return false
        }
    }
    return true
}

func sameApprox(a: Float64, b: Float64, tol: Float64): Bool {
    return same(a, b) || equalWithinAbsOrRel(a, b, tol, tol)
}

let align1 = [0, 1]
let align2 = newIncSet([0, 1])
let align3 = newIncToSet([0, 1])

class IncSet {
    IncSet(
        let x: Int64,
        var y: Int64
    ) {}
}

func newIncSet(inc: Array<Int64>): Array<IncSet> {
    let n = inc.size
    let is_ = Array<IncSet>(n * n, {_ => IncSet(0, 0)})
    for (x in 0..n) {
        for (y in 0..n) {
            is_[x * n + y] = IncSet(inc[x], inc[y])
        }
    }
    return is_
}

class IncToSet {
    IncToSet(
        var dst: Int64,
        let x: Int64,
        let y: Int64
    ) {}
}

func newIncToSet(inc: Array<Int64>): Array<IncToSet> {
    let n = inc.size
    let is_ = Array<IncToSet>(n * n * n, {_ => IncToSet(0, 0, 0)})
    for (i in 0..n) {
        for (x in 0..n) {
            for (y in 0..n) {
                is_[i * n * n + x * n + y] = IncToSet(inc[i], inc[x], inc[y])
            }
        }
    }
    return is_
}

func randomSlice(n: Int64, inc: Int64): Array<Float64> {
    let inc_ = abs(inc)
    let random = Random()
    return Array<Float64>((n - 1) * inc_ + 1, {_ => random.nextFloat64()})
}

func randSlice(n: Int64, inc: Int64, r: Random): Array<Float64> {
    let inc_ = abs(inc)
    return Array<Float64>((n - 1) * inc_ + 1, {_ => r.nextFloat64()})
}

let _ = {=> msgVal}()

let _ = {=> msgGuard}()

let _ = {=> msgReadOnly}()
