package cjnum.internal.asm.f32

import std.unittest.*
import std.unittest.testmacro.*
import cjnum.floats.scalar.same as scalar_same
import cjnum.floats.scalar.equalWithinAbsOrRel
import std.math.*

var nan: Float32 = Float32.NaN
var inf: Float32 = Float32.Inf

func sameApprox(x: Float32, y: Float32, tol: Float32): Bool {
    let a: Float64 = Float64(x)
    let b: Float64 = Float64(y)
    return same(x, y) || equalWithinAbsOrRel(a, b, Float64(tol), Float64(tol))
}

func same(x: Float32, y: Float32): Bool {
    return scalar_same(Float64(x), Float64(y))
}

func sameStrided(ref: Array<Float32>, x: Array<Float32>, inc: Int64): Bool {
    let inc_ = abs(inc)
    for (i in 0..ref.size) {
        if (!same(x[i * inc_], ref[i])) {
            return false
        }
    }
    return true
}

func copy<T>(dst: Array<T>, src: Array<T>): Int64 {
    let minSize = min(src.size, dst.size)
    dst[0..minSize] = src[0..minSize]
    return minSize
}

func guardVector(v: Array<Float32>, g: Float32, gdLn: Int64): Array<Float32> {
    let guarded = Array<Float32>(v.size + gdLn * 2, item: 0.0)
    copy(guarded[gdLn..], v)
    for (i in 0..gdLn) {
        guarded[i] = g
        guarded[guarded.size - 1 - i] = g
    }
    return guarded
}

func isValidGuard(v: Array<Float32>, g: Float32, gdLn: Int64): Bool {
    for (i in 0..gdLn) {
        if (!same(v[i], g) || !same(v[v.size - 1 - i], g)) {
            return false
        }
    }
    return true
}

func guardIncVector(vec: Array<Float32>, gdVal: Float32, inc: Int64, gdLen: Int64): Array<Float32> {
    let inc_ = abs(inc)
    let inrLen = vec.size * inc_
    let guarded = Array<Float32>(inrLen + gdLen * 2, item: 0.0)
    for (i in 0..guarded.size) {
        guarded[i] = gdVal
    }
    for (i in 0..vec.size) {
        guarded[gdLen + i * inc_] = vec[i]
    }
    return guarded
}

func checkValidIncGuard(v: Array<Float32>, g: Float32, inc: Int64, gdLn: Int64): Bool {
    let srcLn = v.size - 2 * gdLn
    for (i in 0..v.size) {
        if (same(v[i], g)) {
        // correct value
        } else if (i < gdLn) {
            return false
        } else if (i > gdLn + srcLn) {
            return false
        } else if ((i - gdLn) % inc == 0 && (i - gdLn) / inc < v.size) {
        // correct value
        } else {
            return false
        }
    }
    return true
}

let align2 = newIncSet(0, 1, 2, 3)
let align3 = newIncToSet(0, 1, 2, 3)

class IncSet {
    IncSet(
        var x: Int64,
        var y: Int64
    ) {}
}

func newIncSet(inc: Array<Int64>): Array<IncSet> {
    let n = inc.size
    var is_ = Array<IncSet>(n * n, item: IncSet(0, 0))
    for (x in 0..inc.size) {
        for (y in 0..inc.size) {
            is_[x * n + y].x = inc[x]
            is_[x * n + y].y = inc[y]
        }
    }
    return is_
}

class IncToSet {
    IncToSet(
        var dst: Int64,
        var x: Int64,
        var y: Int64
    ) {}
}

func newIncToSet(inc: Array<Int64>): Array<IncToSet> {
    let n = inc.size
    var is_ = Array<IncToSet>(n * n * n, item: IncToSet(0, 0, 0))
    for (i in 0..inc.size) {
        for (x in 0..inc.size) {
            for (y in 0..inc.size) {
                is_[i * n * n + x * n + y].dst = inc[i]
                is_[i * n * n + x * n + y].x = inc[x]
                is_[i * n * n + x * n + y].y = inc[y]
            }
        }
    }
    return is_
}

let _ = {=> align2}()

let _ = {=> align3}()
