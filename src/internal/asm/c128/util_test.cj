package cjnum.internal.asm.c128

import cjnum.cmplxs.cscalar as cscalar
import cjnum.complex.*
import cjnum.floats.scalar as scalar

var inf = Float64.Inf
var benchSink: Complex128 = Complex128(0.0, 0.0)

func copy<T>(dst: Array<T>, src: Array<T>): Int64 {
    let minSize = min(src.size, dst.size)
    dst[0..minSize] = src[0..minSize]
    return minSize
}

func guardVector(vec: Array<Complex128>, guard_val: Complex128, guard_len: Int64): Array<Complex128> {
    var guarded = Array<Complex128>(vec.size + guard_len * 2, repeat: Complex128(0.0))
    copy(guarded[guard_len..], vec)
    for (i in 0..guard_len) {
        guarded[i] = guard_val
        guarded[guarded.size - 1 - i] = guard_val
    }
    return guarded
}

func isValidGuard(vec: Array<Complex128>, guard_val: Complex128, guard_len: Int64): Bool {
    for (i in 0..guard_len) {
        if (vec[i] != guard_val || vec[vec.size - 1 - i] != guard_val) {
            return false
        }
    }
    return true
}

func guardIncVector(vec: Array<Complex128>, guard_val: Complex128, inc: Int64, guard_len: Int64): Array<Complex128> {
    var s_ln = vec.size * inc
    if (inc < 0) {
        s_ln = vec.size * -inc
    }
    var guarded = Array<Complex128>(s_ln + guard_len * 2, repeat: Complex128(0.0))
    var cas = 0
    for (i in 0..guarded.size) {
        match {
            case i < guard_len || i > guard_len + s_ln => guarded[i] = guard_val
            case (i - guard_len) % (inc) == 0 && cas < vec.size =>
                guarded[i] = vec[cas]
                cas++
            case _ => guarded[i] = guard_val
        }
    }
    return guarded
}

func checkValidIncGuard(vec: Array<Complex128>, guard_val: Complex128, inc: Int64, guard_len: Int64): Bool {
    var s_ln = vec.size - 2 * guard_len
    if (inc < 0) {
        s_ln = vec.size * -inc
    }
    for (i in 0..vec.size) {
        if (vec[i] == guard_val) {
        // correct value
        } else if (i < guard_len) {
            return false
        } else if (i > guard_len + s_ln) {
            return false
        } else if ((i - guard_len) % inc == 0 && (i - guard_len) / inc < vec.size) {
        // correct value
        } else {
            return false
        }
    }
    return true
}

var align1: Array<Int64> = [0, 1]
var align2 = newIncSet(0, 1)
var align3 = newIncToSet(0, 1)

class IncSet {
    IncSet(
        var x: Int64,
        var y: Int64
    ) {}
}

func newIncSet(inc: Array<Int64>): Array<IncSet> {
    let n = inc.size
    var is_ = Array<IncSet>(n * n, repeat: IncSet(0, 0))
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
    var is_ = Array<IncToSet>(n * n * n, repeat: IncToSet(0, 0, 0))
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

func sameCmplxApprox(a: Complex128, b: Complex128, tol: Float64): Bool {
    return cscalar.same(a, b) || cscalar.equalWithinAbsOrRel(a, b, tol, tol)
}

func sameApprox(a: Float64, b: Float64, tol: Float64): Bool {
    return scalar.same(a, b) || scalar.equalWithinAbsOrRel(a, b, tol, tol)
}

let _ = {=> inf}()

let _ = {=> benchSink}()

let _ = {=> align1}()

let _ = {=> align2}()

let _ = {=> align3}()
