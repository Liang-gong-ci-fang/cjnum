package cjnum.blas.internal.f64

import std.math.*

const msgVal      = "%v: unexpected value at %v Got: %v Expected: %v"
const msgGuard    = "%v: Guard violated in %s vector %v %v"
const msgReadOnly = "%v: modified read-only %v argument"

let nan = Float64.NaN
let inf = Float64.Inf

func newGuardedVector(data: Array<Float64>, inc: Int64): (Array<Float64>, Array<Float64>, Array<Float64>) {
    let iinc = abs(inc)
    let guard = 2 * iinc
    let size = (data.size - 1) * iinc + 1
    let whole = Array<Float64>(size + 2 * guard, item: 0.0)
    let v = whole[guard..(whole.size - guard)]
    for (i in 0..whole.size) {
        whole[i] = Float64.NaN
    }
    for (i in 0..data.size) {
        v[i * iinc] = data[i]
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

// TODO