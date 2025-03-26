package cjnum.floats

import std.math.*

// TODO

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

public func same(s: Array<Float64>, t: Array<Float64>): Bool {
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

// TODO
