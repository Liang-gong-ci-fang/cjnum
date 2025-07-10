package cjnum.mat

import std.core.*

public func mPanics(fn: () -> Unit): (Bool, String) {
    try {
        fn()
    } catch (e: Exception) {
        return (true, e.message)
    }
    return (false, "")
}
