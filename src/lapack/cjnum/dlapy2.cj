package cjnum.lapack.cjnum

import std.math.*

extend LFloat64Implementation {
    public func dlapy2(p: Float64, q: Float64): Float64 {
        var p_ = abs(p)
        var q_ = abs(q)

        if (p_ == Float64.Inf || q_ == Float64.Inf) {
            return Float64.Inf
        }

        if (p_.isNaN() || q_.isNaN()) {
            return Float64.NaN
        }

        if (p_ < q_) {
            (p_, q_) = (q_, p_)
        }

        if (p_ == 0.0) {
            return 0.0
        }

        q_ = q_ / p_
        
        return p_ * sqrt(1.0 + pow(q_, 2))
    }
}