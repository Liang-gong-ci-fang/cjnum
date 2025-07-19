package cjnum.lapack.cjnum

import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlasq3(i0: Int64, n0: Int64, z: Array<Float64>, pp: Int64, dmin: Float64, sigma: Float64,
        desig: Float64, qmax: Float64, nFail: Int64, iter: Int64, nDiv: Int64, ttype: Int64, dmin1: Float64,
        dmin2: Float64, dn: Float64, dn1: Float64, dn2: Float64, g: Float64, tau: Float64): (Int64, Int64, Int64, 
        Float64, Float64, Float64, Float64, Int64, Int64, Int64, Int64, Float64, Float64, Float64, Float64, Float64, 
        Float64, Float64) {
        var i0_ = i0
        var n0_ = n0
        var pp_ = pp
        var dmin_ = dmin
        var sigma_ = sigma
        var desig_ = desig
        var qmax_ = qmax
        var nFail_ = nFail
        var iter_ = iter
        var nDiv_ = nDiv
        var ttype_ = ttype
        var dmin1_ = dmin1
        var dmin2_ = dmin2
        var dn_ = dn
        var dn1_ = dn1
        var dn2_ = dn2
        var g_ = g
        var tau_ = tau

        if (i0_ < 0) {
            throw Exception(i0LT0)
        } else if (n0_ < 0) {
            throw Exception(n0LT0)
        } else if (z.size < 4 * n0_) {
            throw Exception(shortZ)
        } else if (pp_ != 0 && pp_ != 1 && pp_ != 2) {
            throw Exception(badPp)
        }

        let cbias = 1.5

        let n0in = n0_
        let eps = dlamchP
        let tol = eps * 100.0
        let tol2 = tol * tol
        var nn: Int64
        var t: Float64
        while (true) {
            if (n0_ < i0_) {
                return (i0_, n0_, pp_, dmin_, sigma_, desig_, qmax_, nFail_, iter_, nDiv_, ttype_, dmin1_, dmin2_, dn_,
                    dn1_, dn2_, g_, tau_)
            }
            if (n0_ == i0_) {
                z[4 * (n0_ + 1) - 4] = z[4 * (n0_ + 1) + pp_ - 4] + sigma_
                n0_--
                continue
            }
            nn = 4 * (n0_ + 1) + pp_ - 1
            if (n0_ != i0_ + 1) {
                if (z[nn - 5] > tol2 * (sigma_ + z[nn - 3]) && z[nn - 2 * pp_ - 4] > tol2 * z[nn - 7]) {
                    if (z[nn - 9] > tol2 * sigma_ && z[nn - 2 * pp_ - 8] > tol2 * z[nn - 11]) {
                        break
                    }
                } else {
                    z[4 * (n0_ + 1) - 4] = z[4 * (n0_ + 1) + pp_ - 4] + sigma_
                    n0_--
                    continue
                }
            }
            if (z[nn - 3] > z[nn - 7]) {
                (z[nn - 3], z[nn - 7]) = (z[nn - 7], z[nn - 3])
            }
            t = 0.5 * (z[nn - 7] - z[nn - 3] + z[nn - 5])
            if (z[nn - 5] > z[nn - 3] * tol2 && t != 0.0) {
                var s = z[nn - 3] * (z[nn - 5] / t)
                if (s <= t) {
                    s = z[nn - 3] * (z[nn - 5] / (t * (1.0 + sqrt(1.0 + s / t))))
                } else {
                    s = z[nn - 3] * (z[nn - 5] / (t + sqrt(t) * sqrt(t + s)))
                }
                t = z[nn - 7] + (s + z[nn - 5])
                z[nn - 3] *= z[nn - 7] / t
                z[nn - 7] = t
            }
            z[4 * (n0_ + 1) - 8] = z[nn - 7] + sigma_
            z[4 * (n0_ + 1) - 4] = z[nn - 3] + sigma_
            n0_ -= 2
        }
        if (pp_ == 2) {
            pp_ = 0
        }

        if (dmin_ <= 0.0 || n0_ < n0in) {
            if (cbias * z[4 * (i0_ + 1) + pp_ - 4] < z[4 * (n0_ + 1) + pp_ - 4]) {
                let ipn4Out = 4 * (i0_ + n0_ + 2)
                for (j4loop in 4 * (i0_ + 1)..2 * ((i0_ + 1) + (n0_ + 1) - 1) : 4) {
                    let ipn4 = ipn4Out - 1
                    let j4 = j4loop - 1
                    (z[j4 - 3], z[ipn4 - j4 - 4]) = (z[ipn4 - j4 - 4], z[j4 - 3])
                    (z[j4 - 2], z[ipn4 - j4 - 3]) = (z[ipn4 - j4 - 3], z[j4 - 2])
                    (z[j4 - 1], z[ipn4 - j4 - 6]) = (z[ipn4 - j4 - 6], z[j4 - 1])
                    (z[j4], z[ipn4 - j4 - 5]) = (z[ipn4 - j4 - 5], z[j4])
                }
                if (n0_ - i0_ <= 4) {
                    z[4 * (n0_ + 1) + pp_ - 2] = z[4 * (i0_ + 1) + pp_ - 2]
                    z[4 * (n0_ + 1) - pp_ - 1] = z[4 * (i0_ + 1) - pp_ - 1]
                }
                dmin2_ = min(dmin2_, z[4 * (i0_ + 1) - pp_ - 2])
                z[4 * (n0_ + 1) + pp_ - 2] = min(min(z[4 * (n0_ + 1) + pp_ - 2], z[4 * (i0_ + 1) + pp_ - 2]),
                    z[4 * (i0_ + 1) + pp_ + 2])
                z[4 * (n0_ + 1) - pp_ - 1] = min(min(z[4 * (n0_ + 1) - pp_ - 1], z[4 * (i0_ + 1) - pp_ - 1]),
                    z[4 * (i0_ + 1) - pp_ + 3])
                qmax_ = max(max(qmax_, z[4 * (i0_ + 1) + pp_ - 4]), z[4 * (i0_ + 1) + pp_])
                dmin_ = copysign(0.0, -1.0)
            }
        }

        (tau_, ttype_, g_) = dlasq4(i0_, n0_, z, pp_, n0in, dmin_, dmin1_, dmin2_, dn_, dn1_, dn2_, tau_, ttype_, g_)

        var done = true

        while (true) {
            (i0_, n0_, pp_, tau_, sigma_, dmin_, dmin1_, dmin2_, dn_, dn1_, dn2_) = dlasq5(i0_, n0_, z, pp_, tau_,
                sigma_)
            nDiv_ += n0_ - i0_ + 2
            iter_ += 1
            if (dmin_ >= 0.0 && dmin1_ >= 0.0) {
                done = false
                break
            } else if (dmin_ < 0.0 && dmin1_ > 0.0 && z[4 * n0_ - pp_ - 1] < tol * (sigma_ + dn1_) && abs(dn_) < tol *
                sigma_) {
                z[4 * n0_ - pp_ + 1] = 0.0
                dmin_ = 0.0
                done = false
                break
            } else if (dmin_ < 0.0) {
                nFail_++
                if (ttype_ < -22) {
                    tau_ = 0.0
                } else if (dmin1_ > 0.0) {
                    tau_ = (tau_ + dmin_) * (1.0 - 2.0 * eps)
                    ttype_ -= 11
                } else {
                    tau_ = tau_ / 4.0
                    ttype_ -= 12
                }
            } else if (dmin_.isNaN()) {
                if (tau_ == 0.0) {
                    break
                }
                tau_ = 0.0
            } else {
                break
            }
        }

        if (done) {
            (dmin_, dmin1_, dmin2_, dn_, dn1_, dn2_) = dlasq6(i0_, n0_, z, pp_)
            nDiv_ += n0_ - i0_ + 2
            iter_ += 1
            tau_ = 0.0
        }

        if (tau_ < sigma_) {
            desig_ += tau_
            t = sigma_ + desig_
            desig_ -= t - sigma_
        } else {
            t = sigma_ + tau_
            desig_ += sigma_ - (t - tau_)
        }
        sigma_ = t
        return (i0_, n0_, pp_, dmin_, sigma_, desig_, qmax_, nFail_, iter_, nDiv_, ttype_, dmin1_, dmin2_, dn_, dn1_,
            dn2_, g_, tau_)
    }
}
