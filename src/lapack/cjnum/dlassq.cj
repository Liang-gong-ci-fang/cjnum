package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*
import cjnum.blas.*

extend LFloat64Implementation {
    public func dlassq(n: Int64, x: Array<Float64>, incX: Int64, scale: Float64, sumsq: Float64): (Float64, Float64) {
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (incX <= 0) {
            throw Exception(badIncX)
        }
        if (x.size < 1 + (n - 1) * incX) {
            throw Exception(shortX)
        }

        if (scale.isNaN() || sumsq.isNaN()) {
            return (scale, sumsq)
        }

        var scale_ = scale
        var sumsq_ = sumsq
        if (sumsq == 0.0) {
            scale_ = 1.0
        }
        if (scale_ == 0.0) {
            scale_ = 1.0
            sumsq_ = 0.0
        }

        if (n == 0) {
            return (scale_, sumsq_)
        }

        var isBig = false
        var asml = 0.0
        var amed = 0.0
        var abig = 0.0
        var ix = 0
        for (_ in 0..n) {
            var ax = abs(x[ix])
            if (ax > dtbig) {
                ax *= dsbig
                abig += ax * ax
                isBig = true
            } else if (ax < dtsml) {
                if (!isBig) {
                    ax *= dssml
                    asml += ax * ax
                }
            } else {
                amed += ax * ax
            }
            ix += incX
        }

        if (sumsq_ > 0.0) {
            let ax = scale_  * sqrt(sumsq_)
            if (ax > dtbig) {
                if (scale_ > 1.0) {
                    scale_ *= dsbig
                    abig += scale_ * (scale_ * sumsq_)
                } else {
                    abig += scale_ * (scale_ * (dsbig * (dsbig * sumsq_)))
                }
            } else if (ax < dtsml) {
                if (!isBig) {
                    if (scale_ < 1.0) {
                        scale_ *= dssml
                        asml += scale_ * (scale_ * sumsq_)
                    } else {
                        asml += scale_ * (scale_ * (dssml * (dssml * sumsq_)))
                    }
                }
            } else {
                amed += scale_ * (scale_ * sumsq_)
            }
        }

        if (abig > 0.0) {
            if (amed > 0.0 || amed.isNaN()) {
                abig += (amed * dsbig) * dsbig
            }
            scale_ = 1.0 / dsbig
            sumsq_ = abig
        } else if (asml > 0.0) {
            if (amed > 0.0 || amed.isNaN()) {
                amed = sqrt(amed)
                asml = sqrt(asml) / dssml
                let (ymin, ymax) = if (asml > amed) {
                    (amed, asml)
                } else {
                    (asml, amed)
                }
                scale_ = 1.0
                sumsq_ = ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax))
            } else {
                scale_ = 1.0 / dssml
                sumsq_ = asml
            }
        } else {
            scale_ = 1.0
            sumsq_ = amed
        }

        return (scale_, sumsq_)
    }
}