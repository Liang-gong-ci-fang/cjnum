package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*

func iparmq(ispec: Int64, name: String, opts: String, n: Int64, ilo: Int64, ihi: Int64, lwork: Int64): Int64 {
    let nh = ihi - ilo + 1
    var ns = match {
        case nh >= 30 => 4
        case nh >= 60 => 10
        case nh >= 150 => max(10, nh / Int64(log(Float64(nh) / log(2.0))))
        case nh >= 590 => 64
        case nh >= 3000 => 128
        case nh >= 6000 => 256
        case _ => 2
    }
    ns = max(2, ns - (ns % 2))

    match (ispec) {
        case 12 => 
            const nmin = 75
            return nmin
        case 13 =>
            const knwswp = 500
            if (nh <= knwswp) {
                return ns
            }
            return 3 * ns / 2
        case 14 =>
            const nibble = 14
            return nibble
        case 15 =>
            return ns
        case 16 =>
            if (name.size != 6) {
                throw Exception(badName)
            }
            const k22min = 14
            const kacmin = 14
            var acc22 = 0
            if (name[1..] == "GGHRD" || name[1..] == "GGHD3") {
                acc22 = 1
                if (nh >= k22min) {
                    acc22 = 2
                }
            } else if (name[3..] == "EXC") {
                if (nh >= kacmin) {
                    acc22 = 1
                }
                if (nh >= k22min) {
                    acc22 = 2
                }
            } else if (name[1..] == "HSEQR" || name[1..5] == "LAQR") {
                if (ns >= kacmin) {
                    acc22 = 1
                }
                if (ns >= k22min) {
                    acc22 = 2
                }
            }
            return acc22
        case _ => throw Exception(badIspec)
    }
    throw Exception("Mismatched exception")
}