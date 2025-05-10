package cjnum.lapack.cjnum

import std.math.*
import cjnum.lapack.*

func ilaenv(ispec: Int64, name: String, opts: String, n1: Int64, n2: Int64, n3: Int64, n4: Int64): Int64 {
    let sname = name[0] == b'S' || name[0] == b'D'
    let cname = name[0] == b'C' || name[0] == b'Z'
    if (!sname && !cname) {
        throw Exception(badName)
    }

    let c2 = name[1..3]
    let c3 = name[3..6]
    let c4 = c3[1..3]

    match (ispec) {
        case 1 =>
            match (c2) {
                case "GE" =>
                    match (c3) {
                        case "TRF" =>
                            if (sname) {
                                return 64
                            }
                            return 64
                        case "QRF" | "RQF"| "LQF" | "QLF" =>
                            if (sname) {
                                return 32
                            }
                            return 32
                        case "HRD" =>
                            if (sname) {
                                return 32
                            }
                            return 32
                        case "BRD" =>
                            if (sname) {
                                return 32
                            }
                            return 32
                        case "TRI" => 
                            if (sname) {
                                return 64
                            }
                            return 64
                        case _ => throw Exception(badName)
                    }
                case "PO" => 
                    match (c3) {
                        case "TRF" =>
                            if (sname) {
                                return 64
                            }
                            return 64
                        case _ => throw Exception(badName)
                    }
                case "SY" => 
                    match (c3) {
                        case "TRF" => 
                            if (sname) {
                                return 64
                            }
                            return 64
                        case "TRD" =>
                            return 32
                        case "GST" =>
                            return 64
                        case _ => throw Exception(badName)
                    }
                case "HE" =>
                    match (c3) {
                        case "TRF" => return 64
                        case "TRD" => return 32
                        case "GST" => return 64
                        case _ => throw Exception(badName)
                    }
                case "OR" =>
                    match (c3[0]) {
                        case b'G' =>
                            match (c3[1..]) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 32
                                case _ => throw Exception(badName)
                            }
                        case b'M' =>
                            match (c3[1..]) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 32
                                case _ => throw Exception(badName)
                            }
                        case _ => throw Exception(badName)
                    }
                case "UN" =>
                    match (c3[0]) {
                        case b'G' =>
                            match (c3[1..]) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 32
                                case _ => throw Exception(badName)
                            }
                        case b'M' =>
                            match (c3[1..]) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 32
                                case _ => throw Exception(badName)
                            }
                        case _ => throw Exception(badName)
                    }
                case "GB" => 
                    match (c3) {
                        case "TRF" => 
                            if (sname) {
                                if (n4 <= 64) {
                                    return 1
                                }
                                return 32
                            }
                            if (n4 <= 64) {
                                return 1
                            }
                            return 32
                        case _ => throw Exception(badName)
                    }
                case "PB" =>
                    match (c3) {
                        case "TRF" => 
                            if (sname) {
                                if (n2 <= 64) {
                                    return 1
                                }
                                return 32
                            }
                            if (n2 <= 64) {
                                return 1
                            }
                            return 32
                        case _ => throw Exception(badName)
                    }
                case "PT" =>
                    match (c3) {
                        case "TRS" => return 1
                        case _ => throw Exception(badName)
                    }
                case "TR" =>
                    match (c3) {
                        case "TRI" =>
                            if (sname) {
                                return 64
                            }
                            return 64
                        case "EVC" =>
                            if (sname) {
                                return 64
                            }
                        case _ => throw Exception(badName)
                    }
                case "LA" =>
                    match (c3) {
                        case "UUM" =>
                            if (sname) {
                                return 64
                            }
                        case _ => throw Exception(badName)
                    }
                case "ST" =>
                    if (sname && c3 == "EBZ") {
                        return 1
                    }
                    throw Exception(badName)
                case _ => throw Exception(badName)
            }
        case 2 =>
            match (c2) {
                case "GE" =>
                    match (c3) {
                        case "QRF" | "RQF" | "LQF" | "QLF" =>
                            if (sname) {
                                return 2
                            }
                            return 2
                        case "HRD" =>
                            if (sname) {
                                return 2
                            }
                            return 2
                        case "BRD" =>
                            if (sname) {
                                return 2
                            }
                            return 2
                        case "TRI" =>
                            if (sname) {
                                return 2
                            }
                            return 2
                        case _ => throw Exception(badName)
                    }
                case "SY" =>
                    match (c3) {
                        case "TRF" =>
                            if (sname) {
                                return 8
                            }
                            return 8
                        case "TRD" =>
                            if (sname) {
                                return 2
                            }
                            throw Exception(badName)
                        case _ => throw Exception(badName)
                    }
                case "HE" =>
                    if (c3 == "TRD") {
                        return 2
                    }
                    throw Exception(badName)
                case "OR" =>
                    if (!sname) {
                        throw Exception(badName)
                    }
                    match (c3[0]) {
                        case b'G' =>
                            match (c4) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 2
                                case _ => throw Exception(badName)
                            }
                        case b'M' =>
                            match (c4) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 2
                                case _ => throw Exception(badName)
                            }
                        case _ => throw Exception(badName)
                    }
                case "UN" =>
                    if (!sname) {
                        throw Exception(badName)
                    }
                    match (c3[0]) {
                        case b'G' =>
                            match (c4) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 2
                                case _ => throw Exception(badName)
                            }
                        case b'M' =>
                            match (c4) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 2
                                case _ => throw Exception(badName)
                            }
                        case _ => throw Exception(badName)
                    }
                case _ => throw Exception(badName)
            }
        case 3 =>
            match (c2) {
                case "GE" =>
                    match (c3) {
                        case "QRF" | "RQF" | "LQF" | "QLF" =>
                            if (sname) {
                                return 128
                            }
                            return 128
                        case "HRD" =>
                            if (sname) {
                                return 128
                            }
                            return 128
                        case "BRD" =>
                            if (sname) {
                                return 128
                            }
                            return 128
                        case _ => throw Exception(badName)
                    }
                case "SY" =>
                    if (sname && c3 == "TRD") {
                        return 32
                    }
                    throw Exception(badName)
                case "HE" =>
                    if (c3 == "TRD") {
                        return 32
                    }
                    throw Exception(badName)
                case "OR" =>
                    if (!sname) {
                        throw Exception(badName)
                    }
                    match (c3[0]) {
                        case b'G' =>
                            match (c4) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 128
                                case _ => throw Exception(badName)
                            }
                        case _ => throw Exception(badName)
                    }
                case "UN" =>
                    if (!sname) {
                        throw Exception(badName)
                    }
                    match (c3[0]) {
                        case b'G' =>
                            match (c4) {
                                case "QR" | "RQ" | "LQ" | "QL" | "HR" | "TR" | "BR" => return 128
                                case _ => throw Exception(badName)
                            }
                        case _ => throw Exception(badName)
                    }
                case _ => throw Exception(badName)
            }
        case 4 =>
            return 6
        case 5 => 
            return 2
        case 6 =>
            return Int64(Float64(min(n1, n2)) * 1.6)
        case 7 =>
            return 1
        case 8 =>
            return 50
        case 9 =>
            return 25
        case 10 =>
            return 1
        case 11 =>
            return 1
        case 12 | 13 | 14 | 15 | 16 =>
            iparmq(ispec, name, opts, n1, n2, n3, n4)
        case _ => throw Exception(badIspec)
    }
    throw Exception("Mismatched exception")
}