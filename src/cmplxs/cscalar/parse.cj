package cjnum.cmplxs.cscalar

import cjnum.complex.*
import std.core.*
import std.math.*
import std.convert.*

func cParse(s: String): Complex128 {

    if (s.isEmpty()) { 
        throw ParseExpection(information:"quat: failed to parse: ${""}")
    }

    let orig = s;
    var copyS = s;
    let wantClose = s[0] == 40

    if (wantClose) {
        if (s[s.size - 1] != 41) {
            throw ParseExpection(information:"quat: failed to parse: ${orig}")
        }

        let strArray = Array(s);
        let strSlice = strArray[1..(s.size - 1)]
        let newStrArray = Array(strSlice)
        copyS = String.fromUtf8(newStrArray)
    }

    if (copyS.size == 0) {
        throw ParseExpection(information:"quat: failed to parse: ${orig}")
    }

    match (copyS[0]) {
        case 110 | 78 => 
        if (copyS.toAsciiLower().compare("nan") == Ordering.EQ) {
                return Complex128(Float64.NaN, Float64.NaN)
        }
        
        case 105 | 73 => 
        if (copyS.toAsciiLower().compare("inf") == Ordering.EQ) {
                return Complex128(Float64.Inf, Float64.Inf)
        }
        
        case _ => ()
    }

    var q = Complex128(0.0) 
    var parts: Byte = 0
    for (i in 0..4) {
    let get = floatPart(copyS)
    let beg = get[0]
    let end = get[1]
    let p = get[2]

    if ((parts & (1 << p)) != 0 ) {
        throw ParseExpection(information:"quat: failed to parse: ${orig}")
    }

    parts = parts | (1 << p)
    var v = 0.0

    match {
        case end == 1 && copyS[0] == 43 =>
            if(end > copyS.size - 1) {
                throw ParseExpection(information:"quat: failed to parse: ${orig}")
            }
            v = 1.0
        case end == 1 && copyS[0] == 45 => 
            if(end > copyS.size - 1) {
                throw ParseExpection(information:"quat: failed to parse: ${orig}")
            }
            v = -1.0
        case _ => 
        v = Float64.parse(String.fromUtf8(Array(copyS)[beg .. end]))
    }

    copyS = String.fromUtf8(Array(copyS)[end ..= (copyS.size - 1)])
    
    match {
        case p == 0 =>
        q = Complex128((q.real + v), q.imag)
        case p == 1 =>
        q = Complex128(q.real, (q.imag + v))
        copyS = String.fromUtf8(Array(copyS)[1 ..= (copyS.size - 1)])
        case _ => ()
    }

    if (copyS.size == 0) {
        return q
    }

    if (!isSign(copyS[0])) {
        throw ParseExpection(information:"quat: failed to parse: ${orig}") 
    }
    }

    throw ParseExpection(information:"quat: failed to parse: ${orig}")
}

func floatPart(s: String): (Int64, Int64, UInt64){
    const wantMantSign = 0
    const wantMantIntInit = 1
    const wantMantInt = 2
    const wantMantFrac = 3
    const wantExpSign = 4
    const wantExpInt = 5
    const wantInfN = 6
    const wantInfF = 7
    const wantCloseInf = 8
    const wantNaNA = 9
    const wantNaNN = 10
    const wantCloseNaN = 11

    var state = 0
    var i = 0
    var beg: Int64 = 0
    for(rune in s){
        match  {    
            case state == wantMantSign =>
                state = match {
                    case isSign(rune) => wantMantIntInit
                    case isDigit(rune) => wantMantInt
                    case isDot(rune) => wantMantFrac
                    case (rune == 105 || rune == 73) => wantInfN
                    case (rune == 110 || rune == 78) => wantNaNA
                    case _ => throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }

            case state == wantMantIntInit => 
                state = match {
                    case isDigit(rune) => wantMantInt
                    case isDot(rune) => wantMantFrac
                    case (rune == 105) =>
                        if(i == s.size - 1 || s[i + 1] == 45 || s[i + 1] == 43){
                            return (0, i, 1)
                        }
                        wantInfN
                    case (rune == 73) => wantInfN
                    case (rune == 110 || rune == 78) => wantNaNA
                    case _ => throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }

            case state == wantMantInt =>
                state = match {
                    case isDigit(rune) => state
                    case isDot(rune) => wantMantFrac
                    case isExponent(rune) => wantExpSign
                    case isSign(rune) => return (0, i, 0)
                    case (rune == 105) => return (0, i, 1)
                    case _ =>  throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }

            case state == wantMantFrac =>
                state = match {
                    case isDigit(rune) => state
                    case isExponent(rune) => wantExpSign
                    case isSign(rune) => return (0, i, 0)
                    case (rune == 105) => return (0, i, 1)
                    case _ => throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }

            case state == wantExpSign =>
                state = match {
                    case (isSign(rune) || isDigit(rune)) => wantExpInt
                    case _ => throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }

            case state == wantExpInt =>
                state = match {
                    case isDigit(rune) => state 
                    case isSign(rune) => return (0, i, 0)
                    case (rune == 105) => return (0, i, 1)
                    case _ => throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }
            
            case state == wantInfN =>
                if(rune != 110 && rune != 78) {
                    throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }
                state = wantInfF

            case state == wantInfF =>
                if(rune != 70 && rune != 102) {
                    throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }
                state = wantCloseInf

            case state == wantCloseInf => 
                match{
                    case isSign(rune) => return (0, i, 0)
                    case (rune == 105) => return (0, i, 1)
                    case _ => throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }

            case state == wantNaNA => 
                if(rune != 65 && rune != 97) {
                    throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }
                state = wantNaNN
            
            case state == wantNaNN =>
                if (rune != 110 && rune != 78) {
                    throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }
                state = wantCloseNaN
            
            case state == wantCloseNaN =>
                if(isSign(s[0])) {
                    beg = 1
                }
                match {
                    case isSign(rune) => return (beg, i, 0)
                    case (rune == 105) => return (beg, i, 1)
                    case _ => throw ParseExpection(information:"quat: failed to parse in state ${state} with ${rune}: ${s}")
                }

            case _ => () 
        }
        i ++ ;
    }

    if(state == wantMantSign || state == wantExpSign || state == wantExpInt) {
        if(state != wantExpInt || !isDigit(s[s.size - 1])) {
            throw ParseExpection(information:"quat: failed to parse in state ${state} with ${Rune(s[s.size - 1])}: ${s}")
        }

    }

    return (0, s.size, 0)
}

func isSign(r: UInt8): Bool {
    return r == 43 || r == 45
}

func isDigit(r: UInt8): Bool {
    return r >= 48 && r <= 57
}

func isExponent(r: UInt8): Bool {
    return r == 101 || r == 69
}

func isDot(r: UInt8): Bool {
    return r == 46
}