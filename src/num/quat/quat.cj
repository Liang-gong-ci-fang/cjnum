package cjnum.num.quat

import std.math.*
import std.convert.*

public var zero = NqNumber(0.0, 0.0, 0.0, 0.0)

// NqNumber is a float64 precision quaternion.
public class NqNumber {
    /*
    因为Number类的对象属性需要修改，所以这里把对象属性改为了可变的var类型
     */
    public var real: Float64
    public var imag: Float64
    public var jmag: Float64
    public var kmag: Float64
    public init(real: Float64, imag: Float64, jmag: Float64, kmag: Float64) {
        this.real = real
        this.imag = imag
        this.jmag = jmag
        this.kmag = kmag
    }
    /*
    因为仓颉语言不能像go语言一样可以直接比较类的值相等，所以这里
    提供了equal函数用于仓颉比较类的值相等
     */
    public func nqEqual(right: NqNumber): Bool {
        if (this.real == right.real && this.imag == right.imag && this.jmag == right.jmag && this.kmag == right.kmag) {
            return true
        }
        return false
    }
}

// Add returns the sum of x and y.
public func nqAdd(x: NqNumber, y: NqNumber): NqNumber {
    return NqNumber(
        x.real + y.real,
        x.imag + y.imag,
        x.jmag + y.jmag,
        x.kmag + y.kmag,
    )
}

// Sub returns the difference of x and y, x-y.
public func nqSub(x: NqNumber, y: NqNumber): NqNumber {
    return NqNumber(
        x.real - y.real,
        x.imag - y.imag,
        x.jmag - y.jmag,
        x.kmag - y.kmag,
    )
}

// Mul returns the Hamiltonian product of x and y.
public func nqMul(x: NqNumber, y: NqNumber): NqNumber {
    return NqNumber(
        x.real * y.real - x.imag * y.imag - x.jmag * y.jmag - x.kmag * y.kmag,
        x.real * y.imag + x.imag * y.real + x.jmag * y.kmag - x.kmag * y.jmag,
        x.real * y.jmag - x.imag * y.kmag + x.jmag * y.real + x.kmag * y.imag,
        x.real * y.kmag + x.imag * y.jmag - x.jmag * y.imag + x.kmag * y.real,
    )
}

// Scale returns q scaled by f.
public func nqScale(f: Float64, q: NqNumber): NqNumber {
    return NqNumber(f * q.real, f * q.imag, f * q.jmag, f * q.kmag)
}

// Parse converts the string s to a NqNumber. The string may be parenthesized and
// has the format [±]N±Ni±Nj±Nk. The order of the components is not strict.
public func nqParse(s_: String): NqNumber {
    if (s_.size == 0) {
        throw Exception("使用nqParse方法时输入字符串为空")
        // return NqNumber(0.0,0.0,0.0,0.0)
    }
    let orig = s_
    var s: Array<Rune> = s_.toRuneArray()
    let wantClose = (s[0] == r'(')
    if (wantClose) {
        if (s[s.size - 1] != r')') {
            throw Exception("输入字符串没有被表达式完整包围")
            // return NqNumber(0.0,0.0,0.0,0.0)
        }
        s = s.slice(1, s.size - 2)
    }
    if (s.size == 0) {
        throw Exception("输入字符串为空")
        // return NqNumber(0.0,0.0,0.0,0.0)
    }
    var judge = String(s)
    if (s[0] == r'n' || s[0] == r'N') {
        if (judge.toAsciiLower() == "nan") {
            return nqNaN()
        }
    } else if (s[0] == r'i' || s[0] == r'I') {
        if (judge.toAsciiLower() == "inf") {
            return nqInf()
        }
    }

    var q: NqNumber = NqNumber(0.0, 0.0, 0.0, 0.0)
    var parts: Byte = 0
    for (_ in 0..4) {
        var (beg, end, p) = floatPart(judge)
        if ((UInt64(parts) & (1 << p)) != 0) {
            throw Exception("四个部分(实部和三个虚部)存在重复解析")
            // return q
        }
        parts |= 1 << p
        var v: Float64 = 0.0
        if (judge[..end] == "-") {
            if (judge[end..].toRuneArray().size == 0) {
                throw Exception("负号后面是无效内容")
                // return q
            }
            v = -1.0
        } else if (judge[..end] == "+") {
            if (judge[end..].toRuneArray().size == 0) {
                throw Exception("正号后面是无效内容")
                // return q
            }
            v = 1.0
        } else {
            v = Float64.parse(judge[beg..end])
        }
        judge = judge[end..]
        if (p == 0) {
            q.real = v
        } else if (p == 1) {
            q.imag = v
            judge = judge[1..]
        } else if (p == 2) {
            q.jmag = v
            judge = judge[1..]
        } else if (p == 3) {
            q.kmag = v
            judge = judge[1..]
        }
        if (judge.size == 0) {
            return q
        }
        if (!isSign(Rune(judge[0]))) {
            throw Exception("输入字符串的第一个字符${judge[0]}不是有效符号")
            // return q
        }
    }
    throw Exception("四元数解析过程中出现格式错误")
    // return q
}

private func floatPart(s_: String): (Int64, Int64, UInt) {
    var s: Array<Rune> = s_.toRuneArray()

    let wantMantSign = 0
    let wantMantIntInit = 1
    let wantMantInt = 2
    let wantMantFrac = 3
    let wantExpSign = 4
    let wantExpInt = 5

    let wantInfN = 6
    let wantInfF = 7
    let wantCloseInf = 8

    let wantNaNA = 9
    let wantNaNN = 10
    let wantCloseNaN = 11

    var i: Int64 = 0
    var state: Int64 = 0
    var r: Rune = r'0'
    for (r in s) {
        if (state == wantMantSign) {
            if (isSign(r)) {
                state = wantMantIntInit
            } else if (isDigit(r)) {
                state = wantMantInt
            } else if (isDot(r)) {
                state = wantMantFrac
            } else if (r == r'i' || r == r'I') {
                state = wantInfN
            } else if (r == r'n' || r == r'N') {
                state = wantNaNA
            } else {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)
            }
        } else if (state == wantMantIntInit) {
            if (isDigit(r)) {
                state = wantMantInt
            } else if (isDot(r)) {
                state = wantMantFrac
            } else if (r == r'i' || r == r'I') {
                // We need to sneak a look-ahead here.
                if (r == r'i' && (i == s.size - 1 || s[i + 1] == r'-' || s[i + 1] == r'+')) {
                    return (0, i, 1)
                }
                state = wantInfN
            } else if (r == r'n' || r == r'N') {
                state = wantNaNA
            } else {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)
            }
        } else if (state == wantMantInt) {
            if (isDigit(r)) {
            // Do nothing
            } else if (isDot(r)) {
                state = wantMantFrac
            } else if (isExponent(r)) {
                state = wantExpSign
            } else if (isSign(r)) {
                return (0, i, 0)
            } else if (r == r'i') {
                return (0, i, 1)
            } else if (r == r'j') {
                return (0, i, 2)
            } else if (r == r'k') {
                return (0, i, 3)
            } else {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
        } else if (state == wantMantFrac) {
            if (isDigit(r)) {
            // Do nothing
            } else if (isExponent(r)) {
                state = wantExpSign
            } else if (isSign(r)) {
                return (0, i, 0)
            } else if (r == r'i') {
                return (0, i, 1)
            } else if (r == r'j') {
                return (0, i, 2)
            } else if (r == r'k') {
                return (0, i, 3)
            } else {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
        } else if (state == wantExpSign) {
            if (isSign(r) || isDigit(r)) {
                state = wantExpInt
            } else {
                Exception("quat: failed to parse in state ${state} with ${s_}")
                return (0, i, 0)
            }
        } else if (state == wantExpInt) {
            if (isDigit(r)) {
            // Do nothing
            } else if (isSign(r)) {
                return (0, i, 0)
            } else if (r == r'i') {
                return (0, i, 1)
            } else if (r == r'j') {
                return (0, i, 2)
            } else if (r == r'k') {
                return (0, i, 3)
            } else {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
        } else if (state == wantInfN) {
            if (r != r'n' && r != r'N') {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
            state = wantInfF
        } else if (state == wantInfF) {
            if (r != r'f' && r != r'F') {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
            state = wantCloseInf
        } else if (state == wantCloseInf) {
            if (isSign(r)) {
                return (0, i, 0)
            } else if (r == r'i') {
                return (0, i, 1)
            } else if (r == r'j') {
                return (0, i, 2)
            } else if (r == r'k') {
                return (0, i, 3)
            } else {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
        } else if (state == wantNaNA) {
            if (r != r'a' && r != r'A') {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
            state = wantNaNN
        } else if (state == wantNaNN) {
            if (r != r'n' && r != r'N') {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (0,i,0)                
            }
            state = wantCloseNaN
        } else if (state == wantCloseNaN) {
            var beg: Int64 = 0
            if (isSign(s[0])) {
                beg = 1
            }
            if (isSign(r)) {
                return (beg, i, 0)
            } else if (r == r'i') {
                return (beg, i, 1)
            } else if (r == r'j') {
                return (beg, i, 2)
            } else if (r == r'k') {
                return (beg, i, 3)
            } else {
                throw Exception("quat: failed to parse in state ${state} with ${s_}")
                // return (beg,i,0)                
            }
        }
        i += 1
    }
    if (state == wantMantSign || state == wantExpSign || state == wantExpInt) {
        if (!(state == wantExpInt && isDigit(r))) {
            throw Exception("quat: failed to parse in state ${state} with ${s_}")
            // return (0,i,0)  
        }
    }
    return (0, s.size, 0)
}

private func isSign(r: Rune): Bool {
    return r == r'+' || r == r'-'
}

private func isDigit(r: Rune): Bool {
    return r'0' <= r && r <= r'9'
}

private func isExponent(r: Rune): Bool {
    return r == r'e' || r == r'E'
}

private func isDot(r: Rune): Bool {
    return r == r'.'
}

public class ParseError {
    let string: String
    let state: Int64
    let rune: Rune
    public init(string: String, state: Int64, rune: Rune) {
        this.string = string
        this.state = state
        this.rune = rune
    }
    public func error(): String {
        if (this.state < 0) {
            throw Exception("quat: failed to parse: ${this.string}")
            // return this.string
        }
        throw Exception("quat: failed to parse in state ${this.state} with ${this.rune}: ${this.string}")
        // return this.string
    }
}
