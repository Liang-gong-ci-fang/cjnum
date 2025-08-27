package cjnum.num.dualquat

import cjnum.num.quat.*
import cjnum.num.dual.*

// Number is a float64 precision dual quaternion. A dual quaternion
// is a hypercomplex number composed of two quaternions, q₀+q₂ϵ,
// where ϵ²=0, but ϵ≠0. Here, q₀ is termed the real and q₂ the dual.
public class NdqNumber {
    let real: NqNumber
    let dual: NqNumber
    public init(real: NqNumber, dual: NqNumber) {
        this.real = real
        this.dual = dual
    }
}

var zeroQuat = NqNumber(0.0, 0.0, 0.0, 0.0)

/*
以下的TODO部分因为缺乏相关的测试代码支持所以没有实现
*/
// TODO func (d Number) Format(fs fmt.State, c rune) {

// TODO func fmtString(fs fmt.State, c rune, prec, width int, wantPlus bool) string {

// TODO func Add(x, y Number) Number {

// TODO func Sub(x, y Number) Number {

// Mul returns the dual product of x and y.
public func ndqMul(x:NdqNumber, y: NdqNumber): NdqNumber {
	return NdqNumber(
		nqMul(x.real, y.real),
		nqAdd(nqMul(x.real, y.dual), nqMul(x.dual, y.real)),
    )
}

// Inv returns the dual inverse of d.
public func ndqInv(d: NdqNumber): NdqNumber {
	return NdqNumber(
		nqInv(d.real),
		nqScale(-1.0, nqMul(d.dual, nqInv(nqMul(d.real, d.real)))),
    )
}

/*
以下的TODO部分因为缺乏相关的测试代码支持所以没有实现
*/
// TODO func Conj(d Number) Number {

// TODO func ConjDual(d Number) Number {

// TODO func ConjQuat(d Number) Number {

// TODO func Scale(f float64, d Number) Number {

// Abs returns the absolute value of d.
func ndqAbs(d: NdqNumber): NdNumber {
	return NdNumber(
		nqAbs(d.real),
		nqAbs(d.dual),
	)
}