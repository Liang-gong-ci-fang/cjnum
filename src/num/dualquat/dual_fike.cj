package cjnum.num.dualquat

import cjnum.num.quat.*

// PowReal returns d**p, the base-d exponential of p.
//
// Special cases are (in order):
//
//	PowReal(NaN+xϵ, ±0) = 1+NaNϵ for any x
//	PowReal(x, ±0) = 1 for any x
//	PowReal(1+xϵ, y) = 1+xyϵ for any y
//	PowReal(x, 1) = x for any x
//	PowReal(NaN+xϵ, y) = NaN+NaNϵ
//	PowReal(x, NaN) = NaN+NaNϵ
//	PowReal(±0, y) = ±Inf for y an odd integer < 0
//	PowReal(±0, -Inf) = +Inf
//	PowReal(±0, +Inf) = +0
//	PowReal(±0, y) = +Inf for finite y < 0 and not an odd integer
//	PowReal(±0, y) = ±0 for y an odd integer > 0
//	PowReal(±0, y) = +0 for finite y > 0 and not an odd integer
//	PowReal(-1, ±Inf) = 1
//	PowReal(x+0ϵ, +Inf) = +Inf+NaNϵ for |x| > 1
//	PowReal(x+yϵ, +Inf) = +Inf for |x| > 1
//	PowReal(x, -Inf) = +0+NaNϵ for |x| > 1
//	PowReal(x, +Inf) = +0+NaNϵ for |x| < 1
//	PowReal(x+0ϵ, -Inf) = +Inf+NaNϵ for |x| < 1
//	PowReal(x, -Inf) = +Inf-Infϵ for |x| < 1
//	PowReal(+Inf, y) = +Inf for y > 0
//	PowReal(+Inf, y) = +0 for y < 0
//	PowReal(-Inf, y) = Pow(-0, -y)
public func ndqPowReal(d: NdqNumber, p: Float64): NdqNumber {
	if(p==0.0){
		if(nqIsNaN(d.real)){
			return NdqNumber(NqNumber(1.0,0.0,0.0,0.0),nqNaN())
		}
		else if(d.real.nqEqual(zeroQuat)||nqIsInf(d.real)){
			return NdqNumber(NqNumber(1.0,0.0,0.0,0.0),zeroQuat)
		}
	}
	else if(p==1.0){
		return d
	}
	else if(p.isInf()&&p>0.0){
		if(ndqAbs(d).real>1.0){
            if(d.dual.nqEqual(zeroQuat)){
                return NdqNumber(nqInf(),nqNaN())
            }
            return NdqNumber(nqInf(),nqInf())
        }
        return NdqNumber(zeroQuat,nqNaN())
	}
    else if(p.isInf()&&p<0.0){
        if(ndqAbs(d).real>1.0){
            return NdqNumber(zeroQuat,nqNaN())
        }
        if(d.dual.nqEqual(zeroQuat)){
            return NdqNumber(nqInf(),nqNaN())
        }
        return NdqNumber(nqInf(),nqInf())
    }
	let deriv = nqMul(NqNumber(p,0.0,0.0,0.0), nqPow(d.real, NqNumber(p-1.0,0.0,0.0,0.0)))
	return NdqNumber(
		nqPow(d.real, NqNumber(p,0.0,0.0,0.0)),
		nqMul(d.dual, deriv),
    )
}

/*
以下的TODO部分因为缺乏相关的测试代码支持所以没有实现
*/
//TODO func Pow(d, p Number) Number {

// Sqrt returns the square root of d
//
// Special cases are:
//
//	Sqrt(+Inf) = +Inf
//	Sqrt(±0) = (±0+Infϵ)
//	Sqrt(x < 0) = NaN
//	Sqrt(NaN) = NaN
public func ndqSqrt(d: NdqNumber): NdqNumber {
	return ndqPowReal(d, 0.5)
}

// Exp returns e**d, the base-e exponential of d.
//
// Special cases are:
//
//	Exp(+Inf) = +Inf
//	Exp(NaN) = NaN
//
// Very large values overflow to 0 or +Inf.
// Very small values underflow to 1.
public func ndqExp(d: NdqNumber): NdqNumber {
    let fnDeriv = nqExp(d.real)
    return NdqNumber(
        fnDeriv,
        nqMul(fnDeriv, d.dual),
    )
}

// Log returns the natural logarithm of d.
//
// Special cases are:
//
//	Log(+Inf) = (+Inf+0ϵ)
//	Log(0) = (-Inf±Infϵ)
//	Log(x < 0) = NaN
//	Log(NaN) = NaN
public func ndqLog(d: NdqNumber): NdqNumber {
    if (d.real.nqEqual(zeroQuat)) {
        return NdqNumber(nqLog(d.real), nqInf())
    } else if (nqIsInf(d.real)) {
        return NdqNumber(nqLog(d.real), zeroQuat)
    }
    return NdqNumber(
        nqLog(d.real),
        nqMul(d.dual, nqInv(d.real))
    )
}
