package cjnum.internal.math32

import cjnum.floats.scalar.*
import std.math.*
import std.math.numeric.*
import std.random.*
import std.unittest.*
import std.unittest.testmacro.*

const tol:Float64=1e-7

let absTests: Array<Float32> = [0.0, -0.0, 1.0, -1.0, Float32(unan), 123.4, -123.4]
let copySignTests: Array<Array<Float32>> = [[1.0, 2.0, 1.0], [-1.0, -2.0, -1.0], [1.0, -2.0, -1.0], [-1.0, -2.0, -1.0],
    [0.0, 2.0, 0.0], [0.0, -2.0, -0.0], [1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]]
let hypotTest: Array<Array<Float32>> = [[3.0, 4.0, 5.0], [-5.0, -12.0, 13.0], [1.0, 1.0, sqrt(2.0)], [0.0, 0.0, 0.0],
    [0.0, 5.0, 5.0], [5.0, 0.0, 5.0], [-3.0, 4.0, 5.0], [3.0, -4.0, 5.0], [-5.0, -12.0, 13.0],
    [1e6, 1e6, sqrt(1e6 * 1e6 + 1e6 * 1e6)], [Float32(uinf), 1.0, Float32(uinf)], [1.0, Float32(uinf), Float32(uinf)],
    [Float32(unan), 1.0, Float32(unan)], [1.0, Float32(unan), Float32(unan)], [Float32.Max, Float32.Max, Float32.Max]]
let isNaNTests: Array<Float32> = [0.0, 1.0, Float32(unan), -1.0, Float32(uinf), Float32(uneginf)]
let signBitTests: Array<Float32> = [1.0, -1.0, 0.0, -0.0, Float32.Max, Float32.Min]
let signResult: Array<Bool> = [false, true, false, true, false, true]
let sqrtTests: Array<Float32> = [4.0, 0.0, 1234.5, Float32.Max, Float32(uneginf), Float32(uneginf)]
let sqrtTests2: Array<Float32> = [-1.0, Float32.Min]

@Test
public func TestAbs(){ 
    for(test in absTests){
        let result = Abs(test)
        @Expect(result, abs(test))
    }
}

@Test
public func TestCopySign(){
    for(test in copySignTests){
        let result = Copysign(test[0], test[1])
        @Expect(result, test[2])
    }
}

@Test
public func TestHypot(){
    for(test in hypotTest){
        let result = Hypot(test[0], test[1])
        if (result > Float32.Max){
            return 
        }
        @Expect(equalWithinRel(Float64(result), Float64(test[2]), tol))
        
    }
}

@Test
func TestInf(){
        @Expect(Inf(1), Float32(uinf))
        @Expect(Inf(-1), Float32(uneginf))
}

@Test
func TsetIsInf(){
    let posInf = Float32(uinf)
    let negInf = Float32(uneginf)

    @Expect(!IsInf(posInf, 0) || !IsInf(negInf, 0) ||
		!IsInf(posInf, 1) || !IsInf(negInf, -1) ||
		IsInf(posInf, -1) || IsInf(negInf, 1), false)
    @Expect(IsInf(1.0, 0)||IsInf(1.0, 0), false)
}

@Test
func TestIsNaN(){
    for(test in isNaNTests){
        let result = IsNaN(test)
        @Expect(result, IsNaN(test))     
    }
}

@Test
func TestSignbit(){
    var i=0
    for(test in signBitTests){
        let result = Signbit(test)  
        @Expect(result, signResult[i]) 
        i+=1
    }
}

@Test
func TestSqrt(){
    for(test in sqrtTests){
        let y = sqrt(test)
        if (IsNaN(y) && IsNaN(Sqrt(test))) {
			return
		}
        @Expect(equalWithinRel(Float64(Sqrt(test)), Float64(y), tol))    
    }
    for(test2 in sqrtTests2){
            let result: Float32 = Sqrt(test2)
            @Expect(result, Float32(unan))  
        }
}

func Sqrt(x: Float32) : Float32 {
    if(x == 0.0 || IsNaN(x) || IsInf(x, 1)){
        return x
    }
    if(x < 0.0){
        return Float32(unan)
    }
    var ix = x.toBits()
    var exp = Int64(ix >> shift) & mask
    if(exp == 0){
        while((ix & (1 << shift)) == 0) {
            ix = ix << 1
            exp -= 1
        }
        exp+=1
    }
    exp -= bias 
    ix = ix & UInt32((mask << shift) ^ 0xFFFFFFFF)
    ix = ix | (1 << shift)
    if ((exp&1) == 1) { 
		ix = ix << 1
	}
    exp >>= 1 
    ix <<= 1
	var q : UInt32 = 0
    var s : UInt32 = 0
    var r : UInt32= 1 << (shift + 1)
    while (r != 0) {
        let t = s + r
        if (t <= ix) {
            s = t + r
            ix -= t
            q += r
        }
        ix <<= 1
        r >>= 1
    }
    if (ix != 0) { 
		q += q & 1 
	}
	ix = (q >> 1) + (UInt32(exp - 1 + bias) << shift)
    return Float32.fromBits(ix)
}

func alike(a: Float32,b :Float32) : Bool {
    if(IsNaN(a) && IsNaN(b)){
        return true
    }
    if(a == b){
        return Signbit(a) == Signbit(b)
    }
	return false
}