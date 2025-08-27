package cjnum.internal.math32

import cjnum.floats.scalar.*
import std.math.*
import std.math.numeric.*
import std.unittest.*
import std.unittest.testmacro.*

const tol: Float64 = 1e-7

@Test
public func testAbs(){ 
    let absTests: Array<Float32> = [0.0, -0.0, 1.0, -1.0, Float32.NaN, 123.4, -123.4]
    for(test in absTests){
        let result = nAbs(test)
        @Expect(same(Float64(result), Float64(abs(test))))
    }
}

@Test
public func testCopySign(){
    let copySignTests: Array<Array<Float32>> = [[1.0, 2.0, 1.0], [-1.0, -2.0, -1.0], [1.0, -2.0, -1.0], [-1.0, -2.0, -1.0],
    [0.0, 2.0, 0.0], [0.0, -2.0, -0.0], [1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]]
    for(test in copySignTests){
        let result = nCopysign(test[0], test[1])
        @Expect(result, test[2])
    }
}

@Test
public func testHypot(){
    let hypotTest: Array<Array<Float32>> = [[3.0, 4.0, 5.0], [-5.0, -12.0, 13.0], [1.0, 1.0, sqrt(2.0)], [0.0, 0.0, 0.0],
        [0.0, 5.0, 5.0], [5.0, 0.0, 5.0], [-3.0, 4.0, 5.0], [3.0, -4.0, 5.0], [-5.0, -12.0, 13.0],
        [1e6, 1e6, sqrt(1e6 * 1e6 + 1e6 * 1e6)], [Float32.Inf, 1.0, Float32.Inf], [1.0, Float32.Inf, Float32.Inf],
        [Float32.NaN, 1.0, Float32.NaN], [1.0, Float32.NaN, Float32.NaN], [Float32.Max, Float32.Max, Float32.Max]]
    for(test in hypotTest){
        let result = nHypot(test[0], test[1])
        if (result > Float32.Max){
            return 
        }
        @Expect(equalWithinRel(Float64(result), Float64(test[2]), tol))  
    }
}

@Test
public func testInf(){
        @Expect(nInf(1), Float32.Inf)
        @Expect(nInf(-1), -Float32.Inf)
}

@Test
public func tsetIsInf(){
    let posInf = Float32.Inf
    let negInf = -Float32.Inf

    @Expect(!nIsInf(posInf, 0) || !nIsInf(negInf, 0) ||
		!nIsInf(posInf, 1) || !nIsInf(negInf, -1) ||
		nIsInf(posInf, -1) || nIsInf(negInf, 1), false)
    @Expect(nIsInf(1.0, 0)||nIsInf(1.0, 0), false)
}

@Test
public func testIsNaN(){
    let isNaNTests: Array<Float32> = [0.0, 1.0, Float32.NaN, -1.0, Float32.Inf, -Float32.Inf]
    for(test in isNaNTests){
        let result = nIsNaN(test)
        @Expect(result, nIsNaN(test))     
    }
}

@Test
public func testSignbit(){
    let signBitTests: Array<Float32> = [1.0, -1.0, 0.0, -0.0, Float32.Max, Float32.Min]
    let signResult: Array<Bool> = [false, true, false, true, false, true]
    var i = 0
    for(test in signBitTests){
        let result = nSignbit(test)  
        @Expect(result, signResult[i]) 
        i += 1
    }
}

@Test
public func testSqrt(){
    let sqrtTests: Array<Float32> = [4.0, 0.0, 1234.5, Float32.Max, Float32.Inf]
    let sqrtTests2: Array<Float32> = [-1.0, Float32.Min, -Float32.Inf]
    for(test in sqrtTests){
        let y = sqrt(test)
        if (nIsNaN(y) && nIsNaN(nSqrt(test))) {
			return
		}
        @Expect(equalWithinRel(Float64(nSqrt(test)), Float64(y), tol))    
    }
    for(test2 in sqrtTests2){
            let result: Float32 = nSqrt(test2)
            @Expect(same(Float64(result), Float64.NaN))  
        }
}

func nSqrt(x: Float32) : Float32 {
    if(x == 0.0 || nIsNaN(x) || nIsInf(x, 1)){
        return x
    }
    if(x < 0.0){
        return Float32.NaN
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

// func alike(a: Float32,b :Float32) : Bool {
//     if(nIsNaN(a) && nIsNaN(b)){
//         return true
//     }
//     if(a == b){
//         return nSignbit(a) == nSignbit(b)
//     }
// 	return false
// }