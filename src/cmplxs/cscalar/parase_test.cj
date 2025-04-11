package cjnum.cmplxs.cscalar

import std.math.*
import std.math.*
import cjnum.complex.*
import std.unittest.*
import std.unittest.testmacro.*

class parseTestStruct {
    public parseTestStruct(
        let s!: String = "",
        let complex!: Complex128 = Complex128(0.0),
        let expection! : Exception = Exception()
    ) {}
}

let testData1 = [
    parseTestStruct(s:"", expection:ParseExpection(information: "quat: failed to parse: ${""}")),
    parseTestStruct(s:"()", expection:ParseExpection(information: "quat: failed to parse: ()")),
    parseTestStruct(s:"(1", expection:ParseExpection(information: "quat: failed to parse: (1")),
    parseTestStruct(s:"1)", expection:ParseExpection(information: "quat: failed to parse: 1)")),
    parseTestStruct(s:"1+2i+3i", expection:ParseExpection(information: "quat: failed to parse: 1+2i+3i)")),
    parseTestStruct(s:"1e-4i+", expection:ParseExpection(information: "quat: failed to parse: 1e-4i+)")),
    parseTestStruct(s:"1e-4i-", expection:ParseExpection(information: "quat: failed to parse: 1e-4i-)"))
    ]
    
let testData2 = [
    parseTestStruct(s:"1+4i", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"4i+1", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"+1+4i", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"+4i+1", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:".1+.4i", complex:Complex128(0.1, 0.4)),
    parseTestStruct(s:".4i+.1", complex:Complex128(0.1, 0.4)),
    parseTestStruct(s:"+.1+.4i", complex:Complex128(0.1, 0.4)),
    parseTestStruct(s:"+.4i+.1", complex:Complex128(0.1, 0.4)),
    parseTestStruct(s:"1.+4.i", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"4.i+1.", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"+1.+4.i", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"+4.i+1.", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"1.0+4.0i", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"4.0i+1.0", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"+1.0+4.0i", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"+4.0i+1.0", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"1.0e-4+1i", complex:Complex128(1e-4, 1.0)),
    parseTestStruct(s:"1.0e-4+i", complex:Complex128(1e-4, 1.0)),
    parseTestStruct(s:"1.0e-4-i", complex:Complex128(1e-4, -1.0)),
    parseTestStruct(s:"1.0e-4i-1", complex:Complex128(-1.0, 1e-4)),
    parseTestStruct(s:"1.0e-4i+1", complex:Complex128(1.0, 1e-4)),
    parseTestStruct(s:"1e-4+1i", complex:Complex128(1e-4, 1.0)),
    parseTestStruct(s:"1e-4+i", complex:Complex128(1e-4, 1.0)),
    parseTestStruct(s:"1e-4-i", complex:Complex128(1e-4, -1.0)),
    parseTestStruct(s:"1e-4i-1", complex:Complex128(-1.0, 1e-4)),
    parseTestStruct(s:"1e-4i+1", complex:Complex128(1.0, 1e-4)),
    parseTestStruct(s:"(1+4i)", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"(4i+1)", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"(+1+4i)", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"(+4i+1)", complex:Complex128(1.0, 4.0)),
    parseTestStruct(s:"(1e-4+1i)", complex:Complex128(1e-4, 1.0)),
    parseTestStruct(s:"(1e-4+i)", complex:Complex128(1e-4, 1.0)),
    parseTestStruct(s:"(1e-4-i)", complex:Complex128(1e-4, -1.0)),
    parseTestStruct(s:"(1e-4i-1)", complex:Complex128(-1.0, 1e-4)),
    parseTestStruct(s:"(1e-4i+1)", complex:Complex128(1.0, 1e-4)),
    parseTestStruct(s:"Inf", complex:Complex128(Float64.Inf, Float64.Inf)),
    parseTestStruct(s:"inf", complex:Complex128(Float64.Inf, Float64.Inf)),
    parseTestStruct(s:"(Inf+Infi)", complex:Complex128(Float64.Inf, Float64.Inf)),
    parseTestStruct(s:"(-Inf+Infi)", complex:Complex128(-1.0 * Float64.Inf, Float64.Inf)),
    parseTestStruct(s:"(+Inf-Infi)", complex:Complex128(Float64.Inf, -1.0 * Float64.Inf)),
    parseTestStruct(s:"(inf+infi)", complex:Complex128(Float64.Inf, Float64.Inf)),
    parseTestStruct(s:"(-inf+infi)", complex:Complex128(-1.0 * Float64.Inf, Float64.Inf)),
    parseTestStruct(s:"(+inf-infi)", complex:Complex128(Float64.Inf, -1.0 * Float64.Inf)),
    parseTestStruct(s:"NaN", complex:Complex128(Float64.NaN, Float64.NaN)),
    parseTestStruct(s:"nan", complex:Complex128(Float64.NaN, Float64.NaN)),
    parseTestStruct(s:"(nan+nani)", complex:Complex128(Float64.NaN, Float64.NaN)),
    parseTestStruct(s:"(nan-nani)", complex:Complex128(Float64.NaN, Float64.NaN))
    ]

@Test

public func testParse(): Unit {
    for(it in testData1) {
        let _ = @ExpectThrows(cParse(it.s))[it.expection]
    }

    for(it in testData2) {
        let result = cParse(it.s)
        @Expect(same(result, it.complex))
    }
}
