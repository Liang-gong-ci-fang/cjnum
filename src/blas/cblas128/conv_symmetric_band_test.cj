package cjnum.blas.cblas128

import cjnum.blas.*
import cjnum.complex.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

public func newSymmetricBandFrom(a: SymmetricBandCols): SymmetricBand {
    let t = SymmetricBand(n: a.n, k: a.k, stride: a.k + 1, data: Array<Complex128>(a.n * (a.k + 1), repeat: Complex128(0.0)), uplo: a.uplo)
    for (i in 0..t.data.size) {
        t.data[i] = Complex128.NaN
    }
    t.from(a)
    return t
}

extend SymmetricBand <: symmetricBand {
    public func tn(): Int64 {
        return this.n
    }

    public func at(i: Int64, j: Int64): Complex128 {
        var ii = i
        var jj = j
        let b = Band(rows: this.n, cols: this.n, stride: this.stride, data: this.data)
        if (this.uplo == Upper) {
            b.ku = this.k
            if (i > j) {
                ii = j
                jj = i
            }
        } else if (this.uplo == Lower) {
            b.kl = this.k
            if (i < j) {
                ii = j
                jj = i
            }
        } else {
            throw Exception("cblas128: bad BLAS uplo")
        }
        return b.at(ii, jj)
    }

    public func bandwidth(): Int64 { 
        return this.k
    }

    public func tuplo(): Uplo { 
        return this.uplo
    }
}

public func newSymmetricBandColsFrom(a: SymmetricBand): SymmetricBandCols {
    let t = SymmetricBandCols(n: a.n, k: a.k, stride: a.k + 1, data: Array<Complex128>(a.n * (a.k + 1), repeat: Complex128(0.0)), uplo: a.uplo)
    for (i in 0..t.data.size) {
        t.data[i] = Complex128.NaN
    }
    t.from(a)
    return t
}

extend SymmetricBandCols <: symmetricBand {
    public func tn(): Int64 {
        return this.n
    }

    public func at(i: Int64, j: Int64): Complex128 {
        var ii = i
        var jj = j
        let b = BandCols(rows: this.n, cols: this.n, stride: this.stride, data: this.data)
        if (this.uplo == Upper) {
            b.ku = this.k
            if (i > j) {
                ii = j
                jj = i
            }
        } else if (this.uplo == Lower) {
            b.kl = this.k
            if (i < j) {
                ii = j
                jj = i
            }
        } else {
            throw Exception("cblas128: bad BLAS uplo")
        }
        return b.at(ii, jj)
    }

    public func bandwidth(): Int64 { 
        return this.k
    }

    public func tuplo(): Uplo { 
        return this.uplo
    }
}

sealed interface symmetricBand {
	func tn(): Int64
	func at(i: Int64, j: Int64): Complex128
	func bandwidth(): Int64
	func tuplo(): Uplo
}

public func sameSymmetricBand(a: symmetricBand, b: symmetricBand): Bool {
    let an = a.tn()
    let bn = b.tn()
    if (an != bn) {
        return false
    }
    let ak = a.bandwidth()
    let bk = b.bandwidth()
    if (ak != bk) {
        return false
    }
    for (i in 0..an) {
        for (j in 0..an) {
            if (a.at(i, j) != b.at(i, j) || a.at(i, j).isNaN() != b.at(i, j).isNaN()) {
                return false
            }
        }
    }
    return true
}

let symmetricBandTests = [
        SymmetricBand(n: 3, k: 0, stride: 1, uplo: Upper, data: [
            Complex128(1.0), 
            Complex128(2.0), 
            Complex128(3.0)
        ]
    ),
        SymmetricBand(n: 3, k: 0, stride: 1, uplo: Lower, data: [
            Complex128(1.0), 
            Complex128(2.0), 
            Complex128(3.0)
        ]
    ),
        SymmetricBand(n: 3, k: 1, stride: 2, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), 
            Complex128(3.0), Complex128(4.0), 
            Complex128(5.0), Complex128(-1.0)
        ]
    ),
        SymmetricBand(n: 3, k: 1, stride: 2, uplo: Lower, data: [
            Complex128(-1.0), Complex128(1.0), 
            Complex128(2.0), Complex128(3.0), 
            Complex128(4.0), Complex128(5.0)
        ]
    ),
        SymmetricBand(n: 3, k: 2, stride: 3, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0), 
            Complex128(4.0), Complex128(5.0), Complex128(-1.0), 
            Complex128(6.0), Complex128(-2.0), Complex128(-3.0)
        ]
    ),
        SymmetricBand(n: 3, k: 2, stride: 3, uplo: Lower, data: [
            Complex128(-2.0), Complex128(-1.0), Complex128(1.0), 
            Complex128(-3.0), Complex128(2.0), Complex128(4.0), 
            Complex128(3.0), Complex128(5.0), Complex128(6.0)
        ]
    ),
        SymmetricBand(n: 3, k: 0, stride: 5, uplo: Upper, data: [
            Complex128(1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(2.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(3.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
        SymmetricBand(n: 3, k: 0, stride: 5, uplo: Lower, data: [
            Complex128(1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(2.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(3.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
        SymmetricBand(n: 3, k: 1, stride: 5, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(3.0), Complex128(4.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(5.0), Complex128(-1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
        SymmetricBand(n: 3, k: 1, stride: 5, uplo: Lower, data: [
            Complex128(-1.0), Complex128(1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(2.0), Complex128(3.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(4.0), Complex128(5.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
        SymmetricBand(n: 3, k: 2, stride: 5, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0), Complex128(0.0), Complex128(0.0), 
            Complex128(4.0), Complex128(5.0), Complex128(-1.0), Complex128(0.0), Complex128(0.0),  
            Complex128(6.0), Complex128(-2.0), Complex128(-3.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
        SymmetricBand(n: 3, k: 2, stride: 5, uplo: Lower, data: [
            Complex128(-2.0), Complex128(-1.0), Complex128(1.0), Complex128(0.0), Complex128(0.0), 
            Complex128(-3.0), Complex128(2.0), Complex128(4.0), Complex128(0.0), Complex128(0.0),  
            Complex128(3.0), Complex128(5.0), Complex128(6.0), Complex128(0.0), Complex128(0.0)
        ]
    )
]

@Test
public func testConvertSymBand() {
    for (test in symmetricBandTests) {
        let colmajor = newSymmetricBandColsFrom(test)
        @Expect(sameSymmetricBand(colmajor, test))
        let rowmajor = newSymmetricBandFrom(colmajor)
        @Expect(sameSymmetricBand(rowmajor, test))
    }
}