package cjnum.blas.cblas128

import cjnum.complex.*
import cjnum.blas.*
import std.unittest.*
import std.unittest.testmacro.*

public func newTriangularBandFrom(a: TriangularBandCols): TriangularBand {
    let t = TriangularBand(n: a.n, k: a.k, stride: a.k + 1, data: Array<Complex128>(a.n * (a.k +1 ), repeat: Complex128(0.0)), 
        uplo: a.uplo, diag: a.diag)
    for (i in 0..t.data.size) {
        t.data[i] = Complex128.NaN
    }
    t.from(a)
    return t
}

extend TriangularBand <: triangularBand {
    public func tn(): Int64 { 
        return this.n
    }

    public func at(i: Int64, j: Int64): Complex128 {
        if (this.diag == UnitDiag && i == j) {
            return Complex128(1.0)
        }
        let b = Band(rows: this.n, cols: this.n, stride: this.stride, data: this.data)
        if (this.uplo == Upper) {
            if (i > j) {
                return Complex128(0.0)
            }
            b.ku = this.k
        } else if (this.uplo == Lower) {
            if (i < j) {
                return Complex128(0.0)
            }
            b.kl = this.k
        } else {
            throw Exception("cblas128: bad BLAS uplo")
        }
        return b.at(i, j)
    }

    public func bandwidth(): Int64 { 
        return this.k
    }

    public func tuplo(): Uplo { 
        return this.uplo
    }

    public func tdiag(): Diag { 
        return this.diag
    }
}

public func newTriangularBandColsFrom(a: TriangularBand): TriangularBandCols {
    let t = TriangularBandCols(n: a.n, k: a.k, stride: a.k + 1, data: Array<Complex128>(a.n * (a.k + 1), repeat: Complex128(0.0)), 
        uplo: a.uplo, diag: a.diag)
    for (i in 0..t.data.size) {
        t.data[i] = Complex128.NaN
    }
    t.from(a)
    return t
}

extend TriangularBandCols <: triangularBand {
    public func tn(): Int64 { 
        return this.n
    }

    public func at(i: Int64, j: Int64): Complex128 {
        if (this.diag == UnitDiag && i == j) {
            return Complex128(1.0)
        }
        let b = BandCols(rows: this.n, cols: this.n, stride: this.stride, data: this.data)
        if (this.uplo == Upper) {
            if (i > j) {
                return Complex128(0.0)
            }
            b.ku = this.k
        } else if (this.uplo == Lower) {
            if (i < j) {
                return Complex128(0.0)
            }
            b.kl = this.k
        } else {
            throw Exception("cblas128: bad BLAS uplo")
        }
        return b.at(i, j)
    }

    public func bandwidth(): Int64 { 
        return this.k
    }

    public func tuplo(): Uplo { 
        return this.uplo
    }

    public func tdiag(): Diag { 
        return this.diag
    }
}

sealed interface triangularBand {
    func tn(): Int64
    func at(i: Int64, j: Int64): Complex128
    func bandwidth(): Int64
    func tuplo(): Uplo
    func tdiag(): Diag
}

public func sameTriangularBand(a: triangularBand, b: triangularBand): Bool {
    let an = a.tn()
    let bn = b.tn()
    if (an != bn) {
        return false
    }
    if (a.tuplo() != b.tuplo()) {
        return false
    }
    if (a.tdiag() != b.tdiag()) {
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

let triangularBandTests = [
    TriangularBand(n: 3, k: 0, stride: 1, uplo: Upper, data: [
            Complex128(1.0), 
            Complex128(2.0), 
            Complex128(3.0)
        ]
    ),
    TriangularBand(n: 3, k: 0, stride: 1, uplo: Lower, data: [
            Complex128(1.0), 
            Complex128(2.0), 
            Complex128(3.0)
        ]
    ),
    TriangularBand(n: 3, k: 1, stride: 2, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), 
            Complex128(3.0), Complex128(4.0), 
            Complex128(5.0), Complex128(-1.0)
        ]
    ),
    TriangularBand(n: 3, k: 1, stride: 2, uplo: Lower, data: [
            Complex128(-1.0), Complex128(1.0), 
            Complex128(2.0), Complex128(3.0), 
            Complex128(4.0), Complex128(5.0)
        ]
    ),
    TriangularBand(n: 3, k: 2, stride: 3, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0), 
            Complex128(4.0), Complex128(5.0), Complex128(-1.0), 
            Complex128(6.0), Complex128(-2.0), Complex128(-3.0)
        ]
    ),
    TriangularBand(n: 3, k: 2, stride: 3, uplo: Lower, data: [
            Complex128(-2.0), Complex128(-1.0), Complex128(1.0), 
            Complex128(-3.0), Complex128(2.0), Complex128(4.0), 
            Complex128(3.0), Complex128(5.0), Complex128(6.0)
        ]
    ),
    TriangularBand(n: 3, k: 0, stride: 5, uplo: Upper, data: [
            Complex128(1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(2.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(3.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
    TriangularBand(n: 3, k: 0, stride: 5, uplo: Lower, data: [
            Complex128(1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(2.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(3.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
    TriangularBand(n: 3, k: 1, stride: 5, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(3.0), Complex128(4.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(5.0), Complex128(-1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
    TriangularBand(n: 3, k: 1, stride: 5, uplo: Lower, data: [
            Complex128(-1.0), Complex128(1.0), Complex128(0.0), Complex128(0.0), Complex128(0.0), 
            Complex128(2.0), Complex128(3.0), Complex128(0.0), Complex128(0.0), Complex128(0.0),  
            Complex128(4.0), Complex128(5.0), Complex128(0.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
    TriangularBand(n: 3, k: 2, stride: 5, uplo: Upper, data: [
            Complex128(1.0), Complex128(2.0), Complex128(3.0), Complex128(0.0), Complex128(0.0), 
            Complex128(4.0), Complex128(5.0), Complex128(-1.0), Complex128(0.0), Complex128(0.0),  
            Complex128(6.0), Complex128(-2.0), Complex128(-3.0), Complex128(0.0), Complex128(0.0)
        ]
    ),
    TriangularBand(n: 3, k: 2, stride: 5, uplo: Lower, data: [
            Complex128(-2.0), Complex128(-1.0), Complex128(1.0), Complex128(0.0), Complex128(0.0), 
            Complex128(-3.0), Complex128(2.0), Complex128(4.0), Complex128(0.0), Complex128(0.0),  
            Complex128(3.0), Complex128(5.0), Complex128(6.0), Complex128(0.0), Complex128(0.0)
        ]
    )
]

@Test
public func testConvertTriBand() {
    for (test in triangularBandTests) {
        let colmajor = newTriangularBandColsFrom(test)
        @Expect(sameTriangularBand(colmajor, test))
        let rowmajor = newTriangularBandFrom(colmajor)
        @Expect(sameTriangularBand(rowmajor, test))
    }
}