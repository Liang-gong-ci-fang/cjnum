package cjnum.blas.blas64

import cjnum.blas.*
import std.math.*
import std.unittest.*
import std.unittest.testmacro.*

public func newSymmetricFrom(a: SymmetricCols): Symmetric {
    let t = Symmetric(n: a.n, stride: a.n, data: Array<Float64>(a.n * a.n, item: 0.0), uplo: a.uplo)
    t.from(a)
    return t
}

extend Symmetric <: symmetric {
    public func _n(): Int64 { 
        return this.n
    }

    public func at(i: Int64, j: Int64): Float64 {
        var ii = i
        var jj = j
        if (this.uplo == Lower && i < j && j < this.n) {
            ii = j
            jj = i
        }
        if (this.uplo == Upper && i > j) {
            let temp = ii
            ii = jj
            jj = temp
        }
        return this.data[ii * this.stride + jj]
    }

    public func _uplo(): Uplo { 
        return this.uplo
    }
}

public func newSymmetricColsFrom(a: Symmetric): SymmetricCols {
    let t = SymmetricCols(n: a.n, stride: a.n, data: Array<Float64>(a.n * a.n, item: 0.0), uplo: a.uplo)
    t.from(a)
    return t
}

extend SymmetricCols <: symmetric {
    public func _n(): Int64 { 
        return this.n
    }

    public func at(i: Int64, j: Int64): Float64 {
        var ii = i
        var jj = j
        if (this.uplo == Lower && i < j) {
            ii = j
            jj = i
        }
        if (this.uplo == Upper && i > j && i < this.n) {
            let temp = ii
            ii = jj
            jj = temp
        }
        return this.data[ii + jj * this.stride]
    }

    public func _uplo(): Uplo { 
        return this.uplo
    }
}

sealed interface symmetric {
    func _n(): Int64
    func at(i: Int64, j: Int64): Float64
    func _uplo(): Uplo
}

public func sameSymmetric(a: symmetric, b: symmetric): Bool {
    let an = a._n()
    let bn = b._n()
    if (an != bn) {
        return false
    }
    if (a._uplo() != b._uplo()) {
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

let symmetricTests = [
    Symmetric(n: 3, stride: 3, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
    Symmetric(n: 3, stride: 5, data: [1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 0.0, 7.0, 8.0, 9.0, 0.0, 0.0])
]

@Test
public func testConvertSymmetric() {
    for (test in symmetricTests) {
        for (uplo in [Upper, Lower]) {
            test.uplo = uplo
            let colmajor = newSymmetricColsFrom(test)
            @Expect(sameSymmetric(colmajor, test))
            let rowjajor = newSymmetricFrom(colmajor)
            @Expect(sameSymmetric(rowjajor, test))
        }
    }
}