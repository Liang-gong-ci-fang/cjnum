package cjnum.blas.blas32

import cjnum.blas.*
import std.unittest.*
import std.unittest.testmacro.*

public func newTriangularFrom(a: TriangularCols): Triangular {
    let t = Triangular(n: a.n, stride: a.n, data: Array<Float32>(a.n * a.n, repeat: 0.0), diag: a.diag, uplo: a.uplo)
    t.from(a)
    return t
}

extend Triangular <: triangular {
    public func tn(): Int64 {
        return this.n
    }

    public func at(i: Int64, j: Int64): Float32 {
        if (this.diag == UnitDiag && i == j) {
            return 1.0
        }

        if (this.uplo == Lower && i < j && j < this.n) {
            return 0.0
        }

        if (this.uplo == Upper && i > j) {
            return 0.0
        }

        this.data[i * this.stride + j]
    }

    public func tuplo(): Uplo {
        return this.uplo
    }

    public func tdiag(): Diag {
        return this.diag
    }
}

public func newTriangularColsFrom(a: Triangular): TriangularCols {
    let t = TriangularCols(n: a.n, stride: a.n, data: Array<Float32>(a.n * a.n, repeat: 0.0), diag: a.diag, uplo: a.uplo)
    t.from(a)
    return t
}

extend TriangularCols <: triangular {
    public func tn(): Int64 {
        return this.n
    }

    public func at(i: Int64, j: Int64): Float32 {
        if (this.diag == UnitDiag && i == j) {
            return 1.0
        }

        if (this.uplo == Lower && i < j) {
            return 0.0
        }

        if (this.uplo == Upper && i > j && i < this.n) {
            return 0.0
        }

        this.data[i + j * this.stride]
    }

    public func tuplo(): Uplo {
        return this.uplo
    }

    public func tdiag(): Diag {
        return this.diag
    }
}

sealed interface triangular {
    func tn(): Int64
    func at(i: Int64, j: Int64): Float32
    func tuplo(): Uplo
    func tdiag(): Diag
}

func sameTriangular(a: triangular, b: triangular): Bool {
    let an = a.tn()
    let bn = b.tn()
    if (an != bn) {
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

let triangularTests = [Triangular(n: 3, stride: 3, data: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
    Triangular(n: 3, stride: 5, data: [1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0, 0.0, 7.0, 8.0, 9.0, 0.0, 0.0])]

@Test
public func testConvertTriangular() {
    for (test in triangularTests) {
        for (uplo in [Upper, Lower, All]) {
            for (diag in [UnitDiag, NonUnitDiag]) {
                test.uplo = uplo
                test.diag = diag
                let colmajor = newTriangularColsFrom(test)
                @Expect(sameTriangular(colmajor, test))
                let rowmajor = newTriangularFrom(colmajor)
                @Expect(sameTriangular(rowmajor, test))
            }
        }
    }
}
