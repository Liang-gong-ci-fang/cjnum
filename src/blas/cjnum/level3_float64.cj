package cjnum.blas.cjnum

import cjnum.blas.*

extend Float64Implementation {
    public func dgemm(tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dsymm(s: Side, ul: Uplo, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, b: Array<Float64>,
        ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dsyrk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, beta: Float64,
        c: Array<Float64>, ldc: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dsyr2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtrmm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtrsm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
}
