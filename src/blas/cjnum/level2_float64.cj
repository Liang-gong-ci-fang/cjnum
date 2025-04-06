package cjnum.blas.cjnum

import cjnum.blas.*

extend NFloat64Implementation {
    public func dgemv(tA: Transpose, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dgbmv(tA: Transpose, m: Int64, n: Int64, kL: Int64, kU: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        x: Array<Float64>, incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtrmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtbmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtpmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtrsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtbsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dtpsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dsymv(ul: Uplo, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64,
        beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dsbmv(ul: Uplo, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dspmv(ul: Uplo, n: Int64, alpha: Float64, ap: Array<Float64>, x: Array<Float64>, incX: Int64, beta: Float64,
        y: Array<Float64>, incY: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dger(m: Int64, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dsyr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, a: Array<Float64>, lda: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dspr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, ap: Array<Float64>): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dsyr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit {
        //TODO
        throw Exception("Not implemented")
    }
    
    public func dspr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>): Unit {
        //TODO
        throw Exception("Not implemented")
    }
}