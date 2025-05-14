# cjnum 库

### 介绍

cjnum是一个用于 Cangjie 语言的数值计算库，提供了广泛的数学、科学计算和数值分析功能。

### 1 主要接口

#### 1.1 NFloat64
```cangjie
// Float64类型的blas接口
public interface NFloat64 <: Float64Level1 & Float64Level2 & Float64Level3 {}
```

#### 1.2 Float64Level1 & Float64Level2 & Float64Level3
```cangjie
// 向量-向量运算，处理单向量操作或双向量的逐元素运算
public interface Float64Level1 {
    func ddot(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Float64
    func dnrm2(n: Int64, x: Array<Float64>, incX: Int64): Float64
    func dasum(n: Int64, x: Array<Float64>, incX: Int64): Float64
    func idamax(n: Int64, x: Array<Float64>, incX: Int64): Int64
    func dswap(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit
    func dcopy(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit
    func daxpy(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64): Unit
    func drotg(a: Float64, b: Float64): (Float64, Float64, Float64, Float64)
    func drotmg(d1: Float64, d2: Float64, x1: Float64, y1: Float64): (DrotmParams, Float64, Float64, Float64)
    func drot(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64, c: Float64, s: Float64): Unit
    func drotm(n: Int64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64, p: DrotmParams): Unit
    func dscal(n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64): Unit
}

// 矩阵-向量运算，实现矩阵与向量的线性变换
public interface Float64Level2 {
    func dgemv(tA: Transpose, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dgbmv(tA: Transpose, m: Int64, n: Int64, kL: Int64, kU: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        x: Array<Float64>, incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dtrmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit
    func dtbmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit
    func dtpmv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit
    func dtrsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64): Unit
    func dtbsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64): Unit
    func dtpsv(ul: Uplo, tA: Transpose, d: Diag, n: Int64, ap: Array<Float64>, x: Array<Float64>, incX: Int64): Unit
    func dsymv(ul: Uplo, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>, incX: Int64,
        beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dsbmv(ul: Uplo, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, x: Array<Float64>,
        incX: Int64, beta: Float64, y: Array<Float64>, incY: Int64): Unit
    func dspmv(ul: Uplo, n: Int64, alpha: Float64, ap: Array<Float64>, x: Array<Float64>, incX: Int64, beta: Float64,
        y: Array<Float64>, incY: Int64): Unit
    func dger(m: Int64, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit
    func dsyr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, a: Array<Float64>, lda: Int64): Unit
    func dspr(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, ap: Array<Float64>): Unit
    func dsyr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>, lda: Int64): Unit
    func dspr2(ul: Uplo, n: Int64, alpha: Float64, x: Array<Float64>, incX: Int64, y: Array<Float64>, incY: Int64,
        a: Array<Float64>): Unit
}

// 矩阵-矩阵运算，稠密矩阵乘法及其变体
public interface Float64Level3 {
    func dgemm(tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit
    func dsymm(s: Side, ul: Uplo, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, b: Array<Float64>,
        ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit
    func dsyrk(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, beta: Float64,
        c: Array<Float64>, ldc: Int64): Unit
    func dsyr2k(ul: Uplo, t: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit
    func dtrmm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit
    func dtrsm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit
}
```