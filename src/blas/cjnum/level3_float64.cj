package cjnum.blas.cjnum

import cjnum.internal.asm.f64.*
import cjnum.blas.*
import std.collection.*
import std.math.*
import std.sync.*

extend NFloat64Implementation {
    public func dtrsm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit {
        if (s != Left && s != Right) {
            throw Exception(badSide)
        }
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }

        let k = if (s == Left) {
            m
        } else {
            n
        }

        if (lda < max(1, k)) {
            throw Exception(badLdA)
        }
        if (ldb < max(1, n)) {
            throw Exception(badLdB)
        }

        if (m == 0 || n == 0) {
            return 
        }

        if (a.size < lda * (k - 1) + k) {
            throw Exception(shortA)
        }
        if (b.size < ldb * (m - 1) + n) {
            throw Exception(shortB)
        }

        if (alpha == 0.0) {
            for (i in 0..m) {
                let btmp = b[i * ldb..i * ldb + n]
                for (j in 0..btmp.size) {
                    btmp[j] = 0.0
                }
            }
            return 
        }

        let nonUnitDiag = d == NonUnitDiag
        if (s == Left) {
            if (tA == NoTrans) {
                if (ul == Upper) {
                    var i = m - 1
                    while (i >= 0) {
                        let btmp = b[i * ldb..i * ldb + n]
                        if (alpha != 1.0) {
                            scalUnitary(alpha, btmp)
                        }
                        let atmp = a[i * lda + i + 1..i * lda + m]
                        for (ka in 0..atmp.size) {
                            let va = atmp[ka]
                            if (va != 0.0) {
                                let k = ka + i + 1
                                axpyUnitary(-va, b[k * ldb..k * ldb + n], btmp)
                            }
                        }
                        if (nonUnitDiag) {
                            let tmp = 1.0 / a[i * lda + i]
                            scalUnitary(tmp, btmp)
                        }
                        i--
                    }
                    return 
                }

                for (i in 0..m) {
                    let btmp = b[i * ldb..i * ldb + n]
                    if (alpha != 1.0) {
                        scalUnitary(alpha, btmp)
                    }
                    let atmp = a[i * lda..i * lda + i]
                    for (k in 0..atmp.size) {
                        if (atmp[k] != 0.0) {
                            axpyUnitary(-atmp[k], b[k * ldb..k * ldb + n], btmp)
                        }
                    }
                    if (nonUnitDiag) {
                        let tmp = 1.0 / a[i * lda + i]
                        scalUnitary(tmp, btmp)
                    }
                }
                return 
            }

            if (ul == Upper) {
                for (k in 0..m) {
                    let btmpk = b[k * ldb..k * ldb + n]
                    if (nonUnitDiag) {
                        scalUnitary(1.0 / a[k * lda + k], btmpk)
                    }
                    let atmpk = a[k * lda + k + 1..k * lda + m]
                    for (ia in 0..atmpk.size) {
                        let va = atmpk[ia]
                        if (va != 0.0) {
                            let i = ia + k + 1
                            axpyUnitary(-va, btmpk, b[i * ldb..i * ldb + n])
                        }
                    }
                    if (alpha != 1.0) {
                        scalUnitary(alpha, btmpk)
                    }
                }
                return 
            }

            var k = m - 1
            while (k >= 0) {
                let btmpk = b[k * ldb..k * ldb + n]
                if (nonUnitDiag) {
                    scalUnitary(1.0 / a[k * lda + k], btmpk)
                }
                let atmpk = a[k * lda..k * lda + k]
                for (i in 0..atmpk.size) { 
                    let va = atmpk[i]
                    if (va != 0.0) {
                        axpyUnitary(-va, btmpk, b[i * ldb..i * ldb + n])
                    }
                }
                if (alpha != 1.0) {
                    scalUnitary(alpha, btmpk)
                }
                k--
            }
            return 
        }

        if (tA == NoTrans) {
            if (ul == Upper) {
                for (i in 0..m) {
                    let btmp = b[i * ldb..i * ldb + n]
                    if (alpha != 1.0) {
                        scalUnitary(alpha, btmp)
                    }
                    for (k in 0..btmp.size) {
                        let vb = btmp[k]
                        if (vb == 0.0) {
                            continue
                        }
                        if (nonUnitDiag) {
                            btmp[k] /= a[k * lda + k]
                        }
                        axpyUnitary(-btmp[k], a[k * lda + k + 1..k * lda + n], btmp[k + 1..n])
                    }
                }
                return 
            }

            for (i in 0..m) {
                let btmp = b[i * ldb..i * ldb + n]
                if (alpha != 1.0) {
                    scalUnitary(alpha, btmp)
                }
                var k = n - 1
                while (k >= 0) {
                    if (btmp[k] == 0.0) {
                        continue
                    }
                    if (nonUnitDiag) {
                        btmp[k] /= a[k * lda + k]
                    }
                    axpyUnitary(-btmp[k], a[k * lda..k * lda + k], btmp[..k])
                    k--
                }
            }
            return 
        }

        if (ul == Upper) {
            for (i in 0..m) {
                let btmp = b[i * ldb..i * ldb + n]
                var j = n - 1
                while (j >= 0) {
                    var tmp = alpha * btmp[j] - dotUnitary(a[j * lda + j + 1..j * lda + n], btmp[j + 1..])
                    if (nonUnitDiag) {
                        tmp /= a[j * lda + j]
                    }
                    btmp[j] = tmp
                    j--
                }
            }
            return
        }
        for (i in 0..m) {
            let btmp = b[i * ldb..i * ldb + n]
            for (j in 0..n) {
                var tmp = alpha * btmp[j] - dotUnitary(a[j * lda..j * lda + j], btmp[..j])
                if (nonUnitDiag) {
                    tmp /= a[j * lda + j]
                }
                btmp[j] = tmp
            }
        }
    }

    public func dsymm(s: Side, ul: Uplo, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, b: Array<Float64>,
        ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit {
        if (s != Left && s != Right) {
            throw Exception(badSide)
        }
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }

        let k = if (s == Left) {
            m
        } else {
            n
        }

        if (lda < max(1, k)) {
            throw Exception(badLdA)
        }
        if (ldb < max(1, n)) {
            throw Exception(badLdB)
        }
        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (m == 0 || n == 0) {
            return 
        }

        if (a.size < lda * (k - 1) + k) {
            throw Exception(shortA)
        }
        if (b.size < ldb * (m - 1) + n) {
            throw Exception(shortB)
        }
        if (c.size < ldc * (m - 1) + n) {
            throw Exception(shortC)
        }

        if (alpha == 0.0 && beta == 1.0) {
            return 
        }

        if (beta == 0.0) {
            for (i in 0..m) {
                let ctmp = c[i * ldc..i * ldc + n]
                for (j in 0..ctmp.size) {
                    ctmp[j] = 0.0
                }
            }
        }

        if (alpha == 0.0) {
            if (beta != 0.0) {
                for (i in 0..m) {
                    let ctmp = c[i * ldc..i * ldc + n]
                    for (j in 0..n) {
                        ctmp[j] *= beta
                    }
                }
            }
            return
        }

        let isUpper = ul == Upper
        if (s == Left) {
            for (i in 0..m) {
                let atmp = alpha * a[i * lda + i]
                let btmp = b[i * ldb..i * ldb + n]
                let ctmp = c[i * ldc..i * ldc + n]
                for (j in 0..btmp.size) {
                    ctmp[j] *= beta
                    ctmp[j] += atmp * btmp[j]
                }

                for (k in 0..i) {
                    var atmp = if (isUpper) {
                        a[k * lda + i]
                    } else {
                        a[i * lda + k]
                    }
                    atmp *= alpha
                    axpyUnitary(atmp, b[k * ldb..k * ldb + n], ctmp)
                }
                for (k in i + 1..m) {
                    var atmp = if (isUpper) {
                        a[i * lda + k]
                    } else {
                        a[k * lda + i]
                    }
                    atmp *= alpha
                    axpyUnitary(atmp, b[k * ldb..k * ldb + n], ctmp)
                }
            }
            return 
        }

        if (isUpper) {
            for (i in 0..m) {
                var j = n - 1
                while (j >= 0) {
                    let tmp = alpha * b[i * ldb + j]
                    var tmp2 = 0.0
                    let atmp = a[j * lda + j + 1..j * lda + n]
                    let btmp = b[i * ldb + j + 1..i * ldb + n]
                    let ctmp = c[i * ldc + j + 1..i * ldc + n]
                    for (k in 0..atmp.size) {
                        ctmp[k] += tmp * atmp[k]
                        tmp2 += btmp[k] * atmp[k]
                    }
                    c[i * ldc + j] *= beta
                    c[i * ldc + j] += tmp * a[j * lda + j] + alpha * tmp2
                    j--
                }
            }
            return
        }

        for (i in 0..m) {
            for (j in 0..n) {
                let tmp = alpha * b[i * ldb + j]
                var tmp2 = 0.0
                let atmp = a[j * lda..j * lda + j]
                let btmp = b[i * ldb..i * ldb + j]
                let ctmp = c[i * ldc..i * ldc + j]
                for (k in 0..atmp.size) {
                    ctmp[k] += tmp * atmp[k]
                    tmp2 += btmp[k] * atmp[k]
                }
                c[i * ldc + j] *= beta
                c[i * ldc + j] += tmp * a[j * lda + j] + alpha * tmp2
            }
        }
    }

    public func dsyrk(ul: Uplo, tA: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64, beta: Float64,
        c: Array<Float64>, ldc: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != Trans && tA != NoTrans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kLT0)
        }

        let row = if (tA == NoTrans) {
            n
        } else {
            k
        }
        let col = if (tA == NoTrans) {
            k
        } else {
            n
        }

        if (lda < max(1, col)) {
            throw Exception(badLdA)
        }
        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (n == 0) {
            return 
        }

        if (a.size < lda * (row - 1) + col) {
            throw Exception(shortA)
        }
        if (c.size < ldc * (n - 1) + n) {
            throw Exception(shortC)
        }

        if (alpha == 0.0) {
            if (beta == 0.0) {
                if (ul == Upper) {
                    for (i in 0..n) {
                        let ctmp = c[i * ldc + i..i * ldc + n]
                        for (j in 0..ctmp.size) {
                            ctmp[j] = 0.0
                        }
                    }
                    return 
                }

                for (i in 0..n) {
                    let ctmp = c[i * ldc..i * ldc + i + 1]
                    for (j in 0..ctmp.size) {
                        ctmp[j] = 0.0
                    }
                }
                return 
            }

            if (ul == Upper) {
                for (i in 0..n) {
                    let ctmp = c[i * ldc + i..i * ldc + n]
                    for (j in 0..ctmp.size) {
                        ctmp[j] *= beta
                    }
                }
                return 
            }

            for (i in 0..n) {
                let ctmp = c[i * ldc..i * ldc + i + 1]
                for (j in 0..ctmp.size) {
                    ctmp[j] *= beta
                }
            }
            return 
        }

        if (tA == NoTrans) {
            if (ul == Upper) {
                for (i in 0..n) {
                    let ctmp = c[i * ldc + i..i * ldc + n]
                    let atmp = a[i * lda..i * lda + k]
                    if (beta == 0.0) {
                        for (jc in 0..ctmp.size) {
                            let j = jc + i
                            ctmp[jc] = alpha * dotUnitary(atmp, a[j * lda..j * lda + k])
                        }
                    } else {
                        for (jc in 0..ctmp.size) {
                            let vc = ctmp[jc]
                            let j = jc + i
                            ctmp[jc] = vc * beta + alpha * dotUnitary(atmp, a[j * lda..j * lda + k])
                        }
                    }
                }
                return 
            }

            for (i in 0..n) {
                let ctmp = c[i * ldc..i * ldc + i + 1]
                let atmp = a[i * lda..i * lda + k]
                if (beta == 0.0) {
                    for (j in 0..ctmp.size) {
                        ctmp[j] = alpha * dotUnitary(a[j * lda..j * lda + k], atmp)
                    }
                } else {
                    for (j in 0..ctmp.size) {
                        let vc = ctmp[j]
                        ctmp[j] = vc * beta + alpha * dotUnitary(a[j * lda..j * lda + k], atmp)
                    }
                }
            }
            return 
        }

        if (ul == Upper) {
            for (i in 0..n) {
                let ctmp = c[i * ldc + i..i * ldc + n]
                if (beta == 0.0) {
                    for (j in 0..ctmp.size) {
                        ctmp[j] = 0.0
                    }
                } else if (beta != 1.0) {
                    for (j in 0..ctmp.size) {
                        ctmp[j] *= beta
                    }
                }
                for (l in 0..k) {
                    let tmp = alpha * a[l * lda + i]
                    if (tmp != 0.0) {
                        axpyUnitary(tmp, a[l * lda + i..l * lda + n], ctmp)
                    }
                }
            }
            return 
        }

        for (i in 0..n) {
            let ctmp = c[i * ldc..i * ldc + i + 1]
            if (beta != 1.0) {
                for (j in 0..ctmp.size) {
                    ctmp[j] *= beta
                }
            }
            for (l in 0..k) {
                let tmp = alpha * a[l * lda + i]
                if (tmp != 0.0) {
                    axpyUnitary(tmp, a[l * lda..l * lda + i + 1], ctmp)
                }
            }
        }
    }
    
    public func dsyr2k(ul: Uplo, tA: Transpose, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>, lda: Int64,
        b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit {
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != Trans && tA != NoTrans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kLT0)
        }

        let row = if (tA == NoTrans) {
            n  
        } else {
            k
        }
        let col = if (tA == NoTrans) {
            k  
        } else {
            n
        }

        if (lda < max(1, col)) {
            throw Exception(badLdA)
        }
        if (ldb < max(1, col)) {
            throw Exception(badLdB)
        }
        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (n == 0) {
            return 
        }

        if (a.size < lda * (row - 1) + col) {
            throw Exception(shortA)
        }
        if (b.size < ldb * (row - 1) + col) {
            throw Exception(shortB)
        }
        if (c.size < ldc * (n - 1) + n) {
            throw Exception(shortC)
        }

        if (alpha == 0.0) {
            if (beta == 0.0) {
                if (ul == Upper) {
                    for (i in 0..n) {
                        let ctmp = c[i * ldc + i..i * ldc + n]
                        for (j in 0..ctmp.size) {
                            ctmp[j] = 0.0
                        }
                    }
                    return 
                }

                for (i in 0..n) {
                    let ctmp = c[i * ldc..i * ldc + i + 1]
                    for (j in 0..ctmp.size) {
                        ctmp[j] = 0.0
                    }
                }
                return 
            }

            if (ul == Upper) {
                for (i in 0..n) {
                    let ctmp = c[i * ldc + i..i * ldc + n]
                    for (j in 0..ctmp.size) {
                        ctmp[j] *= beta
                    }
                }
                return 
            }

            for (i in 0..n) {
                let ctmp = c[i * ldc..i * ldc + i + 1]
                for (j in 0..ctmp.size) {
                    ctmp[j] *= beta
                }
            }
            return 
        }

        if (tA == NoTrans) {
            if (ul == Upper) {
                for (i in 0..n) {
                    let atmp = a[i * lda..i * lda + k]
                    let btmp = b[i * ldb..i * ldb + k]
                    let ctmp = c[i * ldc + i..i * ldc + n]
                    if (beta == 0.0) {
                        for (jc in 0..ctmp.size) {
                            let j = i + jc
                            var tmp1 = 0.0
                            var tmp2 = 0.0
                            let binner = b[j * ldb..j * ldb + k]
                            let ainner = a[j * lda..j * lda + k]
                            for(l in 0..ainner.size) {
                                tmp1 += ainner[l] * btmp[l]
                                tmp2 += atmp[l] * binner[l]
                            }
                            ctmp[jc] = alpha * (tmp1 + tmp2)
                        }
                    } else {
                        for (jc in 0..ctmp.size) {
                            let j = i + jc
                            var tmp1 = 0.0
                            var tmp2 = 0.0
                            let binner = b[j * ldb..j * ldb + k]
                            let ainner = a[j * lda..j * lda + k]
                            for(l in 0..ainner.size) {
                                tmp1 += ainner[l] * btmp[l]
                                tmp2 += atmp[l] * binner[l]
                            }
                            ctmp[jc] *= beta
                            ctmp[jc] += alpha * (tmp1 + tmp2)
                        }
                    }
                }
                return 
            }

            for (i in 0..n) {
                let atmp = a[i * lda..i * lda + k]
                let btmp = b[i * ldb..i * ldb + k]
                let ctmp = c[i * ldc..i * ldc + i + 1]
                if (beta == 0.0) {
                    for (j in 0..=i) {
                        var tmp1 = 0.0
                        var tmp2 = 0.0
                        let binner = b[j * ldb..j * ldb + k]
                        let ainner = a[j * lda..j * lda + k]
                        for(l in 0..ainner.size) {
                            tmp1 += ainner[l] * btmp[l]
                            tmp2 += atmp[l] * binner[l]
                        }
                        ctmp[j] = alpha * (tmp1 + tmp2)
                    }
                } else {
                    for (j in 0..=i) {
                        var tmp1 = 0.0
                        var tmp2 = 0.0
                        let binner = b[j * ldb..j * ldb + k]
                        let ainner = a[j * lda..j * lda + k]
                        for(l in 0..ainner.size) {
                            tmp1 += ainner[l] * btmp[l]
                            tmp2 += atmp[l] * binner[l]
                        }
                        ctmp[j] *= beta
                        ctmp[j] += alpha * (tmp1 + tmp2)
                    }
                }
            }
            return 
        }

        if (ul == Upper) {
            for (i in 0..n) {
                let ctmp = c[i * ldc + i..i * ldc + n]
                if (beta == 0.0) {
                    for (j in 0..ctmp.size) {
                        ctmp[j] = 0.0
                    }
                } else if (beta == 1.0) {

                } else {
                    for (j in 0..ctmp.size) {
                        ctmp[j] *= beta
                    }
                }
                for (l in 0..k) {
                    let tmp1 = alpha * b[l * ldb + i]
                    let tmp2 = alpha * a[l * lda + i]
                    let btmp = b[l * ldb + i..l * ldb + n]
                    if (tmp1 != 0.0 || tmp2 != 0.0) {
                        let atmp = a[l * lda + i..l * lda + n]
                        for (j in 0..atmp.size) {
                            ctmp[j] += atmp[j] * tmp1 + btmp[j] * tmp2
                        }
                    }
                }
            }
            return 
        }

        for (i in 0..n) {
            let ctmp = c[i * ldc..i * ldc + i + 1]
            if (beta == 0.0) {
                for (j in 0..ctmp.size) {
                    ctmp[j] = 0.0
                }
            } else if (beta == 1.0) {

            } else {
                for (j in 0..ctmp.size) {
                    ctmp[j] *= beta
                }
            }
            for (l in 0..k) {
                let tmp1 = alpha * b[l * ldb + i]
                let tmp2 = alpha * a[l * lda + i]
                let btmp = b[l * ldb..l * ldb + i + 1]
                if (tmp1 != 0.0 || tmp2 != 0.0) {
                    let atmp = a[l * lda..l * lda + i + 1]
                    for (j in 0..atmp.size) {
                        ctmp[j] += atmp[j] * tmp1 + btmp[j] * tmp2
                    }
                }
            }
        }
    }
    
    public func dtrmm(s: Side, ul: Uplo, tA: Transpose, d: Diag, m: Int64, n: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64): Unit {
        if (s != Left && s != Right) {
            throw Exception(badSide)
        }
        if (ul != Lower && ul != Upper) {
            throw Exception(badUplo)
        }
        if (tA != Trans && tA != NoTrans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (d != NonUnitDiag && d != UnitDiag) {
            throw Exception(badDiag)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }

        let k = if (s == Left) {
            m 
        } else {
            n
        }

        if (lda < max(1, k)) {
            throw Exception(badLdA)
        }
        if (ldb < max(1, n)) {
            throw Exception(badLdB)
        }

        if (m == 0 || n == 0) {
            return 
        }

        if (a.size < lda * (k - 1) + k) {
            throw Exception(shortA)
        }
        if (b.size < ldb * (m - 1) + n) {
            throw Exception(shortB)
        }

        if (alpha == 0.0) {
            for (i in 0..m) {
                let btmp = b[i * ldb..i * ldb + n]
                for (j in 0..btmp.size) {
                    btmp[j] = 0.0
                }
            }
            return 
        }

        let nonUnitDiag = d == NonUnitDiag
        if (s == Left) {
            if (tA == NoTrans) {
                if (ul == Upper) {
                    for (i in 0..m) {
                        var tmp = alpha
                        if (nonUnitDiag) {
                            tmp *= a[i * lda + i]
                        }
                        let btmp = b[i * ldb..i * ldb + n]
                        scalUnitary(tmp, btmp)
                        let atmp = a[i * lda + i + 1..i * lda + m]
                        for (ka in 0..atmp.size) {
                            let va = atmp[ka]
                            let k = ka + i + 1
                            if (va != 0.0) {
                                axpyUnitary(alpha * va, b[k * ldb..k * ldb + n], btmp)
                            }
                        }
                    }
                    return 
                }

                var i = m - 1
                while (i >= 0) {
                    var tmp = alpha
                    if (nonUnitDiag) {
                        tmp *= a[i * lda + i]
                    }
                    let btmp = b[i * ldb..i * ldb + n]
                    scalUnitary(tmp, btmp)
                    let atmp = a[i * lda..i * lda + i]
                    for (k in 0..atmp.size) {
                        let va = atmp[k]
                        if (va != 0.0) {
                            axpyUnitary(alpha * va, b[k * ldb..k * ldb + n], btmp)
                        }
                    }
                    i--
                }
                return 
            }

            if (ul == Upper) {
                var k = m - 1
                while (k >= 0) {
                    let btmpk = b[k * ldb..k * ldb + n]
                    let atmpk = a[k * lda + k + 1..k * lda + m]
                    for (ia in 0..atmpk.size) {
                        let va = atmpk[ia]
                        let i = ia + k + 1
                        let btmp = b[i * ldb..i * ldb + n]
                        if (va != 0.0) {
                            axpyUnitary(alpha * va, btmpk, btmp)
                        }
                    }
                    var tmp = alpha
                    if (nonUnitDiag) {
                        tmp *= a[k * lda + k]
                    }
                    if (tmp != 1.0) {
                        scalUnitary(tmp, btmpk)
                    }
                    k--
                }
                return 
            }

            for (k in 0..m) {
                let btmpk = b[k * ldb..k * ldb + n]
                let atmpk = a[k * lda..k * lda + k]
                for (i in 0..atmpk.size) {
                    let va = atmpk[i]
                    let btmp = b[i * ldb..i * ldb + n]
                    if (va != 0.0) {
                        axpyUnitary(alpha * va, btmpk, btmp)
                    }
                }
                var tmp = alpha
                if (nonUnitDiag) {
                    tmp *= a[k * lda + k]
                }
                if (tmp != 1.0) {
                    scalUnitary(tmp, btmpk)
                }
            }
            return 
        }

        if (tA == NoTrans) {
            if (ul == Upper) {
                for (i in 0..m) {
                    let btmp = b[i * ldb..i * ldb + n]
                    var k = n - 1
                    while (k >= 0) {
                        let tmp = alpha * btmp[k]
                        if (tmp == 0.0) {
                            k--
                            continue
                        }
                        btmp[k] = tmp
                        if (nonUnitDiag) {
                            btmp[k] *= a[k * lda + k]
                        }
                        axpyUnitary(tmp, a[k * lda + k + 1..k * lda + n], btmp[k + 1..n])
                        k--
                    }
                }
                return 
            }

            for (i in 0..m) {
                let btmp = b[i * ldb..i * ldb + n]
                for (k in 0..n) {
                    let tmp = alpha * btmp[k]
                    if (tmp == 0.0) {
                        continue
                    }
                    btmp[k] = tmp
                    if (nonUnitDiag) {
                        btmp[k] *= a[k * lda + k]
                    }
                    axpyUnitary(tmp, a[k * lda..k * lda + k], btmp[..k])
                }
            }
            return 
        }

        if (ul == Upper) {
            for (i in 0..m) {
                let btmp = b[i * ldb..i *ldb + n]
                for (j in 0..btmp.size) {
                    var tmp = btmp[j]
                    if (nonUnitDiag) {
                        tmp *= a[j * lda + j]
                    }
                    tmp += dotUnitary(a[j * lda + j + 1..j * lda + n], btmp[j + 1..n])
                    btmp[j] = alpha * tmp
                }
            }
            return 
        }

        for (i in 0..m) {
            let btmp = b[i * ldb..i * ldb + n]
            var j = n - 1
            while (j >= 0) {
                var tmp = btmp[j]
                if (nonUnitDiag) {
                    tmp *= a[j * lda + j]
                }
                tmp += dotUnitary(a[j * lda..j * lda + j], btmp[..j])
                btmp[j] = alpha * tmp
                j--
            }
        }
    }

    public func dgemm(tA: Transpose, tB: Transpose, m: Int64, n: Int64, k: Int64, alpha: Float64, a: Array<Float64>,
        lda: Int64, b: Array<Float64>, ldb: Int64, beta: Float64, c: Array<Float64>, ldc: Int64): Unit {
        if (tA != NoTrans && tA != Trans && tA != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (tB != NoTrans && tB != Trans && tB != ConjTrans) {
            throw Exception(badTranspose)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (k < 0) {
            throw Exception(kLT0)
        }

        let aTrans = tA == Trans || tA == ConjTrans
        if (aTrans) {
            if (lda < max(1, m)) {
                throw Exception(badLdA)
            }
        } else {
            if (lda < max(1, k)) {
                throw Exception(badLdA)
            }
        }

        let bTrans = tB == Trans || tB == ConjTrans
        if (bTrans) {
            if (ldb < max(1, k)) {
                throw Exception(badLdB)
            }
        } else {
            if (ldb < max(1, n)) {
                throw Exception(badLdB)
            }
        }

        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (m == 0 || n == 0) {
            return 
        }

        if (aTrans) {
            if (a.size < (k - 1) * lda + m) {
                throw Exception(shortA)
            }
        } else {
            if (a.size < (m - 1) * lda + k) {
                throw Exception(shortA)
            }
        }

        if (bTrans) {
            if (b.size < (n - 1) * ldb + k) {
                throw Exception(shortB)
            }
        } else {
            if (b.size < (k - 1) * ldb + n) {
                throw Exception(shortB)
            }
        }

        if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }

        if ((alpha == 0.0 || k == 0) && beta == 1.0) {
            return 
        }

        if (beta != 1.0) {
            if (beta == 0.0) {
                for (i in 0..m) {
                    let ctmp = c[i * ldc..i * ldc + n]
                    for (j in 0..ctmp.size) {
                        ctmp[j] = 0.0
                    }
                }
            } else {
                for (i in 0..m) {
                    let ctmp = c[i * ldc..i * ldc + n]
                    for (j in 0..ctmp.size) {
                        ctmp[j] *= beta
                    }
                }
            }
        }

        dgemmParallel(aTrans, bTrans, m, n, k, a, lda, b, ldb, c, ldc, alpha)
    }
}

func dgemmParallel(aTrans: Bool, bTrans: Bool, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64, c: Array<Float64>, ldc: Int64, alpha: Float64): Unit {
    let maxKLen = k
    let parBlocks = blocks(m, blockSize) * blocks(n, blockSize)
    if (parBlocks < minParBlock) {
		dgemmSerial(aTrans, bTrans, m, n, k, a, lda, b, ldb, c, ldc, alpha)
		return
    }

    let wg = ArrayList<Future<Unit>>()

    var i = 0
    var j = 0
    while (i < m) {
        while (j < n) {
            let i_ = i
            let j_ = j
            let fut: Future<Unit> = spawn {
                var leni = blockSize
                if (i_ + leni > m) {
                    leni = m - i_
                }

                var lenj = blockSize
                if (j_ + lenj > n) {
                    lenj = n - j_
                }

                let cSub = sliceView64(c, ldc, i_, j_, leni, lenj)

                var k = 0
                while (k < maxKLen) {
                    var lenk = blockSize
                    if (k + lenk > maxKLen) {
                        lenk = maxKLen - k
                    }
                    let aSub = if (aTrans) {
                        sliceView64(a, lda, k, i_, lenk, leni)
                    } else {
                        sliceView64(a, lda, i_, k, leni, lenk)
                    }
                    let bSub = if (bTrans) {
                        sliceView64(b, ldb, j_, k, lenj, lenk)
                    } else {
                        sliceView64(b, ldb, k, j_, lenk, lenj)
                    }
                    dgemmSerial(aTrans, bTrans, leni, lenj, lenk, aSub, lda, bSub, ldb, cSub, ldc, alpha)
                    k += blockSize
                }
            }
            
            wg.add(fut)
            j += blockSize
        }
        i += blockSize
    }

    for (fut in wg) {
        fut.get()
    }
}

func dgemmSerial(aTrans: Bool, bTrans: Bool, m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64, c: Array<Float64>, ldc: Int64, alpha: Float64): Unit {
    match {
        case !aTrans && !bTrans => dgemmSerialNotNot(m, n, k, a, lda, b, ldb, c, ldc, alpha)
        case aTrans && !bTrans  => dgemmSerialTransNot(m, n, k, a, lda, b, ldb, c, ldc, alpha)
        case !aTrans && bTrans  => dgemmSerialNotTrans(m, n, k, a, lda, b, ldb, c, ldc, alpha)
        case aTrans && bTrans   => dgemmSerialTransTrans(m, n, k, a, lda, b, ldb, c, ldc, alpha)
        case _ => throw Exception("dgemmSerial unreachable")
    }
}

func dgemmSerialNotNot(m: Int64, n: Int64, k: Int64,  a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64, c: Array<Float64>, ldc: Int64, alpha: Float64): Unit {
    for (i in 0..m) {
        let atmp = a[i * lda..i * lda + k]
        let ctmp = c[i * ldc..i * ldc + n]
        for (l in 0..atmp.size) {
            let tmp = alpha * atmp[l]
            if (tmp != 0.0) {
                axpyUnitary(tmp, b[l * ldb..l * ldb + n], ctmp)
            }
        }
    }
}

func dgemmSerialTransNot(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64, c: Array<Float64>, ldc: Int64, alpha: Float64): Unit {
    for (l in 0..k) {
        let atmp = a[l * lda..l * lda + m]
        let btmp = b[l * ldb..l * ldb + n]
        for (i in 0..atmp.size) {
            let tmp = alpha * atmp[i]
            if (tmp != 0.0) {
                let ctmp = c[i * ldc..i * ldc + n]
                axpyUnitary(tmp, btmp, ctmp)
            }
        }
    }
}

func dgemmSerialNotTrans(m: Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64, c: Array<Float64>, ldc: Int64, alpha: Float64): Unit {
    for (i in 0..m) {
        let atmp = a[i * lda..i * lda + k]
        let ctmp = c[i * ldc..i * ldc + n]
        for (j in 0..n) {
            ctmp[j] += alpha * dotUnitary(atmp, b[j * ldb..j * ldb + k])
        }
    }
}

func dgemmSerialTransTrans(m :Int64, n: Int64, k: Int64, a: Array<Float64>, lda: Int64, b: Array<Float64>, ldb: Int64, c: Array<Float64>, ldc: Int64, alpha: Float64): Unit {
    for (l in 0..k) {
        let atmp = a[l * lda..l * lda + m]
        for (i in 0..atmp.size) {
            let tmp = alpha * atmp[i]
            if (tmp != 0.0) {
                let ctmp = c[i * ldc..i * ldc + n]
                axpyInc(tmp, b[l..], ctmp, n, ldb, 1, 0, 0)
            }
        }
    }
}

func sliceView64(a: Array<Float64>, lda: Int64, i: Int64, j: Int64, r: Int64, c: Int64): Array<Float64> {
    return a[i * lda + j..(i + r - 1) * lda + j + c]
}
