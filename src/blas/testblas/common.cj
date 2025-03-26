package cjnum.blas.testblas

import std.math.*
import cjnum.complex.*
import cjnum.blas.*
import std.collection.*

const throwException = true

let znan = Complex128(Float64.NaN, Float64.NaN)

func dTolEqual(a: Float64, b: Float64): Bool {
    if (a.isNaN() && b.isNaN()) {
        return true
    }
    if (a == b) {
        return true
    }
    var a_ = a
    var b_ = b
    let m = max(abs(a), abs(b))
    if (m > 1.0) {
        a_ = a / m
        b_ = b / m
    }
    if (abs(a_ - b_) < 1e-14) {
        return true
    }
    return false
}

func dSliceTolEqual(a: Array<Float64>, b: Array<Float64>): Bool {
    if (a.size != b.size) {
        return false
    }
    for (i in 0..a.size) {
        if (!dTolEqual(a[i], b[i])) {
            return false
        }
    }
    return true
}

func dStridedSliceTolEqual(n: Int64, a: Array<Float64>, inca: Int64, b: Array<Float64>, incb: Int64): Bool {
    var ia = 0
    var ib = 0
    if (inca <= 0) {
        ia = -(n - 1) * inca
    }
    if (incb <= 0) {
        ib = -(n - 1) * incb
    }
    for (i in 0..n) {
        if (!dTolEqual(a[ia], b[ib])) {
            return false
        }
        ia += inca
        ib += incb
    }
    return true
}

func dSliceEqual(a: Array<Float64>, b: Array<Float64>): Bool {
    if (a.size != b.size) {
        return false
    }
    for (i in 0..a.size) {
        if (!dTolEqual(a[i], b[i])) {
            return false
        }
    }
    return true
}

func dCopyTwoTmp(x: Array<Float64>, xTmp: Array<Float64>, y: Array<Float64>, yTmp: Array<Float64>): Unit {
    if (x.size != xTmp.size) {
        throw Exception("x size mismatch")
    }
    if (y.size != yTmp.size) {
        throw Exception("y size mismatch")
    }
    copy(xTmp, x)
    copy(yTmp, y)
}

func copy<T>(dst: Array<T>, src: Array<T>): Int64 {
    let minSize = min(src.size, dst.size)
    dst[0..minSize] = src[0..minSize]
    return minSize
}

func exception(f: () -> Unit): Bool {
    try {
        f()
    } catch (e: Exception) {
        return true
    }
    return false
}

func testpanics(f: () -> Unit, name: String): Bool {
    if (!exception(f)) {
        throw Exception("${name} should panic and does not")
    }
    return true
}

func sliceOfSliceCopy(a: Array<Array<Float64>>): Array<Array<Float64>> {
    return Array<Array<Float64>>(a.size, { i => 
        Array<Float64>(a[i].size, { j => 
            a[i][j]
        })
    })
}

func sliceCopy(a: Array<Float64>): Array<Float64> {
    return Array<Float64>(a.size, { i => a[i] })
}

func flatten(a: Array<Array<Float64>>): ?Array<Float64> {
    if (a.size == 0) {
        return None
    }
    let m = a.size
    let n = a[0].size
    let s = Array<Float64>(m * n, item: 0.0)
    for (i in 0..m) {
        for (j in 0..n) {
            s[i * n + j] = a[i][j]
        }
    }
    return s
}

func unflatten(a: Array<Float64>, m: Int64, n: Int64): Array<Array<Float64>> {
    return Array<Array<Float64>>(m, { i => 
        Array<Float64>(n, { j => 
            a[i * n + j]
        })
    })
}

func flattenTriangular(a: Array<Array<Float64>>, ul: Uplo): Array<Float64> {
    let m = a.size
    let aFlat = Array<Float64>(m * (m + 1) / 2, item: 0.0)
    var k: Int64 = 0
    if (ul == Upper) {
        for (i in 0..m) {
            k += copy(aFlat[k..], a[i][i..])
        }
        return aFlat
    }
    for (i in 0..m) {
        k += copy(aFlat[k..], a[i][..i + 1])
    }
    return aFlat
}

func flattenBanded(a: Array<Array<Float64>>, ku: Int64, kl: Int64): Array<Float64> {
    let m = a.size
    let n = a[0].size
    if (ku < 0 || kl < 0) {
        throw Exception("testblas: negative band length")
    }
    let nRows = m
    let nCols = (ku + kl + 1)
    let aflat = Array<Float64>(nRows * nCols, item: Float64.NaN)

    for (i in 0..nRows) {
        var min_ = -kl
        if (i - kl < 0) {
            min_ = -i
        }
        var max_ = ku
        if (i + ku >= n) {
            max_ = n - i - 1
        }
        for (j in min_..=max_) {
            let col = kl + j
            aflat[i * nCols + col] = a[i][i + j]
        }
    }
    return aflat
}

func makeIncremented(x: Array<Float64>, inc: Int64, extra: Int64): Array<Float64> {
    if (inc == 0) {
        throw Exception("zero inc")
    }
    let absinc = abs(inc)
    let xcopy = Array<Float64>(x.size, item: 0.0)
    if (inc > 0) {
        copy(xcopy, x)
    } else {
        for (i in 0..x.size) {
            xcopy[i] = x[x.size - i - 1]
        }
    }

    var counter: Float64 = 100.0
    let xnew = ArrayList<Float64>()
    for (i in 0..xcopy.size) {
        xnew.append(xcopy[i])
        if (i != x.size - 1) {
            for (j in 0..absinc - 1) {
                xnew.append(counter)
                counter += 1.0
            }
        }
    }
    for (i in 0..extra) {
        xnew.append(counter)
        counter += 1.0
    }
    return xnew.toArray()
}

func makeIncremented32(x: Array<Float32>, inc: Int64, extra: Int64): Array<Float32> {
    if (inc == 0) {
        throw Exception("zero inc")
    }
    let absinc = abs(inc)
    let xcopy = Array<Float32>(x.size, item: 0.0)
    if (inc > 0) {
        copy(xcopy, x)
    } else {
        for (i in 0..x.size) {
            xcopy[i] = x[x.size - i - 1]
        }
    }

    var counter: Float32 = 100.0
    let xnew = ArrayList<Float32>()
    for (i in 0..xcopy.size) {
        xnew.append(xcopy[i])
        if (i != x.size - 1) {
            for (j in 0..absinc - 1) {
                xnew.append(counter)
                counter += 1.0
            }
        }
    }
    for (i in 0..extra) {
        xnew.append(counter)
        counter += 1.0
    }
    return xnew.toArray()
}

func allPairs(x: Array<Int64>, y: Array<Int64>): Array<Array<Int64>> {
    let p = ArrayList<Array<Int64>>()
    for (v0 in x) {
        for (v1 in y) {
            p.append([v0, v1])
        }
    }
    return p.toArray()
}

func sameFloat64(a: Float64, b: Float64): Bool {
    return a == b || a.isNaN() && b.isNaN()
}

// TODO Complex

func transString(t: Transpose): String {
    return if (t == NoTrans) {
        "NoTrans"
    } else if (t == Trans) {
        "Trans"
    } else if (t == ConjTrans) {
        "ConjTrans"
    } else {
        "unknown trans"
    }
}

func uploString(uplo: Uplo): String {
    return if (uplo == Lower) {
        "Lower"
    } else if (uplo == Upper) {
        "Upper"
    } else {
        "unknown uplo"
    }
}

func sideString(side: Side): String {
    return if (side == Left) {
        "Left"
    } else if (side == Right) {
        "Right"
    } else {
        "unknown side"
    }
}

func diagString(diag: Diag): String {
    return if (diag == UnitDiag) {
        "UnitDiag"
    } else if (diag == NonUnitDiag) {
        "NonUnitDiag"
    } else {
        "unknown diag"
    }
}

// TODO complex