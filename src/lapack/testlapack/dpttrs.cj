package cjnum.lapack.testlapack

import std.random.*
import cjnum.lapack.*
import cjnum.blas.blas64.General as General64

public func dpttrsTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (n in [0, 1, 2, 3, 4, 5, 10, 20, 50, 51, 52, 53, 54, 100]) {
        for (nrhs in [0, 1, 2, 3, 4, 5, 10, 20, 50]) {
            for (ldb in [max(1, nrhs), nrhs + 3]) {
                _dpttrsTest(lFloat64, rnd, n, nrhs, ldb)
            }
        }
    }
    return true
}

// func _dpttrsTest(lFloat64: LFloat64, rnd: Random, n: Int64, nrhs: Int64, ldb: Int64): Unit {
//     const tol = 1e-15

//     let name = "n=${n}, nrhs=${nrhs}, ldb = ${ldb}"
//     // let (d, e) = newRandomSymTridiag(n, rnd)

//     let itd = StringReader(BufferedInputStream(
//         ByteBuffer(File.readFrom("/home/yanami/code/tmp/d.txt")))).lines()
	
// 	let d = Array<Float64>(n, repeat: 0.0)

//     var i = 0
//     while (let Some(num) <- itd.next()) {
//         d[i] = Decimal.parse(num).toFloat64()
//         i++
//     }

//     let ite = StringReader(BufferedInputStream(
//         ByteBuffer(File.readFrom("/home/yanami/code/tmp/e.txt")))).lines()
	
// 	let e = Array<Float64>(n - 1, repeat: 0.0)

//     var j = 0
//     while (let Some(num) <- ite.next()) {
//         if (num == "NaN") {
//             e[j] = Float64.NaN
//         } else {
//             e[j] = Decimal.parse(num).toFloat64()
//         }
//         j++
//     }

//     let dFac = d.clone()
//     let eFac = e.clone()

//     let ok = lFloat64.dpttrf(n, dFac, eFac)
//     if (!ok) {
//         throw Exception("${name}: bad test matrix, Dpttrf failed")
//     }

//     let xWant = randomGeneral(n, nrhs, ldb, rnd)

//     let itx = StringReader(BufferedInputStream(
//         ByteBuffer(File.readFrom("/home/yanami/code/tmp/x.txt")))).lines()
	
// 	let x = Array<Float64>(23, repeat: 0.0)

//     var k = 0
//     while (let Some(num) <- itx.next()) {
//         if (num == "NaN") {
//             x[k] = Float64.NaN
//         } else {
//             x[k] = Decimal.parse(num).toFloat64()
//         }
//         k++
//     }

//     xWant.data = x

//     let b = zeros(n, nrhs, ldb)
//     dstmm(n, nrhs, d, e, xWant.data, xWant.stride, b.data, b.stride)

//     lFloat64.dpttrs(n, nrhs, dFac, eFac, b.data, b.stride)

//     let resid = dpttrsResidual(b, xWant)
//     if (resid > tol) {
//         throw Exception("${name}: unexpected solution: |diff| = ${resid}, want <= ${tol}")
//     }
// }


func _dpttrsTest(lFloat64: LFloat64, rnd: Random, n: Int64, nrhs: Int64, ldb: Int64): Unit {
    const tol = 1e-15

    let name = "n=${n}, nrhs=${nrhs}, ldb = ${ldb}"
    let (d, e) = newRandomSymTridiag(n, rnd)
    let dFac = d.clone()
    let eFac = e.clone()

    let ok = lFloat64.dpttrf(n, dFac, eFac)
    if (!ok) {
        throw Exception("${name}: bad test matrix, Dpttrf failed")
    }

    let xWant = randomGeneral(n, nrhs, ldb, rnd)

    let b = zeros(n, nrhs, ldb)
    dstmm(n, nrhs, d, e, xWant.data, xWant.stride, b.data, b.stride)

    lFloat64.dpttrs(n, nrhs, dFac, eFac, b.data, b.stride)

    let resid = dpttrsResidual(b, xWant)
    if (resid > tol) {
        throw Exception("${name}: unexpected solution: |diff| = ${resid}, want <= ${tol}")
    }
}

func dstmm(m: Int64, n: Int64, d: Array<Float64>, e: Array<Float64>, b: Array<Float64>, ldb: Int64, c: Array<Float64>, ldc: Int64): Unit {
    if (m == 0 || n == 0) {
        return 
    }
    if (m == 1) {
        let d0 = d[0]
        for (j in 0..n) {
            c[j] = d0 * b[j]
        }
        return 
    }
    for (j in 0..n) {
        c[j] = d[0] * b[j] + e[0] * b[ldb + j]
    }
    for (i in 1..m-1) {
        for (j in 0..n) {
            c[i * ldc + j] = e[i - 1] * b[(i - 1) * ldb + j] + d[i] * b[i * ldb + j] + e[i] * b[(i + 1) * ldb + j]
        }
    }
    for (j in 0..n) {
        c[(m - 1) * ldc + j] = e[m - 2] * b[(m - 2) * ldb + j] + d[m - 1] * b[(m - 1) * ldb + j]
    }
}

func dpttrsResidual(xGot: General64, xWant: General64): Float64 {
    let (n, nrhs) = (xGot.rows, xGot.cols)
    let d = zeros(n, nrhs, nrhs)
    for (i in 0..n) {
        for (j in 0..nrhs) {
            d.data[i * d.stride + j] = xGot.data[i * xGot.stride + j] - xWant.data[i * xWant.stride + j]
            // let a = xGot.data[i * xGot.stride + j]
            // let b = xWant.data[i * xWant.stride + j]
            // var comp = 0.0 // 补偿项
            // let diff = a - b
            // // 补偿计算（减少抵消误差）
            // let tmp = diff - comp
            // let sum = d.data[i * d.stride + j] + tmp
            // comp = (sum - d.data[i * d.stride + j]) - tmp
            // d.data[i * d.stride + j] = sum
        }
    }
    return (dlange(MaxColumnSum, n, nrhs, d.data, d.stride)) / Float64(n)
}