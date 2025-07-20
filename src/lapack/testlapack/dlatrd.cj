package cjnum.lapack.testlapack

import cjnum.blas.blas64.General as General64
import cjnum.blas.blas64.*
import cjnum.blas.*
import cjnum.lapack.*
import std.math.*
import std.random.*

class dlatrdTestStruct {
    dlatrdTestStruct(
        let n: Int64,
        let nb: Int64,
        let lda: Int64,
        let ldw: Int64
    ) {}
}

public func dlatrdTest(lFloat64: LFloat64): Bool {
    const tol = 1e-14

    let rnd = Random(1)
    for (uplo in [Upper, Lower]) {
        for (test in [
                dlatrdTestStruct(5, 2, 0, 0),
                dlatrdTestStruct(5, 5, 0, 0),
                
                dlatrdTestStruct(5, 3, 10, 11),
                dlatrdTestStruct(5, 5, 10, 11)
            ]) {
            let n = test.n
            let nb = test.nb
            let lda = if (test.lda == 0) {
                n
            } else {
                test.lda
            }
            let ldw = if (test.ldw == 0) {
                nb
            } else {
                test.ldw
            }

            let a = Array<Float64>(n * lda, repeat: 0.0)
            for (i in 0..a.size) {
                a[i] = rnd.nextGaussianFloat64()
            }

            let e = Array<Float64>(n - 1, repeat: 0.0)
            for (i in 0..e.size) {
                e[i] = Float64.NaN
            }
            let tau = Array<Float64>(n - 1, repeat: 0.0)
            for (i in 0..tau.size) {
                tau[i] = Float64.NaN
            }
            let w = Array<Float64>(n * ldw, repeat: 0.0)
            for (i in 0..w.size) {
                w[i] = Float64.NaN
            }

            let aCopy = Array<Float64>(a.size, repeat: 0.0)
            copy(aCopy, a)

            lFloat64.dlatrd(uplo, n, nb, a, lda, e, tau, w, ldw)

            let q = General64(
                rows: n, 
                cols: n, 
                stride: n,
                data: Array<Float64>(n * n, repeat: 0.0)
            )

            for (i in 0..n) {
                q.data[i * q.stride + i] = 1.0
            }
            if (uplo == Upper) {
                for (i in n - 1..=n - nb : -1) {
                    if (i == 0) {
                        continue
                    }
                    let v = Vector(
                        inc: 1, 
                        data: Array<Float64>(n, repeat: 0.0)
                    )
                    for (j in 0..i - 1) {
                        v.data[j] = a[j * lda + i]
                    }
                    v.data[i - 1] = 1.0

                    let h = General64(
                        rows: n, cols: n, stride: n, data: Array<Float64>(n * n, repeat: 0.0)
                    )
                    for (j in 0..n) {
                        h.data[j * n + j] = 1.0
                    }
                    ger(-tau[i - 1], v, v, h)

                    let qTmp = General64(
                        rows: n, cols: n, stride: n, data: Array<Float64>(n * n, repeat: 0.0)
                    )
                    copy(qTmp.data, q.data)
                    gemm(NoTrans, NoTrans, 1.0, qTmp, h, 0.0, q)
                }
            } else {
                for (i in 0..nb) {
                    if (i == n - 1) {
                        continue
                    }

                    let v = Vector(
                        inc: 1, 
                        data: Array<Float64>(n, repeat: 0.0)
                    )
                    v.data[i + 1] = 1.0
                    for (j in i + 2..n) {
                        v.data[j] = a[j * lda + i]
                    }

                    let h = General64(
                        rows: n, cols: n, stride: n, data: Array<Float64>(n * n, repeat: 0.0)
                    )
                    for (j in 0..n) {
                        h.data[j * n + j] = 1.0
                    }
                    ger(-tau[i], v, v, h)

                    let qTmp = General64(
                        rows: n, cols: n, stride: n, data: Array<Float64>(n * n, repeat: 0.0)
                    )
                    for (k in 0..q.data.size) {
                        qTmp.data[k] = q.data[k]
                    }
                    gemm(NoTrans, NoTrans, 1.0, qTmp, h, 0.0, q)
                }
            }

            let name = "uplo=${uplo},n=${n},nb=${nb}"
            if (residualOrthogonal(q, false) > tol) {
                throw Exception("Case ${name}: Q is not orthogonal")
            }

            let aGen = genFromSym(Symmetric(n: n, stride: lda, uplo: uplo, data: aCopy))
            if (!dlatrdCheckDecomposition(uplo, n, nb, e, a, lda, aGen, q, tol)) {
                throw Exception("Case ${name}: Decomposition mismatch")
            }
        }
    }
    return true
}

func dlatrdCheckDecomposition(uplo: Uplo, n: Int64, nb: Int64, e: Array<Float64>, a: Array<Float64>,
    lda: Int64, aGen: General64, q: General64, tol: Float64): Bool {
    let tmp = General64(
        rows: n,
        cols: n,
        stride: n,
        data: Array<Float64>(n * n, repeat: 0.0)
    )
    let ans = General64(
        rows: n,
        cols: n,
        stride: n,
        data: Array<Float64>(n * n, repeat: 0.0)
    )
    gemm(Trans, NoTrans, 1.0, q, aGen, 0.0, tmp)
    gemm(NoTrans, NoTrans, 1.0, tmp, q, 0.0, ans)

    if (uplo == Upper) {
        for (i in n - nb..n) {
            for (j in 0..n) {
                let v = ans.data[i * ans.stride + j]
                if (i == j) {
                    if (abs(v - a[i * lda + j]) > tol) {
                        return false
                    }
                } else if (i == j - 1) {
                    if (abs(a[i * lda + j] - 1.0) > tol) {
                        return false
                    }
                    if (abs(v - e[i]) > tol) {
                        return false
                    }
                } else if (i == j + 1) {
                } else {
                    if (abs(v) > tol) {
                        return false
                    }
                }
            }
        }
    } else {
        for (i in 0..nb) {
            for (j in 0..n) {
                let v = ans.data[i * ans.stride + j]
                if (i == j) {
                    if (abs(v - a[i * lda + j]) > tol) {
                        return false
                    }
                } else if (i == j - 1) {
                } else if (i == j + 1) {
                    if (abs(a[i * lda + j] - 1.0) > tol) {
                        return false
                    }
                    if (abs(v - e[i - 1]) > tol) {
                        return false
                    }
                } else {
                    if (abs(v) > tol) {
                        return false
                    }
                }
            }
        }
    }
    return true
}

func genFromSym(a: Symmetric): General64 {
    let n = a.n
    let lda = a.stride
    let uplo = a.uplo
    let b = General64(
        rows: n,
        cols: n,
        stride: n,
        data: Array<Float64>(n * n, repeat: 0.0)
    )

    for (i in 0..n) {
        for (j in i..n) {
            var v = a.data[i * lda + j]
            if (uplo == Lower) {
                v = a.data[j * lda + i]
            }
            b.data[i * n + j] = v
            b.data[j * n + i] = v
        }
    }
    return b
}
