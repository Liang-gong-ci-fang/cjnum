package cjnum.lapack.testlapack

import std.random.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.lapack.*
import cjnum.blas.blas64.*
import cjnum.blas.blas64.General as General64

class DgetrsTestStruct {
    DgetrsTestStruct(
        let n: Int64,
        let nrhs: Int64,
        let lda: Int64,
        let ldb: Int64,
        let tol: Float64
    ) {}
}

public func dgetrsTest(lFloat64: LFloat64): Bool {
    let rnd = Random(1)
    for (trans in [NoTrans, Trans]) {
        for (test in [
            DgetrsTestStruct(3, 3, 0, 0, 1e-12),
            DgetrsTestStruct(3, 5, 0, 0, 1e-12),
            DgetrsTestStruct(5, 3, 0, 0, 1e-12),

            DgetrsTestStruct(3, 3, 8, 10, 1e-12),
            DgetrsTestStruct(3, 5, 8, 10, 1e-12),
            DgetrsTestStruct(5, 3, 8, 10, 1e-12),

            DgetrsTestStruct(300, 300, 0, 0, 1e-8),
            DgetrsTestStruct(300, 500, 0, 0, 1e-8),
            DgetrsTestStruct(500, 300, 0, 0, 1e-6),

            DgetrsTestStruct(300, 300, 700, 600, 1e-8),
            DgetrsTestStruct(300, 500, 700, 600, 1e-8),
            DgetrsTestStruct(500, 300, 700, 600, 1e-6)
        ]) {
            let n = test.n
            let nrhs = test.nrhs
            let lda = if (test.lda == 0) {
                n
            } else {
                test.lda
            }
            let ldb = if (test.ldb == 0) {
                nrhs
            } else {
                test.ldb
            }
            let a = Array<Float64>(n * lda, { _ => rnd.nextFloat64() })
            let b = Array<Float64>(n * ldb, { _ => rnd.nextFloat64() })
            let aCopy = a.clone()
            let bCopy = b.clone()

            let ipiv = Array<Int64>(n, { _=> rnd.nextInt64() })

            lFloat64.dgetrf(n, n, a, lda, ipiv)
            lFloat64.dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb)

            let A = General64(
                rows: n,
                cols: n,
                stride: lda,
                data: aCopy
            )
            let B =General64(
                rows: n,
                cols: nrhs,
                stride: ldb,
                data: bCopy
            )
            let X = General64(
                rows: n,
                cols: nrhs,
                stride: ldb,
                data: b
            )
            let tmp = General64(
                rows: n,
                cols: nrhs,
                stride: ldb,
                data: Array<Float64>(n * ldb, repeat: 0.0)
            )
            copy(tmp.data, bCopy)
            gemm(trans, NoTrans, 1.0, A, X, 0.0, B)
            if (!fEqualApprox(tmp.data, bCopy, test.tol)) {
                throw Exception("Linear solve mismatch. trans = ${trans}, n = ${n}, nrhs = ${nrhs}, lda = ${lda}, ldb = ${ldb}")
            }
        }
    }
    return true
}