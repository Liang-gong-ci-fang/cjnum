package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlarfx(side: Side, m: Int64, n: Int64, v: Array<Float64>, tau: Float64, c: Array<Float64>, ldc: Int64,
        work: Array<Float64>): Unit {
        if (side != Left && side != Right) {
            throw Exception(badSide)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (ldc < max(1, n)) {
            throw Exception(badLdC)
        }

        if (m == 0 || n == 0) {
            return
        }

        var nh = m
        var lwork = n
        if (side == Right) {
            nh = n
            lwork = m
        }
        if (v.size < nh) {
            throw Exception(shortV)
        }
        if (c.size < (m - 1) * ldc + n) {
            throw Exception(shortC)
        }
        if (nh > 10 && work.size < lwork) {
            throw Exception(shortWork)
        }

        if (tau == 0.0) {
            return
        }

        if (side == Left) {
            match (m) {
                case 0 => return

                case 1 =>
                    let t0 = 1.0 - tau * v[0] * v[0]
                    for (j in 0..n) {
                        c[j] *= t0
                    }
                    return

                case 2 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                    }
                    return

                case 3 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                    }
                    return

                case 4 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    let v3 = v[3]
                    let t3 = tau * v3
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j] + v3 * c[3 * ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                        c[3 * ldc + j] -= sum * t3
                    }
                    return

                case 5 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    let v3 = v[3]
                    let t3 = tau * v3
                    let v4 = v[4]
                    let t4 = tau * v4
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j] + v3 * c[3 * ldc + j] + v4 * c[4 *
                                ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                        c[3 * ldc + j] -= sum * t3
                        c[4 * ldc + j] -= sum * t4
                    }
                    return

                case 6 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    let v3 = v[3]
                    let t3 = tau * v3
                    let v4 = v[4]
                    let t4 = tau * v4
                    let v5 = v[5]
                    let t5 = tau * v5
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j] + v3 * c[3 * ldc + j] + v4 * c[4 *
                                ldc + j] + v5 * c[5 * ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                        c[3 * ldc + j] -= sum * t3
                        c[4 * ldc + j] -= sum * t4
                        c[5 * ldc + j] -= sum * t5
                    }
                    return

                case 7 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    let v3 = v[3]
                    let t3 = tau * v3
                    let v4 = v[4]
                    let t4 = tau * v4
                    let v5 = v[5]
                    let t5 = tau * v5
                    let v6 = v[6]
                    let t6 = tau * v6
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j] + v3 * c[3 * ldc + j] + v4 * c[4 *
                                ldc + j] + v5 * c[5 * ldc + j] + v6 * c[6 * ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                        c[3 * ldc + j] -= sum * t3
                        c[4 * ldc + j] -= sum * t4
                        c[5 * ldc + j] -= sum * t5
                        c[6 * ldc + j] -= sum * t6
                    }
                    return

                case 8 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    let v3 = v[3]
                    let t3 = tau * v3
                    let v4 = v[4]
                    let t4 = tau * v4
                    let v5 = v[5]
                    let t5 = tau * v5
                    let v6 = v[6]
                    let t6 = tau * v6
                    let v7 = v[7]
                    let t7 = tau * v7
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j] + v3 * c[3 * ldc + j] + v4 * c[4 *
                                ldc + j] + v5 * c[5 * ldc + j] + v6 * c[6 * ldc + j] + v7 * c[7 * ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                        c[3 * ldc + j] -= sum * t3
                        c[4 * ldc + j] -= sum * t4
                        c[5 * ldc + j] -= sum * t5
                        c[6 * ldc + j] -= sum * t6
                        c[7 * ldc + j] -= sum * t7
                    }
                    return

                case 9 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    let v3 = v[3]
                    let t3 = tau * v3
                    let v4 = v[4]
                    let t4 = tau * v4
                    let v5 = v[5]
                    let t5 = tau * v5
                    let v6 = v[6]
                    let t6 = tau * v6
                    let v7 = v[7]
                    let t7 = tau * v7
                    let v8 = v[8]
                    let t8 = tau * v8
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j] + v3 * c[3 * ldc + j] + v4 * c[4 *
                                ldc + j] + v5 * c[5 * ldc + j] + v6 * c[6 * ldc + j] + v7 * c[7 * ldc + j] + v8 * c[8 *
                                ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                        c[3 * ldc + j] -= sum * t3
                        c[4 * ldc + j] -= sum * t4
                        c[5 * ldc + j] -= sum * t5
                        c[6 * ldc + j] -= sum * t6
                        c[7 * ldc + j] -= sum * t7
                        c[8 * ldc + j] -= sum * t8
                    }
                    return

                case 10 =>
                    let v0 = v[0]
                    let t0 = tau * v0
                    let v1 = v[1]
                    let t1 = tau * v1
                    let v2 = v[2]
                    let t2 = tau * v2
                    let v3 = v[3]
                    let t3 = tau * v3
                    let v4 = v[4]
                    let t4 = tau * v4
                    let v5 = v[5]
                    let t5 = tau * v5
                    let v6 = v[6]
                    let t6 = tau * v6
                    let v7 = v[7]
                    let t7 = tau * v7
                    let v8 = v[8]
                    let t8 = tau * v8
                    let v9 = v[9]
                    let t9 = tau * v9
                    for (j in 0..n) {
                        let sum = v0 * c[j] + v1 * c[ldc + j] + v2 * c[2 * ldc + j] + v3 * c[3 * ldc + j] + v4 * c[4 *
                                ldc + j] + v5 * c[5 * ldc + j] + v6 * c[6 * ldc + j] + v7 * c[7 * ldc + j] + v8 * c[8 *
                                ldc + j] + v9 * c[9 * ldc + j]
                        c[j] -= sum * t0
                        c[ldc + j] -= sum * t1
                        c[2 * ldc + j] -= sum * t2
                        c[3 * ldc + j] -= sum * t3
                        c[4 * ldc + j] -= sum * t4
                        c[5 * ldc + j] -= sum * t5
                        c[6 * ldc + j] -= sum * t6
                        c[7 * ldc + j] -= sum * t7
                        c[8 * ldc + j] -= sum * t8
                        c[9 * ldc + j] -= sum * t9
                    }
                    return

                case _ =>
                    dlarf(side, m, n, v, 1, tau, c, ldc, work)
                    return
            }
        }

        match (n) {
            case 0 => return

            case 1 =>
                let t0 = 1.0 - tau * v[0] * v[0]
                for (j in 0..m) {
                    c[j * ldc] *= t0
                }
                return

            case 2 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                }
                return

            case 3 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                }
                return

            case 4 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                let v3 = v[3]
                let t3 = tau * v3
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2] + v3 * cs[3]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                    cs[3] -= sum * t3
                }
                return

            case 5 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                let v3 = v[3]
                let t3 = tau * v3
                let v4 = v[4]
                let t4 = tau * v4
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2] + v3 * cs[3] + v4 * cs[4]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                    cs[3] -= sum * t3
                    cs[4] -= sum * t4
                }
                return

            case 6 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                let v3 = v[3]
                let t3 = tau * v3
                let v4 = v[4]
                let t4 = tau * v4
                let v5 = v[5]
                let t5 = tau * v5
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2] + v3 * cs[3] + v4 * cs[4] + v5 * cs[5]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                    cs[3] -= sum * t3
                    cs[4] -= sum * t4
                    cs[5] -= sum * t5
                }
                return

            case 7 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                let v3 = v[3]
                let t3 = tau * v3
                let v4 = v[4]
                let t4 = tau * v4
                let v5 = v[5]
                let t5 = tau * v5
                let v6 = v[6]
                let t6 = tau * v6
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2] + v3 * cs[3] + v4 * cs[4] + v5 * cs[5] + v6 * cs[6]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                    cs[3] -= sum * t3
                    cs[4] -= sum * t4
                    cs[5] -= sum * t5
                    cs[6] -= sum * t6
                }
                return

            case 8 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                let v3 = v[3]
                let t3 = tau * v3
                let v4 = v[4]
                let t4 = tau * v4
                let v5 = v[5]
                let t5 = tau * v5
                let v6 = v[6]
                let t6 = tau * v6
                let v7 = v[7]
                let t7 = tau * v7
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2] + v3 * cs[3] + v4 * cs[4] + v5 * cs[5] + v6 * cs[6] +
                        v7 * cs[7]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                    cs[3] -= sum * t3
                    cs[4] -= sum * t4
                    cs[5] -= sum * t5
                    cs[6] -= sum * t6
                    cs[7] -= sum * t7
                }
                return

            case 9 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                let v3 = v[3]
                let t3 = tau * v3
                let v4 = v[4]
                let t4 = tau * v4
                let v5 = v[5]
                let t5 = tau * v5
                let v6 = v[6]
                let t6 = tau * v6
                let v7 = v[7]
                let t7 = tau * v7
                let v8 = v[8]
                let t8 = tau * v8
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2] + v3 * cs[3] + v4 * cs[4] + v5 * cs[5] + v6 * cs[6] +
                        v7 * cs[7] + v8 * cs[8]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                    cs[3] -= sum * t3
                    cs[4] -= sum * t4
                    cs[5] -= sum * t5
                    cs[6] -= sum * t6
                    cs[7] -= sum * t7
                    cs[8] -= sum * t8
                }
                return

            case 10 =>
                let v0 = v[0]
                let t0 = tau * v0
                let v1 = v[1]
                let t1 = tau * v1
                let v2 = v[2]
                let t2 = tau * v2
                let v3 = v[3]
                let t3 = tau * v3
                let v4 = v[4]
                let t4 = tau * v4
                let v5 = v[5]
                let t5 = tau * v5
                let v6 = v[6]
                let t6 = tau * v6
                let v7 = v[7]
                let t7 = tau * v7
                let v8 = v[8]
                let t8 = tau * v8
                let v9 = v[9]
                let t9 = tau * v9
                for (j in 0..m) {
                    let cs = c[j * ldc..]
                    let sum = v0 * cs[0] + v1 * cs[1] + v2 * cs[2] + v3 * cs[3] + v4 * cs[4] + v5 * cs[5] + v6 * cs[6] +
                        v7 * cs[7] + v8 * cs[8] + v9 * cs[9]
                    cs[0] -= sum * t0
                    cs[1] -= sum * t1
                    cs[2] -= sum * t2
                    cs[3] -= sum * t3
                    cs[4] -= sum * t4
                    cs[5] -= sum * t5
                    cs[6] -= sum * t6
                    cs[7] -= sum * t7
                    cs[8] -= sum * t8
                    cs[9] -= sum * t9
                }
                return

            case _ =>
                dlarf(side, m, n, v, 1, tau, c, ldc, work)
                return
        }
    }
}
