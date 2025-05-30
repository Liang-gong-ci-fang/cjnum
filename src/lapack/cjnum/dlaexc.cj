package cjnum.lapack.cjnum

import cjnum.blas.*
import cjnum.blas.blas64.*
import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dlaexc(wantq: Bool, n: Int64, t: Array<Float64>, ldt: Int64, q: Array<Float64>, ldq: Int64, j1: Int64, n1: Int64,
        n2: Int64, work: Array<Float64>): Bool {
        if (n < 0) {
            throw Exception(nLT0)
        } else if (ldt < max(1, n)) {
            throw Exception(badLdT)
        } else if (wantq && ldt < max(1, n)) {
            throw Exception(badLdQ)
        } else if (j1 < 0 || n <= j1) {
            throw Exception(badJ1)
        } else if (work.size < n) {
            throw Exception(shortWork)
        } else if (n1 < 0 || n1 > 2) {
            throw Exception(badN1)
        } else if (n2 < 0 || n2 > 2) {
            throw Exception(badN2)
        }

        if (n == 0 || n1 == 0 || n2 == 0) {
            return true
        }

        if (t.size < (n - 1) * ldt + n) {
            throw Exception(shortT)
        } else if (wantq && q.size < (n - 1) * ldq + n) {
            throw Exception(shortQ)
        }

        if (j1 + n1 >= n) {
            return true
        }

        let j2 = j1 + 1
        var j3 = j1 + 2

        let nFloat64 = nFloat64Implementation()

        if (n1 == 1 && n2 == 1) {
            let t11 = t[j1 * ldt + j1]
            let t22 = t[j2 * ldt + j2]

            let (cs, sn, _) = dlartg(t[j1 * ldt + j2], t22 - t11)

            if (n - j3 > 0) {
                nFloat64.drot(n - j3, t[j1 * ldt + j3..], 1, t[j2 * ldt + j3..], 1, cs, sn)
            }
            if (j1 > 0) {
                nFloat64.drot(j1, t[j1..], ldt, t[j2..], ldt, cs, sn)
            }

            t[j1 * ldt + j1] = t22
            t[j2 * ldt + j2] = t11

            if (wantq) {
                nFloat64.drot(n, q[j1..], ldq, q[j2..], ldq, cs, sn)
            }

            return true
        }

        let nd = n1 + n2
        var d = Array<Float64>(16, item: 0.0)
        const ldd = 4
        dlacpy(All, nd, nd, t[j1 * ldt + j1..], ldt, d[..], ldd)
        let dnorm = dlange(MaxAbs, nd, nd, d[..], ldd, work)

        let eps = dlamchP
        let thresh = max(10.0 * eps * dnorm, dlamchS / eps)

        var x = Array<Float64>(4, item: 0.0)
        const ldx = 2
        let (scale, _, _) = dlasy2(false, false, -1, n1, n2, d[..], ldd, d[n1 * ldd + n1..], ldd, d[n1..], ldd, x[..],
            ldx)

        match {
            case n1 == 1 && n2 == 2 =>
                let u = Array<Float64>([scale, x[0], 1.0])
                let (_, tau) = dlarfg(3, x[1], u[..2], 1)
                let t11 = t[j1 * ldt + j1]

                dlarfx(Left, 3, 3, u[..], tau, d[..], ldd, work)
                dlarfx(Right, 3, 3, u[..], tau, d[..], ldd, work)

                if (max(abs(d[2 * ldd]), max(abs(d[2 * ldd + 1]), abs(d[2 * ldd + 2] - t11))) > thresh) {
                    return false
                }

                dlarfx(Left, 3, n - j1, u[..], tau, t[j1 * ldt + j1..], ldt, work)
                dlarfx(Right, j2 + 1, 3, u[..], tau, t[j1..], ldt, work)

                t[j3 * ldt + j1] = 0.0
                t[j3 * ldt + j2] = 0.0
                t[j3 * ldt + j3] = t11

                if (wantq) {
                    dlarfx(Right, n, 3, u[..], tau, q[j1..], ldq, work)
                }

            case n1 == 2 && n2 == 1 =>
                let u = Array<Float64>([1.0, -x[ldx], scale])
                let (_, tau) = dlarfg(3, -x[0], u[1..], 1)
                let t33 = t[j3 * ldt + j3]

                dlarfx(Left, 3, 3, u[..], tau, d[..], ldd, work)
                dlarfx(Right, 3, 3, u[..], tau, d[..], ldd, work)

                if (max(abs(d[ldd]), max(abs(d[2 * ldd]), abs(d[0] - t33))) > thresh) {
                    return false
                }

                dlarfx(Right, j3 + 1, 3, u[..], tau, t[j1..], ldt, work)
                dlarfx(Left, 3, n - j1 - 1, u[..], tau, t[j1 * ldt + j2..], ldt, work)

                t[j1 * ldt + j1] = t33
                t[j2 * ldt + j1] = 0.0
                t[j3 * ldt + j1] = 0.0

                if (wantq) {
                    dlarfx(Right, n, 3, u[..], tau, q[j1..], ldq, work)
                }

            case _ =>
                let u1 = Array<Float64>([1.0, -x[ldx], scale])
                let (_, tau1) = dlarfg(3, -x[0], u1[1..], 1)

                let temp = -tau1 * (x[1] + Float64(u1[1]) * x[ldx + 1])
                let u2 = Array<Float64>([1.0, -temp * u1[2], scale])
                let (_, tau2) = dlarfg(3, -temp * u1[1] - x[ldx + 1], u2[1..], 1)

                dlarfx(Left, 3, 4, u1[..], tau1, d[..], ldd, work)
                dlarfx(Right, 4, 3, u1[..], tau1, d[..], ldd, work)
                dlarfx(Left, 3, 4, u2[..], tau2, d[ldd..], ldd, work)
                dlarfx(Right, 4, 3, u2[..], tau2, d[1..], ldd, work)

                let m1 = max(abs(d[2 * ldd]), abs(d[2 * ldd + 1]))
                let m2 = max(abs(d[3 * ldd]), abs(d[3 * ldd + 1]))
                if (max(m1, m2) > thresh) {
                    return false
                }

                let j4 = j1 + 3
                dlarfx(Left, 3, n - j1, u1[..], tau1, t[j1 * ldt + j1..], ldt, work)
                dlarfx(Right, j4 + 1, 3, u1[..], tau1, t[j1..], ldt, work)
                dlarfx(Left, 3, n - j1, u2[..], tau2, t[j2 * ldt + j1..], ldt, work)
                dlarfx(Right, j4 + 1, 3, u2[..], tau2, t[j2..], ldt, work)

                t[j3 * ldt + j1] = 0.0
                t[j3 * ldt + j2] = 0.0
                t[j4 * ldt + j1] = 0.0
                t[j4 * ldt + j2] = 0.0

                if (wantq) {
                    dlarfx(Right, n, 3, u1[..], tau1, q[j1..], ldq, work)
                    dlarfx(Right, n, 3, u2[..], tau2, q[j2..], ldq, work)
                }
        }

        if (n2 == 2) {
            let (a, b) = (t[j1 * ldt + j1], t[j1 * ldt + j2])
            let (c, d) = (t[j2 * ldt + j1], t[j2 * ldt + j2])
            var cs: Float64
            var sn: Float64
            (t[j1 * ldt + j1], t[j1 * ldt + j2], t[j2 * ldt + j1], t[j2 * ldt + j2], _, _, _, _, cs, sn) = dlanv2(a, b,
                c, d)
            if (n - j1 - 2 > 0) {
                nFloat64.drot(n - j1 - 2, t[j1 * ldt + j1 + 2..], 1, t[j2 * ldt + j1 + 2..], 1, cs, sn)
            }
            if (j1 > 0) {
                nFloat64.drot(j1, t[j1..], ldt, t[j2..], ldt, cs, sn)
            }
            if (wantq) {
                nFloat64.drot(n, q[j1..], ldq, q[j2..], ldq, cs, sn)
            }
        }
        if (n1 == 2) {
            let j3 = j1 + n2
            let j4 = j3 + 1
            let (a, b) = (t[j3 * ldt + j3], t[j3 * ldt + j4])
            let (c, d) = (t[j4 * ldt + j3], t[j4 * ldt + j4])
            var cs: Float64
            var sn: Float64
            (t[j3 * ldt + j3], t[j3 * ldt + j4], t[j4 * ldt + j3], t[j4 * ldt + j4], _, _, _, _, cs, sn) = dlanv2(a, b,
                c, d)
            if (n - j3 - 2 > 0) {
                nFloat64.drot(n - j3 - 2, t[j3 * ldt + j3 + 2..], 1, t[j4 * ldt + j3 + 2..], 1, cs, sn)
            }
            nFloat64.drot(j3, t[j3..], ldt, t[j4..], ldt, cs, sn)
            if (wantq) {
                nFloat64.drot(n, q[j3..], ldq, q[j4..], ldq, cs, sn)
            }
        }

        return true
    }
}
