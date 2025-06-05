package cjnum.lapack.cjnum

import cjnum.lapack.*
import std.math.*

extend LFloat64Implementation {
    public func dtrexc(compq: UpdateSchurComp, n: Int64, t: Array<Float64>, ldt: Int64, q: Array<Float64>, ldq: Int64,
        ifst: Int64, ilst: Int64, work: Array<Float64>): (Int64, Int64, Bool) {
        var ifst_: Int64 = ifst
        var ilst_: Int64 = ilst
        if (compq != UpdateSchur && compq != UpdateSchurNone) {
            throw Exception(badUpdateSchurComp)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (ldt < max(1, n)) {
            throw Exception(badLdT)
        }
        if (ldq < 1 || (compq == UpdateSchur && ldq < n)) {
            throw Exception(badLdQ)
        }
        if ((ifst_ < 0 || n <= ifst_) && n > 0) {
            throw Exception(badIfst)
        }
        if ((ilst_ < 0 || n <= ilst_) && n > 0) {
            throw Exception(badIlst)
        }

        if (n == 0) {
            return (ifst_, ilst_, true)
        }

        if (t.size < (n - 1) * ldt + n) {
            throw Exception(shortT)
        }
        if (compq == UpdateSchur && q.size < (n - 1) * ldq + n) {
            throw Exception(shortQ)
        }
        if (work.size < n) {
            throw Exception(shortWork)
        }

        if (n == 1) {
            return (ifst_, ilst_, true)
        }

        if (ifst_ > 0 && t[ifst_ * ldt + ifst_ - 1] != 0.0) {
            ifst_--
        }
        var nbf = 1
        if (ifst_ + 1 < n && t[(ifst_ + 1) * ldt + ifst_] != 0.0) {
            nbf = 2
        }

        if (ilst_ > 0 && t[ilst_ * ldt + ilst_ - 1] != 0.0) {
            ilst_--
        }
        var nbl = 1
        if (ilst_ + 1 < n && t[(ilst_ + 1) * ldt + ilst_] != 0.0) {
            nbl = 2
        }

        var ok = true
        let wantq = compq == UpdateSchur

        if (ifst_ == ilst_) {
            return (ifst_, ilst_, true)
        } else if (ifst_ < ilst_) {
            if (nbf == 2 && nbl == 1) {
                ilst_--
            } else if (nbf == 1 && nbl == 2) {
                ilst_++
            }
            var here = ifst_
            while (here < ilst_) {
                if (nbf == 1 || nbf == 2) {
                    var nbnext = 1
                    if (here + nbf + 1 < n && t[(here + nbf + 1) * ldt + here + nbf] != 0.0) {
                        nbnext = 2
                    }
                    ok = dlaexc(wantq, n, t, ldt, q, ldq, here, nbf, nbnext, work)
                    if (!ok) {
                        return (ifst_, here, false)
                    }
                    here += nbnext
                    if (nbf == 2 && t[(here + 1) * ldt + here] == 0.0) {
                        nbf = 3
                    }
                    continue
                }

                var nbnext = 1
                if (here + 3 < n && t[(here + 3) * ldt + here + 2] != 0.0) {
                    nbnext = 2
                }
                ok = dlaexc(wantq, n, t, ldt, q, ldq, here + 1, 1, nbnext, work)
                if (!ok) {
                    return (ifst_, here, false)
                }
                if (nbnext == 1) {
                    dlaexc(wantq, n, t, ldt, q, ldq, here, 1, nbnext, work)
                    here++
                    continue
                }

                if (t[(here + 2) * ldt + here + 1] == 0.0) {
                    nbnext = 1
                }
                if (nbnext == 2) {
                    ok = dlaexc(wantq, n, t, ldt, q, ldq, here, 1, nbnext, work)
                    if (!ok) {
                        return (ifst_, here, false)
                    }
                } else {
                    dlaexc(wantq, n, t, ldt, q, ldq, here, 1, 1, work)
                    dlaexc(wantq, n, t, ldt, q, ldq, here + 1, 1, 1, work)
                }
                here += 2
            }
            return (ifst_, here, true)
        } else {
            var here = ifst_
            while (here > ilst_) {
                var nbnext = 1
                if (here >= 2 && t[(here - 1) * ldt + here - 2] != 0.0) {
                    nbnext = 2
                }
                if (nbf == 1 || nbf == 2) {
                    ok = dlaexc(wantq, n, t, ldt, q, ldq, here - nbnext, nbnext, nbf, work)
                    if (!ok) {
                        return (ifst_, here, false)
                    }
                    here -= nbnext
                    if (nbf == 2 && t[(here + 1) * ldt + here] == 0.0) {
                        nbf = 3
                    }
                    continue
                }

                ok = dlaexc(wantq, n, t, ldt, q, ldq, here - nbnext, nbnext, 1, work)
                if (!ok) {
                    return (ifst_, here, false)
                }
                if (nbnext == 1) {
                    dlaexc(wantq, n, t, ldt, q, ldq, here, nbnext, 1, work)
                    here--
                    continue
                }

                if (t[here * ldt + here - 1] == 0.0) {
                    nbnext = 1
                }
                if (nbnext == 2) {
                    ok = dlaexc(wantq, n, t, ldt, q, ldq, here - 1, 2, 1, work)
                    if (!ok) {
                        return (ifst_, here, false)
                    }
                } else {
                    dlaexc(wantq, n, t, ldt, q, ldq, here, 1, 1, work)
                    dlaexc(wantq, n, t, ldt, q, ldq, here - 1, 1, 1, work)
                }
                here -= 2
            }
            return (ifst_, here, true)
        }
    }
}
