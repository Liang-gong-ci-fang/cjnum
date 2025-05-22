package cjnum.lapack.cjnum

import std.math.*
import cjnum.blas.*
import cjnum.lapack.*

extend LFloat64Implementation {
    public func dlasr(side: Side, pivot: Pivot, direct: Direct, m: Int64, n: Int64, c: Array<Float64>,
        s: Array<Float64>, a: Array<Float64>, lda: Int64): Unit {
        if (side != Left && side != Right) {
            throw Exception(badSide)
        }
        if (pivot != Variable && pivot != Top && pivot != Bottom) {
            throw Exception(badPivot)
        }
        if (direct != Forward && direct != Backward) {
            throw Exception(badDirect)
        }
        if (m < 0) {
            throw Exception(mLT0)
        }
        if (n < 0) {
            throw Exception(nLT0)
        }
        if (lda < max(1, n)) {
            throw Exception(badLdA)
        }

        if (m == 0 || n == 0) {
            return 
        }

        if (side == Left) {
            if (c.size < m - 1) {
                throw Exception(shortC)
            }
            if (s.size < m - 1) {
                throw Exception(shortS)
            }
        } else {
            if (c.size < n - 1) {
                throw Exception(shortC)
            }
            if (s.size < n - 1) {
                throw Exception(shortS)
            }
        }
        if (a.size < (m - 1) * lda + n) {
            throw Exception(shortA)
        }

        if (side == Left) {
            if (pivot == Variable) {
                if (direct == Forward) {
                    for (j in 0..m - 1) {
                        let ctmp = c[j]
                        let stmp = s[j]
                        if (ctmp != 1.0 || stmp != 0.0) {
                            for (i in 0..n) {
                                let tmp2 = a[j * lda + i]
                                let tmp = a[(j + 1) * lda + i]
                                a[(j + 1) * lda + i] = ctmp * tmp - stmp * tmp2
                                a[j * lda + i] = stmp * tmp + ctmp * tmp2
                            }
                        }
                    }
                    return
                }
                var j = m - 2
                while (j >= 0) {
                    let ctmp = c[j]
                    let stmp = s[j]
                    if (ctmp != 1.0 || stmp != 0.0) {
                        for (i in 0..n) {
                            let tmp2 = a[j * lda + i]
                            let tmp = a[(j + 1) * lda + i]
                            a[(j + 1) * lda + i] = ctmp * tmp - stmp * tmp2
                            a[j * lda + i] = stmp * tmp + ctmp * tmp2
                        }
                    }
                    j--
                }
                return
            } else if (pivot == Top) {
                if (direct == Forward) {
                    for (j in 1..m) {
                        let ctmp = c[j - 1]
                        let stmp = s[j - 1]
                        if (ctmp != 1.0 || stmp != 0.0) {
                            for (i in 0..n) {
                                let tmp = a[j * lda + i]
                                let tmp2 = a[i]
                                a[j * lda + i] = ctmp * tmp - stmp * tmp2
                                a[i] = stmp * tmp + ctmp * tmp2
                            }
                        }
                    }
                    return
                }
                var j = m - 1
                while (j >= 1) {
                    let ctmp = c[j - 1]
                    let stmp = s[j - 1]
                    if (ctmp != 1.0 || stmp != 0.0) {
                        for (i in 0..n) {
                            let ctmp = c[j - 1]
                            let stmp = s[j - 1]
                            if (ctmp != 1.0 || stmp != 0.0) {
                                for (i in 0..n) {
                                    let tmp = a[j * lda + i]
                                    let tmp2 = a[i]
                                    a[j * lda + i] = ctmp * tmp - stmp * tmp2
                                    a[i] = stmp * tmp + ctmp * tmp2
                                }
                            }
                        }
                    }
                    j--
                }
                return
            }
            if (direct == Forward) {
                for (j in 0..m - 1) {
                    let ctmp = c[j]
                    let stmp = s[j]
                    if (ctmp != 1.0 || stmp != 0.0) {
                        for (i in 0..n) {
                            let tmp = a[j * lda + i]
                            let tmp2 = a[(m - 1) * lda + i]
                            a[j * lda + i] = stmp * tmp2 + ctmp * tmp
                            a[(m - 1) * lda + i] = ctmp * tmp2 - stmp * tmp
                        }
                    }
                }
                return
            }
            var j = m - 2
            while (j >= 0) {
                let ctmp = c[j]
                let stmp = s[j]
                if (ctmp != 1.0 || stmp != 0.0) {
                    for (i in 0..n) {
                        let tmp = a[j * lda + i]
                        let tmp2 = a[(m - 1) * lda + i]
                        a[j * lda + i] = stmp * tmp2 + ctmp * tmp
                        a[(m - 1) * lda + i] = ctmp * tmp2 - stmp * tmp
                    }
                }
                j--
            }
            return
        }
        if (pivot == Variable) {
            if (direct == Forward) {
                for (j in 0..n - 1) {
                    let ctmp = c[j]
                    let stmp = s[j]
                    if (ctmp != 1.0 || stmp != 0.0) {
                        for (i in 0..m) {
                            let tmp = a[i * lda + j + 1]
                            let tmp2 = a[i * lda + j]
                            a[i * lda + j + 1] = ctmp * tmp - stmp * tmp2
                            a[i * lda + j] = stmp * tmp + ctmp * tmp2
                        }
                    }
                }
                return
            }
            var j = n - 2
            while (j >= 0) {
                let ctmp = c[j]
                let stmp = s[j]
                if (ctmp != 1.0 || stmp != 0.0) {
                    for (i in 0..m) {
                        let tmp = a[i * lda + j + 1]
                        let tmp2 = a[i * lda + j]
                        a[i * lda + j + 1] = ctmp * tmp - stmp * tmp2
                        a[i * lda + j] = stmp * tmp + ctmp * tmp2
                    }
                }
                j--
            }
            return
        } else if (pivot == Top) {
            if (direct == Forward) {
                for (j in 1..n) {
                    let ctmp = c[j - 1]
                    let stmp = s[j - 1]
                    if (ctmp != 1.0 || stmp != 0.0) {
                        for (i in 0..m) {
                            let tmp = a[i * lda + j]
                            let tmp2 = a[i * lda]
                            a[i * lda + j] = ctmp * tmp - stmp * tmp2
                            a[i * lda] = stmp * tmp + ctmp * tmp2
                        }
                    }
                }
                return
            }
            var j = n - 1
            while (j >= 1) {
                let ctmp = c[j - 1]
                let stmp = s[j - 1]
                if (ctmp != 1.0 || stmp != 0.0) {
                    for (i in 0..m) {
                        let tmp = a[i * lda + j]
                        let tmp2 = a[i * lda]
                        a[i * lda + j] = ctmp * tmp - stmp * tmp2
                        a[i * lda] = stmp * tmp + ctmp * tmp2
                    }
                }
                j--
            }
            return
        }
        if (direct == Forward) {
            for (j in 0..n - 1) {
                let ctmp = c[j]
                let stmp = s[j]
                if (ctmp != 1.0 || stmp != 0.0) {
                    for (i in 0..m) {
                        let tmp = a[i * lda + j]
                        let tmp2 = a[i * lda + n - 1]
                        a[i * lda + j] = stmp * tmp2 + ctmp * tmp
                        a[i * lda + n - 1] = ctmp * tmp2 - stmp * tmp
                    }
                }
            }
            return
        }
        var j = n - 2
        while (j >= 0) {
            let ctmp = c[j]
            let stmp = s[j]
            if (ctmp != 1.0 || stmp != 0.0) {
                for (i in 0..m) {
                    let tmp = a[i * lda + j]
                    let tmp2 = a[i * lda + n - 1]
                    a[i * lda + j] = stmp * tmp2 + ctmp * tmp
                    a[i * lda + n - 1] = ctmp * tmp2 - stmp * tmp
                }
            }
            j--
        }
    }
}
