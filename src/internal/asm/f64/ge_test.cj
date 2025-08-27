package cjnum.internal.asm.f64

import std.unittest.*
import std.unittest.testmacro.*

class gerStruct {
    gerStruct(
        let x: Array<Float64>,
        let y: Array<Float64>,
        let a: Array<Float64>,
        let want: Array<Float64>
    ) {}
}

let gerTest = [
    gerStruct([2.0], [4.4], [10.0], [18.8]),
    gerStruct([-2.0, -3.0, 0.0], [-1.1, 5.0],
        [1.3, 2.4, 2.6, 2.8, -1.3, -4.3], [3.5, -7.6, 5.9, -12.2, -1.3, -4.3]),
    gerStruct([-2.0, 7.0, 12.0], [-1.1, 0.0, 6.0],
        [1.3, 2.4, 3.5, 2.6, 2.8, 3.3, -1.3, -4.3, -9.7],
        [3.5, 2.4, -8.5, -5.1, 2.8, 45.3, -14.5, -4.3, 62.3]),
    gerStruct([-2.0, -3.0, 0.0, 1.0, 2.0], [-1.1, 5.0, 0.0],
        [1.3, 2.4, 3.5, 2.6, 2.8, 3.3, -1.3, -4.3, -9.7, 8.0, 9.0, -10.0, -12.0, -14.0, -6.0],
        [3.5, -7.6, 3.5, 5.9, -12.2, 3.3, -1.3, -4.3, -9.7, 6.9, 14.0, -10.0, -14.2, -4.0, -6.0]),
    gerStruct([-2.0, -3.0, 0.0], [-1.1, 5.0, 0.0, 9.0, 19.0, 22.0],
        [1.3, 2.4, 3.5, 4.8, 1.11, -9.0, 2.6, 2.8, 3.3, -3.4, 6.2, -8.7, -1.3, -4.3, -9.7, -3.1, 8.9, 8.9]
        , [3.5, -7.6, 3.5, -13.2, -36.89, -53.0, 5.9, -12.2, 3.3, -30.4, -50.8, -74.7, -1.3, -4.3, -9.7,
            -3.1, 8.9, 8.9]),
    gerStruct([-2.0, 0.0, 2.0, 0.0, 7.0], [-1.1, 8.0, 7.0, 3.0, 5.0],
            [1.3, 2.4, 3.5, 2.2, 8.3, 2.6, 2.8, 3.3, 4.4, -1.5, -1.3, -4.3, -9.7, -8.8, 6.2, 8.0, 9.0, -10.0, -11.0,
                12.0, -12.0, -14.0, -6.0, -2.0, 4.0],
            [3.5, -13.6, -10.5, -3.8, -1.7, 2.6, 2.8, 3.3, 4.4, -1.5, -3.5, 11.7, 4.3, -2.8, 16.2, 8.0, 9.0, -10.0,
                -11.0, 12.0, -19.700000000000003, 42.0, 43.0, 19.0, 39.0]),
    gerStruct([-2.0, 8.0, 9.0, -3.0, -1.2, 5.0, 4.5], [-1.1, nan, 19.0, 11.0, -9.22, 7.0,
        3.3],
            [1.3, 2.4, 3.5, 4.8, 1.11, -9.0, 2.2, 2.6, 2.8, 3.3, -3.4, 6.2, -8.7, 5.1, -1.3, -4.3, -9.7, -3.1, 8.9, 8.9,
                8.0, 5.0, -2.5, 1.8, -3.6, 2.8, 4.9, 7.0, -1.3, -4.3, -9.7, -3.1, 8.9, 8.9, 8.0, 2.6, 2.8, 3.3, -3.4,
                6.2, -8.7, 5.1, 1.3, 2.4, 3.5, 4.8, 1.11, -9.0, 2.2],
            [3.5, nan, -34.5, -17.2, 19.55, -23.0, -4.4, -6.2, nan, 155.3, 84.6, -67.56, 47.3, 31.5, -11.2, nan, 161.3,
                95.9, -74.08, 71.9, 37.7, 8.3, nan, -55.2, -36.6, 30.46, -16.1, -2.9, 0.02, nan, -32.5, -16.3, 19.964,
                0.5, 4.04, -2.9, nan, 98.3, 51.6, -39.9, 26.3, 21.6, -3.65, nan, 89.0, 54.3, -40.38, 22.5, 17.05])
]

@Test
public func testGer() {
    const tol = 1e-15
    const xGdVal = -0.5
    const yGdVal = 1.5
    const aGdVal = 10.0
    const gdLn = 4

    for (i in 0..gerTest.size) {
        let m = gerTest[i].x.size
        let n = gerTest[i].y.size
        for (align in align2) {
            let xgLn = gdLn + align.x
            let ygLn = gdLn + align.y
            let agLn = gdLn + align.x ^ align.y

            let xg = guardVector(gerTest[i].x, xGdVal, xgLn)
            let yg = guardVector(gerTest[i].y, yGdVal, ygLn)

            let x = xg[xgLn..xg.size - xgLn]
            let y = yg[ygLn..yg.size - ygLn]

            let ag = guardVector(gerTest[i].a, aGdVal, agLn)
            let a = ag[agLn..ag.size - agLn]

            let alpha = 1.0
            ger(m, n, alpha, x, 1, y, 1, a, n)
            for (j in 0..gerTest[i].want.size) {
                @Expect(sameApprox(a[j], gerTest[i].want[j], tol))
                return
            }

            @Expect(isValidGuard(xg, xGdVal, xgLn))
            @Expect(isValidGuard(yg, yGdVal, ygLn))
            @Expect(isValidGuard(ag, aGdVal, agLn))
            @Expect(equalStrided(gerTest[i].x, x, 1))
            @Expect(equalStrided(gerTest[i].y, y, 1))
        }

        for (inc in newIncSet(1, 2)) {
            let xg = guardIncVector(gerTest[i].x, xGdVal, inc.x, gdLn)
            let yg = guardIncVector(gerTest[i].y, yGdVal, inc.y, gdLn)
            let x = xg[gdLn..xg.size - gdLn]
            let y = yg[gdLn..yg.size - gdLn]
            let ag = guardVector(gerTest[i].a, aGdVal, gdLn)
            let a = ag[gdLn..ag.size - gdLn]

            let alpha = 3.5
            ger(m, n, alpha, x, inc.x, y, inc.y, a, n)
            for (j in 0..gerTest[i].want.size) {
                let want = alpha * (gerTest[i].want[j] - gerTest[i].a[j]) + gerTest[i].a[j]
                @Expect(sameApprox(a[j], want, tol))
            }

            checkValidIncGuard(xg, xGdVal, inc.x, gdLn)
            checkValidIncGuard(yg, yGdVal, inc.y, gdLn)

            @Expect(isValidGuard(ag, aGdVal, gdLn))
            @Expect(equalStrided(gerTest[i].x, x, inc.x))
            @Expect(equalStrided(gerTest[i].x, x, inc.y))
        }
    }
}
