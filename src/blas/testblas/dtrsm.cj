package cjnum.blas.testblas

import std.random.*
import cjnum.blas.*
import cjnum.floats.*

class dtrsmTestStruct {
    dtrsmTestStruct(
        let s!: Side,
        let ul!: Uplo,
        let tA!: Transpose,
        let d!: Diag,
        let m!: Int64,
        let n!: Int64,
        let alpha!: Float64,
        let a!: Array<Array<Float64>>,
        let b!: Array<Array<Float64>>,
        let want!: Array<Array<Float64>>
    ) {}
}

public func dtrsmTest(nFloat64: NFloat64): Bool {
    let tests = [
        dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: NoTrans,
            d: NonUnitDiag,
            m: 3,
            n: 2,
            alpha: 2.0,
            a: [
                [1.0, 2.0, 3.0],
                [0.0, 4.0, 5.0],
                [0.0, 0.0, 5.0]
            ],
            b: [
                [3.0, 6.0],
                [4.0, 7.0],
                [5.0, 8.0]
            ],
            want: [
                [1.0, 3.4],
                [-0.5, -0.5],
                [2.0, 3.2]
            ]
        ),
		dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: NoTrans,
            d: UnitDiag,
			m: 3,
            n: 2,
            alpha: 2.0,
            a: [
                [1.0, 2.0, 3.0],
				[0.0, 4.0, 5.0],
				[0.0, 0.0, 5.0]
            ],
            b: [
                [3.0, 6.0],
                [4.0, 7.0],
				[5.0, 8.0]
            ],
            want: [
                [60.0, 96.0],
                [-42.0, -66.0],
				[10.0, 16.0]
			]
		),
		dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: NoTrans,
            d: NonUnitDiag,
			m: 3,
            n: 4,
            alpha: 2.0,
            a: [
                [1.0, 2.0, 3.0],
				[0.0, 4.0, 5.0],
				[0.0, 0.0, 5.0]
			],
			b: [
                [3.0, 6.0, 2.0, 9.0],
                [4.0, 7.0, 1.0, 3.0],
                [5.0, 8.0, 9.0, 10.0]
            ],
			want: [
                [1.0, 3.4, 1.2, 13.0],
                [-0.5, -0.5, -4.0, -3.5],
                [2.0, 3.2, 3.6, 4.0]
            ]
		),
		dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: NoTrans,
            d: UnitDiag,
			m: 3,
			n: 4,
			alpha: 2.0,
            a: [
                [1.0, 2.0, 3.0],
                [0.0, 4.0, 5.0],
                [0.0, 0.0, 5.0]
			],
			b: [
                [3.0, 6.0, 2.0, 9.0],
                [4.0, 7.0, 1.0, 3.0],
                [5.0, 8.0, 9.0, 10.0]
            ],
			want: [
                [60.0, 96.0, 126.0, 146.0],
                [-42.0, -66.0, -88.0, -94.0],
                [10.0, 16.0, 18.0, 20.0]
            ]
		),
		dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: NoTrans,
            d: NonUnitDiag,
			m: 3,
			n: 2,
			alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 4.0, 0.0],
                [5.0, 6.0, 7.0]
			],
			b: [
                [3.0, 6.0],
                [4.0, 7.0],
                [5.0, 8.0]
            ],
			want: [
                [4.5, 9.0],
                [-0.375, -1.5],
                [-0.75, -12.0 / 7.0]
            ]
		),
		dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: NoTrans,
            d: UnitDiag,
			m: 3,
			n: 2,
			alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 4.0, 0.0],
                [5.0, 6.0, 7.0]
			],
			b: [
                [3.0, 6.0],
                [4.0, 7.0],
                [5.0, 8.0]
            ],
			want: [
				[9.0, 18.0],
				[-15.0, -33.0],
                [60.0, 132.0]
			]
		),
		dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: NoTrans,
            d: NonUnitDiag,
			m: 3,
            n: 4,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
				[3.0, 4.0, 0.0],
				[5.0, 6.0, 7.0]
			],
			b: [
                [3.0, 6.0, 2.0, 9.0],
                [4.0, 7.0, 1.0, 3.0],
                [5.0, 8.0, 9.0, 10.0]
            ],
			want: [
                [4.5, 9.0, 3.0, 13.5],
                [-0.375, -1.5, -1.5, -63.0 / 8.0],
                [-0.75, -12.0 / 7.0, 3.0, 39.0 / 28.0]
            ]
		),
		dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: NoTrans,
            d: UnitDiag,
			m: 3,
            n: 4,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
				[3.0, 4.0, 0.0],
				[5.0, 6.0, 7.0]
			],
			b: [
                [3.0, 6.0, 2.0, 9.0],
                [4.0, 7.0, 1.0, 3.0],
                [5.0, 8.0, 9.0, 10.0]
            ],
			want: [
                [9.0, 18.0, 6.0, 27.0],
				[-15.0, -33.0, -15.0, -72.0],
                [60.0, 132.0, 87.0, 327.0]
            ]
		),
        dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: Trans,
            d: NonUnitDiag,
            m: 3,
            n: 2,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 7.0]
            ],
            b: [
                [3.0, 6.0],
                [4.0, 7.0],
                [5.0, 8.0]
            ],
            want: [
                [4.5, 9.0],
                [-0.3, -1.2],
                [-6.0 / 35.0, -24.0 / 35.0]
            ]
        ),
        dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: Trans,
            d: UnitDiag,
            m: 3,
            n: 2,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 7.0]
            ],
            b: [
                [3.0, 6.0],
                [4.0, 7.0],
                [5.0, 8.0]
            ],
            want: [
                [9.0, 18.0],
                [-15.0, -33.0],
                [69.0, 150.0]
            ]
        ),
        dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: Trans,
            d: NonUnitDiag,
            m: 3,
            n: 4,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 7.0]
            ],
            b: [
                [3.0, 6.0, 6.0, 7.0],
                [4.0, 7.0, 8.0, 9.0],
                [5.0, 8.0, 10.0, 11.0]
            ],
            want: [
                [4.5, 9.0, 9.0, 10.5],
                [-0.3, -1.2, -0.6, -0.9],
                [-6.0 / 35.0, -24.0 / 35.0, -12.0 / 35.0, -18.0 / 35.0]
            ]
        ),
        dtrsmTestStruct(
            s: Left,
            ul: Upper,
            tA: Trans,
            d: UnitDiag,
            m: 3,
            n: 4,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 7.0]
            ],
            b: [
                [3.0, 6.0, 6.0, 7.0],
                [4.0, 7.0, 8.0, 9.0],
                [5.0, 8.0, 10.0, 11.0]
            ],
            want: [
                [9.0, 18.0, 18.0, 21.0],
                [-15.0, -33.0, -30.0, -36.0],
                [69.0, 150.0, 138.0, 165.0]
            ]
        ),
        dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: Trans,
            d: NonUnitDiag,
            m: 3,
            n: 2,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 4.0, 0.0],
                [5.0, 6.0, 8.0]
            ],
            b: [
                [3.0, 6.0],
                [4.0, 7.0],
                [5.0, 8.0]
            ],
            want: [
                [-0.46875, 0.375],
                [0.1875, 0.75],
                [1.875, 3.0]
            ]
        ),
        dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: Trans,
            d: UnitDiag,
            m: 3,
            n: 2,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 4.0, 0.0],
                [5.0, 6.0, 8.0]
            ],
            b: [
                [3.0, 6.0],
                [4.0, 7.0],
                [5.0, 8.0]
            ],
            want: [
                [168.0, 267.0],
                [-78.0, -123.0],
                [15.0, 24.0]
            ]
        ),
        dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: Trans,
            d: NonUnitDiag,
            m: 3,
            n: 4,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 4.0, 0.0],
                [5.0, 6.0, 8.0]
            ],
            b: [
                [3.0, 6.0, 2.0, 3.0],
                [4.0, 7.0, 4.0, 5.0],
                [5.0, 8.0, 6.0, 7.0]
            ],
            want: [
                [-0.46875, 0.375, -2.0625, -1.78125],
                [0.1875, 0.75, -0.375, -0.1875],
                [1.875, 3.0, 2.25, 2.625]
            ],
        ),
		dtrsmTestStruct(
            s: Left,
            ul: Lower,
            tA: Trans,
            d: UnitDiag,
			m: 3,
            n: 4,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
				[3.0, 4.0, 0.0],
				[5.0, 6.0, 8.0]
			],
			b: [
                [3.0, 6.0, 2.0, 3.0],
                [4.0, 7.0, 4.0, 5.0],
                [5.0, 8.0, 6.0, 7.0]
            ],
			want: [
                [168.0, 267.0, 204.0, 237.0],
                [-78.0, -123.0, -96.0, -111.0],
                [15.0, 24.0, 18.0, 21.0]
            ]
		),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: NoTrans,
            d: NonUnitDiag,
            m: 4,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 7.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
                [16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
            want: [
                [15.0, -2.4, -48.0 / 35.0],
                [19.5, -3.3, -66.0 / 35.0],
                [24.0, -4.2, -2.4],
                [28.5, -5.1, -102.0 / 35.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: NoTrans,
            d: UnitDiag,
            m: 4,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
                [16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
            want: [
                [30.0, -57.0, 258.0],
                [39.0, -75.0, 339.0],
                [48.0, -93.0, 420.0],
                [57.0, -111.0, 501.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: NoTrans,
            d: NonUnitDiag,
            m: 2,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 7.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
            want: [
                [15.0, -2.4, -48.0 / 35.0],
                [19.5, -3.3, -66.0 / 35.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: NoTrans,
            d: UnitDiag,
            m: 2,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
            want: [
                [30.0, -57.0, 258.0],
                [39.0, -75.0, 339.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: NoTrans,
            d: NonUnitDiag,
            m: 4,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
                [16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
            want: [
                [4.2, 1.2, 4.5],
                [5.775, 1.65, 5.625],
                [7.35, 2.1, 6.75],
                [8.925, 2.55, 7.875]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: NoTrans,
            d: UnitDiag,
            m: 4,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
                [16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
            want: [
                [435.0, -183.0, 36.0],
                [543.0, -228.0, 45.0],
                [651.0, -273.0, 54.0],
                [759.0, -318.0, 63.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: NoTrans,
            d: NonUnitDiag,
            m: 2,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
            want: [
                [4.2, 1.2, 4.5],
                [5.775, 1.65, 5.625]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: NoTrans,
            d: UnitDiag,
            m: 2,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
            want: [
                [435.0, -183.0, 36.0],
                [543.0, -228.0, 45.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: Trans,
            d: NonUnitDiag,
            m: 4,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
				[16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
            want: [
                [4.2, 1.2, 4.5],
                [5.775, 1.65, 5.625],
				[7.35, 2.1, 6.75],
				[8.925, 2.55, 7.875]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: Trans,
            d: UnitDiag,
            m: 4,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
				[16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
            want: [
				[435.0, -183.0, 36.0],
				[543.0, -228.0, 45.0],
                [651.0, -273.0, 54.0],
                [759.0, -318.0, 63.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: Trans,
            d: NonUnitDiag,
            m: 2,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
            want: [
                [4.2, 1.2, 4.5],
                [5.775, 1.65, 5.625]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Upper,
            tA: Trans,
            d: UnitDiag,
            m: 2,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 3.0, 4.0],
                [0.0, 5.0, 6.0],
                [0.0, 0.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
            want: [
				[435.0, -183.0, 36.0],
				[543.0, -228.0, 45.0]
            ]
        ),
        dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: Trans,
            d: NonUnitDiag,
            m: 4,
            n: 3,
            alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
				[3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
            ],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
				[16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
            want: [
                [15.0, -2.4, -1.2],
				[19.5, -3.3, -1.65],
				[24.0, -4.2, -2.1],
				[28.5, -5.1, -2.55]
            ]
        ),
		dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: Trans,
            d: UnitDiag,
			m: 4,
			n: 3,
			alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
			],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0],
				[16.0, 17.0, 18.0],
                [19.0, 20.0, 21.0]
            ],
			want: [
				[30.0, -57.0, 258.0],
				[39.0, -75.0, 339.0],
                [48.0, -93.0, 420.0],
                [57.0, -111.0, 501.0]
			]
		),
		dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: Trans,
            d: NonUnitDiag,
			m: 2,
			n: 3,
			alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
			],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
			want: [
				[15.0, -2.4, -1.2],
				[19.5, -3.3, -1.65]
			]
		),
		dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: Trans,
            d: UnitDiag,
			m: 2,
			n: 3,
			alpha: 3.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
			],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
			want: [
				[30.0, -57.0, 258.0],
				[39.0, -75.0, 339.0]
			]
		),
		dtrsmTestStruct(
            s: Right,
            ul: Lower,
            tA: Trans,
            d: UnitDiag,
			m: 2,
			n: 3,
			alpha: 0.0,
            a: [
                [2.0, 0.0, 0.0],
                [3.0, 5.0, 0.0],
                [4.0, 6.0, 8.0]
			],
            b: [
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ],
			want: [
				[0.0, 0.0, 0.0],
				[0.0, 0.0, 0.0]
			]
		)
    ]

    let r = Random(1)
    for (i in 0..tests.size) {
        let test = tests[i]
        let m = test.m
        let n = test.n
        let na = if (test.s == Right) {
            n
        } else {
            m
        }
        for (lda in [na, na + 3]) {
            for (ldb in [n, n + 5]) {
                let a = Array<Float64>(lda * na, { _ =>
                    r.nextFloat64()
                })
                for (i in 0..na) {
                    for (j in 0..na) {
                        a[i * lda + j] = test.a[i][j]
                    }
                }

                let b = Array<Float64>(ldb * m, { _ =>
                    r.nextFloat64()
                })
                for (i in 0..m) {
                    for (j in 0..n) {
                        b[i * ldb + j] = test.b[i][j]
                    }
                }

                nFloat64.dtrsm(test.s, test.ul, test.tA, test.d, test.m, test.n, test.alpha, a, lda, b, ldb)

                let want = b.clone()
                for (i in 0..m) {
                    for (j in 0..n) {
                        want[i * ldb + j] = test.want[i][j]
                    }
                }

                if (!fEqualApprox(want, b, 1e-13)) {
                    throw Exception("Case ${i}: Want ${want}, got ${b}.")
                }
            }
        }
    }
    return true
}