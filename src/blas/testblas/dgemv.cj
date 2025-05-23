package cjnum.blas.testblas

import std.math.*
import cjnum.blas.*
import cjnum.floats.*

class DgemvTestStruct {
    DgemvTestStruct(
        let name!: String,
        let m!: Int64,
        let n!: Int64,
        let a!: Array<Array<Float64>>,
        let tA!: Transpose,
        let x!: Array<Float64>,
        let incX!: Int64,
        let y!: Array<Float64>,
        let incY!: Int64,
        let subcases!: Array<DgemvSubcase>
    ) {}
}

class DgemvSubcase {
    DgemvSubcase(
        let mulXNeg1!: Bool = false,
        let mulYNeg1!: Bool = false,
        let alpha!:    Float64,
        let beta!:     Float64,
        let ans!:      Array<Float64>
    ) {}
}

let dgemvTests = [
    DgemvTestStruct(
        name: "M_gt_N_Inc1_NoTrans",
        tA: NoTrans,
        m: 5,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0],
            [1.0, 1.0, 2.0],
            [9.0, 2.0, 5.0]
        ],
        incX: 1,
        incY: 1,
        x: [1.0, 2.0, 3.0],
        y: [7.0, 8.0, 9.0, 10.0, 11.0],
        subcases: [
            DgemvSubcase(
                alpha: 0.0,
                beta: 0.0,
                ans: [0.0, 0.0, 0.0, 0.0, 0.0]
            ),
            DgemvSubcase(
                alpha: 0.0,
                beta: 1.0,
                ans: [7.0, 8.0, 9.0, 10.0, 11.0]
            ),
            DgemvSubcase(
                alpha: 1.0,
                beta: 0.0,
                ans: [40.8, 43.9, 33.0, 9.0, 28.0]
            ),
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [284.4, 303.2, 210.0, 12.0, 158.0]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_gt_N_Inc1_Trans",
        tA: Trans,
        m: 5,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0],
            [1.0, 1.0, 2.0],
            [9.0, 2.0, 5.0]
        ],
        incX: 1,
        incY: 1,
        x: [1.0, 2.0, 3.0, -4.0, 5.0],
        y: [7.0, 8.0, 9.0],
        subcases: [
            DgemvSubcase(
                alpha: 0.0,
                beta: 0.0,
                ans: [0.0, 0.0, 0.0]
            ),
            DgemvSubcase(
                alpha: 0.0,
                beta: 1.0,
                ans: [7.0, 8.0, 9.0]
            ),
            DgemvSubcase(
                alpha: 1.0,
                beta: 0.0,
                ans: [94.3, 40.2, 52.3]
            ),
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [712.4, 273.6, 364.4]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_eq_N_Inc1_NoTrans",
        tA: NoTrans,
        m: 3,
        n: 3,
        a: [
			[4.1, 6.2, 8.1],
			[9.6, 3.5, 9.1],
			[10.0, 7.0, 3.0]
        ],
        incX: 1,
        incY: 1,
        x: [1.0, 2.0, 3.0],
        y: [7.0, 2.0, 2.0],
        subcases: [
            DgemvSubcase(
                alpha: 0.0,
                beta: 0.0,
                ans: [0.0, 0.0, 0.0]
            ),
            DgemvSubcase(
                alpha: 0.0,
                beta: 1.0,
                ans: [7.0, 2.0, 2.0]
            ),
            DgemvSubcase(
                alpha: 1.0,
                beta: 0.0,
                ans: [40.8, 43.9, 33.0]
            ),
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [40.8 * 8.0 - 6.0 * 7.0, 43.9 * 8.0 - 6.0 * 2.0, 33.0 * 8.0 - 6.0 * 2.0]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_eq_N_Inc1_Trans",
        tA: Trans,
        m: 3,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0]
        ],
        incX: 1,
        incY: 1,
        x: [1.0, 2.0, 3.0],
        y: [7.0, 2.0, 2.0],
        subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [384.4, 261.6, 270.4]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_lt_N_Inc1_NoTrans",
        tA: NoTrans,
        m: 3,
        n: 5,
        a: [
            [4.1, 6.2, 8.1, 10.0, 7.0],
            [9.6, 3.5, 9.1, -2.0, 9.0],
            [10.0, 7.0, 3.0, 1.0, -5.0]
        ],
        incX: 1,
        incY: 1,
        x: [1.0, 2.0, 3.0, -7.6, 8.1],
        y: [7.0, 2.0, 2.0],
        subcases: [
            DgemvSubcase(
                alpha: 0.0,
                beta: 0.0,
                ans: [0.0, 0.0, 0.0]
            ),
            DgemvSubcase(
                alpha: 0.0,
                beta: 1.0,
                ans: [7.0, 2.0, 2.0]
            ),
            DgemvSubcase(
                alpha: 1.0,
                beta: 0.0,
                ans: [21.5, 132.0, -15.1]
            ),
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [21.5 * 8.0 - 6.0 * 7.0, 132.0 * 8.0 - 6.0 * 2.0, -15.1 * 8.0 - 6.0 * 2.0]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_lt_N_Inc1_Trans",
        tA: Trans,
        m: 3,
        n: 5,
        a: [
			[4.1, 6.2, 8.1, 10.0, 7.0],
			[9.6, 3.5, 9.1, -2.0, 9.0],
			[10.0, 7.0, 3.0, 1.0, -5.0]
        ],
        incX: 1,
        incY: 1,
        x: [1.0, 2.0, 3.0],
        y: [7.0, 2.0, 2.0, -3.0, 5.0],
        subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [384.4, 261.6, 270.4, 90.0, 50.0]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_gt_N_Part1_NoTrans",
        tA: NoTrans,
        m: 5,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0],
            [1.0, 1.0, 2.0],
            [9.0, 2.0, 5.0]
        ],
        incX: 1,
        incY: 2,
        x: [1.0, 2.0, 3.0],
        y: [7.0, 100.0, 8.0, 101.0, 9.0, 102.0, 10.0, 103.0, 11.0],
        subcases: [
            DgemvSubcase(
                alpha: 0.0,
                beta: 0.0,
                ans: [0.0, 100.0, 0.0, 101.0, 0.0, 102.0, 0.0, 103.0, 0.0]
            ),
            DgemvSubcase(
                alpha: 0.0,
                beta: 1.0,
                ans: [7.0, 100.0, 8.0, 101.0, 9.0, 102.0, 10.0, 103.0, 11.0]
            ),
            DgemvSubcase(
                alpha: 1.0,
                beta: 0.0,
                ans: [40.8, 100.0, 43.9, 101.0, 33.0, 102.0, 9.0, 103.0, 28.0]
            ),
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [284.4, 100.0, 303.2, 101.0, 210.0, 102.0, 12.0, 103.0, 158.0]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_gt_N_Part1_Trans",
        tA: Trans,
        m: 5,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0],
            [1.0, 1.0, 2.0],
            [9.0, 2.0, 5.0]
        ],
        incX: 1,
        incY: 2,
        x: [1.0, 2.0, 3.0, -4.0, 5.0],
        y: [7.0, 100.0, 8.0, 101.0, 9.0],
        subcases: [
            DgemvSubcase(
                alpha: 0.0,
                beta: 0.0,
                ans: [0.0, 100.0, 0.0, 101.0, 0.0]
            ),
            DgemvSubcase(
                alpha: 0.0,
                beta: 1.0,
                ans: [7.0, 100.0, 8.0, 101.0, 9.0]
            ),
            DgemvSubcase(
                alpha: 1.0,
                beta: 0.0,
                ans: [94.3, 100.0, 40.2, 101.0, 52.3]
            ),
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [712.4, 100.0, 273.6, 101.0, 364.4]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_gt_N_IncNot1_NoTrans",
        tA: NoTrans,
        m: 5,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0],
            [1.0, 1.0, 2.0],
            [9.0, 2.0, 5.0]
        ],
        incX: 2,
        incY: 3,
		x: [1.0, 15.0, 2.0, 150.0, 3.0],
		y: [7.0, 2.0, 6.0, 8.0, -4.0, -5.0, 9.0, 1.0, 1.0, 10.0, 19.0, 22.0, 11.0],
        subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [284.4, 2.0, 6.0, 303.2, -4.0, -5.0, 210.0, 1.0, 1.0, 12.0, 19.0, 22.0, 158.0]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [220.4, 2.0, 6.0, 311.2, -4.0, -5.0, 322.0, 1.0, 1.0, -4.0, 19.0, 22.0, 222.0]
            ),
            DgemvSubcase(
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [182.0, 2.0, 6.0, 24.0, -4.0, -5.0, 210.0, 1.0, 1.0, 291.2, 19.0, 22.0, 260.4]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [246.0, 2.0, 6.0, 8.0, -4.0, -5.0, 322.0, 1.0, 1.0, 299.2, 19.0, 22.0, 196.4]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_gt_N_IncNot1_Trans",
        tA: Trans,
        m: 5,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0],
            [1.0, 1.0, 2.0],
            [9.0, 2.0, 5.0]
        ],
        incX: 2,
        incY: 3,
		x: [1.0, 15.0, 2.0, 150.0, 3.0, 8.0, -3.0, 6.0, 5.0],
		y: [7.0, 2.0, 6.0, 8.0, -4.0, -5.0, 9.0],
        subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [720.4, 2.0, 6.0, 281.6, -4.0, -5.0, 380.4]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [219.6, 2.0, 6.0, 316.0, -4.0, -5.0, 195.6]
            ),
            DgemvSubcase(
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [392.4, 2.0, 6.0, 281.6, -4.0, -5.0, 708.4]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [207.6, 2.0, 6.0, 316.0, -4.0, -5.0, 207.6]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_eq_N_IncNot1_NoTrans",
        tA: NoTrans,
        m: 3,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0]
        ],
        incX: 2,
        incY: 3,
        x: [1.0, 15.0, 2.0, 150.0, 3.0],
        y: [7.0, 2.0, 6.0, 8.0, -4.0, -5.0, 9.0],
        subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [284.4, 2.0, 6.0, 303.2, -4.0, -5.0, 210.0]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [220.4, 2.0, 6.0, 311.2, -4.0, -5.0, 322.0]
            ),
            DgemvSubcase(
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [222.0, 2.0, 6.0, 303.2, -4.0, -5.0, 272.4]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [334.0, 2.0, 6.0, 311.2, -4.0, -5.0, 208.4]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_eq_N_IncNot1_Trans",
        tA: Trans,
        m: 3,
        n: 3,
        a: [
            [4.1, 6.2, 8.1],
            [9.6, 3.5, 9.1],
            [10.0, 7.0, 3.0]
        ],
        incX: 2,
        incY: 3,
        x: [1.0, 15.0, 2.0, 150.0, 3.0],
        y: [7.0, 2.0, 6.0, 8.0, -4.0, -5.0, 9.0],
        subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [384.4, 2.0, 6.0, 225.6, -4.0, -5.0, 228.4]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [290.0, 2.0, 6.0, 212.8, -4.0, -5.0, 310.0]
            ),
            DgemvSubcase(
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [240.4, 2.0, 6.0, 225.6, -4.0, -5.0, 372.4]
            ),
            DgemvSubcase(
                mulXNeg1: true,
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [322.0, 2.0, 6.0, 212.8, -4.0, -5.0, 278.0]
            )
        ]
    ),
    DgemvTestStruct(
        name: "M_lt_N_IncNot1_NoTrans",
        tA: NoTrans,
        m: 3,
        n: 5,
        a: [
			[4.1, 6.2, 8.1, 10.0, 11.0],
			[9.6, 3.5, 9.1, -3.0, -2.0],
            [10.0, 7.0, 3.0, -7.0, -4.0]
        ],
        incX: 2,
        incY: 3,
		x: [1.0, 15.0, 2.0, 150.0, 3.0, -2.0, -4.0, 8.0, -9.0],
		y: [7.0, 2.0, 6.0, 8.0, -4.0, -5.0, 9.0],
        subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [-827.6, 2.0, 6.0, 543.2, -4.0, -5.0, 722.0]
            ),
			DgemvSubcase(
                mulXNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [-93.2, 2.0, 6.0, -696.8, -4.0, -5.0, -1070.0]
			),
			DgemvSubcase(
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [734.0, 2.0, 6.0, 543.2, -4.0, -5.0, -839.6]
            ),
			DgemvSubcase(
                mulXNeg1: true,
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
				ans: [-1058.0, 2.0, 6.0, -696.8, -4.0, -5.0, -105.2]
			)
        ]
    ),
	DgemvTestStruct(
        name: "M_lt_N_IncNot1_Trans",
        tA: Trans,
        m: 3,
        n: 5,
        a: [
			[4.1, 6.2, 8.1, 10.0, 11.0],
			[9.6, 3.5, 9.1, -3.0, -2.0],
            [10.0, 7.0, 3.0, -7.0, -4.0]
		],
		incX: 2,
		incY: 3,
		x: [1.0, 15.0, 2.0, 150.0, 3.0],
        y: [7.0, 2.0, 6.0, 8.0, -4.0, -5.0, 9.0, -4.0, -1.0, -9.0, 1.0, 1.0, 2.0],
		subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [384.4, 2.0, 6.0, 225.6, -4.0, -5.0, 228.4, -4.0, -1.0, -82.0, 1.0, 1.0, -52.0]
            ),
			DgemvSubcase(
                mulXNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [290.0, 2.0, 6.0, 212.8, -4.0, -5.0, 310.0, -4.0, -1.0, 190.0, 1.0, 1.0, 188.0]
            ),
			DgemvSubcase(
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [-82.0, 2.0, 6.0, -184.0, -4.0, -5.0, 228.4, -4.0, -1.0, 327.6, 1.0, 1.0, 414.4]
            ),
			DgemvSubcase(
                mulXNeg1: true,
				mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [158.0, 2.0, 6.0, 88.0, -4.0, -5.0, 310.0, -4.0, -1.0, 314.8, 1.0, 1.0, 320.0]
            )
		]
	),
	DgemvTestStruct(
        name: "M_eq_N_Lg_IncNot1_Trans",
        tA: Trans,
        m: 7,
        n: 7,
        a: [
			[4.1, 6.2, 8.1, 2.5, 3.3, 7.4, 9.3],
			[9.6, 3.5, 9.1, 1.2, 5.4, 4.8, 8.7],
			[10.0, 7.0, 3.0, 2.0, 4.0, 1.0, 12.0],
			[9.6, 3.5, 9.1, 1.2, 5.4, 4.8, 8.7],
			[4.1, 6.2, 8.1, 2.5, 3.3, 7.4, 9.3],
			[10.0, 7.0, 3.0, 2.0, 4.0, 1.0, 12.0],
			[9.6, 3.5, 9.1, 1.2, 5.4, 4.8, 8.7]
		],
		incX: 2,
		incY: 3,
		x: [1.0, 105.0, 2.0, 150.0, 3.0, 200.0, 4.0, 300.0, 5.0, 600.0, 6.0, 700.0, 7.0],
        y: [7.0, 200.0, 600.0, 8.0, -400.0, -500.0, 9.0, 600.0, 700.0, 10.0, 500.0, 400.0, 11.0, 200.0, 300.0, 12.0, 100.0, 200.0, 13.0, 300.0, 400.0, 14.0],
		subcases: [
            DgemvSubcase(
                alpha: 8.0,
                beta: -6.0,
                ans: [1873.1999999999998, 200.0, 600.0, 1117.6, -400.0, -500.0, 1497.1999999999998, 600.0, 700.0, 328.8, 500.0, 400.0, 942.0, 200.0, 300.0, 854.4000000000001, 100.0, 200.0, 2137.2, 300.0, 400.0, 14.0]
            ),
			DgemvSubcase(
                mulXNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [1690.8, 200.0, 600.0, 1148.0, -400.0, -500.0, 1562.8, 600.0, 700.0, 357.6, 500.0, 400.0, 897.2, 200.0, 300.0, 998.4, 100.0, 200.0, 2103.6000000000004, 300.0, 400.0, 14.0]
            ),
			DgemvSubcase(
                mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [2173.2, 200.0, 600.0, 878.4000000000001, -400.0, -500.0, 954.0, 600.0, 700.0, 328.8, 500.0, 400.0, 1485.1999999999998, 200.0, 300.0, 1093.6, 100.0, 200.0, 1837.1999999999998, 300.0, 400.0, 14.0]
            ),
			DgemvSubcase(
                mulXNeg1: true,
				mulYNeg1: true,
                alpha: 8.0,
                beta: -6.0,
                ans: [2139.6, 200.0, 600.0, 1022.4, -400.0, -500.0, 909.2, 600.0, 700.0, 357.6, 500.0, 400.0, 1550.8, 200.0, 300.0, 1124.0, 100.0, 200.0, 1654.8, 300.0, 400.0, 14.0]
            )
		]
	)
]

public func dgemvTest(nFloat64: NFloat64): Bool {
    for (test in dgemvTests) {
        for (i in 0..test.subcases.size) {
            let cas = test.subcases[i]
            dgemvcomp(test, cas, i, nFloat64)
            dgemvbad(test, cas, i, nFloat64)
        }
    }
    return true
}

func dgemvcomp(test: DgemvTestStruct, cas: DgemvSubcase, i: Int64, nFloat64: NFloat64): Unit {
    let x = sliceCopy(test.x)
    let y = sliceCopy(test.y)
    let a = sliceOfSliceCopy(test.a)
    let aFlat = flatten(a)
    let lda = test.n

    let incX = if (cas.mulXNeg1) {
        -test.incX
    } else {
        test.incX
    }
    let incY = if (cas.mulYNeg1) {
        -test.incY
    } else {
        test.incY
    }

    let f = { => 
        nFloat64.dgemv(test.tA, test.m, test.n, cas.alpha, aFlat, lda, x, incX, cas.beta, y, incY)
    }
    if (exception(f)) {
        throw Exception("Test ${test.name} case ${i}: unexpected panic")
        // TODO fix where these codes are here?
        // if (throwException) {
        //     nFloat64.dgemv(test.tA, test.m, test.n, cas.alpha, aFlat, lda, x, incX, cas.beta, y, incY)
        // }
        // return 
    }

    if (!dSliceEqual(x, test.x)) {
        throw Exception("Test ${test.name}, case ${i}: x modified during call")
    }
    let aFlat2 = flatten(sliceOfSliceCopy(test.a))
    if (!dSliceEqual(aFlat2, aFlat)) {
        throw Exception("Test ${test.name}, case ${i}: a modified during call")
    }

    if (!dSliceTolEqual(cas.ans, y)) {
        throw Exception("Test ${test.name}, case ${i}: answer mismatch: Expected ${cas.ans}, Found ${y}")
    }
}

func dgemvbad(test: DgemvTestStruct, cas: DgemvSubcase, i: Int64, nFloat64: NFloat64): Unit {
    let x = sliceCopy(test.x)
    let y = sliceCopy(test.y)
    let a = sliceOfSliceCopy(test.a)
    let aFlatRow = flatten(a)
    let ldaRow = test.n

    var f = { => 
        nFloat64.dgemv(b'X', test.m, test.n, cas.alpha, aFlatRow, ldaRow, x, test.incX, cas.beta, y, test.incY)
    }
    if (!exception(f)) {
        throw Exception("Test ${test.name} case ${i}: no panic for bad transpose")
    }

    f = { => 
        nFloat64.dgemv(test.tA, -2, test.n, cas.alpha, aFlatRow, ldaRow, x, test.incX, cas.beta, y, test.incY)
    }
    if (!exception(f)) {
        throw Exception("Test ${test.name} case ${i}: no panic for m negative")
    }

    f = { => 
        nFloat64.dgemv(test.tA, test.m, -4, cas.alpha, aFlatRow, ldaRow, x, test.incX, cas.beta, y, test.incY)
    }
    if (!exception(f)) {
        throw Exception("Test ${test.name} case ${i}: no panic for n negative")
    }

    f = { =>
        nFloat64.dgemv(test.tA, test.m, test.n, cas.alpha, aFlatRow, ldaRow, x, 0, cas.beta, y, test.incY)
    }
    if (!exception(f)) {
        throw Exception("Test ${test.name} case ${i}: no panic for incX zero")
    }

    f = { =>
        nFloat64.dgemv(test.tA, test.m, test.n, cas.alpha, aFlatRow, ldaRow, x, test.incX, cas.beta, y, 0)
    }
    if (!exception(f)) {
        throw Exception("Test ${test.name} case ${i}: no panic for incY zero")
    }

    f = { =>
        nFloat64.dgemv(test.tA, test.m, test.n, cas.alpha, aFlatRow, ldaRow - 1, x, test.incX, cas.beta, y, test.incY)
    }
    if (!exception(f)) {
        throw Exception("Test ${test.name} case ${i}: no panic for lda too small row major")
    }
}