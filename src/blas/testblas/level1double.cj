package cjnum.blas.testblas

import std.math.*
import cjnum.blas.*
import cjnum.floats.*
import cjnum.floats.scalar.*

const pi = 3.14159265358979323846264338327950288419716939937510582097494459

public class DoubleOneVectorCase {
    DoubleOneVectorCase(
        public let name!: String,
        public let x!: Array<Float64>,
        public let incx!: Int64,
        public let n!: Int64,
        public let exception!: Bool,
        public let dasum!: Float64 = 0.0,
        public let dnrm2!: Float64 = 0.0,
        public let idamax!: Int64 = 0,
        public let dscalCases!: Array<DScalCase>
    ) {}
}

public class DScalCase {
    DScalCase(
        public let alpha!: Float64,
        public let ans!: Array<Float64>,
        public let name!: String = ""
    ) {}
}

public let doubleOneVectorCases = [
    DoubleOneVectorCase(
        name: "AllPositive",
        x: [6.0, 5.0, 4.0, 2.0, 6.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: 23.0,
        dnrm2: 10.81665382639196787935766380241148783875388972153573863813135,
        idamax: 0,
        dscalCases: [
            DScalCase(
                alpha: 0.0,
                ans: [0.0, 0.0, 0.0, 0.0, 0.0],
            ),
            DScalCase(
                alpha: 1.0,
                ans: [6.0, 5.0, 4.0, 2.0, 6.0],
            ),
            DScalCase(
                alpha: -2.0,
                ans: [-12.0, -10.0, -8.0, -4.0, -12.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "LeadingZero",
        x: [0.0, 1.0],
        incx: 1,
        n: 2,
        exception: false,
        dasum: 1.0,
        dnrm2: 1.0,
        idamax: 1,
        dscalCases: [
            DScalCase(
                alpha: 0.0,
                ans: [0.0, 0.0],
            ),
            DScalCase(
                alpha: 1.0,
                ans: [0.0, 1.0]
            ),
            DScalCase(
                alpha: -2.0,
                ans: [0.0, -2.0]
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "MaxInMiddle",
        x: [6.0, 5.0, 9.0, 0.0, 6.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: 26.0,
        dnrm2: 13.34166406412633371248943627250846646911846482744007727141318,
        idamax: 2,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-12.0, -10.0, -18.0, 0.0, -12.0]
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "MaxAtEnd",
        x: [6.0, 5.0, -9.0, 0.0, 10.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: 30.0,
        dnrm2: 15.55634918610404553681857596630667886426639062914642880494347,
        idamax: 4,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-12.0, -10.0, 18.0, 0.0, -20.0]
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "AllNegative",
        x: [-6.0, -5.0, -4.0, -2.0, -6.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: 23.0,
        dnrm2: 10.81665382639196787935766380241148783875388972153573863813135,
        idamax: 0,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [12.0, 10.0, 8.0, 4.0, 12.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "AllMixed",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: 23.0,
        dnrm2: 10.81665382639196787935766380241148783875388972153573863813135,
        idamax: 0,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [12.0, -10.0, -8.0, 4.0, 12.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "ZeroN",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: 1,
        n: 0,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "OneN",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: 1,
        n: 1,
        exception: false,
        dasum: 6.0,
        dnrm2: 6.0,
        idamax: 0,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [12.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "PositiveExactInc",
        x: [-6.0, 5.0, 10.0, -2.0, -5.0],
        incx: 2,
        n: 3,
        exception: false,
        dasum: 21.0,
        dnrm2: 12.68857754044952038019377274608948979173952662752515253090272,
        idamax: 1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [12.0, 5.0, -20.0, -2.0, 10.0]
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "PositiveOffInc",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0, 8.0, 10.0, 11.0],
        incx: 3,
        n: 3,
        exception: false,
        dasum: 18.0,
        dnrm2: 11.83215956619923208513465658312323409683100246158868064575943,
        idamax: 2,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [12.0, 5.0, 4.0, 4.0, -6.0, 8.0, -20.0, 11.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "PositiveShortInc",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0, 8.0, 10.0, 11.0],
        incx: 3,
        n: 2,
        exception: false,
        dasum: 8.0,
        dnrm2: 6.324555320336758663997787088865437067439110278650433653715009,
        idamax: 0,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [12.0, 5.0, 4.0, 4.0, -6.0, 8.0, 10.0, 11.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NegativeInc",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: -1,
        n: 5,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NegativeExactInc",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: -2,
        n: 3,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NegativeOffInc",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0, 8.0, 10.0, 11.0],
        incx: -3,
        n: 2,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0, 8.0, 10.0, 11.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NegativeShortInc",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0, 8.0, 10.0, 11.0],
        incx: -3,
        n: 2,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0, 8.0, 10.0, 11.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NegativeN",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: 2,
        n: -5,
        exception: true,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "ZeroInc",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: 0,
        n: 5,
        exception: true,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "OutOfBounds",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: 2,
        n: 6,
        exception: true,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NegativeOutOfBounds",
        x: [-6.0, 5.0, 4.0, -2.0, -6.0],
        incx: -2,
        n: 6,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-6.0, 5.0, 4.0, -2.0, -6.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NaN",
        x: [Float64.NaN, 2.0],
        incx: 1,
        n: 2,
        exception: false,
        dasum: Float64.NaN,
        dnrm2: Float64.NaN,
        idamax: 0,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [Float64.NaN, -4.0],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [0.0, 0.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NaNInc",
        x: [Float64.NaN, Float64.NaN, 2.0],
        incx: 2,
        n: 2,
        exception: false,
        dasum: Float64.NaN,
        dnrm2: Float64.NaN,
        idamax: 0,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [Float64.NaN, Float64.NaN, -4.0],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [0.0, Float64.NaN, 0.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "Empty",
        x: [],
        incx: 1,
        n: 0,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "EmptyZeroInc",
        x: [],
        incx: 0,
        n: 0,
        exception: true,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "EmptyReverse",
        x: [],
        incx: -1,
        n: 0,
        exception: false,
        dasum: 0.0,
        dnrm2: 0.0,
        idamax: -1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "MultiInf",
        x: [5.0, Float64.Inf, -Float64.Inf, 8.0, 9.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: Float64.Inf,
        dnrm2: Float64.Inf,
        idamax: 1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-10.0, -Float64.Inf, Float64.Inf, -16.0, -18.0],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [0.0, 0.0, 0.0, 0.0, 0.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "NaNInf",
        x: [5.0, Float64.NaN, -Float64.Inf, 8.0, 9.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: Float64.NaN,
        dnrm2: Float64.NaN,
        idamax: 2,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-10.0, Float64.NaN, Float64.Inf, -16.0, -18.0],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [0.0, 0.0, 0.0, 0.0, 0.0],
            )
        ]
    ),
    DoubleOneVectorCase(
        name: "InfNaN",
        x: [5.0, Float64.Inf, Float64.NaN, 8.0, 9.0],
        incx: 1,
        n: 5,
        exception: false,
        dasum: Float64.NaN,
        dnrm2: Float64.NaN,
        idamax: 1,
        dscalCases: [
            DScalCase(
                alpha: -2.0,
                ans: [-10.0, -Float64.Inf, Float64.NaN, -16.0, -18.0],
            ),
            DScalCase(
                alpha: 0.0,
                ans: [0.0, 0.0, 0.0, 0.0, 0.0],
            )
        ]
    )
]

public class DoubleTwoVectorCase {
    DoubleTwoVectorCase(
        public let name!: String,
        public let x!: Array<Float64>,
        public let y!: Array<Float64>,
        public var xTmp!: Array<Float64> = [],
        public let yTmp!: Array<Float64> = [],
        public let incx!: Int64,
        public let incy!: Int64,
        public let n!: Int64,
        public let exception!: Bool,
        public let daxpyCases!: Array<DaxpyCase>,
        public let ddotAns!: Float64 = 0.0,
        public let dswapAns!: DTwoVecAnswer = DTwoVecAnswer(),
        public let dcopyAns!: DTwoVecAnswer = DTwoVecAnswer(),
        public let drotCases!: Array<DrotCase>,
        public let drotmCases!: Array<DrotmCase>
    ) {}
}

public class DaxpyCase {
    DaxpyCase(
        public let alpha!: Float64,
        public let ans!: Array<Float64>
    ) {}
}

public class DrotCase {
    DrotCase(
        public let c!: Float64,
        public let s!: Float64,
        public let xAns!: Array<Float64>,
        public let yAns!: Array<Float64>
    ) {}
}

public class DrotmCase {
    DrotmCase(
        public let p!: DrotmParams,
        public let xAns!: Array<Float64>,
        public let yAns!: Array<Float64>,
        public let name!: String = ""
    ) {}
}

public class DTwoVecAnswer {
    DTwoVecAnswer(
        public let x!: Array<Float64> = [],
        public let y!: Array<Float64> = []
    ) {}
}

public let doubleTwoVectorCases = [
    DoubleTwoVectorCase(
        name: "UnitaryInc",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 1,
        incy: 1,
        n: 6,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 1.0,
                ans: [18.0, 13.0, -2.0, 10.0, 20.0, 4.0],
            ),
            DaxpyCase(
                alpha: 2.0,
                ans: [28.0, 28.0, -8.0, 13.0, 34.0, 11.0],
            ),
            DaxpyCase(
                alpha: -3.0,
                ans: [-22.0, -47.0, 22.0, -2.0, -36.0, -24.0],
            ),
            DaxpyCase(
                alpha: 0.0,
                ans: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0],
            )
        ],
        ddotAns: 110.0,
        dswapAns: DTwoVecAnswer(
            x: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0],
            y: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        ),
        dcopyAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
            y: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        ),
        drotCases: [
            DrotCase(
                c: cos(0.0),
                s: sin(0.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0],
            ),
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [12.444023964292095, 12.749380282068351, -3.7473736752571014, 5.677251193294846, 15.224018588957296, 5.076299724034451],
                yAns: [3.024279678886205, -8.151889500183792, 6.160940718590796, 5.076299724034451, -0.4788089421498931, -5.677251193294846],
            ),
            DrotCase(
                c: cos(0.5 * pi),
                s: sin(0.5 * pi),
                xAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0],
                yAns: [-10.0, -15.0, 6.0, -3.0, -14.0, -7.0],
            ),
            DrotCase(
                c: cos(pi),
                s: sin(pi),
                xAns: [-10.0, -15.0, 6.0, -3.0, -14.0, -7.0],
                yAns: [-8.0, 2.0, -4.0, -7.0, -6.0, 3.0],
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Identity,
                    h: [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0],
                name: "Neg2Flag"
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h: [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [8.2, 13.7, -5.8, 2.0, 12.0, 6.6],
                yAns: [5.0, 0.5, 1.4, 3.8, 4.4, -0.8],
                name: "Neg1Flag"
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: OffDiagonal,
                    h: [1.0, 0.1, -0.1, 1.0]
                ),
                xAns: [9.2, 15.2, -6.4, 2.3, 13.4, 7.3],
                yAns: [9.0, -0.5, 3.4, 7.3, 7.4, -2.3],
                name: "ZeroFlag"
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: Diagonal,
                    h: [0.5, -1.0, 1.0, 0.7]
                ),
                xAns: [13.0, 5.5, 1.0, 8.5, 13.0, 0.5],
                yAns: [-4.4, -16.4, 8.8, 1.9, -9.8, -9.1],
                name: "OneFlag"
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "UnitaryIncLong",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0, 8.0, -9.0, 10.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, 7.0, -6.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 1,
        incy: 1,
        n: 6,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 1.0,
                ans: [18.0, 13.0, -2.0, 10.0, 20.0, 4.0, 7.0, -6.0]
            ),
            DaxpyCase(
                alpha: 2.0,
                ans: [28.0, 28.0, -8.0, 13.0, 34.0, 11.0, 7.0, -6.0]
            ),
            DaxpyCase(
                alpha: -3.0,
                ans: [-22.0, -47.0, 22.0, -2.0, -36.0, -24.0, 7.0, -6.0]
            ),
            DaxpyCase(
                alpha: 0.0,
                ans: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, 7.0, -6.0]
            )
        ],
        ddotAns: 110.0,
        dswapAns: DTwoVecAnswer(
            x: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, 8.0, -9.0, 10.0],
            y: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0, 7.0, -6.0]
		),
        dcopyAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0, 8.0, -9.0, 10.0],
            y: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0, 7.0, -6.0]
		),
        drotCases: [
            DrotCase(
                c: cos(0.0),
                s: sin(0.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0, 8.0, -9.0, 10.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, 7.0, -6.0]
            ),
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [12.444023964292095, 12.749380282068351, -3.7473736752571014, 5.677251193294846, 15.224018588957296, 5.076299724034451, 8.0, -9.0, 10.0],
                yAns: [3.024279678886205, -8.151889500183792, 6.160940718590796, 5.076299724034451, -0.4788089421498931, -5.677251193294846, 7.0, -6.0]
            ),
            DrotCase(
                c: cos(0.5 * pi),
                s: sin(0.5 * pi),
                xAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, 8.0, -9.0, 10.0],
                yAns: [-10.0, -15.0, 6.0, -3.0, -14.0, -7.0, 7.0, -6.0]
            ),
            DrotCase(
                c: cos(pi),
                s: sin(pi),
                xAns: [-10.0, -15.0, 6.0, -3.0, -14.0, -7.0, 8.0, -9.0, 10.0],
                yAns: [-8.0, 2.0, -4.0, -7.0, -6.0, 3.0, 7.0, -6.0]
            )
        ],
		drotmCases: [
			DrotmCase(
				p: DrotmParams(
					flag: Identity,
					h: [0.9, 0.1, -0.1, 0.5]
				),
				xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0, 8.0, -9.0, 10.0],
				yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, 7.0, -6.0],
				name: "Neg2Flag"
			),
			DrotmCase(
				p: DrotmParams(
					flag: Rescaling,
					h: [0.9, 0.1, -0.1, 0.5]
				),
				xAns: [8.2, 13.7, -5.8, 2.0, 12.0, 6.6, 8.0, -9.0, 10.0],
				yAns: [5.0, 0.5, 1.4, 3.8, 4.4, -0.8, 7.0, -6.0],
				name: "Neg1Flag"
			),
			DrotmCase(
				p: DrotmParams(
					flag: OffDiagonal,
					h: [1.0, 0.1, -0.1, 1.0]
				),
				xAns: [9.2, 15.2, -6.4, 2.3, 13.4, 7.3, 8.0, -9.0, 10.0],
				yAns: [9.0, -0.5, 3.4, 7.3, 7.4, -2.3, 7.0, -6.0],
				name: "ZeroFlag"
			),
			DrotmCase(
				p: DrotmParams(
					flag: Diagonal,
					h: [0.5, -1.0, 1.0, 0.7]
				),
				xAns: [13.0, 5.5, 1.0, 8.5, 13.0, 0.5, 8.0, -9.0, 10.0],
				yAns: [-4.4, -16.4, 8.8, 1.9, -9.8, -9.1, 7.0, -6.0],
				name: "OneFlag"
			)
		]
    ),
    DoubleTwoVectorCase(
        name: "PositiveInc",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 2,
        incy: 3,
        n: 3,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [28.0, -2.0, 4.0, -5.0, 6.0, -3.0, 24.0, 10.0]
            )
        ],
        ddotAns: -18.0,
        dswapAns: DTwoVecAnswer(
            x: [8.0, 15.0, 7.0, 3.0, -4.0, 7.0],
            y: [10.0, -2.0, 4.0, -6.0, 6.0, -3.0, 14.0, 10.0]
        ),
        dcopyAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
            y: [10.0, -2.0, 4.0, -6.0, 6.0, -3.0, 14.0, 10.0]
        ),
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [12.444023964292095, 15.0, -2.479518890035003, 3.0, 10.997835971550302, 7.0],
                yAns: [3.024279678886205, -2.0, 4.0, 8.879864079700745, 6.0, -3.0, -9.541886812516392, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h: [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [8.2, 15.0, -6.1, 3.0, 13.0, 7.0],
                yAns: [5.0, -2.0, 4.0, 2.9, 6.0, -3.0, -0.6, 10.0]
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: OffDiagonal,
                    h: [1.0, 0.1, -0.1, 1.0]
                ),
                xAns: [9.2, 15.0, -6.7, 3.0, 14.4, 7.0],
                yAns: [9.0, -2.0, 4.0, 6.4, 6.0, -3.0, -2.6, 10.0],
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: Diagonal,
                    h: [0.5, -1.0, 1.0, 0.7]
                ),
                xAns: [13.0, 15.0, 4.0, 3.0, 3.0, 7.0],
                yAns: [-4.4, -2.0, 4.0, 10.9, 6.0, -3.0, -16.8, 10.0]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "NegativeInc",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: -2,
        incy: -3,
        n: 3,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [28.0, -2.0, 4.0, -5.0, 6.0, -3.0, 24.0, 10.0],
            )
        ],
        ddotAns: -18.0,
        dswapAns: DTwoVecAnswer(
            x: [8.0, 15.0, 7.0, 3.0, -4.0, 7.0],
            y: [10.0, -2.0, 4.0, -6.0, 6.0, -3.0, 14.0, 10.0],
        ),
        dcopyAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
            y: [10.0, -2.0, 4.0, -6.0, 6.0, -3.0, 14.0, 10.0],
        ),
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [12.444023964292095, 15.0, -2.479518890035003, 3.0, 10.997835971550302, 7.0],
                yAns: [3.024279678886205, -2.0, 4.0, 8.879864079700745, 6.0, -3.0, -9.541886812516392, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:  [0.9, 0.1, -0.1, 0.5],
                ),
                xAns: [8.2, 15.0, -6.1, 3.0, 13.0, 7.0],
                yAns: [5.0, -2.0, 4.0, 2.9, 6.0, -3.0, -0.6, 10.0]
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: OffDiagonal,
                    h:  [1.0, 0.1, -0.1, 1.0],
                ),
                xAns: [9.2, 15.0, -6.7, 3.0, 14.4, 7.0],
                yAns: [9.0, -2.0, 4.0, 6.4, 6.0, -3.0, -2.6, 10.0]
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: Diagonal,
                    h:  [0.5, -1.0, 1.0, 0.7],
                ),
                xAns: [13.0, 15.0, 4.0, 3.0, 3.0, 7.0],
                yAns: [-4.4, -2.0, 4.0, 10.9, 6.0, -3.0, -16.8, 10.0]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "MixedInc1",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 2,
        incy: -3,
        n: 3,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [36.0, -2.0, 4.0, -5.0, 6.0, -3.0, 16.0, 10.0],
            )
        ],
        ddotAns: 30.0,
        dswapAns: DTwoVecAnswer(
            x: [-4.0, 15.0, 7.0, 3.0, 8.0, 7.0],
            y: [14.0, -2.0, 4.0, -6.0, 6.0, -3.0, 10.0, 10.0],
        ),
        dcopyAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
            y: [14.0, -2.0, 4.0, -6.0, 6.0, -3.0, 10.0, 10.0],
        ),
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [7.372604823403701, 15.0, -2.479518890035003, 3.0, 16.069255112438693, 7.0],
                yAns: [1.333806631923407, -2.0, 4.0, 8.879864079700745, 6.0, -3.0, -7.851413765553595, 10.0],
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [9.4, 15.0, -6.1, 3.0, 11.8, 7.0],
                yAns: [5.4, -2.0, 4.0, 2.9, 6.0, -3.0, -1.0, 10.0],
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: OffDiagonal,
                    h:    [1.0, 0.1, -0.1, 1.0]
                ),
                xAns: [10.4, 15.0, -6.7, 3.0, 13.2, 7.0],
                yAns: [9.4, -2.0, 4.0, 6.4, 6.0, -3.0, -3.0, 10.0]
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: Diagonal,
                    h:    [0.5, -1.0, 1.0, 0.7]
                ),
                xAns: [1.0, 15.0, 4.0, 3.0, 15.0, 7.0],
                yAns: [-8.4, -2.0, 4.0, 10.9, 6.0, -3.0, -12.8, 10.0]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "MixedInc2",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: -2,
        incy: 3,
        n: 3,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [36.0, -2.0, 4.0, -5.0, 6.0, -3.0, 16.0, 10.0],
            )
        ],
        ddotAns: 30.0,
        dswapAns: DTwoVecAnswer(
            x: [-4.0, 15.0, 7.0, 3.0, 8.0, 7.0],
            y: [14.0, -2.0, 4.0, -6.0, 6.0, -3.0, 10.0, 10.0]
        ),
        dcopyAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
            y: [14.0, -2.0, 4.0, -6.0, 6.0, -3.0, 10.0, 10.0]
        ),
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [7.372604823403701, 15.0, -2.479518890035003, 3.0, 16.069255112438693, 7.0],
                yAns: [1.333806631923407, -2.0, 4.0, 8.879864079700745, 6.0, -3.0, -7.851413765553595, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [9.4, 15.0, -6.1, 3.0, 11.8, 7.0],
                yAns: [5.4, -2.0, 4.0, 2.9, 6.0, -3.0, -1.0, 10.0],
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: OffDiagonal,
                    h:    [1.0, 0.1, -0.1, 1.0]
                ),
                xAns: [10.4, 15.0, -6.7, 3.0, 13.2, 7.0],
                yAns: [9.4, -2.0, 4.0, 6.4, 6.0, -3.0, -3.0, 10.0]
            ),
            DrotmCase(
                p: DrotmParams(
                    flag: Diagonal,
                    h:    [0.5, -1.0, 1.0, 0.7]
                ),
                xAns: [1.0, 15.0, 4.0, 3.0, 15.0, 7.0],
                yAns: [-8.4, -2.0, 4.0, 10.9, 6.0, -3.0, -12.8, 10.0]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "ZeroN",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: -2,
        incy: 3,
        n: 0,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ],
        dswapAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
            y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
        ),
        dcopyAns: DTwoVecAnswer(
            x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
            y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
        ),
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "NegativeN",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: -2,
        incy: 3,
        n: -3,
        exception: true,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [36.0, -2.0, 4.0, -5.0, 6.0, -3.0, 16.0, 10.0]
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [8.2, 13.7, -5.8, 2.0, 12.0, 6.6],
                yAns: [5.0, 0.5, 1.4, 3.8, 4.4, -0.8]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "ZeroIncX",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 0,
        incy: 3,
        n: 2,
        exception: true,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [36.0, -2.0, 4.0, -5.0, 6.0, -3.0, 16.0, 10.0]
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [8.2, 13.7, -5.8, 2.0, 12.0, 6.6],
                yAns: [5.0, 0.5, 1.4, 3.8, 4.4, -0.8]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "ZeroIncY",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 1,
        incy: 0,
        n: 2,
        exception: true,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [36.0, -2.0, 4.0, -5.0, 6.0, -3.0, 16.0, 10.0]
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [8.2, 13.7, -5.8, 2.0, 12.0, 6.6],
                yAns: [5.0, 0.5, 1.4, 3.8, 4.4, -0.8]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "OutOfBoundsX",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 8,
        incy: 2,
        n: 2,
        exception: true,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [36.0, -2.0, 4.0, -5.0, 6.0, -3.0, 16.0, 10.0]
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [8.2, 13.7, -5.8, 2.0, 12.0, 6.6],
                yAns: [5.0, 0.5, 1.4, 3.8, 4.4, -0.8]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "OutOfBoundsY",
        x: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
        y: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0],
        xTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        yTmp: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        incx: 2,
        incy: 8,
        n: 2,
        exception: true,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: [36.0, -2.0, 4.0, -5.0, 6.0, -3.0, 16.0, 10.0]
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [10.0, 15.0, -6.0, 3.0, 14.0, 7.0],
                yAns: [8.0, -2.0, 4.0, 7.0, 6.0, -3.0, -4.0, 10.0]
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "Empty",
        x: [],
        y: [],
        incx: 1,
        incy: 1,
        n: 0,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: []
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [],
                yAns: []
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [],
                yAns: []
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "EmptyZeroIncX",
        x: [],
        y: [],
        incx: 0,
        incy: 1,
        n: 0,
        exception: true,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: []
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [],
                yAns: []
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [],
                yAns: []
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "EmptyZeroIncY",
        x: [],
        y: [],
        incx: 1,
        incy: 0,
        n: 0,
        exception: true,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: []
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [],
                yAns: []
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [],
                yAns: []
            )
        ]
    ),
    DoubleTwoVectorCase(
        name: "EmptyReverse",
        x: [],
        y: [],
        incx: -1,
        incy: -1,
        n: 0,
        exception: false,
        daxpyCases: [
            DaxpyCase(
                alpha: 2.0,
                ans: []
            )
        ],
        drotCases: [
            DrotCase(
                c: cos(25.0 * pi / 180.0),
                s: sin(25.0 * pi / 180.0),
                xAns: [],
                yAns: []
            )
        ],
        drotmCases: [
            DrotmCase(
                p: DrotmParams(
                    flag: Rescaling,
                    h:    [0.9, 0.1, -0.1, 0.5]
                ),
                xAns: [],
                yAns: []
            )
        ]
    )
]

public func ddotTest(nFloat64: NFloat64): Bool {
    for (c in doubleTwoVectorCases) {
        dCopyTwoTmp(c.x, c.xTmp, c.y, c.yTmp)
        if (c.exception) {
            let f = { =>
                nFloat64.ddot(c.n, c.xTmp, c.incx, c.yTmp, c.incy)
                return 
            }
            if (!exception(f)) {
                throw Exception("${c.name} should panic and does not")
            }
            continue
        }
        let dot = nFloat64.ddot(c.n, c.xTmp, c.incx, c.yTmp, c.incy)
        if (!dTolEqual(dot, c.ddotAns)) {
            throw Exception("ddot: mismatch ${c.name}: expected ${c.ddotAns}, found ${dot}")
        }
    }
    return true
}

public func dnrm2Test(nFloat64: NFloat64): Bool {
    for (c in doubleOneVectorCases) {
        if (c.exception) {
            let f = { => 
                nFloat64.dnrm2(c.n, c.x, c.incx)
                return 
            }
            if (!exception(f)) {
                throw Exception("${c.name} should panic and does not")
            }
            continue
        }
        let v = nFloat64.dnrm2(c.n, c.x, c.incx)
        if (!dTolEqual(v, c.dnrm2)) {
            throw Exception("dnrm2: mismatch ${c.name}: expected ${c.dnrm2}, found ${v}")
        }
    }
    return true
}

public func dasumTest(nFloat64: NFloat64): Bool {
    for (c in doubleOneVectorCases) {
        if (c.exception) {
            let f = { =>
                nFloat64.dasum(c.n, c.x, c.incx)
                return 
            }
            if (!exception(f)) {
                throw Exception("${c.name} should panic and does not")
                return false
            }
            continue
        }
        let v = nFloat64.dasum(c.n, c.x, c.incx)
        if (!dTolEqual(v, c.dasum)) {
            throw Exception("dasum: mismatch ${c.name}: expected ${c.dnrm2}, found ${v}")
            return false
        }
    }
    return true
}

public func idamaxTest(nFloat64: NFloat64): Bool {
    for (c in doubleOneVectorCases) {
        if (c.exception) {
            let f = { =>
                nFloat64.idamax(c.n, c.x, c.incx)
                return 
            }
            if (!exception(f)) {
                throw Exception("${c.name} should panic and does not")
                return false
            }
            continue
        }
        let v = nFloat64.idamax(c.n, c.x, c.incx)
        if (v != c.idamax) {
            if (fHasNaN(c.x)) {
                println("idamax: mismatch ${c.name}: expected ${c.idamax}, found ${v}")
            } else {
                throw Exception("idamax: mismatch ${c.name}: expected ${c.idamax}, found ${v}")
            }
        }
    }
    return true
}

public func dswapTest(nFloat64: NFloat64): Bool {
    for (c in doubleTwoVectorCases) {
        dCopyTwoTmp(c.x, c.xTmp, c.y, c.yTmp)
        if (c.exception) {
            let f = { => 
                nFloat64.dswap(c.n, c.xTmp, c.incx, c.yTmp, c.incy)
                return 
            }
            if (!exception(f)) {
                throw Exception("${c.name} should panic and does not")
            }
            continue
        }
        nFloat64.dswap(c.n, c.xTmp, c.incx, c.yTmp, c.incy)
        if (!dSliceTolEqual(c.xTmp, c.dswapAns.x)) {
            throw Exception("dswap: mismatch ${c.name} (x): expected ${c.dswapAns.x}, found ${c.xTmp}")
        }
        if (!dSliceTolEqual(c.yTmp, c.dswapAns.y)) {
            throw Exception("dswap: mismatch ${c.name} (y): expected ${c.dswapAns.y}, found ${c.yTmp}")
        }
    }
    return true
}

public func dcopyTest(nFloat64: NFloat64): Bool {
    for (c in doubleTwoVectorCases) {
        dCopyTwoTmp(c.x, c.xTmp, c.y, c.yTmp)
        if (c.exception) {
            let f = { => 
                nFloat64.dcopy(c.n, c.xTmp, c.incx, c.yTmp, c.incy)
            }
            if (!exception(f)) {
                throw Exception("${c.name} should panic and does not")
            }
            continue
        }
        nFloat64.dcopy(c.n, c.xTmp, c.incx, c.yTmp, c.incy)
        if (!dSliceTolEqual(c.xTmp, c.dcopyAns.x)) {
            throw Exception("dcopy: mismatch ${c.name} (x): expected ${c.dcopyAns.x}, found ${c.xTmp}")
            return false
        }
        if (!dSliceTolEqual(c.yTmp, c.dcopyAns.y)) {
            throw Exception("dcopy: mismatch ${c.name} (y): expected ${c.dcopyAns.y}, found ${c.yTmp}")
            return false
        }
    }
    return true
}

public func daxpyTest(nFloat64: NFloat64): Bool {
    for (c in doubleTwoVectorCases) {
        for (kind in c.daxpyCases) {
            dCopyTwoTmp(c.x, c.xTmp, c.y, c.yTmp)
            if (c.exception) {
                let f = { => 
                    nFloat64.daxpy(c.n, kind.alpha, c.xTmp, c.incx, c.yTmp, c.incy)
                    return 
                }
                if (!exception(f)) {
                    throw Exception("${c.name} should panic and does not")
                }
                continue
            }
            nFloat64.daxpy(c.n, kind.alpha, c.xTmp, c.incx, c.yTmp, c.incy)
            if (!dSliceTolEqual(c.yTmp, kind.ans)) {
                throw Exception("daxpy: mismatch ${c.name}: expected ${kind.ans}, found ${c.yTmp}")
            }
        }
    }
    return true
}

public class DrotgTestStruct {
    DrotgTestStruct(
        public let name!: String,
        public let a!: Float64 = 0.0,
        public let b!: Float64 = 0.0,
        public let c!: Float64 = 0.0,
        public let s!: Float64 = 0.0,
        public let r!: Float64 = 0.0,
        public let z!: Float64 = 0.0
    ) {}
}

public let drotgTests = [
    DrotgTestStruct(
        name: "ZeroAB",
        c: 1.0
    ),
    DrotgTestStruct(
        name: "PosA_ZeroB",
        a: 0.5,
        c: 1.0,
        r: 0.5
    ),
    DrotgTestStruct(
        name: "NegA_ZeroB",
        a: -4.6,
        c: 1.0,
        r: -4.6
    ),
    DrotgTestStruct(
        name: "ZeroA_PosB",
        b: 3.0,
        s: 1.0,
        r: 3.0,
        z: 1.0
    ),
    DrotgTestStruct(
        name: "ZeroA_NegB",
        b: -0.3,
        s: 1.0,
        r: -0.3,
        z: 1.0
    ),
    DrotgTestStruct(
		name: "PosA_PosB_AGTB",
		a:    5.0,
		b:    0.3,
		c:    0.99820484546577868593549038000,
		s:    0.05989229072794672115612942280,
		r:    5.00899191454727744602429072688,
		z:    0.05989229072794672115612942280
    ),
    DrotgTestStruct(
		name: "PosA_PosB_ALTB",
		a:    3.0,
		b:    4.0,
		c:    3.0 / 5.0,
		s:    4.0 / 5.0,
		r:    5.0,
		z:    5.0 / 3.0
    ),
    DrotgTestStruct(
		name: "PosA_NegB_AGTB",
		a:    2.6,
		b:    -0.9,
		c:    0.94498607344025815971847507095,
		s:    -0.32711056388316628605639521686,
		r:    2.751363298439520872718790879655,
		z:    -0.3271105638831662860563952168
    ),
    DrotgTestStruct(
		name: "PosA_NegB_ALTB",
		a:    2.6,
		b:    -2.9,
		c:    -0.6675450157520258540548049558,
		s:    0.7445694406464903756765132200,
		r:    -3.8948684188300893100043812234,
		z:    1.0 / -0.6675450157520258540548049558
    ),
    DrotgTestStruct(
		name: "NegA_PosB_AGTB",
		a:    -11.4,
		b:    10.3,
		c:    0.7419981952497362418487847947,
		s:    -0.6704018781642353764072353847,
		r:    -15.363918770938617534070671122,
		z:    -0.6704018781642353764072353847
    ),
    DrotgTestStruct(
		name: "NegA_PosB_ALTB",
		a:    -1.4,
		b:    10.3,
		c:    -0.1346838895922121112404717523,
		s:    0.9908886162855605326977564640,
		r:    10.394710193170370442523552032,
		z:    1.0 / -0.1346838895922121112404717523
    ),
    DrotgTestStruct(
		name: "NegA_NegB_AGTB",
		a:    -11.4,
		b:    10.3,
		c:    0.7419981952497362418487847947,
		s:    -0.6704018781642353764072353847,
		r:    -15.363918770938617534070671122,
		z:    -0.6704018781642353764072353847,
    ),
    DrotgTestStruct(
		name: "NegA_NegB_ALTB",
		a:    -1.4,
		b:    -10.3,
		c:    0.1346838895922121112404717523,
		s:    0.9908886162855605326977564640,
		r:    -10.394710193170370442523552032,
		z:    1.0 / 0.1346838895922121112404717523,
    )
]

public func drotgTest(nFloat64: NFloat64, skipExtreme: Bool): Bool {
    for (test in drotgTests) {
        let (c, s, r, z) = nFloat64.drotg(test.a, test.b)
        if (!dTolEqual(c, test.c)) {
            throw Exception("drotg: c mismatch ${test.name}: expected ${test.c}, found ${c}")
        }
        if (!dTolEqual(s, test.s)) {
            throw Exception("drotg: s mismatch ${test.name}: expected ${test.s}, found ${s}")
        }
        if (!dTolEqual(r, test.r)) {
            throw Exception("drotg: r mismatch ${test.name}: expected ${test.r}, found ${r}")
        }
        if (!dTolEqual(z, test.z)) {
            throw Exception("drotg: z mismatch ${test.name}: expected ${test.z}, found ${z}")
        }
    }

	const ulp    = 0x1p-52
	const safmin = 0x1p-1022
	const safmax = 1.0 / safmin
	const tol    = 20.0 * ulp

    var values = [
		-1.0 / ulp,
		-1.0,
		-1.0 / 3.0,
		-ulp,
		0.0,
		ulp,
		1.0 / 3.0,
		1.0,
		1.0 / ulp,
		-Float64.Inf,
		Float64.Inf,
		Float64.NaN
    ]
    if (!skipExtreme) {
        values = [
            -1.0 / ulp,
            -1.0,
            -1.0 / 3.0,
            -ulp,
            0.0,
            ulp,
            1.0 / 3.0,
            1.0,
            1.0 / ulp,
            -Float64.Inf,
            Float64.Inf,
            Float64.NaN,
            -safmax,
            -safmin,
            safmin,
            safmax
        ]
    }

    for (f in values) {
        for (g in values) {
            let name = "Case f=${f}, g=${g}"

            let (cs, sn, r, _) = nFloat64.drotg(f, g)
            
            if (f.isNaN() || g.isNaN()) {
                if (!r.isNaN()) {
                    throw Exception("${name}: unexpected r=${r}; want NaN")
                }
            } else if (f.isInf() || g.isInf()) {
                if (!r.isNaN() && !r.isInf()) {
                    throw Exception("${name}: unexpected r=${r}; want NaN or Inf")
                }
            } else {
                var d = max(abs(f), abs(g))
                d = min(max(safmin, d), safmax)
                let fs = f / d
                let gs = g / d
                let rs = r / d

                var rnorm = abs(rs)
                if (rnorm == 0.0) {
                    rnorm = max(abs(fs), abs(gs))
                    if (rnorm == 0.0) {
                        rnorm = 1.0
                    }
                }
                var resid = abs(rs - (cs * fs + sn * gs)) / rnorm
                if (resid > tol) {
                    throw Exception("${name}: cs*f + sn*g != r; resid=${resid}")
                }

                resid = abs(-sn * fs + cs * gs)
                if (resid > tol) {
                    throw Exception("${name}: cs*f + sn*g != 0; resid=${resid}")
                }

                resid = abs(1.0 - (cs * cs + sn * sn))
                if (resid > tol) {
                    throw Exception("${name}: cs*f + sn*g != 1; resid=${resid}")
                }

                if (abs(f) > abs(g)) {
                    if (cs < 0.0) {
                        throw Exception("${name}: cs is negative; cs=${cs}")
                    }
                } else {
                    if (cs * copysign(1.0, f) * copysign(1.0, g) < 0.0) {
                        throw Exception("${name}: sign of cs doesn't match sign of f and g; cs=${cs}, sign(f)=${copysign(1.0, f)}, sign(g)=${copysign(1.0, g)}")
                    }
                }
            }
        }
    }
    return true
}

func copysign(f: Float64, sign: Float64): Float64 {
    const signBit: UInt64 = 1 << 63
    return Float64.fromBits((f.toBits() & (!signBit)) | sign.toBits() & signBit)
}

public class drotmgTestStruct {
    drotmgTestStruct(
        let name!: String,
        let d1!: Float64,
        let d2!: Float64,
        let x1!: Float64,
        let y1!: Float64,
        let d1Want!: Float64,
        let d2Want!: Float64,
        let x1Want!: Float64,
        let pWant!: DrotmParams
    ) {}
}

public func drotmgTest(nFloat64: NFloat64): Bool {
    const tol = 1e-14

    let tests = [
        drotmgTestStruct(
			name:   "NegD1",
			d1:     -4.0,
			d2:     6.0,
			x1:     8.0,
			y1:     -4.0,
			d1Want: 0.0,
			d2Want: 0.0,
			x1Want: 0.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [0.0, 0.0, 0.0, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "ZeroD2",
			d1:     4.0,
			d2:     0.0,
			x1:     8.0,
			y1:     -5.0,
			d1Want: 4.0,
			d2Want: 0.0,
			x1Want: 8.0,
			pWant: DrotmParams(
				flag: Identity,
				h:    [0.0, 0.0, 0.0, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "ZeroY1",
			d1:     4.0,
			d2:     -6.0,
			x1:     8.0,
			y1:     0.0,
			d1Want: 4.0,
			d2Want: -6.0,
			x1Want: 8.0,
			pWant: DrotmParams(
				flag: Identity,
				h:    [0.0, 0.0, 0.0, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "NegQ2_and_AQ1_LT_AQ2",
			d1:     8.0,
			d2:     -6.0,
			x1:     4.0,
			y1:     8.0,
			d1Want: 0.0,
			d2Want: 0.0,
			x1Want: 0.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [0.0, 0.0, 0.0, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "ZeroD1",
			d1:     0.0,
			d2:     2.0,
			x1:     8.0,
			y1:     4.0,
			d1Want: 2.0,
			d2Want: 0.0,
			x1Want: 4.0,
			pWant: DrotmParams(
				flag: Diagonal,
				h:    [0.0, 0.0, 0.0, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "AbsQ1_GT_AbsQU__D2_Pos",
			d1:     2.0,
			d2:     3.0,
			x1:     8.0,
			y1:     5.0,
			d1Want: 1.2610837438423645,
			d2Want: 1.8916256157635467,
			x1Want: 12.6875,
			pWant: DrotmParams(
				flag: OffDiagonal,
				h:    [0.0, -0.625, 0.9375, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "AbsQ1_GT_AbsQU__D2_Neg",
			d1:     2.0,
			d2:     -3.0,
			x1:     8.0,
			y1:     5.0,
			d1Want: 4.830188679245283,
			d2Want: -7.245283018867925,
			x1Want: 3.3125,
			pWant: DrotmParams(
				flag: OffDiagonal,
				h:    [0.0, -0.625, -0.9375, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "AbsQ1_LT_AbsQU__D2_Pos",
			d1:     2.0,
			d2:     3.0,
			x1:     5.0,
			y1:     8.0,
			d1Want: 2.3801652892561984,
			d2Want: 1.586776859504132,
			x1Want: 121.0 / 12.0,
			pWant: DrotmParams(
				flag: Diagonal,
				h:    [5.0 / 12.0, 0.0, 0.0, 0.625],
            )
        ),
        drotmgTestStruct(
			name:   "D1=D2_X1=X2",
			d1:     2.0,
			d2:     2.0,
			x1:     8.0,
			y1:     8.0,
			d1Want: 1.0,
			d2Want: 1.0,
			x1Want: 16.0,
			pWant: DrotmParams(
				flag: Diagonal,
				h:    [1.0, 0.0, 0.0, 1.0],
            )
        ),
        drotmgTestStruct(
			name:   "RD1_Big_RD2_Big_Flag_0",
			d1:     1600000000.0,
			d2:     800000000.0,
			x1:     8.0,
			y1:     7.0,
			d1Want: 68.96627824858757,
			d2Want: 34.483139124293785,
			x1Want: 45312.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [4096.0, -3584.0, 1792.0, 4096.0],
            )
        ),
        drotmgTestStruct(
			name:   "RD1_Big_RD2_Big_Flag_1",
			d1:     800000000.0,
			d2:     1600000000.0,
			x1:     8.0,
			y1:     7.0,
			d1Want: 57.6914092640818,
			d2Want: 28.8457046320409,
			x1Want: 47396.57142857142,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [2340.5714285714284, -4096.0, 4096.0, 4681.142857142857],
            )
        ),
        drotmgTestStruct(
			name:   "RD1_Big_RD2_Med_Flag_0",
			d1:     20000000.0,
			d2:     2.0,
			x1:     8.0,
			y1:     8.0,
			d1Want: 1.1920927762985347,
			d2Want: 1.9999998000000199,
			x1Want: 32768.0032768,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [4096.0, -1.0, 0.0004096, 1.0],
            )
        ),
        drotmgTestStruct(
			name:   "RD1_Big_RD2_Med_Flag_1",
			d1:     2.0,
			d2:     20000000000.0,
			x1:     8.0,
			y1:     80000000000.0,
			d1Want: 1192.0928955078125,
			d2Want: 2.0,
			x1Want: 327680000000000.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [4.096e-17, -1.0, 4096.0, 1e-10],
            )
        ),
        drotmgTestStruct(
			name:   "D1_Big_D2_Small_Flag_1",
			d1:     0.000000014,
			d2:     2000000000.0,
			x1:     0.000008,
			y1:     8000000.0,
			d1Want: 119.20928955078125,
			d2Want: 0.234881024,
			x1Want: 32768000000.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [2.8671999999999997e-26, -0.000244140625, 4096.0, 2.44140625e-16],
            )
        ),
        drotmgTestStruct(
			name:   "RD1_Med_RD2_Big_Flag_0",
			d1:     2.0,
			d2:     20000000000.0,
			x1:     80000000.0,
			y1:     8.0,
			d1Want: 1.9998000199980002,
			d2Want: 1191.9736981379988,
			x1Want: 80008000.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [1.0, -0.0004096, 1000.0, 4096.0],
            )
        ),
        drotmgTestStruct(
			name:   "D1_Med_D2_Big_Flag_1",
			d1:     20000000000.0,
			d2:     0.4,
			x1:     80000000.0,
			y1:     80000000000000000.0,
			d1Want: 0.39999998000000103,
			d2Want: 1192.092835903171,
			x1Want: 80000004000000000.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [50.0, -4096.0, 1.0, 4.096e-06],
            )  
        ),
        drotmgTestStruct(
			name:   "RD1_Med_RD2_Small_Flag_0",
			d1:     1.2,
			d2:     0.000000000045,
			x1:     2.7,
			y1:     8.0,
			d1Want: 1.1999999996049382,
			d2Want: 0.0007549747197514486,
			x1Want: 2.700000000888889,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [1.0, -0.0007233796296296296, 1.1111111111111111e-10, 0.000244140625],
            )  
        ),
        drotmgTestStruct(
			name:   "RD1_Med_RD2_Small_Flag_1",
			d1:     1.2,
			d2:     0.000000000045,
			x1:     2.7,
			y1:     80000000000.0,
			d1Want: 0.0007549747199770676,
			d2Want: 1.19999999996355,
			x1Want: 19531250.000593264,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [0.0002197265625, -1.0, 0.000244140625, 3.375e-11],
            )  
        ),
        drotmgTestStruct(
			name:   "D1_Small_D2_Big_Flag_1",
			d1:     120000000000000000.0,
			d2:     0.000000000012345,
			x1:     0.08,
			y1:     8000000000000.0,
			d1Want: 0.00010502490698765249,
			d2Want: 216.1836123957717,
			x1Want: 3851666919.8055897,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [23731773997.569866, -16777216.0, 0.000244140625, 1.6777216e-07],
            )
        ),
        drotmgTestStruct(
			name:   "RD1_Small_RD2_Med_Flag_0",
			d1:     0.0000000002,
			d2:     20.0,
			x1:     0.8,
			y1:     0.000000008,
			d1Want: 0.003355409645903541,
			d2Want: 19.99980000199998,
			x1Want: 0.000195314453125,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [0.000244140625, -1e-08, 0.24414062499999997, 1.0],
            )
        ),
        drotmgTestStruct(
			name:   "RD1_Small_RD2_Med_Flag_1",
			d1:     0.02,
			d2:     0.000000000004,
			x1:     0.008,
			y1:     8000000.0,
			d1Want: 6.710886366445568e-05,
			d2Want: 0.019999999900000003,
			x1Want: 1953.125009765625,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [0.0012207031250000002, -1.0, 0.000244140625, 1e-09],
            )
        ),
        drotmgTestStruct(
			name:   "OpenBLAS#1452",
			d1:     5.9e-8,
			d2:     5.960464e-8,
			x1:     1.0,
			y1:     150.0,
			d1Want: 0.9999559282289687,
			d2Want: 0.9898121986058326,
			x1Want: 0.03662270484346241,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [1.6110934624105326e-06, -0.000244140625, 0.000244140625, 1.6276041666666668e-06],
			)
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#1",
			d1:     0.1,
			d2:     0.3,
			x1:     1.2,
			y1:     0.2,
			d1Want: 12.0 / 130.0,
			d2Want: 36.0 / 130.0,
			x1Want: 1.3,
			pWant: DrotmParams(
				flag: OffDiagonal,
				h:    [0.0, -1.0 / 6.0, 0.5, 0.0],
			)
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#2",
			d1:     0.7,
			d2:     0.2,
			x1:     0.6,
			y1:     4.2,
			d1Want: 14.0 / 75.0,
			d2Want: 49.0 / 75.0,
			x1Want: 4.5,
			pWant: DrotmParams(
				flag: Diagonal,
				h:    [0.5, 0.0, 0.0, 1.0 / 7.0],
            )
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#3",
			d1:     0.0,
			d2:     0.0,
			x1:     0.0,
			y1:     0.0,
			d1Want: 0.0,
			d2Want: 0.0,
			x1Want: 0.0,
			pWant: DrotmParams(
				flag: Identity,
				h:    [0.0, 0.0, 0.0, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#4",
			d1:     4.0,
			d2:     -1.0,
			x1:     2.0,
			y1:     4.0,
			d1Want: 0.0,
			d2Want: 0.0,
			x1Want: 0.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [0.0, 0.0, 0.0, 0.0],
            )
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#5",
			d1:     6e-10,
			d2:     2e-2,
			x1:     1e5,
			y1:     10.0,
			d1Want: 45e-11 * 4096.0 * 4096.0,
			d2Want: 15e-3,
			x1Want: 4e5 / (3.0 * 4096.0),
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [1.0 / 4096.0, -1e-4, 1e4 / (3.0 * 4096.0), 1.0],
            )
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#6",
			d1:     4e10,
			d2:     2e-2,
			x1:     1e-5,
			y1:     10.0,
			d1Want: 4e10 / (1.5 * 4096.0 * 4096.0),
			d2Want: 2e-2 / 1.5,
			x1Want: 6.144e-2,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [4096.0, -1e6, 5e-7 * 4096.0, 1.0],
            )
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#7",
			d1:     2e-10,
			d2:     4e-2,
			x1:     1e5,
			y1:     10.0,
			d1Want: 4.0 / 150.0,
			d2Want: (2e-10 / 1.5) * 4096.0 * 4096.0,
			x1Want: 15.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [5e-5, -1.0 / 4096.0, 1.0, 1e4 / 4096.0],
            )
        ),
        drotmgTestStruct(
			name:   "netlib/BLAS/TESTING#8",
			d1:     2e10,
			d2:     4e-2,
			x1:     1e-5,
			y1:     10.0,
			d1Want: 4.0 / 150.0,
			d2Want: 2e10 / (1.5 * 4096.0 * 4096.0),
			x1Want: 15.0,
			pWant: DrotmParams(
				flag: Rescaling,
				h:    [5e5, -4096.0, 1.0, 4.096e-3],
            )
        )
    ]

    for (test in tests) {
        let (p, d1, d2, x1) = nFloat64.drotmg(test.d1, test.d2, test.x1, test.y1)

        if (p.flag != test.pWant.flag) {
            throw Exception("${test.name}: unexpected flag, want ${test.pWant.flag}, got ${p.flag}")
        }
        for (i in 0..p.h.size) {
            if (!equalWithinAbsOrRel(test.pWant.h[i], p.h[i], tol, tol)) {
                throw Exception("${test.name}: unexpected H, want ${test.pWant.h}, got ${p.h}")
            }
        }
        if (!equalWithinAbsOrRel(d1, test.d1Want, tol, tol)) {
            throw Exception("${test.name}: unexpected d1, want ${test.d1Want}, got ${d1}")
        }
        if (!equalWithinAbsOrRel(d2, test.d2Want, tol, tol)) {
            throw Exception("${test.name}: unexpected d2, want ${test.d2Want}, got ${d2}")
        }
        if (!equalWithinAbsOrRel(x1, test.x1Want, tol, tol)) {
            throw Exception("${test.name}: unexpected x1, want ${test.x1Want}, got ${x1}")
        }

        let y = [test.y1]
        nFloat64.drotm(1, [test.x1], 1, y, 1, p)
        if (d2 >= 0.0) {
            y[0] *= sqrt(d2)
        }
        if (!equalWithinAbs(y[0], 0.0, 1e-10)) {
            throw Exception("${test.name}: y1 not zeroed out, got ${y[0]}")
        }
    } 
    return true
}

public func drotTest(nFloat64: NFloat64): Bool {
    for (c in doubleTwoVectorCases) {
        for (kind in c.drotCases) {
            dCopyTwoTmp(c.x, c.xTmp, c.y, c.yTmp)
            if (c.exception) {
                let f = { => 
                    nFloat64.drot(c.n, c.xTmp, c.incx, c.yTmp, c.incy, kind.c, kind.s)
                    return 
                }
                if (!exception(f)) {
                    throw Exception("${c.name} should panic and does not")
                }
                continue
            }
            nFloat64.drot(c.n, c.xTmp, c.incx, c.yTmp, c.incy, kind.c, kind.s)
            if (!dSliceTolEqual(c.xTmp, kind.xAns)) {
                throw Exception("drot: x mismatch ${c.name}: expected ${kind.xAns}, found ${c.xTmp}")
            }
            if (!dSliceTolEqual(c.yTmp, kind.yAns)) {
                throw Exception("drot: y mismatch ${c.name}: expected ${kind.yAns}, found ${c.yTmp}")
            }
        }
    }
    return true
}

public func drotmTest(nFloat64: NFloat64): Bool {
    for (c in doubleTwoVectorCases) {
        for (kind in c.drotmCases) {
            dCopyTwoTmp(c.x, c.xTmp, c.y, c.yTmp)
            if (c.exception) {
                let f = { => 
                    nFloat64.drotm(c.n, c.xTmp, c.incx, c.yTmp, c.incy, kind.p)
                    return 
                }
                if (!exception(f)) {
                    throw Exception("${c.name} should panic and does not")
                }
                continue
            }
            nFloat64.drotm(c.n, c.xTmp, c.incx, c.yTmp, c.incy, kind.p)
            if (!dSliceTolEqual(c.xTmp, kind.xAns)) {
                throw Exception("drotm: mismatch ${c.name} ${kind.name}: expected ${kind.xAns}, found ${c.xTmp}")
            }
            if (!dSliceTolEqual(c.yTmp, kind.yAns)) {
                throw Exception("drotm: mismatch ${c.name} ${kind.name}: expected ${kind.yAns}, found ${c.yTmp}")
            }
        }
    }
    return true
}

public func dscalTest(nFloat64: NFloat64): Bool {
    for (c in doubleOneVectorCases) {
        for (kind in c.dscalCases) {
            var xTmp = c.x.clone()
            if (c.exception) {
                let f = { => 
                    nFloat64.dscal(c.n, kind.alpha, c.x.clone(), c.incx)
                    return 
                }
                if (!exception(f)) {
                    throw Exception("${c.name} should panic and does not")
                }
                continue
            }
            nFloat64.dscal(c.n, kind.alpha, xTmp, c.incx)
            if (!dSliceTolEqual(xTmp, kind.ans)) {
                throw Exception("dscal: mismatch ${c.name}, ${kind.name}: expected ${kind.ans}, found ${xTmp}")
            }
        }
    }
    return true
}