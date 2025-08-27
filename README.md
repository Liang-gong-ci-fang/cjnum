<div align="center">
<img alt="" src="https://raw.gitcode.com/naxida/cjnum/attachment/uploads/0454254b-8d77-4615-bf1c-ba30a7ad44f4/cjnum.jpg" />
</div>

<p></p>

<p align="center">
<img alt="" src="https://img.shields.io/badge/release-v0.1.0-brightgreen" style="display: inline-block;" />
<img alt="" src="https://img.shields.io/badge/cjc-v1.0.1-brightgreen" style="display: inline-block;" />
<!-- <img alt="" src="https://img.shields.io/badge/cjcov-1.0%25-brightgreen" style="display: inline-block;" /> -->
<img alt="" src="https://img.shields.io/badge/state-孵化-brightgreen" style="display: inline-block;" />
<!-- <img alt="" src="https://img.shields.io/badge/domain-HOS/Cloud-brightgreen" style="display: inline-block;" /> -->
</p>

<p></p>

<div align="center">
<img alt="" src="https://raw.gitcode.com/Cangjie-SIG/cjgrapht/attachment/uploads/f9ffe377-e653-4fa3-b03e-e6ade35594a3/MEGAPROJECT.jpg" />
</div>

<p></p>

## <img alt="" src="./doc/readme-image/readme-icon-introduction.png" style="display: inline-block;" width=3%/> 1 介绍

### 1.1 项目特性

1. 这是一个用于 Cangjie 语言的数值计算库，提供了广泛的数学、科学计算和数值分析功能。

2. 本项目基于 go 版本的 gonum(https://github.com/gonum/gonum)进行开发和改进，原项目是一个用于数值计算的开源库。本项目对其进行了移植，以适应仓颉环境的需求。

### 1.2 项目计划

1. 逐步实现更多的数值计算功能，包括矩阵运算，统计计算，数值积分，优化算法等。

### 1.3 已实现的功能

1. package blas: 基本线性代数子程序。已实现blas功能的类型有Float64。
2. package floats: 基本浮点运算。包括求切片，求点积，求最值等。

## <img alt="" src="./doc/readme-image/readme-icon-framework.png" style="display: inline-block;" width=3%/> 2 架构

### 2.1 项目结构

```shell
.
├── README.md
├── LICENSE
├── CHANGELOG
├── cjpm.toml
├── doc
│   └── readme-image
└── src                             
    ├─ blas
    │  ├─ blas32
    │  ├─ blas64
    │  ├─ cblas128
    │  ├─ cjnum
    │  ├─ testblas
    │  └─ blas.cj                                 
    ├─ cmplxs
    │  └─ cacalar   
    ├─ complex
    ├─ floats
    │  └─ scalar
    ├─ internal
    │  ├─ asm
    │  │  ├─ f32
    │  │  └─ f64
    │  └─ math32 
    └─ lapack
        ├─ cjnum
        └─ lapack.cj 
```

### 2.2 接口说明
#### 2.2.1 blas接口

```cangjie
// Float64类型的blas接口
public interface NFloat64 <: Float64Level1 & Float64Level2 & Float64Level3 {}

// 向量-向量运算，处理单向量操作或双向量的逐元素运算
public interface Float64Level1

// 矩阵-向量运算，实现矩阵与向量的线性变换
public interface Float64Level2

// 矩阵-矩阵运算，稠密矩阵乘法及其变体
public interface Float64Level3
```

### 3.1 编译构建（Win/Linux/Mac）

```shell
cjpm build
```

### 3.2 功能示例
#### 3.2.1 求向量的点积

```cangjie
import cjnum.blas.*
import cjnum.blas.blas64.*

let nFloat64 = nFloat64Implementation()

let x = [10.0, 15.0, -6.0, 3.0, 14.0, 7.0]
let y = [8.0, -2.0, 4.0, 7.0, 6.0, -3.0]
let incX = 1
let incX = 1
let n = 6

let dot = nFloat64.ddot(n, x, incX, y, incY)
// result: dot = 110.0
```

#### 3.2.2 三角矩阵乘向量

```cangjie
import cjnum.blas.*
import cjnum.blas.blas64.*

let nFloat64 = nFloat64Implementation()

let n = 3
let a = [
    [5.0, 0.0, 0.0], 
    [6.0, 9.0, 0.0], 
    [7.0, 10.0, 13.0]
]
let x = [3.0, 4.0, 5.0]
let d = NonUnitDiag // 对角元素不全为1
let ul = Lower      // 下三角矩阵
let tA = Trans      // 取矩阵 a 的转置
let incX = 1   // 步长

let aFlat = flatten(a)  // 将矩阵 a 扁平化，得到[5, 0, 0, 6, 9, 0, 7, 10, 13]
nFloat64.dtrmv(ul, tA, d, n, aFlat, n, x, incX)
// result: x = [74.0, 86.0, 65.0]
```

#### 3.2.3 通用矩阵乘法

```cangjie
import cjnum.blas.*
import cjnum.blas.blas64.*

let nFloat64 = NFloat64Implementation()

// m, n, k矩阵维度, a是m * k, b是k * n, c是m * n
let m = 4
let n = 3
let k = 2
let alpha = 2.0  // 标量系数
let beta =  0.5  // 标量系数
let a = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
    [10.0, 11.0]
]
let b = [
    [1.0, 5.0, 6.0],
    [5.0, -8.0, 8.0]
]
let c = [
    [4.0, 8.0, -9.0],
    [12.0, 16.0, -8.0],
    [1.0, 5.0, 15.0],
    [-3.0, -4.0, 7.0]
]

let tA = NoTrans // 非转置
let tB = NoTrans 
let aFlat = flatten(a) // 扁平化矩阵
let bFlat = flatten(b)
let cFlat = flatten(c)
let lda = a[0].size
let ldb = b[0].size
let ldc = c[0].size

nFloat64.dgemm(tA, tB, m, n, k, alpha, aFlat, lda, bFlat, ldb, beta, cFlat, ldc)
// result: cFlat = [24.0, -18.0, 39.5, 64.0, -32.0, 124.0, 94.5, -55.5, 219.5, 128.5, -78.5, 299.5]
```

## <img alt="" src="./doc/readme-image/readme-icon-contribute.png" style="display: inline-block;" width=3%/> 4 参与贡献

本项目由 [SIGCANGJIE / 仓颉兴趣组](https://gitcode.com/SIGCANGJIE) 实现并维护。技术支持和意见反馈请提Issue。

本项目基于  Apache License 2.0，欢迎给我们提交PR，欢迎参与任何形式的贡献。

本项目committer：[@naxida](https://gitcode.com/naxida)

This project is supervised by [@zhangyin_gitcode](https://gitcode.com/zhangyin_gitcode) (HUAWEI Developer Advocate).

![](https://raw.gitcode.com/SIGCANGJIE/homepage/attachment/uploads/9b648c07-efc2-4eb3-b02f-eab18c77beea/devadvocate.png)