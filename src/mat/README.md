
### 1.1 已实现的测试

1. src\mat\band_test.cj: 实现了所有的测试,同时添加了BandDense类型对象的比较相等的方法。
2. src\mat\diagonal_test.cj: 实现了所有的测试，因为缺乏对randDiagDense方法的测试，所以这里没有实现该方法，同时添加了DiagDense类型对象的比较相等的方法。

### 1.2 未完成的工作

1. src\mat\list_test.cj：整个文件的函数没有实现。
2. src\mat\matrix_test.cj: func TestCol，func TestRow，func TestCond，func TestDet，func TestDot，func TestEqual，func TestMax...func TestMulVecToer没有实现。
3. 与src\mat\band_test.cj和src\mat\diagonal_test.cj测试函数无关的类，接口和对应方法没有实现。