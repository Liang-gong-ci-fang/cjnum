package cjnum.mat

import cjnum.blas.blas64.*
import std.random.*
import std.unittest.*
import std.unittest.testmacro.*

// TODO 
// // legalDims returns whether {m,n} is a valid dimension of the given matrix type.
// private func legalDims(a: Matrix, m:Int64, n: Int64): Bool {
// 	if(a is Untransposer){
//         let a_=a as Untransposer
//         return legalDims(a_.getOrThrow().mUntranspose(), n, m)
//     }
//     else if(a is Dense||a is BasicMatrix||a is BandDense||a is BasicMatrix){
//         if (m < 0 || n < 0) {
// 			return false
// 		}
// 		return true
//     }
//     else{
//         throw Exception("legal dims type not coded")
//     }
//     // TODO
// 	// case *SymDense, *TriDense, *basicSymmetric, *basicTriangular,
// 	// 	*SymBandDense, *basicSymBanded, *TriBandDense, *basicTriBanded,
// 	// 	*basicDiagonal, *DiagDense, *Tridiag:
// 	// 	if m < 0 || n < 0 || m != n {
// 	// 		return false
// 	// 	}
// 	// 	return true
// 	// case *VecDense, *basicVector:
// 	// 	if m < 0 || n < 0 {
// 	// 		return false
// 	// 	}
// 	// 	return n == 1
	
// }

// // testMatrices is a list of matrix types to test.
// // This test relies on the fact that the implementations of Triangle do not
// // corrupt the value of Uplo when they are empty. This test will fail
// // if that changes (and some mechanism will need to be used to force the
// // correct TriKind to be read).
// private let testMatrices = Array<Matrix>(
//     [
//         Dense(General(),0,0),
// 	    BasicMatrix(General(),0,0),
// 	    Transpose(Dense(General(),0,0))

//     //,
// 	// &VecDense{mat: blas64.Vector{Inc: 1}},
// 	// &VecDense{mat: blas64.Vector{Inc: 10}},
// 	// &basicVector{},
// 	// Transpose{&VecDense{mat: blas64.Vector{Inc: 1}}},
// 	// Transpose{&VecDense{mat: blas64.Vector{Inc: 10}}},
// 	// Transpose{&basicVector{}},

// 	// &BandDense{mat: blas64.Band{KL: 2, KU: 1}},
// 	// &BandDense{mat: blas64.Band{KL: 1, KU: 2}},
// 	// Transpose{&BandDense{mat: blas64.Band{KL: 2, KU: 1}}},
// 	// Transpose{&BandDense{mat: blas64.Band{KL: 1, KU: 2}}},
// 	// TransposeBand{&BandDense{mat: blas64.Band{KL: 2, KU: 1}}},
// 	// TransposeBand{&BandDense{mat: blas64.Band{KL: 1, KU: 2}}},

// 	// &SymDense{},
// 	// &basicSymmetric{},
// 	// Transpose{&basicSymmetric{}},

// 	// &TriDense{mat: blas64.Triangular{Uplo: blas.Upper}},
// 	// &TriDense{mat: blas64.Triangular{Uplo: blas.Lower}},
// 	// &basicTriangular{mat: blas64.Triangular{Uplo: blas.Upper}},
// 	// &basicTriangular{mat: blas64.Triangular{Uplo: blas.Lower}},
// 	// Transpose{&TriDense{mat: blas64.Triangular{Uplo: blas.Upper}}},
// 	// Transpose{&TriDense{mat: blas64.Triangular{Uplo: blas.Lower}}},
// 	// TransposeTri{&TriDense{mat: blas64.Triangular{Uplo: blas.Upper}}},
// 	// TransposeTri{&TriDense{mat: blas64.Triangular{Uplo: blas.Lower}}},
// 	// Transpose{&basicTriangular{mat: blas64.Triangular{Uplo: blas.Upper}}},
// 	// Transpose{&basicTriangular{mat: blas64.Triangular{Uplo: blas.Lower}}},
// 	// TransposeTri{&basicTriangular{mat: blas64.Triangular{Uplo: blas.Upper}}},
// 	// TransposeTri{&basicTriangular{mat: blas64.Triangular{Uplo: blas.Lower}}},

// 	// &SymBandDense{},
// 	// &basicSymBanded{},
// 	// Transpose{&basicSymBanded{}},

// 	// &SymBandDense{mat: blas64.SymmetricBand{K: 2}},
// 	// &basicSymBanded{mat: blas64.SymmetricBand{K: 2}},
// 	// Transpose{&basicSymBanded{mat: blas64.SymmetricBand{K: 2}}},
// 	// TransposeBand{&basicSymBanded{mat: blas64.SymmetricBand{K: 2}}},

// 	// &TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}},
// 	// &TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}},
// 	// &basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}},
// 	// &basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}},
// 	// Transpose{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// Transpose{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},
// 	// Transpose{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// Transpose{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},
// 	// TransposeTri{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// TransposeTri{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},
// 	// TransposeTri{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// TransposeTri{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},
// 	// TransposeBand{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// TransposeBand{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},
// 	// TransposeBand{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// TransposeBand{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},
// 	// TransposeTriBand{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// TransposeTriBand{&TriBandDense{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},
// 	// TransposeTriBand{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Upper}}},
// 	// TransposeTriBand{&basicTriBanded{mat: blas64.TriangularBand{K: 2, Uplo: blas.Lower}}},

// 	// &DiagDense{},
// 	// &DiagDense{mat: blas64.Vector{Inc: 10}},
// 	// Transpose{&DiagDense{}},
// 	// Transpose{&DiagDense{mat: blas64.Vector{Inc: 10}}},
// 	// TransposeTri{&DiagDense{}},
// 	// TransposeTri{&DiagDense{mat: blas64.Vector{Inc: 10}}},
// 	// TransposeBand{&DiagDense{}},
// 	// TransposeBand{&DiagDense{mat: blas64.Vector{Inc: 10}}},
// 	// TransposeTriBand{&DiagDense{}},
// 	// TransposeTriBand{&DiagDense{mat: blas64.Vector{Inc: 10}}},
// 	// &basicDiagonal{},
// 	// Transpose{&basicDiagonal{}},
// 	// TransposeTri{&basicDiagonal{}},
// 	// TransposeBand{&basicDiagonal{}},
// 	// TransposeTriBand{&basicDiagonal{}},

// 	// &Tridiag{},
// 	// Transpose{&Tridiag{}},
// 	// TransposeBand{&Tridiag{}},
//     ]
// )

// private let sizes = Array<(Int64,Int64)>([
// 	(1, 1),
// 	(1, 3),
// 	(3, 1),

// 	(6, 6),
// 	(6, 11),
// 	(11, 6)
// ])

// // @Test
// public func testOneInputFunc(
// 	// name is the name of the function being tested.
// 	name: String,

// 	// f is the function being tested.
// 	f: (a: Matrix)->Any ,

// 	// denseComparison performs the same operation, but using Dense matrices for
// 	// comparison.
// 	denseComparison: (a: Dense)->Any ,

// 	// sameAnswer compares the result from two different evaluations of the function
// 	// and returns true if they are the same. The specific function being tested
// 	// determines the definition of "same". It may mean identical or it may mean
// 	// approximately equal.
// 	sameAnswer: (a:Any, b: Any)-> Bool,

// 	// legalType returns true if the type of the input is a legal type for the
// 	// input of the function.
// 	legalType: (a: Matrix)-> Bool,

// 	// legalSize returns true if the size is valid for the function.
// 	legalSize: (r:Int64, c: Int64)-> Bool
// ) {
// 	let src = Random(1)
// 	for (aMat in testMatrices) {
// 		for (test in sizes) {
// 			// Skip the test if the argument would not be assignable to the
// 			// method's corresponding input parameter or it is not possible
// 			// to construct an argument of the requested size.
// 			if (!legalType(aMat)) {
// 				continue
// 			}
// 			if (!legalDims(aMat, test.ar, test.ac)) {
// 				continue
// 			}
// 			a := makeRandOf(aMat, test.ar, test.ac, src)

// 			// Compute the true answer if the sizes are legal.
// 			dimsOK := legalSize(test.ar, test.ac)
// 			var want interface{}
// 			if dimsOK {
// 				var aDense Dense
// 				aDense.CloneFrom(a)
// 				want = denseComparison(&aDense)
// 			}
// 			aCopy := makeCopyOf(a)
// 			// Test the method for a zero-value of the receiver.
// 			aType, aTrans := untranspose(a)
// 			errStr := fmt.Sprintf("%v(%T), size: %#v, atrans %t", name, aType, test, aTrans)
// 			var got interface{}
// 			panicked, err := panics(func() { got = f(a) })
// 			if !dimsOK && !panicked {
// 				t.Errorf("Did not panic with illegal size: %s", errStr)
// 				continue
// 			}
// 			if dimsOK && panicked {
// 				t.Errorf("Panicked with legal size: %s: %v", errStr, err)
// 				continue
// 			}
// 			if !equal(a, aCopy) {
// 				t.Errorf("First input argument changed in call: %s", errStr)
// 			}
// 			if !dimsOK {
// 				continue
// 			}
// 			if !sameAnswer(want, got) {
// 				t.Errorf("Answer mismatch: %s; got %v, want %v", errStr, got, want)
// 			}
// 		}
// 	}
// }
