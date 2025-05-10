package cjnum.lapack.cjnum

import cjnum.lapack.*

public class LFloat64Implementation {}

let lFloat64Implementation = LFloat64Implementation()

let dlamchE = 0x1p-53

let dlamchB = 2

let dlamchP = Float64(dlamchB) * dlamchE

let dlamchS = 0x1p-1022

let dtsml = 0x1p-511
let dtbig = 0x1p486
let dssml = 0x1p537
let dsbig = 0x1p-538

func copysign(f: Float64, sign: Float64): Float64 {
    const signBit: UInt64 = 1 << 63
    return Float64.fromBits((f.toBits() & (!signBit)) | sign.toBits() & signBit)
}
