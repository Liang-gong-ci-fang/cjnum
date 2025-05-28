package cjnum.lapack.cjnum

import cjnum.lapack.*

public class LFloat64Implementation {}

let lFloat64Implementation = LFloat64Implementation()

const dlamchE = 0x1p-53

const dlamchB = 2

const dlamchP = 2.0 * dlamchE

const dlamchS = 0x1p-1022

const dtsml = 0x1p-511
const dtbig = 0x1p486
const dssml = 0x1p537
const dsbig = 0x1p-538

func copysign(f: Float64, sign: Float64): Float64 {
    const signBit: UInt64 = 1 << 63
    return Float64.fromBits((f.toBits() & (!signBit)) | sign.toBits() & signBit)
}

let _ = {=>lFloat64Implementation}()