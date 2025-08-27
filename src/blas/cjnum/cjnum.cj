package cjnum.blas.cjnum

const blockSize   = 64 
const minParBlock = 4

public class NFloat64Implementation {}

let nFloat64Implementation = NFloat64Implementation()

func blocks(dim: Int64, bsize: Int64): Int64 {
    return (dim + bsize - 1) / bsize
}

// TODO