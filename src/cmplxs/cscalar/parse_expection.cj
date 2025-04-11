package cjnum.cmplxs.cscalar


public class ParseExpection <: Exception {
    public ParseExpection(
        information!: String = ""
    ) {
        super(information)
    }
}