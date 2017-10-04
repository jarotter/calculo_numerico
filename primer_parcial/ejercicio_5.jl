#Pkg.add("MatrixDepot")
#Pkg.update("MatrixDepot")
include("lib/qr.jl")
using(MatrixDepot)

A = matrixdepot("fiedler", 25);
qr_simple(A, 1000, 1e-10)
qr_dynamic(A, 20, 1e-10)
eigfact(A)[:values];
