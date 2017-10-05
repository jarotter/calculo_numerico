#Pkg.add("MatrixDepot")
#Pkg.update("MatrixDepot")
include("lib/qr.jl")
using(MatrixDepot)

A = matrixdepot("fiedler", 25);
qr_simple(A, 2000, 1e-10)
qr_dynamic(A, 25, 1e-10)
eigfact(A)[:values];
