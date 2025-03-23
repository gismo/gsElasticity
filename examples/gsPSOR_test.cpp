/** @file sparseSolvers_example.cpp

    @brief Testing the use of sparse linear solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>
#include <gsElasticity/gsPSOR.h>

using namespace gismo;


int main(int argc, char** argv)
{
    std::string fn("");
    index_t mat_size = 10;

    gsCmdLine cmd("Testing the use of sparse linear solvers.");
    cmd.addPlainString("try", "Name of the solver to try", fn);
    cmd.addInt("n", "size", "Size of the matrices", mat_size);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Construct the geometry
    index_t N = 8;
    index_t p = 2;
    gsMultiPatch<> mp;
    gsKnotVector<> kv(0,1,N-1,p+1,1,p); // start, end, num. interior knots (nelements-1), multiplicity of end knots, degree
    gsTensorBSplineBasis<2,real_t> basis(kv,kv);
    gsMatrix<> coefs = basis.anchors().transpose();
    mp.addPatch(basis.makeGeometry(coefs));
    // mp.degreeIncrease(numElev);

    // Construct the basis
    gsMultiBasis<> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";

    gsMatrix<> cSupp(2,2);
    cSupp.row(0)<<0.0,0.5; // x
    cSupp.row(1)<<0.4,0.6; // y
    gsFunctionExpr<> c_ini("if((x>="+util::to_string(cSupp(0,0))+") and \
                               (x<="+util::to_string(cSupp(0,1))+") and \
                               (y>="+util::to_string(cSupp(1,0))+") and \
                               (y<="+util::to_string(cSupp(1,1))+"), 1, .9)",2);

    /// L2 projection
    gsMatrix<> x;
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(mb);
    auto w = A.getSpace(mb);
    auto c = A.getCoeff(c_ini);
    auto G = A.getMap(mp);
    auto sol = A.getSolution(w,x);
    w.setup();
    A.initSystem();
    A.assemble( w * w.tr() * meas(G), w * c * meas(G) );

    gsSparseMatrix<> M = A.matrix();
    gsMatrix<> F = A.rhs();

    gsInfo<<"M.minCoeff() = "<<M.coeffs().minCoeff()<<"\n";
    gsInfo<<"M.maxCoeff() = "<<M.coeffs().maxCoeff()<<"\n";
    gsInfo<<"F.minCoeff() = "<<F.minCoeff()<<"\n";
    gsInfo<<"F.maxCoeff() = "<<F.maxCoeff()<<"\n";

    gsPSOR<real_t> PSOR(M);
    PSOR.options().setInt("MaxIterations",1000);
    PSOR.options().setReal("tolU",1e-10);
    PSOR.options().setReal("tolPos",1e-10);
    PSOR.options().setReal("tolNeg",1e-10);
    PSOR.solve(-F,x);
    gsInfo<<"PSOR iterations = "<<PSOR.numIter()<<"\n";

    gsDebug<<"x^T (M x + q) = "<<x.transpose()*(M*x - F)<<"\n";
    gsDebug<<"M x + b       = "<<(M*x - F).transpose().minCoeff()<<"\n";
    gsDebug<<"x             =  "<<x.transpose().minCoeff()<<"\n";

    gsExprEvaluator<> ev(A);
    gsInfo<<"x.minCoeff() = "<<x.minCoeff()<<"\n";
    gsInfo<<"x.maxCoeff() = "<<x.maxCoeff()<<"\n";
    gsInfo<<"error = "<<ev.integral((sol-c)*(sol-c)*meas(G))<<"\n";

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
#endif
    solver.compute(M);
    x = solver.solve(F);
    gsInfo<<"x.minCoeff() = "<<x.minCoeff()<<"\n";
    gsInfo<<"x.maxCoeff() = "<<x.maxCoeff()<<"\n";
    gsInfo<<"error = "<<ev.integral((sol-c)*(sol-c)*meas(G))<<"\n";



    return EXIT_SUCCESS;
}
