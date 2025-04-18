/** @file thbSplineBasis_example.cpp

    @brief Tutorial on gsTHBSplineBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

// Look also in basis_example and bSplineBasis_example.

//! [Include namespace]
#include <string>
#include <gismo.h>
#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsElasticity/gsPSOR.h>
#include <gsModeling/gsRBFCurve.h>
//! [Include namespace]

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numHRef = 0;
    index_t numUHRef = 0;
    index_t numElev = 0;
    std::string output;
    std::string input, geometry;


    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Degree elevation",numElev);
    cmd.addInt("r", "numHRef","Number of elements in the crack size", numHRef);
    cmd.addInt("R", "numUHRef","Number of pre-refinements", numUHRef);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("i", "input", "Input XML file", input);
    cmd.addString("g", "geometry", "Geometry XML file", geometry);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(!input.empty(),"Input file not provided");
    GISMO_ASSERT(!geometry.empty(),"Geometry file not provided");


    if (output.empty())
        output = "./output/";

    std::string outputdir = output + gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);

    //! [Parse command line]

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// Read the input files
    //////////////////////////////////////////////////////////////////////////////////////////////////

    gsFileData<> fd(input);
    gsInfo << "Input file "<< fd.lastPath() <<"\n";

    GISMO_ASSERT(fd.hasLabel("crack"), "Crack not found in the input file.");
    GISMO_ASSERT(fd.hasLabel("material"), "Material parameters not found in the input file.");

    gsMultiPatch<> crack;
    fd.getLabel("crack", crack);

    gsOptionList materialParameters;
    fd.getLabel("material", materialParameters);

    real_t l0 = materialParameters.getReal("l0");
    real_t Gc = materialParameters.getReal("Gc");
    real_t beta = materialParameters.getReal("beta");
    index_t order = materialParameters.getInt("order");
    index_t AT = materialParameters.getInt("AT");
    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");

    fd.clear();
    fd.read(geometry);
    gsInfo << "Geometry file "<< fd.lastPath() <<"\n";
    // GISMO_ASSERT(fd.hasLabel("geometry"), "Material parameters not found in the input file.");
    // fd.getLabel("geometry", mp);
    gsMultiPatch<> mp;
    fd.getFirst(mp);
    gsMultiBasis<> mb(mp);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// Initial basis
    //////////////////////////////////////////////////////////////////////////////////////////////////

    //// Compute curve length
    gsMultiBasis<> curveBasis(crack.basis(0));
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(curveBasis);
    auto C = ev.getMap(crack);
    real_t length = ev.integral(1.0 * meas(C));
    gsInfo<<"Curve length = "<<length<<"\n";

    gsRBFCurve<real_t,Constant> fun(crack,beta,beta);
    if(plot) gsWriteParaview(mp,fun,outputdir+"initial",100000);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// Local projection
    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// L2 projection
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(mb);
    auto w = A.getSpace(mb);
    auto c = A.getCoeff(fun);
    auto G = A.getMap(mp);
    w.setup();
    A.initSystem();
    A.assemble( w * w.tr() * meas(G), w * c * meas(G) );
    gsMatrix<> mm = A.matrix() * gsMatrix<>::Ones(A.matrix().rows(),1);
    // count nonzeros in the rhs
    index_t nnz = (A.rhs().array() > 0).count();
    gsMatrix<> f(nnz,1), m(nnz,1);
    for (index_t i = 0, j = 0; i < A.rhs().rows(); ++i)
        if (A.rhs()(i,0) > 0)
        {
            f(j,0) = A.rhs()(i,0);
            m(j,0) = mm(i,0);
            ++j;
        }

    gsMatrix<> tmpCoefs = m.cwiseInverse().cwiseProduct(f);
    gsMatrix<> coefs(w.mapper().freeSize(),1);
    coefs.setZero();
    for (index_t i = 0, j = 0; i < w.mapper().size(); ++i)
        if (A.rhs()(i,0) > 0)
            coefs(i,0) = tmpCoefs(j++,0);

    gsMultiPatch<> damage;
    damage.addPatch(mb.basis(0).makeGeometry(give(coefs)));
    if(plot) gsWriteParaview(mp,damage,outputdir+"L2",10000);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// PF assembler
    //////////////////////////////////////////////////////////////////////////////////////////////////
    gsBoundaryConditions<> bc_d;
    bc_d.setGeoMap(mp);
    gsPhaseFieldAssemblerBase<real_t> * pfAssembler;
    if      (order == 2 && AT == 1)
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Second,PFmode::AT1>(mp,mb,bc_d);
    else if (order == 4 && AT == 1)
    {
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Fourth,PFmode::AT1>(mp,mb,bc_d);
        pfAssembler->options().setReal("cw",4.44847);
    }
    else if (order == 2 && AT == 2)
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Second,PFmode::AT2>(mp,mb,bc_d);
    else if (order == 4 && AT == 2)
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Fourth,PFmode::AT2>(mp,mb,bc_d);
    else
        GISMO_ERROR("Invalid order and/or AT model");

    pfAssembler->options().setReal("l0",l0);
    pfAssembler->options().setReal("Gc",Gc);
    pfAssembler->initialize();

    gsInfo<<"Assembling phase-field problem"<<std::flush;
    gsMatrix<> D, deltaD;
    pfAssembler->constructSolution(damage,D);
    gsSparseMatrix<> Q, QPhi;
    gsMatrix<> q;
    gsMatrix<> R;
    pfAssembler->assembleMatrix();
    pfAssembler->matrix_into(QPhi);
    Q = QPhi;
    pfAssembler->assembleVector();
    pfAssembler->rhs_into(q);
    gsInfo<<". Done\n";

    R = Q * D + q;

    real_t energy;
    for (index_t it = 0; it!=10000; it++)
    {
        gsInfo<<"Iteration "<<it<<":"<<std::flush;
        gsPSOR<real_t> PSORsolver(Q);
        PSORsolver.options().setInt("MaxIterations",3000);
        PSORsolver.options().setSwitch("Verbose",false);
        PSORsolver.options().setReal("tolU",1e-4);
        PSORsolver.options().setReal("tolNeg",1e-6);
        PSORsolver.options().setReal("tolPos",1e-6);
        PSORsolver.solve(R,deltaD); // deltaD = Q \ R
        D += deltaD;

        R = Q * D + q;
        gsInfo<<"||dD|| = "<<deltaD.norm()<<", ||E|| = "<<R.norm()<<"\n";

        if (deltaD.norm()/D.norm() < 1e-10)
        {
            gsInfo<<"Converged\n";
            break;
        }
    }

    energy = (0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value();
    gsInfo<<"E_d   = "<<energy<<"\n";

    pfAssembler->constructSolution(D,damage);

    if(plot) gsWriteParaview(mp,damage,outputdir+"L2_after",100000);

    gsFileData<> fd_out;
    fd_out.addWithLabel(damage,outputdir+"initial");
    fd_out.save(outputdir+"initial");

    delete pfAssembler;
    return EXIT_SUCCESS;
}

