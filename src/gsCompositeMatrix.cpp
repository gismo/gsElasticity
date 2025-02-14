/** @file gsCompositeMatrix.h

    @brief Provides a simple class to defined a composite matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

// #pragma once

#include <gismo.h>

using namespace gismo;

template<class T>
gsMatrix<T> gsCompositeMatrix(	const T Exx,
	                            const T Eyy,
	                            const T Ezz,
	                            const T Gxy,
	                            const T Gyz,
	                            const T Gxz,
	                            const T nuxy,
	                            const T nuyz,
	                            const T nuxz
	                          )
{
	gsMatrix<T> G(6,6);
	G.setZero();
	T D = (1-nuxy*nuxy-nuyz*nuyz-nuxz*nuxz-2*nuxy*nuyz*nuxz) / (Exx*Eyy*Ezz);
	G(0,0) = (1-nuyz*nuyz) / (Eyy*Ezz); // Gzz
	G(1,1) = (1-nuxz*nuxz) / (Exx*Ezz); // Gyy
	G(2,2) = (1-nuxy*nuxy) / (Exx*Eyy); // Gzz

	// Gij = nuij + nuik*nukj / Ejj*Ekk
	G(0,1) = G(1,0) = (nuxy+nuxz*nuyz) / (Eyy*Ezz); // Gxy
	G(0,2) = G(2,0) = (nuxz+nuxy*nuyz) / (Eyy*Ezz); // Gxz
	G(1,2) = G(2,1) = (nuyz+nuxy*nuxz) / (Exx*Ezz); // Gyz

	G *= 1.0/D;

	G(3,3) = Gxy; // Factor 2?
	G(4,4) = Gyz; // Factor 2?
	G(5,5) = Gxz; // Factor 2?
	return G;
}

template<class T>
gsMatrix<T> gsCompositeMatrixRaw(	const T G11,const T G12,const T G13,const T G14,const T G15,const T G16,
												const T G22,const T G23,const T G24,const T G25,const T G26,
	                            							const T G33,const T G34,const T G35,const T G36,
	                            										const T G44,const T G45,const T G46,
	                            													const T G55,const T G56,
	                            																const T G66
	                          	)
{
	gsMatrix<T> G(6,6);
	G.setZero();
	G(0,1) = G12;
	G(0,2) = G13;
	G(0,3) = G14;
	G(0,4) = G15;
	G(0,5) = G16;

	G(1,2) = G23;
	G(1,3) = G24;
	G(1,4) = G25;
	G(1,5) = G26;

	G(2,3) = G34;
	G(2,4) = G35;
	G(2,5) = G36;

	G(3,4) = G45;
	G(3,5) = G46;

	G(4,5) = G56;

	G = G + G.transpose();
	G(0,0) = G11;
	G(1,1) = G22;
	G(2,2) = G33;
	G(3,3) = G44;
	G(4,4) = G55;
	G(5,5) = G66;

	return G;
}

template<class T>
gsMatrix<T> gsCompositeMatrix(	const T Exx,
	                            const T Eyy,
	                            const T Gxy,
	                            const T nuxy
	                          )
{
	gsMatrix<T> G(3,3);
	G.setZero();

	G(0,0) = Exx*nuxy / ((1+nuxy)*(1-2*nuxy)) + Exx / ((1+nuxy)); // Gxx
	G(1,1) = Eyy*nuxy / ((1+nuxy)*(1-2*nuxy)) + Eyy / ((1+nuxy)); // Gxx
	G(0,1) = G(1,0) = (nuxy*Exx) / ((1+nuxy)*(1-2*nuxy)); // Gxy

	G(2,2) = Gxy;
	return G;
}

template<class T>
gsMatrix<T> gsCompositeMatrixRaw(	const T G11,const T G12,const T G13,
												const T G22,const T G23,
	                            							const T G33
	                          	)
{
	gsMatrix<T> G(3,3);
	G.setZero();
	G(0,1) = G12;
	G(0,2) = G13;

	G(1,2) = G23;

	G = G + G.transpose();
	G(0,0) = G11;
	G(1,1) = G22;
	G(2,2) = G33;

	return G;
}