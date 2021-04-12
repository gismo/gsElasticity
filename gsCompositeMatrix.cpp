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

#pragma once

namespace gismo
{
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
}
