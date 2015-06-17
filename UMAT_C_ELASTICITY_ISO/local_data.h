////////////////////////////////////////////////////////////////////////
///                                                                  ///
/// Created by Martin Genet, 2008-2015                               ///
///                                                                  ///
/// Laboratoire de Mécanique et de Technologie (LMT), Cachan, France ///
/// Lawrence Berkeley National Laboratory, California, USA           ///
/// University of California at San Francisco, USA                   ///
/// Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ///
///                                                                  ///
////////////////////////////////////////////////////////////////////////

#ifndef elasticity_iso_local_data_cpp
#define elasticity_iso_local_data_cpp

#include "containers/vec_mat_tools.h"

#include "UTILS/ELASTICITY/elasticity_isotrope_H.cpp"

#include "UTILS/ELASTICITY/elasticity_isotrans_H.cpp"
#include "UTILS/ELASTICITY/elasticity_isotrope_H.cpp"
#include "UTILS/ELASTICITY/elasticity_isotrope_S.cpp"
#include "UTILS/ELASTICITY/elasticity_orthoeq_H.cpp"
#include "UTILS/ELASTICITY/elasticity_orthoeq_S.cpp"
#include "UTILS/ELASTICITY/elasticity_orthotrope_H.cpp"
#include "UTILS/ELASTICITY/elasticity_orthotrope_S.cpp"

/**
 *
 * This structure contains all necessary data to compute the stress, energies, and search direction.
 * Intelligible names are given for better readability.
 * If the variable already exists 'as is' in the UMatData structure, it is defined here as a reference on the other one. If not, it is defined normally.
 *
 */

namespace umat_c_elasticity_iso {

template <class TUmatData>
struct LocalData
{
    typedef typename TUmatData::TVec_nvec    TVec_nvec;
    typedef typename TUmatData::TMatGen_nvec TMatGen_nvec;

    /// fields
    TVec_nvec epsilon, &sigma;

    /// elasticity parameters
    double &E0, &N0;

    /// energy
    double &ener;

    /// research direction
    TMatGen_nvec &H0;

    LocalData (const TUmatData &umat_data) : sigma(umat_data.stress),
                                             E0(umat_data.props[0]),
                                             N0(umat_data.props[1]),
                                             H0(umat_data.ddsdde),
                                             ener(umat_data.sse)
    {
        epsilon = umat_data.stran + umat_data.dstran;
//         PRINT(epsilon);
    }

}; // struct LocalData

} // namespace umat_c_elasticity_iso

#endif // #ifndef elasticity_iso_local_data_cpp
