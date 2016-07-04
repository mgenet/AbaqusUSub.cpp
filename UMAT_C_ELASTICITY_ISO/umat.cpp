////////////////////////////////////////////////////////////////////////
///                                                                  ///
/// Created by Martin Genet, 2008-2016                               ///
///                                                                  ///
/// Laboratoire de Mécanique et de Technologie (LMT), Cachan, France ///
/// Lawrence Berkeley National Laboratory, California, USA           ///
/// University of California at San Francisco, USA                   ///
/// Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ///
/// École Polytechnique, Palaiseau, France                           ///
///                                                                  ///
////////////////////////////////////////////////////////////////////////

#ifndef elasticity_iso_umat_cpp
#define elasticity_iso_umat_cpp

#include "options.cpp"
#include "local_data.cpp"

/**
 *
 * This is the specific UMat function for isotropic elasticity (2D and 3D).
 * First, a LocalData structure is created, that contains all necessary data to compute the stress, energies, and search direction.
 *
 */

namespace umat_c_elasticity_iso {

template <class TUmatData>
void umat (TUmatData &umat_data)
{
    typedef LocalData<TUmatData> TLocalData;

    /// localdata
    TLocalData local_data(umat_data);
//     PRINT(local_data.epsilon);

    /// stiffness
    local_data.H0 = elasticity_isotrope_H<ndim> (
        local_data.E0,
        local_data.N0);
//     PRINT(local_data.H0);

    /// stress
    local_data.sigma = local_data.H0 * local_data.epsilon;
//     PRINT(local_data.sigma);

    /// ener
    local_data.ener = dot_vec_col(local_data.sigma, local_data.epsilon)/2;
//     PRINT(local_data.ener);

} // void umat

} // namespace umat_c_elasticity_iso

#endif // #ifndef elasticity_iso_umat_cpp
