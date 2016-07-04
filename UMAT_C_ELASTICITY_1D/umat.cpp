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

#ifndef elasticity_1D_umat_cpp
#define elasticity_1D_umat_cpp

#include "options.cpp"
#include "local_data.cpp"

namespace umat_c_elasticity_1d {

template <class TUmatData>
void umat(TUmatData &umat_data)
{
    typedef LocalData<TUmatData> TLocalData;

    /// localdata
    TLocalData local_data(umat_data);
//     PRINT(local_data.epsilon);

    /// stiffness
    local_data.H0(0, 0) = local_data.E0;
//     PRINT(local_data.H0);

    /// stress
    local_data.sigma = local_data.H0 * local_data.epsilon;
//     PRINT(local_data.sigma);

    /// ener
    local_data.ener = dot_vec_col(local_data.sigma, local_data.epsilon)/2;
//     PRINT(local_data.ener);

} // void umat

} // namespace umat_c_elasticity_1d

#endif // #ifndef elasticity_1D_umat_cpp
