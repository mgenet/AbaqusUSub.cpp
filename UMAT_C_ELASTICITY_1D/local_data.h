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

#ifndef elasticity_1D_local_data_h
#define elasticity_1D_local_data_h

#include "containers/vec_mat_tools.h"

namespace umat_c_elasticity_1d {

template <class TUmatData>
struct LocalData
{
    typedef typename TUmatData::TVec_nvec    TVec_nvec;
    typedef typename TUmatData::TMatGen_nvec TMatGen_nvec;

    /// fields
    TVec_nvec epsilon, &sigma;

    /// elasticity parameters
    double &E0;

    /// energy
    double &ener;

    /// research direction
    TMatGen_nvec &H0;

    LocalData(const TUmatData &umat_data) : sigma(umat_data.stress),
                                            E0(umat_data.props[0]),
                                            H0(umat_data.ddsdde),
                                            ener(umat_data.sse)
    {
        epsilon = umat_data.stran + umat_data.dstran;
//         PRINT(epsilon);
    }

}; // struct LocalData

} // namespace umat_c_elasticity_1d

#endif // #ifndef elasticity_1D_local_data_h
