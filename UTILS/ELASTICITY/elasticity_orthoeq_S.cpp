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

#ifndef elasticity_orthoeq_S_cpp
#define elasticity_orthoeq_S_cpp

#include "containers/vec_mat_tools.h"

#include "elasticity_orthotrope_S.cpp"

template<unsigned int ndim>
inline LMT::Mat<double, LMT::Gen<get_nb_vec(ndim)>, LMT::Dense<LMT::Row> > elasticity_orthoeq_S(
    const double &E,
    const double &N,
    const double &G)
{
    if (ndim == 2)
    {
        return elasticity_orthotrope_2D_S(E, E, N, G);
    }
    else if (ndim == 3)
    {
        return elasticity_orthotrope_3D_S(E, E, E, N, N, N, G, G, G);
    }
    else
    {
        assert(0);
    }

} // elasticity_orthoeq_S

#endif // #ifndef elasticity_orthoeq_S_cpp
