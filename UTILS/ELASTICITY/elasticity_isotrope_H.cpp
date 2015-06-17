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

#ifndef elasticity_isotrope_H_cpp
#define elasticity_isotrope_H_cpp

#include "containers/vec_mat_tools.h"

template<unsigned int ndim>
inline LMT::Mat<double, LMT::Gen<get_nb_vec(ndim)>, LMT::Dense<LMT::Row> > elasticity_isotrope_H(
    const double      &E,
    const double      &N,
    const std::string &type_stress_2D="plane_stress")
{
    double l,m;

    if ((ndim == 2) and (type_stress_2D == "plane_stress"))
    {
        l = (E * N) / ((1. + N) * (1. - N));
        m = E / 2. / (1. + N);
    }
    else if ((ndim == 3) or  ((ndim == 2) and (type_stress_2D == "plane_strain")))
    {
        l = (E * N) / ((1. + N) * (1. - 2. * N));
        m = E / 2. / (1. + N);
    }
    else
    {
        assert (0);
    }

    return l * self_tens_prod_vec_col(I2_vec<ndim>()) + 2 * m * self_sym_tens_prod_vec_col(I2_vec<ndim>());

} // elasticity_isotrope_H

#endif // #ifndef elasticity_isotrope_H_cpp
