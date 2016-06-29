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

#ifndef elasticity_isotrope_S_cpp
#define elasticity_isotrope_S_cpp

#include "containers/vec_mat_tools.h"

template<unsigned int ndim>
inline LMT::Mat<double, LMT::Gen<get_nb_vec(ndim)>, LMT::Dense<LMT::Row> > elasticity_isotrope_S(
    const double &E,
    const double &N,
    const std::string &type_stress_2D="plane_stress")
{
    if ((ndim == 3) or ((ndim == 2) and (type_stress_2D == "plane_stress")))
    {
        return ((1+N)/E) * self_sym_tens_prod_vec_col(I2_vec<2>()) - (N/E) * self_tens_prod_vec_col(I2_vec<2>());
    }
    else
    {
        assert (0);
    }

} // elasticity_isotrope_S

#endif // #ifndef elasticity_isotrope_S_cpp
