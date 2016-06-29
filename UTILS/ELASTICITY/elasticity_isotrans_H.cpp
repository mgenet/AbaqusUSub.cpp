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

#ifndef elasticity_isotrans_H_cpp
#define elasticity_isotrans_H_cpp

#include "elasticity_orthotrope_H.cpp"

template<unsigned int ndim>
inline LMT::Mat<double, LMT::Gen<get_nb_vec(ndim)>, LMT::Dense<LMT::Row> > elasticity_isotrans_H(
    const double &E_out,
    const double &E_in,
    const double &N_out,
    const double &N_in,
    const double &G_out,
    const double &G_in)
{
    if (ndim == 3)
    {
        return elasticity_orthotrope_3D_H(E_out, E_in, E_in, N_out, N_out, N_in, G_out, G_out, G_in);
    }
    else
    {
        assert(0);
    }

} // elasticity_isotrans_H

#endif // #ifndef elasticity_isotrans_H_cpp
