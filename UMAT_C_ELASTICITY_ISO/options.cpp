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

#ifndef elasticity_iso_options_cpp
#define elasticity_iso_options_cpp

namespace umat_c_elasticity_iso {

constexpr unsigned int ndim = SPACE_NDIM; // can be SPACE_NDIM, 2 or 3

constexpr unsigned int nvec = ndim*(ndim+1)/2;
constexpr unsigned int npro = 2;
constexpr unsigned int nsta = 0;

} // namespace umat_c_elasticity_iso

#endif // #ifndef elasticity_iso_options_cpp
