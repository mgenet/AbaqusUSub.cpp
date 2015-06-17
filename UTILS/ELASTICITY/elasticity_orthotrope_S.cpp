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

#ifndef elasticity_orthotrope_S_cpp
#define elasticity_orthotrope_S_cpp

inline LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > elasticity_orthotrope_2D_S(
    const double &E1,
    const double &E2,
    const double &N12,
    const double &G12)
{
    LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > S;
    S.set(0.);

    S(0,0) = 1./E1;
    S(1,1) = 1./E2;

    S(2,2) = 1./G12;

    S(0,1) = -N12/E1; S(1,0) = S(0,1);

    return S;

} // elasticity_orthotrope_2D_S

inline LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > elasticity_orthotrope_2D_S(
    const LMT::Vec<double, 2> &E,
    const LMT::Vec<double, 1> &N,
    const LMT::Vec<double, 1> &G)
{
    return elasticity_orthotrope_2D_S(E[0], E[1], N[0], G[0]);

} // elasticity_orthotrope_2D_S

inline LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > elasticity_orthotrope_3D_S(
    const double &E1,
    const double &E2,
    const double &E3,
    const double &N12,
    const double &N13,
    const double &N23,
    const double &G12,
    const double &G13,
    const double &G23)
{
    LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > S;
    S.set(0.);

    S(0,0) = 1./E1;
    S(1,1) = 1./E2;
    S(2,2) = 1./E3;

    S(3,3) = 1./2./G12;
    S(4,4) = 1./2./G13;
    S(5,5) = 1./2./G23;

    S(0,1) = -N12/E1; S(1,0) = S(0,1);
    S(2,0) = -N13/E1; S(0,2) = S(2,0);
    S(1,2) = -N23/E2; S(2,1) = S(1,2);

    return S;

} // elasticity_orthotrope_3D_S

inline LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > elasticity_orthotrope_3D_S(
    const LMT::Vec<double, 3> &E,
    const LMT::Vec<double, 3> &N,
    const LMT::Vec<double, 3> &G)
{
    return elasticity_orthotrope_3D_S(E[0], E[1], E[2], N[0], N[1], N[2], G[0], G[1], G[2]);
} // elasticity_orthotrope_3D_S

inline LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > elasticity_orthotrope_S(
    const double &E1,
    const double &E2,
    const double &N12,
    const double &G12)
{
    return elasticity_orthotrope_2D_S(E1, E2, N12, G12);
} // elasticity_orthotrope_S

inline LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > elasticity_orthotrope_S(
    const double &E1,
    const double &E2,
    const double &E3,
    const double &N12,
    const double &N13,
    const double &N23,
    const double &G12,
    const double &G13,
    const double &G23)
{
    return elasticity_orthotrope_3D_S(E1, E2, E3, N12, N13, N23, G12, G13, G23);
} // elasticity_orthotrope_S

#endif // #ifndef elasticity_orthotrope_S_cpp
