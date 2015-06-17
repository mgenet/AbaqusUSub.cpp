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

#ifndef elasticity_orthotrope_H_cpp
#define elasticity_orthotrope_H_cpp

inline LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > elasticity_orthotrope_2D_H(
    const double &E1,
    const double &E2,
    const double &N12,
    const double &G12)
{
    LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > H;
    H.set(0.);

    H(0,0) = pow(E1, 2) / (E1 - pow(N12, 2) * E2);
    H(1,1) = (E1 * E2)  / (E1 - pow(N12, 2) * E2);

    H(2,2) = 2 * G12;

    H(0,1) = (N12 * E1 * E2) / (E1 - pow(N12, 2) * E2); H(1,0) = H(0,1);

    return H;

} // elasticity_orthotrope_2D_H

inline LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > elasticity_orthotrope_2D_H(
    const LMT::Vec<double, 2> &E,
    const LMT::Vec<double, 1> &N,
    const LMT::Vec<double, 1> &G)
{
    return elasticity_orthotrope_2D_H(E[0], E[1], N[0], G[0]);

} // elasticity_orthotrope_2D_H

inline LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > elasticity_orthotrope_3D_H(
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
    LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > H;
    H.set(0.);

    H(0,0) = (-E2 + N23 * N23 * E3) * E1 * E1 / (-E1 * E2 + E1 * N23 * N23 * E3 + N12 * N12 * E2 * E2 + 2 * N12 * E2 * N23 * N13 * E3 + N13 * N13 * E2 * E3);
    H(1,1) = (-E1 + N13 * N13 * E3) * E2 * E2 / (-E1 * E2 + E1 * N23 * N23 * E3 + N12 * N12 * E2 * E2 + 2 * N12 * E2 * N23 * N13 * E3 + N13 * N13 * E2 * E3);
    H(2,2) = -(E1 - N12 * N12 * E2) * E2 * E3 / (-E1 * E2 + E1 * N23 * N23 * E3 + N12 * N12 * E2 * E2 + 2 * N12 * E2 * N23 * N13 * E3 + N13 * N13 * E2 * E3);

    H(3,3) = 2 * G12;
    H(4,4) = 2 * G13;
    H(5,5) = 2 * G23;

    H(0,1) = -(N12 * E2 + N23 * N13 * E3) * E1 * E2 / (-E1 * E2 + E1 * N23 * N23 * E3 + N12 * N12 * E2 * E2 + 2 * N12 * E2 * N23 * N13 * E3 + N13 * N13 * E2 * E3);
    H(0,2) = -(N12 * N23 + N13) * E2 * E1 * E3 / (-E1 * E2 + E1 * N23 * N23 * E3 + N12 * N12 * E2 * E2 + 2 * N12 * E2 * N23 * N13 * E3 + N13 * N13 * E2 * E3);
    H(1,2) = -(N23 * E1 + N13 * N12 * E2) * E2 * E3 / (-E1 * E2 + E1 * N23 * N23 * E3 + N12 * N12 * E2 * E2 + 2 * N12 * E2 * N23 * N13 * E3 + N13 * N13 * E2 * E3);

    H(1,0) = H(0,1);
    H(2,0) = H(0,2);
    H(2,1) = H(1,2);

    return H;

} // elasticity_orthotrope_3D_H

inline LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > elasticity_orthotrope_3D_H(
    const LMT::Vec<double, 3> &E,
    const LMT::Vec<double, 3> &N,
    const LMT::Vec<double, 3> &G)
{
    return elasticity_orthotrope_3D_H(E[0], E[1], E[2], N[0], N[1], N[2], G[0], G[1], G[2]);

} // elasticity_orthotrope_3D_H

inline LMT::Mat<double, LMT::Gen<3>, LMT::Dense<LMT::Row> > elasticity_orthotrope_H(
    const double &E1,
    const double &E2,
    const double &N12,
    const double &G12)
{
    return elasticity_orthotrope_2D_H(E1, E2, N12, G12);

} // elasticity_orthotrope_H

inline LMT::Mat<double, LMT::Gen<6>, LMT::Dense<LMT::Row> > elasticity_orthotrope_H(
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
    return elasticity_orthotrope_3D_H(E1, E2, E3, N12, N13, N23, G12, G13, G23);

} // elasticity_orthotrope_H

#endif // #ifndef elasticity_orthotrope_H_cpp
