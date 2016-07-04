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

#ifndef umat_data_h
#define umat_data_h

#include "containers/vec.h"
#include "containers/mat.h"

#include "containers/vec_mat_tools.h"

/**
 *
*/
template<unsigned int ndim,
         unsigned int nvec,
         unsigned int npro,
         unsigned int nsta>
struct UMatData
{
    typedef LMT::Vec<double      > TVec;
    typedef LMT::Vec<double,   2 > TVec_2;
    typedef LMT::Vec<double,   3 > TVec_3;
    typedef LMT::Vec<double, ndim> TVec_ndim;
    typedef LMT::Vec<double, nvec> TVec_nvec;
    typedef LMT::Vec<double, nsta> TVec_nsta;
    typedef LMT::Vec<double, npro> TVec_npro;

    typedef LMT::Vec<TVec_nvec> TVec_TVec_nvec;

    typedef LMT::Mat<double, LMT::Gen<ndim> , LMT::Dense<LMT::Row> > TMatGen_ndim; // Need Row storage because of fortran
    typedef LMT::Mat<double, LMT::Sym<ndim> , LMT::Dense<LMT::Row> > TMatSym_ndim; // Need Row storage because of fortran
    typedef LMT::Mat<double, LMT::Diag<ndim>, LMT::Dense<LMT::Row> > TMatDia_ndim; // Need Row storage because of fortran
    typedef LMT::Mat<double, LMT::Gen<nvec> , LMT::Dense<LMT::Row> > TMatGen_nvec; // Need Row storage because of fortran
    typedef LMT::Mat<double, LMT::Sym<nvec> , LMT::Dense<LMT::Row> > TMatSym_nvec; // Need Row storage because of fortran
    typedef LMT::Mat<double, LMT::Diag<nvec>, LMT::Dense<LMT::Row> > TMatDia_nvec; // Need Row storage because of fortran

    char         &cmname; /// char stays char, because we don't know its size, so we can't reinterpret it as a string
    int          &noel;
    int          &npt;
    TVec_ndim    &coords;

    int          &nprops;
    TVec_npro    &props;

    TVec_2       &time;
    double       &dtime;
    double       &pnewdt;
    int          &kstep;
    int          &kinc;

    int          &nstatv;
    TVec_nsta    &statev;

    int          &ntens;
    TMatGen_ndim &dfgrd0;
    TMatGen_ndim &dfgrd1;
    TMatGen_ndim &drot;
    TVec_nvec    &stran;
    TVec_nvec    &dstran;
    TMatGen_nvec &ddsdde;
    TVec_nvec    &stress;

    double &temp;
    double &dtemp;

    double &sse;
    double &spd;
    double &scd;

    /// constructor
    UMatData(
        double *stress_,
        double *statev_,
        double *ddsdde_,
        double *sse_,
        double *spd_,
        double *scd_,
        double *rpl_,
        double *ddsddt_,
        double *drplde_,
        double *drpldt_,
        double *stran_,
        double *dstran_,
        double *time_,
        double *dtime_,
        double *temp_,
        double *dtemp_,
        double *predef_,
        double *dpred_,
        char   *cmname_,
        int    *ndi_,
        int    *nshr_,
        int    *ntens_,
        int    *nstatv_,
        double *props_,
        int    *nprops_,
        double *coords_,
        double *drot_,
        double *pnewdt_,
        double *celent_,
        double *dfgrd0_,
        double *dfgrd1_,
        int    *noel_,
        int    *npt_,
        int    *layer_,
        int    *kspt_,
        int    *kstep_,
        int    *kinc_,
        short   cmname_len_) : cmname(*cmname_),
                               noel(*noel_),
                               npt(*npt_),
                               coords(*reinterpret_cast<TVec_ndim*> (coords_)),
                               nprops(*nprops_),
                               props(*reinterpret_cast<TVec_npro*> (props_)),
                               time(*reinterpret_cast<TVec_2*> (time_)),
                               dtime(*dtime_),
                               pnewdt(*pnewdt_),
                               kstep(*kstep_),
                               kinc(*kinc_),
                               nstatv(*nstatv_),
                               statev(*reinterpret_cast<TVec_nsta*> (statev_)),
                               ntens(*ntens_),
                               dfgrd0(*reinterpret_cast<TMatGen_ndim*> (dfgrd0_)),
                               dfgrd1(*reinterpret_cast<TMatGen_ndim*> (dfgrd1_)),
                               drot(*reinterpret_cast<TMatGen_ndim*> (drot_)),
                               stran(*reinterpret_cast<TVec_nvec*> (stran_)),
                               dstran(*reinterpret_cast<TVec_nvec*> (dstran_)),
                               ddsdde(*reinterpret_cast<TMatGen_nvec*> (ddsdde_)),
                               stress(*reinterpret_cast<TVec_nvec*> (stress_)),
                               temp(*temp_),
                               dtemp(*dtemp_),
                               sse(*sse_),
                               spd(*spd_),
                               scd(*scd_)
    {
        rem_sqr_vec_col(stran);
        rem_sqr_vec_col(dstran);
        add_sqr_vec_col(stress);
    }

    /// destructor
    ~UMatData()
    {
        add_sqr_vec_col(stran);
        add_sqr_vec_col(dstran);
        rem_sqr_vec_col(stress);
        rem_sqr_mat_col_gen(ddsdde);
    }

    /// print
    void print() const
    {
//         std::cout << "name: "                 << &cmname << std::endl;
        std::cout << "element number: "       << noel    << std::endl;
        std::cout << "integration point number: "   << npt     << std::endl;
//         std::cout << "integration point position: " << coords  << std::endl;

//         std::cout << "number of material constants: "       << nprops       << std::endl;
//         std::cout << "size of list of material constants: " << props.size() << std::endl;
//         std::cout << "material constants: "                 << props        << std::endl;

        std::cout << "time: "             << time   << std::endl;
//         std::cout << "dtime: "            << dtime  << std::endl;
//         std::cout << "pnewdt: "           << pnewdt << std::endl;
        std::cout << "step number: "      << kstep  << std::endl;
        std::cout << "increment number: " << kinc   << std::endl;

//         std::cout << "number of state variables: "       << nstatv        << std::endl;
//         std::cout << "size of list of state variables: " << statev.size() << std::endl;
//         std::cout << "state variables: "                 << statev        << std::endl;

//         std::cout << "ntens: "  << ntens  << std::endl;
//         std::cout << "dfgrd0: " << dfgrd0 << std::endl;
//         std::cout << "dfgrd1: " << dfgrd1 << std::endl;
//         std::cout << "stran: "  << stran  << std::endl;
//         std::cout << "dstran: " << dstran << std::endl;
//         std::cout << "ddsdde: " << ddsdde << std::endl;
//         std::cout << "stress: " << stress << std::endl;

//         std::cout << "temp: "  << temp  << std::endl;
//         std::cout << "dtemp: " << dtemp << std::endl;

//         std::cout << "specific strain energy: "       << sse << std::endl;
//         std::cout << "specific plastic dissipation: " << spd << std::endl;
//         std::cout << "specific creep dissipation: "   << scd << std::endl;
    }
};

#endif
