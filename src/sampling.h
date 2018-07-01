/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
#ifndef _ROME_SAMPLING
#define _ROME_SAMPLING

#include <stdio.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <assert.h>
#include <fstream>
#include "mkl.h"

#include "./mpi.h"
#include "./math.h"
#include "./symmetries.h"
#include "./metadata.h"
#include "./array_vector.h"

typedef std::vector< std::vector<FDOUBLE> > Vector2d;
typedef std::vector<FDOUBLE> 				Vector1d;

// turn off sampling 3D
#define SAMPLING3D

#ifdef SAMPLING3D
#include "../Healpix_2.15a/healpix_base.h"
#endif

#define NOPRIOR 0
#define PRIOR_ROTTILT_PSI 1

class HealpixSampler
{
    // The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
    static inline void checkDirection(FDOUBLE& rot,FDOUBLE& tilt){
        // The following was incorrect?!
        if (tilt < 0.)
        {
            tilt = -tilt;
            rot += 180.;
        }
        bool is_ok = false;
        while (!is_ok)
        {
            if (rot > 180.)
                rot -= 360.;
            else if (rot < -180.)
                rot += 360.;
            else
                is_ok = true;
        }
    }
    
public:
    // translation
    FDOUBLE offset_step,offset_range;
    std::vector<FDOUBLE> translations_x;
    std::vector<FDOUBLE> translations_y;
    // shift image twice to reduce memory traffic
    std::vector<std::pair<int, int>> single_translations_index;
    std::vector<FDOUBLE> single_translations;
    int nrTrans;
    // rotation
    FDOUBLE psi_step;
    std::vector<FDOUBLE> psi_angles;
    int nrPsi;
    // 3D direction
    // vector with sampling points described by angles
    std::vector<FDOUBLE> rot_angles;
    std::vector<FDOUBLE> tilt_angles;
    // vector with the original pixel number in the healpix object
    std::vector<int> directions_ipix;
    // adaptive_oversampling = 1,2,3,4,5
#ifdef SAMPLING3D
    Healpix_Base healpix_base;
#endif
    int healpix_order;
    bool is_3D;
    int nrPix;
    std::string fn_sym;
    // List of symmetry operators
    std::vector <Matrix2D<FDOUBLE> > R_repository, L_repository;
    // Two numbers that describe the symmetry group
    int pgGroup,pgOrder;
public:
    //----------- Random perturbation ---------
    FDOUBLE random_perturbation;
    FDOUBLE perturbation_factor;
    
    HealpixSampler();
    ~HealpixSampler();
    
    void initialize(FDOUBLE _offset_step = -1.,FDOUBLE _offset_range = -1.,FDOUBLE _psi_step = -1,
                    int _order = -1,std::string _fn_sym = "C1");
    
    void resetRandomlyPerturbedSampling();
    
    size_t NrDir(int oversampling_order = 0) const;
    
    size_t NrPsi(int oversampling_order = 0) const;
    
    size_t NrTrans(int oversampling_order = 0) const;
    
    size_t NrPoints(int oversampling_order) const;
    
    double getAngularSampling(int adaptive_oversampling = 0) const;
    
    // How often is each orientation oversampled?
    size_t oversamplingFactorOrientations(int oversampling_order) const;
    
    // How often is each translation oversampled?
    int oversamplingFactorTranslations(int oversampling_order) const;

    void getTranslation(int itrans, FDOUBLE &trans_x,FDOUBLE &trans_y) const
    {
        trans_x = translations_x[itrans];
        trans_y = translations_y[itrans];
    }
    
    void getTranslations(int itrans,int oversampling_order,FDOUBLE* over_trans_x,FDOUBLE* over_trans_y) const;
    
    //
    void getAllTranslationsAndOverTrans(int oversampling_order,Vector1d& trans_x,Vector1d& trans_y,
                                        Vector1d& trans_x_over,Vector1d& trans_y_over) const;
    
    //
    void getOrientations(int idir,int ipsi,int oversampling_order,FDOUBLE* over_psi,
                         FDOUBLE* over_rot = nullptr,FDOUBLE* over_tilt = nullptr) const;
    
    //
    void getAllTranslations(int oversampling_order, FDOUBLE* my_translations_x, FDOUBLE* my_translations_y) const;
    
    //
    void getAllTranslations(int oversampling_order,
                            Vector2d& my_translations_x,
                            Vector2d& my_translations_y) const;
    
    // 2D
    void getAllOrientations(int oversampling_order,FDOUBLE* my_psi) const;
    
    // 3D
    void getAllOrientations(int oversampling_order,
                            Vector2d& my_psi,
                            Vector2d& my_rot,
                            Vector2d& my_tilt) const;
    
    // Eliminate symmetry-equivalent points from the sampling_points_vector and sampling_points_angles vectors
    // This function first calls removeSymmetryEquivalentPointsGeometric,
    // and then checks each point versus all others to calculate an angular distance
    // If this distance is less than 0.8 times the angular sampling, the point is deleted
    // This cares care of sampling points near the edge of the geometrical considerations
    void removeSymmetryEquivalentPoints(FDOUBLE max_ang);
    
    // eliminate symmetry-related points based on simple geometrical considerations,
    // symmetry group, symmetry order
    void removeSymmetryEquivalentPointsGeometric(const int symmetry, int sym_order,
                                                 std::vector <Matrix1D<FDOUBLE> >  &sampling_points_vector);
    
    // Write all orientations as a sphere in a bild file
    // Mainly useful for debugging
    void writeAllOrientationsToBild(std::string fn_bild, std::string rgb, FDOUBLE size);
    
    // Write a BILD file describing the angular distribution
    // R determines the radius of the sphere on which cylinders will be placed
    // Rmax_frac determines the length of the longest cylinder (relative to R, 0.2 + +20%)
    // width_frac determines how broad each cylinder is. frac=1 means they touch each other
    void writeBildFileOrientationalDistribution(VectorOfFDOUBLE &pdf_direction,std::string &fn_bild, FDOUBLE R,
                                                FDOUBLE offset = 0., FDOUBLE Rmax_frac = 0.3, FDOUBLE width_frac = 0.5);
    
    //
    void writeOutSampling(std::string fn_sampling);
    
    //
    void readFromSampling(std::string fn_sampling);
};


class SamplingGrid
{
public:
    Vector2d exp_over_trans_x;
    Vector2d exp_over_trans_y;
    Vector2d exp_over_psi;
    Vector2d exp_over_tilt;
    Vector2d exp_over_rot;
    Vector1d exp_trans_x;
    Vector1d exp_trans_y;
    Vector1d exp_trans_x_over;
    Vector1d exp_trans_y_over;
    // single translation(~exp_trans_x[itrans] = exp_single_trans[exp_single_trans_index[itrans].x])
    Vector1d exp_single_trans;
    std::vector<std::pair<int, int>> exp_single_trans_index;
    //
    int exp_nr_trans,exp_nr_dir,exp_nr_psi,exp_nr_over_trans,exp_nr_over_rot;
    int adaptive_oversampling,current_oversampling;
    SamplingGrid():	exp_nr_trans(0),exp_nr_dir(0),exp_nr_psi(0),exp_nr_over_trans(0),exp_nr_over_rot(0),
    adaptive_oversampling(0),current_oversampling(0){}
    ~SamplingGrid(){}
    void initialize(HealpixSampler& sampler3d,int _adaptive_oversampling);
    void finalize();
    void computeGrid(HealpixSampler& sampler3d,int _current_oversampling);
    inline void getProjectMatrix(FDOUBLE A[][3],int idir,int ipsi,int iover_rot)
    {
        assert(iover_rot < exp_nr_over_rot);
        assert(idir < exp_nr_dir);
        assert(ipsi < exp_nr_psi);
        int iorient = idir*exp_nr_psi+ipsi;
        Euler_angles2matrix(exp_over_rot[iorient][iover_rot],exp_over_tilt[iorient][iover_rot],
                            exp_over_psi[iorient][iover_rot], A);
    }
    inline void getShiftxy(FDOUBLE& shiftx,FDOUBLE& shifty,int itrans,int iover_trans)
    {
        assert(itrans < exp_nr_trans);
        assert(iover_trans < exp_nr_over_trans);
        shiftx = exp_over_trans_x[itrans][iover_trans];
        shifty = exp_over_trans_y[itrans][iover_trans];
    }
    inline void getShiftxy(FDOUBLE& shiftx,FDOUBLE& shifty,
                           FDOUBLE& shiftxOver,FDOUBLE& shiftyOver,
                           int itrans,int iover_trans)
    {
        assert(itrans < exp_nr_trans);
        assert(iover_trans < exp_nr_over_trans);
        shiftx = exp_trans_x[itrans];
        shifty = exp_trans_y[itrans];
        shiftxOver = exp_trans_x_over[iover_trans];
        shiftyOver = exp_trans_y_over[iover_trans];
    }
    void testGetShift();
};

// info is var_num*3 matrix each row store basic info(start,end,order(indexBits))
// notice ,if we declare an pointer out of this function and then allocate the memmory for X in this function
// the allocated memory will be free....
namespace GTMSampling {
    
    double* uniformSampling(const double *infos,int var_num,int &K,int &L);
    
    double* sample1L(const double *infos,int &K,int &L);
    
    double* sample2L(const double *infos,int &K,int &L,bool inv = false);
    
    double* sample3L(const double *infos,int &K,int &L,bool inv = false);
    
    double* sample4L(const double *infos,int &K,int &L,bool inv = false);
    
    double* sample5L(const double *infos,int &K,int &L,bool inv = false);
    
    // double *sampleHealpix(int order,int &K,int &L);
    
}


#endif /* defined(_ROME_SAMPLING) */
