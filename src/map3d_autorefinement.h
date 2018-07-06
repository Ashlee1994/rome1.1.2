/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Bevin Brett
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

#ifndef MAP3D_AUTOREFINEMENT_H_
#define MAP3D_AUTOREFINEMENT_H_

#include "./map3d_optimizer_new.h"

#define MAX_NR_ITER_WO_RESOL_GAIN 1
#define MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES 1
//
// --- node assignment example
// node: (0,1,2,3,4), (5,6,7,8),   9
//        half_one ,  half_two,  master
#define IF_HALF_ONE_MASTER_NODE if(node == 0)
#define IF_HALF_TWO_MASTER_NODE if(split_half_world_node == 0 && node != 0)
#define IF_HALF_ONE_TWO_MASTER_NODE if(split_half_world_node == 0)
#define IF_MASTER_NODE if(node == nodes-1)
#define IF_ALL_MASTER_NODE if(node == nodes-1 || split_half_world_node == 0)
#define IF_HALF_ONE_TWO_NODE if(node != nodes-1)
#define MPI_HALF_ONE_MASTER_NODE 0
#define MPI_HALF_TWO_MASTER_NODE (nodes-1)/2
//
namespace Map3dAutoRefinement
{
    using namespace Map3dOptimizer_new;
    
    // Global parameters to store accuracy on rot and trans
    extern double acc_rot;
    extern double acc_trans;
    
    // HiddenVariableMonitor hiddenVarMonitor
#define SEP
#define ELT(T,N,V) static auto& N = hiddenVarMonitor.N;
    MAPOPTIMIZER_OLD_HIDDEN_VAR_ELT
#undef ELT
#undef SEP
    
#define SEP
#define ELT(T,N,V) extern T N;
    AUTOREFINE_BASE_VARS
#undef SEP
#undef ELT
    //
    void setupMLoptimizer();
    
    void prepare();
    
    void iterate(bool verb = false);
    
    void maximizationBcastData();
    
    void destroyMLoptimizer();
    
    // Adjust angular sampling based on the expected angular accuracies for auto-refine procedure
    void updateAngularSampling(bool verb);
    
    // Calculate expected error in orientational assignments
    // Based on comparing projections of the model and see how many degrees apart gives rise to difference of power > 3*sigma^ of the noise
    void calculateExpectedAngularErrors(const MetaDataTable& exp_metadata,int my_first_particle, int my_last_particle);
    
    // Check convergence for auto-refine procedure
    void checkConvergence();
    
    // Print convergence information to screen for auto-refine procedure
    void printConvergenceStats();
    
    // Join two independent reconstructions ate the lowest frequencies to avoid convergence in distinct orientations
    void joinTwoHalvesAtLowResolution();
    
    // Join the sums from two random halves
    void combineWeightedSumsTwoRandomHalves();
    
    // When refining two random halves separately, the master receives both models, calculates FSC and the power of their difference
    // and sends these curves, together with new tau2_class estimates to all slaves...
    void compareTwoHalves();
    
    // Write temporary data and weight arrays from the backprojector to disc to allow unregularized reconstructions
    void writeTemporaryDataAndWeightArrays();
    
    //  Read temporary data and weight arrays from disc and perform unregularized reconstructions
    //  Also write the unregularized reconstructions to disc.
    void readTemporaryDataAndWeightArraysAndReconstruct();
} // end namespace Map3dAutoRefinement

#endif /* defined(MAP3D_AUTOREFINEMENT_H_) */
