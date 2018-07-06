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

#include "./map3d_autorefinement.h"

namespace Map3dAutoRefinement
{
    using namespace Map3dOptimizer_new;
    
    // Global parameters to store accuracy on rot and trans
    double acc_rot;
    double acc_trans;
    
#define SEP
#define ELT(T,N,V) T N;
    AUTOREFINE_BASE_VARS
#undef SEP
#undef ELT
    
#ifdef USEMPI
    MPI::Intracomm split_half_world;
#endif
    int split_half_world_node = -1;
    
void setupMLoptimizer()
{
    Map3dOptimizer_new::setupMLoptimizer();
#ifdef USEMPI
    split_half_world = getNewIntracomm();
    IF_HALF_ONE_TWO_NODE split_half_world_node = split_half_world.Get_rank();
#endif
}

void prepare()
{
    Map3dOptimizer_new::prepare();
}

void iterate(bool verb)
{
#ifdef USEMPI
    if (do_split_random_halves && node == nodes-1){
        assert(split_half_world_node == -1);
    }
#endif
    // Update the current resolution and image sizes, and precalculate resolution pointers
    // The rest of the time this will be done after maximization and before writing output files,
    // so that current resolution is in the output files of the current iteration
    IF_HALF_ONE_TWO_NODE {
        bool set_by_ini_high = ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc) );
        mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
    }
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif
    
    // setupStatusTracer();
    
    // TODO
    // continue
    // if (continue_fn!="NULL") {
    //     readResult();
    //     // After the first iteration the references are always CTF-corrected
    //     if (do_ctf_correction)
    //         refs_are_ctf_corrected = true;
    // }
    
    bool has_already_reached_convergence = false;
    for (iter = iter + 1; iter <= nr_iter; iter++)
    {
        IF_HALF_ONE_MASTER_NODE std::cout<<"current_resolution = "<<mapModel.current_resolution<<std::endl;
        NODE0ONLY showPopulation();
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(iter, "iterate()_iter", __FILE__, __LINE__);
#endif
        
        if (do_auto_refine && node == 0)
            printConvergenceStats();
        
        const double starttime = dtime();
        // update coarse_size,current_size,Npix_per_shell,Mresol_coarse,Mresol_fine
        IF_HALF_ONE_TWO_NODE {
            double angularSampler = sampler3d.getAngularSampling();
            mapModel.updateImageSizeAndResolutionPointers(	Npix_per_shell, Mresol_coarse, Mresol_fine,
                                                          	coarse_size, current_size,adaptive_oversampling, angularSampler,
                                                          	mlModel.ave_Pmax, has_high_fsc_at_limit, do_use_all_data);
        }
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(current_size, "updateImageSizeAndResolutionPointers()_current_size", __FILE__, __LINE__);
        global_data_stream.foutInt(coarse_size, "updateImageSizeAndResolutionPointers()_coarse_size", __FILE__, __LINE__);
        global_data_stream.foutInt(Npix_per_shell.wptrAll(), (ori_size/2+1), "updateImageSizeAndResolutionPointers()_Npix_per_shell", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_fine.wptrAll(), current_size*(current_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_fine", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_coarse.wptrAll(), coarse_size*(coarse_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_coarse", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        IF_HALF_ONE_MASTER_NODE std::cout<<"current_size = "<<current_size<<" coarse_size = "<<coarse_size<<std::endl;
        
        // set the exp_metadata and exp_image_data and some other data
        // allocate the maximum data
        IF_HALF_ONE_TWO_NODE prepareExpMap();
        
        // Calculate expected error in orientational assignments
        if (!(iter==1 && do_firstiter_cc))
        {
            int n_trials_acc = 100;
            NODE0ONLY {
                // check error
                n_trials_acc = std::min(n_trials_acc, metadata.numberOfParticles(1));
                calculateExpectedAngularErrors(metadata, 0, n_trials_acc-1);
            }
            // The reconstructing slave Bcast acc_rottilt, acc_psi, acc_trans to all other nodes!
            MPI::COMM_WORLD.Bcast(&acc_rot, 1, MPI_FDOUBLE, 0);
            MPI::COMM_WORLD.Bcast(&acc_trans, 1, MPI_FDOUBLE, 0);
            //
#ifdef DATA_STREAM
            IF_HALF_ONE_TWO_NODE
            {
                global_data_stream.foutInt(&n_trials_acc, 1, "expectation()_n_trials_acc", __FILE__, __LINE__);
                for (int iclass = 0; iclass < nr_classes; iclass++)
                {
                    global_data_stream.foutDouble(mlModel.acc_rot.wptrAll()[iclass], "expectation()_mymodel_acc_rot_iclass", __FILE__, __LINE__);
                    global_data_stream.foutDouble(mlModel.acc_trans.wptrAll()[iclass], "expectation()_mymodel_acc_trans_iclass", __FILE__, __LINE__);
                    global_data_stream.foutDouble(mlModel.orientability_contrib[iclass].wptrAll(), mlModel.orientability_contrib[iclass].size(),
                                                  "expectation()_mymodel_orientability_contrib_iclass", __FILE__, __LINE__);
                    global_data_stream.check();global_data_stream.flush();
                }
                global_data_stream.foutDouble(acc_rot, "expectation()_acc_rot", __FILE__, __LINE__);
                global_data_stream.foutDouble(acc_trans, "expectation()_acc_trans", __FILE__, __LINE__);
                global_data_stream.check();global_data_stream.flush();
            }
#endif
        }
        
        IF_HALF_ONE_TWO_NODE
        {
            if(do_auto_refine && iter > 1)
            	updateAngularSampling(node==0);
        }
        
        IF_HALF_ONE_TWO_NODE prepareExpData();
        
//        if (iter == 12) {
//            //
//            MPI::COMM_WORLD.Barrier();
//            std::cerr<<"debug here exit...."<<std::endl;
//            EXIT_ABNORMALLY;
//        }
        
        IF_HALF_ONE_TWO_NODE expectation();
        
#ifdef USEMPI
        
        IF_HALF_ONE_TWO_NODE mlModel.reduceData(split_half_world);
        
        //#define DO_RECONSTRUCT_EACH_NODE
#ifdef DO_RECONSTRUCT_EACH_NODE
        IF_HALF_ONE_TWO_NODE mapModel.reduceData(split_half_world,false);
#else
        IF_HALF_ONE_TWO_NODE mapModel.reduceData(split_half_world); // TODO : may only need to reduce to master node.
#endif
        gatherMetaDataToMaster(metadata, nodes-1, do_split_random_halves);
        
        MPI::COMM_WORLD.Barrier();
        
#endif
        
        // Write out data and weight arrays to disc in order to also do an unregularized reconstruction ?????
        if (do_auto_refine && has_converged)
            writeTemporaryDataAndWeightArrays();
        
        // Inside iterative refinement: do FSC-calculation BEFORE the solvent flattening, otherwise over-estimation of resolution
        // anyway, now that this is done inside BPref, there would be no other way...
        if (do_split_random_halves)
        {
            // For asymmetric molecules, join 2 half-reconstructions at the lowest resolutions to prevent them from diverging orientations
            if (low_resol_join_halves > 0.)
                joinTwoHalvesAtLowResolution();
            
            // Calculate gold-standard FSC curve
            compareTwoHalves();

            // For automated sampling procedure
            IF_HALF_ONE_TWO_NODE // the master does not have the correct mymodel.current_size, it only handles metadata!
            {
                // Check that incr_size is at least the number of shells as between FSC=0.5 and FSC=0.143
                int fsc05   = -1;
                int fsc0143 = -1;
                auto fsc_halves_class_0class = mlModel.fsc_halves_class[0].rptr(ori_size/2+1);
                for (int i = 0; i < ori_size/2+1; i++)
                {
                    if (fsc_halves_class_0class[i] < 0.5 && fsc05 < 0)
                        fsc05 = i;
                    if (fsc_halves_class_0class[i] < 0.143 && fsc0143 < 0)
                        fsc0143 = i;
                }
                // At least fsc05 - fsc0143 + 5 shells as incr_size
                mapModel.incr_size = std::max(mapModel.incr_size, fsc0143 - fsc05 + 5);
                has_high_fsc_at_limit = (fsc_halves_class_0class[current_size/2 - 1] > 0.2);
            }
            
#ifdef DATA_STREAM
            global_data_stream.foutDouble(mapModel.incr_size, "iterate()_incr_size", __FILE__, __LINE__);
            global_data_stream.foutDouble(has_high_fsc_at_limit, "iterate()_has_high_fsc_at_limit", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            // Upon convergence join the two random halves
            if (do_join_random_halves || do_always_join_random_halves)
            {
                combineWeightedSumsTwoRandomHalves();
            }
        }
        
        //
        hiddenVarMonitor.reduceData();
        IF_HALF_ONE_TWO_NODE
        {
            IF_HALF_ONE_MASTER_NODE
            {
                maximization(do_split_random_halves,(do_join_random_halves || do_always_join_random_halves));
            }
            IF_HALF_TWO_MASTER_NODE
            {
                maximization(do_split_random_halves,do_join_random_halves);
            }
            maximizationBcastData();
        }
        hiddenVarMonitor.bcastData(0);// TODO : may not need
        
        // Also perform the unregularized reconstruction
        if (do_auto_refine && has_converged)
            readTemporaryDataAndWeightArraysAndReconstruct();
        
        // Make sure all nodes have the same resolution, set the data_vs_prior array from half1 also for half2
        if (do_split_random_halves)
        {
            std::vector<FDOUBLE> data_vs_prior_class;
            IF_HALF_ONE_TWO_NODE assert(mlModel.data_vs_prior_class[0].size()==ori_size/2+1);
            data_vs_prior_class.resize(ori_size/2+1);
            for (int iclass = 0; iclass < nr_classes; iclass++)
            {
                NODE0ONLY copy(data_vs_prior_class.data(), data_vs_prior_class.size(), mlModel.data_vs_prior_class[iclass].rptrAll(), mlModel.data_vs_prior_class[iclass].size());
                MPI::COMM_WORLD.Bcast(data_vs_prior_class.data(), data_vs_prior_class.size(), MPI_FDOUBLE, 0);
                IF_HALF_ONE_TWO_NODE copy(mlModel.data_vs_prior_class[iclass].wptrAll(), mlModel.data_vs_prior_class[iclass].size(), data_vs_prior_class.data(), data_vs_prior_class.size());
            }
        }
        
        // Apply masks to the reference images
        // At the last iteration, do not mask the map for validation purposes
        if(do_solvent && !has_converged){
            IF_HALF_ONE_TWO_NODE mapModel.applySolventFlatten(mask_fn);
        }
        
#ifdef DATA_STREAM
        IF_HALF_ONE_TWO_NODE
        {
            global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), mapModel.Irefs[0].dimzyx, "iterate()_Iref[0]", __FILE__, __LINE__);
            global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), mapModel.Irefs[nr_classes-1].dimzyx, "iterate()_Iref[nr_classes-1]", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
        }
#endif
        
        IF_HALF_ONE_TWO_NODE
        {
            // Re-calculate the current resolution, do this before writing to get the correct values in the output files
            bool set_by_ini_high = ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc) );
            mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
            mapModel.printResolution(mlModel.data_vs_prior_class,set_by_ini_high);
        }
        
#ifdef DATA_STREAM
        global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        // If we are joining random halves, then do not write an optimiser file so that it cannot be restarted!
        bool do_write_optimiser = !do_join_random_halves;
        // Write out final map without iteration number in the filename
        // if (do_join_random_halves)
        //     iter = -1;
        
        // Write output files
        IF_HALF_ONE_TWO_MASTER_NODE
        {
            bool do_write_sampling = false;bool do_write_data = false;bool do_write_optimiser = false;
            bool do_write_model = true;bool do_write_mrc = true;int random_subset = node==0?1:2;
           	writeResult(do_write_sampling, do_write_data, do_write_optimiser, do_write_model, do_write_mrc, random_subset);
        }
        IF_MASTER_NODE
        {
            bool do_write_sampling = true;bool do_write_data = true;
            bool do_write_model = false;bool do_write_mrc = false;int random_subset = -1;
           	writeResult(do_write_sampling, do_write_data, do_write_optimiser, do_write_model, do_write_mrc, random_subset);
        }
        
        if (do_auto_refine && has_converged)
        {
            if (verb > 0)
            {
                std::cout << " Auto-refine: Refinement has converged, stopping now... " << std::endl;
                // std::cout << " Auto-refine: + Final reconstruction from all particles is saved as: " <<  fn_out << "_class001.mrc" << std::endl;
                // std::cout << " Auto-refine: + Final model parameters are stored in: " << fn_out << "_model.star" << std::endl;
                // std::cout << " Auto-refine: + Final data parameters are stored in: " << fn_out << "_data.star" << std::endl;
                std::cout << " Auto-refine: + Final resolution (without masking) is: " << 1./mapModel.current_resolution << std::endl;
                if (acc_rot < 10.)
                    std::cout << " Auto-refine: + But you may want to run relion_postprocess to mask the unfil.mrc maps and calculate a higher resolution FSC" << std::endl;
                else
                {
                    std::cout << " Auto-refine: + WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles!" << std::endl;
                    std::cout << " Auto-refine: + WARNING: This has been observed to lead to spurious FSC curves, so be VERY wary of inflated resolution estimates..." << std::endl;
                    std::cout << " Auto-refine: + WARNING: You most probably do NOT want to publish these results!" << std::endl;
                    std::cout << " Auto-refine: + WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
                }
                // if (do_use_reconstruct_images)
                //     std::cout << " Auto-refine: + Used rlnReconstructImageName images for final reconstruction. Ignore filtered map, and only assess the unfiltered half-reconstructions!" << std::endl;
            }
            break;
        }
        
        // Show more stats
        NODE0ONLY PerformanceCounter::showAll(std::cout, false);
        NODE0ONLY PerformanceCounter::showAll(std::cerr, true);
        
        // update the metadata,free exp_metadata and exp_image_data
        IF_HALF_ONE_TWO_NODE endExpData();
        
        IF_HALF_ONE_TWO_NODE
        {
            if (do_auto_refine) checkConvergence();
        }
        
        MPI::COMM_WORLD.Bcast(&has_converged, 1, MPI::BOOL, 0);
        MPI::COMM_WORLD.Bcast(&do_join_random_halves, 1, MPI::BOOL, 0);
        
        const double endtime = dtime();
        NODE0ONLY std::cout<<"**** iteration "<<iter<<" completed in "<<endtime-starttime<<" seconds ****"<<std::endl<<std::flush;
    }
    
    NODE0ONLY showPopulation();
}

void destroyMLoptimizer()
{
    Map3dOptimizer_new::destroyMLoptimizer();
}

void maximizationBcastData()
{
#ifdef USEMPI
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        split_half_world.Bcast(mapModel.Irefs[iclass].wptr(), mapModel.Irefs[iclass].dimzyx, MPI_FDOUBLE, 0);
        split_half_world.Bcast(mlModel.data_vs_prior_class[iclass].wptrAll(), mlModel.data_vs_prior_class[iclass].size(), MPI_FDOUBLE, 0);
        split_half_world.Bcast(mlModel.sigma2_class[iclass].wptrAll(), mlModel.sigma2_class[iclass].size(), MPI_FDOUBLE, 0);
        // ??
        split_half_world.Bcast(mlModel.tau2_class[iclass].wptrAll(), mlModel.tau2_class[iclass].size(), MPI_FDOUBLE, 0);
    }
#endif
}
    
//
void setupAutorefine()
{
    if (!do_auto_refine) return;
    
    if (nr_classes > 1 && do_split_random_halves)
        ERROR_REPORT("ERROR: One cannot use --split_random_halves with more than 1 reference... You could first classify, and then refine each class separately using --random_halves.");
    
    if (do_join_random_halves && !do_split_random_halves)
        ERROR_REPORT("ERROR: cannot join random halves because they were not split in the previous run");
    
    if (do_always_join_random_halves)
        std::cout << " Joining half-reconstructions at each iteration: this is a developmental option to test sub-optimal FSC usage only! " << std::endl;
    
    // if (do_auto_refine)
    {
        nr_iter = 999;
        has_fine_enough_angular_sampling = false;
        
        if (iter == 0 && sampler3d_healpix_order >= autosampling_hporder_local_searches)
        {
            ERROR_REPORT("sampler3d_healpix_order shuld be smaller than autosampling_hporder_local_searches...");
            // mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
            // sampling.orientational_prior_mode = PRIOR_ROTTILT_PSI;
            // DOUBLE rottilt_step = sampling.getAngularSampling(adaptive_oversampling);
            // mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 2. * 2. * rottilt_step * rottilt_step;
        }
    }
    //
    // if (do_split_random_halves && node->size <= 2)
    //     REPORT_ERROR("MlOptimiserMpi::initialiseWorkLoad: at least 3 MPI processes are required when splitting data into random halves");
    // else if(node->size <= 1)
    //     REPORT_ERROR("MlOptimiserMpi::initialiseWorkLoad: at least 2 MPI processes are required, otherwise use the sequential program");
    
    assert(random_seed!=0);
    // Get the same random number generator seed for all mpi nodes
    // if (random_seed == -1)
    // {
    //     if (node->isMaster())
    //     {
    //         random_seed = time(NULL);
    //         for (int slave = 1; slave < node->size; slave++)
    //             node->relion_MPI_Send(&random_seed, 1, MPI_INT, slave, MPITAG_RANDOMSEED, MPI_COMM_WORLD);
    //     }
    //     else
    //     {
    //         MPI_Status status;
    //         node->relion_MPI_Recv(&random_seed, 1, MPI_INT, 0, MPITAG_RANDOMSEED, MPI_COMM_WORLD, status);
    //     }
    // }
    
    // First split the data into two random halves and then randomise the particle order
    // if (do_split_random_halves)
    //     mydata.divideOriginalParticlesInRandomHalves(random_seed);
    
    // Randomise the order of the particles
    // mydata.randomiseOriginalParticlesOrder(random_seed, do_split_random_halves);
    
    // Also randomize random-number-generator for perturbations on the angles
    // init_random_generator(random_seed);
    
    
    // if (node->isMaster())
    // {
    // The master never participates in any actual work
    //     my_first_ori_particle_id = 0;
    //     my_last_ori_particle_id = -1;
    // }
    // else
    // {
    if (do_split_random_halves)
    {
        // int nr_slaves_subset1 = (node->size - 1) / 2;
        // int nr_slaves_subset2 = nr_slaves_subset1;
        // if ( (node->size - 1) % 2 != 0)
        //     nr_slaves_subset1 += 1;
        // if (node->myRandomSubset() == 1)
        // {
        //     // Divide first half of the images
        //     divide_equally(mydata.numberOfOriginalParticles(1), nr_slaves_subset1, node->rank / 2, my_first_ori_particle_id, my_last_ori_particle_id);
        // }
        // else
        // {
        //     // Divide second half of the images
        //     divide_equally(mydata.numberOfOriginalParticles(2), nr_slaves_subset2, node->rank / 2 - 1, my_first_ori_particle_id, my_last_ori_particle_id);
        //     my_first_ori_particle_id += mydata.numberOfOriginalParticles(1);
        //     my_last_ori_particle_id += mydata.numberOfOriginalParticles(1);
        // }
        //
        // testing node1
        nr_global_images = metadata.numberOfParticles(1);
        nr_local_images = divide_equally_pool(nr_global_images, nr_pool, nodes, node, first_local_image, last_local_image);
    }
    else
    {
        // int nr_slaves = (node->size - 1);
        // divide_equally(mydata.numberOfOriginalParticles(), nr_slaves, node->rank - 1, my_first_ori_particle_id, my_last_ori_particle_id);
    }
    
    // }
    
}
    
void updateAngularSampling(bool verb)
{
//    if (iter >= 3)
//    {
//        has_fine_enough_angular_sampling = true;
//        return;
//    }
    
    if (!do_split_random_halves)
        ERROR_REPORT("Map3dAutoRefinement::updateAngularSampling: BUG! updating of angular sampling should only happen for gold-standard (auto-) refinements.");
    
    // Only change the sampling if the resolution has not improved during the last 2 iterations
    // AND the hidden variables have not changed during the last 2 iterations
    double old_rottilt_step = sampler3d.getAngularSampling(adaptive_oversampling);
    
    //
#ifdef DATA_STREAM
    global_data_stream.foutDouble(old_rottilt_step, "updateAngularSampling()_old_rottilt_step", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif

    // Only use a finer angular sampling is the angular accuracy is still above 75% of the estimated accuracy
    // If it is already below, nothing will change and eventually nr_iter_wo_resol_gain or nr_iter_wo_large_hidden_variable_changes will go above MAX_NR_ITER_WO_RESOL_GAIN
    if (nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN && nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES)
    {
        // Old rottilt step is already below 75% of estimated accuracy: have to stop refinement
        if (old_rottilt_step < 0.75 * acc_rot)
        {
            // don't change angular sampling, as it is already fine enough
            has_fine_enough_angular_sampling = true;
            
        }
        else
        {
            has_fine_enough_angular_sampling = false;
            
            // A. Use translational sampling as suggested by acc_trans
            
            // Prevent very coarse translational samplings: max 1.5
            // Also stay a bit on the safe side with the translational sampling: 75% of estimated accuracy
            double new_step = std::min(1.5, 0.75 * acc_trans) * std::pow(2., adaptive_oversampling);
            //
#ifdef DATA_STREAM
            global_data_stream.foutDouble(new_step, "updateAngularSampling()_new_step", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            // Search ranges are five times the last observed changes in offsets
            double new_range = 5. * current_changes_optimal_offsets;
            // New range can only become 30% bigger than the previous range (to prevent very slow iterations in the beginning)
            new_range = std::min(1.3*sampler3d.offset_range, new_range);
            // Prevent too narrow searches: always at least 3x3 pixels in the coarse search
            if (new_range < 1.5 * new_step)
                new_range = 1.5 * new_step;
            // Also prevent too wide searches: that will lead to memory problems:
            // If steps size < 1/4th of search range, then decrease search range by 50%
            if (new_range > 4. * new_step)
                new_range /= 2.;
            
            //If even that was not enough: use coarser step size and hope things will settle down later...
            if (new_range > 4. * new_step)
                new_step = new_range / 4.;
            sampler3d.setTranslations(new_step, new_range);
            //
#ifdef DATA_STREAM
            global_data_stream.foutDouble(new_step, "updateAngularSampling()_new_step", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            // B. Use twice as fine angular sampling
            int new_hp_order;
            double new_rottilt_step, new_psi_step;
            if (mapModel.ref_dim == 3)
            {
                new_hp_order = sampler3d.healpix_order + 1;
                new_rottilt_step = new_psi_step = 360. / (6 * round(std::pow(2., new_hp_order + adaptive_oversampling)));
            }
            else
            {
                ERROR_REPORT("Map3dAutoRefinement::updateAngularSampling BUG: ref_dim should be three");
            }

            // Set the new sampling in the sampling-object
            sampler3d.setOrientations(new_hp_order, new_psi_step * std::pow(2., adaptive_oversampling));
            //
#ifdef DATA_STREAM
            global_data_stream.foutInt(new_hp_order, "updateAngularSampling()_new_hp_order", __FILE__, __LINE__);
            global_data_stream.foutDouble(new_psi_step, "updateAngularSampling()_new_psi_step", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            // Resize the pdf_direction arrays to the correct size and fill with an even distribution
            mlModel.initialisePdfDirection(sampler3d.NrDir());
            //
#ifdef DATA_STREAM
            global_data_stream.foutInt(mlModel.nr_directions, "updateAngularSampling()_nr_directions", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            // Reset iteration counters
            nr_iter_wo_resol_gain = 0;
            nr_iter_wo_large_hidden_variable_changes = 0;

            // Reset smallest changes hidden variables
            smallest_changes_optimal_classes = 9999999;
            smallest_changes_optimal_offsets = 999.;
            smallest_changes_optimal_orientations = 999.;

            // If the angular sampling is smaller than autosampling_hporder_local_searches, then use local searches of +/- 6 times the angular sampling
            if (new_hp_order >= autosampling_hporder_local_searches)
            {
                // Switch ON local angular searches
                mlModel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
                sampler3d.orientational_prior_mode = PRIOR_ROTTILT_PSI;
                sigma2_angle = 2. * new_rottilt_step;
                mlModel.sigma2_rot = mlModel.sigma2_tilt = mlModel.sigma2_psi = 2. * 2. * new_rottilt_step * new_rottilt_step;
                do_local_searching = true;
            }
            //
#ifdef DATA_STREAM
            global_data_stream.foutInt(mlModel.orientational_prior_mode, "updateAngularSampling()_mymodel.orientational_prior_mode", __FILE__, __LINE__);
            global_data_stream.foutInt(sampler3d.orientational_prior_mode, "updateAngularSampling()_mymodel.orientational_prior_mode", __FILE__, __LINE__);
            global_data_stream.foutDouble(mlModel.sigma2_rot, "updateAngularSampling()_mymodel.sigma2_rot", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
        }
    }
    
    // Print to screen
    if (verb)
    {
        std::cout <<" Auto-refine: updateAngularSampling() "<<std::endl;
		std::cout << " Auto-refine: Angular step= " << sampler3d.getAngularSampling(adaptive_oversampling) << " degrees; local searches= ";
        if (sampler3d.orientational_prior_mode == NOPRIOR)
            std:: cout << "false" << std::endl;
        else
            std:: cout << "true" << std::endl;
        std::cout << " Auto-refine: Offset search range= " << sampler3d.offset_range << " pixels; offset step= " << sampler3d.getTranslationalSampling(adaptive_oversampling) << " pixels"<<std::endl;
        std::cout << " Auto-refine: offset step2= "<<sampler3d.offset_step<<",Healpix order= "<<sampler3d.healpix_order<<std::endl;
        std::cout << " Auto-refine: nr_iter_wo_resol_gain = "<<nr_iter_wo_resol_gain<<",nr_iter_wo_large_hidden_variable_changes = "<<nr_iter_wo_large_hidden_variable_changes<<std::endl;
        std::cout << " Auto-refine: (old_rottilt_step < 0.75 * acc_rot?) ("<<sampler3d.getAngularSampling(adaptive_oversampling)<<" < "<<0.75*acc_rot<<")??"<<std::endl;
        std::cout << " Auto-refine: --------------------------------------------------------------------------------- "<<std::endl;
    }
}

// Calculate expected error in orientational assignments
// Based on comparing projections of the model and see how many degrees apart gives rise to difference of power > 3*sigma^ of the noise
void calculateExpectedAngularErrors(const MetaDataTable& metadata,int my_first_particle, int my_last_particle)
{
    int n_trials = my_last_particle-my_first_particle+1;
    
    // Use smaller images in the first pass, but larger ones in the second pass
    int exp_current_image_size = current_size;
    
    // Separate angular error estimate for each of the classes
    acc_rot = acc_trans = 999.; // later XMIPP_MIN will be taken to find the best class...
    
    // P(X | X_1) / P(X | X_2) = exp ( |F_1 - F_2|^2 / (-2 sigma2) )
    // exp(-4.60517) = 0.01
    double pvalue = 4.60517;
    
    int exp_current_image_Fsize2 = exp_current_image_size*(exp_current_image_size/2+1);
    std::vector<double> Fctf(exp_current_image_Fsize2);
    auto Fctf_wptr = Fctf.data();
    std::vector<double> F1Real(exp_current_image_Fsize2),F1Imag(exp_current_image_Fsize2);
    std::vector<double> F2Real(exp_current_image_Fsize2),F2Imag(exp_current_image_Fsize2);
    auto F1Real_wptr = F1Real.data();auto F1Imag_wptr = F1Imag.data();
    auto F2Real_wptr = F2Real.data();auto F2Imag_wptr = F2Imag.data();
    double A1[3][3], A2[3][3];
    
    std::cout << " Estimating accuracies in the orientational assignment ... " << std::endl;
    // init_progress_bar(n_trials * mymodel.nr_classes);
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        // Don't do this for (almost) empty classes
        if (mlModel.pdf_class.rptrAll()[iclass] < 0.01)
        {
            mlModel.acc_rot.wptrAll()[iclass]   = 999.;
            mlModel.acc_trans.wptrAll()[iclass] = 999.;
            continue;
        }
        
        double acc_rot_class = 0.;
        double acc_trans_class = 0.;
        // Particles are already in random order, so just move from 0 to n_trials
        for (int iimage = 0; iimage < n_trials; iimage++)
        {
            int igroup = metadata.accessAll(iimage).GROUP_NO - 1;
            // Get CTF for this particle
            if (do_ctf_correction)
            {
                CTF ctf;
                ctf.setValues(metadata.accessAll(iimage).CTF_DEFOCUS_U,
                              metadata.accessAll(iimage).CTF_DEFOCUS_V,
                              metadata.accessAll(iimage).CTF_DEFOCUS_ANGLE,
                              metadata.accessAll(iimage).CTF_VOLTAGE,
                              metadata.accessAll(iimage).CTF_CS,
                              metadata.accessAll(iimage).CTF_Q0,
                              metadata.accessAll(iimage).CTF_BFAC);
                
                ctf.getFftwImage(Fctf_wptr, exp_current_image_size, ori_size, ori_size, pixel_size,
                                 particleModel.ctf_phase_flipped, particleModel.only_flip_phases, particleModel.intact_ctf_first_peak, true);
            }
            
            // Search 2 times: ang and off
            // Don't estimate rotational accuracies if we're doing do_skip_rotate (for faster movie-frame alignment)
            bool do_skip_rotate = false;
            int imode_start = (do_skip_rotate) ? 1 : 0;
            for (int imode = imode_start; imode < 2; imode++)
            {
                double ang_error = 0.;
                double sh_error = 0.;
                double ang_step;
                double sh_step;
                double my_snr = 0.;
                
                // Search for ang_error and sh_error where there are at least 3-sigma differences!
                // 13feb12: change for explicit probability at P=0.01
                while (my_snr <= pvalue)
                {
                    // Graduallly increase the step size
                    if (ang_error < 0.2)
                        ang_step = 0.05;
                    else if (ang_error < 1.)
                        ang_step = 0.1;
                    else if (ang_error < 2.)
                        ang_step = 0.2;
                    else if (ang_error < 5.)
                        ang_step = 0.5;
                    else if (ang_error < 10.)
                        ang_step = 1.0;
                    else if (ang_error < 20.)
                        ang_step = 2;
                    else
                        ang_step = 5.0;
                    
                    if (sh_error < 0.2)
                        sh_step = 0.05;
                    else if (sh_error < 1.)
                        sh_step = 0.1;
                    else if (sh_error < 2.)
                        sh_step = 0.2;
                    else if (sh_error < 5.)
                        sh_step = 0.5;
                    else if (sh_error < 10.)
                        sh_step = 1.0;
                    else
                        sh_step = 2.0;
                    
                    ang_error += ang_step;
                    sh_error += sh_step;
                    
                    // Prevent an endless while by putting boundaries on ang_error and sh_error
                    if ( (imode == 0 && ang_error > 30.) || (imode == 1 && sh_error > 10.) )
                        break;
                    
                    dontShare_Random_generator.init(random_seed + iimage);
                    
                    double rot1 = metadata.accessAll(iimage).ROT;
                    double tilt1 = metadata.accessAll(iimage).TILT;
                    double psi1 = metadata.accessAll(iimage).PSI;
                    double xoff1 = 0.;
                    double yoff1 = 0.;
                    double zoff1 = 0.;
                    
                    // Get the FT of the first image
                    Euler_angles2matrix(rot1, tilt1, psi1, A1);
                    mapModel.get2DFourierTransform(iclass, F1Real_wptr, F1Imag_wptr, exp_current_image_size, A1, IS_NOT_INV);
                    
                    // Apply the angular or shift error
                    double rot2 = rot1;
                    double tilt2 = tilt1;
                    double psi2 = psi1;
                    double xshift = xoff1;
                    double yshift = yoff1;
                    double zshift = zoff1;
                    
                    // Perturb psi or xoff , depending on the mode
                    if (imode == 0)
                    {
                        assert(mapModel.ref_dim==3);
                        if (mapModel.ref_dim == 3)
                        {
                            // Randomly change rot, tilt or psi
                            double ran = dontShare_Random_generator.rnd_unif();
                            if (ran < 0.3333)
                                rot2 = rot1 + ang_error;
                            else if (ran < 0.6667)
                                tilt2 = tilt1 + ang_error;
                            else
                                psi2  = psi1 + ang_error;
                        }
                        else
                        {
                            psi2  = psi1 + ang_error;
                        }
                    }
                    else
                    {
                        // Randomly change xoff or yoff
                        double ran = dontShare_Random_generator.rnd_unif();
                        if (ran < 0.5)
                            xshift = xoff1 + sh_error;
                        else
                            yshift = yoff1 + sh_error;
                    }
                    
                    if (imode == 0)
                    {
                        // Get new rotated version of reference
                        Euler_angles2matrix(rot2, tilt2, psi2, A2);
                        mapModel.get2DFourierTransform(iclass, F2Real_wptr, F2Imag_wptr, exp_current_image_size, A2, IS_NOT_INV);
                    }
                    else
                    {
                        // Get shifted version
                        shiftImageInFourierTransformNoTab(F1Real_wptr, F1Imag_wptr, F2Real_wptr, F2Imag_wptr, exp_current_image_size, -xshift, -yshift, ori_size);
                    }
                    
                    // Apply CTF to F1 and F2 if necessary
                    if (do_ctf_correction)
                    {
                        for (int n = 0; n < exp_current_image_Fsize2; n++) {
                            F1Real_wptr[n] *= Fctf_wptr[n];
                            F1Imag_wptr[n] *= Fctf_wptr[n];
                            F2Real_wptr[n] *= Fctf_wptr[n];
                            F2Imag_wptr[n] *= Fctf_wptr[n];
                        }
                    }
                    
                    auto myMresol = (exp_current_image_size == coarse_size) ? Mresol_coarse.rptrAll() : Mresol_fine.rptrAll();
                    auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].rptrAll();
                    my_snr = 0.;
                    for (int n = 0; n < exp_current_image_Fsize2; n++)
                    {
                        int ires = myMresol[n];
                        if (ires > 0)
                        {
                            auto diff2 = (F1Real_wptr[n]-F2Real_wptr[n])*(F1Real_wptr[n]-F2Real_wptr[n])+(F1Imag_wptr[n]-F2Imag_wptr[n])*(F1Imag_wptr[n]-F2Imag_wptr[n]);
                            my_snr += diff2 / (2 * sigma2_fudge * sigma2_noise_igroup[ires]);
                        }
                    }

// #define DEBUG_ANGACC
#ifdef DEBUG_ANGACC
                    if (imode==0)
                    {
                        std::cerr << " ang_error= " << ang_error << std::endl;
                        std::cerr << " rot1= " << rot1 << " tilt1= " << tilt1 << " psi1= " << psi1 << std::endl;
                        std::cerr << " rot2= " << rot2 << " tilt2= " << tilt2 << " psi2= " << psi2 << std::endl;
                        
                    }
                    else
                    {
                        std::cerr << " xshift= " << xshift << " yshift= " << yshift << " zshift= " << zshift << std::endl;
                        std::cerr << " sh_error= " << sh_error << std::endl;
                    }
                    std::cerr << " my_snr= " << my_snr << std::endl;
                    FourierTransformer transformer;
                    MultidimArray<DOUBLE> spec_img(mymodel.ori_size), spec_diff(mymodel.ori_size), count(mymodel.ori_size);
                    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F1)
                    {
                        long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
                        spec_img(idx) += norm(dAkij(F1, k, i, j));
                        spec_diff(idx) += norm(dAkij(F1, k, i, j) - dAkij(F2, k, i, j));
                        count(idx) += 1.;
                    }
                    spec_img /= count;
                    spec_diff /= count;
                    for (int i=0; i < XSIZE(F1); i++)
                        std::cerr << " i= " << i << " spec_img(i)= "<< spec_img(i) << " spec_diff(i)= "<< spec_diff(i) << " sigma2_noise(i)= "<<  mymodel.sigma2_noise[group_id](i)
                        << " count(i)= " << count(i) << " sum-diff-norm=" << count(i)*spec_diff(i)<< " sum-diff-norm/sigma2= " << count(i)*(spec_diff(i)/mymodel.sigma2_noise[group_id](i))<< std::endl;
                    Image<DOUBLE> tt;
                    if (mymodel.data_dim == 3)
                        tt().resize(YSIZE(F1), YSIZE(F1), YSIZE(F1));
                    else
                        tt().resize(YSIZE(F1), YSIZE(F1));
                    transformer.inverseFourierTransform(F1, tt());
                    CenterFFT(tt(),false);
                    tt.write("F1.spi");
                    transformer.inverseFourierTransform(F2, tt());
                    CenterFFT(tt(),false);
                    tt.write("F2.spi");
                    std::cerr << "Written F1.spi and F2.spi. Press any key to continue... "<<std::endl;
                    char c;
                    std::cin >> c;
#endif
                    
                    // Only for the psi-angle and the translations, and only when my_prob < 0.01 calculate a histogram of the contributions at each resolution shell
                    if (my_snr > pvalue && imode == 0)
                    {
                        for (int n = 0; n < exp_current_image_Fsize2; n++) {
                            int ires = myMresol[n];
                            if (ires > 0) {
                                auto diff2 = (F1Real_wptr[n]-F2Real_wptr[n])*(F1Real_wptr[n]-F2Real_wptr[n])+(F1Imag_wptr[n]-F2Imag_wptr[n])*(F1Imag_wptr[n]-F2Imag_wptr[n]);
                                mlModel.orientability_contrib[iclass].wptrAll()[ires] += diff2 / (2 * sigma2_fudge * sigma2_noise_igroup[ires]);
                            }
                        }
                    }
                    
                } // end while my_snr >= pvalue
                if (imode == 0)
                    acc_rot_class += ang_error;
                else if (imode == 1)
                    acc_trans_class += sh_error;
            } // end for imode
        }// end for iimage
        
        // progress_bar(n_trials*iclass + my_metadata_entry);
        
        mlModel.acc_rot.wptrAll()[iclass]   = acc_rot_class / (double)n_trials;
        mlModel.acc_trans.wptrAll()[iclass] = acc_trans_class / (double)n_trials;
        
        // Store normalised spectral contributions to orientability
        double orientability_contrib_sum = sumVec(mlModel.orientability_contrib[iclass].wptrAll(), mlModel.ori_Fsize);
        if (orientability_contrib_sum > 0.)
            mlModel.orientability_contrib[iclass] /= orientability_contrib_sum;
        
        // Keep the orientational accuracy of the best class for the auto-sampling approach
        acc_rot     = std::min(mlModel.acc_rot.rptrAll()[iclass], acc_rot);
        acc_trans   = std::min(mlModel.acc_trans.rptrAll()[iclass], acc_trans);
        
        
        // Richard's formula with Greg's constant
        //DOUBLE b_orient = (acc_rot_class*acc_rot_class* particle_diameter*particle_diameter) / 3000.;
        //std::cout << " + expected B-factor from the orientational errors = "
        //		<< b_orient<<std::endl;
        // B=8 PI^2 U^2
        //std::cout << " + expected B-factor from the translational errors = "
        //		<< 8 * PI * PI * mymodel.pixel_size * mymodel.pixel_size * acc_trans_class * acc_trans_class << std::endl;
        
    } // end loop iclass
    // progress_bar(n_trials * mymodel.nr_classes);
    
    std::cout << " Auto-refine: calculateExpectedAngularErrors() "<<std::endl;
    std::cout << " Auto-refine: Estimated accuracy angles (acc_rot) = " << acc_rot<< " degrees; offsets= " << acc_trans << " pixels" << std::endl;
    // Warn for inflated resolution estimates
    if (acc_rot > 10. && do_auto_refine)
    {
        std::cout << " Auto-refine: WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles (yet)!" << std::endl;
        std::cout << " Auto-refine: WARNING: You probably need not worry if the accuracy improves during the next few iterations." << std::endl;
        std::cout << " Auto-refine: WARNING: However, if the problem persists it may lead to spurious FSC curves, so be wary of inflated resolution estimates..." << std::endl;
        std::cout << " Auto-refine: WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
    }
}

void joinTwoHalvesAtLowResolution()
{
    if (!do_split_random_halves)
        ERROR_REPORT("BUG: you should not be in joinTwoHalvesAtLowResolution!");
    
    // Loop over all classes (this will be just one class for now...)
    FDOUBLE myres = std::max(low_resol_join_halves, 1./mapModel.current_resolution);
    int lowres_r_max = ceil(ori_size * mapModel.pixel_size / myres);
    
    for (int iclass = 0; iclass < nr_classes; iclass++ )
    {
        IF_HALF_ONE_TWO_MASTER_NODE
        {
#ifdef DATA_STREAM
            global_data_stream.foutInt(lowres_r_max, "joinTwoHalvesAtLowResolution()_lowres_r_max", __FILE__, __LINE__);
            global_data_stream.foutDoubleComplex((std::complex<double>*)mapModel.backprojector[iclass].data.wptr(),mapModel.backprojector[iclass].data.dimzyx,
                                                 "joinTwoHalvesAtLowResolution()_wsum_model.BPref[iclass].data_before", __FILE__, __LINE__);
            global_data_stream.foutDouble(mapModel.backprojector[iclass].weight.wptr(), mapModel.backprojector[iclass].weight.dimzyx,
                                          "joinTwoHalvesAtLowResolution()_wsum_model.BPref[iclass].weight_before", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            Vol<MKL_Complex > lowres_data;
            Vol<FDOUBLE > lowres_weight;
            mapModel.backprojector[iclass].getLowResDataAndWeight(lowres_data, lowres_weight, lowres_r_max);
            
#ifdef DATA_STREAM
            global_data_stream.foutDoubleComplex((std::complex<double>*)lowres_data.wptr(), lowres_data.dimzyx, "joinTwoHalvesAtLowResolution()_lowres_data_before", __FILE__, __LINE__);
            global_data_stream.foutDouble(lowres_weight.wptr(), lowres_weight.dimzyx, "joinTwoHalvesAtLowResolution()_lowres_weight_before", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            IF_HALF_TWO_MASTER_NODE
            {
                // The second slave sends its lowres_data and lowres_weight to the first slave
                MPI::COMM_WORLD.Send((FDOUBLE*)lowres_data.rptr(), 2*lowres_data.dimzyx, MPI_FDOUBLE, MPI_HALF_ONE_MASTER_NODE, 188);
                MPI::COMM_WORLD.Send(lowres_weight.rptr(), lowres_weight.dimzyx, MPI_FDOUBLE, MPI_HALF_ONE_MASTER_NODE, 189);
                
                // Now the first slave is calculating the average....
                
                // Then the second slave receives the average back from the first slave
                MPI::COMM_WORLD.Recv((FDOUBLE*)lowres_data.wptr(), 2*lowres_data.dimzyx, MPI_FDOUBLE, MPI_HALF_ONE_MASTER_NODE, 190);
                MPI::COMM_WORLD.Recv(lowres_weight.wptr(), lowres_weight.dimzyx, MPI_FDOUBLE, MPI_HALF_ONE_MASTER_NODE, 191);
            }
            IF_HALF_ONE_MASTER_NODE
            {
                std::cout << " Averaging half-reconstructions up to " << myres << " Angstrom resolution to prevent diverging orientations ..." << std::endl;
                std::cout << " Note that only for higher resolutions the FSC-values are according to the gold-standard!" << std::endl;
                
                Vol<MKL_Complex> lowres_data_half2;
                Vol<FDOUBLE > lowres_weight_half2;
                lowres_data_half2.init(lowres_data.dimz, lowres_data.dimy, lowres_data.dimx);
                lowres_weight_half2.init(lowres_weight.dimz, lowres_data.dimy, lowres_weight.dimx);

                // The first slave receives the average from the second slave
                MPI::COMM_WORLD.Recv((FDOUBLE*)lowres_data_half2.wptr(), 2*lowres_data_half2.dimzyx, MPI_FDOUBLE, MPI_HALF_TWO_MASTER_NODE, 188);
                MPI::COMM_WORLD.Recv((FDOUBLE*)lowres_weight_half2.wptr(), lowres_weight_half2.dimzyx, MPI_FDOUBLE, MPI_HALF_TWO_MASTER_NODE, 189);
                
                // The first slave calculates the average of the two lowres_data and lowres_weight arrays
                for (size_t n = 0; n < lowres_data.dimzyx; n++)
                {
                    ACCESS(lowres_data, 0, 0, n).real += ACCESS(lowres_data_half2, 0, 0, n).real;
                    ACCESS(lowres_data, 0, 0, n).imag += ACCESS(lowres_data_half2, 0, 0, n).imag;
                    ACCESS(lowres_data, 0, 0, n).real   /= 2.;
                    ACCESS(lowres_data, 0, 0, n).imag   /= 2.;
                    ACCESS(lowres_weight, 0, 0, n) += ACCESS(lowres_weight_half2, 0, 0, n) ;
                    ACCESS(lowres_weight, 0, 0, n) /= 2.;
                }
                
                // The first slave sends the average lowres_data and lowres_weight also back to the second slave
                MPI::COMM_WORLD.Send((FDOUBLE*)lowres_data.rptr(), 2*lowres_data.dimzyx, MPI_FDOUBLE, MPI_HALF_TWO_MASTER_NODE, 190);
                MPI::COMM_WORLD.Send(lowres_weight.rptr(), lowres_weight.dimzyx, MPI_FDOUBLE, MPI_HALF_TWO_MASTER_NODE, 191);
            }
#ifdef DATA_STREAM
            global_data_stream.foutDoubleComplex((std::complex<double>*)lowres_data.wptr(), lowres_data.dimzyx, "joinTwoHalvesAtLowResolution()_lowres_data_after", __FILE__, __LINE__);
            global_data_stream.foutDouble(lowres_weight.wptr(), lowres_weight.dimzyx, "joinTwoHalvesAtLowResolution()_lowres_weight_after", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
            // Now that both slaves have the average lowres arrays, set them back into the backprojector
            mapModel.backprojector[iclass].setLowResDataAndWeight(lowres_data, lowres_weight, lowres_r_max);
            //
#ifdef DATA_STREAM
            global_data_stream.foutDoubleComplex((std::complex<double>*)mapModel.backprojector[iclass].data.wptr(),mapModel.backprojector[iclass].data.dimzyx,
                                                 "joinTwoHalvesAtLowResolution()_wsum_model.BPref[iclass].data_after", __FILE__, __LINE__);
            global_data_stream.foutDouble(mapModel.backprojector[iclass].weight.wptr(), mapModel.backprojector[iclass].weight.dimzyx,
                                          "joinTwoHalvesAtLowResolution()_wsum_model.BPref[iclass].weight_after", __FILE__, __LINE__);
            global_data_stream.check();global_data_stream.flush();
#endif
        }
    }
}
    
void compareTwoHalves()
{
    //
    if (!do_split_random_halves)
        ERROR_REPORT("ERROR: you should not be in MlOptimiserMpi::compareTwoHalves!");
    
    std::vector<FDOUBLE> fsc(ori_size/2+1,0);
    // Loop over all classes
    for (int iclass = 0; iclass < nr_classes; iclass++ )
    {
        IF_HALF_ONE_TWO_MASTER_NODE
        {
            // The first two slaves calculate the downsampled average
            Vol<MKL_Complex > avg1;
            mapModel.backprojector[iclass].getDownsampledAverage(avg1);
            IF_HALF_TWO_MASTER_NODE
            {
                // The second slave sends its average to the first slave
                MPI::COMM_WORLD.Send((FDOUBLE*)avg1.rptr(), 2*avg1.dimzyx, MPI_FDOUBLE, MPI_HALF_ONE_MASTER_NODE, 888);
            }
            IF_HALF_ONE_MASTER_NODE
            {
                std::cout << " Calculating gold-standard FSC ..."<< std::endl;
                // The first slave receives the average from the second slave and calculates the FSC between them
                Vol<MKL_Complex > avg2;
                avg2.init(avg1.dimz, avg1.dimy, avg1.dimx);
                MPI::COMM_WORLD.Recv((FDOUBLE*)avg2.wptr(), 2*avg2.dimzyx, MPI_FDOUBLE, MPI_HALF_TWO_MASTER_NODE, 888);
                mapModel.backprojector[iclass].calculateDownSampledFourierShellCorrelation(avg1, avg2, fsc.data());
            }
        }
        // Now slave 1 sends the fsc curve to everyone else
        MPI::COMM_WORLD.Bcast(fsc.data(), ori_size/2+1, MPI_FDOUBLE, 0);
        //
        IF_HALF_ONE_TWO_NODE
            copy(mlModel.fsc_halves_class[iclass].wptrAll(), mlModel.fsc_halves_class[iclass].size(), fsc.data(), fsc.size());

#ifdef DATA_STREAM
        global_data_stream.foutDouble(fsc.data(), ori_size/2+1, "mymodel.fsc_halves_class[iclass]", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    }
}

void combineWeightedSumsTwoRandomHalves()
{
    int color;
    IF_HALF_ONE_TWO_MASTER_NODE color = 0;
    else color = 1;
    // put random two half master node to master_world
    MPI::Intracomm master_world;
    master_world = MPI::COMM_WORLD.Split(color, node);
    // reduce wsum for two random half
    IF_HALF_ONE_TWO_MASTER_NODE mlModel.reduceData(master_world);
    //
    // reduce wsum Iref
    IF_HALF_ONE_TWO_MASTER_NODE mapModel.reduceData(master_world);
}
    
void checkConvergence()
{
//    if (iter >= 3)
//    {
//        assert(has_fine_enough_angular_sampling);
//        has_converged = true;
//        do_join_random_halves = true;
//        // In the last iteration, include all data until Nyquist
//        do_use_all_data = true;
//        return;
//    }
    has_converged = false;
    if ( has_fine_enough_angular_sampling && nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN && nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES )
    {
        has_converged = true;
        do_join_random_halves = true;
        // In the last iteration, include all data until Nyquist
        do_use_all_data = true;
    }
}

void writeTemporaryDataAndWeightArrays()
{
    IF_HALF_ONE_TWO_MASTER_NODE
    {
        int random_subset = node==0?1:2;
        std::string fn = write_path+write_fn+"_half"+num2str(random_subset,1);
        // Write out temporary arrays for all classes
        for (int iclass = 0; iclass < mlModel.nr_classes; iclass++)
        {
            if (mlModel.pdf_class.rptrAll()[iclass] > 0.)
            {
                std::string fn_tmp = fn+"_class"+num2str(iclass+1,3);
                mapModel.backprojector[iclass].data.writeToDisk(fn_tmp+"_data",false);
                mapModel.backprojector[iclass].weight.writeToDisk(fn_tmp+"_weight",false);
            }
        }
    }
}

void readTemporaryDataAndWeightArraysAndReconstruct()
{
    IF_HALF_ONE_TWO_MASTER_NODE
    {
        int random_subset = node==0?1:2;
        std::string fn = write_path+write_fn+"_half"+num2str(random_subset,1);
        //
        for (int iclass = 0; iclass < mlModel.nr_classes; iclass++)
        {
            std::vector<FDOUBLE> dummy(ori_size/2+1,0);
            Vol<FDOUBLE> Iunreg;
            Iunreg.init(ori_size, ori_size, ori_size, true);
            if (mlModel.pdf_class.rptrAll()[iclass] > 0.)
            {
                std::string fn_tmp = fn+"_class"+num2str(iclass+1,3);
                mapModel.backprojector[iclass].data.readFromDisk(fn_tmp+"_data",false);
                mapModel.backprojector[iclass].weight.readFromDisk(fn_tmp+"_weight",false);
                
                //
                // Now perform the unregularized reconstruction
                mapModel.backprojector[iclass].reconstruct(Iunreg, gridding_nr_iter, false, 1, dummy.data(), dummy.data(), dummy.data(),
                                                           dummy.data(), 1, false, true, -1, omp_get_max_threads());

//                // Update header information
//                DOUBLE avg, stddev, minval, maxval;
//                Iunreg().computeStats(avg, stddev, minval, maxval);
//                Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
//                Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
//                Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
//                Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
//                Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, mymodel.pixel_size);
//                Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, mymodel.pixel_size);
//                Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, mymodel.pixel_size);
                
                // And write the resulting model to disc
                std::string fn_mrc = fn_tmp + "_unfil.mrc";
                Mrcs::MrcVolume oneImage(Iunreg.wptr(),ori_size,pixel_size);
                oneImage.write(fn_mrc,mapModel.refHead);
                
                // remove temporary arrays from the disc
                remove((fn_tmp+"_data.vol").c_str());
                remove((fn_tmp+"_weight.vol").c_str());
            }
            Iunreg.fini();
        }
    }
    
}
    
void printConvergenceStats()
{
    std::cout << " Auto-refine: --------------------------------------------------------------------------------- "<<std::endl;
    std::cout << " Auto-refine: printConvergenceStats() "<<std::endl;
    std::cout << " Auto-refine: Iteration= "<< iter<< std::endl;
    std::cout << " Auto-refine: Resolution= "<< 1./mapModel.current_resolution<< " (no gain for " << nr_iter_wo_resol_gain << " iter) "<< std::endl;
    std::cout << " Auto-refine: Changes in angles= " << current_changes_optimal_orientations
              << " degrees; and in offsets= " << current_changes_optimal_offsets
   			  << " pixels (no gain for " << nr_iter_wo_large_hidden_variable_changes << " iter) "<< std::endl;
    //std::cout << " offset_range  = "<<sampler3d.offset_range<<" ,offset_step = "<<sampler3d.offset_step<<" ,healpix order = "<<sampler3d.healpix_order<<std::endl;
    //std::cout << " do local searching = "<<mlModel.orientational_prior_mode<<" "<<sampler3d.orientational_prior_mode<<std::endl;
    std::cout << " has_fine_enough_angular_sampling = "<<has_fine_enough_angular_sampling<<std::endl;
    if (has_converged)
    {
        std::cout << " Auto-refine: Refinement has converged, entering last iteration where two halves will be combined..."<<std::endl;
    }
    std::cout << "------------------------------------------------------"<<std::endl;
    
#ifdef DATA_STREAM
    global_data_stream.foutInt(iter, "printConvergenceStats()_iter", __FILE__, __LINE__);
    global_data_stream.foutDouble(1./mapModel.current_resolution, "printConvergenceStats()_Resolution", __FILE__, __LINE__);
    global_data_stream.foutInt(nr_iter_wo_resol_gain, "printConvergenceStats()_nr_iter_wo_resol_gain", __FILE__, __LINE__);
    global_data_stream.foutDouble(current_changes_optimal_orientations, "printConvergenceStats()_current_changes_optimal_orientations", __FILE__, __LINE__);
    global_data_stream.foutDouble(current_changes_optimal_offsets, "printConvergenceStats()_current_changes_optimal_offsets", __FILE__, __LINE__);
    global_data_stream.foutInt(nr_iter_wo_large_hidden_variable_changes, "printConvergenceStats()_nr_iter_wo_large_hidden_variable_changes", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}
    
}// end namespace Map3dAutoRefinement
