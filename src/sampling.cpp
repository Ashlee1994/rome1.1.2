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
#include "sampling.h"

HealpixSampler::HealpixSampler()
{
    is_3D = false;
    healpix_order = -1;
    fn_sym = "C1";
    pgGroup = pgOrder = 0;
    random_perturbation = 0;
    perturbation_factor = 0.5;
}

HealpixSampler::~HealpixSampler()
{
    translations_x  .resize(0);
    translations_y  .resize(0);
    psi_angles      .resize(0);
    rot_angles      .resize(0);
    tilt_angles     .resize(0);
    healpix_order   = -1;
    single_translations			.resize(0);
    single_translations_index	.resize(0);
}


void HealpixSampler::initialize(FDOUBLE _offset_step /*= -1.*/,FDOUBLE _offset_range /*= -1.*/,
                                FDOUBLE _psi_step /*= -1*/,int _order /*= -1*/,std::string _fn_sym /*= "C1"*/)
{
    // --------- initialize offset,rotate,healpix step -----------
    if (_offset_step > 0. && _offset_range >= 0.)
    {
        offset_step = _offset_step;
        offset_range = _offset_range;
    }
    else{
        offset_step = 2;
        offset_range = 10;
    }
    
    if(_order > 0){
        healpix_order = _order;
        assert(healpix_order < 5);
        is_3D = true;
    }
    
    if(_psi_step == -1){
        if(is_3D)
            psi_step = 360. / (6 * round(std::pow(2., healpix_order)));
        else
            psi_step = 10;
    }
    else
        psi_step = _psi_step;
    
    // -------- initialize translations  ---------------
    nrTrans = 0;
    int maxr = ceil(offset_range / offset_step);
    // alloc enough space for translations
    translations_x  .resize((2*maxr+1)*(2*maxr+1));
    translations_y  .resize((2*maxr+1)*(2*maxr+1));
    int trans_ind_x = 0;
    for (int ix = -maxr; ix <= maxr; ix++,trans_ind_x++)
    {
        FDOUBLE xoff = ix * offset_step;
        single_translations.push_back(xoff);
        int trans_ind_y = 0;
        for (int iy = -maxr; iy <= maxr; iy++,trans_ind_y++)
        {
            FDOUBLE yoff = iy * offset_step;
            if (xoff*xoff + yoff*yoff <= offset_range * offset_range){
                translations_x[nrTrans] = xoff;
                translations_y[nrTrans] = yoff;
                nrTrans++;
                //
                single_translations_index.push_back(std::make_pair(trans_ind_x,trans_ind_y));
            }
        }
        
    }
    //
//#define DEBUG_SAMPLING
#ifdef DEBUG_SAMPLING
    std::cout<<std::setw(10)<<"itrans"<<std::setw(15)<<"translations_x"<<std::setw(15)<<"translations_y"<<std::endl;
    for(int itrans = 0;itrans < nrTrans;itrans++){
        std::cout<<std::setw(10)<<itrans<<std::setw(15)<<translations_x[itrans]<<std::setw(15)<<translations_y[itrans]<<std::endl;
        double trans_x = single_translations[single_translations_index[itrans].first];
        double trans_y = single_translations[single_translations_index[itrans].second];
        if(trans_x!=translations_x[itrans] || trans_y!=translations_y[itrans]){
            std::cout<<"!!!!"<<std::setw(10)<<itrans<<std::setw(15)<<trans_x<<std::setw(15)<<trans_y<<" "<<std::endl;
            ERROR_REPORT("diff..sampling");
        }
    }
#endif
    // ---------- initialize psi_angles   -----------------
    nrPsi = 0;
    int nr_psi = ceil(360./psi_step);
    psi_step = 360./(double)nr_psi;
    
    psi_angles  .resize(nr_psi);
    for (int ipsi = 0; ipsi < nr_psi; ipsi++)
    {
        psi_angles[nrPsi] = ipsi * psi_step;
        nrPsi++;
    }
    
    // ----------- initialize rot_angles,tilt_angles  --------------
    if (is_3D) {
#ifdef SAMPLING3D
        healpix_base.Set(healpix_order, NEST);
        fn_sym = _fn_sym;
        
        { // symmetry
            // Set up symmetry
            SymList SL;
            SL.isSymmetryGroup(fn_sym, pgGroup, pgOrder);
            SL.read_sym_file(fn_sym);
            
            // Precalculate (3x3) symmetry matrices
            Matrix2D<FDOUBLE>  L(4, 4), R(4, 4);
            Matrix2D<FDOUBLE>  Identity(3,3);
            Identity.initIdentity();
            R_repository.clear();
            L_repository.clear();
            R_repository.push_back(Identity);
            L_repository.push_back(Identity);
            for (int isym = 0; isym < SL.SymsNo(); isym++)
            {
                SL.get_matrices(isym, L, R);
                R.resize(3, 3);
                L.resize(3, 3);
                R_repository.push_back(R);
                L_repository.push_back(L);
            }
        }
        //
        nrPix = healpix_base.Npix();
        double zz, phi;
        
        rot_angles  	.resize(nrPix);
        tilt_angles 	.resize(nrPix);
        directions_ipix	.resize(nrPix);
        for (int ipix = 0; ipix < nrPix; ipix++)
        {
            healpix_base.pix2ang_z_phi(ipix, zz, phi);
            rot_angles[ipix] = rad2deg(phi);
            tilt_angles[ipix] = acosd(zz);
            checkDirection(rot_angles[ipix], tilt_angles[ipix]);
            directions_ipix[ipix] = ipix;
        }
        
        {// symmetry
//            #define DEBUG_SAMPLING
#ifdef  DEBUG_SAMPLING
            writeAllOrientationsToBild("./orients_all.bild", "1 0 0 ", 0.020);
#endif
            // Now remove symmetry-related pixels
            // TODO check size of healpix_base.max_pixrad
            removeSymmetryEquivalentPoints(0.5 * RAD2DEG(healpix_base.max_pixrad()));
            
#ifdef  DEBUG_SAMPLING
            writeAllOrientationsToBild("./orients_sym.bild", "0 1 0 ", 0.021);
#endif
            nrPix = directions_ipix.size();
        }
        
#endif
    }
    else{
        nrPix = 1;
        assert(_fn_sym=="C1");
        fn_sym = "C1"; // This may not be set yet if restarting a 2D run....
    }
    // initialize random number
    
    resetRandomlyPerturbedSampling();
}

void HealpixSampler::resetRandomlyPerturbedSampling()
{
    // Actual instance of random perturbation
    // Add to the random perturbation from the last iteration, so it keeps changing strongly...
    random_perturbation += dontShare_Random_generator.rnd_unif(0.5*perturbation_factor, perturbation_factor);
    random_perturbation = realWrap(random_perturbation, -perturbation_factor, perturbation_factor);
}

size_t HealpixSampler::NrDir(int oversampling_order /*= 0*/) const
{
    if (oversampling_order == 0)
        return nrPix;
    else
        return round(std::pow(2., oversampling_order * 2))*nrPix;
}

size_t HealpixSampler::NrPsi(int oversampling_order /*= 0*/) const
{
    if (oversampling_order == 0)
        return nrPsi;
    else
        return round(pow(2., oversampling_order)) * nrPsi;
}

size_t HealpixSampler::NrTrans(int oversampling_order /*= 0*/) const
{
    if (oversampling_order == 0)
        return nrTrans;// translations_x.size();
    else
        return round(pow(2., oversampling_order * 2)) * nrTrans;//translations_x.size();
    
}

size_t HealpixSampler::NrPoints(int oversampling_order) const
{
    if (is_3D)
        return NrPsi(oversampling_order) * NrTrans(oversampling_order) * NrDir(oversampling_order);
    else
        return NrPsi(oversampling_order) * NrTrans(oversampling_order);
}

double HealpixSampler::getAngularSampling(int adaptive_oversampling /*= 0*/) const
{
    if (is_3D) {
        int order =  healpix_order + adaptive_oversampling;
        return 360. / (6 * round(std::pow(2., order)));
    }
    else
        return psi_step / pow(2., adaptive_oversampling);
}

size_t HealpixSampler::oversamplingFactorOrientations(int oversampling_order) const
{
	assert(oversampling_order < 4);
    if (is_3D)
        return round(std::pow(2., oversampling_order * 3));
    else
        return round(std::pow(2., oversampling_order));
}

int HealpixSampler::oversamplingFactorTranslations(int oversampling_order) const
{
	assert(oversampling_order < 4);
	return round(pow(2., oversampling_order * 2));
}

void HealpixSampler::getTranslations(int itrans,int oversampling_order,FDOUBLE* over_trans_x,FDOUBLE* over_trans_y) const
{
    assert(oversampling_order < 4);
    
    if (oversampling_order == 0)
    {
        over_trans_x[0] = translations_x[itrans];
        over_trans_y[0] = translations_y[itrans];
    }
    else
    {
        int nr_oversamples = round(std::pow(2., oversampling_order));
        for (int iover_trans_y = 0; iover_trans_y < nr_oversamples; iover_trans_y++)
        {
            for (int iover_trans_x = 0; iover_trans_x < nr_oversamples; iover_trans_x++)
            {
                FDOUBLE over_yoff = translations_y[itrans] - 0.5 * offset_step + (0.5 + iover_trans_y) * offset_step / nr_oversamples;
                FDOUBLE over_xoff = translations_x[itrans] - 0.5 * offset_step + (0.5 + iover_trans_x) * offset_step / nr_oversamples;
                
                int iover_trans = iover_trans_y*nr_oversamples+iover_trans_x;
                
                over_trans_x[iover_trans] = over_xoff;
                over_trans_y[iover_trans] = over_yoff;
            }
        }
    }
    
    if (fabs(random_perturbation) > 0.)
    {
        double myperturb = random_perturbation * offset_step;
        int nr_over_trans = oversamplingFactorTranslations(oversampling_order);
        for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++)
        {
            over_trans_x[iover_trans] += myperturb;
            over_trans_y[iover_trans] += myperturb;
        }
    }
}

void HealpixSampler::getAllTranslationsAndOverTrans(int oversampling_order,Vector1d& trans_x,Vector1d& trans_y,
                                                    Vector1d& trans_x_over,Vector1d& trans_y_over) const
{
    assert(oversampling_order < 4);
    
    for (int itrans = 0; itrans < nrTrans; itrans++)
    {
        trans_x[itrans] = translations_x[itrans];
        trans_y[itrans] = translations_y[itrans];
    }
    //
    if (oversampling_order > 0)
    {
        int nr_oversamples = round(std::pow(2., oversampling_order));
        for (int iover_trans_y = 0; iover_trans_y < nr_oversamples; iover_trans_y++)
        {
            for (int iover_trans_x = 0; iover_trans_x < nr_oversamples; iover_trans_x++)
            {
                double over_yoff = - 0.5 * offset_step + (0.5 + iover_trans_y) * offset_step / nr_oversamples;
                double over_xoff = - 0.5 * offset_step + (0.5 + iover_trans_x) * offset_step / nr_oversamples;
                
                int iover_trans = iover_trans_y*nr_oversamples+iover_trans_x;
                
                trans_x_over[iover_trans] = over_xoff;
                trans_y_over[iover_trans] = over_yoff;
            }
        }
    }
    else
    {
        trans_x_over[0] = 0;
        trans_y_over[0] = 0;
    }
    //
    if (fabs(random_perturbation) > 0.)
    {
        double myperturb = random_perturbation * offset_step;
#ifdef TTTT // where to add perturbation
        for (int itrans = 0; itrans < nrTrans; itrans++) {
            trans_x[itrans] += myperturb;
            trans_y[itrans] += myperturb;
        }
#else
        if (oversampling_order > 0){
            int nr_oversamples = round(std::pow(2., oversampling_order));
            int exp_nr_over_trans = nr_oversamples*nr_oversamples;
            for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
                trans_x_over[iover_trans] += myperturb;
                trans_y_over[iover_trans] += myperturb;
            }
        }
        else{
            trans_x_over[0] += myperturb;
            trans_y_over[0] += myperturb;
        }
#endif
    }
}

//void HealpixSampler::getAllSingleTranslationsAndOverTrans(int oversampling_order,Vector1d& trans_x,Vector1d& trans_y,
//                                                    	  Vector1d& trans_x_over,Vector1d& trans_y_over) const
//{
//    
//}

void HealpixSampler::getOrientations(int idir,int ipsi,int oversampling_order,FDOUBLE* over_psi,
                                     FDOUBLE* over_rot /*= nullptr*/,FDOUBLE* over_tilt /*= nullptr*/) const
{
    assert(oversampling_order < 4);
    if (!is_3D) { // 2D case
        if (oversampling_order == 0)
        {
            over_psi[0] = psi_angles[ipsi];
        }
        else
        {
            // for 2D sampling, only push back oversampled psi rotations
            int nr_ipsi_over = round(pow(2., oversampling_order));
            for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
            {
                double overpsi = psi_angles[ipsi] - 0.5 * psi_step + (0.5 + ipsi_over) * psi_step / nr_ipsi_over;
                over_psi[ipsi_over] = overpsi;
            }
        }
    }
    else{ // 3D case
#ifdef SAMPLING3D
        if (oversampling_order == 0)
        {
            over_psi[0] = psi_angles[ipsi];
            over_tilt[0] = tilt_angles[idir];
            over_rot[0] = rot_angles[idir];
        }
        else{
            Healpix_Base HealPixOver(oversampling_order + healpix_order, NEST);
            int fact = HealPixOver.Nside()/healpix_base.Nside();
            int nr_ipsi_over = round(std::pow(2., oversampling_order));
            int nr_idir_over = round(std::pow(2., oversampling_order*2));
            int x, y, face;
            FDOUBLE rot, tilt;
            int ipix = directions_ipix[idir];
            int idir_over = 0;
            // Get x, y and face for the original, coarse grid
            healpix_base.nest2xyf(ipix, x, y, face);
            // Loop over the oversampled Healpix pixels on the fine grid
            for (int j = fact * y; j < fact * (y+1); ++j)
            {
                for (int i = fact * x; i < fact * (x+1); ++i)
                {
                    int overpix = HealPixOver.xyf2nest(i, j, face);
                    // this one always has to be double (also for SINGLE_PRECISION CALCULATIONS) for call to external library
                    double zz, phi;
                    HealPixOver.pix2ang_z_phi(overpix, zz, phi);
                    rot = rad2deg(phi);
                    tilt = acosd(zz);
                    
                    // The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
                    checkDirection(rot, tilt);
                    
                    for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
                    {
                        double overpsi = psi_angles[ipsi] - 0.5 * psi_step + (0.5 + ipsi_over) * psi_step / nr_ipsi_over;
                        int iover = idir_over*nr_ipsi_over+ipsi_over;
                        assert(iover >= 0);
                        over_rot[iover] = rot;
                        over_tilt[iover] = tilt;
                        over_psi[iover] = overpsi;
                    }// end loop ipsi_over
                    idir_over++;
                }// end loop i
            } // end loop j
            assert(idir_over == nr_idir_over);
        }// oversampling_order != 0
#endif
    }// 3D case
    
    
    // Random perturbation
    if (fabs(random_perturbation) > 0.)
    {
        double myperturb = random_perturbation * getAngularSampling();
        
        if (!is_3D) { // 2D case
            int nr_ipsi_over = round(pow(2., oversampling_order));
            for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
                over_psi[ipsi_over] += myperturb;
        }
        else{ // 3D case
            FDOUBLE A[3][3],R[3][3],AR[3][3];
            int nr_over = oversamplingFactorOrientations(oversampling_order);
            for (int iover = 0; iover < nr_over; iover++) {
                Euler_angles2matrix(over_rot[iover],over_tilt[iover],over_psi[iover],A);
                Euler_angles2matrix(myperturb, myperturb, myperturb, R);
#define A_MULTIPLY_R(i,j) AR[i][j] = A[i][0]*R[0][j]+A[i][1]*R[1][j]+A[i][2]*R[2][j];
                A_MULTIPLY_R(0,0) A_MULTIPLY_R(0,1) A_MULTIPLY_R(0,2)
                A_MULTIPLY_R(1,0) A_MULTIPLY_R(1,1) A_MULTIPLY_R(1,2)
                A_MULTIPLY_R(2,0) A_MULTIPLY_R(2,1) A_MULTIPLY_R(2,2)
                Euler_matrix2angles(AR,over_rot[iover],over_tilt[iover],over_psi[iover]);
            }// end loop iover
        }// 3D case
    }
}

void HealpixSampler::getAllTranslations(int oversampling_order, FDOUBLE* my_translations_x, FDOUBLE* my_translations_y) const
{
    assert(oversampling_order < 4);
    //
    int nr_over_trans = oversamplingFactorTranslations(oversampling_order);
    for (int itrans = 0; itrans < nrTrans; itrans++)
        getTranslations(itrans,oversampling_order,
                        my_translations_x+itrans*nr_over_trans,
                        my_translations_y+itrans*nr_over_trans);
}

void HealpixSampler::getAllTranslations(int oversampling_order,Vector2d& my_translations_x,Vector2d& my_translations_y) const
{
    assert(oversampling_order < 4);
    //
    int nr_over_trans = oversamplingFactorTranslations(oversampling_order);
    for (int itrans = 0; itrans < nrTrans; itrans++)
        getTranslations(itrans,oversampling_order,
                        my_translations_x[itrans].data(),
                        my_translations_y[itrans].data());
}

// 2D
void HealpixSampler::getAllOrientations(int oversampling_order,FDOUBLE* my_psi) const
{
    assert(oversampling_order < 4);
    // TODO
    int nr_over = oversamplingFactorOrientations(oversampling_order);
    int idir = 0;
    for (int ipsi = 0; ipsi < nrPsi; ipsi++) {
        int offset = ipsi*nr_over;
        getOrientations(idir,ipsi,oversampling_order,my_psi+offset);
    }
}

// 3D
void HealpixSampler::getAllOrientations(int oversampling_order,Vector2d& my_psi,Vector2d& my_rot,Vector2d& my_tilt) const
{
    assert(oversampling_order < 4);
    int nr_over = oversamplingFactorOrientations(oversampling_order);
	#pragma omp parallel for collapse(2)
    for (int idir = 0; idir < nrPix; idir++) {
        for (int ipsi = 0; ipsi < nrPsi; ipsi++) {
            getOrientations(idir,ipsi,oversampling_order,
                            my_psi[idir*nrPsi+ipsi].data(),
                            my_rot[idir*nrPsi+ipsi].data(),
                            my_tilt[idir*nrPsi+ipsi].data());
        }
    }
}

//
void HealpixSampler::writeOutSampling(std::string fn_sampling)
{
    std::ofstream  samplingFile;
    samplingFile.open((fn_sampling+".star").c_str(), std::ios::out);
    {
        samplingFile << std::endl;
        samplingFile << "data_sampling_general" <<std::endl;
        samplingFile << std::endl;
#define COUTMETADATA(v1,v2) samplingFile << "_rln" << std::left<<std::setw(30) << v1 << std::right<<std::setw(18) << v2 <<std::endl;
        //
        COUTMETADATA("Is3DSampling"				, is_3D					)
        COUTMETADATA("Is3DTranslationalSampling", 0						)
        COUTMETADATA("HealpixOrder"				, healpix_order			)
        COUTMETADATA("SymmetryGroup"			, fn_sym				)
        COUTMETADATA("TiltAngleLimit"			, -91					)
        COUTMETADATA("PsiStep"					, psi_step				)
        COUTMETADATA("OffsetRange"				, offset_range			)
        COUTMETADATA("OffsetStep"				, offset_step			)
        COUTMETADATA("SamplingPerturbInstance"	, random_perturbation	)
        COUTMETADATA("SamplingPerturbFactor"	, 0.5					)
        //
#undef COUTMETADATA
        samplingFile << std::endl << std::endl;
    }
    // data_sampling_directions
    {
        SmallMetataDataTable data_sampling_directions("data_sampling_directions");
        data_sampling_directions.appendName({"AngleRot","AngleTilt"});
        data_sampling_directions.appendType({ElemTypeDouble,ElemTypeDouble});
        data_sampling_directions.appendElem({rot_angles.data(),tilt_angles.data()}, nrPix);
        data_sampling_directions.print(samplingFile);
    }
    samplingFile.close();
}

void HealpixSampler::readFromSampling(std::string fn_sampling)
{
    ifstreamCheckingExistence samplingFile(fn_sampling.c_str());
    {
        bool startingRead = false;
        std::string line;
        while (true) {
            if (startingRead) {
                double doubleTmp;std::string stringTmp;
#define CINMETADATADOUBLE(V) samplingFile >> line;samplingFile >> doubleTmp;MASTERNODE std::cout<<std::setw(30)<<line<<" "<<doubleTmp<<std::endl;
#define CINMETADATASTR(V) samplingFile >> line;samplingFile >> stringTmp;MASTERNODE std::cout<<std::setw(30)<<line<<" "<<stringTmp<<std::endl;
                //
                CINMETADATADOUBLE(	"Is3DSampling"				)
                CINMETADATADOUBLE(	"Is3DTranslationalSampling"	)
                CINMETADATADOUBLE(	"HealpixOrder"				)
                CINMETADATASTR(		"SymmetryGroup"				)
                CINMETADATADOUBLE(	"TiltAngleLimit"			)
                CINMETADATADOUBLE(	"PsiStep"					)
                CINMETADATADOUBLE(	"OffsetRange"				)
                CINMETADATADOUBLE(	"OffsetStep"				)
                CINMETADATADOUBLE(	"SamplingPerturbInstance"	);//random_perturbation = doubleTmp;
                CINMETADATADOUBLE(	"SamplingPerturbFactor"		);ERROR_CHECK(0.5!=doubleTmp, "Set sampling SamplingPerturbFactor.");
                //
#undef CINMETADATADOUBLE
#undef CINMETADATASTR
                break;
            }
            else{
                getline(samplingFile,line);
                MASTERNODE std::cout<<line<<std::endl;
                if ((line.find("data_sampling_general") !=std::string::npos) ){
                    startingRead = true;
                    getline(samplingFile,line);assert(line=="");// escape a empty line
                }
                ERROR_CHECK(samplingFile.eof(), "end of sampling file,can not find data_sampling_general.");
            }
        }
    }
    samplingFile.close();
}

// --------------------------------------------------------------------------------------------------- //

void SamplingGrid::initialize(HealpixSampler& sampler3d,int _adaptive_oversampling)
{
    exp_nr_trans = sampler3d.NrTrans();
    exp_nr_dir = sampler3d.NrDir();
    exp_nr_psi = sampler3d.NrPsi();
    size_t nr_orientation = exp_nr_dir*exp_nr_psi;
    adaptive_oversampling = _adaptive_oversampling;
    int exp_nr_over_trans_max = sampler3d.oversamplingFactorTranslations(adaptive_oversampling);
    int exp_nr_over_rot_max = sampler3d.oversamplingFactorOrientations(adaptive_oversampling);
    exp_over_psi	.resize(nr_orientation, std::vector<FDOUBLE>(exp_nr_over_rot_max,0));
    exp_over_rot	.resize(nr_orientation, std::vector<FDOUBLE>(exp_nr_over_rot_max,0));
    exp_over_tilt	.resize(nr_orientation, std::vector<FDOUBLE>(exp_nr_over_rot_max,0));
    exp_over_trans_x.resize(exp_nr_trans, std::vector<FDOUBLE>(exp_nr_over_trans_max,0));
    exp_over_trans_y.resize(exp_nr_trans, std::vector<FDOUBLE>(exp_nr_over_trans_max,0));
    exp_trans_x		.resize(exp_nr_trans);
    exp_trans_y		.resize(exp_nr_trans);
    exp_trans_x_over.resize(exp_nr_over_trans_max);
    exp_trans_y_over.resize(exp_nr_over_trans_max);
    exp_nr_over_trans = 0;
    exp_nr_over_rot = 0;
}

void SamplingGrid::finalize()
{
    exp_over_trans_x.resize(0);
    exp_over_trans_y.resize(0);
    exp_over_psi	.resize(0);
    exp_over_tilt	.resize(0);
    exp_over_rot	.resize(0);
    exp_trans_x		.resize(0);
    exp_trans_y		.resize(0);
    exp_trans_x_over.resize(0);
    exp_trans_y_over.resize(0);
    exp_nr_trans = 0;exp_nr_dir = 0;exp_nr_psi = 0;
    exp_nr_over_trans = 0;exp_nr_over_rot = 0;
    adaptive_oversampling = 0;current_oversampling = 0;
}

void SamplingGrid::computeGrid(HealpixSampler& sampler3d,int _current_oversampling)
{
    assert(current_oversampling <= adaptive_oversampling);
    current_oversampling = _current_oversampling;
    sampler3d.getAllOrientations(current_oversampling, exp_over_psi, exp_over_rot, exp_over_tilt);
    sampler3d.getAllTranslations(current_oversampling, exp_over_trans_x, exp_over_trans_y);
    sampler3d.getAllTranslationsAndOverTrans(current_oversampling, exp_trans_x, exp_trans_y, exp_trans_x_over, exp_trans_y_over);
    exp_nr_over_rot = sampler3d.oversamplingFactorOrientations(current_oversampling);
    exp_nr_over_trans = sampler3d.oversamplingFactorTranslations(current_oversampling);
}

void SamplingGrid::testGetShift()
{
    //std::cout<<"Starting compare get shift.."<<std::endl;
    for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
            FDOUBLE shiftx1,shifty1;
            FDOUBLE shiftx2,shifty2,shiftx2Over,shifty2Over;
            getShiftxy(shiftx1, shifty1, itrans, iover_trans);
            getShiftxy(shiftx2, shifty2, shiftx2Over, shifty2Over, itrans, iover_trans);
            //std::cout<<"itrans = "<<itrans<<",iover_trans = "<<iover_trans<<std::endl;
            if (fabs(shiftx1-shiftx2-shiftx2Over) > 1e-20 || fabs(shifty1-shifty2-shifty2Over) > 1e-20) {
                ERROR_REPORT("Wrong getShiftxy..");
                //std::cout<<shiftx1<<" != "<<shiftx2<<" + "<<shiftx2Over<<","<<shifty1<<" != "<<shifty2<<" + "<<shifty2Over<<std::endl;
            }
            else{
                //std::cout<<shiftx1<<" = "<<shiftx2<<" + "<<shiftx2Over<<","<<shifty1<<" = "<<shifty2<<" + "<<shifty2Over<<std::endl;
            }
        }
    }
}

///

namespace GTMSampling{
    // info is var_num*3 matrix each row store basic info(start,end,order(indexBits))
    // notice ,if we declare an pointer out of this function and then allocate the memmory for X in this function
    // the allocated memory will be free...
    double* uniformSampling(const double *infos,int var_num,int &K,int &L){
        
        L = var_num;
        
        K = 1;
        int bit_sum = 0;
        for(int var = 0;var < L;var++){
            bit_sum += infos[var*3+2];
            K *= (1 << (int)infos[var*3+2]);
        }
        // K = pow(2,bit_sum)
        // std::cout<<"K = "<<K<<",L = "<<L<<std::endl;
        // std::cout<<"need space:"<<(double)K*L*8./1024./1024.<<" MB"<<std::endl;
        
        ERROR_CHECK(bit_sum > 32, "bit sum out of range(32).");
        
        double *X = new double [K*L];
        int subIndex = 0;
        int tempIndex;
        
        for(int index = 0;index < K;index++){
            tempIndex = index;
            
            for(int var = 0;var < L;var++){
                subIndex = tempIndex & ((1 << (int)infos[var*3+2]) - 1);
                
                X[index*L+var] = infos[var*3+0] + (infos[var*3+1] - infos[var*3+0])/((1 << (int)infos[var*3+2]) - 1)*subIndex;
                tempIndex = tempIndex >> (int)infos[var*3+2];
            }
        }
        
        std::cout<<"end sampling."<<std::endl;
        return X;
    }
    
    
    double* sample1L(const double *infos,int &K,int &L){
        
        L = 1;
        K = (int)infos[2];
        
        double *X = new double [K*L];
        
        for(int i = 0;i < K;i++)
            X[i] = infos[0] + (infos[1] - infos[0])/(K-1)*i;
        
        return X;
    }
    
    
    //   example:infos[6] = {1,3,3,1,2,2}
    //   result:
    //   1 1
    //   2 1
    //   3 1
    //   1 2
    //   2 2
    //   3 2
    //
    double* sample2L(const double *infos,int &K,int &L,bool inv){
        
        L = 2;
        K = (int)infos[2]*(int)infos[5];
        
        double *X = new double [K*L];
        int i,j;
        
        if(inv){
            for(int k = 0;k < K;k++){
                i = k / (int)infos[5];
                j = k % (int)infos[5];
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*j;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i = k / (int)infos[2];
                j = k % (int)infos[2];
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*j;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i;
            }
        }
        
        return X;
    }
    
    
    //   example:infos[9] = {1,3,3,1,2,2,1,2,2}
    //   result:
    //   1               1               1
    //   1               1               2
    //   1               2               1
    //   1               2               2
    //   2               1               1
    //   2               1               2
    //   2               2               1
    //   2               2               2
    //   3               1               1
    //   3               1               2
    //   3               2               1
    //   3               2               2
    //
    double* sample3L(const double *infos,int &K,int &L,bool inv){
        
        L = 3;
        K = (int)infos[2]*(int)infos[5]*(int)infos[8];
        
        double *X = new double [K*L];
        int i1,i2,i3,i;
        if(inv){
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[8]);
                
                i = k % ((int)infos[5]*(int)infos[8]);
                i2 = i / (int)infos[8];
                i3 = i % (int)infos[8];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[2]);
                
                i = k % ((int)infos[5]*(int)infos[2]);
                i2 = i / (int)infos[2];
                i3 = i % (int)infos[2];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i3;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i1;
            }
        }
        
        return X;
    }
    
    double* sample4L(const double *infos,int &K,int &L,bool inv){
        
        L = 4;
        K = (int)infos[2]*(int)infos[5]*(int)infos[8]*(int)infos[11];
        
        double *X = new double [K*L];
        int i1,i2,i3,i4,i;
        if(inv){
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[8]*(int)infos[11]);
                
                i = k % ((int)infos[5]*(int)infos[8]*(int)infos[11]);
                i2 = i / ((int)infos[8]*(int)infos[11]);
                
                i = i % ((int)infos[8]*(int)infos[11]);
                i3 = i / (int)infos[11];
                i4 = i % (int)infos[11];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i4 = k / ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                
                i = k % ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                i3 = i / ((int)infos[5]*(int)infos[2]);
                
                i = i % ((int)infos[5]*(int)infos[2]);
                i2 = i / (int)infos[2];
                i1 = i % (int)infos[2];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
            }
        }
        
        return X;
    }
    
    double* sample5L(const double *infos,int &K,int &L,bool inv){
        
        L = 5;
        K = (int)infos[2]*(int)infos[5]*(int)infos[8]*(int)infos[11]*(int)infos[14];
        
        double *X = new double [K*L];
        int i1,i2,i3,i4,i5,i;
        if(inv){
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[8]*(int)infos[11]*(int)infos[14]);
                
                i = k % ((int)infos[5]*(int)infos[8]*(int)infos[11]*(int)infos[14]);
                i2 = i / ((int)infos[8]*(int)infos[11]*(int)infos[14]);
                
                i = i % ((int)infos[8]*(int)infos[11]*(int)infos[14]);
                i3 = i / ((int)infos[11]*(int)infos[14]);
                
                i = i % ((int)infos[11]*(int)infos[14]);
                i4 = i / (int)infos[14];
                i5 = i % (int)infos[14];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
                X[k*L+4] = infos[12] + (infos[13] - infos[12])/((int)infos[14]-1)*i5;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i5 = k / ((int)infos[11]*(int)infos[8]*(int)infos[5]*(int)infos[2]);
                
                i = k % ((int)infos[11]*(int)infos[8]*(int)infos[5]*(int)infos[2]);
                i4 = i / ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                
                i = i % ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                i3 = i / ((int)infos[5]*(int)infos[2]);
                
                i = i % ((int)infos[5]*(int)infos[2]);
                i2 = i / (int)infos[2];
                i1 = i % (int)infos[2];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
                X[k*L+4] = infos[12] + (infos[13] - infos[12])/((int)infos[14]-1)*i5;
            }
        }
        
        return X;
    }
    
    // #include "./Healpix_2.15a/healpix_base.h"
    // double *sampleHealpix(int order,int &K,int &L){
    
    // 	std::cout<<"staring s sampling."<<std::endl;
    
    //     const double Pi = 3.1415926535897;
    
    //     Healpix_Base healpix_base;
    //     int healpix_order = order;   //N_side = 2^order
    //     healpix_base.Set(healpix_order, RING);   //initialize,RING//NEST
    //     int Npix = healpix_base.Npix();
    
    //     K = Npix;
    //     L = 2;
    
    //     ///////////////////some unnecessary test codes/////////////////////////////////
    //     std::cout<<"Order = "<<healpix_order<<std::endl;
    //     std::cout<<"Nside = "<<healpix_base.Nside()<<std::endl;
    //     std::cout<<"Npix = "<<Npix<<std::endl;
    
    //     double *X = (double*)_mm_malloc(sizeof(double)*2*Npix,64);
    
    //     double zz,phi;
    //     for(int ipix = 0;ipix < Npix;ipix++){
    //       healpix_base.pix2ang_z_phi(ipix, zz, phi);
    //       X[2*ipix] = (acos(zz)/*-Pi/2*/);
    //       X[2*ipix+1] = phi;
    //     }
    
    //     std::cout<<"end sampling."<<std::endl;
    
    //     return X;
    // }
}

/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
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
// The way symmetry is handled was copied from Xmipp.
// The original disclaimer is copied below
/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
void HealpixSampler::writeAllOrientationsToBild(std::string fn_bild, std::string rgb, FDOUBLE size)
{
    std::ofstream out;
    out.open (fn_bild.c_str());
    if (!out)
        ERROR_REPORT( (std::string)"HealpixSampling::writeAllOrientationsToBild: Cannot write file: " + fn_bild);
    
    
    out << ".color 1 0 0 \n";
    out << ".arrow 0 0 0 1 0 0 0.01 \n";
    out << ".color 0 1 0 \n";
    out << ".arrow 0 0 0  0 1 0 0.01 \n";
    out << ".color 0 0 1 \n";
    out << ".arrow 0 0 0 0 0 1 0.01 \n";
    
    
    Matrix1D<FDOUBLE> v(3);
    out << ".color " << rgb << std::endl;
    
    for (unsigned long int ipix = 0; ipix < rot_angles.size(); ipix++)
    {
        Euler_angles2direction(rot_angles[ipix], tilt_angles[ipix], v.data().data());
        out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v) <<  " " << std::to_string((long long)size) << std::endl;
    }
    
    out.close();
}

void HealpixSampler::writeBildFileOrientationalDistribution(VectorOfFDOUBLE &pdf_direction,
                                                             std::string &fn_bild, FDOUBLE R, FDOUBLE offset, FDOUBLE Rmax_frac, FDOUBLE width_frac)
{
    if (!is_3D)
        return;

    if (pdf_direction.size() != rot_angles.size())
        ERROR_REPORT("HealpixSampling::writeBildFileOrientationalDistribution XSIZE(pdf_direction) != rot_angles.size()!");
    
    
    FDOUBLE pdfmax, pdfmin, pdfmean, pdfsigma;
    computeStats(pdf_direction.rptrAll(), pdf_direction.size(), pdfmean, pdfsigma, pdfmin, pdfmax);
    
    std::ofstream fh_bild;
    fh_bild.open(fn_bild.c_str(), std::ios::out);
    if (!fh_bild)
        ERROR_REPORT("HealpixSampling::writeBildFileOrientationalDistribution: cannot open " + fn_bild);
    
    // 2 * PI * R = 360 degrees, 2*radius should cover angular sampling at width_frac=1
    FDOUBLE width = width_frac * PI*R*(getAngularSampling()/360.);
    Matrix1D<FDOUBLE> v(3);
    
    for (int iang = 0; iang < rot_angles.size(); iang++)
    {
        FDOUBLE pdf = pdf_direction[iang];
        
        // Don't make a cylinder for pdf==0
        if (pdf > 0.)
        {
            // Colour from blue to red according to deviations from sigma_pdf
            FDOUBLE colscale = (pdf - pdfmean) / pdfsigma;
            colscale = std::min(colscale, (FDOUBLE)5.);
            colscale = std::max(colscale, (FDOUBLE)-1.);
            colscale /= 6.;
            colscale += 1./6.; // colscale ranges from 0 (-5 sigma) to 1 (+5 sigma)
            
            // The length of the cylinder will depend on the pdf_direction
            FDOUBLE Rp = R + Rmax_frac * R * pdf / pdfmax;
            
            Euler_angles2direction(rot_angles[iang], tilt_angles[iang], v.data().data());
            
            // Don't include cylinders with zero length, as chimera will complain about that....
            if (fabs((R - Rp) * XX(v)) > 0.01 ||
                fabs((R - Rp) * YY(v)) > 0.01 ||
                fabs((R - Rp) * ZZ(v)) > 0.01)
            {
                // The width of the cylinders will be determined by the sampling:
                fh_bild << ".color " << colscale << " 0 " << 1. - colscale << std::endl;
                fh_bild << ".cylinder "
                << R  * XX(v) + offset << " "
                << R  * YY(v) + offset << " "
                << R  * ZZ(v) + offset << " "
                << Rp * XX(v) + offset << " "
                << Rp * YY(v) + offset << " "
                << Rp * ZZ(v) + offset << " "
                << width
                <<"\n";
            }
        }
        
    }
    
    // Close and write file to disc
    fh_bild.close();
}

void HealpixSampler::removeSymmetryEquivalentPoints(FDOUBLE max_ang)
{
    // Maximum distance
    FDOUBLE cos_max_ang = cos(DEG2RAD(max_ang));
    FDOUBLE my_dotProduct;
    Matrix1D<FDOUBLE>  direction(3), direction1(3);
    std::vector< Matrix1D<FDOUBLE> > directions_vector;
    
    // Calculate all vectors and fill directions_vector
    for (long int i = 0; i < rot_angles.size(); i++)
    {
        Euler_angles2direction(rot_angles[i], tilt_angles[i], direction.data().data());
        directions_vector.push_back(direction);
    }
    
    // First call to conventional remove_redundant_points
    removeSymmetryEquivalentPointsGeometric(pgGroup, pgOrder, directions_vector);
    
#ifdef  DEBUG_SAMPLING
    writeAllOrientationsToBild("orients_sym0.bild", "0 1 0", 0.021);
#endif
    
    // Only correct the seams (i.e. the borders of the asymmetric units) for small numbers of directions
    // For large numbers, the sampling is very fine and the probability distributions are probably delta functions anyway
    // Large numbers take long times to calculate...
    // Only a small fraction of the points at the border of the AU is thrown away anyway...
    if (rot_angles.size() < 4000)
    {
        // Create no_redundant vectors
        std::vector <Matrix1D<FDOUBLE> > no_redundant_directions_vector;
        std::vector <FDOUBLE> no_redundant_rot_angles;
        std::vector <FDOUBLE> no_redundant_tilt_angles;
        std::vector <int> no_redundant_directions_ipix;
        
        // Then check all points versus each other
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            
            direction1=directions_vector[i];
            bool uniq = true;
            
            //for (long int k = 0; k < no_redundant_directions_vector.size(); k++)
            // i is probably closer to latest additions: loop backwards over k....
            for (long int k = no_redundant_directions_vector.size() -1; k >= 0; k--)
            {
                for (int j = 0; j < R_repository.size(); j++)
                {
                    auto dir_mul_R = no_redundant_directions_vector[k].transpose() * R_repository[j];
                    direction =  L_repository[j] * ( (dir_mul_R).transpose() );
                    //Calculate distance
                    my_dotProduct = dotProduct(direction,direction1);
                    if (my_dotProduct > cos_max_ang)
                    {
                        uniq = false;
                        break;
                    }
                }// for j
                if (!uniq) break;
            } // for k
            
            if (uniq)
            {
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        } // for i
        
        // Now overwrite the rot/tilt_angles and directions_vectors with their no_redundant counterparts
        rot_angles = no_redundant_rot_angles;
        tilt_angles = no_redundant_tilt_angles;
        directions_ipix = no_redundant_directions_ipix;
    }
}

void HealpixSampler::removeSymmetryEquivalentPointsGeometric(const int symmetry,int sym_order,
                                                             std::vector <Matrix1D<FDOUBLE> >  &directions_vector)
{
    Matrix2D<FDOUBLE>  L(4, 4), R(4, 4);
    Matrix2D<FDOUBLE>  aux(3, 3);
    Matrix1D<FDOUBLE>  row1(3), row2(3), row(3);
    
    std::vector <Matrix1D<FDOUBLE> > no_redundant_directions_vector;
    std::vector <FDOUBLE> no_redundant_rot_angles;
    std::vector <FDOUBLE> no_redundant_tilt_angles;
    std::vector <int> no_redundant_directions_ipix;
    
    FDOUBLE my_dotProduct;
    if (symmetry == pg_CN)
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= (-180. / sym_order) &&
                rot_angles[i] <= (180. / sym_order))
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry == pg_CI  ||
             symmetry == pg_CS )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (tilt_angles[i] <= 90)
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNV )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=    0. / sym_order &&
                rot_angles[i] <=  180. / sym_order)
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNH )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180. / sym_order &&
                rot_angles[i] <=  180. / sym_order &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_SN )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180.*2. / sym_order &&
                rot_angles[i] <=  180.*2. / sym_order &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DN )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180. / (sym_order) + 90. &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DNV )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90.  &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DNH )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90. &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=   90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_T )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)-0.942809, (FDOUBLE)0., (FDOUBLE)0.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_3_fold_axis_3(3);
        _3_fold_axis_2_by_3_fold_axis_3 = vectorR3((FDOUBLE)0.471405, (FDOUBLE)0.272165, (FDOUBLE)0.7698);
        _3_fold_axis_2_by_3_fold_axis_3.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_3_by_3_fold_axis_1(3);
        _3_fold_axis_3_by_3_fold_axis_1 = vectorR3((FDOUBLE)0.471404, (FDOUBLE)0.816497, (FDOUBLE)0.);
        _3_fold_axis_3_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=     90. &&
                rot_angles[i] <=   150. ||
                rot_angles[i] ==     0
                )
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_3_fold_axis_3) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_3_by_3_fold_axis_1) >= 0
                    )
                {
                    no_redundant_rot_angles.push_back(rot_angles[i]);
                    no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_TD )
    {//OK
        Matrix1D<FDOUBLE>  _2_fold_axis_1_by_3_fold_axis_2(3);
        _2_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)-0.942809, (FDOUBLE)0., (FDOUBLE)0.);
        _2_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_3_fold_axis_5(3);
        _3_fold_axis_2_by_3_fold_axis_5 = vectorR3((FDOUBLE)0.471405, (FDOUBLE)0.272165, (FDOUBLE)0.7698);
        _3_fold_axis_2_by_3_fold_axis_5.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_5_by_2_fold_axis_1(3);
        _3_fold_axis_5_by_2_fold_axis_1 = vectorR3((FDOUBLE)0., (FDOUBLE)0.471405, (FDOUBLE)-0.666667);
        _3_fold_axis_5_by_2_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            //           if ( rot_angles[i]>=     120. &&
            //                 rot_angles[i]<=   150. ||
            //                 rot_angles[i]==     0
            //              )
            if (
                dotProduct(directions_vector[i], _2_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_2_by_3_fold_axis_5) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_5_by_2_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_TH )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_2_fold_axis_1(3);
        _3_fold_axis_1_by_2_fold_axis_1 = vectorR3((FDOUBLE)-0.816496, (FDOUBLE)0., (FDOUBLE)0.);
        _3_fold_axis_1_by_2_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _2_fold_axis_1_by_2_fold_axis_2(3);
        _2_fold_axis_1_by_2_fold_axis_2 = vectorR3((FDOUBLE)0.707107, (FDOUBLE)0.408248, (FDOUBLE)-0.57735);
        _2_fold_axis_1_by_2_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _2_fold_axis_2_by_3_fold_axis_1(3);
        _2_fold_axis_2_by_3_fold_axis_1 = vectorR3((FDOUBLE)-0.408248, (FDOUBLE)-0.707107, (FDOUBLE)0.);
        _2_fold_axis_2_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            //           if ( rot_angles[i]>=     120. &&
            //                 rot_angles[i]<=   150. ||
            //                 rot_angles[i]==     0
            //              )
            if (
                dotProduct(directions_vector[i], _3_fold_axis_1_by_2_fold_axis_1) >= 0 &&
                dotProduct(directions_vector[i], _2_fold_axis_1_by_2_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _2_fold_axis_2_by_3_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_O )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)-1., (FDOUBLE)1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3((FDOUBLE)1., (FDOUBLE)1., (FDOUBLE)0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3((FDOUBLE)-1., (FDOUBLE)1., (FDOUBLE)0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if ((rot_angles[i] >=   45. &&
                 rot_angles[i] <=  135. &&
                 tilt_angles[i] <=  90.) ||
                rot_angles[i] ==  0.
                )
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                    )
                {
                    no_redundant_rot_angles.push_back(rot_angles[i]);
                    no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_OH )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)-1., (FDOUBLE)1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3((FDOUBLE)1., (FDOUBLE)1., (FDOUBLE)0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3((FDOUBLE)-1., (FDOUBLE)1., (FDOUBLE)0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90. &&
                rot_angles[i] <=  135. &&
                tilt_angles[i] <=  90.)
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                    )
                {
                    no_redundant_rot_angles.push_back(rot_angles[i]);
                    no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_I || symmetry  == pg_I2)
    {//OK
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3((FDOUBLE)-0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3((FDOUBLE)0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I1)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0, 90, 0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)-0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);// TODO : check this??double to float????
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I3)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)-0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,-31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)0., (FDOUBLE)1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)-0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
                dotProduct(directions_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) <= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5)
    {//OK
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
        exit(0);
    }
    else if (symmetry  == pg_IH || symmetry  == pg_I2H)
    {//OK
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3((FDOUBLE)-0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3((FDOUBLE)0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  vectorR3((FDOUBLE)1.,(FDOUBLE)0.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I1H)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0, 90, 0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)-0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  A * vectorR3((FDOUBLE)1.,(FDOUBLE)0.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I3H)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)0., (FDOUBLE)1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)-0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3((FDOUBLE)0.,(FDOUBLE)1.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   >= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   >= 0 &&
                dotProduct(directions_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_2_fold_axis)     >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4H)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,-31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)0., (FDOUBLE)1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)-0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3((FDOUBLE)0.,(FDOUBLE)1.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
                dotProduct(directions_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) <= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_2_fold_axis)     >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5H)
    {//OK
        std::cerr << "ERROR: pg_I5H Symmetry not implemented" << std::endl;
        exit(0);
    }
    else
    {
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }
    
    
    // Now overwrite the rot/tilt_angles and directions_vectors with their no_redundant counterparts
    rot_angles = no_redundant_rot_angles;
    tilt_angles = no_redundant_tilt_angles;
    directions_vector = no_redundant_directions_vector;
    directions_ipix = no_redundant_directions_ipix;
    
    
}

