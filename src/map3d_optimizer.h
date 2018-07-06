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
#if defined(ROME_MAP2D)
#error "map3d used in ROME2D project"
#endif
#ifndef MAP3D_OPTIMIZER_H_
#define MAP3D_OPTIMIZER_H_

#include "./sampling.h"
#include "./mrcs.h"
#include "./image.h"
#include "./fft_fftw3.h"
#include "./ctf.h"
#include "./time.h"
#include "./mpi.h"
#include "./metadata.h"
#include "./initialize.h"
#include "./string.h"
#include "./exp_model.h"
#include "./map_model.h"
#include "./ml_model.h"
#include "./progressbar.h"
#include "./statusTracer.h"
#include "./map_optimizer_base_old.h"

namespace Map3dOptimizer
{

    
};


#endif
