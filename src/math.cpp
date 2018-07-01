/***************************************************************************
 *
 * Authors: "Yongbei(Galow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
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

#include "math.h"
#include "time.h"
#include <random>

// Compute statistics.
//
// The average, standard deviation, minimum and maximum value are
// returned.
void computeStats(const FDOUBLE* data, int size, FDOUBLE& avg, FDOUBLE& stddev, FDOUBLE& minval, FDOUBLE& maxval)
{
    avg = 0;
    stddev = 0;
    
    minval = maxval = data[0];
    
    for (int i = 0; i < size; i++)
    {
        FDOUBLE Tval=data[i];
        FDOUBLE val=Tval;//static_cast< DOUBLE >(Tval);
        avg += val;
        stddev += val * val;
        
        if (Tval > maxval)
            maxval = Tval;
        else if (Tval < minval)
            minval = Tval;
    }
    
    avg /= size;
    
    if (size > 1)
    {
        stddev = stddev / size - avg * avg;
        stddev *= size / (size - 1);
        
        // Foreseeing numerical instabilities
        stddev = sqrt( /*static_cast< DOUBLE >*/(fabs(stddev)) );
    }
    else
        stddev = 0;
}

//  -----------  euler function  --------------
void Euler_angles2matrix(FDOUBLE alpha, FDOUBLE beta, FDOUBLE gamma,
                         FDOUBLE A[][3])
{
    FDOUBLE ca, sa, cb, sb, cg, sg;
    FDOUBLE cc, cs, sc, ss;
    
    alpha = deg2rad(alpha);
    beta  = deg2rad(beta);
    gamma = deg2rad(gamma);
    
    ca = cos(alpha);
    cb = cos(beta);
    cg = cos(gamma);
    sa = sin(alpha);
    sb = sin(beta);
    sg = sin(gamma);
    cc = cb * ca;
    cs = cb * sa;
    sc = sb * ca;
    ss = sb * sa;
    
    A[0][0] =  cg * cc - sg * sa;
    A[0][1] =  cg * cs + sg * ca;
    A[0][2] = -cg * sb;
    A[1][0] = -sg * cc - cg * sa;
    A[1][1] = -sg * cs + cg * ca;
    A[1][2] = sg * sb;
    A[2][0] =  sc;
    A[2][1] =  ss;
    A[2][2] = cb;
}

// Matrix --> Euler angles
//#define DEBUG_EULER

#if !defined(FLT_EPSILON)
#define FLT_EPSILON 1.19209e-07
#endif

void Euler_matrix2angles(const FDOUBLE A[][3],
                         FDOUBLE &alpha,FDOUBLE &beta, FDOUBLE &gamma)
{
    FDOUBLE abs_sb, sign_sb;
    
    abs_sb = sqrt(A[0][2] * A[0][2] + A[1][2] * A[1][2]);
    if (abs_sb > 16*FLT_EPSILON)
    {
        gamma = atan2(A[1][2], -A[0][2]);
        alpha = atan2(A[2][1], A[2][0]);
        
        if (abs(sin(gamma)) < FLT_EPSILON)
            sign_sb = sgn(-A[0][2] / cos(gamma));
        // if (sin(alpha)<FLT_EPSILON) sign_sb=sgn(-A[0][2]/cos(gamma));
        // else sign_sb=(sin(alpha)>0) ? sgn(A[2][1]):-sgn(A[2][1]);
        else
            sign_sb = (sin(gamma) > 0) ? sgn(A[1][2]) : -sgn(A[1][2]);
        beta  = atan2(sign_sb * abs_sb, A[2][2]);
    }
    else
    {
        if (sgn(A[2][2]) > 0)
        {
            // Let's consider the matrix as a rotation around Z
            alpha = 0;
            beta  = 0;
            gamma = atan2(-A[1][0], A[0][0]);
        }
        else
        {
            alpha = 0;
            beta  = rome_pi;
            gamma = atan2(A[1][0], -A[0][0]);
        }
    }
    
    gamma = rad2deg(gamma);
    beta  = rad2deg(beta);
    alpha = rad2deg(alpha);
    
#ifdef DEBUG_EULER
    std::cout << "abs_sb " << abs_sb << std::endl;
    std::cout << "A(1,2) " << A(1, 2) << " A(0,2) " << A(0, 2) << " gamma "
    << gamma << std::endl;
    std::cout << "A(2,1) " << A(2, 1) << " A(2,0) " << A(2, 0) << " alpha "
    << alpha << std::endl;
    std::cout << "sign sb " << sign_sb << " A(2,2) " << A(2, 2)
    << " beta " << beta << std::endl;
#endif
}
#undef FLT_EPSILON

void Euler_angles2direction(FDOUBLE alpha, FDOUBLE beta,FDOUBLE v[])
{
    FDOUBLE ca, sa, cb, sb;
    FDOUBLE sc, ss;
    
    alpha = deg2rad(alpha);
    beta  = deg2rad(beta);
    
    ca = cos(alpha);
    cb = cos(beta);
    sa = sin(alpha);
    sb = sin(beta);
    sc = sb * ca;
    ss = sb * sa;
    
    v[0] = sc;
    v[1] = ss;
    v[2] = cb;
}

#if defined(MAP3D_OLD_BASELINE) // YongBei baseline for Map2d and Map3d

class Random_generator::Internals {
public:
    float X1, X2;
    int call;
    Internals() : call(0) {
    }
    ~Internals() {
        fini();
    }
    void init(int seed)
    {
        fini();
        if (seed < 0)
            std::srand(static_cast <unsigned> (dtime()) );
        else
            std::srand(static_cast <unsigned> (seed) );
    }
    void fini(){
        return;
    }
    float randUniform01() {
        return static_cast<float>(std::rand())/(static_cast <float>(RAND_MAX));
    }
private:
};

#else // Bevin

class Random_generator::Internals {
public:
    float X1, X2;
    int call;
    Internals() : call(0), re(new std::default_random_engine()), ud(0.0, 1.0) {
    }
    ~Internals() {
        fini();
    }
    void init(int seed)
    {
        fini();
        if (seed < 0)
            re = new std::default_random_engine(rd());
        else
            re = new std::default_random_engine(seed);
    }
    void fini(){
        delete re; re = NULL;
    }
    float randUniform01() {
        return ud(*re);
    }
private:
    std::random_device          rd;
    std::default_random_engine *re;
    std::uniform_real_distribution<float> ud;
};

#endif

Random_generator::Random_generator() {
    _internals = new Internals;
}

Random_generator::~Random_generator() {
    delete _internals; _internals = NULL;
}

void Random_generator::init(int seed)
{
    _internals->init(seed);
}

float Random_generator::rnd_unif(float a, float b)
{
    if (a == b)
        return a;
    else
        return a + _internals->randUniform01()*(b-a);
}

// Gaussian distribution ---------------
float Random_generator::rnd_gaus(float mu, float sigma)
{
    // TODO - replace this with the C++ random standard one

    float U1, U2, W, mult;
    
    if (sigma == 0)
        return mu;
    
    auto& call = _internals->call;
    auto& X1   = _internals->X1;
    auto& X2   = _internals->X2;

    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (float) X2);
    }
    
    do
    {
        U1 = -1 + _internals->randUniform01() * 2;
        U2 = -1 + _internals->randUniform01() * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mu + sigma * (float) X1);
    
}


Random_generator dontShare_Random_generator;
