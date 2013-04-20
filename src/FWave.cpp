#include <cmath>
#include "FWave.h"

using namespace std;

const double FWave::GRAVITY = 9.81;

void FWave::solve(double h_l, double hu_l, double h_r, double hu_r, double &netUpdateLeft_h, double &netUpdateLeft_hu, double &netUpdateRight_h, double &netUpdateRight_hu, double &waveSpeedLeft, double &waveSpeedRight){
	double lambda_1;
	double lambda_2;
	computeEigenvalues (h_l, hu_l, h_r, hu_r, lambda_1, lambda_2);
} 

double FWave::computeHeight (double h_l, double h_r){
	return 0.5 * ( h_l + h_r );
}

double FWave::computeParticleVelocity(double h_l, double hu_l, double h_r, double hu_r){
	return ( hu_l / sqrt(h_l) + ( hu_r / sqrt(h_r) ) ) / ( sqrt(h_l) + sqrt(h_r) );
}

void FWave::computeEigenvalues (double h_l, double hu_l, double h_r, double hu_r, double &lambda_1, double &lambda_2){
	lambda_1 = ( computeParticleVelocity(h_l, hu_l, h_r, hu_r) - sqrt( GRAVITY * computeHeight(h_l, h_r) ) );
	lambda_2 = ( computeParticleVelocity(h_l, hu_l, h_r, hu_r) + sqrt( GRAVITY * computeHeight(h_l, h_r) ) );
}

/**
 * **Compute flux jump**
 * 
 * \f{align*}{ \Delta f := f(q_r) - f(q_l)
 * &=
 * \begin{bmatrix}
 * hu_r \\
 * hu^2_r + \frac{1}{2} g h_r
 * \end{bmatrix}   
 * -
 * \begin{bmatrix}
 * hu_l \\
 * hu^2_l + \frac{1}{2} g h_l
 * \end{bmatrix} 
 * \\
 * &=
 * \begin{bmatrix}
 * hu_r - hu_l \\
 * hu^2_r - hu^2_l + \frac{1}{2} g (h_r - h_l)
 * \end{bmatrix}
 * \\
 * &=
 * \begin{bmatrix}
 * \Delta f_1 \\
 * \Delta f_2
 * \end{bmatrix}\f}
 * where \f$g\f$ denotes the gravity constant.
 */
void FWave::computeFluxJump(double h_l, double hu_l, double h_r, double hu_r, double &delta_f_1, double &delta_f_2) {
    delta_f_1 = hu_r - hu_l;
    delta_f_2 = (hu_r * (hu_r / h_r)) - (hu_l * (hu_l / h_l)) + 0.5 * GRAVITY * (h_r*h_r - h_l*h_l);
}
