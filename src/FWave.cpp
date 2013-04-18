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
