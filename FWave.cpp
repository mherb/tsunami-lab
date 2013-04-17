#include "FWave.h"
#include <math.h>

using namespace std;


GRAVITY = 9.81;

static void solve(double h_l, double hu_l, double h_r, double hu_r, double &netUpdateLeft_h, double &netUpdateLeft_hu, double &netUpdateRight_h, double &netUpdateRight_hu, double &waveSpeedLeft, double &waveSpeedRight){
	double lambda_1;
	double lambda_2;
	computeEigenvalues (h_l, hu_l, h_r, hu_r, lambda_1, lambda_2);
} 


static double computeHeight (double h_l, double h_r){
	return 0.5*(h_l+h_r);
}


static double computeParticleVelocity(double h_l, double hu_l, double h_r, double hu_r){
	return (hu_l/sqrt(h_l) + (hu_r/sqrt(h_r)))/(sqrt(h_l)+sqrt(h_r));
}


static void computeEigenvalues (double h_l, double hu_l, double h_r, double hu_r, double &lambda_1, double &lambda_2){
	lambda_1 = (computeParticleVelocity(h_l, hu_l, h_r, hu_r) - sqrt(GRAVITY*computeHeight(h_l, h_r)));
	lambda_2 = (computeParticleVelocity(h_l, hu_l, h_r, hu_r) + sqrt(GRAVITY*computeHeight(h_l, h_r)));
}
