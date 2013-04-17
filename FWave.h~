#ifndef FWAVE_H
#define FWAVE_H

class FWave {
    private: 
        static const double GRAVITY;
        static void computeEigenvalues(
            double h_l,
            double hu_l,
            double h_r,
            double hu_r,
            double &lambda_1,
            double &lambda_2
        );
        static void computeNetUpdates(
            double h_l,
            double hu_l,
            double h_r,
            double hu_r,
            double &netUpdateLeft_h,
            double &netUpdateLeft_hu,
            double &netUpdateRight_h,
            double &netUpdateRight_hu
        );
    public:
        static void solve(
            double h_l,
            double hu_l,
            double h_r,
            double hu_r,
            double &netUpdateLeft_h,
            double &netUpdateLeft_hu,
            double &netUpdateRight_h,
            double &netUpdateRight_hu,
            double &waveSpeedLeft,
            double &waveSpeedRight
        );
};

#endif