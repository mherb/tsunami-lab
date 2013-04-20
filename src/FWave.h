#ifndef FWAVE_H
#define FWAVE_H

class FWave {
    private:
        /**
         * Gravity constant in m/s^2
         */
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
        static double computeParticleVelocity(
            double h_l,
            double hu_l,
            double h_r,
            double hu_r
        );
        static double computeHeight (
            double h_l,
            double h_r
        );
                
        /**
         * Computes the flux jump as a function of left and right water height and momentum.
         *
         * @param[in]   h_l         Left water height
         * @param[in]   hu_l        Left momentum
         * @param[in]   h_r         Right water height
         * @param[in]   hu_l        Right momentum
         * @param[out]  delta_f_1   First component of flux jump
         * @param[out]  delta_f_2   Second component of flux jump
         */
        static void computeFluxJump(
            double h_l,
            double hu_l,
            double h_r,
            double hu_r,
            double &delta_f_1,
            double &delta_f_2
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
