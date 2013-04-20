
#ifndef FWAVETEST_H
#define FWAVETEST_H

#include <cxxtest/TestSuite.h>

#define private public
#define protected public
#include <FWave.h>

class FWaveTest : public CxxTest::TestSuite {
    private:
        /**
         * Numerical tolerance for assertions
         */
        const static double TOLERANCE = 1e-10;
    public:
        void testParticleVelocity() {
            // h_l = 16
            // u_l = 1.5
            // h_r = 9
            // u_r = -0.5
            TS_ASSERT_DELTA(
                FWave::computeParticleVelocity(16.0, 24.0, 9.0, -4.5),
                0.6428571428571428571,
                TOLERANCE
            );
        }

        void testComputeFluxJump() {
            double delta_f_1, delta_f_2;
            
            // h_l = 16
            // hu_l = 24
            // h_r = 9
            // hu_r = -4.5
            FWave::computeFluxJump(16.0, 24.0, 9.0, -4.5, delta_f_1, delta_f_2);
        
            TSM_ASSERT_DELTA("First component", delta_f_1, -28.5, TOLERANCE);
            TSM_ASSERT_DELTA("Second component", delta_f_2, -892.125, TOLERANCE);
        
            // h_l = 12
            // hu_l = -15
            // h_r = 6
            // hu_r = 5.5
            FWave::computeFluxJump(12.0, -15.0, 6.0, 5.5, delta_f_1, delta_f_2);
        
            TSM_ASSERT_DELTA("First component", delta_f_1, 20.5, TOLERANCE);
            TSM_ASSERT_DELTA("Second component", delta_f_2, -543.44833333333333333333, TOLERANCE);
        }
        
        void testComputeEigencoefficients() {
            // h_l = 6
            // hu_l = 13
            // h_r = 6.5
            // hu_r = -3.5
            
            double lambda_1, lambda_2;
            double alpha_1, alpha_2;
            lambda_1 = -7.04318942880952746985;
            lambda_2 = 8.61727033455629779978;
            
            FWave::computeEigencoefficients(6.0, 13.0, 6.5, -3.5, lambda_1, lambda_2, alpha_1, alpha_2);
            
            TSM_ASSERT_DELTA("First eigencoefficient", alpha_1, -9.35854767054606546293, TOLERANCE);
            TSM_ASSERT_DELTA("Second eigencoefficient", alpha_2, -7.14145232945393453707, TOLERANCE);
        }
};

#endif
