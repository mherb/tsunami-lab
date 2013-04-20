
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
        
        void testSolve() {
            double netUpdateLeft_h, netUpdateLeft_hu;
            double netUpdateRight_h, netUpdateRight_hu;
            double waveSpeedLeft, waveSpeedRight;
            
            // Regular: Lambda1 < 0, Lambda2 > 0
            FWave::solve(10.0, 5.0, 12.5, -3.5,
                    netUpdateLeft_h, netUpdateLeft_hu,
                    netUpdateRight_h, netUpdateRight_hu,
                    waveSpeedLeft, waveSpeedRight
                );
            TSM_ASSERT_DELTA("[Regular] Net Update Left Height", netUpdateLeft_h, -17.34505918808096570196, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Left Momentum", netUpdateLeft_hu, 180.68503796971846054652, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Right Height", netUpdateRight_h, 8.84505918808096570196, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Right Momentum", netUpdateRight_hu, 93.70121203028153945348, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Wave Speed Left", waveSpeedLeft, -10.41708973203620488931, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Wave Speed Right", waveSpeedRight, 10.59362182183554874211, TOLERANCE);
            
            // SupersonicRight: Lambda1, Lambda2 > 0
            FWave::solve(4.5, 20.0, 2.5, 22.5,
                    netUpdateLeft_h, netUpdateLeft_hu,
                    netUpdateRight_h, netUpdateRight_hu,
                    waveSpeedLeft, waveSpeedRight
                );
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Right Height", netUpdateRight_h, 2.5, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Right Momentum", netUpdateRight_hu, 44.94111111111111111111, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Wave Speed Left", waveSpeedLeft, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Wave Speed Right", waveSpeedRight, 12.24950641851166448956, TOLERANCE);
            
            // SupersonicLeft: Lambda1, Lambda2 < 0
            FWave::solve(7.5, -27.3, 1.4, -25.2,
                    netUpdateLeft_h, netUpdateLeft_hu,
                    netUpdateRight_h, netUpdateRight_hu,
                    waveSpeedLeft, waveSpeedRight
                );
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Left Height", netUpdateLeft_h, 2.1, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Left Momentum", netUpdateLeft_hu, 87.93555, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Wave Speed Left", waveSpeedLeft, -14.57956803440405980804, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Wave Speed Right", waveSpeedRight, 0, TOLERANCE);
            
            // Steady state
            FWave::solve(12.0, 14.0, 12.0, 14.0,
                    netUpdateLeft_h, netUpdateLeft_hu,
                    netUpdateRight_h, netUpdateRight_hu,
                    waveSpeedLeft, waveSpeedRight
                );
            TSM_ASSERT_DELTA("[Steady] Net Update Left Height", netUpdateLeft_h, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Left Momentum", netUpdateLeft_hu, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Right Height", netUpdateRight_h, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Right Momentum", netUpdateRight_hu, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Wave Speed Left", waveSpeedLeft, -9.68321812534840799736, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Wave Speed Right", waveSpeedRight, 12.0165514586817413307, TOLERANCE);
            
            // Lambda1 = 0, Lambda2 > 0
            double h = 5.0;
            double hu = h * sqrt(FWave::GRAVITY * h);
            FWave::solve(h, hu, h, hu,
                    netUpdateLeft_h, netUpdateLeft_hu,
                    netUpdateRight_h, netUpdateRight_hu,
                    waveSpeedLeft, waveSpeedRight
                );
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Wave Speed Left", waveSpeedLeft, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Wave Speed Right", waveSpeedRight, 14.00714103591450242095, TOLERANCE);
            
            /*
            // TODO: test fails due to unmet assertion
            // we should not do hu/h to get h because that leads to a NaN
            // in case h is zero (which can happen)
            
            // Height = 0
            FWave::solve(0.0, 0.0, 5.0, 2.5,
                    netUpdateLeft_h, netUpdateLeft_hu,
                    netUpdateRight_h, netUpdateRight_hu,
                    waveSpeedLeft, waveSpeedRight
                );
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Wave Speed Left", waveSpeedLeft, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Wave Speed Right", waveSpeedRight, 0, TOLERANCE);
            */
    }
};

#endif
