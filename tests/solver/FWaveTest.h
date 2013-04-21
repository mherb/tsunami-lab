
#ifndef SOLVER_FWAVETEST_H
#define SOLVER_FWAVETEST_H

#include <cxxtest/TestSuite.h>

#define private public
#define protected public
#include <solver/FWave.hpp>

typedef double T;

class FWaveTest : public CxxTest::TestSuite {
    private:
        /**
         * Numerical tolerance for assertions
         */
        const static T TOLERANCE = 1e-10;
        
        solver::FWave<T> fwave;
        
    public:
        
        void testStoreParameters() {
            fwave.storeParameters(16.0, 9.0, 24.0, -4.5);
            TSM_ASSERT_DELTA("[Regular] Left water heigt", fwave.h_l, 16.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Right water heigt", fwave.h_r, 9.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Left velocity", fwave.u_l, 1.5, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Right velocity", fwave.u_r, -0.5, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Left momentum", fwave.hu_l, 24.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Right momentum", fwave.hu_r, -4.5, TOLERANCE);
            
            fwave.storeParameters(0.0, 9.0, 0.0, -4.5);
            TSM_ASSERT_DELTA("[ZeroHeight] Left water heigt", fwave.h_l, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Right water heigt", fwave.h_r, 9.0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Left velocity", fwave.u_l, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Right velocity", fwave.u_r, -0.5, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Left momentum", fwave.hu_l, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Right momentum", fwave.hu_r, -4.5, TOLERANCE);
        }
        
        void testParticleVelocity() {
            // h_l = 16
            // u_l = 1.5
            // h_r = 9
            // u_r = -0.5
            TS_ASSERT_DELTA(
                fwave.computeParticleVelocity(16.0, 9.0, 1.5, -0.5, 24.0, -4.5),
                0.6428571428571428571,
                TOLERANCE
            );
        }
        
        void testComputeFluxJump() {
            // h_l = 16
            // hu_l = 24
            // h_r = 9
            // hu_r = -4.5
            fwave.storeParameters(16.0, 9.0, 24.0, -4.5);
            fwave.computeFluxJump();
        
            TSM_ASSERT_DELTA("First component", fwave.delta_f_1, -28.5, TOLERANCE);
            TSM_ASSERT_DELTA("Second component", fwave.delta_f_2, -892.125, TOLERANCE);
        
            // h_l = 12
            // hu_l = -15
            // h_r = 6
            // hu_r = 5.5
            fwave.storeParameters(12.0, 6.0, -15.0, 5.5);
            fwave.computeFluxJump();
        
            TSM_ASSERT_DELTA("First component", fwave.delta_f_1, 20.5, TOLERANCE);
            TSM_ASSERT_DELTA("Second component", fwave.delta_f_2, -543.44833333333333333333, TOLERANCE);
        }
        
        void testComputeEigencoefficients() {
            // h_l = 6
            // hu_l = 13
            // h_r = 6.5
            // hu_r = -3.5
            
            fwave.storeParameters(6.0, 6.5, 13.0, -3.5);
            fwave.lambda_1 = -7.04318942880952746985;
            fwave.lambda_2 = 8.61727033455629779978;
            
            fwave.computeEigencoefficients();
            
            TSM_ASSERT_DELTA("First eigencoefficient", fwave.alpha_1, -9.35854767054606546293, TOLERANCE);
            TSM_ASSERT_DELTA("Second eigencoefficient", fwave.alpha_2, -7.14145232945393453707, TOLERANCE);
        }
        
        void testComputeNetUpdates() {
            T netUpdateLeft_h, netUpdateLeft_hu;
            T netUpdateRight_h, netUpdateRight_hu;
            T waveSpeedLeft, waveSpeedRight;
            T maxWaveSpeed;
            
            // Regular: Lambda1 < 0, Lambda2 > 0
            fwave.computeNetUpdates(10.0, 12.5, 5.0, -3.5, 0.0, 0.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[Regular] Net Update Left Height", netUpdateLeft_h, -17.34505918808096570196, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Left Momentum", netUpdateLeft_hu, 180.68503796971846054652, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Right Height", netUpdateRight_h, 8.84505918808096570196, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Right Momentum", netUpdateRight_hu, 93.70121203028153945348, TOLERANCE);
            // TSM_ASSERT_DELTA("[Regular] Wave Speed Left", waveSpeedLeft, -10.41708973203620488931, TOLERANCE);
            // TSM_ASSERT_DELTA("[Regular] Wave Speed Right", waveSpeedRight, 10.59362182183554874211, TOLERANCE);
            
            // SupersonicRight: Lambda1, Lambda2 > 0
            fwave.computeNetUpdates(4.5, 2.5, 20.0, 22.5, 0.0, 0.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Right Height", netUpdateRight_h, 2.5, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Right Momentum", netUpdateRight_hu, 44.94111111111111111111, TOLERANCE);
            // TSM_ASSERT_DELTA("[SupersonicRight] Wave Speed Left", waveSpeedLeft, 0, TOLERANCE);
            // TSM_ASSERT_DELTA("[SupersonicRight] Wave Speed Right", waveSpeedRight, 12.24950641851166448956, TOLERANCE);
            
            // SupersonicLeft: Lambda1, Lambda2 < 0
            fwave.computeNetUpdates(7.5, 1.4, -27.3, -25.2, 0.0, 0.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Left Height", netUpdateLeft_h, 2.1, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Left Momentum", netUpdateLeft_hu, 87.93555, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            // TSM_ASSERT_DELTA("[SupersonicLeft] Wave Speed Left", waveSpeedLeft, -14.57956803440405980804, TOLERANCE);
            // TSM_ASSERT_DELTA("[SupersonicLeft] Wave Speed Right", waveSpeedRight, 0, TOLERANCE);
            
            // Steady state
            fwave.computeNetUpdates(12.0, 12.0, 14.0, 14.0, 0.0, 0.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[Steady] Net Update Left Height", netUpdateLeft_h, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Left Momentum", netUpdateLeft_hu, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Right Height", netUpdateRight_h, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Right Momentum", netUpdateRight_hu, 0.0, TOLERANCE);
            // TSM_ASSERT_DELTA("[Steady] Wave Speed Left", waveSpeedLeft, -9.68321812534840799736, TOLERANCE);
            // TSM_ASSERT_DELTA("[Steady] Wave Speed Right", waveSpeedRight, 12.0165514586817413307, TOLERANCE);
            
            // Lambda1 = 0, Lambda2 > 0
            T h = 5.0;
            T hu = h * sqrt(fwave.gravity * h);
            fwave.computeNetUpdates(h, h, hu, hu, 0.0, 0.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            // TSM_ASSERT_DELTA("[LambdaZero] Wave Speed Left", waveSpeedLeft, 0, TOLERANCE);
            // TSM_ASSERT_DELTA("[LambdaZero] Wave Speed Right", waveSpeedRight, 14.00714103591450242095, TOLERANCE);
            
            // Height = 0
            fwave.computeNetUpdates(0.0, 5.0, 0.0, 2.5, 0.0, 0.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[ZeroHeight] Net Update Left Height", netUpdateLeft_h, -11.13068051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Net Update Left Momentum", netUpdateLeft_hu, 49.55681948558561664715, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Net Update Right Height", netUpdateRight_h, 13.63068051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Net Update Right Momentum", netUpdateRight_hu, 74.31818051441438335285, TOLERANCE);
            // TSM_ASSERT_DELTA("[ZeroHeight] Wave Speed Left", waveSpeedLeft, -4.45227220576575334114, TOLERANCE);
            // TSM_ASSERT_DELTA("[ZeroHeight] Wave Speed Right", waveSpeedRight, 5.45227220576575334114, TOLERANCE);
    }
};

#endif
