
#ifndef SOLVER_FWAVETEST_H
#define SOLVER_FWAVETEST_H

#include <cxxtest/TestSuite.h>

#define private public
#define protected public
#include <solver/FWave.hpp>

typedef double T;

/**
 * Unit Testing class for FWave
 */
class FWaveTest : public CxxTest::TestSuite {
    private:
        /**
         * Numerical tolerance for assertions
         */
        const static T TOLERANCE = 1e-10;
        
        solver::FWave<T> fwave;
        
    public:
        
        /**
         * Test private storeParameters()
         */
        void testStoreParameters() {
            fwave.storeParameters(16.0, 9.0, 24.0, -4.5, -16.0, -16.0);
            TSM_ASSERT_DELTA("[Regular] Left water heigt", fwave.h_l, 16.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Right water heigt", fwave.h_r, 9.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Left velocity", fwave.u_l, 1.5, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Right velocity", fwave.u_r, -0.5, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Left momentum", fwave.hu_l, 24.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Right momentum", fwave.hu_r, -4.5, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Left bathymetry", fwave.b_l, -16.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Right bathymetry", fwave.b_r, -16.0, TOLERANCE);
            
            fwave.storeParameters(0.0, 9.0, 0.0, -4.5, 0.0, -9.0);
            TSM_ASSERT_DELTA("[ZeroHeight] Left water heigt", fwave.h_l, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Right water heigt", fwave.h_r, 9.0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Left velocity", fwave.u_l, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Right velocity", fwave.u_r, -0.5, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Left momentum", fwave.hu_l, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Right momentum", fwave.hu_r, -4.5, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Left bathymetry", fwave.b_l, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeight] Right bathymetry", fwave.b_r, -9.0, TOLERANCE);
        }
        
       /**
        * Test numerical flux jump computation
        */
        void testComputeFluxJump() {
            // h_l = 16
            // hu_l = 24
            // h_r = 9
            // hu_r = -4.5
            fwave.storeParameters(16.0, 9.0, 24.0, -4.5, -16.0, -16.0);
            fwave.computeFluxJump();
        
            TSM_ASSERT_DELTA("First component", fwave.delta_f_1, -28.5, TOLERANCE);
            TSM_ASSERT_DELTA("Second component", fwave.delta_f_2, -892.125, TOLERANCE);
        
            // h_l = 12
            // hu_l = -15
            // h_r = 6
            // hu_r = 5.5
            fwave.storeParameters(12.0, 6.0, -15.0, 5.5, -12.0, -12.0);
            fwave.computeFluxJump();
        
            TSM_ASSERT_DELTA("First component", fwave.delta_f_1, 20.5, TOLERANCE);
            TSM_ASSERT_DELTA("Second component", fwave.delta_f_2, -543.44833333333333333333, TOLERANCE);
        }
        
       /**
        * Test computation of eigencoefficients
        */
        void testComputeEigencoefficients() {
            // h_l = 6
            // hu_l = 13
            // h_r = 6.5
            // hu_r = -3.5
            
            fwave.storeParameters(6.0, 6.5, 13.0, -3.5, -6.0, -6.0);
            fwave.lambda_1 = -7.04318942880952746985;
            fwave.lambda_2 = 8.61727033455629779978;
            
            fwave.computeEigencoefficients();
            
            TSM_ASSERT_DELTA("First eigencoefficient", fwave.alpha_1, -9.35854767054606546293, TOLERANCE);
            TSM_ASSERT_DELTA("Second eigencoefficient", fwave.alpha_2, -7.14145232945393453707, TOLERANCE);
        }
        
        /**
         * Test computation of net updates
         */
        void testComputeNetUpdates() {
            T netUpdateLeft_h, netUpdateLeft_hu;
            T netUpdateRight_h, netUpdateRight_hu;
            T waveSpeedLeft, waveSpeedRight;
            T maxWaveSpeed;
            
            // Regular: Lambda1 < 0, Lambda2 > 0
            fwave.computeNetUpdates(10.0, 12.5, 5.0, -3.5, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[Regular] Net Update Left Height", netUpdateLeft_h, -17.34505918808096570196, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Left Momentum", netUpdateLeft_hu, 180.68503796971846054652, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Right Height", netUpdateRight_h, 8.84505918808096570196, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Net Update Right Momentum", netUpdateRight_hu, 93.70121203028153945348, TOLERANCE);
            TSM_ASSERT_DELTA("[Regular] Max Wave Speed", maxWaveSpeed, 10.59362182183554874211, TOLERANCE);
            
            // SupersonicRight: Lambda1, Lambda2 > 0
            fwave.computeNetUpdates(4.5, 2.5, 20.0, 22.5, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Right Height", netUpdateRight_h, 2.5, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Net Update Right Momentum", netUpdateRight_hu, 44.94111111111111111111, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicRight] Max Wave Speed", maxWaveSpeed, 12.24950641851166448956, TOLERANCE);
            
            // SupersonicLeft: Lambda1, Lambda2 < 0
            fwave.computeNetUpdates(7.5, 1.4, -27.3, -25.2, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Left Height", netUpdateLeft_h, 2.1, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Left Momentum", netUpdateLeft_hu, 87.93555, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[SupersonicLeft] Max Wave Speed", maxWaveSpeed, 14.57956803440405980804, TOLERANCE);
            
            // Steady state
            fwave.computeNetUpdates(12.0, 12.0, 14.0, 14.0, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[Steady] Net Update Left Height", netUpdateLeft_h, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Left Momentum", netUpdateLeft_hu, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Right Height", netUpdateRight_h, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Net Update Right Momentum", netUpdateRight_hu, 0.0, TOLERANCE);
            TSM_ASSERT_DELTA("[Steady] Max Wave Speed", maxWaveSpeed, 12.0165514586817413307, TOLERANCE);
            
            // Lambda1 = 0, Lambda2 > 0
            T h = 5.0;
            T hu = h * sqrt(fwave.gravity * h);
            fwave.computeNetUpdates(h, h, hu, hu, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[LambdaZero] Max Wave Speed", maxWaveSpeed, 14.00714103591450242095, TOLERANCE);
            
            // Left Height = 0
            fwave.computeNetUpdates(0.0, 5.0, 0.0, 2.5, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[ZeroHeightLeft] Net Update Left Height", netUpdateLeft_h, -11.13068051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightLeft] Net Update Left Momentum", netUpdateLeft_hu, 49.55681948558561664715, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightLeft] Net Update Right Height", netUpdateRight_h, 13.63068051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightLeft] Net Update Right Momentum", netUpdateRight_hu, 74.31818051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightLeft] Max Wave Speed", maxWaveSpeed, 5.45227220576575334114, TOLERANCE);
            
            // Right Height = 0
            fwave.computeNetUpdates(5.0, 0.0, 2.5, 0.0, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[ZeroHeightRight] Net Update Left Height", netUpdateLeft_h, 11.13068051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightRight] Net Update Left Momentum", netUpdateLeft_hu, -49.55681948558561664715, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightRight] Net Update Right Height", netUpdateRight_h, -13.63068051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightRight] Net Update Right Momentum", netUpdateRight_hu, -74.31818051441438335285, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightRight] Max Wave Speed", maxWaveSpeed, 5.45227220576575334114, TOLERANCE);
            
            // Both Height = 0
            fwave.computeNetUpdates(0.0, 0.0, 2.5, 1.5, -50.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[ZeroHeightBoth] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightBoth] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightBoth] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightBoth] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[ZeroHeightBoth] Max Wave Speed", maxWaveSpeed, 0, TOLERANCE);
            
            // Left Bathymetry >= 0 (Reflecting Boundary)
            // h_l, hu_l is ignored
            fwave.computeNetUpdates(10.0, 5.0, 10.0, -2.5, 0.0, -50.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[DryLeft] Net Update Left Height", netUpdateLeft_h, -2.5, TOLERANCE);
            TSM_ASSERT_DELTA("[DryLeft] Net Update Left Momentum", netUpdateLeft_hu, 17.50892629489312802619, TOLERANCE);
            TSM_ASSERT_DELTA("[DryLeft] Net Update Right Height", netUpdateRight_h, -2.5, TOLERANCE);
            TSM_ASSERT_DELTA("[DryLeft] Net Update Right Momentum", netUpdateRight_hu, -17.50892629489312802619, TOLERANCE);
            TSM_ASSERT_DELTA("[DryLeft] Max Wave Speed", maxWaveSpeed, 7.00357051795725121047, TOLERANCE);
            
            // Right Bathymetry >= 0 (Reflecting Boundary)
            // h_r, hu_r is ignored
            fwave.computeNetUpdates(12.5, 5.0, 6.5, 10.0, -50.0, 1.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[DryRight] Net Update Left Height", netUpdateLeft_h, -6.5, TOLERANCE);
            TSM_ASSERT_DELTA("[DryRight] Net Update Left Momentum", netUpdateLeft_hu, 71.97851241863782780576, TOLERANCE);
            TSM_ASSERT_DELTA("[DryRight] Net Update Right Height", netUpdateRight_h, -6.5, TOLERANCE);
            TSM_ASSERT_DELTA("[DryRight] Net Update Right Momentum", netUpdateRight_hu, -71.97851241863782780576, TOLERANCE);
            TSM_ASSERT_DELTA("[DryRight] Max Wave Speed", maxWaveSpeed, 11.07361729517505043166, TOLERANCE);
            
            // Both Bathymetry >= 0
            fwave.computeNetUpdates(4.5, 3.5, 2.5, 1.5, 0.0, 0.0,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    maxWaveSpeed
                );
            TSM_ASSERT_DELTA("[DryBoth] Net Update Left Height", netUpdateLeft_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[DryBoth] Net Update Left Momentum", netUpdateLeft_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[DryBoth] Net Update Right Height", netUpdateRight_h, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[DryBoth] Net Update Right Momentum", netUpdateRight_hu, 0, TOLERANCE);
            TSM_ASSERT_DELTA("[DryBoth] Max Wave Speed", maxWaveSpeed, 0, TOLERANCE);
    }
};

#endif
