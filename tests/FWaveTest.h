
#ifndef FWAVETEST_H
#define FWAVETEST_h

#include <cxxtest/TestSuite.h>

#define private public
#define protected public
#include <FWave.h>

class FWaveTest : public CxxTest::TestSuite {
    public:
        void testParticleVelocity() {
            // h_l = 16
            // u_l = 1.5
            // h_r = 9
            // u_r = -0.5
            TS_ASSERT_EQUALS(
                FWave::computeParticleVelocity(16.0, 24.0, 9.0, -4.5),
                (4.5 / 7.0)
            );
        }
};

#endif
