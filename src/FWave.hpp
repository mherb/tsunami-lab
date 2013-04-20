#ifndef FWAVE_HPP
#define FWAVE_HPP

#include <cmath>
#include <cassert>

/**
 * fWave solver
 */
class FWave {
    private:
        /**
         * Gravity in m/s^2
         */
        const static double GRAVITY;
        
        /**
         * Numerical tolerance for comparisons
         */
        const static double TOLERANCE;
        
        static void computeEigenvalues(double h_l, double hu_l, double h_r, double hu_r, double &lambda_1, double &lambda_2) {
            lambda_1 = ( computeParticleVelocity(h_l, hu_l, h_r, hu_r) - sqrt( GRAVITY * computeHeight(h_l, h_r) ) );
            lambda_2 = ( computeParticleVelocity(h_l, hu_l, h_r, hu_r) + sqrt( GRAVITY * computeHeight(h_l, h_r) ) );
        }
        
        static double computeParticleVelocity( double h_l, double hu_l, double h_r, double hu_r) {
            return ( hu_l / sqrt(h_l) + ( hu_r / sqrt(h_r) ) ) / ( sqrt(h_l) + sqrt(h_r) );
        }
        
        static double computeHeight(double h_l, double h_r) {
            return 0.5 * ( h_l + h_r );
        }
             
        /**
         * Computes the flux jump as a function of left and right water height and momentum.
         *
         * @param[in]   h_l         Left water height
         * @param[in]   hu_l        Left momentum
         * @param[in]   h_r         Right water height
         * @param[in]   hu_r        Right momentum
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
        ) {
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
            delta_f_1 = hu_r - hu_l;
            delta_f_2 = (hu_r * (hu_r / h_r)) - (hu_l * (hu_l / h_l)) + 0.5 * GRAVITY * (h_r*h_r - h_l*h_l);
        }
        
        /**
         * Computes the eigencoefficients needed to compute the resulting waves
         *
         * @param[in]   h_l         Left water height
         * @param[in]   hu_l        Left momentum
         * @param[in]   h_r         Right water height
         * @param[in]   hu_r        Right momentum
         * @param[in]   lambda_1    First eigenvalue
         * @param[in]   lambda_2    Second eigenvalue
         * @param[out]  alpha_1     First eigencoefficient
         * @param[out]  alpha_2     Second eigencoefficient
         */
        static void computeEigencoefficients(
            double h_l,
            double hu_l,
            double h_r,
            double hu_r,
            double lambda_1,
            double lambda_2,
            double &alpha_1,
            double &alpha_2
        ) {
           /**
            * **Compute Eigencoefficients**
            * \f[
            * \begin{bmatrix}
            *     \alpha_1 \\
            *     \alpha_2
            * \end{bmatrix} =
            * R^{-1} \cdot \Delta f
            * \text{ with }
            * R = \begin{bmatrix}
            *     1 & 1 \\
            *     \lambda_1 & \lambda_2
            * \end{bmatrix},
            * \Delta f = \begin{bmatrix}
            *     \Delta f_1 \\
            *     \Delta f_2
            * \end{bmatrix}
            * \f]
            * 
            * With 
            * \f[
            * R^{-1} = 
            * \frac{1}{\lambda_2 - \lambda_1} \begin{bmatrix}
            *     \lambda_2 & -1 \\
            *     -\lambda_1 & 1
            * \end{bmatrix}
            * \f]
            * we end up with
            * \f{align*}{
            * \begin{bmatrix}
            *     \alpha_1 \\
            *     \alpha_2
            * \end{bmatrix}
            * &= 
            * \frac{1}{\lambda_2 - \lambda_1} \begin{bmatrix}
            *     \lambda_2 & -1 \\
            *     -\lambda_1 & 1
            * \end{bmatrix}
            * \cdot
            * \begin{bmatrix}
            *     \Delta f_1 \\
            *     \Delta f_2
            * \end{bmatrix}
            * \\
            * &=
            * \frac{1}{\lambda_{\text{diff}}} \begin{bmatrix}
            *     \lambda_2 \Delta f_1 - \Delta f_2 \\
            *     -\lambda_1 \Delta f_1 + \Delta f_2
            * \end{bmatrix},
            * \qquad \lambda_{\text{diff}} = \lambda_2 - \lambda_1
            * \f}
            */
         
            double delta_f_1, delta_f_2;
            computeFluxJump(h_l, hu_l, h_r, hu_r, delta_f_1, delta_f_2);
        
            double lambda_diff = lambda_2 - lambda_1;
            assert(lambda_diff > TOLERANCE || lambda_diff < -TOLERANCE); // lambda_diff must not be zero
        
            alpha_1 = (lambda_2 * delta_f_1 - delta_f_2) / lambda_diff;
            alpha_2 = (-lambda_1 * delta_f_1 + delta_f_2) / lambda_diff;
        }
    public:
       
        /**
         * Compute net updates and resulting wavespeeds at an edge
         *
         * @param[in]   h_l                 Left water height
         * @param[in]   hu_l                Left momentum
         * @param[in]   h_r                 Right water height
         * @param[in]   hu_r                Right momentum
         * @param[out]  netUpdateLeft_h     Left net update (water height)
         * @param[out]  netUpdateLeft_hu    Left net update (momentum)
         * @param[out]  netUpdateRight_h    Right net update (water height)
         * @param[out]  netUpdateRight_hu   Right net update (momentum)
         * @param[out]  waveSpeedLeft       Left wave speed
         * @param[out]  waveSpeedRight      Right wave speed
         */
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
        ) {
            // Height cannot be negative
            assert(h_l >= 0);
            assert(h_r >= 0);
        
            double lambda_1, lambda_2;
            double alpha_1, alpha_2;
            computeEigenvalues(h_l, hu_l, h_r, hu_r, lambda_1, lambda_2);
            computeEigencoefficients(h_l, hu_l, h_r, hu_r, lambda_1, lambda_2, alpha_1, alpha_2);
        
            netUpdateLeft_h = 0.0;
            netUpdateLeft_hu = 0.0;
            netUpdateRight_h = 0.0;
            netUpdateRight_hu = 0.0;
            waveSpeedLeft = 0.0;
            waveSpeedRight = 0.0;
        
            /**
             * **Compute Wavespeeds**
             *
             * \f[
             * \begin{bmatrix} \lambda_l \\ \lambda_r \end{bmatrix} =
             * \begin{cases}
             *     \begin{bmatrix} \lambda_1 & 0 \end{bmatrix}^{\mathsf T} & \lambda_1 > 0 \land \lambda_2 > 0 \\
             *     \begin{bmatrix} 0 & \lambda_2 \end{bmatrix}^{\mathsf T} & \lambda_1 < 0 \land \lambda_2 < 0 \\
             *     \begin{bmatrix} \lambda_1 & \lambda_2 \end{bmatrix}^{\mathsf T} & \text{else} \\
             * \end{cases}
             * \f]
             */
            if(std::signbit(lambda_1) == std::signbit(lambda_2)) {
                // same sign
                if(std::signbit(lambda_1)) {
                    // both negative
                    waveSpeedLeft = lambda_1;
                } else {
                    // both positive
                    waveSpeedRight = lambda_2;
                }
            } else {
                // opposite sign
                waveSpeedLeft = lambda_1;
                waveSpeedRight = lambda_2;
            }
        
            /**
             * **Compute net updates**
             *
             * \f{align*}{
             * A^-\Delta Q &=
             * \sum\limits_{p:\lambda_p < 0} Z_p =
             * \sum\limits_{p:\lambda_p < 0} \begin{bmatrix} \alpha_p \\ \alpha_p\lambda_p\end{bmatrix} =
             * \begin{bmatrix} \Delta h_l \\ \Delta hu_l \end{bmatrix} \\
             * A^+\Delta Q &=
             * \sum\limits_{p:\lambda_p > 0} Z_p =
             * \sum\limits_{p:\lambda_p > 0} \begin{bmatrix} \alpha_p \\ \alpha_p\lambda_p\end{bmatrix} =
             * \begin{bmatrix} \Delta h_r \\ \Delta hu_r \end{bmatrix} \\
             * \f}
             */
            if(lambda_1 > TOLERANCE) {
                // to right
                netUpdateRight_h += alpha_1;
                netUpdateRight_hu += alpha_1 * lambda_1;
            } else if(lambda_1 < -TOLERANCE) {
                // to left
                netUpdateLeft_h += alpha_1;
                netUpdateLeft_hu += alpha_1 * lambda_1;
            } else {
                // TODO: Fix behavior
                // lambda is (numerically) zero, which implies hu is zero
                // however, alpha may be not be zero, so we need to handle that case
                // in some sensible way
                assert(alpha_1 < TOLERANCE && alpha_1 > -TOLERANCE);
            }
        
            if(lambda_2 > TOLERANCE) {
                // to right
                netUpdateRight_h += alpha_2;
                netUpdateRight_hu += alpha_2 * lambda_2;
            } else if(lambda_2 < -TOLERANCE) {
                // to left
                netUpdateLeft_h += alpha_2;
                netUpdateLeft_hu += alpha_2 * lambda_2;
            } else {
                // TODO: Fix behavior
                // lambda is (numerically) zero, which implies hu is zero
                // however, alpha may be not be zero, so we need to handle that case
                // in some sensible way
                assert(alpha_2 < TOLERANCE && alpha_2 > -TOLERANCE);
            }
        }
};

const double FWave::GRAVITY = 9.81;
const double FWave::TOLERANCE = 1e-10;

#endif
