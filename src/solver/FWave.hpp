#ifndef SOLVER_FWAVE_HPP
#define SOLVER_FWAVE_HPP

#include <cmath>
#include <cassert>

namespace solver {
    /**
     * fWave solver
     */
    template <typename T> class FWave {
        private:
            /**
             * Gravity in m/s^2
             */
            T gravity;
         
            /**
             * Numerical tolerance for comparisons
             */
            T tolerance;
         
            /**
             * Left-side water height
             */
            T h_l;
         
            /**
             * Right-side water height
             */
            T h_r;
         
            /**
             * Left-side particle velocity
             */
            T u_l;
         
            /**
             * Right-side particle velocity
             */
            T u_r;
         
            /**
             * Left-side momentum
             */
            T hu_l;
         
            /**
             * Right-side momentum
             */
            T hu_r;
         
            /**
             * First Roe eigenvalue
             */
            T lambda_1;
         
            /**
             * Second Roe eigenvalue
             */
            T lambda_2;
         
            /**
             * First component of flux jump
             */
            T delta_f_1;
         
            /**
             * Second component of flux jump
             */
            T delta_f_2;
         
            /**
             * First eigencoefficient in wave decomposition corresponding to first eigenvalue
             */
            T alpha_1;
         
            /**
             * Second eigencoefficient in wave decomposition corresponding to second eigenvalue
             */
            T alpha_2;
         
            /**
             * Store input parameters as instance variables
             * 
             * @param[in] h_l   Left-side water height
             * @param[in] h_r   Right-side water height
             * @param[in] hu_l  Left-side momentum
             * @param[in] hu_r  Right-side momentum
             */
            void storeParameters(T h_l, T h_r, T hu_l, T hu_r) {
                this->h_l = h_l;
                this->h_r = h_r;
                this->hu_l = hu_l;
                this->hu_r = hu_r;
             
                if(h_l > tolerance)
                    this->u_l = hu_l / h_l;
                else
                    this->u_l = 0.0;
             
                if(h_r > tolerance)
                    this->u_r = hu_r / h_r;
                else
                    this->u_r = 0.0;
            }
         
           /**
            * Compute eigenvalues used in wave decomposition and store in instance variables
            */
            void computeEigenvalues() {
                T velocity = computeParticleVelocity(h_l, h_r, u_l, u_r, hu_l, hu_r);
                T height = sqrt( gravity * computeHeight(h_l, h_r) );
                lambda_1 = velocity - height;
                lambda_2 = velocity + height;
            }
         
           /**
            * Compute particle velocity used in eigenvalue computation
            *
            * @param[in]   h_l         Left-side water height
            * @param[in]   h_r         Right-side water height
            * @param[in]   u_l         Left-side velocity
            * @param[in]   u_r         Right-side velocity
            * @param[in]   hu_l        Left-side momentum
            * @param[in]   hu_r        Right-side momentum
            * @return                  Particle velocity
            */
            static T computeParticleVelocity( T h_l, T h_r, T u_l, T u_r, T hu_l, T hu_r) {
                return ( u_l * sqrt(h_l) + u_r * sqrt(h_r) ) / ( sqrt(h_l) + sqrt(h_r) );
            }
         
           /**
            * Computes water height used in eigenvalue computation
            *
            * @param[in]   h_l         Left water height
            * @param[in]   h_r         Right water height
            * @return                  Water height
            */
            static T computeHeight(T h_l, T h_r) {
                return 0.5 * ( h_l + h_r );
            }
              
            /**
             * Compute the flux jump as a function of left and right water height and momentum.
             * Results are stored in instance variables delta_f_1 and delta_f_2
             */
            void computeFluxJump() {
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
                delta_f_2 = (hu_r * u_r) - (hu_l * u_l) + 0.5 * gravity * (h_r*h_r - h_l*h_l);
            }
         
            /**
             * Compute the eigencoefficients needed to compute the resulting waves.
             * Results are stored in instance variables alpha_1 and alpha_2
             */
            void computeEigencoefficients() {
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
             
                computeFluxJump();
             
                T lambda_diff = lambda_2 - lambda_1;
                assert(lambda_diff > tolerance || lambda_diff < -tolerance); // lambda_diff must not be zero
             
                alpha_1 = (lambda_2 * delta_f_1 - delta_f_2) / lambda_diff;
                alpha_2 = (-lambda_1 * delta_f_1 + delta_f_2) / lambda_diff;
            }
        public:
            /**
             * Constructor
             *
             * @param[in]   _gravity            Gravity constant in m/s^2
             * @param[in]   _tolerance          Numerical tolerance
             */
            FWave(T _gravity = 9.81, T _tolerance = 1e-10):
                gravity(_gravity), tolerance(_tolerance) {}
         
            /**
             * Compute net updates and maximum wavespeed at an edge
             *
             * @param[in]   h_l                 Left-side water height
             * @param[in]   h_r                 Right-side water height
             * @param[in]   hu_l                Left-side momentum
             * @param[in]   hu_r                Right-side momentum
             * @param[in]   b_l                 Left-side bathymetry
             * @param[in]   b_r                 Right-side bathymetry
             * @param[out]  netUpdateLeft_h     Left net update (water height)
             * @param[out]  netUpdateRight_h    Right net update (water height)
             * @param[out]  netUpdateLeft_hu    Left net update (momentum)
             * @param[out]  netUpdateRight_hu   Right net update (momentum)
             * @param[out]  maxWaveSpeed        Maximum wave speed
             */
            void computeNetUpdates(
                T h_l,
                T h_r,
                T hu_l,
                T hu_r,
                T b_l,
                T b_r,
                T &netUpdateLeft_h,
                T &netUpdateRight_h,
                T &netUpdateLeft_hu,
                T &netUpdateRight_hu,
                T &maxWaveSpeed
            ) {
                // Height cannot be negative
                assert(h_l >= 0);
                assert(h_r >= 0);
             
                // Store params
                // TODO: store and obey bathymetry
                storeParameters(h_l, h_r, hu_l, hu_r);
             
                computeEigenvalues();
                computeEigencoefficients();
             
                netUpdateLeft_h = 0.0;
                netUpdateLeft_hu = 0.0;
                netUpdateRight_h = 0.0;
                netUpdateRight_hu = 0.0;
                T waveSpeedLeft = 0.0;
                T waveSpeedRight = 0.0;
         
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
             
                // TODO: set maxWaveSpeed correctly
         
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
                if(lambda_1 > tolerance) {
                    // to right
                    netUpdateRight_h += alpha_1;
                    netUpdateRight_hu += alpha_1 * lambda_1;
                } else if(lambda_1 < -tolerance) {
                    // to left
                    netUpdateLeft_h += alpha_1;
                    netUpdateLeft_hu += alpha_1 * lambda_1;
                } else {
                    // TODO: Fix behavior
                    // lambda is (numerically) zero, which implies hu is zero
                    // however, alpha may be not be zero, so we need to handle that case
                    // in some sensible way
                    assert(alpha_1 < tolerance && alpha_1 > -tolerance);
                }
         
                if(lambda_2 > tolerance) {
                    // to right
                    netUpdateRight_h += alpha_2;
                    netUpdateRight_hu += alpha_2 * lambda_2;
                } else if(lambda_2 < -tolerance) {
                    // to left
                    netUpdateLeft_h += alpha_2;
                    netUpdateLeft_hu += alpha_2 * lambda_2;
                } else {
                    // TODO: Fix behavior
                    // lambda is (numerically) zero, which implies hu is zero
                    // however, alpha may be not be zero, so we need to handle that case
                    // in some sensible way
                    assert(alpha_2 < tolerance && alpha_2 > -tolerance);
                }
            }
    };
};

#endif
