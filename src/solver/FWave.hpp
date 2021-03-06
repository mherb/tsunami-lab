#ifndef SOLVER_FWAVE_HPP
#define SOLVER_FWAVE_HPP

#include <cmath>
#include <cassert>
#include <algorithm>

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
             * Left-side bathymetry
             */
            T b_l;
         
            /**
             * Right-side bathymetry
             */
            T b_r;
         
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
             * @param[in] b_l  Left-side bathymetry
             * @param[in] b_r  Right-side bathymetry
             */
            inline void storeParameters(T h_l, T h_r, T hu_l, T hu_r, T b_l, T b_r) {
                this->h_l = h_l;
                this->h_r = h_r;
                this->hu_l = hu_l;
                this->hu_r = hu_r;
                
                this->b_l = b_l;
                this->b_r = b_r;
             
                if(h_l > tolerance)
                    this->u_l = hu_l / h_l;
                else
                    this->u_l = (T)0.0;
             
                if(h_r > tolerance)
                    this->u_r = hu_r / h_r;
                else
                    this->u_r = (T)0.0;
            }
         
           /**
            * Compute eigenvalues used in wave decomposition and store in instance variables
            */
            inline void computeEigenvalues() {
                T sqrt_h_l = std::sqrt(h_l);
                T sqrt_h_r = std::sqrt(h_r);
                T velocity = ( u_l * sqrt_h_l + u_r * sqrt_h_r ) / ( sqrt_h_l + sqrt_h_r );
                T phaseVelocity = std::sqrt( (T)0.5 * gravity * ( h_l + h_r ) );
                lambda_1 = velocity - phaseVelocity;
                lambda_2 = velocity + phaseVelocity;
            }
            
            /**
             * Compute the flux jump as a function of left and right water height and momentum.
             * Results are stored in instance variables delta_f_1 and delta_f_2
             */
            inline void computeFluxJump() {
               /**
                * **Compute flux jump**
                * 
                * \f{align*}{ \Delta f := f(q_r) - f(q_l) - \Delta x \Psi_{i-1/2}
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
                * -
                * \begin{bmatrix}
                * 0 \\
                * -g(b_r-b_l)\frac{h_l+h_r}{2}
                * \end{bmatrix}
                * \\
                * &=
                * \begin{bmatrix}
                * hu_r - hu_l \\
                * hu^2_r - hu^2_l + \frac{1}{2} g (h_r^2 - h_l^2) + \frac{1}{2} g(b_r-b_l)(h_l+h_r)
                * \end{bmatrix}
                * \\
                * &=
                * \begin{bmatrix}
                * hu_r - hu_l \\
                * hu^2_r - hu^2_l + \frac{1}{2} g (h_r (h_r + b_r - b_l) + h_l (-h_l + b_r - b_l))
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
                delta_f_2 = (hu_r * u_r) - (hu_l * u_l)
                    + (T)0.5 * gravity * (h_r * (h_r + b_r - b_l) + h_l * (-h_l + b_r - b_l));
            }
            
            /**
             * Compute the eigencoefficients needed to compute the resulting waves.
             * Results are stored in instance variables alpha_1 and alpha_2
             */
            inline void computeEigencoefficients() {
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
            FWave(T _gravity = (T)9.81, T _tolerance = (T)1e-10):
                gravity(_gravity), tolerance(_tolerance) {}
            
            /** Destructor */
            ~FWave() {}
            
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
                T waveSpeedLeft, waveSpeedRight;
                computeNetUpdates(h_l, h_r,
                    hu_l, hu_r,
                    b_l, b_r,
                    netUpdateLeft_h, netUpdateRight_h,
                    netUpdateLeft_hu, netUpdateRight_hu,
                    waveSpeedLeft, waveSpeedRight
                );
                
                /**
                 * **Compute Maximum wavespeeds**
                 *
                 * \f[\lambda_{\text{max}} = \max\left\{ \left| \lambda_1 \right|, \left| \lambda_2 \right| \right\}\f]
                 */
                // Note: We're not using the raw lambda_1/2 here, but the returned left and 
                // right-going wave speeds, since it doesn't make any difference if we take
                // the maximum abolute value of both
                maxWaveSpeed = std::max(std::fabs(waveSpeedLeft), std::fabs(waveSpeedRight));
            }
            /**
             * Compute net updates and left and right-going wave speeds
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
             * @param[out]  waveSpeedLeft       Left-going wave speed
             * @param[out]  waveSpeedRight      Right-going wave speed
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
                T &waveSpeedLeft,
                T &waveSpeedRight
            ) {
                // Height cannot be negative
                assert(h_l >= (T)0.0);
                assert(h_r >= (T)0.0);
                
                // init values
                netUpdateLeft_h = (T)0.0;
                netUpdateLeft_hu = (T)0.0;
                netUpdateRight_h = (T)0.0;
                netUpdateRight_hu = (T)0.0;
                waveSpeedLeft = (T)0.0;
                waveSpeedRight = (T)0.0;
                                
                // handle edge case "both heights numerically zero"
                // in that case, just return zero for everything
                if(h_l < tolerance && h_r < tolerance)
                    return;
                
                // Check for dry-wet / wet-dry cases
                // Handle bathymetry edge-cases (wet-dry / dry-wet)
                if(b_l <= -tolerance && b_r <= -tolerance) {
                    // left and right cell are both wet
                    storeParameters(h_l, h_r, hu_l, hu_r, b_l, b_r);
                } else if(b_l <= -tolerance) {
                    // right one is a dry cell: reflecting boundary
                    // right cell should be the same height and bathymetry
                    // but opposite momentum
                    storeParameters(h_l, h_l, hu_l, -hu_l, b_l, b_l);
                } else if(b_r <= -tolerance) {
                    // left one is a dry cell: reflecting boundary
                    // left cell should be the same height and bathymetry
                    // but opposite momentum
                    storeParameters(h_r, h_r, -hu_r, hu_r, b_r, b_r);
                } else {
                    // left and right cell are both dry
                    // nothing to do here
                    return;
                }
                
                computeEigenvalues();
                computeEigencoefficients();
                
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
                if(lambda_2 >= (T)0.0)
                    waveSpeedRight = lambda_2;
                else
                    waveSpeedRight = (T)0.0;
                
                if(lambda_1 <= (T)0.0)
                    waveSpeedLeft = lambda_1;
                else
                    waveSpeedLeft = (T)0.0;
                
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
                if(lambda_1 < (T)0.0) {
                    // to left
                    if(b_l < -tolerance) {
                        netUpdateLeft_h += alpha_1;
                        netUpdateLeft_hu += alpha_1 * lambda_1;
                    }
                } else {
                    // to right
                    if(b_r < -tolerance) {
                        netUpdateRight_h += alpha_1;
                        netUpdateRight_hu += alpha_1 * lambda_1;
                    }
                }
         
                if(lambda_2 >= (T)0.0) {
                    // to right
                    if(b_r < -tolerance) {
                        netUpdateRight_h += alpha_2;
                        netUpdateRight_hu += alpha_2 * lambda_2;
                    }
                } else {
                    // to left
                    if(b_l < -tolerance) {
                        netUpdateLeft_h += alpha_2;
                        netUpdateLeft_hu += alpha_2 * lambda_2;
                    }
                }
            }
    };
};

#endif
