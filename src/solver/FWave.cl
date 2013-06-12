#include <math.h> 


            /**
             * Gravity in m/s^2
             */
            __constant float gravity = 9.81;
         
            /**
             * Numerical tolerance for comparisons
             */
            __constant float tolerance = 1e-10;


            /**
             * OpenCl Kernel function
			 * 
			 * Compute net updates and wavespeed of each cell in one row
             *
             * @param[in]   h                 water height
             * @param[in]   hu                momentum
             * @param[in]   b                 bathymetry
             * @param[in]   sizex             grid size in x-dimension
             * @param[in]   sizey             grid size in y-dimension
             * @param[out]  netUpdate_h       net updates (water height)
             * @param[out]  netUpdate_hu      net updates (momentum)
             * @param[out]  WaveSpeed         wave speeds
             */
		__kernel void computeNetUpdates_opencl(
			__global float* h,
            __global float* hu,
            __global float* b,
			__global size_t sizex,
			__global size_t sizey,
            __global float* netUpdate_h,
            __global float* netUpdate_hu,
            __global float* WaveSpeed
            ) {

			//Thread IDs in x and y dimension
			size_t xId = get_global_id(0);
			size_t yId = get_global_id(1);

			//make sure we stay inside the array indices
			if (xId < sizex-1 && yId < sizey){
				//left water height
				float h_l=h[xId][yId];
				//right water height
				float h_r=h[xId+1][yId];
				//left momentum
				float hu_l=hu[xId][yId];
				//right momentum
				float hu_r=hu[xId+1][yId];
				//left bathymetry
				float b_l=b[xId][yId];
				//right bathymetry
				float b_r=b[xId+1][yId];

                // init values
                netUpdate_h[xId][yId] = 0.0;
                netUpdate_h[xId+1][yId] = 0.0;
				netUpdate_hu[xId][yId] = 0.0;
                netUpdate_hu[xId+1][yId] = 0.0;
                WaveSpeed[xId][yId] = 0.0;
                WaveSpeed[xId+1][yId] = 0.0;
                                
                // handle edge case "both heights numerically zero"
                // in that case, just return zero for everything
                if(h_l < tolerance && h_r < tolerance)
                    return;
                
                // Check for dry-wet / wet-dry cases
                // Handle bathymetry edge-cases (wet-dry / dry-wet)
                if(b_l <= -tolerance && b_r <= -tolerance) {
                    // left and right cell are both wet
					//no need to change values
                } else if(b_l <= -tolerance) {
                    // right one is a dry cell: reflecting boundary
                    // right cell should be the same height and bathymetry
                    // but opposite momentum
					h_r=h_l;
					hu_r=-hu_l;
					b_r=b_l;
                } else if(b_r <= -tolerance) {
                    // left one is a dry cell: reflecting boundary
                    // left cell should be the same height and bathymetry
                    // but opposite momentum
					h_l=h_r;
					hu_l=-hu_r;
					b_l=b_r;
                } else {
                    // left and right cell are both dry
                    // nothing to do here
                    return;
                }

        	    /**
				 * Compute eigenvalues lambda_1 and lambda_2
         		 */
                float sqrt_h_l = sqrt(h_l);
                float sqrt_h_r = sqrt(h_r);
                float velocity = ( u_l * sqrt_h_l + u_r * sqrt_h_r ) / ( sqrt_h_l + sqrt_h_r );
                float phaseVelocity = sqrt(0.5 * gravity * ( h_l + h_r ) );
                float lambda_1 = velocity - phaseVelocity;
                float lambda_2 = velocity + phaseVelocity;
    

                /**
         	     * Compute the eigencoefficients alpha_1 and alpha_2 needed to compute the resulting waves.
        	     */
                float delta_f_1 = hu_r - hu_l;
                float delta_f_2 = (hu_r * u_r) - (hu_l * u_l)
                    + 0.5 * gravity * (h_r * (h_r + b_r - b_l) + h_l * (-h_l + b_r - b_l));

                float lambda_diff = lambda_2 - lambda_1;
             
                float alpha_1 = (lambda_2 * delta_f_1 - delta_f_2) / lambda_diff;
                float alpha_2 = (-lambda_1 * delta_f_1 + delta_f_2) / lambda_diff;
            
                
                /**
                 * **Compute Wavespeeds**
                 */
				//WaveSpeed right
                if(lambda_2 >= 0.0)
                    WaveSpeed[xId+1][yId] = lambda_2;
                else
                    WaveSpeed[xId+1][yId] = 0.0;

				//WaveSpeed left
                if(lambda_1 <= 0.0)
                    WaveSpeed[xId][yId] = lambda_1;
                else
                    WaveSpeed[xId][yId] = 0.0;
                
                /**
                 * **Compute net updates**
                 */
                if(lambda_1 < 0.0) {
                    // to left
                    if(b_l < -tolerance) {
                        netUpdate_h[xId][yId] += alpha_1;
                        netUpdate_hu[xId][yId] += alpha_1 * lambda_1;
                    }
                } else {
                    // to right
                    if(b_r < -tolerance) {
                        netUpdate_h[xId+1][yId] += alpha_1;
                        netUpdate_hu[xId+1][yId] += alpha_1 * lambda_1;
                    }
                }
         
                if(lambda_2 >= 0.0) {
                    // to right
                    if(b_r < -tolerance) {
                        netUpdate_h[xId+1][yId] += alpha_2;
                        netUpdate_hu[xId+1][yId] += alpha_2 * lambda_2;
                    }
                } else {
                    // to left
                    if(b_l < -tolerance) {
                        netUpdate_h[xId][yId] += alpha_2;
                        netUpdate_hu[xId][yId] += alpha_2 * lambda_2;
                    }
                }
		}
	};

