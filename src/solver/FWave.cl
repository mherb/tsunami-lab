
// Gravity in m/s^2
__constant float gravity = 9.81f;

// Numerical tolerance for comparisons
__constant float tolerance = 1e-10;

/**
 * OpenCl FWave solver Kernel function
 * 
 * Compute net updates and maximum wavespeed for a single edge
 *
 * @param[in]   h_l                 Left-side water height
 * @param[in]   h_r                 Right-side water height
 * @param[in]   hu_l                Left-side momentum
 * @param[in]   hu_r                Right-side momentum
 * @param[in]   b_l                 Left-side bathymetry
 * @param[in]   b_r                 Right-side bathymetry
 * @param[out]  net_update_h_l      Left net update (water height)
 * @param[out]  net_update_h_r      Right net update (water height)
 * @param[out]  net_update_hu_l     Left net update (momentum)
 * @param[out]  net_update_hu_r     Right net update (momentum)
 * @param[out]  max_wave_speed      Maximum wave speed
 */
__kernel void computeNetUpdates(
    float h_l, float h_r,
    float hu_l, float hu_r,
    float b_l, float b_r,
    __global float* net_update_h_l, __global float* net_update_h_r,
    __global float* net_update_hu_l, __global float* net_update_hu_r,
    __global float* max_wave_speed) {
    
    // Hint: maybe use native_* functions (Table 6.9) for roots and division?
    
    const float u_l = (h_l > tolerance) ? (hu_l / h_l) : 0.f;
    const float u_r = (h_r > tolerance) ? (hu_r / h_r) : 0.f;
    
    // TODO: init the buffers with zero using clEnqueueFillBuffer
    // so we don't have to do that in the kernel
    *net_update_h_l = 0.f;
    *net_update_h_r = 0.f;
    *net_update_hu_l = 0.f;
    *net_update_hu_r = 0.f;
    *max_wave_speed = 0.f;
    
    // handle edge case "both heights numerically zero"
    // in that case, just return zero for everything
    if(h_l < tolerance && h_r < tolerance)
        return;
    
    // Boundary type: Regular (0), Reflecting Left Boundary (1), Reflecting Right Boundary (2)
    // TODO: maybe use typedef or #defines here
    unsigned char boundary_type = 0;
    
    // Check for dry-wet / wet-dry cases
    // Handle bathymetry edge-cases (wet-dry / dry-wet)
    if(b_l <= -tolerance && b_r <= -tolerance) {
        // left and right cell are both wet
		// no need to change values
    } else if(b_l <= -tolerance) {
        // right one is a dry cell: reflecting boundary
        // right cell should be the same height and bathymetry
        // but opposite momentum
		h_r = h_l;
		hu_r = -hu_l;
		b_r = b_l;
        boundary_type = 2;
    } else if(b_r <= -tolerance) {
        // left one is a dry cell: reflecting boundary
        // left cell should be the same height and bathymetry
        // but opposite momentum
		h_l = h_r;
		hu_l = -hu_r;
		b_l = b_r;
        boundary_type = 1;
    } else {
        // left and right cell are both dry
        // nothing to do here
        return;
    }
    
    // Compute eigenvalues lambda_1 and lambda_2
    float sqrt_h_l = sqrt(h_l);
    float sqrt_h_r = sqrt(h_r);
    // Hint: use dot product here
    float velocity = ( u_l * sqrt_h_l + u_r * sqrt_h_r ) / ( sqrt_h_l + sqrt_h_r );
    // Hint: maybe use dot product here too?
    float phase_velocity = sqrt(0.5f * gravity * ( h_l + h_r ) );
    float lambda_1 = velocity - phase_velocity;
    float lambda_2 = velocity + phase_velocity;

    // Compute the eigencoefficients alpha_1 and alpha_2 needed
    // to compute the resulting waves.
    float delta_f_1 = hu_r - hu_l;
    // Hint: use dot product(s) here
    float delta_f_2 = (hu_r * u_r) - (hu_l * u_l)
        + 0.5f * gravity * (h_r * (h_r + b_r - b_l) + h_l * (-h_l + b_r - b_l));

    float lambda_diff = lambda_2 - lambda_1;

    // Hint: use dot product(s) here
    float alpha_1 = (lambda_2 * delta_f_1 - delta_f_2) / lambda_diff;
    float alpha_2 = (-lambda_1 * delta_f_1 + delta_f_2) / lambda_diff;
    
    // Compute max wavespeed
    *max_wave_speed = fmax(fabs(lambda_1), fabs(lambda_2));
    
    // Compute net updates
    // Hint: use dot products here
    if(lambda_1 < 0.f) {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1) {
            *net_update_h_l += alpha_1;
            *net_update_hu_l += alpha_1 * lambda_1;
        }
    } else {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2) {
            *net_update_h_r += alpha_1;
            *net_update_hu_r += alpha_1 * lambda_1;
        }
    }

    if(lambda_2 >= 0.f) {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2) {
            *net_update_h_r += alpha_2;
            *net_update_hu_r += alpha_2 * lambda_2;
        }
    } else {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1) {
            *net_update_h_l += alpha_2;
            *net_update_hu_l += alpha_2 * lambda_2;
        }
    }
}
