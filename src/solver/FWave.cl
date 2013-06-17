
// Gravity in m/s^2
__constant float gravity = 9.81f;

// Numerical tolerance for comparisons
__constant float tolerance = 1e-10;

/**
 * OpenCl FWave solver Kernel function
 * 
 * Compute net updates and maximum wavespeed for a single edge
 *
 * @param[in]   h                   Left and right water height
 * @param[in]   hu                  Left and right momentum
 * @param[in]   b                   Left and right bathymetry
 * @param[out]  net_update_h        Left and right net update (water height)
 * @param[out]  net_update_hu       Left and right net update (momentum)
 * @param[out]  max_wave_speed      Maximum wave speed
 */
__kernel void computeNetUpdates(
    float2 h,
    float2 hu,
    float2 b,
    __global float2* net_update_h,
    __global float2* net_update_hu,
    __global float* max_wave_speed) {
    
    // Hint: maybe use native_* functions (Table 6.9) for roots and division?
    
    // TODO: init the buffers with zero using clEnqueueFillBuffer
    // so we don't have to do that in the kernel
    *net_update_h = (float2)(0.f, 0.f);
    *net_update_hu = (float2)(0.f, 0.f);
    *max_wave_speed = 0.f;
    
    // handle edge case "both heights numerically zero"
    // in that case, just return zero for everything
    if(h.x < tolerance && h.y < tolerance)
        return;
    
    // Boundary type: Regular (0), Reflecting Left Boundary (1), Reflecting Right Boundary (2)
    // TODO: maybe use typedef or #defines here
    unsigned char boundary_type = 0;
    
    // Check for dry-wet / wet-dry cases
    // Handle bathymetry edge-cases (wet-dry / dry-wet)
    if(b.x <= -tolerance && b.y <= -tolerance) {
        // left and right cell are both wet
		// no need to change values
    } else if(b.x <= -tolerance) {
        // right one is a dry cell: reflecting boundary
        // right cell should be the same height and bathymetry
        // but opposite momentum
		h.y = h.x;
		hu.y = -hu.x;
		b.y = b.x;
        boundary_type = 2;
    } else if(b.y <= -tolerance) {
        // left one is a dry cell: reflecting boundary
        // left cell should be the same height and bathymetry
        // but opposite momentum
		h.x = h.y;
		hu.x = -hu.y;
		b.x = b.y;
        boundary_type = 1;
    } else {
        // left and right cell are both dry
        // nothing to do here
        return;
    }
    
    const float u_l = (h.x > tolerance) ? (hu.x / h.x) : 0.f;
    const float u_r = (h.y > tolerance) ? (hu.y / h.y) : 0.f;
	const float2 u = (float2)(u_l, u_r);
    
    // Compute eigenvalues lambda.x and lambda.y
    float sqrt_h = (float2)(sqrt(h.x), sqrt(h.y));
    // Hint: use dot product here
    float velocity = ( u_l * sqrt_h.x + u_r * sqrt_h.y ) / ( sqrt_h.x + sqrt_h.y );
    // Hint: maybe use dot product here too?
    float phase_velocity = sqrt(0.5f * gravity * ( h.x + h.y ) );
    float2 lambda = (float2)(velocity - phase_velocity, velocity + phase_velocity);

    // Compute the eigencoefficients alpha.x and alpha.y needed
    // to compute the resulting waves.
    float2 delta_f = (float2)(hu.y - hu.x,
		(hu.y * u_r) - (hu.x * u_l)
        + 0.5f * gravity * (h.y * (h.y + b.y - b.x) + h.x * (-h.x + b.y - b.x)));

    float lambda_diff = lambda.y - lambda.x;

    // Hint: use dot product(s) here
    float2 alpha = (float2)((lambda.y * delta_f.x - delta_f.y) / lambda_diff, (-lambda.x * delta_f.x + delta_f.y) / lambda_diff);
    
    // Compute max wavespeed
    *max_wave_speed = fmax(fabs(lambda.x), fabs(lambda.y));
    
    // Compute net updates
    // Hint: use dot products here
    if(lambda.x < 0.f) {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1) {
            *net_update_h.x += alpha.x;
            *net_update_hu.x += alpha.x * lambda.x;
        }
    } else {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2) {
            *net_update_h.y += alpha.x;
            *net_update_hu.y += alpha.x * lambda.x;
        }
    }

    if(lambda.y >= 0.f) {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2) {
            *net_update_h.y += alpha.y;
            *net_update_hu.y += alpha.y * lambda.y;
        }
    } else {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1) {
            *net_update_h.x += alpha.y;
            *net_update_hu.x += alpha.y * lambda.y;
        }
    }
}
