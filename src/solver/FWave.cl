
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


	float2 h = (float2)(h_l, h_r);
    float2 hu = (float2)(hu_l, hu_r);
    float2 b= (float2)(b_l, b_r);

    
    // Hint: maybe use native_* functions (Table 6.9) for roots and division?
    
    // TODO: init the buffers with zero using clEnqueueFillBuffer
    // so we don't have to do that in the kernel
    float2 net_update_h = (float2)(0.f, 0.f);
    float2 net_update_hu = (float2)(0.f, 0.f);
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
    

    float2 sqrt_h = (float2)(sqrt(h.x), sqrt(h.y));
    // compute velocity
    float velocity = (dot(u, sqrt_h)) / ( sqrt_h.x + sqrt_h.y );
    // Hint: maybe use dot product here too?
    float phase_velocity = sqrt(0.5f * gravity * ( h.x + h.y ) );
    // Compute eigenvalues lambda
    float2 lambda_1 = (float2)(velocity - phase_velocity, -1); // (lambda.x, -1)
	float2 lambda_2 = (float2)(velocity + phase_velocity, -1); // (lambda.y, -1)

	// Compute auxiliary vector delta_f_h to be able to use dot product when computing delta_f
	float2 delta_f_h = (float2)(-h.x + b.y - b.x, h.y + b.y - b.x);
    float2 delta_f = (float2)(hu.y - hu.x,
		(hu.y * u.y) - (hu.x * u.x)
        + 0.5f * gravity * (dot(h, delta_f_h)));

    float lambda_diff = lambda_2.x - lambda_1.x;

    // Compute the eigencoefficients alpha needed
    // to compute the resulting waves.
    float2 alpha = (float2)(dot(lambda_1, delta_f), -dot(lambda_2, delta_f)) / lambda_diff;
    
    // Compute max wavespeed
    *max_wave_speed = fmax(fabs(lambda_1.x), fabs(lambda_2.x));

    // Compute net updates
    if(lambda.x < 0.f) {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1) {
            net_update_h.x += alpha.x;
            net_update_hu.x += alpha.x * lambda_1.x;
        }
    } else {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2) {
            net_update_h.y += alpha.x;
            net_update_hu.y += alpha.x * lambda_1.x;
        }
    }
    if(lambda.y >= 0.f) {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2) {
            net_update_h.y += alpha.y;
            net_update_hu.y += alpha.y * lambda_2.x;
        }
    } else {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1) {
            net_update_h.x += alpha.y;
            net_update_hu.x += alpha.y * lambda_2.x;
        }
    }

	*net_update_h_l = net_update_h.x;
	*net_update_h_r = net_update_h.y;
	*net_update_hu_l = net_update_hu.x;
	*net_update_hu_r = net_update_hu.y;
}
