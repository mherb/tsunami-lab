// Use local or global memory according to preprocessor definition
#ifdef MEM_LOCAL
#define SOLVER_MEM __local
#else
#define SOLVER_MEM __global
#endif

// Gravity in m/s^2
#define GRAVITY 9.81f

// Numerical tolerance for comparisons
#define TOLERANCE 0.00001

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
    SOLVER_MEM float* net_update_h_l, SOLVER_MEM float* net_update_h_r,
    SOLVER_MEM float* net_update_hu_l, SOLVER_MEM float* net_update_hu_r,
    SOLVER_MEM float* max_wave_speed) {
    
    float2 h  = (float2)( h_l,  h_r  );
    float2 hu = (float2)( hu_l, hu_r );
    float2 b  = (float2)( b_l,  b_r  );
    
    // Hint: maybe use native_* functions (Table 6.9) for roots and division?
    
    // TODO: init the buffers with zero using clEnqueueFillBuffer
    // so we don't have to do that in the kernel
    *net_update_h_l = 0.f;
    *net_update_h_r = 0.f;
    *net_update_hu_l = 0.f;
    *net_update_hu_r = 0.f;
    *max_wave_speed = 0.f;
    
    // handle edge case "both heights numerically zero"
    // in that case, just return zero for everything
    if(h.x < TOLERANCE && h.y < TOLERANCE)
        return;
    
    // Boundary type: Regular (0), Reflecting Left Boundary (1), Reflecting Right Boundary (2)
    // TODO: maybe use typedef or #defines here
    unsigned char boundary_type = 0;
    
    // Check for dry-wet / wet-dry cases
    // Handle bathymetry edge-cases (wet-dry / dry-wet)
    if(b.x <= -TOLERANCE && b.y <= -TOLERANCE) {
        // left and right cell are both wet
        // no need to change values
    } else if(b.x <= -TOLERANCE) {
        // right one is a dry cell: reflecting boundary
        // right cell should be the same height and bathymetry
        // but opposite momentum
        h.y = h.x;
        hu.y = -hu.x;
        b.y = b.x;
        boundary_type = 2;
    } else if(b.y <= -TOLERANCE) {
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
    
    const float u_l = (h.x > TOLERANCE) ? (hu.x / h.x) : 0.f;
    const float u_r = (h.y > TOLERANCE) ? (hu.y / h.y) : 0.f;
    float2 u = (float2)(u_l, u_r);
    
    
    float2 sqrt_h = sqrt(h);
    // compute velocity
    float velocity = dot(u, sqrt_h) / ( sqrt_h.x + sqrt_h.y );
    //compute phase velocity
    float phase_velocity = sqrt(0.5f * GRAVITY * ( h.x + h.y ) );
    // Compute eigenvalues lambda
    float2 lambda = (float2)(velocity - phase_velocity, velocity + phase_velocity);
    
    // Compute flux jump
    u.x = -u.x; // negate first component of u to use dot product
    float2 delta_f = (float2)(hu.y - hu.x,
        dot(u, hu) + 0.5f * GRAVITY * dot(h, (float2)(-h.x + b.y - b.x, h.y + b.y - b.x)));
    
    // Compute the eigencoefficients alpha needed
    // to compute the resulting waves.
    float2 alpha = (float2)(dot((float2)(lambda.y, -1.f), delta_f), dot((float2)(-lambda.x, 1.f), delta_f));
    alpha /= (lambda.y - lambda.x);
    
    // Compute max wavespeed
    *max_wave_speed = fmax(fabs(lambda.x), fabs(lambda.y));
    
    float2 net_update_left =  (float2)( 0.f, 0.f );
    float2 net_update_right = (float2)( 0.f, 0.f );
    // Compute net updates
    if(lambda.x < 0.f) {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1)
            net_update_left += (float2)(alpha.x, alpha.x * lambda.x);
    } else {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2)
            net_update_right += (float2)(alpha.x, alpha.x * lambda.x);
    }
    if(lambda.y >= 0.f) {
        // to right: if right bound is not reflecting => update
        if(boundary_type != 2)
            net_update_right += (float2)(alpha.y, alpha.y * lambda.y);
    } else {
        // to left: if left bound is not reflecting => update
        if(boundary_type != 1)
            net_update_left += (float2)(alpha.y, alpha.y * lambda.y);
    }
    
    *net_update_h_l = net_update_left.x;
    *net_update_h_r = net_update_right.x;
    *net_update_hu_l = net_update_left.y;
    *net_update_hu_r = net_update_right.y;
}
