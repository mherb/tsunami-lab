
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
    if(b.x > -tolerance && b.y > -tolerance) {
        // left and right cell are both dry
        // nothing to do here
        return;
	}else{

	// b.x <= -tolerance && b.y > -tolerance
	int rightCellDry = (b.x <= -tolerance)*(b.y > -tolerance);
        // right one is a dry cell: reflecting boundary
        // right cell should be the same height and bathymetry
        // but opposite momentum
		h.y = rightCellDry * h.x + !rightCellDry * h.y;
		hu.y = -rightCellDry * hu.x + !rightCellDry * hu.y;
		b.y = rightCellDry * b.x + !rightCellDry * b.y;
        boundary_type = rightCellDry * 2 + !rightCellDry * boundary_type;

	// b.y <= -tolerance && b.x > -tolerance
	int leftCellDry = (b.y <= -tolerance)*(b.x > -tolerance);
        // left one is a dry cell: reflecting boundary
        // left cell should be the same height and bathymetry
        // but opposite momentum
		h.x = leftCellDry * h.y + !leftCellDry * h.x;
		hu.x = - leftCellDry * hu.y + !leftCellDry * hu.x;
		b.x = leftCellDry * b.y + !leftCellDry * b.x;
        boundary_type = leftCellDry * 1 + !leftCellDry * boundary_type;
    }
    
    const float u_l = (h.x > tolerance) ? (hu.x / h.x) : 0.f;
    const float u_r = (h.y > tolerance) ? (hu.y / h.y) : 0.f;
	const float2 u = (float2)(u_l, u_r);
    

    float2 sqrt_h = (float2)(sqrt(h.x), sqrt(h.y));
    // compute velocity
    float velocity = (dot(u, sqrt_h)) / ( sqrt_h.x + sqrt_h.y );
	//compute phase velocity
    float phase_velocity = sqrt(0.5f * gravity * ( h.x + h.y ) );
    // Compute eigenvalues lambda
 	float2 lambda = (float2)(velocity - phase_velocity, velocity + phase_velocity);

	// Compute auxiliary vector delta_f_h to be able to use dot product when computing delta_f
	float2 delta_f_h = (float2)(-h.x + b.y - b.x, h.y + b.y - b.x);
    float2 delta_f = (float2)(hu.y - hu.x,
		(hu.y * u.y) - (hu.x * u.x)
        + 0.5f * gravity * (dot(h, delta_f_h)));

    float lambda_diff = lambda.y - lambda.x;

    // Compute the eigencoefficients alpha needed
    // to compute the resulting waves.
    float2 alpha = (float2)((lambda.y * delta_f.x - delta_f.y), (-lambda.x * delta_f.x + delta_f.y))/ lambda_diff;
    
    // Compute max wavespeed
    *max_wave_speed = fmax(fabs(lambda.x), fabs(lambda.y));

    // Compute net updates
        //if(lambda.x < 0.f && boundary_type != 1)
        net_update_h.x += (lambda.x < 0.f)*(boundary_type != 1)*(alpha.x);
        net_update_hu.x += (lambda.x < 0.f)*(boundary_type != 1)*(alpha.x * lambda.x);
        
		//if(lambda.x >= 0.f && boundary_type != 2)
		net_update_h.y += (lambda.x >= 0.f)*(boundary_type != 2)*(alpha.x);
        net_update_hu.y += (lambda.x >= 0.f)*(boundary_type != 2)*(alpha.x * lambda.x);
    
   		//if(lambda.y >= 0.f && boundary_type != 2)
        net_update_h.y += (lambda.y >= 0.f)*(boundary_type != 2)*(alpha.y);
        net_update_hu.y += (lambda.y >= 0.f)*(boundary_type != 2)*(alpha.y * lambda.y);

       	//if(lambda.y < 0.f && boundary_type != 1)
        net_update_h.x += (lambda.y < 0.f)*(boundary_type != 1)*(alpha.y);
        net_update_hu.x += (lambda.y < 0.f)*(boundary_type != 1)*(alpha.y * lambda.y);

	*net_update_h_l = net_update_h.x;
	*net_update_h_r = net_update_h.y;
	*net_update_hu_l = net_update_hu.x;
	*net_update_hu_r = net_update_hu.y;
}
