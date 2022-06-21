
struct Ray
{
	float3 origin;
	float3 direction;
	float t;
};




__kernel void raytrace(__global float4* colorBuffer,
						int image_width,
						int image_height,
						float3 vertical,
						float3 horizontal,
						float3 lower_left_corner,
						float3 origin)
{
	int id = get_global_id(0) / 32;
	int row = id % image_width; //This will depend on how the memory is laid out in the 2d array. 
    int col = id / image_width; //If it's not row-major, then we'll need to flip these two statements.

	float u =  (float) row  / ((float)image_width - 1.0);
	float v = (float)col / ((float)image_height - 1.0);

	//printf("col %d, row %d , colorBuffer %d\n", col, row, id);

	struct Ray ray;
	ray.origin = origin;
	ray.direction = lower_left_corner + u * horizontal + v * vertical - origin;
	ray.t = 0.0001;

	float3 unit_direction = normalize(ray.direction);
	float t = 0.5 * (unit_direction.y + 1.0);
	__local float4 _sharedMemory[32];
	_sharedMemory[get_local_id(0)] =  (1.0 - t) * (float4)(1.0, 1.0, 1.0, 0.0) + t * (float4)(0.5, 0.7, 1.0, 0.0);
	barrier(CLK_LOCAL_MEM_FENCE);
	float4 color = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
	if (get_local_id(0) == 0)
	{
		for (int i = 0; i < 32; i++)
			color += _sharedMemory[i];
		colorBuffer[id] = color;
	}
		 
	// printf("%f\n", color);
	
}