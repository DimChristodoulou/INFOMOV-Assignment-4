
struct Ray
{
	float3 origin;
	float3 direction;
	float t;
};




__kernel void raytrace(__global float3* colorBuffer,
						int image_width,
						int image_height,
						float3 vertical,
						float3 horizontal,
						float3 lower_left_corner,
						float3 origin)
{
	int j = image_height - 1 - get_global_id(1);
	int i = get_global_id(0);
	float u =  (float) i  / ((float)image_width - 1.0);
	float v = (float)j / ((float)image_height - 1.0);

	struct Ray ray;
	ray.origin = origin;
	ray.direction = lower_left_corner + u * horizontal + v * vertical - origin;
	ray.t = 0.0001;
	float3 unit_direction = normalize(ray.direction);
	float t = 0.5 * (unit_direction.y + 1.0);
	float3 color =  (1.0 - t) * (float3)(1.0, 1.0, 1.0) + t * (float3)(0.5, 0.7, 1.0);
	colorBuffer[j * image_width + i] = color;
}