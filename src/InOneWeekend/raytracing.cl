
typedef struct cl_Ray
{
	float3 origin;
	float3 direction;
	float t;
} Ray;

typedef struct cl_hit_record {
	float3 point;
	float3 normal;
} hit_record;

typedef struct cl_material {
	float3 albedo;
	double fuzz;
	double ir;
	int materialType;	// 0 for lambertian, 1 for metal, 2 for dielectric
} Material;

// RNG - Marsaglia's xor32
uint RandomUInt(uint *seed)
{
	*seed ^= *seed << 13;
	*seed ^= *seed >> 17;
	*seed ^= *seed << 5;
	return *seed;
}

float3 UnitVector(float3 v) {
	return v / length(v);
}

float length_squared(float3 vec) {
    return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}

float random_float(uint *seed) {
    // Returns a random real in [0,1).
    return RandomUInt(seed) / (0x7fff + 1.0);
}

float random_float_in_range(uint *seed, float min, float max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_float(seed);
}

// Returns a vec between min and max
float3 randomVec(uint *seed, float min, float max){
	return (float3)(random_float_in_range(seed,min,max), random_float_in_range(seed,min,max), random_float_in_range(seed,min,max));
}

float3 random_in_unit_sphere(uint *seed) {
    while (true) {
        float3 p = randomVec(seed, -1,1);
        if (length_squared(p) >= 1) continue;
        return p;
    }
}

float SquaredLength(float3 vec) {
	return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

bool scatter(Material *mat, Ray ray, hit_record *hit, float3 *attenuation, Ray *scattered){
	int id = get_global_id(0);

	// Lambertian
	if(mat->materialType == 0){
		
	}
	// Metal
	else if(mat->materialType == 1){
	}
	// Dielectric
	else if(mat->materialType == 2){
	}
}

float4 ray_color(Ray* r)
{
	hit_record rec;
	if (hit(0.001, 10000000000.0, rec))
	{

	}
}

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

	uint seed = 0x12345678 + id;
	if(get_global_id(0)==0){
		printf("%u %u\n", RandomUInt(&seed), RandomUInt(&seed));
	}

	// Ray ray;
	// ray.origin = origin;
	// ray.direction = lower_left_corner + u * horizontal + v * vertical - origin;
	// ray.t = 0.0001;

	// float3 unit_direction = normalize(ray.direction);
	// float t = 0.5 * (unit_direction.y + 1.0);
	// __local float4 _sharedMemory[32];
	// _sharedMemory[get_local_id(0)] =  (1.0 - t) * (float4)(1.0, 1.0, 1.0, 0.0) + t * (float4)(0.5, 0.7, 1.0, 0.0);
	// barrier(CLK_LOCAL_MEM_FENCE);
	// float4 color = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
	// if (get_local_id(0) == 0)
	// {
	// 	for (int i = 0; i < 32; i++)
	// 		color += _sharedMemory[i];
	// 	colorBuffer[id] = color;
	// }
		 
	// printf("%f\n", color);
	
}