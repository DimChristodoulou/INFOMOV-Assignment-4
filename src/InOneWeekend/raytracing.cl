#define PI 3.141592653589793238462643383279 

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


 uint Next(uint u)
{
	uint v = u * 3935559000370003845 + 2691343689449507681;

	v ^= v >> 21;
	v ^= v << 37;
	v ^= v >> 4;

	v *= 4768777513237032717;

	v ^= v << 20;
	v ^= v >> 41;
	v ^= v << 5;

	return v;
}
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
     return *seed * 2.3283064365387e-10f;
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
	//printf("%f \n", random_float(seed));
	float theta = random_float_in_range(seed, 0, 2.0 * PI);
	//printf("%f \n", theta);
	float z = random_float_in_range(seed, -1.0, 1.0);
	//printf("%f \n", z);

	float z2 = z * z;
	float sqrOneMinusZ2 = sqrt(1 - z2);
	return (float3)(sqrOneMinusZ2 * cos(theta), sqrOneMinusZ2 * sin(theta), z);
    /*while (true) {
        float3 p = randomVec(seed, -1,1);
        if (length_squared(p) >= 1) continue;
        return p;
    }*/
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
	uint id = get_global_id(0);/*get_global_id(0) / 32;*/
	int row = id % image_width; //This will depend on how the memory is laid out in the 2d array. 
    int col = id / image_width; //If it's not row-major, then we'll need to flip these two statements.

	float u =  (float) row  / ((float)image_width - 1.0);
	float v = (float)col / ((float)image_height - 1.0);

	//printf("col %d, row %d , colorBuffer %d\n", col, row, id);

	uint seed = Next(0x12345678 + id);
	float3 dir = (float3)(0.0f, 0.0f, 0.0f);
	if(get_local_id(0)==0){
		float3 dir = random_in_unit_sphere(&seed);
		//dir = (float3)(1, 1, 1);
		printf(" %u, %f %f, %f\n",get_group_id(0), dir.x, dir.y, dir.z);
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