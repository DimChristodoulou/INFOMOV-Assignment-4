#define PI 3.141592653589793238462643383279 

typedef struct cl_Ray
{
	float3 origin;
	float3 direction;
	float t;
} Ray;

typedef struct cl_material {
	float3 albedo;
	float fuzz;
	float ir;
	float3 padding0;
	int materialType;
	int3 padding1;	// 0 for lambertian, 1 for metal, 2 for dielectric
} Material;


typedef struct cl_hit_record {
	float3 p;
	float3 normal;
	float t;
	float padding0;
	Material mat;
	bool front_face;
	bool padding1;
	bool padding2;
	bool padding3;
} hit_record;

typedef struct cl_sphere
{
	float3 center;
	float radius;
	Material mat;
} Sphere;


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

void set_face_normal(hit_record* rec, Ray* ray, float3 outward_normal)
{
	rec->front_face = dot(ray->direction, outward_normal) < 0;
	rec->normal = rec->front_face ? outward_normal : -outward_normal;
}
float3 random_unit_vector(uint *seed) {
    return UnitVector(random_in_unit_sphere(seed));
}

float SquaredLength(float3 vec) {
	return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

bool near_zero(float3 vec) {
	// Return true if the vector is close to zero in all dimensions.
	const float s = 1e-8f;
	return (fabs(vec.x) < s) && (fabs(vec.y) < s) && (fabs(vec.z) < s);
}

bool scatter(Ray *ray, hit_record *hit, float3 *attenuation, Ray *scattered, uint *seed){
	int id = get_global_id(0);

	// Lambertian
	if(hit->mat.materialType == 0){
		float3 scatter_direction = hit->normal + random_unit_vector(seed);

		// Catch degenerate scatter direction
		if (near_zero(scatter_direction)){
			scatter_direction = hit->normal;
		}

		scattered->origin = hit->p;
		scattered->direction = scatter_direction;
		*attenuation = hit->mat.albedo;
		return true;
	}
	// Metal
	else if(hit->mat.materialType == 1){
	}
	// Dielectric
	else if(hit->mat.materialType == 2){
	}
}

bool sphere_hit(Sphere* sphere,
	Ray* ray, float t_min, float t_max, hit_record* rec)
{
	float3 center = sphere->center;
	float radius = sphere->radius;
	float3 oc = ray->origin - center;
	float a = SquaredLength(ray->direction);
	float half_b = dot(oc, ray->direction);
	float c = SquaredLength(oc) - radius * radius;

	float discriminant = half_b * half_b - a * c;
	if (discriminant < 0) return false;
	float sqrtd = sqrt(discriminant);


	float root = (-half_b - sqrtd) / a;
	if (root < t_min || t_max < root)
	{
		root = (-half_b + sqrtd) / a;
		if (root < t_min || t_max < root)
			return false;
	}

	rec->t = root;
	rec->p = ray->origin + root * ray->direction;
	float3 outward_normal = (rec->p - center) / radius;
	set_face_normal(rec, ray, outward_normal);
	rec->mat = sphere->mat;
	return true;


}


bool world_hit(Ray* ray, float t_min, float t_max, hit_record* rec, int nSpheres, Sphere* sphereBuffer)
{
	hit_record temp_rec;
	bool hit_anything = false;
	float closest_so_far = t_max;
	for (int i = 0; i < nSpheres; i++)
	{
		Sphere sphere = sphereBuffer[i];
		if (sphere_hit(&sphere, ray, t_min, t_max, rec))
		{
			hit_anything = true;
			closest_so_far = temp_rec.t;
			(*rec) = temp_rec;
		}
	}

	return hit_anything;
}


float4 ray_color(Ray* r, int nSpheres, Sphere* sphereBuffer, uint *seed)
{
	hit_record rec;
	if (world_hit(r, 0.001, 10000000000.0, &rec, nSpheres, sphereBuffer))
	{
		float4 attenuation;
		Ray scattered;
		if (scatter(r, &rec, &attenuation, &scattered, seed))
			return attenuation;
		return (float4)(0.0, 0.0, 0.0, 0.0);
	}

	float3 unit_direction = normalize(r->direction);
	float t = 0.5 * (unit_direction.y + 1.0);
	return (1.0 - t)* (float4)(1.0, 1.0, 1.0, 0.0) + t * (float4)(0.5, 0.7, 1.0, 0.0);
}

__kernel void raytrace(	__global float4* colorBuffer,
						__global Sphere* sphereBuffer,
						int image_width,
						int image_height,
						float3 vertical,
						float3 horizontal,
						float3 lower_left_corner,
						float3 origin,
						int nSpheres)
{
	uint id = get_global_id(0);/*get_global_id(0) / 32;*/
	int row = id % image_width; //This will depend on how the memory is laid out in the 2d array. 
    int col = id / image_width; //If it's not row-major, then we'll need to flip these two statements.

	float u =  (float) row  / ((float)image_width - 1.0);
	float v = (float)col / ((float)image_height - 1.0);

	//printf("col %d, row %d , colorBuffer %d\n", col, row, id);

	uint seed = Next(0x12345678 + id);
	// float3 dir = (float3)(0.0f, 0.0f, 0.0f);
	// if(get_local_id(0)==0){
	// 	float3 dir = random_in_unit_sphere(&seed);
	// 	//dir = (float3)(1, 1, 1);
	// 	printf(" %u, %f %f, %f\n",get_group_id(0), dir.x, dir.y, dir.z);
	// }

	Ray ray;
	ray.origin = origin;
	ray.direction = lower_left_corner + u * horizontal + v * vertical - origin;
	ray.t = 0.0001;

	//float3 unit_direction = normalize(ray.direction);
	//float t = 0.5 * (unit_direction.y + 1.0);
	float4 color = (0.7f, 0.65f, 0.9f, 0.0f);
	
	__local float4 _sharedMemory[32];
	_sharedMemory[get_local_id(0)] = ray_color(&ray, nSpheres, sphereBuffer, &seed);
	barrier(CLK_LOCAL_MEM_FENCE);
	
	if (get_local_id(0) == 0)
	{
		for (int i = 0; i < 32; i++)
			color += _sharedMemory[i];
		colorBuffer[id] = color;
	}
		 
	// printf("%f\n", color);
	
}