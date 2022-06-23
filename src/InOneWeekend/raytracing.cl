#define PI 3.141592653589793238462643383279 
#define NSAMPLES 32

typedef struct cl_Ray
{
	float3 origin;
	float3 direction;
	float t;
} Ray;

typedef struct __attribute__ ((__packed__)) cl_material {
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

typedef struct __attribute__ ((__packed__)) cl_sphere
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

void set_face_normal(hit_record* rec, Ray* ray, float3 outward_normal) {
	rec->front_face = (dot(ray->direction, outward_normal) < 0);

	if(rec->front_face){
		rec->normal = outward_normal;
	}
	else {
		rec->normal = -outward_normal;
	}
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

float3 reflect(float3 v, float3 n) {
    return v - 2*dot(v, n)*n;
}

float reflectance(float cosine, float ref_idx) {
	// Use Schlick's approximation for reflectance.
	float r0 = (1-ref_idx) / (1+ref_idx);
	r0 = r0*r0;
	return r0 + (1-r0)*pow((1 - cosine),5);
}

float3 refract(float3 uv, float3 n, float etai_over_etat) {
    float cos_theta = fmin(dot(-uv, n), 1.0);
    float3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    float3 r_out_parallel = -sqrt(fabs(1.0 - SquaredLength(r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
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
		scattered->direction = normalize(scatter_direction);
		scattered->t = 0.0001;
		*attenuation = hit->mat.albedo;
		return true;
	}
	// Metal
	else if(hit->mat.materialType == 1){
		// vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
		// scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
		// attenuation = albedo;
		// return (dot(scattered.direction(), rec.normal) > 0);
		float3 reflected = reflect(UnitVector(ray->direction), hit->normal);
		scattered->origin = hit->p;
		scattered->direction = normalize(reflected + (hit->mat.fuzz*random_in_unit_sphere(seed)));
		scattered->t = 0.0001;
		*attenuation = hit->mat.albedo;
		return (dot(scattered->direction, hit->normal) > 0);
	}
	// Dielectric
	else if(hit->mat.materialType == 2){
		*attenuation = (float3)(1.0, 1.0, 1.0);
		float refraction_ratio = hit->front_face ? (1.0/hit->mat.ir) : hit->mat.ir;
		
		float3 unit_direction = normalize(ray->direction);
		float cos_theta = fmin(dot(-unit_direction, hit->normal), 1.0);
		float sin_theta = sqrt(1.0 - cos_theta*cos_theta);

		bool cannot_refract = refraction_ratio * sin_theta > 1.0;

		if(cannot_refract || reflectance(cos_theta, refraction_ratio) > random_float(seed)){
			scattered->direction = normalize(reflect(unit_direction, hit->normal));
		}
		else{
			scattered->direction = normalize(refract(unit_direction, hit->normal, refraction_ratio));
		}

		scattered->t = 0.0001;
		scattered->origin = hit->p;
		return true;
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
	if (root < t_min || t_max < root) {
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
		if (sphere_hit(&sphere, ray, t_min, closest_so_far, &temp_rec))
		{
			hit_anything = true;
			closest_so_far = temp_rec.t;
			(*rec) = temp_rec;
		}
	}

	return hit_anything;
}


float4 ray_color(Ray* r, int nSpheres, Sphere* sphereBuffer, uint *seed, bool* hit_skybox)
{
	hit_record rec;
	if (world_hit(r, 0.0001, 10000000000.0, &rec, nSpheres, sphereBuffer))
	{
		float4 attenuation;
		Ray scattered;
		if (scatter(r, &rec, &attenuation, &scattered, seed))
		{
			(*r) = scattered;
			return attenuation;

		}

		return (float4)(0.0, 0.0, 0.0, 0.0);
	}

	float3 unit_direction = normalize(r->direction);
	float t = 0.5 * (unit_direction.y + 1.0);
	*hit_skybox = true;
	return (1.0 - t)* (float4)(1.0, 1.0, 1.0, 0.0) + t * (float4)(1.0, 0.7, 0.5, 0.0);
}

__kernel void raytrace(	__global float4* colorBuffer,
						__global Sphere* sphereBuffer,
						int image_width,
						int image_height,
						float3 vertical,
						float3 horizontal,
						float3 lower_left_corner,
						float3 origin,
						int maxDepth,
						int nSpheres,
						int nSamples)
{
	uint id = get_group_id(0);/*get_global_id(0) / 32;*/
	int row = id % image_width; //This will depend on how the memory is laid out in the 2d array. 
    int col = id / image_width; //If it's not row-major, then we'll need to flip these two statements.

	float u =  (float) row  / ((float)image_width - 1.0);
	float v = (float)col / ((float)image_height - 1.0);

	//printf("col %d, row %d , colorBuffer %d\n", col, row, id);

	uint seed = Next(0x12345678 + id);
	// float3 dir = (float3)(0.0f, 0.0f, 0.0f);
	if(get_global_id(0)==0){
		// for(size_t i =0; i< nSpheres; i++){
		// 	if(sphereBuffer[i].mat.materialType == 2){
		// 		printf("IR: %f \n", sphereBuffer[i].mat.ir);
		// 	}
		// }
	// 	printf("%d %f %f %f %f %f %f %f %d\n", sizeof(Sphere), s.center.x, s.center.y, s.center.z, s.radius, s.mat.albedo.x, s.mat.albedo.y, s.mat.albedo.z, s.mat.materialType);
	//	printf("%d\n", nSpheres);
	}

	Ray ray;
	ray.origin = origin;
	ray.direction = lower_left_corner + u * horizontal + v * vertical - origin;
	ray.t = 0.0001;

	//float3 unit_direction = normalize(ray.direction);
	//float t = 0.5 * (unit_direction.y + 1.0);
	//float4 color = (0.7f, 0.65f, 0.9f, 0.0f);
	
	__local float4 _sharedMemory[NSAMPLES];
	float4 rtCol = (float4)(1.0f, 1.0f, 1.0f, 1.0f);
	bool hit_skybox = false;
	for (int i = 0; i < maxDepth; i++)
	{
		if (!hit_skybox)
			rtCol *= ray_color(&ray, nSpheres, sphereBuffer, &seed, &hit_skybox);
		else if	(!hit_skybox && (i == maxDepth - 1))
			rtCol = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
		/*else
			rtCol = rtCol;*/

	}
	_sharedMemory[get_local_id(0)] = rtCol;
	barrier(CLK_LOCAL_MEM_FENCE);
	
	if (get_group_id(0) == 0)
	{
		//printf("%f\n", _sharedMemory[get_local_id(0)]);
	}

	if (get_local_id(0) == 0)
	{
		float4 color = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
		for (int i = 0; i < NSAMPLES; i++)
			color += _sharedMemory[i];
		colorBuffer[id] = color;
		//if (get_group_id(0) == 0)
			//printf("%f\n", colorBuffer[id]);
	}
	//if (get_global_id(0) == 0)
	//{
	//	printf("%d\n", nSpheres);
	//	for (int i = 0; i < nSpheres; i++)
	//	{
	//		//printf("%f\n", sphereBuffer[i].center);
	//	}
	//}

	// printf("%f\n", color);
	
}