//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "precomp.h"
#include "rtweekend.h"

#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include <chrono>
#include <iostream>


typedef struct cl_material {
    float3 albedo;
    float fuzz_ir;
    float padding1, padding2;   // TODO: MAYBE REMOVE THE EXTRA PADDING
    uint materialType; // 0 for lambertian, 1 for metal, 2 for dielectric
} Material;

typedef struct cl_sphere
{
    float3 center;
    float radius;
    Material mat;
} CL_Sphere;

color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }
    
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(1.0, 0.7, 0.5);
}


hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    return world;
}
TheApp* CreateApp() { return new MyApp(); }


int main(){
    
    Kernel::InitCL();
    Kernel kernel("src/InOneWeekend/raytracing.cl", "raytrace");

    // Image
    const float aspect_ratio = 16.0 / 9.0;
    const int image_width = 1080;
    const int image_height = 607;
    const int samples_per_pixel = 32;
    const int max_depth = 50;
    static int nPixels = image_width * image_height;

    // Camera
    point3 lookfrom(13, 2, 3);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    vec3 w = unit_vector(lookfrom - lookat);
    vec3 u = unit_vector(cross(vup, w));
    vec3 v = cross(w, u);

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    
    // World
    auto world = random_scene();

    CL_Sphere* spheres = new CL_Sphere[world.objects.size()];

    for (size_t i = 0; i < world.objects.size(); i++)
    {
        sphere *obj = dynamic_cast<sphere*>(world.objects[i].get());
        spheres[i].center = float3(obj->center.x(), obj->center.y(), obj->center.z());
        spheres[i].radius = obj->radius;
        
        material* generic_mat = dynamic_cast<material*>(obj->mat_ptr.get());
        
        if (generic_mat->get_mat_type() == 0) {
            lambertian* lamb = dynamic_cast<lambertian*>(obj->mat_ptr.get());
            spheres[i].mat.albedo = float3(lamb->albedo.x(), lamb->albedo.y(), lamb->albedo.z());
            spheres[i].mat.materialType = 0;
        }
        else if (generic_mat->get_mat_type() == 1) {
            metal* met = dynamic_cast<metal*>(obj->mat_ptr.get());
            spheres[i].mat.albedo = float3(met->albedo.x(), met->albedo.y(), met->albedo.z());
            spheres[i].mat.fuzz_ir = met->fuzz;
            spheres[i].mat.materialType = 1;
        }
        else if (generic_mat->get_mat_type() == 2){
            dielectric* diel = dynamic_cast<dielectric*>(obj->mat_ptr.get());
            spheres[i].mat.albedo = float3(0, 0, 0);
            spheres[i].mat.fuzz_ir = diel->ir;
            spheres[i].mat.materialType = 2;
        }
    }

    cl_int error;
    float4* pixel_color = new float4[nPixels];

    // GPU Buffers
    cl_mem colorBuffer = clCreateBuffer(Kernel::GetContext(), CL_MEM_READ_WRITE, nPixels * sizeof(float4), 0, 0);
    cl_mem sphereBuffer = clCreateBuffer(Kernel::GetContext(), CL_MEM_READ_WRITE, world.objects.size() * sizeof(CL_Sphere), 0, 0);

    error = clEnqueueWriteBuffer(Kernel::GetQueue(), sphereBuffer, true, 0, world.objects.size() * sizeof(CL_Sphere), spheres, 0, 0, 0);

    int numOfSpheres = world.objects.size();
    float lens_radius = aperture / 2;

    kernel.SetArgument(0, &colorBuffer);
    kernel.SetArgument(1, &sphereBuffer);
    kernel.SetArgument(2, image_width);
    kernel.SetArgument(3, image_height);
    kernel.SetArgument(4, cam.get_vertical());
    kernel.SetArgument(5, cam.get_horizontal());
    kernel.SetArgument(6, cam.get_lower_left_corner());
    kernel.SetArgument(7, cam.get_origin());
    kernel.SetArgument(8, max_depth);
    kernel.SetArgument(9, numOfSpheres);     
    kernel.SetArgument(10, samples_per_pixel);
    kernel.SetArgument(11, lens_radius);
    float3 fw = float3(w.x(), w.y(), w.z());
    kernel.SetArgument(12, fw);
    fw = float3(u.x(), u.y(), u.z());
    kernel.SetArgument(13, fw);
    fw = float3(v.x(), v.y(), v.z());
    kernel.SetArgument(14, fw);

    Timer t;
    kernel.Run(image_width * image_height * samples_per_pixel, samples_per_pixel);
    clFinish(kernel.GetQueue());
    
    error = clEnqueueReadBuffer(Kernel::GetQueue(), colorBuffer, true, 0, nPixels * sizeof(float4), pixel_color, 0, 0, 0);
    std::cerr << error << ' ' << std::flush;

    // Render   
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            // RUN ON CPU
            /*color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);*/

            // RUN ON GPU
            int index = j * image_width + i;
            write_color(std::cout, color(pixel_color[index].x, pixel_color[index].y, pixel_color[index].z), samples_per_pixel);
        }
    }
       
    std::cerr << "frame time:" << t.elapsed() * 1000 << " msec" << std::endl;
    std::cerr << "\nDone.\n";
   
    system("pause");

}
