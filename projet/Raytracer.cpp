//Raytracer.cpp : Defines the entry point for the console application 
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <stdio.h>
#include <vector>
#include "Vector.h"
#include <iostream>

#define M_PI 3.1415926535897

void save_img(const char* filename, const unsigned char* pixels, int W, int H){
    FILE* F = fopen(filename,"w");
    if (!F)
        return ;

    fprintf(F,"P3\n%d %d\n255\n", W, H);
   for(int i = 0; i < H; i ++){
       for(int j = 0; j < W; j++){
           fprintf(F, "%d %d %d\n", pixels[(i*W +j) * 3 + 0], pixels[(i*W +j) * 3 + 1], pixels[(i*W +j) * 3 + 2]);
       }
   }
   fclose(F);
   


}
bool intersectionSphere(const Ray d , const Sphere& s, Vector& P, Vector& N) {
    // resout a*t² + b*t + c = 0 
    double a = 1 ; 
    double b = 2*dot(d.direction, d.origin - s.O);
    double c = (d.origin - s.O).getNorm2() - s.R*s.R;

    double delta = b*b - 4 * a*c;
    if(delta < 0 ) return false;
    double t1 = (-b - sqrt(delta)) / (2 * a);
    double t2 = (-b + sqrt(delta)) / (2 * a);

    if (t2 < 0) return false; 

    double t; 
    if (t1 > 0)
        t = t1;
    else 
        t = t2;

    P = d.origin + t*d.direction;
    N = (P - s.O).getNormalized();

    return true;
}

bool intersectionRect(const Ray d , const Rectangle& rect, Vector& P, Vector& N) {
    double tmp1 = d.origin[0] + d.direction[0] * (rect.z - d.origin[2]) / d.direction[2];
    double tmp2 = d.origin[1] + d.direction[1] * (rect.z - d.origin[2]) / d.direction[2];

    if(rect.x1 <= tmp1 && tmp1 <= rect.x2 && rect.y1 <= tmp2 && tmp2 <= rect.y2){
        double t = (rect.z - d.origin[2]) / d.direction[2];
        P = d.origin + t*d.direction;
        Vector v(0, 0, 1);
        N = v;
        return true;
    }
   /* double t = (rect.z - d.origin[2]) / d.direction[2];
    P = d.origin + t*d.direction;
    Vector v(0, 0, -1);
    N = v;
    */
    return false;
}

int main()
{
    int W = 1024; 
    int H = 1024; 
    double fov = 60 * M_PI / 180;

    Sphere s(Vector(0,0,-55), 20, Vector(1, 0, 0));
    Rectangle rect(-20, -10, 40, 20, -65, Vector(1, 0, 0));

    Vector position_lumiere(-15, -70, -40);   
    double intensite_lumiere = 1500000;

    std::vector<unsigned char> image(W*H *3);

    /*
    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction(j - W / 2, i - H / 2, -W / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(Vector(0,0,0), direction);
            Vector P, N;
            bool has_inter = intersectionSphere(r,s,P,N);

            double intensite_pixel = 0; 
            if (has_inter) {
                intensite_pixel = intensite_lumiere * std::max(0., dot((position_lumiere-P).getNormalized(), N)) / (position_lumiere -P).getNorm2();
            }

            image[((H - i -1)*W +j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel)); // rouge 
            image[((H - i -1)*W +j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel));  // vert
            image[((H - i -1)*W +j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel));  // bleu
        }

    */

    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction(j - W / 2, i - H / 2, -W / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(Vector(0,0,0), direction);
            Vector P, N;
           bool has_inter = intersectionRect(r,rect,P,N);
         //   bool has_inter = intersectionSphere(r,s,P,N);

            Vector intensite_pixel(0,0,0);
            //double intensite_pixel = 0;
            if (has_inter) {
                intensite_pixel = s.albedo * intensite_lumiere * std::max(0., dot((position_lumiere-P).getNormalized(), N)) / (position_lumiere -P).getNorm2();
              //  intensite_pixel = intensite_lumiere * std::max(0., dot((position_lumiere-P).getNormalized(), N)) / (position_lumiere -P).getNorm2();
               //  std::cout << "intensité pixel: " << intensite_pixel << std::endl;
            }

            image[((H - i -1)*W +j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0])); // rouge 
            image[((H - i -1)*W +j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1]));  // vert
            image[((H - i -1)*W +j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2]));  // bleu

            //image[((H - i -1)*W +j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel)); // rouge 
            //image[((H - i -1)*W +j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel));  // vert
            //image[((H - i -1)*W +j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel));  // bleu
        }
    
    }

    save_img("out.ppm" , &image[0], W, H);
    return 0;

}
