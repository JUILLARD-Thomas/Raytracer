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

int main()
{
    int W = 1024; 
    int H = 1024; 
    double fov = 60 * M_PI / 180;



    Scene s; 
    s.addSphere(Vector(30,30,-1), 21, Vector(1, 0, 0)); //boule rouge
    s.addSphere(Vector(1,5,-70), 30, Vector(1, 1, 0)); // boule jaune
    s.addSphere(Vector(0,-2000 -20, 0), 2000, Vector(1, 1, 1)); // sol
    s.addSphere(Vector(0,2000 + 100,0), 2000, Vector(1, 1, 1)); // plafond
    s.addSphere(Vector(-100 - 50,0,0), 100, Vector(0, 1, 0)); //mur droit
    s.addSphere(Vector(2000 + 50,0,0), 2000, Vector(0, 0, 1)); // mur gauche
    s.addSphere(Vector(0,0,-100 - 50), 100, Vector(0, 1, 1));
    
    s.addRect(0, 0, 30, 30, -30, Vector(0, 1, 0));

    Vector position_lumiere(15, 70, 30);   
    double intensite_lumiere = 1500000;

    std::vector<unsigned char> image(W*H *3);

    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction(j - W / 2, i - H / 2, -W / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(Vector(0,0,80), direction);
            Vector P, N;
           int shape_id;
         //  std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;
           bool has_inter = s.intersection(r,P,N,shape_id);
          // std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;

            Vector intensite_pixel(0,0,0);
            if (has_inter) {
                intensite_pixel = s.shapes[shape_id]->getAlbedo() * intensite_lumiere * std::max(0., dot((position_lumiere-P).getNormalized(), N)) / (position_lumiere -P).getNorm2();
            }

            image[((H - i -1)*W +j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0])); // rouge 
            image[((H - i -1)*W +j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1]));  // vert
            image[((H - i -1)*W +j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2]));  // bleu
        }
    
    }

    save_img("out.ppm" , &image[0], W, H);
    return 0;

}
