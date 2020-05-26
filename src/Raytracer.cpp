//Raytracer.cpp : Defines the entry point for the console application 
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <jsoncpp/json/json.h>
#include <fstream>
#include <iostream>

#include <string.h>
#include <stdio.h>
#include <vector>
#include "Vector.h"


using namespace std;

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

bool parseCommandLine(int argc, char** argv, int &level, char* fileIn, char* fileOut){
    if(argc != 7){
        return false;
    }
 

    for(int i = 1; (i+1) < argc; i+= 2){
        if(0 == strcmp("-i", argv[i])){
            strcpy(fileIn, argv[i+1]);
        }
         else if(0 == strcmp("-n", argv[i])){
            level = atoi(argv[i+1]);

        } else if(0 == strcmp("-o", argv[i])){
             strcpy(fileOut, argv[i+1]);
        }
    }
    return true;
}

void parseFile(const char* filename, Scene &scene, Vector &position_l, double &intensite_l ){
    std::ifstream ifs(filename);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj); // reader can also read strings

    const Json::Value& spheres = obj["spheres"]; // array of characters
    const Json::Value& rectangles = obj["rectangles"]; // array of characters
    const Json::Value& triangles = obj["triangles"]; // array of characters
    const Json::Value& cylindres = obj["cylindres"]; // array of characters
    
    cout << "nbsphere" << spheres.size() << endl;
    cout << "nbrect" << rectangles.size() << endl;
    for (int i = 0; i < spheres.size(); i++){
        Vector axe(Vector(spheres[i]["axeX"].asInt(), spheres[i]["axeY"].asInt(), spheres[i]["axeZ"].asInt()));
        int rayon =  spheres[i]["rayon"].asInt();
        Vector couleur(spheres[i]["couleur"][0].asInt(), spheres[i]["couleur"][1].asInt(), spheres[i]["couleur"][2].asInt());
        scene.addSphere(axe, rayon, couleur);
    }
    for(int i = 0; i < rectangles.size(); i++){
        double x = rectangles[i]["axeX"].asDouble();
        double y = rectangles[i]["axeY"].asDouble();
        double z = rectangles[i]["axeZ"].asDouble();
        double width = rectangles[i]["width"].asDouble();
        double height = rectangles[i]["height"].asDouble();
        Vector couleur(rectangles[i]["couleur"][0].asInt(), rectangles[i]["couleur"][1].asInt(), rectangles[i]["couleur"][2].asInt());
    
        scene.addRect(x,y,width,height,z,couleur);
    }

    intensite_l = obj["lumiere"]["intensite"].asDouble();
    position_l = Vector(obj["lumiere"]["x"].asDouble(), obj["lumiere"]["y"].asDouble(), obj["lumiere"]["z"].asDouble());
    
    cout << "nbforme" << scene.shapes.size() << endl;
    // s.addRect(0, 0, 30, 30, -30, Vector(0, 1, 0));
}

int main(int argc, char* argv[]){

    int level = 1;
    char fileIn[20];
    char fileOut[20];
    if(! parseCommandLine(argc, argv, level,fileIn, fileOut)){
        return -1;
    }
    
    int W = 1024; 
    int H = 1024; 
    double fov = 60 * M_PI / 180;
    Scene s; 

    Vector position_lumiere;
    double intensite_lumiere;
    parseFile(fileIn, s, position_lumiere, intensite_lumiere);

    

// diffenrent en fonction level
    switch (level)
    {
    case 2 :
        /* level2(H, w, SCENE, fileOUT, position_l, intensite_l, ) */
        break;
    
    default:
        break;
    }
    
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

    save_img(fileOut , &image[0], W, H);

    
    return 0;

}
