//Raytracer.cpp : Defines the entry point for the console application 
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <GL/glut.h>


#include <jsoncpp/json/json.h>
#include <fstream>
#include <iostream>

#include <string.h>
#include <stdio.h>
#include <vector>
#include "Vector.h"


using namespace std;

#define M_PI 3.1415926535897
int H = 1024;
int W = 1024;

int level = 1;
char fileIn[150];
char fileOut[150];
double fov = 60 * M_PI / 180;
Scene s;
Vector camPosition(0,0,80);
Vector camDirect(0,0,0);


//Vector position_lumiere;
//double intensite_lumiere; 



unsigned char data[1024*1024*3];
 
 /* l'image dans la fenêtre openGL */ 
void displayMe()
{
    cout << "go disPLAYYYYYY" <<endl;
    glClearColor( 0, 0, 0, 1 );
    glClear(GL_COLOR_BUFFER_BIT);

    glDrawPixels( W, H, GL_RGB, GL_UNSIGNED_BYTE, data );

    glutSwapBuffers();
    
}

/* Sauvegarde de l'image sous le format ppm */ 

void save_img_ppm(const char* filename, const unsigned char* pixels){
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

/* Défini les arguments à mettre dans les lignes de commandes */ 
bool parseCommandLine(int argc, char** argv, int &level, char* fileIn, char* fileOut){
    if(argc != 7){
        return false;
    }
 

    for(int i = 1; (i+1) < argc; i+= 2){
        cout << "avant cmp " << argc << endl;
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
/* fonction permettant de parser grace aux .json */
void parseFile(const char* filename, Scene &scene){
    std::ifstream ifs(filename);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);

    const Json::Value& spheres = obj["spheres"]; 
    const Json::Value& rectangles = obj["rectangles"];
    const Json::Value& triangles = obj["triangles"];
    const Json::Value& cylindres = obj["cylindres"]; 
    
    cout << "nbsphere" << spheres.size() << endl;
    cout << "nbrect" << rectangles.size() << endl;
    cout << "nbTRITRI ->" << triangles.size() << endl;

    for (int i = 0; i < spheres.size(); i++){
        bool mirror = (i == 0) ? true : false;
        Vector axe(Vector(spheres[i]["axeX"].asInt(), spheres[i]["axeY"].asInt(), spheres[i]["axeZ"].asInt()));
        int rayon =  spheres[i]["rayon"].asInt();
        Vector couleur(spheres[i]["couleur"][0].asInt(), spheres[i]["couleur"][1].asInt(), spheres[i]["couleur"][2].asInt());
        scene.addSphere(axe, rayon, (i == 0) ? true : false, couleur);
    }
    for(int i = 0; i < rectangles.size(); i++){
        
        Vector a(rectangles[i]["A"][0].asInt(),rectangles[i]["A"][1].asInt(), rectangles[i]["A"][2].asInt());
        Vector b(rectangles[i]["B"][0].asInt(),rectangles[i]["B"][1].asInt(), rectangles[i]["B"][2].asInt());
        Vector c(rectangles[i]["C"][0].asInt(),rectangles[i]["C"][1].asInt(), rectangles[i]["C"][2].asInt());
        Vector d(rectangles[i]["D"][0].asInt(),rectangles[i]["D"][1].asInt(), rectangles[i]["D"][2].asInt());
        Vector couleur(rectangles[i]["couleur"][0].asInt(), rectangles[i]["couleur"][1].asInt(), rectangles[i]["couleur"][2].asInt());
    
        scene.addRect(a,b,c,d, false,couleur);
    }

    for(int i = 0; i < triangles.size(); i++){
        Vector x(triangles[i]["x"][0].asInt(), triangles[i]["x"][1].asInt(), triangles[i]["x"][2].asInt());
        Vector y(triangles[i]["y"][0].asInt(), triangles[i]["y"][1].asInt(), triangles[i]["y"][2].asInt());
        Vector z(triangles[i]["z"][0].asInt(), triangles[i]["z"][1].asInt(), triangles[i]["z"][2].asInt());
        Vector couleur(triangles[i]["couleur"][0].asInt(), triangles[i]["couleur"][1].asInt(), triangles[i]["couleur"][2].asInt());
    
        scene.addTriangle(x,y,z,false, couleur);
    }

    for (int i = 0; i < cylindres.size(); i++){

        Vector pointA(cylindres[i]["pointA"][0].asInt(),cylindres[i]["pointA"][1].asInt(),cylindres[i]["pointA"][2].asInt());
        double rayon = cylindres[i]["rayon"].asDouble();
        Vector vectV(cylindres[i]["V"][0].asInt(),cylindres[i]["V"][1].asInt(),cylindres[i]["V"][2].asInt());
        double hauteur = cylindres[i]["hauteur"].asDouble();
        
        Vector couleur(cylindres[i]["couleur"][0].asInt(),cylindres[i]["couleur"][1].asInt(),cylindres[i]["couleur"][2].asInt());
         cout << "A: " << pointA[0] << " " << pointA[1] << " " << pointA[2] << endl;
         cout << "V: " << vectV[0] << " " << vectV[1] << " " << vectV[2] << endl;
         cout << "couleur: " << couleur[0] << " " << couleur[1] << " " << couleur[2] << endl;
         cout << "h: " << hauteur << endl;
         cout << "r: " << rayon << endl;
        scene.addCylindre(pointA, rayon, vectV, hauteur, false, couleur);

    }
    s.intensite_lumiere = obj["lumiere"]["intensite"].asDouble();
    s.position_lumiere = Vector(obj["lumiere"]["x"].asDouble(), obj["lumiere"]["y"].asDouble(), obj["lumiere"]["z"].asDouble());
    //intensite_l = obj["lumiere"]["intensite"].asDouble();
    //position_l = Vector(obj["lumiere"]["x"].asDouble(), obj["lumiere"]["y"].asDouble(), obj["lumiere"]["z"].asDouble());
    
    std::cout << "nbforme" << scene.shapes.size() << endl;
    // s.addRect(0, 0, 30, 30, -30, Vector(0, 1, 0));
}



Vector getColor(const Ray &r, Scene &s, int nb_rebonds){
     if(nb_rebonds == 0){
        return Vector(0,0,0);
    }

     Vector P, N;
    int shape_id;
    double t;
    bool has_inter = s.intersection(r,P,N,shape_id, t);

    Vector intensite_pixel(0,0,0);
    if (has_inter) {
         if(s.shapes[shape_id]->isMirror){
           // cout << "mirror" << endl;
            Vector direction_mirror = r.direction - 2 * dot(N, r.direction) * N;
            Ray rayon_mirror(P + 0.001*N, direction_mirror);
            intensite_pixel = getColor(rayon_mirror, s, nb_rebonds - 1);
        } else{
            Ray ray_light(P + 0.01 * N, (s.position_lumiere - P).getNormalized());
            Vector P_light,N_light;

            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
            double d_light2 = (s.position_lumiere - P ).getNorm2();
            if (has_inter_light && t_light * t_light < d_light2){
                intensite_pixel = Vector(0,0,0);
            }
            else{
                intensite_pixel = s.shapes[shape_id]-> albedo * s.intensite_lumiere * std::max(0., dot((s.position_lumiere-P).getNormalized(), N)) /d_light2;
            }
        }
        
        
    }
    return intensite_pixel;
/*

    if (has_inter) {

        if(s.shapes[shape_id]->isMirror){
           // cout << "mirror" << endl;
            Vector direction_mirror = r.direction - 2 * dot(N, r.direction) * N;
            Ray rayon_mirror(P + 0.001*N, direction_mirror);
            intensite_pixel = getColor(rayon_mirror, s, nb_rebonds - 1);
        } else{
          //  cout << "normal " << endl;
            Ray ray_light(P + 0.01 * N, (s.position_lumiere - P).getNormalized());
            Vector P_light,N_light;

            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
            double d_light2 = (s.position_lumiere - P ).getNorm2();
            if (has_inter_light && t_light * t_light < d_light2)
            {
                intensite_pixel = Vector(0,0,0);
            }
            else
            {
                intensite_pixel = s.shapes[shape_id]-> albedo * s.intensite_lumiere * std::max(0., dot((s.position_lumiere-P).getNormalized(), N)) /d_light2;
            }
        }
    }


return intensite_pixel;*/
    
}




/* définir le niveau 1 d'exécution du projet */ 
void levelOne(int H,int W,Scene &s,char* fileOut, double fov){
    std::vector<unsigned char> image(W*H *3);

    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction(j - W / 2, i - H / 2, -W / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(Vector(0,0,80), direction);
            Vector P, N;
            int shape_id;
            double t;
         //  std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;
            bool has_inter = s.intersection(r,P,N,shape_id, t);
          // std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;

            Vector intensite_pixel(0,0,0);
            if (has_inter) {
                intensite_pixel = s.shapes[shape_id]->albedo;
                intensite_pixel[0] = (intensite_pixel[0] > 0) ? 255 : 0;
                intensite_pixel[1] = (intensite_pixel[1] > 0) ? 255 : 0;
                intensite_pixel[2] = (intensite_pixel[2] > 0) ? 255 : 0; 
            }

            image[((H - i -1)*W +j) * 3 + 0] =  intensite_pixel[0]; // rouge 
            image[((H - i -1)*W +j) * 3 + 1] = intensite_pixel[1];  // vert
            image[((H - i -1)*W +j) * 3 + 2] = intensite_pixel[2];  // bleu
        }
    }

    save_img_ppm(fileOut , &image[0]);
}
/* définir le niveau 2 d'exécution du projet */ 
void levelTwo(int H,int W,Scene &s,char* fileOut, double fov){
    std::vector<unsigned char> image(W*H *3);


    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction((j + camDirect[0] - W  / 2), i + camDirect[1] - H / 2, camDirect[0] -W  / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(Vector(0,0,80), direction);
            Vector P, N;
            int shape_id;
            double t;
         //  std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;
            bool has_inter = s.intersection(r,P,N,shape_id, t);
          // std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;

            Vector intensite_pixel(0,0,0);
            if (has_inter) {
               // cout << "INTERSECTION !!!!" << endl;
              //  std::cout << "intersection " << N[0] << " " << N[1] << " " << N[2] << std::endl;
                intensite_pixel = s.shapes[shape_id]-> albedo * s.intensite_lumiere * std::max(0., dot((s.position_lumiere-P).getNormalized(), N)) / (s.position_lumiere -P).getNorm2();
            }

            image[((H - i -1)*W +j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0])); // rouge 
            image[((H - i -1)*W +j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1]));  // vert
            image[((H - i -1)*W +j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2]));  // bleu
        }
    }

    save_img_ppm(fileOut , &image[0]);
}


void level333(int H,int W,Scene &s,char* fileOut, double fov){
    std::vector<unsigned char> image(W*H *3);
    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction((j + camDirect[0] - W  / 2), i + camDirect[1] - H / 2, camDirect[0] -W  / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(Vector(0,0,80), direction);
           

            Vector intensite_pixel = getColor(r, s, 5);

            image[((H - i -1)*W +j) * 3 + 0] = std::min(255., std::max(0., std::pow(intensite_pixel[0], 1/2.2))); // rouge 
            image[((H - i -1)*W +j) * 3 + 1] = std::min(255., std::max(0., std::pow(intensite_pixel[1], 1/2.2)));  // vert
            image[((H - i -1)*W +j) * 3 + 2] = std::min(255., std::max(0., std::pow(intensite_pixel[2], 1/2.2)));  // bleu
        }
    }

    save_img_ppm(fileOut , &image[0]);
}


/* définir le niveau 2 bis d'exécution du projet */ 
void levelThree3(int H,int W,Scene &s,char* fileOut, double fov){
    std::vector<unsigned char> image(W*H *3);
    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction((j + camDirect[0] - W  / 2), i + camDirect[1] - H / 2, camDirect[0] -W  / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(Vector(0,0,80), direction);
            Vector P, N;
            int shape_id;
            double t;
            bool has_inter = s.intersection(r,P,N,shape_id, t);

            Vector intensite_pixel(0,0,0);
            if (has_inter) {
              
                Ray ray_light(P + 0.01 * N, (s.position_lumiere - P).getNormalized());
                Vector P_light,N_light;

                int sphere_id_light;
                double t_light;
                bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);
                double d_light2 = (s.position_lumiere - P ).getNorm2();
                if (has_inter_light && t_light * t_light < d_light2){
                    intensite_pixel = Vector(0,0,0);
                }
                else{
                    intensite_pixel = s.shapes[shape_id]-> albedo * s.intensite_lumiere * std::max(0., dot((s.position_lumiere-P).getNormalized(), N)) /d_light2;
                }
            }

            

            image[((H - i -1)*W +j) * 3 + 0] = std::min(255., std::max(0., std::pow(intensite_pixel[0], 1/2.2))); // rouge 
            image[((H - i -1)*W +j) * 3 + 1] = std::min(255., std::max(0., std::pow(intensite_pixel[1], 1/2.2)));  // vert
            image[((H - i -1)*W +j) * 3 + 2] = std::min(255., std::max(0., std::pow(intensite_pixel[2], 1/2.2)));  // bleu
        }
    }

    save_img_ppm(fileOut , &image[0]);
}

/* définir le niveau 3 d'exécution du projet */ 
void levelThree(int H,int W,Scene &s,char* fileOut, double fov){
   
   // std::vector<unsigned char> image(W*H *3);
    std::cout << "On recalcule TOUUUUUUUUUUUT " << camPosition[0] << endl;
    for (int i =0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector direction((j + camDirect[0] - W  / 2), i + camDirect[1] - H / 2, camDirect[0] -W  / (2 * tan(fov / 2)));
            direction.normalize();

            Ray r(camPosition, direction);

            Vector P, N;
            int shape_id;
            double t;
         //  std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;
            bool has_inter = s.intersection(r,P,N,shape_id, t);
          // std::cout << "test intersection: " << i << "," << j << "\n" << std::endl;

            Vector intensite_pixel(0,0,0);
            if (has_inter) {
                intensite_pixel = s.shapes[shape_id]-> albedo * s.intensite_lumiere * std::max(0., dot((s.position_lumiere-P).getNormalized(), N)) / (s.position_lumiere -P).getNorm2();
            }

            data[(i*W +j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0])); // rouge 
            data[(i*W +j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1]));  // vert
            data[(i*W +j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2]));  // bleu
        }
    }

}
/* fonctions de détections des touches pour le déplacement de la caméra */
void vClavier(unsigned char key, int x, int y){
    bool move = false;
	switch (key){
	
	/* Mode de déplacement vers la gauche */
		case 'q' :
		case 'Q' :
            move = true;
            cout << "Un peu a gauche" << endl;
            camPosition[0] += -10;
            break;

    /* Mode de déplacement vers la droite */
		case 'd' :
		case 'D' :
            move = true;
            cout << "Un peu a doite" << endl;
            camPosition[0] += 10;
            break;

    /* Mode de déplacement en haut */
		case 'z' :
		case 'Z' :
            move = true;
            cout << "Un peu a haut" << endl;
            camPosition[1] += 10;
            break;

    /* Mode de déplacement en bas */
		case 's' :
		case 'S' :
            move = true;
            cout << "Un peu a BAS" << endl;
            camPosition[1] += -10;
            break;

    /* Mode de déplacement en profondeur -  */
		case 'w' :
		case 'W' :
            move = true;
            cout << "Un peu a -profond" << endl;
            camPosition[2] += 10;
            break;
    
    /* Mode de déplacement en profondeur + */
		case 'x' :
		case 'X' :
            move = true;
            cout << "Un peu a +profond" << endl;
            camPosition[2] -= 10;
            break;

    /* inclinaison de la caméra en haut */
        case 'r' :
		case 'R' :
            move = true;
            cout << "Tema e haut" << endl;
            camDirect[1] += 40;
            break;
    /* inclinaison de la caméra en bas */ 
        case 'f' :
		case 'F' :
            move = true;
            cout << "Tema en bas" << endl;
            camDirect[1] -= 40;
            break;

    /* inclinaison de la caméra vers la droite */
        case 'c' :
		case 'C' :
            move = true;
            cout << "Tema a droite" << endl;
            camDirect[0] -= 20;
            break;

    /* inclinaison de la caméra vers la gauche */ 
        case 'v' :
		case 'V' :
            move = true;
            cout << "Tema a gauche" << endl;
            camDirect[0] += 20;
            break;

		default :
			printf(" Touche: %c\n Souris a: %d %d \n",key,x,y); 
			break;
	}
    if(move){
        levelThree( H, W, s, fileOut,fov );
        glutPostRedisplay();
    }
}



/* fonction de main pour lancer le projet */ 
int main(int argc, char* argv[]){

    
    if(! parseCommandLine(argc, argv, level,fileIn, fileOut)) return -1;

    cout << "level: " << level << "\nfileIn: " << fileIn << "\nfileOut: " << fileOut << endl;


    for(int i = 0; i < argc; i++){
        std::cout << i << " " << argv[i] << std::endl;
    }    

    
    parseFile(fileIn, s);
    cout << "nbforme " << s.shapes.size() << endl;
   

   

    switch (level)
    {
    case 1 :
        /* Fonction pour le niveau 1 */
        levelOne( H, W, s, fileOut,fov );
        break;
    
    case 2 :
        /* Fonction pour le niveau 2 */
        levelTwo( H, W, s, fileOut, fov );
        break;

    case 3 :
        /* Fonction pour le niveau 3 */
        levelThree( H, W, s, fileOut, fov );
        glutInit(&argc, argv); //init la lib glut
        glutInitDisplayMode(GLUT_SINGLE); //mask on touche pas
        glutInitWindowSize(1024, 1024); //size
        glutCreateWindow("Hello world!");

        glutDisplayFunc(displayMe);
        glutKeyboardFunc(vClavier); //poulouLou  
        
        glutMainLoop(); //sert de wait
        break;
    case 4 :
       level333( H, W, s, fileOut, fov );
        break;
    
    default:
        cout << "Level: " << level << "doesn't exist" << endl;
        break;
    }

    
    return 0;

}
