#pragma once
#include <math.h>
#include <vector>
#include <iostream>
#include <memory>

class Vector {
    public: 
        Vector(double x = 0, double y = 0, double z = 0) {
            coord[0] = x;
            coord[1] = y;
            coord[2] = z;
        }
        const double& operator[](int i) const {
            return coord[i];
        }

        double& operator[] (int i) {
            return coord[i];
        }

        double getNorm2() {
            return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
        }
        void normalize() {
            double norm = sqrt(getNorm2());
            coord[0] /= norm;
            coord[1] /= norm;
            coord[2] /= norm;
        }
        Vector getNormalized() {
            Vector result(*this);
            result.normalize();
            return result;
        }




    private:
        double coord[3];
};

Vector operator+(const Vector& a, const Vector &b);
Vector operator-(const Vector& a, const Vector &b);
Vector operator*(double a, const Vector &b);
Vector operator*( const Vector &b, double a);
Vector operator/(const Vector& a, double b);
double dot(const Vector&a, const Vector& b);

class Ray {
    public : 
        Ray(const Vector& o, const Vector&  d) : origin(o) , direction(d) {};
        Vector origin, direction;
};

class Shape {
    public:
        virtual bool intersection(const Ray d , Vector& P, Vector& N, double &t) = 0;
        virtual ~Shape(){};
        
        virtual Vector getAlbedo() = 0;

};

class Sphere : public Shape {
    public:

    Sphere(const Vector &origin, double rayon, const Vector &couleur): O(origin), R(rayon), albedo(couleur) {};

    bool intersection(const Ray d , Vector& P, Vector& N, double &t) {
    // resout a*t² + b*t + c = 0
        double a = 1 ; 
        double b = 2*dot(d.direction, d.origin - O);
        double c = (d.origin - O).getNorm2() - R*R;

        double delta = b*b - 4 * a*c;

        if(delta < 0 ) return false;
        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);

        if (t2 < 0) return false; 

        if (t1 > 0)
            t = t1;
        else 
            t = t2;

        P = d.origin + t*d.direction;
        N = (P - O).getNormalized();
       // std::cout << "sphere: t -> " << t << std::endl;

        return true;
}

    Vector getAlbedo(){
        return albedo;
    }

    Vector O;
    double R;
    Vector albedo;


};

class Rectangle : public Shape {
    public:
    Rectangle(double x, double y, double width, double height, double Z, const Vector &couleur) {
        x1 = x;
        x2 = x + width;
        y1 = y;
        y2 = y + height;
        z = Z;
        albedo = couleur;
    };

    bool intersection(const Ray d , Vector& P, Vector& N, double &t) {
        // std::cout << "rect ok" <<  std::endl;
        double tmp1 = d.origin[0] + d.direction[0] * (z - d.origin[2]) / d.direction[2];
        double tmp2 = d.origin[1] + d.direction[1] * (z - d.origin[2]) / d.direction[2];


        if(x1 <= tmp1 && tmp1 <= x2 && y1 <= tmp2 && tmp2 <= y2){
            t = (z - d.origin[2]) / d.direction[2];
            P = d.origin + t*d.direction;
            Vector v(0, 0, 1);
            N = v;
         //   std::cout << "rect: t -> " << t << std::endl;
            return true;
        }
    return false;
}

Vector getAlbedo(){
        return albedo;
    }

    double x1;
    double x2;
    double y1;
    double y2;
    double z;
    Vector albedo;
};

class Scene {
    public : 
        Scene() {};
        void addSphere(const Vector &origin, double rayon, const Vector &couleur) { shapes.push_back(std::unique_ptr<Shape>(new Sphere(origin, rayon, couleur))); }

        void addRect(double x, double y, double width, double height, double Z, const Vector &couleur) { shapes.push_back(std::unique_ptr<Shape>(new Rectangle(x, y, width, height, Z, couleur))); }


        bool intersection(const Ray d , Vector& P, Vector& N, int &shape_id) {

            bool has_inter = false;
            double min_t = 1E99;
           // std::cout << "lancer d'un rayon:" << std::endl;

            for (int i =0; i < shapes.size(); i++) {
                Vector localP, localN;
                double t;
               bool local_has_inter = shapes[i]->intersection(d, localP, localN, t);
              

               if(local_has_inter){
                   if(i >= 2){
                   //    std::cout << "POU\nLOU\n\n\nintersection:  " << i << "\n\n_n" << std::endl;
                   }
                    has_inter =true;
             //   std::cout << "intersection: " << i << std::endl;

                   if(t < min_t) {
                  //  std::cout << "lMais il faut choisir:" << t << " < " << min_t << std::endl;

                       min_t = t;
                       P = localP;
                       N = localN;
                       shape_id = i;
                   }else{
                      //  std::cout << "Mais pas ASSEZ fort " << t << " >" << min_t <<  std::endl;
                   }
               }
            }

         //   std::cout << "fin d'un rayon:" << std::endl;


            return has_inter;
        }

        std::vector<std::unique_ptr<Shape>> shapes; 
};