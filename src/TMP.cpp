#pragma once
#include <math.h>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>

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

Vector cross(const Vector&a, const Vector&b);

/*
Vector random_cos(const Vector &N) {
    int threadid = omp_get_thread_num();
    double r1 = uniform(engine[threadid]);
    double r2 = uniform(engine[threadid]);
    double sr2 = sqrt(1. - r2);
    Vector direction_aleatoire_repere_local(cos(2 * M_PIr1)sr2, sin(2 * M_PIr1)sr2, sqrt(r2));
    /Vector aleatoire(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);
    Vector tangent1 = cross(N, aleatoire); tangent1.normalize();/
    Vector tangent1;
    Vector absN(abs(N[0]), abs(N[1]), abs(N[2]));
    if (absN[0] <= absN[1] && absN[0]<=absN[2]) {
        tangent1 = Vector(0, -N[2], N[1]);
    }else
        if (absN[1] <= absN[0] && absN[1]<=absN[2]) {
            tangent1 = Vector(-N[2], 0, N[0]);
        } else
            tangent1 = Vector(-N[1], N[0], 0);
    tangent1.normalize();
    Vector tangent2 = cross(tangent1, N);

    return direction_aleatoire_repere_local[2] * N + direction_aleatoire_repere_local[0] * tangent1 + direction_aleatoire_repere_local[1] * tangent2;
}*/


class Ray {
    public : 
        Ray(const Vector& o, const Vector&  d) : origin(o) , direction(d) {};
        Vector origin, direction;
};

class Shape {
    public:
        Shape(Vector coul, bool trans) : albedo(coul), transp(trans) {};
        virtual bool intersection(const Ray d , Vector& P, Vector& N, double &t) = 0;
        
        
        Vector albedo;
        bool transp;

};

class Sphere : public Shape {
    public:

    Sphere(const Vector &origin, double rayon, const Vector &couleur): O(origin), R(rayon), Shape(couleur, false) {
    };

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

    Vector O;
    double R;


};

class Triangle : public Shape {
    public: 
    Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector &couleur, bool mirror = false, bool transp = false) : A(A), B(B), C(C), Shape(couleur, false){};

    bool intersection(const Ray d , Vector& P, Vector& N, double &t) {
        N = cross(B- A, C- A).getNormalized();
        t = dot(C - d.origin, N ) /dot(d.direction, N);
        if (t < 0) return false;

        P = d.origin + t*d.direction;
        Vector u = B - A ; 
        Vector v = C - A ; 
        Vector w = P - A; 

        double m11 = u.getNorm2();
        double m12 = dot(u,v);
        double m22 = v.getNorm2();
        double detm = m11 *m22 - m12*m12;

        double b11 = dot(w,u);
        double b21 = dot(w,v);
        double detb = b11*m22 - b21*m12;
        double beta = detb /detm; //coord barycentrique w.r.t à B

        double g12 = b11;
        double g22 = b21; 
        double detg = m11*g22 - m12*g12;
        double gamma = detg / detm; //coord barycentrique w.r.t à C

        double alpha = 1 - beta - gamma;
        if(alpha < 0 || alpha > 1) return false; 
        if(beta < 0 || beta > 1) return false; 
        if(gamma < 0 || gamma > 1) return false; 
        if(alpha + beta + gamma > 1) return false; 


        return true;
    }


    Vector A, B, C, albedo;
    bool mirror, transp;    

};


class Rectangle : public Shape {
    public:
    Rectangle(const Vector& A, const Vector& B, const Vector& C, const Vector& D, const Vector &couleur): A(A), B(B), C(C), D(D), Shape(couleur, false) {

    };
    
    bool intersection(const Ray d , Vector& P, Vector& N, double &t) {
        N = cross(B- A, C- A).getNormalized();
        t = dot(C - d.origin, N ) /dot(d.direction, N);
        if (t < 0) return false;

        P = d.origin + t*d.direction;

        
       Vector V1 = (B - A).getNormalized();
       Vector V2 = (C - B).getNormalized();
       Vector V3 = (D - C).getNormalized();
       Vector V4 = (A - D).getNormalized();
       Vector V5 = (P - A).getNormalized();
       Vector V6 = (P - B).getNormalized();
       Vector V7 = (P - C).getNormalized();
       Vector V8 = (P - D).getNormalized();

       if (dot(V1, V5) < 0.0) return false;
       if (dot(V2, V6) < 0.0) return false;
       if (dot(V3, V7) < 0.0) return false;
       if (dot(V4, V8) < 0.0) return false;
        return true;
    }
    
/*
        if(!((A[0] <= P[0] && P[0] <= B[0]) || (A[0] >= P[0] && P[0] >= B[0]))){
            return false;
        } else if(!((A[1] <= P[1] && P[1] <= C[1]) || (A[1] >= P[1] && P[1] >= C[1]))){
            return false;
        } else if(!((A[2] <= P[2] && P[2] <= D[2]) || (A[2] >= P[2] && P[2] >= D[2]))){
            return false;
        }
        
         
    return true;
    


    float t = -(-Origin.y + 20) / (-Dir.y);
    if (t > 0) {
        Vector hitPoint = Origin + Dir * t;
        Vec V1 = (P2 - P1).norm();
        Vec V3 = (P4 - P3).norm();
        Vec V4 = (hitPoint - P1).norm();
        Vec V5 = (hitPoint - P3).norm();
        float V1dotV4 = V1.dot(V4);
        float V3dotV5 = V3.dot(V5);
        if (V1dotV4 > 0 && V3dotV5 > 0) {
            return true;
        }
    }
    return false;
}
*/



    Vector A, B, C, D;
};

class Scene {
    public : 
        Scene() {};
        void addSphere(const Vector &origin, double rayon, const Vector &couleur) { shapes.push_back(std::unique_ptr<Shape>(new Sphere(origin, rayon, couleur))); }

        void addRect(const Vector& A, const Vector& B, const Vector& C,const Vector& D, const Vector &couleur) { shapes.push_back(std::unique_ptr<Shape>(new Rectangle(A, B, C, D, couleur))); }
        
        void addTriangle(const Vector& A, const Vector& B, const Vector& C, const Vector &couleur, bool mirror = false, bool transp = false){(shapes.push_back(std::unique_ptr<Shape>(new Triangle(A, B, C, couleur, mirror, transp)))); }
        
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
        Sphere* lumiere;
        double intensite_lumiere;
};