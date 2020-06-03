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
        /* Retourne une coordonée :  x ou y ou z */
        const double& operator[](int i) const {
            return coord[i];
        }
        /* Retourne une coordonée :  x ou y ou z */
        double& operator[] (int i) {
            return coord[i];
        }
        /* renvoie la norme au carré */
        double getNorm2() {
            return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
        }
        /* Sert à normaliser  x + y + z = 1 */
        void normalize() {
            double norm = sqrt(getNorm2());
            coord[0] /= norm;
            coord[1] /= norm;
            coord[2] /= norm;
        }
        /* renvoie la Norme */ 
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


class Ray {
    public : 
        Ray(const Vector& o, const Vector&  d) : origin(o) , direction(d) {};
        Vector origin, direction;
};

class Shape {
    public:
        Shape(Vector coul, bool trans, bool isMirror) : albedo(coul), transp(trans), isMirror(isMirror) {};
        virtual bool intersection(const Ray d , Vector& P, Vector& N, double &t) = 0;
        
        
        Vector albedo;
        bool transp;
        bool isMirror;

};

class Sphere : public Shape {
    public:

    Sphere(const Vector &origin, double rayon, bool isMirror, const Vector &couleur): O(origin), R(rayon), Shape(couleur, false, isMirror) {
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
    Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector &couleur, bool isMirror, bool transp = false) : A(A), B(B), C(C), Shape(couleur, false, isMirror){};

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
    Rectangle(const Vector& A, const Vector& B, const Vector& C, const Vector& D, bool isMirror, const Vector &couleur): A(A), B(B), C(C), D(D), Shape(couleur, false, isMirror) {

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

    Vector A, B, C, D;
};

class Cylindre : public Shape {
    public:
    Cylindre(const Vector& A, double r, const Vector& V, double h, bool isMirror, const Vector &couleur): C(A), r(r), V(V), h(h), Shape(couleur, false, isMirror) {
        Vector v(V[0], V[1], V[2]);
        B = h * v.getNormalized() + A; 
        V2 = Vector(-V[0],-V[1],-V[2]);
    };
    
    bool intersection(const Ray d , Vector& P, Vector& N, double &t) {
        //std::cout << "TEST" << std::endl;
        Vector L= d.origin - C;
        Vector w = cross(d.direction, V);
        double w2 = w.getNorm2();

        if(w2 == 0){
            double a = dot(L, V);
            Vector D = L - a * V;
            double d2 = D.getNorm2();
            if(d2 > (r * r)){
                return false;
            }
        }
        Vector wn = w.getNormalized();
        double R = abs(dot(L,wn));
        if(R > r){
            return false;
        }
        Vector E = cross(L, V);
        t = -(dot(E, wn))/sqrt(w2); 
        Vector F = cross(wn, V);
        Vector Fn = F.getNormalized();
        double s = sqrt(r*r - R*R)/ abs(dot(d.direction, Fn));
      //  Vector P1 = d.origin + dot(t - s, d.direction);
     //   Vector P2 = d.origin + dot(t + s, d.direction);
         Vector P1 = d.origin + (t - s) * d.direction;
        Vector P2 = d.origin + (t + s) * d.direction;
        if(dot(L,V) < r){
            P = P2;
        } else{
            P = P1;
        }
       // P = ((P1 - d.origin).getNorm2() < (P2 - d.origin).getNorm2()) ? P1 : P2;
      //  std::cout << "intersection " << t <<  std::endl;
   //   P[2] = -P[2]; 
        Vector CP = P - C;
        double CQ = dot(CP, V);
        Vector QP = CP - CQ * V;
        N = QP/r;
     //   N[2] = -N[2];
        Vector MYP1 = P - C;
        Vector MYP2 = P - B;
       // std::cout << "\nB " << B[0] << " " << B[1] << " " << B[2] <<  std::endl;
       // std::cout << "Vect BC " << MYP2[0] << " " << MYP2[1] << " " << MYP2[2] <<  std::endl;
       // std::cout << "V2 " << V2[0] << " " << V2[1] << " " << V2[2] <<  std::endl;
        
        if(dot(MYP1,V) >= 0 && dot(MYP2, V2) >= 0){
            return true;
        }
       
      //  N = Vector(0,0,1);

        return false;

        
    }

    Vector C, V,V2, B;
    double r,h;
};

class Scene {
    public : 
        Scene() {};
        void addSphere(const Vector &origin, double rayon, bool isMirror, const Vector &couleur) { shapes.push_back(std::unique_ptr<Shape>(new Sphere(origin, rayon, isMirror, couleur))); }

        void addRect(const Vector& A, const Vector& B, const Vector& C,const Vector& D, bool isMirror, const Vector &couleur) { shapes.push_back(std::unique_ptr<Shape>(new Rectangle(A, B, C, D, isMirror, couleur))); }
        
        void addTriangle(const Vector& A, const Vector& B, const Vector& C, bool isMirror, const Vector &couleur, bool transp = false){(shapes.push_back(std::unique_ptr<Shape>(new Triangle(A, B, C, couleur, isMirror, transp)))); }

        void addCylindre(const Vector& A, double r, const Vector& V, double h, bool isMirror, const Vector &couleur){shapes.push_back(std::unique_ptr<Shape>(new Cylindre(A,r,V,h, isMirror,couleur)));}
        
        bool intersection(const Ray d , Vector& P, Vector& N, int &shape_id, double &min_t) const {

            bool has_inter = false;
            min_t = 1E99;
           // std::cout << "lancer d'un rayon:" << std::endl;

            for (int i =0; i < shapes.size(); i++) {
                Vector localP, localN;
                double t;
               bool local_has_inter = shapes[i]->intersection(d, localP, localN, t);
              

               if(local_has_inter){

                    has_inter =true;

                   if(t < min_t) {

                       min_t = t;
                       P = localP;
                       N = localN;
                       shape_id = i;
                   }

               }
            }

         //   std::cout << "fin d'un rayon:" << std::endl;

            return has_inter;
        }

        std::vector<std::unique_ptr<Shape>> shapes;
        Sphere* lumiere;
        double intensite_lumiere;
        Vector position_lumiere;
};