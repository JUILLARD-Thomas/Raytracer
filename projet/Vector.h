#pragma once
#include <math.h>

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

class Sphere {
    public:
    Sphere(const Vector &origin, double rayon, const Vector &couleur) : O(origin), R(rayon), albedo(couleur) {};

    Vector O;
    double R;
    Vector albedo;


};

class Rectangle {
    public:
    Rectangle(double x, double y, double width, double height, double Z, const Vector &couleur) : x1(x), x2(x + width), y1(y), y2(y+height), z(Z), albedo(couleur) {};

    double x1;
    double x2;
    double y1;
    double y2;
    double z;
    Vector albedo;
};