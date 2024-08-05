// ------------------------------------------ //
//                                            //
//              Batardeau Engine              //
//                   *---*                    //
//                                            //
// Eng: NSI                                   //
// Ent: Terrasol (Setec)                      //
// Ver: 0.2                                   //
// Dat: 2024/08/03                            //
//                                            //
// ------------------------------------------ //

// Ce code constitue une première implémentation
// d'un code de calcul par éléments finis mettant
// en jeu des éléments de coque mince rectangulaires.

// Ce code est un miroir de celui développé en Python,
// conçu dans l'objectif de mesurer le gain de performances
// rendu possible par l'utilisation de C++.

// Pour cette première version, il n'intègre pas d'éléments
// annexes, tels que des poutres et des ancrages.

// Versions:
// +++++++++
//      v0.1 - Première implémentation du code pour coque seule, retranscription du code Python.
//      v0.2 - Amélioration de l'étape d'assemblage, en internalisant la multiplication de rotation.
//             On passe de 25% à 6-7% en part de temps pour l'assemblage.



// I - Importation des librairies nécessaires
// ------------------------------------------

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>

using namespace Eigen;
using namespace std;

#include <chrono>
using namespace std::chrono;

// II - Création des classes adaptées
// ----------------------------------

// a) Classe "vecteur" customisée

class vecteur{
    public:
    double dx,dy,dz;
    vecteur(){};
    vecteur(double i, double j, double k){
        dx = i;
        dy = j;
        dz = k;
    };
    double norm(){
        return sqrt(dx*dx + dy*dy + dz*dz);
    };
    vecteur normalize(){
        vecteur vec;
        double norme = sqrt(dx*dx + dy*dy + dz*dz);
        vec.dx=dx/norme;vec.dy=dy/norme;vec.dz=dz/norme;
        return vec;
    };
    double Pscalar(vecteur vec){
        return vec.dx*dx + vec.dy*dy + vec.dz*dz;
    };
    vecteur Pvector(vecteur vec){
        return vecteur(dy*vec.dz-dz*vec.dy,dz*vec.dx-dx*vec.dz,dx*vec.dy-dy*vec.dx);
    };
};

vecteur eX(1,0,0);
vecteur eY(0,1,0);
vecteur eZ(0,0,1);

// b) Classe "Point" customisée

class point{
    public:
    double x,y,z;
    int no, ktz;
    point(){};
    point(double i, double j, double k, int no_in, int ktz_in){
        x=i;y=j;z=k;no=no_in;ktz=ktz_in;
    }
    point(double i, double j, double k){
        x=i;y=j;z=k;no=-1;ktz=0;
    }
    vecteur vectorTo(point pt){
        return vecteur(pt.x-x,pt.y-y,pt.z-z);
    };
    vecteur vectorFrom(point pt){
        return vecteur(x-pt.x,y-pt.y,z-pt.z);
    };
    point interpolTo(point pt, double taux, int no_in, int ktz_in){
        return point(x+taux*(pt.x-x),y+taux*(pt.y-y),z+taux*(pt.z-z),no_in,ktz_in);
    }
};

typedef Triplet<double> T;

// III - Création des fonctions d'assemblage appropriées
// -----------------------------------------------------

Matrix3d m3_0;
Matrix3d m3_1;
Matrix3d z3 = Matrix3d::Constant(0);

// Cette fonction prend en entrée les points 1, 2 et 4 d'un élément
// Elle renvoie alors la matrice de passage de la base locale dans la base globale
Matrix3d local2global(point pt1,point pt2,point pt4){
    vecteur ex = (pt1.vectorTo(pt2)).normalize();
    vecteur ey = (pt1.vectorTo(pt4)).normalize();
    vecteur ez = ex.Pvector(ey);

    m3_0    << ex.Pscalar(eX), ey.Pscalar(eX), ez.Pscalar(eX),
               ex.Pscalar(eY), ey.Pscalar(eY), ez.Pscalar(eY),
               ex.Pscalar(eZ), ey.Pscalar(eZ), ez.Pscalar(eZ);
    return m3_0;
};

MatrixXd m24_0(24,24);
MatrixXd m24_1(24,24);
MatrixXd m24_2(24,24);

double tpsMultiplication=0.;
double tpsRemplissage_K=0.;
double tpsR3 = 0.;

// Cette fonction prend en entrée les quatre points de l'élément pt1 à pt4, le module d'Young E, le coefficient de Poisson nu et l'épaisseur t
// Elle ajoute les triplets de la matrice de raideur élémentaire dans la base globale dans la liste totale des triplets.
void global_stiffness_elem(vector<T> &lst, point pt1, point pt2, point pt3, point pt4, double E, double nu, double t){
    double B = E*t/(1-nu*nu);
    double D = E*t*t*t/(12*(1-nu*nu));
    double a = (pt1.vectorTo(pt2)).norm();
    double b = (pt1.vectorTo(pt4)).norm();
    int ktz1=pt1.ktz;int ktz2=pt2.ktz;int ktz3=pt3.ktz;int ktz4=pt4.ktz;

    auto start_m3 = high_resolution_clock::now();
    m3_1 = local2global(pt1, pt2, pt4);
    double r11=m3_1(0,0);double r12=m3_1(0,1);double r13=m3_1(0,2);
    double r21=m3_1(1,0);double r22=m3_1(1,1);double r23=m3_1(1,2);
    double r31=m3_1(2,0);double r32=m3_1(2,1);double r33=m3_1(2,2);
    auto stop_m3 = high_resolution_clock::now();
    auto duration_m3 = duration_cast<microseconds>(stop_m3 - start_m3);
    tpsR3 += duration_m3.count();

    auto start_mat = high_resolution_clock::now();
    m24_2 <<D*(r13*r13)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r11*(nu/8+1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r21*(nu/8+1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r31*(nu/8+1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r13*(D*r11*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),D*(r13*r13)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r11*(1/8-3*nu/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r21*(1/8-3*nu/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r31*(1/8-3*nu/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r13*(D*r11*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a))),D*(r13*r13)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r11*(-nu/8-1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r21*(-nu/8-1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r31*(-nu/8-1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r13*(D*r11*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)-b/(a*a))),D*(r13*r13)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r11*(3*nu/8-1/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r21*(3*nu/8-1/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r31*(3*nu/8-1/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r13*(D*r11*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a))),
            D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r11*(nu/8+1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r23*r23)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r21*(nu/8+1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r31*(nu/8+1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r23*(D*r11*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r11*(1/8-3*nu/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r23*r23)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r21*(1/8-3*nu/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r31*(1/8-3*nu/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r23*(D*r11*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a))),D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r11*(-nu/8-1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r23*r23)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r21*(-nu/8-1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r31*(-nu/8-1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r23*(D*r11*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)-b/(a*a))),D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r11*(3*nu/8-1/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r23*r23)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r21*(3*nu/8-1/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r31*(3*nu/8-1/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r23*(D*r11*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a))),
            D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r11*(nu/8+1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r21*(nu/8+1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r33*r33)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r31*(nu/8+1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r33*(D*r11*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r11*(1/8-3*nu/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r21*(1/8-3*nu/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r33*r33)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r31*(1/8-3*nu/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r33*(D*r11*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a))),D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r11*(-nu/8-1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r21*(-nu/8-1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r33*r33)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r31*(-nu/8-1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r33*(D*r11*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)-b/(a*a))),D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r11*(3*nu/8-1/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r21*(3*nu/8-1/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r33*r33)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r31*(3*nu/8-1/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r33*(D*r11*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a))),
            D*r11*r13*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r13*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r23*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r23*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r33*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r33*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),ktz1*(r13*r13)+r11*(-D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(-D*nu*r11+D*r12*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz1*r13*r23+r11*(-D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(-D*nu*r21+D*r22*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz1*r13*r33+r11*(-D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(-D*nu*r31+D*r32*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),D*r11*r13*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r13*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r23*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r23*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r33*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r33*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*(r11*r11)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r12*r12)*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r11*r13*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r13*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r23*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r23*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r33*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r33*(nu/(5*b)-1/(5*b)+b/(a*a)),D*(r11*r11)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r12*r12)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r13*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r13*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r23*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r23*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r33*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r33*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*(r11*r11)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r12*r12)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),
            D*r13*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r23*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r23*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r33*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r33*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),ktz1*r13*r23+r21*(-D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(-D*nu*r11+D*r12*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz1*(r23*r23)+r21*(-D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(-D*nu*r21+D*r22*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz1*r23*r33+r21*(-D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(-D*nu*r31+D*r32*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),D*r13*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r23*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r23*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r33*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r33*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r22*r22)*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r13*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r22*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r23*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r23*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r33*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r33*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r21*r21)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r22*r22)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r13*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r23*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r23*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r33*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r33*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r22*r22)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),
            D*r13*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r23*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r23*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r31*r33*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*r33*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),ktz1*r13*r33+r31*(-D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(-D*nu*r11+D*r12*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz1*r23*r33+r31*(-D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(-D*nu*r21+D*r22*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz1*(r33*r33)+r31*(-D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(-D*nu*r31+D*r32*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),D*r13*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r23*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r23*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r31*r33*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*r33*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r32*r32)*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r13*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r32*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r23*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r23*r32*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r31*r33*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*r33*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r31*r31)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r32*r32)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r13*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r23*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r23*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r31*r33*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*r33*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r32*r32)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),
            D*(r13*r13)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r11*(3*nu/8-1/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r21*(3*nu/8-1/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r31*(3*nu/8-1/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r13*(D*r11*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),D*(r13*r13)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r11*(-nu/8-1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r21*(-nu/8-1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r31*(-nu/8-1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r13*(D*r11*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),D*(r13*r13)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r11*(1/8-3*nu/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r21*(1/8-3*nu/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r31*(1/8-3*nu/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r13*(D*r11*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),D*(r13*r13)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r11*(nu/8+1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r21*(nu/8+1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r31*(nu/8+1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r13*(D*r11*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)+b/(a*a))),
            D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r11*(3*nu/8-1/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r23*r23)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r21*(3*nu/8-1/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r31*(3*nu/8-1/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r23*(D*r11*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r11*(-nu/8-1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r23*r23)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r21*(-nu/8-1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r31*(-nu/8-1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r23*(D*r11*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r11*(1/8-3*nu/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r23*r23)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r21*(1/8-3*nu/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r31*(1/8-3*nu/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r23*(D*r11*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r11*(nu/8+1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r23*r23)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r21*(nu/8+1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r31*(nu/8+1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r23*(D*r11*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)+b/(a*a))),
            D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r11*(3*nu/8-1/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r21*(3*nu/8-1/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r33*r33)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r31*(3*nu/8-1/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r33*(D*r11*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r11*(-nu/8-1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r21*(-nu/8-1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r33*r33)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r31*(-nu/8-1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r33*(D*r11*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r11*(1/8-3*nu/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r21*(1/8-3*nu/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r33*r33)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r31*(1/8-3*nu/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r33*(D*r11*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r11*(nu/8+1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r21*(nu/8+1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r33*r33)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r31*(nu/8+1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r33*(D*r11*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)+b/(a*a))),
            D*r11*r13*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r13*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r23*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r23*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r33*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r33*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*(r11*r11)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r12*r12)*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r11*r13*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r13*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r23*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r23*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r33*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r33*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),ktz2*(r13*r13)+r11*(D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(D*nu*r11+D*r12*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz2*r13*r23+r11*(D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(D*nu*r21+D*r22*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz2*r13*r33+r11*(D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(D*nu*r31+D*r32*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),D*r11*r13*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r13*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r23*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r23*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r33*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r33*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*(r11*r11)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r12*r12)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r13*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r13*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r23*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r23*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r33*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r33*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*(r11*r11)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r12*r12)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),
            D*r13*r21*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r23*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r23*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r33*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r33*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r22*r22)*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r13*r21*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r23*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r23*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r33*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r33*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),ktz2*r13*r23+r21*(D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(D*nu*r11+D*r12*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz2*(r23*r23)+r21*(D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(D*nu*r21+D*r22*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz2*r23*r33+r21*(D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(D*nu*r31+D*r32*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),D*r13*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r23*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r23*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r33*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r33*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r22*r22)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r22*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r23*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r23*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r33*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r33*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r21*r21)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r22*r22)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),
            D*r13*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r23*r31*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r23*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r31*r33*(a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*r33*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r32*r32)*((a*a)*nu-(a*a)+10*(b*b))/(15*a*b),D*r13*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r23*r31*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r23*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r31*r33*(2*a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*r33*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),ktz2*r13*r33+r31*(D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(D*nu*r11+D*r12*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz2*r23*r33+r31*(D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(D*nu*r21+D*r22*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),ktz2*(r33*r33)+r31*(D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(D*nu*r31+D*r32*(-4*(a*a)*nu+4*(a*a)+20*(b*b))/(15*a*b)),D*r13*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r23*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r23*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r31*r33*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*r33*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r32*r32)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r32*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r23*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r23*r32*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r31*r33*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*r33*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r31*r31)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r32*r32)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),
            D*(r13*r13)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r11*(-nu/8-1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r21*(-nu/8-1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r31*(-nu/8-1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r13*(D*r11*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)+b/(a*a))),D*(r13*r13)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r11*(3*nu/8-1/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r21*(3*nu/8-1/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r31*(3*nu/8-1/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r13*(D*r11*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r13*(D*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),D*(r13*r13)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r11*(nu/8+1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r21*(nu/8+1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r31*(nu/8+1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r13*(D*r11*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),D*(r13*r13)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r11*(1/8-3*nu/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r21*(1/8-3*nu/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r31*(1/8-3*nu/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r13*(D*r11*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r13*(D*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),
            D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r11*(-nu/8-1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r23*r23)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r21*(-nu/8-1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r31*(-nu/8-1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r23*(D*r11*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)+b/(a*a))),D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r11*(3*nu/8-1/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r23*r23)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r21*(3*nu/8-1/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r31*(3*nu/8-1/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r23*(D*r11*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r23*(D*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r11*(nu/8+1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r23*r23)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r21*(nu/8+1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r31*(nu/8+1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r23*(D*r11*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r11*(1/8-3*nu/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r23*r23)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r21*(1/8-3*nu/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r31*(1/8-3*nu/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r23*(D*r11*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r23*(D*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),
            D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r11*(-nu/8-1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r21*(-nu/8-1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r33*r33)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r31*(-nu/8-1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r33*(D*r11*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)+b/(a*a))),D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r11*(3*nu/8-1/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r21*(3*nu/8-1/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r33*r33)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r31*(3*nu/8-1/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r33*(D*r11*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),r33*(D*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a))),D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r11*(nu/8+1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r21*(nu/8+1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r33*r33)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r31*(nu/8+1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r33*(D*r11*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a))),D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r11*(1/8-3*nu/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r21*(1/8-3*nu/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r33*r33)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r31*(1/8-3*nu/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r33*(D*r11*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),r33*(D*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a))),
            D*r11*r13*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r13*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r23*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r23*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r33*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r33*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*(r11*r11)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r12*r12)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r13*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r13*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r23*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r23*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r33*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r33*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*(r11*r11)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r12*r12)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r13*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r13*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r23*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r23*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r33*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r33*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),ktz3*(r13*r13)+r11*(-D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(-D*nu*r11+D*r12*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz3*r13*r23+r11*(-D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(-D*nu*r21+D*r22*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz3*r13*r33+r11*(-D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(-D*nu*r31+D*r32*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),D*r11*r13*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r13*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r23*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r23*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r33*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r33*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*(r11*r11)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r12*r12)*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),
            D*r13*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r22*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r23*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r23*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r33*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r33*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r21*r21)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r22*r22)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r13*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r22*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r23*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r23*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r33*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r33*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r22*r22)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r22*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r23*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r23*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r33*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r33*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),ktz3*r13*r23+r21*(-D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(-D*nu*r11+D*r12*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz3*(r23*r23)+r21*(-D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(-D*nu*r21+D*r22*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz3*r23*r33+r21*(-D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(-D*nu*r31+D*r32*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),D*r13*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r23*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r23*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r33*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r33*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r22*r22)*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),
            D*r13*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r32*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r23*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r23*r32*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r31*r33*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*r33*(-nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r31*r31)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r32*r32)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r13*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r23*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r23*r32*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r31*r33*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*r33*(-4*nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r32*r32)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r23*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r23*r32*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r31*r33*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*r33*(4*nu/(5*b)+1/(5*b)+2*b/(a*a)),ktz3*r13*r33+r31*(-D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(-D*nu*r11+D*r12*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz3*r23*r33+r31*(-D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(-D*nu*r21+D*r22*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz3*(r33*r33)+r31*(-D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(-D*nu*r31+D*r32*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),D*r13*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r23*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r23*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r31*r33*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*r33*(nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r32*r32)*((a*a)*(nu-1)+10*(b*b))/(15*a*b),
            D*(r13*r13)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r11*(1/8-3*nu/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r21*(1/8-3*nu/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r11*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r12*(B*r31*(1/8-3*nu/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r13*(D*r11*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a))),D*(r13*r13)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r11*(nu/8+1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r21*(nu/8+1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r11*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r12*(B*r31*(nu/8+1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r13*(D*r11*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)-b/(a*a))),r13*(D*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)-b/(a*a))),D*(r13*r13)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r11*(3*nu/8-1/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r21*(3*nu/8-1/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r11*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r12*(B*r31*(3*nu/8-1/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r13*(D*r11*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a))),D*(r13*r13)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r11*(-nu/8-1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r21*(-nu/8-1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r11*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r12*(B*r31*(-nu/8-1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r13*(D*r11*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r13*(D*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),
            D*r13*r23*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r11*(1/8-3*nu/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r23*r23)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r21*(1/8-3*nu/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r21*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r22*(B*r31*(1/8-3*nu/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r23*(D*r11*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a))),D*r13*r23*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r11*(nu/8+1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r23*r23)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r21*(nu/8+1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r21*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r22*(B*r31*(nu/8+1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r23*(D*r11*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)-b/(a*a))),r23*(D*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)-b/(a*a))),D*r13*r23*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r11*(3*nu/8-1/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r23*r23)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r21*(3*nu/8-1/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r21*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r22*(B*r31*(3*nu/8-1/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r23*(D*r11*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a))),D*r13*r23*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r11*(-nu/8-1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r23*r23)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r21*(-nu/8-1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r21*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r22*(B*r31*(-nu/8-1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r23*(D*r11*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r23*(D*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),
            D*r13*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r12*(3*nu/8-1/8)+B*r11*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r11*(1/8-3*nu/8)+B*r12*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*r23*r33*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r22*(3*nu/8-1/8)+B*r21*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r21*(1/8-3*nu/8)+B*r22*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),D*(r33*r33)*(-4*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)+2*b/(a*a*a))+r31*(B*r32*(3*nu/8-1/8)+B*r31*((a*a)*(nu-1)+(b*b))/(6*a*b))+r32*(B*r31*(1/8-3*nu/8)+B*r32*(-4*(a*a)+(b*b)*(1-nu))/(12*a*b)),r33*(D*r11*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r21*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r31*(-2*a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a))),D*r13*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r12*(nu/8+1/8)+B*r11*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r11*(nu/8+1/8)+B*r12*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*r23*r33*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r22*(nu/8+1/8)+B*r21*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r21*(nu/8+1/8)+B*r22*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),D*(r33*r33)*(-2*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)-2*b/(a*a*a))+r31*(B*r32*(nu/8+1/8)+B*r31*((a*a)*(nu-1)-2*(b*b))/(12*a*b))+r32*(B*r31*(nu/8+1/8)+B*r32*(-2*(a*a)+(b*b)*(nu-1))/(12*a*b)),r33*(D*r11*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*(-nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r21*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*(-nu/(5*b)+1/(5*b)-b/(a*a))),r33*(D*r31*(-a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*(-nu/(5*b)+1/(5*b)-b/(a*a))),D*r13*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r12*(1/8-3*nu/8)+B*r11*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r11*(3*nu/8-1/8)+B*r12*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*r23*r33*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r22*(1/8-3*nu/8)+B*r21*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r21*(3*nu/8-1/8)+B*r22*((a*a)+(b*b)*(nu-1))/(6*a*b)),D*(r33*r33)*(2*a/(b*b*b)+4*nu/(5*a*b)-14/(5*a*b)-4*b/(a*a*a))+r31*(B*r32*(1/8-3*nu/8)+B*r31*((a*a)*(1-nu)-4*(b*b))/(12*a*b))+r32*(B*r31*(3*nu/8-1/8)+B*r32*((a*a)+(b*b)*(nu-1))/(6*a*b)),r33*(D*r11*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*(nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*(nu/(5*b)-1/(5*b)-2*b/(a*a))),D*r13*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r12*(-nu/8-1/8)+B*r11*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r11*(-nu/8-1/8)+B*r12*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*r23*r33*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r22*(-nu/8-1/8)+B*r21*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r21*(-nu/8-1/8)+B*r22*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),D*(r33*r33)*(4*a/(b*b*b)-4*nu/(5*a*b)+14/(5*a*b)+4*b/(a*a*a))+r31*(B*r32*(-nu/8-1/8)+B*r31*((a*a)*(1-nu)+2*(b*b))/(6*a*b))+r32*(B*r31*(-nu/8-1/8)+B*r32*(2*(a*a)+(b*b)*(1-nu))/(6*a*b)),r33*(D*r11*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),r33*(D*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a))),
            D*r11*r13*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r13*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r23*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r23*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r33*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r12*r33*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*(r11*r11)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r12*r12)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r13*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r13*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r23*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r23*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r33*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r12*r33*(nu/(5*b)-1/(5*b)+b/(a*a)),D*(r11*r11)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r12*r12)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r11*r13*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r13*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r23*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r23*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r33*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r12*r33*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*(r11*r11)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r12*r12)*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r11*r13*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r13*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r23*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r23*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r11*r33*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r12*r33*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),ktz4*(r13*r13)+r11*(D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(D*nu*r11+D*r12*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz4*r13*r23+r11*(D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(D*nu*r21+D*r22*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz4*r13*r33+r11*(D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r12*(D*nu*r31+D*r32*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),
            D*r13*r21*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r22*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r23*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r23*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r21*r33*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r22*r33*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r21*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r22*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r22*r22)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r21*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r22*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r23*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r23*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r21*r33*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r22*r33*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r21*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r22*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r21*r21)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r22*r22)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r13*r21*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r22*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r23*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r23*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r21*r33*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r22*r33*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r21*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r22*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r21*r21)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r22*r22)*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r21*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r22*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r23*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r23*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r21*r33*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r22*r33*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),ktz4*r13*r23+r21*(D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(D*nu*r11+D*r12*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz4*(r23*r23)+r21*(D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(D*nu*r21+D*r22*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz4*r23*r33+r21*(D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r22*(D*nu*r31+D*r32*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),
            D*r13*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r13*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r23*r31*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r23*r32*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r31*r33*(2*a/(b*b)-nu/(5*a)+1/(5*a))+D*r32*r33*(4*nu/(5*b)+1/(5*b)-b/(a*a)),D*r11*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r12*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*r22*r32*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+(b*b)*nu-(b*b))/(15*a*b)+D*(r32*r32)*(4*(a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r13*r32*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r23*r31*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r23*r32*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r31*r33*(a/(b*b)+nu/(5*a)-1/(5*a))+D*r32*r33*(nu/(5*b)-1/(5*b)+b/(a*a)),D*r11*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r12*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r21*r31*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*r22*r32*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*(r31*r31)*(5*(a*a)+(b*b)*(1-nu))/(15*a*b)+D*(r32*r32)*((a*a)*(1-nu)+5*(b*b))/(15*a*b),D*r13*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r13*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r23*r31*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r23*r32*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r31*r33*(-a/(b*b)+4*nu/(5*a)+1/(5*a))+D*r32*r33*(-nu/(5*b)+1/(5*b)+2*b/(a*a)),D*r11*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r12*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r21*r31*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*r22*r32*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*(r31*r31)*(10*(a*a)+4*(b*b)*nu-4*(b*b))/(15*a*b)+D*(r32*r32)*((a*a)*(nu-1)+10*(b*b))/(15*a*b),D*r13*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r13*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r23*r31*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r23*r32*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),D*r31*r33*(-2*a/(b*b)-4*nu/(5*a)-1/(5*a))+D*r32*r33*(-4*nu/(5*b)-1/(5*b)-2*b/(a*a)),ktz4*r13*r33+r31*(D*nu*r12+D*r11*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(D*nu*r11+D*r12*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz4*r23*r33+r31*(D*nu*r22+D*r21*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(D*nu*r21+D*r22*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b)),ktz4*(r33*r33)+r31*(D*nu*r32+D*r31*(20*(a*a)-4*(b*b)*nu+4*(b*b))/(15*a*b))+r32*(D*nu*r31+D*r32*(4*(a*a)*(1-nu)+20*(b*b))/(15*a*b));

    auto stop_mat = high_resolution_clock::now();
    auto duration_mat = duration_cast<microseconds>(stop_mat - start_mat);
    tpsRemplissage_K += duration_mat.count();

    
    auto start_mult = high_resolution_clock::now();

    auto stop_mult = high_resolution_clock::now();
    auto duration_mult = duration_cast<microseconds>(stop_mult - start_mult);
    tpsMultiplication += duration_mult.count();
    
    for (int i=0; i<24; i++){
        for (int j=0; j<=i; j++){
            double mij = m24_2(i,j);
            if ((mij>1e-3) || (mij<-1e-3)){
                int i_loc = i/6;
                int j_loc = j/6;
                int no_i, no_j;
                switch(i_loc) {
                    case 0: no_i = pt1.no; break;
                    case 1: no_i = pt2.no; break;
                    case 2: no_i = pt3.no; break;
                    case 3: no_i = pt4.no; break;
                };
                switch(j_loc) {
                    case 0: no_j = pt1.no; break;
                    case 1: no_j = pt2.no; break;
                    case 2: no_j = pt3.no; break;
                    case 3: no_j = pt4.no; break;
                };
                lst.push_back(T(6*no_i+i%6,6*no_j+j%6,mij));
                if (j<i){
                    lst.push_back(T(6*no_j+j%6,6*no_i+i%6,mij));
                };
                //cout << 6*no_i+i%6 << "," << 6*no_j+j%6 << "," << mij << endl;
            };
        };
    };
};

// Cette fonction prend en entrée les quatre points d'un élément et la force surfacique normale exercée sur celui-ci
// Elle complète la liste de triplets pour le vecteur-force.
void global_load_elem(vector<T> &lst, point pt1, point pt2, point pt3, point pt4, double gz){
    double a = (pt1.vectorTo(pt2)).norm();
    double b = (pt1.vectorTo(pt4)).norm();
    VectorXd F_loc(24);
    F_loc << 0,0,1,b/6.,-a/6.,0,0,0,1,b/6.,a/6.,0,0,0,1,-b/6.,a/6.,0,0,0,1,-b/6.,-a/6.,0;
    F_loc = a*b/4*gz*F_loc;
    m3_1 = local2global(pt1, pt2, pt4);
    m24_1 <<m3_1, z3, z3, z3, z3, z3, z3, z3,
            z3, m3_1, z3, z3, z3, z3, z3, z3,
            z3, z3, m3_1, z3, z3, z3, z3, z3,
            z3, z3, z3, m3_1, z3, z3, z3, z3,
            z3, z3, z3, z3, m3_1, z3, z3, z3,
            z3, z3, z3, z3, z3, m3_1, z3, z3,
            z3, z3, z3, z3, z3, z3, m3_1, z3,
            z3, z3, z3, z3, z3, z3, z3, m3_1;
    VectorXd F_glo(24);
    auto start_mult = high_resolution_clock::now();
    F_glo = m24_1*F_loc;
    auto stop_mult = high_resolution_clock::now();
    auto duration_mult = duration_cast<microseconds>(stop_mult - start_mult);
    tpsMultiplication += duration_mult.count();
    
    for (int i=0; i<24; i++){
        double fi = F_glo(i);
        if ((fi>1e-3) || (fi<-1e-3)){
            int i_loc = i/6;
            int no_i;
            switch(i_loc) {
                case 0: no_i = pt1.no; break;
                case 1: no_i = pt2.no; break;
                case 2: no_i = pt3.no; break;
                case 3: no_i = pt4.no; break;
            };
            lst.push_back(T(6*no_i+i%6,0,fi));
        };
    };
};

// IV - Création de la fonction de maillage
// ----------------------------------------

class carac_mail{
    public:
    int Nv, SNi;
    carac_mail(int Nv_in, int SNi_in){
        Nv = Nv_in;
        SNi = SNi_in;
    };
};

carac_mail mailleur(double lh, double lv, double H, vector<point> &lst, vector<point> lst_coins, int nb_coins){
    int compteur = 0;
    int SNi=0;
    int Nv = H/lv+1;

    for (int i=0; i<nb_coins;i++){
        int Ni = int((lst_coins[i].vectorTo(lst_coins[i+1])).norm()/lh)+1;
        SNi += Ni;
        for (int k=0; k<Ni; k++){
            double taux = double(k)/double(Ni);
            int ktz=1;
            if (k==0){
                ktz = 0;
            };
            lst.push_back(lst_coins[i].interpolTo(lst_coins[i+1],taux,compteur,ktz));
            compteur++;
        };
    };

    lst.reserve((Nv+1)*SNi);

    for (int k=1; k<Nv+1; k++){
        for (int j=0; j<SNi; j++){
            point pt_ori = lst[j];
            lst.push_back(point(pt_ori.x,pt_ori.y,pt_ori.z+H/double(Nv)*k,compteur,pt_ori.ktz));
            compteur += 1;
        };
    };
    return carac_mail(Nv, SNi);
};

// V - Création de la fonction de traitement total
// -----------------------------------------------

void batardeau(double lh, double lv, double H, vector<point> lst_coins, int nb_coins, double E, double nu, double t){
    vector<T> triplet_K;
    vector<T> triplet_F;
    vector<point> vertices;

    carac_mail kar = mailleur(lh, lv, H, vertices, lst_coins, nb_coins);
    int Nv=kar.Nv; int SNi=kar.SNi;
    int nb_total = SNi*(Nv+1);
    cout << "Nombre de nodes : " << nb_total << endl;

    triplet_K.reserve(36*4*nb_total+6*SNi);
    triplet_F.reserve(6*4*nb_total);

    auto start = high_resolution_clock::now();

    for (int k=0; k<Nv; k++){
        for (int l=0; l<SNi; l++){
            point pt1, pt2, pt3, pt4;
            if (l<SNi-1){
                pt1 = vertices[k*SNi+l];
                pt2 = vertices[k*SNi+l+1];
                pt3 = vertices[(k+1)*SNi+l+1];
                pt4 = vertices[(k+1)*SNi+l];
            };
            if (l==SNi-1) {
                pt1 = vertices[k*SNi+l];
                pt2 = vertices[k*SNi];
                pt3 = vertices[(k+1)*SNi];
                pt4 = vertices[(k+1)*SNi+l];
            };

            global_stiffness_elem(triplet_K,pt1,pt2,pt3,pt4,E,nu,t);
            global_load_elem(triplet_F,pt1,pt2,pt3,pt4,10e3*(H-pt1.z));
        };
    };

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Durée remplissage [s] : " << duration.count()/1000000.0 << endl;

    double maxK = 1e20;
    for (int k=0; k<6*SNi; k++){
        triplet_K.push_back(T(k,k,maxK));
    };

    SparseMatrix<double> mat_K(6*nb_total,6*nb_total);
    mat_K.setFromTriplets(triplet_K.begin(),triplet_K.end());

    SparseMatrix<double> mat_F(6*nb_total,1);
    mat_F.setFromTriplets(triplet_F.begin(),triplet_F.end());

    SimplicialCholesky<SparseMatrix<double>> chol(mat_K);
    VectorXd q_sol = chol.solve(mat_F);

    //ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
    //SimplicialLDLT<SparseMatrix<double>> solver;
    //SimplicialLLT<SparseMatrix<double>> solver;
    //solver.compute(mat_K);
    //VectorXd q_sol = solver.solve(mat_F);

    ofstream MyFile("output_batardeau.csv");

    MyFile << "x,y,z,ux,uy,uz" <<endl;
    point stockp;
    for (int k=0; k<nb_total; k++){
        stockp = vertices[k];
        MyFile << stockp.x << "," << stockp.y << "," << stockp.z << "," << q_sol(6*k) << "," << q_sol(6*k+1) << "," << q_sol(6*k+2) << endl;
    };

    MyFile.close();
};

// VI - Fonction main, script du programme
// ---------------------------------------

int main() {

    auto start2 = high_resolution_clock::now();
    double lh, lv, H, E, nu, t;
    int nb_coins;
    string extract;
    vector<point> lst_coins;

    ifstream rFile("input_batardeau.txt");
    
    rFile >> lh >> lv >> H >> E >> nu >> t >> nb_coins;
    
    for (int k=0; k<nb_coins; k++){
        double x0,y0;
        rFile >> x0 >> y0;
        lst_coins.push_back(point(x0,y0,0));
    };
    lst_coins.push_back(lst_coins[0]);


    rFile.close();

    batardeau(lh,lv,H,lst_coins,nb_coins, E, nu, t);

    cout << "Durée création K [s] : " << tpsRemplissage_K/1000000.0 << endl;
    cout << "Durée création R [s] : " << tpsR3/1000000.0 << endl;
    cout << "Durée multiplication [s] : " << tpsMultiplication/1000000.0 << endl;

    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>(stop2 - start2);
    cout << "Durée totale [s] : " << duration2.count()/1000000.0 << endl;

    return 0;
};
