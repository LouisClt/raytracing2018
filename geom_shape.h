
#ifndef GEOM_SHAPE_HEADER
#define GEOM_SHAPE_HEADER

#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <utility>
#include <vector>
#include <iostream>
#include <vector>
#include <float.h>
#include <float.h>
#include <random>





class Vector3D
{
private:
    double V[3];
public:
    Vector3D(double a, double b, double c);
    Vector3D();
    /// Bracket operator.
     double& operator[](unsigned n);
     Vector3D operator+(Vector3D vect1) const ;
     Vector3D(  const Vector3D& vect1)  ;
     Vector3D operator*(  double a) const  ;
     double operator*(  Vector3D a) const  ;
     Vector3D operator-( Vector3D vect1) const;
    Vector3D p_term( Vector3D vect1) const  ;
    Vector3D p_vect(Vector3D vect1);
    double max();
    bool operator!=(Vector3D vect);
    bool operator==(Vector3D vect);
    int min();
    int arg_max()
    {
        double max_abs=0.;
        int ind_max=0;
        for (int i =0;i<3;i++)
        {
            if (abs(V[i])>max_abs)
            {
                max_abs= abs(V[i]);
                ind_max =i;
            }

        }
        return(ind_max);
    }

    void normalize();
};


// Forward declare 
class Ray;
class Scene;

class Object
{
    friend Scene;
    //friend Sphere3D;
    //friend Triangle;

    protected: 
        Vector3D Color_abs;
        Vector3D Emissivity;
        double N_opt;
        bool Is_light;
        bool Is_specular;
        bool Is_transparent;
    public: 
        Object(): Color_abs(Vector3D(1,0,0)),N_opt(1.),Is_light(false),Is_specular(false),Is_transparent(false){};
        Object( Vector3D color_abs, Vector3D emissivity,double n_opt, bool is_specular,bool is_transparent): Color_abs(color_abs), Emissivity(emissivity),N_opt(n_opt),
        Is_specular(is_specular),Is_transparent(is_transparent){};

        virtual ~Object(){};
        virtual bool intersect( Ray ray, double& t, Object* &object)=0;
        //virtual Vector3D pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count)=0;
        virtual Vector3D n(Vector3D X )=0; 
        virtual Vector3D pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count);
        
        virtual Vector3D getAlbedo(Vector3D n, Vector3D p)
        {
            Vector3D color1 (0,0,0);
            //Vector3D color = n.p_term(Vector3D(1,1,1))*(P*P);
            //color= color*(1/sqrt(color*color));
            
            Vector3D offset (0,0,0);
         
            
                color1[0]=(cos(p[0]/5));//*abs(sin(n[i]+offset[i]));
                color1[1]=(cos(p[1]/5));
                color1[2]=(cos(p[0]/10)*sin(p[1]/10));

            color1= n.p_term(color1);
                for(int i=0; i<3; i++)
                {
                    color1[i]= fabs(color1[i]);
                }
            
            return(color1);
        };



};




class Geometry;

class Triangle: public Object 
{
    friend Scene;
    private:
        Vector3D A, B,C,Na,Nb,Nc,UVa,UVb,UVc, N;
        Geometry* Geom_pt;
        int FaceGroup, w,h;


    public:
        static long unsigned  COUNT_OF_INTERSECT;
        virtual ~Triangle(){}; 
        Triangle();
        Triangle(Vector3D a, Vector3D b, Vector3D c): A(a),B(b),C(c)
        {
            Vector3D CA =(A-C);
            CA.normalize();
            Vector3D CB= B-C;
            CB.normalize();
            N= CA.p_vect(CB);
            N.normalize();
        };


        Triangle(Vector3D a, Vector3D b, Vector3D c, 
                Vector3D na, Vector3D nb,Vector3D nc,
                Vector3D uva, Vector3D uvb, Vector3D uvc,
                Geometry* geom_pt,
                int faceGroup, int w1, int h1): A(a),B(b),C(c),Na(na),Nb(nb),Nc(nc),Geom_pt(geom_pt),
                    UVa(uva),UVb(uvb), UVc(uvc),FaceGroup(faceGroup),w(w1),h(h1)
        {
            Vector3D CA =(A-C);
            CA.normalize();
            Vector3D CB= B-C;
            CB.normalize();
            N= CA.p_vect(CB);
            N.normalize();
        };


        // const  ..& ??
         bool intersect( Ray ray, double& t, Object* &object);
        //Vector3D pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count);
        Vector3D n(Vector3D X =Vector3D(0,0,0));
        Vector3D& operator[](unsigned n)
        {
            if (n==0)
            {
                return(A);
            }
            if (n==1)
            {
                return(B);
            }
            if (n==2)
            {
                return(C);
            }
            else 
            {
                return(N);
            }
            
        };

        Vector3D mean()
        {
            return((A+B+C)*(1./3.));
        }
        
        virtual Vector3D getAlbedo(Vector3D n, Vector3D P);

};


class Sphere3D: public Object
{
    friend Scene;
    private:
        Vector3D Center;
        

        double Radius;
        

        
        
    public: 
        Sphere3D();
        Sphere3D( Vector3D center, double radius);
        Sphere3D( Vector3D center, double radius, Vector3D color_abs, Vector3D emissivity);
        Sphere3D( Vector3D center, double radius, 
        Vector3D color_abs, Vector3D emissivity, bool is_specular,bool is_transparent);
        ~Sphere3D(){};
        double radius() const;
        Vector3D center() const ;
        bool intersect( Ray ray, double& t, Object* &object);
        
        Vector3D n(Vector3D X );

};

/*
class TriVect: public Object 
{

     protected:
        std::vector <Triangle> Tri_Vect;

   public:

    TriVect(const char* obj, double scaling=0., const Vector3D& offset=Vector3D(0,0,0));
    )virtual bool intersect( Ray ray, double& t)=0;
    //virtual Vector3D pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count)=0;
    virtual Vector3D n(Vector3D X =Vector3D(0,0,0))=0; 
    virtual Vector3D pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count);

};*/
class Ray
{
    private:
    Vector3D Center;
    Vector3D U;
    public: 
    Ray();
    Ray( Vector3D center,  Vector3D u);
     Vector3D center();
     Vector3D u();
     void normalize();
};

class Geometry; 

class Scene 
{
    private:
     std::vector<Object*> Obj_vect;
    int Light_number;
    
    public: 
        static int MAX_RBDS_CNT;
        Scene();
        ~Scene()
        {
            for(unsigned i=0; i<Obj_vect.size();i++)
            {
                
                //delete Obj_vect[i];
                //Obj_vect[i]=0;
            }
        };

        void add_sphere(Sphere3D sphere);
        void add_triangle(Triangle triangle);
        void add_geom(Geometry* geom1);
        //void add_lum(Sphere3D lum);
        bool intersect( Ray ray,double& t); 
        Vector3D pix_color (Ray ray, int rebonds_count, bool direct_light_op);
        Object* get_direct_light (Vector3D P);

};

#endif