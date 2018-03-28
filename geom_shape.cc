
#include "geom_shape.h"
#include <math.h>
#include <limits.h>
#include <float.h>
#include <random>
#include "geometry.h"
#include <typeinfo>
//std::default_random_engine engine;
//std::uniform_real_distribution <double> u(0,1);

std::default_random_engine engine;
std::uniform_real_distribution <double> u(0,1);

double PI = 3.14159265;
//==========================================
/// 
///                Class Vector3D
///         
//==========================================




Vector3D::Vector3D()
{
    V[0]=0.;
    V[1]=0.;
    V[2]=0.;
}


Vector3D::Vector3D(double a, double b, double c)
        
{
    V[0]=a;
    V[1]=b;
    V[2]=c;
}




/// Bracket operator.
double& Vector3D::operator[](unsigned n)
{
    if ((n<0)||(n>2))
    {
        throw ;
    }
    else
    {
        return(V[n]);
    }
}




 Vector3D Vector3D::operator+( Vector3D vect1) const
{
    return(Vector3D(vect1[0]+V[0],
    vect1[1]+V[1],
    vect1[2]+V[2]));
}

 Vector3D Vector3D::operator*( double a) const
{
    return(Vector3D(V[0]*a,
    V[1]*a,
    V[2]*a));
}



 Vector3D Vector3D::operator-( Vector3D vect1) const 
{
    return(Vector3D(V[0]-vect1[0],
    V[1]-vect1[1],
    V[2]-vect1[2]));
}

double Vector3D::operator*(Vector3D vect1) const 
{
    return(vect1[0]*V[0]+
    vect1[1]*V[1]+
    vect1[2]*V[2]);
}
Vector3D Vector3D::p_term( Vector3D vect1) const  
{
    return(Vector3D(vect1[0]*V[0],
    vect1[1]*V[1],
    vect1[2]*V[2]));
}
 void Vector3D::normalize()
 {
    double a=  sqrt((1./((*this) * (*this))));
    for (unsigned i=0;i<3;i++)
    {
        V[i]=a*V[i];
    }
 }

Vector3D Vector3D::p_vect(Vector3D vect1)
{
    return( Vector3D(
        V[1]*vect1[2]-vect1[1]*V[2],
        -V[0]*vect1[2]+vect1[0]*V[2],
        V[0]*vect1[1]-vect1[0]*V[1]));

}

int Vector3D::min()
{
    double  a= abs(V[0]);
    int i =0;
    if (abs(V[1])<a)
    {
        i=1;
        a = abs(V[1]);
    }
    if (abs(V[2])<a)
    {
        i=2;
        //a = V[1]<a;
    } 
    return(i);
}

Vector3D::Vector3D(  const Vector3D& vect1) 
{
    V[0] = vect1.V[0];
    V[1] = vect1.V[1];
    V[2] = vect1.V[2];
}

bool Vector3D::operator!=(Vector3D vect)
{
    for (int i=0; i<3;i++)
    {
        if (V[i]!=vect[i])
        {
            return(true);
        }
    }
    return(false);
}

bool Vector3D::operator==(Vector3D vect)
{
    for (int i=0; i<3;i++)
    {
        if (V[i]!=vect[i])
        {
            return(false);
        }
    }
    return(true);
}
double Vector3D::max()
{
    double a=0.;
    for (int i=0;i<3;i++)
    {
        if (abs(V[i])>a)
        {
            a=abs(V[i]);
        }
    }
    return(a);
};
//==============================================================================
///
///                             Class Light
/// 
//==============================================================================
/*Light::Light(Vector3D center, double int_dbl)
{
    Vector3D temp_vect(1.,1.,1.);
    Center=center;
    I=temp_vect*int_dbl;
}


Light::Light(Vector3D center, Vector3D i): Center(center),I(i)
{
    //Center(center);
    //I(i);
}
Vector3D Light::center()
{
    return(Center);
}
Vector3D Light::i()
{
    return(I);
}
void Light::light_abs(Vector3D albedo)
{
    I[0] = I[0] * albedo[0];
    I[1] = I[1] * albedo[1];
    I[2] = I[2] * albedo[2];

}
*/


//===============================================================
///
///                         Class Object
///
//===============================================================

// 

Vector3D  Object::pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count)
{
    if (Is_light)
    {
        Vector3D color = Emissivity;
        //double a=color.max();
        //color=color*(1/a);
        return(color);
    }

    
    
    Vector3D n= this->n(P);
    Vector3D i = ray.u();
    double norm_squared_n= n*n;
    n = n*(1/sqrt(norm_squared_n));
    double i_scal_n= i*n;
    double Eps= pow(10,-9);
    
    if (Is_specular)
    {
        if (rebonds_count==0)
        {
            return(Vector3D());
        }
        else
        {
            // creation nouveau rayon
            Vector3D c_corr= P + n*Eps;
            Vector3D i=ray.u();
            i.normalize();
            Vector3D direct =  i - n*(2* (i*n));
            Ray ray1(c_corr,direct);
            ray1.normalize();
            return(Color_abs.p_term(scene->pix_color(ray1, rebonds_count-1,false)));
        }

        
    }
    else if(Is_transparent) 
    {
        if (n*i<0.)// o nentre dans la sphère 
        {
            double temp_var = 1- (1/N_opt)*(1/N_opt)*
                            (1-(i_scal_n)*(i_scal_n));
            if (temp_var>=0)
            {
                Vector3D temp_0;
                Vector3D Tn= temp_0-n*sqrt(temp_var);
                Vector3D Tt= (i-n*(i_scal_n))*(1/N_opt);// 1 pour nair
                Vector3D T= Tn+Tt;
                Vector3D P_corr = P-n*Eps;
                Ray ray1(P_corr,T);
                double t=DBL_MAX;
                Object* obj = 0;
                intersect(ray1,t,obj);
                P_corr= P_corr + T*t;
                Ray ray2 (P_corr,T);
                return(pix_color(ray2,P_corr,scene,rebonds_count));
            }
            // suite impossible
            /*else // Si on est réfléchi 
            {
                // ON reprend les mêmes lignes que précédemment 
                // revoir ..???
                // creation nouveau rayon
                Vector3D c_corr= P + n*Eps;
                Vector3D i=ray.u();
                i.normalize();
                Vector3D direct =  i - n*(2* (i*n));
                Ray ray1(c_corr,direct);
                ray1.normalize();
                scene->pix_color(lum,ray1,colors, rebonds_count-1);
            }*/
        }
        else // ON sort de la sphère 
        {
            double temp_var = 1- (N_opt)*(N_opt)*
                            (1-(i_scal_n)*(i_scal_n));
            if (temp_var>=0)
            {
                Vector3D Tn= n*sqrt(temp_var);
                Vector3D Tt= (i+n*(i_scal_n))*(N_opt);// 1 pour nair
                Vector3D T= Tn+Tt;
                Vector3D P_corr = P+n*Eps;
                Ray ray1(P_corr,T);
                return(scene->pix_color(ray1,rebonds_count,false));
            }
            else
            {
                Vector3D c_corr= P - n*Eps;
                Vector3D direct =  i - n*(2* (i*n));
                Ray ray1(c_corr,direct);
                double t=DBL_MAX;
                Object* obj=0;
                intersect(ray1,t, obj);
                c_corr= c_corr + direct*t;
                Ray ray2 (c_corr,direct);
                return(pix_color(ray2,c_corr,scene,rebonds_count-1));

                //scene->pix_color(lum,ray1,colors, rebonds_count-1);
                // ne rien faire on fait comme si on restait dans la sphère 
            }
            
        }
    }
    else
    {
        // sphere normale
        Vector3D colors (0.,0.,0.);

        Sphere3D* lum= dynamic_cast<Sphere3D*> (scene->get_direct_light(P));
        Vector3D I=lum->Emissivity;// pour l'instant on définit ici l'intensité...
        Vector3D l=P-lum->center();
        
        double norm_squared_l= l*l;
        l.normalize();
        //Vector3D c_l_corr= P + n*Eps;
        
        int m = l.min();
        Vector3D b1;
        if (m==0)
        {
            b1 = Vector3D(0,-l[2],l[1]);
            b1.normalize();
        }
        else if (m==1)
        {
            b1 = Vector3D(-l[2],0,l[0]);
            b1.normalize();
        }
        else 
        {
            b1 = Vector3D(l[1],-l[0],0);
            b1.normalize();
        }


        //std::random_device rd;
        //std::mt19937 mt(rd());
        //std::uniform_real_distribution<double> u(1.0, 10.0);


        Vector3D c1= l.p_vect(b1);
        double r1 = u(engine);
        double r2 = u(engine);
        double x = cos(2*PI*r1)*sqrt(1-r2);
        double y = sin(2*PI*r1)*sqrt(1-r2);
        double z = sqrt(r2);
        Vector3D new_vect = b1*x+ c1*y + l*z;
        Vector3D new_P= new_vect*(lum->radius()+Eps)+lum->center();
        if ( typeid(*this) == typeid( Triangle ))
        {
            if ((i*n)>=0)
            {
                n = Vector3D() -n*(1+Eps);
            }
        }
        //Vector3D color_abs1=Color_abs;
        Vector3D color_abs1= getAlbedo(n,P);

        Vector3D new_vect2 = new_P-P-n*Eps;


        Ray new_ray (P+n*Eps,new_vect2 );
 

        double t=DBL_MAX;
        //bool intersect_only_with_not_transp=-1.;
        
        if (scene->intersect(new_ray,t)&&(t<=1.))
        {
          
        } 
        else 
        {
            Vector3D O_newS=new_P-lum->center();
            Vector3D n_prime= (O_newS)*(1/lum->radius());
            Vector3D OP= (P-lum->center());
            OP.normalize();
            Vector3D PprimeP= (P-new_P);
            double norm_sq_PPprime= PprimeP*PprimeP;
            double one_on_norm= 1/sqrt(norm_sq_PPprime);
            PprimeP= PprimeP*(one_on_norm);
            Vector3D PPprime= (new_P-P);
            PPprime=PPprime*(one_on_norm);
            double cos_theta= (n*PPprime);
            double cos_theta_prime=(n_prime*PprimeP);
            double cos_theta_2prime=n_prime*OP;
            //if (cos_theta<0||cos_theta_prime<0||cos_theta_2prime<0)
            //{
            //    std::cout<<cos_theta<<" "<<cos_theta_prime<<" 
            // "<<cos_theta_2prime<<std::endl;
            //}
            
          
        
            if ( typeid(*this) == typeid( Triangle ))
            {
                //std::cout<<" ok enter"<<std::endl;
                 //cos_theta= abs(n*PPprime);
                //double cos_theta_prime=(n_prime*PprimeP);
                //cos_theta_prime=abs(n_prime*PprimeP);
            }

            double pix_intensity= lum->radius()*lum->radius()*std::max(cos_theta_prime,0.)/cos_theta_2prime*std::max(cos_theta,0.)/(norm_sq_PPprime);// normalkement PI/PI
            for (unsigned i=0; i<3;i++)
            {
                colors[i]=colors[i]+pix_intensity*color_abs1[i]*I[i];
            }
        }
           
        //return(colors);        
        m = n.min();
       
        if (m==0)
        {
            b1 = Vector3D(0,-n[2],n[1]);
            b1.normalize();
        }
        else if (m==1)
        {
            b1 = Vector3D(-n[2],0,n[0]);
            b1.normalize();
        }
        else 
        {
            b1 = Vector3D(n[1],-n[0],0);
            b1.normalize();
        }


        //std::random_device rd;
        //std::mt19937 mt(rd());
        //std::uniform_real_distribution<double> u(1.0, 10.0);


         c1= n.p_vect(b1);
         r1 = u(engine);
        r2 = u(engine);
        x = cos(2*PI*r1)*sqrt(1-r2);
        y = sin(2*PI*r1)*sqrt(1-r2);
        z = sqrt(r2);
        new_vect = b1*x+ c1*y + n*z;
        Ray new_ray4 (P+n*Eps,new_vect );
        //lum.light_abs(Color_abs);
        return(colors+ color_abs1.p_term(scene->pix_color(new_ray4,rebonds_count-1,true)));
    }
}



//==================================================================
/// 
///                     Class Triangle 
///
//==================================================================
long unsigned  Triangle::COUNT_OF_INTERSECT=0;
bool Triangle::intersect( Ray ray, double& t, Object* &object)
{

    Triangle::COUNT_OF_INTERSECT ++;

    // intersection rayon-plan
     t = ((A-ray.center())*N)/(ray.u()*N);
     object = this;
    if (t<0.)
    {
        return(false);
    }
    else
    {
        Vector3D P=ray.center()+ray.u()*t;//intersection
        double APAB = (P-A)*(B-A);
        double APAC = (P-A)*(C-A);
        double ACAC = (C-A)*(C-A);
        double ACAB = (C-A)*(B-A);
        //double APAC = (P-A)*(C-A);
        double ABAB = (B-A)*(B-A);
        
        double det = ABAB*ACAC-ACAB*ACAB;
        double Beta = (APAB*ACAC-APAC*ACAB)/det;
        if (Beta<0)
        {
            return(false);
        }
        else 
        {
            double Gamma = (ABAB*APAC-APAB*ACAB)/det;
            if (Gamma<0.)
            {
                return(false);
            }
            else
            {
                if (Gamma+Beta>1.)
                {
                    return(false);
                }
                else
                {
                    return(true);
                }
            }
        }


    }
}



/*Vector3D Triangle::pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count)
{
    return(Vector3D());
}*/

Vector3D Triangle::n(Vector3D P )
{

    double APAB = (P-A)*(B-A);
    double APAC = (P-A)*(C-A);
    double ACAC = (C-A)*(C-A);
    double ACAB = (C-A)*(B-A);
    //double APAC = (P-A)*(C-A);
    double ABAB = (B-A)*(B-A);
    
    double det = ABAB*ACAC-ACAB*ACAB;
    double Beta = (APAB*ACAC-APAC*ACAB)/det;
    double Gamma = (ABAB*APAC-APAB*ACAB)/det;
    double alpha = 1. -Beta - Gamma;

    
    /*double normerr = sqrt((P-A*alpha - B*Beta - C*Gamma)*(P-A*alpha - B*Beta - C*Gamma));
    if (normerr>0.0001)
    {
        std::cout<<" WTF"<<normerr<<std::endl;
    }//*/
    //Vector3D color = n.p_term(Vector3D(1,1,1))*(P*P);
    //color= color*(1/sqrt(color*color));
    
    //Vector3D offset (0,0,0);
    Vector3D n1 = Na*alpha + Nb*Beta + Nc*Gamma;
    return(n1);
}
Vector3D Triangle::getAlbedo(Vector3D n, Vector3D P)
{

    double APAB = (P-A)*(B-A);
    double APAC = (P-A)*(C-A);
    double ACAC = (C-A)*(C-A);
    double ACAB = (C-A)*(B-A);
    //double APAC = (P-A)*(C-A);
    double ABAB = (B-A)*(B-A);
    
    double det = ABAB*ACAC-ACAB*ACAB;
    double Beta = (APAB*ACAC-APAC*ACAB)/det;
    double Gamma = (ABAB*APAC-APAB*ACAB)/det;
    double alpha = 1. -Beta - Gamma;

    
    /*double normerr = sqrt((P-A*alpha - B*Beta - C*Gamma)*(P-A*alpha - B*Beta - C*Gamma));
    if (normerr>0.0001)
    {
        std::cout<<" WTF"<<normerr<<std::endl;
    }//*/
    //Vector3D color = n.p_term(Vector3D(1,1,1))*(P*P);
    //color= color*(1/sqrt(color*color));
    
    //Vector3D offset (0,0,0);
    Vector3D UV = UVa*alpha + UVb*Beta + UVc*Gamma;

    //int w= (*Textures_pt)[FaceGroup][0].size();
    //int h= (*Textures_pt)[FaceGroup].size();
    Vector3D color1 (0,0,0);
    int i = UV[1]*(h-1);
    int j = UV[0]*(w-1);
    for (int k=0; k<3; k++)
    {
        color1[k] = double((Geom_pt->textures)[FaceGroup][ (i*w+j)*3+k])/255;
    }
    
    
    //color1[0]=(cos(P[0]/5));//*abs(sin(n[i]+offset[i]));
    //color1[1]=(sin(P[1]/5));
    //color1[2]=(cos(P[2]/10)*sin(P[1]/10));

    //color1= n.p_term(color1);
        //for(int i=0; i<3; i++)
        //{
            //   color1[i]= fabs(color1[i]);
        // }
    
    return(color1);
}
//=========================================================
///
///             Class Sphere3D
///
//=========================================================
double N_SPHERE=1.2;
Sphere3D::Sphere3D():  Object(Vector3D(0.95,0.95,0.95),Vector3D(),
            N_SPHERE,false,false), Center(),Radius(1.)
{
     Is_light= false;
}
Sphere3D::Sphere3D( Vector3D center, double radius):
    Object(Vector3D(0.95,0.95,0.95),Vector3D(),N_SPHERE,false,false),
    Center(center),
    Radius(radius)
{
    Is_light= false;
}
Sphere3D::Sphere3D( Vector3D center, double radius, Vector3D color_abs, Vector3D emissivity):
        Object(color_abs,emissivity,N_SPHERE,false,false),
        Center(center),Radius(radius)
{
    if (Emissivity==Vector3D())
    {
         Is_light= false;
    }
    else 
    {
        Is_light=true;
    }
}

Sphere3D::Sphere3D( Vector3D center, double radius, 
        Vector3D color_abs, Vector3D emissivity, bool is_specular, bool is_transparent):
        Object(color_abs,emissivity,N_SPHERE,
            is_specular,is_transparent),
        Center(center),Radius(radius)
{
    if (Emissivity==Vector3D())
    {
        Is_light=false;
    }
    else 
    {
        Is_light=true;
    }
}

double Sphere3D::radius() const
{
    return(Radius);
}
Vector3D Sphere3D::center() const 
{
    return(Center);
}

bool Sphere3D::intersect( Ray ray, double& t, Object* &object)
{
    object =this;
    Vector3D c=ray.center();
    Vector3D u=ray.u();

    double a =(u*u);
    double b = 2*(u*(c-Center));
    double c1 = (c-Center)*(c-Center)-
            Radius*Radius;
    double delta1 = b*b - 4*a*c1;

    if (delta1>=0.)
    {
        //std::cout<<delta<<std::endl;
        double t1=(-b+sqrt(delta1))/(2*a);
        double t2=(-b-sqrt(delta1))/(2*a);
        if (std::max(t1,t2)>=0.)
        {
            if(t1>t2)
            {
                if (t2>=0)
                {
                    t=(t2);
                }
                else
                {
                    t=(t1);
                }
            }
            else
            {
                if (t1>=0)
                {
                    t=(t1);
                }
                else
                {
                    t=(t2);
                }
            }
            return(true);
        }
    }
    t=DBL_MAX;
    return(false);
}



Vector3D Sphere3D::n(Vector3D X )
{
    Vector3D n=(X-Center);
    n.normalize();
    return(n);
}
//=====================================================
///
///             Class Ray
///             
//====================================================

Ray::Ray( Vector3D center,  Vector3D u):
    Center(center),U(u)
{

}

Vector3D Ray::center()
{
    return(Center);
}
Vector3D Ray::u()
{
    return(U);
}

 void Ray::normalize()
 {
     U=U*(1/sqrt(U*U));
 }


//==============================================================================
///
///                         Class Scene
///
//==============================================================================
int Scene::MAX_RBDS_CNT=15;

Scene::Scene()
{
    Light_number=0;
}

void Scene::add_sphere(Sphere3D sphere)
{
    if (sphere.Emissivity!=Vector3D())
    {
        std::vector<Object*>::iterator it;
        it = Obj_vect.begin();
        Sphere3D* sphere_pt = new Sphere3D (sphere);
        Obj_vect.insert(it,sphere_pt);
        Light_number++;
    }
    else 
    {
        Sphere3D* sphere_pt = new Sphere3D (sphere);
        Obj_vect.push_back(sphere_pt);
    }
}

void Scene::add_triangle(Triangle triangle)
{   
    Triangle* triangle_pt = new Triangle (triangle);
    Obj_vect.push_back(triangle_pt);
}
 
void Scene::add_geom(Geometry* geom1)
{
    
    //Geometry* geom_pt = new Geometry (*geom1);
    Obj_vect.push_back(geom1);
}
bool Scene::intersect( Ray ray, double& t)
{
    t = DBL_MAX;
    //unsigned closest_i=0;
    unsigned n = Obj_vect.size();
    for (unsigned i=0; i<n;i++)
    {
        double t_temp(DBL_MAX);
        Object* obj=0;
        bool intersection=Obj_vect[i]->intersect(ray,t_temp,obj);
        if (intersection&&t_temp<t)
        {
            t=t_temp;
            //closest_i=i;
        }
    }
    if (t==DBL_MAX)
    {
        return (false);
    }
    return(true);

}

Vector3D Scene::pix_color (Ray ray, int rebonds_count, bool indirect_light_op)
{
    double t = DBL_MAX;
    unsigned closest_i=0;
    unsigned n = Obj_vect.size();
    Object* closest_obj= 0;
    //int rays_generated =15;
    for (unsigned i=0; i<n;i++)
    {
        double t_temp(DBL_MAX);
        Object* obj=0;
        bool intersection=Obj_vect[i]->intersect(ray,t_temp,obj);
        if (intersection&&t_temp<t)
        {
            t=t_temp;
            closest_obj = obj;
        }
    }
    if (t==DBL_MAX)
    {
        
        return(Vector3D());
    }
    // else ....
    /*if genere_rays
    {
        Vector3D P= ray.center()+ray.u()*t;
        Sphere_vect[closest_i].pix_color(lum,ray,P,colors, this, rebonds_count);
        int N=15;
        for (int i=0; i<N; i++)
        {
            int m= P- Sphere_vect[closest_i]
        }
    }*/

    Vector3D P= ray.center()+ray.u()*t;
    if ((rebonds_count>=0) &&(!(indirect_light_op&&(Obj_vect[closest_i]->Is_light))))
    {
        return(closest_obj->pix_color(ray,P, this, rebonds_count));
    }
    else
    {
        return(Vector3D(0.,0.,0.));
    }      
}

Object* Scene::get_direct_light (Vector3D P)
{
    return(Obj_vect[rand()%Light_number]);
    // D'autres méthodes sont 
    // possibles sur comment échantillonnes les lumières pour réduire le bruit
    // (ie si une sphere de lumière est toute petite etc )
}








































BVH::BVH(int i_debut, int i_fin , Geometry* mother_g ):Mother_geom(mother_g), I_debut(i_debut), I_fin(i_fin)
{
    Vector3D max_vertex (0,0,0);
    Vector3D min_vertex (DBL_MAX,DBL_MAX,DBL_MAX);
    for (int i=I_debut; i<I_fin;i++)
    {
        for (int k=0; k <3; k++)
        {
            for (int m = 0; m<3;m++)
            {
                if ((*Mother_geom)[i][k][m]<min_vertex[m])
                {
                    min_vertex[m]=(*Mother_geom)[i][k][m];
                }
                if ((*Mother_geom)[i][k][m]>max_vertex[m])
                {
                    max_vertex[m]=(*Mother_geom)[i][k][m];
                }
            }
        }
    }
    Boite = Bbox (min_vertex, max_vertex);
    if ((I_fin-I_debut)<=MAX_TRIANGLE_N)
    {
        f_g =0;
        f_d = 0;
        
    }
    else// on doit encore séparer l'abre BVH  
    {
        int split_ind = I_debut ;
        Vector3D diag = max_vertex - min_vertex;
        Dim_of_split = diag.arg_max();
        //Value_of_split=(*Mother_geom)[split_ind].mean()[Dim_of_split]; 
        Value_of_split= (min_vertex + diag * (0.5))[Dim_of_split]; 
        
        for (int i =I_debut ; i<I_fin; i++)
        {
            if ((*Mother_geom)[i].mean()[Dim_of_split]<(Value_of_split))
            {
                if (split_ind==i)
                {
                    split_ind ++;
                }
                else 
                {
                    Triangle temp1 = (*Mother_geom)[i];
                    Triangle temp2 = (*Mother_geom)[split_ind];
                    (*Mother_geom)[i] = temp2;
                    (*Mother_geom)[split_ind] = temp1;
                    split_ind ++;
                }
                

            }
        }
    
        if (split_ind == I_debut || split_ind==I_fin)
        {
            f_g = 0;
            f_d = 0;
        }
        else 
        {
             f_g = new BVH(I_debut, split_ind,Mother_geom);
            f_d = new BVH(split_ind, I_fin,Mother_geom);
        }

    }
}




bool BVH::intersect(Ray ray,double& t, Object* &object )
{
    if (f_g!=0)// pas une feuille 
    {
        double t_g = 0;
        double t_d = 0;
        bool boolg = f_g->boite().intersect(ray, t_g, object  );
        bool boold = f_d->boite().intersect(ray, t_d, object  );
        bool b1=false, b2 =false;
        if (boolg&&boold)
        {
            if (t_g<t_d)
            {
                if (t>t_g)
                {
                    b1=f_g->intersect(ray,t,object);
                }
                if (t>t_d)
                {
                    b2 = f_d->intersect(ray,t,object);
                }
                
                //f_d->intersect(ray,t,object);
            }
            else 
            {
                if (t>t_d)
                {
                    b2 = f_d->intersect(ray,t,object);
                }
                if (t>t_g)
                {
                    b1 = f_g->intersect(ray,t,object);
                }
            }

        }
        else if (boolg)
            {
                if (t>t_g)
                {
                    b1 = f_g->intersect(ray,t,object);
                }	
            }
        else if (boold)
            {
                if (t>t_d)
                {
                    b2 = f_d->intersect(ray,t,object);
                }
            }
       return(b1 || b2);
    }
    else // on est une feuille 
    {
        double t_temp= DBL_MAX;
        int arg_t_min= -1;
        Object* obj=0;
        //if (!(Boite.intersect(ray,t_temp,obj)))
        //{
        //    return(false);
        //}
        t_temp= DBL_MAX;
        bool inter= false;
        for (int i =I_debut;i<I_fin;i++)
        {
            if ((*Mother_geom)[i].intersect(ray,t_temp,obj))
            {
                if (t_temp<t)
                {
                    inter=true;
                    t=t_temp;
                    object= &((*Mother_geom)[i]);
                }	
            }	
        }
        return(inter );
    }
}



