// Raytracer.cpp : Defines the entry point for the console application.
#define _CRT_SECURE_NO_WARNINGS // for Visual Studio 2017 (maybe 2015 as well)

#include <iostream>
#include <vector>
#include "geom_shape.h"
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <math.h> // sqrt
#include <stdio.h>
#include <math.h>
#include <random>
#include "geometry.h"


void save_image(const char* filename, const unsigned char* tableau, int w, int h) { // (0,0) is top-left corner

	FILE *f;

	int filesize = 54 + 3 * w*h;  

	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };
	
	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);
	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	f = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	unsigned char *row = new unsigned char[w * 3];
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++) 
		{
			row[j * 3] = tableau[(w*(h - i - 1) * 3) + j * 3+2];
			row[j * 3+1] = tableau[(w*(h - i - 1) * 3) + j * 3+1];
			row[j * 3+2] = tableau[(w*(h - i - 1) * 3) + j * 3];
		}
		fwrite(row, 3, w, f);
		fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
	}
	fclose(f);
	delete[] row;
}
//std::default_random_engine engine;
//std::uniform_real_distribution <double> u(0,1);


int main()
{
	int max_rbd_cnt=10;
	int N=75;

	Scene::MAX_RBDS_CNT=max_rbd_cnt;
	srand (time(NULL));
	int W = 1024;
	int H = 1024;
	double I= 5000000;
	double PI = 3.14159265;
	double R= 5.;
	double emissivity = I / (PI*R*R);
	std::vector<unsigned char> img(W*H * 3, 0);
	                   
	double alpha = 30.*2.*PI/180.;// =fov

	double tan_alpha_on_2 = tan (alpha/2);
	Vector3D center_lum (15,20,30);
	Vector3D center_lum2 (20,10,8);
	
	//img[(10 * W + 50)*3] = 255;  // pixel at (x, y) = (50, 10) is red
	Vector3D color1 (0.95,0.01,0.1);
	Vector3D color2 (0.05,0.95,0.1);
	Vector3D color3 (0.05,0.1,0.97);
	Vector3D color4 (0.6,1.0,0.2);
	Vector3D color5 (0.95,0.15,0.97);
	Vector3D color6 (0.1,0.95,0.9);
	Vector3D color_mir(0.75,0.75,0.75);

	Vector3D emissivity1(0.,0.,0.);
	Vector3D emissivity2(1.,1.,1.);
	Vector3D emissivity3(1.,0.2,0.);
	emissivity2 = emissivity2*emissivity;
	emissivity3 = emissivity3*(1.5*emissivity);

	//Vector3D color_mir(0.75,0.75,0.75);
	Sphere3D lum(center_lum,R,color1,emissivity2);
	Sphere3D lum2(center_lum2,R/2,color1,emissivity3);
	
	Sphere3D sphere1 (Vector3D(),10.);

	Sphere3D sphere2(Vector3D(1000,0,0),940.,color1,emissivity1);
	Sphere3D sphere3 (Vector3D(0,1000,0),940.,color2,emissivity1);
	Sphere3D sphere4 (Vector3D(-1000,0,0),940.,color3,emissivity1);
	Sphere3D sphere5 (Vector3D(0,-1000,0),940.,color4,emissivity1);
	Sphere3D sphere6 (Vector3D(0,0,1000),940.,color5,emissivity1);
	Sphere3D sphere7 (Vector3D(0,0,-1000),900.,color6,emissivity1);


	Sphere3D sphere8 (Vector3D(-25.,25.,-5.),10.,color_mir,emissivity1,true,false);
	Sphere3D sphere9 (Vector3D(25.,25.,-5.),5.,color1,emissivity1,false,true);
	Scene scene1;

	Triangle triangle1(Vector3D(0,34.6,-15),Vector3D(-20,0,-15),Vector3D(20.,0,-15));

	//Geometry geom1 ("cube.obj",5.,Vector3D(-30.,-30.,-30));

	Geometry* geom_pt2 = new Geometry ("girl.obj",38.,Vector3D(-0.,-50.,-0));
	scene1.add_sphere(lum);
	//scene1.add_sphere(lum2);
	//scene1.add_sphere(sphere1);
	scene1.add_sphere(sphere2);
	scene1.add_sphere(sphere3);
	scene1.add_sphere(sphere4);
	scene1.add_sphere(sphere5);
	scene1.add_sphere(sphere6);
	scene1.add_sphere(sphere7);
	scene1.add_sphere(sphere8);
	scene1.add_sphere(sphere9);
	//scene1.add_triangle(triangle1);
	//scene1.add_geom(&geom1);

	geom_pt2->add_texture("visage.bmp");
	geom_pt2->add_texture("cheveux.bmp");
	geom_pt2->add_texture("corps.bmp");
	geom_pt2->add_texture("pantalon.bmp");
	geom_pt2->add_texture("accessoires.bmp");
	geom_pt2->add_texture("mains.bmp");
	geom_pt2->build_obj();

	scene1.add_geom(geom_pt2);
	std::default_random_engine engine;

	std::uniform_real_distribution <double> u(0,1);
	engine.seed(time(NULL));
#pragma omp parallel for schedule(static,1)
	for (int i=0; i<H;i++)
	{
		std::cout<<i <<std::endl;
		for (int j=0;j<W;j++)
		{
			Vector3D computed_colors=Vector3D(0.,0.,0.) ;


			for (int k=0;k<N; k++)
			{
				
				double r1= u(engine);// antialiasing 
				double r2= u(engine);
				double r3= u(engine);
				double r4 = u(engine);
				double d = 55.; // distance à laquelle les rayons se croisents
				// pour la profondeur de champ
				double sigma = 0.25;
				double sigma2= 0.8;
				double mu=0.;
				 double xi= sqrt(-2*log(r1))*cos(2*PI*r2)*sigma+mu;
        		double yi= sqrt(-2*log(r1))*sin(2*PI*r2)*sigma+mu;
				
				
				double xi1= sqrt(-2*log(r3))*cos(2*PI*r4)*sigma2+mu;
        		double yi1= sqrt(-2*log(r3))*sin(2*PI*r4)*sigma2+mu;
				
				Vector3D u (j-W*0.5+xi+0.5, 0.5+i-H*0.5+yi, -W/(2*(tan_alpha_on_2)));
				u.normalize();
				
				u=u*d+Vector3D (-xi1,-yi1,0);// 1/ sqrt(3) pour tan(alpha/2)
				u.normalize();
				Vector3D c (xi1,yi1,55);// rempalcer par atan ????
				Ray rayon (c,u);
				rayon.normalize();
				Vector3D colors=Vector3D(0.,0.,0.) ;
				colors = scene1.pix_color(rayon,max_rbd_cnt,false);
				computed_colors = computed_colors + colors*(1/static_cast<double>(N));
			}
			double consta= 1./2.2;
			Vector3D computed_colors2 = computed_colors;
			double a = (computed_colors2.max());
			if (a>772.41)
			{
				computed_colors2=computed_colors2*(772.41/a);
			}
			img[((H-i-1) * W + j)*3] =std::min(255., floor(pow(255.*computed_colors2[0],consta)));
			img[((H-i-1) * W + j)*3+1] =std::min(255., floor(pow(255.*computed_colors2[1],consta)));
			img[((H-i-1) * W + j)*3+2] = std::min(255., floor(pow(255.*computed_colors2[2],consta)));
		}
		std::cout<<double(Triangle::COUNT_OF_INTERSECT)/(i+1)/W / (N)/(max_rbd_cnt+1)<<std::endl;
	}
	std::cout<<double(Triangle::COUNT_OF_INTERSECT)/H/W / N/(max_rbd_cnt+1)<<std::endl;
	save_image("test.bmp", &img[0], W, H);

    return 0;
}