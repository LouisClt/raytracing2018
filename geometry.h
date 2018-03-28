#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER

#include <vector>
#include <string>
#include "geom_shape.h"
#include <stdio.h>
#include <cstring>
#include <map>



class Bbox
{
	private:
		Vector3D Min_vertex, Max_vertex;
	public:
		Bbox(){};

		Bbox(Vector3D min_vertex, Vector3D max_vertex): Min_vertex(min_vertex),Max_vertex(max_vertex)
		{
		};
		bool intersect(Ray ray, double& t, Object* object=0)
		{
			Vector3D t_min, t_max;
			t_min[0] = ((Min_vertex-ray.center())*Vector3D(1,0,0))/(ray.u()
			*Vector3D(1,0,0));
			t_min[1] = ((Min_vertex-ray.center())*Vector3D(0,1,0))/(ray.u()
			*Vector3D(0,1,0));
			t_min[2] = ((Min_vertex-ray.center())*Vector3D(0,0,1))/(ray.u()
			*Vector3D(0,0,1));

			t_max[0] = ((Max_vertex-ray.center())*Vector3D(1,0,0))/(ray.u()
			*Vector3D(1,0,0));
			t_max[1] = ((Max_vertex-ray.center())*Vector3D(0,1,0))/(ray.u()
			*Vector3D(0,1,0));
			t_max[2] = ((Max_vertex-ray.center())*Vector3D(0,0,1))/(ray.u()
			*Vector3D(0,0,1));

			for (int i  =0; i <3; i++)
			{
				if (t_max[i]<t_min[i])
				{
					double temp= t_max[i];
					t_max[i]= t_min[i];
					t_min[i]= temp;
				}
			}
			double t_max_dbl=DBL_MAX;
			double t_min_dbl = 0;
			for (int i =0;i<3;i++)
			{
				if (t_max[i]<t_max_dbl)
				{
					t_max_dbl = t_max[i];
				}
				if (t_min[i]>t_min_dbl)
				{
					t_min_dbl = t_min[i];
				}
			}
			t = t_min_dbl;
			return(t_max_dbl >= t_min_dbl);

		};


};


class BVH
{
	private:
		Geometry* Mother_geom;
		BVH* f_g;
		BVH* f_d;
		Bbox Boite;
		double Value_of_split;
		int Dim_of_split;
		int I_debut, I_fin;

	public: 
		int MAX_TRIANGLE_N=20;
		BVH();
		Bbox boite()
		{
			return(Boite);
		};
		BVH(int i_debut, int i_fin , Geometry* mother_g);

		bool intersect(Ray ray,double& t, Object* &object );
	

};


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) {
	};
	int vtxi, vtxj, vtxk;
	int uvi, uvj, uvk;
	int ni, nj, nk;
	int faceGroup;
};

class Geometry: public Object {
public:
	Geometry() {};
	Geometry(const char* obj, double scaling, const Vector3D& offset) {
		readOBJ(obj);

		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * scaling + offset;
		}
		
	}
	void build_obj()
	{
		
		Vector3D max_vertex (0,0,0);
		Vector3D min_vertex (DBL_MAX,DBL_MAX,DBL_MAX);
		for (int i=0; i<indices.size();i++)
		{
			
			Tri_vect.push_back(
				Triangle(
					vertices[indices[i].vtxi],
					vertices[indices[i].vtxj],
					vertices[indices[i].vtxk]
					,normals[indices[i].ni],
					normals[indices[i].nj],
					normals[indices[i].nk], 
					uvs[indices[i].uvi],
					uvs[indices[i].uvj],
					uvs[indices[i].uvk],
					this,
					indices[i].faceGroup,
					w[indices[i].faceGroup],
					h[indices[i].faceGroup]
					));
			for (int k=0; k <3; k++)
			{
				for (int m = 0; m<3;m++)
				{
					if (Tri_vect[i][k][m]<min_vertex[m])
					{
						min_vertex[m]=Tri_vect[i][k][m];
					}
					if (Tri_vect[i][k][m]>max_vertex[m])
					{
						max_vertex[m]=Tri_vect[i][k][m];
					}
				}
			}
		}
		boite_tot = Bbox(min_vertex, max_vertex);
		bvh = new BVH(0,Tri_vect.size(),this);
	}
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");

		std::map<std::string, int> groupNames;
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				if (groupNames.find(std::string(grp)) != groupNames.end()) {
					curGroup = groupNames[std::string(grp)];
				}
				else {
					curGroup = groupNames.size();
					groupNames[std::string(grp)] = curGroup;
				}
			}
			if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
				sscanf(line, "mtllib %[^\n]\n", matfile);
			}
			if (line[0] == 'v' && line[1] == ' ') {
				Vector3D vec;
				Vector3D col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[2], &vec[1], &col[0], &col[1], &col[2]) == 6) {
					vertices.push_back(vec);
					vertexcolors.push_back(col);
				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);  // helmet
																				 //vec[2] = -vec[2]; //car2
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector3D vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); //girl
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector3D vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;

				char* consumedline = line + 1;
				int offset;
				t.faceGroup = curGroup;
				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;

					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
						
					}
				}


				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.faceGroup = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}


		}
		fclose(f);
	}

	bool intersect( Ray ray, double& t, Object* &object)
	{
		t = DBL_MAX;
		int arg_t_min= -1;
		double t_temp = DBL_MAX;
		Object* obj=0;
		if (!(boite_tot.intersect(ray,t_temp,obj)))
		{
			return(false);
		}
		bvh->intersect(ray,t,object);
		/*
		t_temp= DBL_MAX;
        bool inter= false;
        for (int i =0;i<Tri_vect.size();i++)
        {
            if ((Tri_vect)[i].intersect(ray,t_temp,obj))
            {
                if (t_temp<t)
                {
                    inter=true;
                    t=t_temp;
                    object= &((Tri_vect)[i]);
                }	
            }	
        }
        return(inter );//*/
		return(t!=DBL_MAX);
	};

	 Vector3D pix_color (Ray ray, Vector3D P,Scene* scene, int rebonds_count)
	{
	  return(Vector3D());
	};

	void add_texture(const char* filename) {

		textures.resize(textures.size() + 1);
		w.resize(w.size() + 1);
		h.resize(h.size() + 1);

		FILE* f;
		f = fopen(filename, "rb");
		unsigned char info[54];
		fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

		w[w.size() - 1] = *(int*)&info[18]; // extract image height and width from header
		h[h.size() - 1] = *(int*)&info[22];

		int size = 3 * w[w.size() - 1] * h[h.size() - 1];
		textures[textures.size() - 1].resize(size); // allocate 3 bytes per pixel
		fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
		fclose(f);

		for (int i = 0; i < size; i += 3) {
			std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
		}
	}


	Vector3D n(Vector3D X)
	{
		return(Vector3D());
	};
	
	Triangle& operator[](unsigned n)
	{
		return(Tri_vect[n]);
	};
	std::vector<TriangleIndices> indices;
	std::vector<Vector3D> vertices;
	std::vector<Vector3D> normals;
	std::vector<Vector3D> uvs; // Vector en 3D mais on n'utilise que 2 composantes
	std::vector<Vector3D> vertexcolors;
	std::vector <Triangle> Tri_vect;
	BVH* bvh;
	Bbox boite_tot;

	std::vector<std::vector<unsigned char> > textures;
	std::vector<int> w, h;

};





#endif