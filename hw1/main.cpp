#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "math.h"
#include <unistd.h>
#include <fstream>

typedef unsigned char RGB[3];



struct mainTriangle
{
	int material_id;

	parser::Vec3f vertex0;
	parser::Vec3f vertex1;
	parser::Vec3f vertex2;

	parser::Vec3f edge1;
	parser::Vec3f edge2;

	parser::Vec3f normalVector;
};

struct Ray
{
	parser::Vec3f org;
	parser::Vec3f d;
};

struct Color
{
	int r;
	int g;
	int b;
};

double dotProduct(parser::Vec3f v1 , parser::Vec3f v2);
parser::Vec3f doRayTracing(Ray ray, int max_recursion_depth);

parser::Scene scene;
std::vector<mainTriangle> globalTriangles;
std::vector<parser::Sphere> globalSpheres;

// VECTOR ISLEMLERI
	parser::Vec3f addTwoVector(parser::Vec3f v1, parser::Vec3f v2)
	{
		parser:: Vec3f add;
		add.x = v1.x + v2.x;
		add.y = v1.y + v2.y;
		add.z = v1.z + v2.z;
		return add;
	}


	parser::Vec3f vectorDiff(parser::Vec3f v1 , parser::Vec3f v2)
	{
		parser:: Vec3f diff;
		diff.x = v1.x - v2.x;
		diff.y = v1.y - v2.y;
		diff.z = v1.z - v2.z;
		return diff;

	}

	parser::Vec3f scalarVectorMult(float a, parser::Vec3f v1)
	{
		parser::Vec3f v2;
		v2.x = v1.x * a;
		v2.y = v1.y * a;
		v2.z = v1.z * a;

		return v2;
	}

	double dotProduct(parser::Vec3f v1 , parser::Vec3f v2){

		return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
	}


	parser::Vec3f crossProduct(parser::Vec3f v1, parser::Vec3f v2)
	{
		parser::Vec3f cross_P;


		cross_P.x = v1.y * v2.z - v1.z * v2.y; 
	    cross_P.y = -(v1.x * v2.z - v1.z * v2.x); 
	    cross_P.z = v1.x * v2.y - v1.y * v2.x; 
	    return cross_P;
	}
	parser::Vec3f normalize(parser::Vec3f v1)
	{
		// normalize
		double length = sqrt(dotProduct(v1,v1)); 

		parser::Vec3f v2;

		v2.x = v1.x/length;
		v2.y = v1.y/length;
		v2.z = v1.z/length;
		return v2;

	}

Ray findRayInThePixel(parser::Camera camera, parser::Vec3f pixel)
{
	Ray ray;
	ray.org.x = camera.position.x;
	ray.org.y = camera.position.y;
	ray.org.z = camera.position.z;

	// find d

	ray.d.x = pixel.x - ray.org.x;
	ray.d.y = pixel.y - ray.org.y;
	ray.d.z = pixel.z - ray.org.z;

	// normalize
	ray.d = normalize(ray.d);

	return ray;
}

bool raySphereIntersection(
	Ray ray, 
	parser::Sphere sphere, 
	double *t
	)
{
	float EPSILON = scene.shadow_ray_epsilon;

	parser::Vec3f L;
	L.x = (scene.vertex_data[sphere.center_vertex_id - 1].x) - ray.org.x;
	L.y = (scene.vertex_data[sphere.center_vertex_id - 1].y) - ray.org.y;
	L.z = (scene.vertex_data[sphere.center_vertex_id - 1].z) - ray.org.z;


	double tca = dotProduct(L,ray.d);
	
	if(tca < 0)
	//if(tca < EPSILON)
		return false;

	else{
		
		double s2 = (dotProduct(L,L))-(tca*tca);
		double s = sqrt(s2);
		
		if(s>sphere.radius)
			return false;

		else
		{
			double thc = sqrt((sphere.radius*sphere.radius) - s2);
			double t0 = tca - thc;
			*t = t0;
			return true;
		}
	}
}

bool rayTriangleIntersection(
	Ray ray, 
	mainTriangle triangle, 
	double *t
	)
{

	float EPSILON = 0.0000001;
	//float EPSILON = scene.shadow_ray_epsilon * (1.0/10.0);

    parser::Vec3f h, s, q;
    float a,f,u,v;

    //std::cout << "ray.d " << std::endl;

    //std::cout << ray.d.x << " " << ray.d.y << " " << ray.d.z << std::endl;
    h = crossProduct(ray.d, triangle.edge2);
    a = dotProduct(triangle.edge1, h);

    if (a > -EPSILON && a < EPSILON)
        return false;

    f = 1.0/a;

    s = vectorDiff(ray.org, triangle.vertex0);
    u = f * dotProduct(s, h);
    if (u < 0.0 || u > 1.0)
        return false;

    q = crossProduct(s, triangle.edge1);
    v = f * dotProduct(ray.d, q);

    if (v < 0.0 || u + v > 1.0)
        return false;
    
    float t1 = f * dotProduct(triangle.edge2, q);
    
   	//std::cout << "T 1: " << t1 << std::endl; 
    if (t1 > EPSILON) // ray intersection
    {
    	*t = t1;
    	//std::cout << "EPSILON " << std::endl;
        return true;
    }
    else
    	return false;
}	

// finds ambient light
parser::Vec3f findL_a(int hitMaterial_id)
{	
	parser::Vec3f k_a = scene.materials[hitMaterial_id-1].ambient;
	//find I_a;
	//I_a  = ambient light
	parser::Vec3f I_a = scene.ambient_light;
	
	//find L_a
	//L_a = (k_a)(I_a)
	parser::Vec3f L_a;
	L_a.x = k_a.x * I_a.x;
	L_a.y = k_a.y * I_a.y;
	L_a.z = k_a.z * I_a.z;

	return L_a;
}

// finds E(intensity/rsquare)
parser::Vec3f findE(parser::Vec3f wi, parser::Vec3f I)
{
	float rSquare= dotProduct(wi,wi);

	parser::Vec3f E;
	E.x = I.x/rSquare;
	E.y = I.y/rSquare;
	E.z = I.z/rSquare;

	return E;
}

parser::Vec3f findL_d(
	parser::Vec3f wi, 
	parser::Vec3f normalVector, 
	parser::Vec3f E,
	int hitMaterial_id
	)
{
	// find cos_teta
	//cos_teta = max{0,dot(wi,n)}
	parser::Vec3f L_d;
	float cos_teta = dotProduct(wi,normalVector);

	if(cos_teta < 0.0)
	{
		L_d.x = 0;
		L_d.y = 0;
		L_d.z = 0;
	}
	else
	{
		//find L_d
		// L_d = (k_d)(cos_teta)(I/rSquare)
		parser::Vec3f k_d = scene.materials[hitMaterial_id-1].diffuse;
		L_d.x = k_d.x * E.x;
		L_d.y = k_d.y * E.y;
		L_d.z = k_d.z * E.z;

		L_d = scalarVectorMult(cos_teta,L_d);
	}

	return L_d;
}

parser::Vec3f findL_s(
	Ray ray,
	parser::Vec3f intersectionPoint,
	parser::Vec3f wi,
	parser::Vec3f normalVector,
	parser::Vec3f E,
	int hitMaterial_id
	)
{	
	//find wo
	//wo = intersectpoint - camerapos;
	//camerapo = ray.org
	parser::Vec3f wo= normalize(vectorDiff(ray.org,intersectionPoint));
	

	//find h
	//h = (wi+wo) / ||wi+wo||
	parser::Vec3f h = normalize(addTwoVector(wi,wo));

	//find cos_a
	//cos_a = max{0,dot(n,h)}
	float cos_a = dotProduct(normalVector,h);
	parser::Vec3f L_s;

	if(cos_a<0.0)
	{
		L_s.x = 0;
		L_s.y = 0;
		L_s.z = 0;
	}
	else
	{
		parser::Vec3f k_s = scene.materials[hitMaterial_id-1].specular;

		L_s.x = k_s.x * E.x;
		L_s.y = k_s.y * E.y;
		L_s.z = k_s.z * E.z;
		//find ph_Exp
		//ph_Exp = phong exponent

		float ph_Exp = scene.materials[hitMaterial_id-1].phong_exponent;
		while(ph_Exp>0)
		{
			L_s = scalarVectorMult(cos_a,L_s);
			ph_Exp--;
		}
	}
	return L_s;
}

parser::Vec3f doClamping(parser::Vec3f L)
{

	if(L.x > 255)
		L.x = 255;
	if(L.x < 0)
		L.x = 0;
	if(L.y > 255)
		L.y = 255;
	if(L.y < 0)
		L.y = 0;
	if(L.z > 255)
		L.z = 255;
	if(L.z < 0)
		L.z = 0;

	return L;
}

parser::Vec3f findL_m(
	Ray ray, 
	parser::Vec3f intersectionPoint, 
	parser::Vec3f normalVector,
	int hitMaterial_id,
	int max_recursion_depth
	)
{
	parser::Vec3f L_m;
	L_m.x = 0;
	L_m.y = 0;
	L_m.z = 0;

	parser::Vec3f k_m = scene.materials[hitMaterial_id-1].mirror;

	if(k_m.x == 0 && k_m.y == 0 && k_m.z == 0)
	{
		return L_m;
	}

	parser::Vec3f wo= normalize(vectorDiff(ray.org,intersectionPoint));
	
	float n_wo = dotProduct(normalVector, wo);
	n_wo = 2 * n_wo;

	parser::Vec3f wr = addTwoVector(scalarVectorMult(-1, wo), scalarVectorMult(n_wo, normalVector));

	wr = normalize(wr);
	
	Ray reflectionRay;
	/*
	reflectionRay.org.x = intersectionPoint.x;
	reflectionRay.org.y = intersectionPoint.y;
	reflectionRay.org.z = intersectionPoint.z;
	*/
	reflectionRay.org.x = intersectionPoint.x + scene.shadow_ray_epsilon;
	reflectionRay.org.y = intersectionPoint.y + scene.shadow_ray_epsilon;
	reflectionRay.org.z = intersectionPoint.z + scene.shadow_ray_epsilon;
	
	reflectionRay.d = wr;
	/*
	reflectionRay.d.x = wr.x + scene.shadow_ray_epsilon;
	reflectionRay.d.y = wr.y + scene.shadow_ray_epsilon;
	reflectionRay.d.z = wr.z + scene.shadow_ray_epsilon;

	reflectionRay.d = normalize(reflectionRay.d);
	*/

	L_m = doRayTracing(reflectionRay, max_recursion_depth);
	L_m.x = k_m.x * L_m.x;
	L_m.y = k_m.y * L_m.y;
	L_m.z = k_m.z * L_m.z;

	//L_m = addTwoVector(araVector, L_m);
	return L_m;
}

bool isWiIntersect(
	parser::Vec3f wi, 
	parser::Vec3f intersectionPoint,
	parser::Vec3f lightPosition
	)
{
	const float EPSILON = scene.shadow_ray_epsilon;

	parser::Vec3f diff = vectorDiff(lightPosition, intersectionPoint);

	float tToLight = sqrt(dotProduct(diff, diff));

	parser::Vec3f wiEpsilon = scalarVectorMult(EPSILON, wi);
	
	Ray shadowRay;
	shadowRay.org = addTwoVector(intersectionPoint, wiEpsilon);
	shadowRay.d = wi;

	double t0 = 0.0f;

	int sphereSize = globalSpheres.size();
	for (int i = 0; i < sphereSize; i++)
	{
		bool hit = raySphereIntersection(shadowRay, globalSpheres[i], &t0);

		if (hit /*&& t0 > EPSILON*/ && t0 < tToLight)
		{
			return true;
		}
	}
	int triangleSize = globalTriangles.size();
	for(int k = 0; k < triangleSize; k++){

		bool hit = rayTriangleIntersection(shadowRay, globalTriangles[k], &t0);

		if (hit /*&& t0 > EPSILON*/ && t0 < tToLight)
		{
			return true;
		}
	}

	return false;
}

parser::Vec3f doRayTracing(Ray ray, int max_recursion_depth)
{
	if(max_recursion_depth == 0)
	{
		parser::Vec3f L;
		L.x = 0;
		L.y = 0;
		L.z = 0;
		return L;
	}
	max_recursion_depth--;

	double minT = 2000000.0;
	int sphereHit = -1;
	int trianglehit = -1;
	double t0 = 0.0f;

	int hitMaterial_id = 0;

	parser::Vec3f normalVector;
	parser::Vec3f intersectionPoint;

	int sphereSize = globalSpheres.size();
	for(int k = 0; k < sphereSize; k++)
	{
		bool hit = raySphereIntersection(ray, globalSpheres[k], &t0);

		if (hit && t0< minT)
		{
			minT = t0;
			sphereHit = k;

			hitMaterial_id = globalSpheres[k].material_id;

			parser::Vec3f center = scene.vertex_data[globalSpheres[k].center_vertex_id - 1];
			
			intersectionPoint = scalarVectorMult(t0, ray.d);
			intersectionPoint = addTwoVector(intersectionPoint, ray.org);

			normalVector = vectorDiff(intersectionPoint, center);
			normalVector = normalize(normalVector);
		}
	}
	int triangleSize = globalTriangles.size();
	for(int k = 0; k < triangleSize; k++){

		bool hit = rayTriangleIntersection(ray, globalTriangles[k], &t0);

		if (hit && t0<minT)
		{
			//std::cout << "hit" << std::endl;
			minT = t0;
			trianglehit = k;
			sphereHit = -1;
			// which materaial id hits
			hitMaterial_id = globalTriangles[k].material_id;

			normalVector = globalTriangles[k].normalVector;

			intersectionPoint = scalarVectorMult(t0, ray.d);
			intersectionPoint = addTwoVector(intersectionPoint, ray.org);
		}
	}

	parser::Vec3f L;

	if(sphereHit != -1 || trianglehit != -1)
	{
		
		parser::Vec3f L_a = findL_a(hitMaterial_id);

		L = L_a;

		for(int k = 0; k < scene.point_lights.size(); k++)
		{
			parser::Vec3f wi = vectorDiff(scene.point_lights[k].position,intersectionPoint);
			
			parser::Vec3f I = scene.point_lights[k].intensity;

			parser::Vec3f E = findE(wi, I);

			wi = normalize(wi);
			
			// point light I
			if(isWiIntersect(wi, intersectionPoint, scene.point_lights[k].position))
				continue;
			else
			{
				parser::Vec3f L_d = findL_d(wi, normalVector, E, hitMaterial_id);

				parser::Vec3f L_s = findL_s(ray, intersectionPoint, wi, normalVector, E, hitMaterial_id);

				
				L = addTwoVector(L, addTwoVector(L_s,L_d));
				

				//L = addTwoVector(L, addTwoVector(L_s, L_d));
				/*
				parser::Vec3f L_m;

				parser::Vec3f k_m = scene.materials[hitMaterial_id-1].mirror;
				if(k_m.x == 0 && k_m.y == 0 && k_m.z == 0)
				{
					L_m.x = 0;
					L_m.y = 0;
					L_m.z = 0;
				}
				else
				{
					parser::Vec3f wo= normalize(vectorDiff(ray.org,intersectionPoint));
					
					float n_wo = dotProduct(normalVector, wo);
					n_wo = 2 * n_wo;

					parser::Vec3f wr = addTwoVector(scalarVectorMult(-1, wo), scalarVectorMult(n_wo, normalVector));

					wr = normalize(wr);

					Ray reflectionRay;
					reflectionRay.org = intersectionPoint;
					reflectionRay.d = wr;				

					L_m = doRayTracing(reflectionRay, max_recursion_depth);

					L_m.x = k_m.x * L_m.x;
					L_m.y = k_m.y * L_m.y;
					L_m.z = k_m.z * L_m.z;
				}
				
				L = addTwoVector(
					L_m,
					addTwoVector(L, addTwoVector(L_s,L_d))
					);
				*/
			}
		}

		L = addTwoVector(
		findL_m(ray, intersectionPoint, normalVector, hitMaterial_id, max_recursion_depth),
		L);

	} 
	else
	{
		L.x = scene.background_color.x;
		L.y = scene.background_color.y;
		L.z = scene.background_color.z;
	}

	
	return L;
}

void renderForEachCamera(parser::Camera camera)
{
	std::vector<parser::Vec3f> pixels_pos;

	Color** pixelColors = new Color*[camera.image_height];
	for (int i = 0; i < camera.image_height; i++)
	{
		pixelColors[i] = new Color[camera.image_width];
	}

	/*
	// first, find center
	parser::Vec3f center = addTwoVector(camera.position, scalarVectorMult(camera.near_distance,camera.gaze));

	// second, find q
	parser::Vec3f v = normalize(camera.up);
	// calculate u
	parser::Vec3f minusGaze = normalize(scalarVectorMult(-1, camera.gaze));
	parser::Vec3f u = crossProduct(v, minusGaze);
	parser::Vec3f u_left = scalarVectorMult(camera.near_plane.x, u);
	parser::Vec3f top_v = scalarVectorMult(camera.near_plane.w, v);
	parser::Vec3f q;
	q = addTwoVector(center, u_left);
	q = addTwoVector(q, top_v);
	*/
	parser::Vec3f center = addTwoVector(camera.position, scalarVectorMult(camera.near_distance,normalize(camera.gaze)));

	// v is up 
	parser::Vec3f v = camera.up;	

	// w is minus gaze
	parser::Vec3f w = normalize(scalarVectorMult(-1, camera.gaze));

	parser::Vec3f u = normalize(crossProduct(v, w));

	v = crossProduct(w,u);
	

	parser::Vec3f u_left = scalarVectorMult(camera.near_plane.x, u);
	parser::Vec3f top_v = scalarVectorMult(camera.near_plane.w, v);
	parser::Vec3f q;
	q = addTwoVector(center, u_left);
	q = addTwoVector(q, top_v);
	

	for (int i = 0; i < camera.image_height; i++)
	{
		for (int j = 0; j < camera.image_width; j++)
		{

			float s_u = (camera.near_plane.y - camera.near_plane.x) * ((j + 0.5) / (float)camera.image_width);
			float s_v = (camera.near_plane.w - camera.near_plane.z) * ((i + 0.5) / (float)camera.image_height);

			parser::Vec3f su_u = scalarVectorMult(s_u, u);
			parser::Vec3f sv_v = scalarVectorMult(s_v, v);
			sv_v = scalarVectorMult(-1, sv_v);

			parser::Vec3f pixelCoord;
			pixelCoord = addTwoVector(su_u, sv_v);
			pixelCoord = addTwoVector(pixelCoord, q);

			Ray ray = findRayInThePixel(camera ,pixelCoord);

			parser::Vec3f currentPixelColor = doRayTracing(ray, scene.max_recursion_depth + 1);
			
			//std::cout << "i: " << i << " j: " << j << std::endl;
			currentPixelColor = doClamping(currentPixelColor);

			pixelColors[i][j].r = currentPixelColor.x;
			pixelColors[i][j].g = currentPixelColor.y;
			pixelColors[i][j].b = currentPixelColor.z;
		}
	}

	int l = 0;
	unsigned char* image = new unsigned char [camera.image_width * camera.image_height * 3];
	
	for (int i = 0; i < camera.image_height; i++)
	{
		for (int j = 0; j < camera.image_width; j++)
		{
			//std::cout<< "l: " << l << " i: " << i << " j: " << j << std::endl;
            image[l++] = pixelColors[i][j].r;
            image[l++] = pixelColors[i][j].g;
            image[l++] = pixelColors[i][j].b;
		}
	}

	write_ppm(camera.image_name.c_str(), image, camera.image_width, camera.image_height);

}

// store materials in global
	void parseTriangles()
	{
		int triangleSize = scene.triangles.size();
	    for (int i = 0; i < triangleSize; i++)
	    {
	    	mainTriangle newTriangle;

	    	newTriangle.material_id = scene.triangles[i].material_id;

	    	newTriangle.vertex0 = scene.vertex_data[scene.triangles[i].indices.v0_id -1];
	    	newTriangle.vertex1 = scene.vertex_data[scene.triangles[i].indices.v1_id -1];
	    	newTriangle.vertex2 = scene.vertex_data[scene.triangles[i].indices.v2_id -1];

		    newTriangle.edge1 = vectorDiff(newTriangle.vertex1, newTriangle.vertex0);
		    newTriangle.edge2 = vectorDiff(newTriangle.vertex2, newTriangle.vertex0);

			newTriangle.normalVector = crossProduct(newTriangle.edge1, newTriangle.edge2);
			newTriangle.normalVector = normalize(newTriangle.normalVector);

			globalTriangles.push_back(newTriangle);
	    }
	}

	void parseMeshes()
	{
		int meshesSize = scene.meshes.size();
	    for (int i = 0; i < meshesSize; i++)
	    {
	    	int facesSize = scene.meshes[i].faces.size();
	    	for(int l = 0; l < facesSize; l++)
	    	{
	    		mainTriangle newTriangle;

	    		newTriangle.material_id = scene.meshes[i].material_id;

	    		newTriangle.vertex0 = scene.vertex_data[scene.meshes[i].faces[l].v0_id -1];
	    		newTriangle.vertex1 = scene.vertex_data[scene.meshes[i].faces[l].v1_id -1];
	    		newTriangle.vertex2 = scene.vertex_data[scene.meshes[i].faces[l].v2_id -1];


	    		newTriangle.edge1 = vectorDiff(newTriangle.vertex1, newTriangle.vertex0);
			    newTriangle.edge2 = vectorDiff(newTriangle.vertex2, newTriangle.vertex0);

				newTriangle.normalVector = crossProduct(newTriangle.edge1, newTriangle.edge2);
				newTriangle.normalVector = normalize(newTriangle.normalVector);

				globalTriangles.push_back(newTriangle);
			}
	    }
	}

	void parseSpheres()
	{
		int sphereSize = scene.spheres.size();
		for (int i = 0; i < sphereSize; i++)
	    {
	    	parser::Sphere newSphere;
			newSphere.material_id = scene.spheres[i].material_id;
	    	newSphere.center_vertex_id = scene.spheres[i].center_vertex_id;
	    	newSphere.radius = scene.spheres[i].radius;
	    	globalSpheres.push_back(newSphere);
	    }
	}

int main(int argc, char* argv[])
{
    scene.loadFromXml(argv[1]);
    
    parseTriangles();
    parseMeshes();
    parseSpheres();

    for (int i = 0; i < scene.cameras.size(); ++i)
    {
    	renderForEachCamera(scene.cameras[i]);
    }

    return 0;
}
