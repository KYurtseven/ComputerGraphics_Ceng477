#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "hw2_types.h"
#include "hw2_math_ops.h"
#include "hw2_file_ops.h"
#include <iostream>
#include <vector>
#include <math.h>

///////////////// Global variables - constants - structures//////////////

    #define PI 3.14159265358979323846264338327950288

    typedef struct{

        int vertexIds[3];
        double vertex1[4];
        double vertex2[4];
        double vertex3[4];
        double center[4];
    
    } Mytriangle;


    Camera cameras[100];
    int numberOfCameras = 0;

    Model models[1000];
    int numberOfModels = 0;

    Color colors[100000];
    int numberOfColors = 0;

    Translation translations[1000];
    int numberOfTranslations = 0;

    Rotation rotations[1000];
    int numberOfRotations = 0;

    Scaling scalings[1000];
    int numberOfScalings = 0;

    Vec3 vertices[100000];
    int numberOfVertices = 0;

    Color backgroundColor;

    // backface culling setting, default disabled
    int backfaceCullingSetting = 0;

    Color **image;

///////////////////////////////////////////////////////////////////////


void multiply_3_4_MatrixWithVec4d(double r[3], double m[3][4], double v[4]) {
    int i, j;
    double total;
    for (i = 0; i < 3; i++) {
        total = 0;
        for (j = 0; j < 4; j++)
            total += m[i][j] * v[j];
        r[i] = total;
    }

}

Mytriangle findMytriangle(Triangle t){
    int index;

    Mytriangle ct;

    index = t.vertexIds[0];
    ct.vertexIds[0] = index;
    ct.vertex1[0] = vertices[index].x;
    ct.vertex1[1] = vertices[index].y;
    ct.vertex1[2] = vertices[index].z;
    ct.vertex1[3] = 1.0;

    index = t.vertexIds[1];
    ct.vertexIds[1] = index;
    ct.vertex2[0] = vertices[index].x;
    ct.vertex2[1] = vertices[index].y;
    ct.vertex2[2] = vertices[index].z;
    ct.vertex2[3] = 1.0;

    index = t.vertexIds[2];
    ct.vertexIds[2] = index;
    ct.vertex3[0] = vertices[index].x;
    ct.vertex3[1] = vertices[index].y;
    ct.vertex3[2] = vertices[index].z;
    ct.vertex3[3] = 1.0;


    ct.center[0] = (ct.vertex1[0] + ct.vertex2[0] + ct.vertex3[0])/3.0;
    ct.center[1] = (ct.vertex1[1] + ct.vertex2[1] + ct.vertex3[1])/3.0;
    ct.center[2] = (ct.vertex1[2] + ct.vertex2[2] + ct.vertex3[2])/3.0;
    ct.center[3] = 1.0;
    
    return ct;

}

void createMvp(Camera cam,double Mvp[3][4]){

    Mvp[0][0] = cam.sizeX / 2.0;
    Mvp[0][1] = 0.0;
    Mvp[0][2] = 0.0;
    Mvp[0][3] = (cam.sizeX - 1.0) / 2.0 ;
    
    Mvp[1][0] = 0.0;
    Mvp[1][1] = cam.sizeY / 2.0;
    Mvp[1][2] = 0.0;
    Mvp[1][3] = (cam.sizeY - 1.0 ) / 2.0 ;
    
    Mvp[2][0] = 0.0;
    Mvp[2][1] = 0.0;
    Mvp[2][2] = 1.0 / 2.0 ; 
    Mvp[2][3] = 1.0 / 2.0 ;

}

void createMper(Camera cam,double Mper[4][4]){

    Mper[0][0] = 2.0 * cam.n / (cam.r-cam.l);
    Mper[0][1] = 0.0;
    Mper[0][2] = (cam.r + cam.l) / (cam.r - cam.l);
    Mper[0][3] = 0.0;
    
    Mper[1][0] = 0.0;
    Mper[1][1] = 2.0 * cam.n / (cam.t-cam.b);
    Mper[1][2] = (cam.t + cam.b) / (cam.t - cam.b);
    Mper[1][3] = 0.0;
    
    Mper[2][0] = 0.0;
    Mper[2][1] = 0.0;
    Mper[2][2] = -1.0* ((cam.f + cam.n) / (cam.f - cam.n));
    Mper[2][3] = -2.0* cam.f * cam.n / (cam.f-cam.n);
    
    Mper[3][0] = 0.0;
    Mper[3][1] = 0.0;
    Mper[3][2] = -1.0;
    Mper[3][3] = 0.0;

}

void createMcam(Camera cam,double Mcam[4][4]){
    
    double r1 = (cam.u.x * cam.pos.x) + (cam.u.y * cam.pos.y) + (cam.u.z * cam.pos.z);
    double r2 = cam.v.x * cam.pos.x + cam.v.y * cam.pos.y + cam.v.z * cam.pos.z; 
    double r3 = cam.w.x * cam.pos.x + cam.w.y * cam.pos.y + cam.w.z * cam.pos.z;  

    Mcam[0][0] = cam.u.x ;
    Mcam[0][1] = cam.u.y ;
    Mcam[0][2] = cam.u.z ;
    Mcam[0][3] = -1.0 * r1;
    
    Mcam[1][0] = cam.v.x ;
    Mcam[1][1] = cam.v.y ;
    Mcam[1][2] = cam.v.z ;
    Mcam[1][3] = -r2;
    
    Mcam[2][0] = cam.w.x ;
    Mcam[2][1] = cam.w.y ;
    Mcam[2][2] = cam.w.z ;
    Mcam[2][3] = -r3;
    
    Mcam[3][0] = 0.0;
    Mcam[3][1] = 0.0;
    Mcam[3][2] = 0.0;
    Mcam[3][3] = 1.0;

}

void createScaleMatrix(Scaling scal, double S[4][4]){
    
    makeIdentityMatrix(S);
    S[0][0] = scal.sx;
    S[1][1] = scal.sy;
    S[2][2] = scal.sz;

}

void createRotationMatrix(Rotation rot, double R[4][4]){
    
    double angle_radian = (rot.angle * PI ) / 180.0 ;
    Vec3 u;

    u.x = rot.ux;
    u.y = rot.uy;
    u.z = rot.uz;
    
    u = normalizeVec3(u);

    double sin_a = sin(angle_radian);
    double cos_a = cos(angle_radian);

    R[0][0] = cos_a + u.x * u.x * (1.0 - cos_a); 
    R[0][1] = u.x * u.y * (1.0 - cos_a) - (u.z * sin_a);
    R[0][2] = u.x * u.z * (1.0 - cos_a) + (u.y * sin_a);
    R[0][3] = 0.0;
    
    R[1][0] = u.x * u.y * (1.0 - cos_a) + (u.z * sin_a); 
    R[1][1] = cos_a + u.y * u.y * (1.0 - cos_a);
    R[1][2] = u.y * u.z * (1.0 - cos_a) - (u.x * sin_a);
    R[1][3] = 0.0;

    R[2][0] = u.x * u.z * (1.0 - cos_a) - (u.y * sin_a);
    R[2][1] = u.y * u.z * (1.0 - cos_a) + (u.x * sin_a);
    R[2][2] = cos_a + u.z * u.z * (1.0 - cos_a);
    R[2][3] = 0.0;

    R[3][0] = 0.0;
    R[3][1] = 0.0;
    R[3][2] = 0.0;
    R[3][3] = 1.0;

}

void createTranslationMatrix(Translation trans, double T[4][4]){
    
    makeIdentityMatrix(T);

    T[0][3] = trans.tx;
    T[1][3] = trans.ty;
    T[2][3] = trans.tz;

}

void perspectiveDivide(double m[4]){
    
    double w = m[3];

    m[0] = m[0] / w;
    m[1] = m[1] / w;
    m[2] = m[2] / w;
    m[3] = m[3] / w;

}

void draw(double x, double y, Color c){

    image[(int)x][(int)y].r = (int) c.r;
    image[(int)x][(int)y].g = (int) c.g;
    image[(int)x][(int)y].b = (int) c.b; 

}

void draw1(double x0 , double x1 , double y0 , double y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(x1-x0);
    double colorConstantG = (color1.g-color0.g)/(x1-x0);
    double colorConstantB = (color1.b-color0.b)/(x1-x0);

    Color c = color0;

    double y = y0;
    double midPoint = (x0+1) * (y0 - y1) + (y0+0.5)*(x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = x0 ; i <= x1 ; i++){


        draw(i,y,c);

        c.r += colorConstantR;
        c.b += colorConstantB;
        c.g += colorConstantG;
        
        //below the line
        if(midPoint < 0.0){
            //go NE
            y++;
            midPoint = midPoint + (y0 - y1) + (x1 - x0);
        }

        //on the line or above the line
        else
            //go E
            midPoint = midPoint + (y0 - y1);
    }

}

void draw2(double x0 , double x1, double y0, double y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(y1-y0);
    double colorConstantG = (color1.g-color0.g)/(y1-y0);
    double colorConstantB = (color1.b-color0.b)/(y1-y0);

    Color c = color0;

    double x = x0;
    double midPoint = (x0 + 0.5)*(y0 - y1) + (y0 + 1) *(x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = y0 ; i <= y1 ; i++){

        draw(x,i,c);

        c.r += colorConstantR;
        c.b += colorConstantB;
        c.g += colorConstantG;
        
        if(midPoint > 0.0){
            //go NE
            x++;
            midPoint = midPoint + (y0 - y1) + (x1 - x0);
        }

        else
            //go N
            midPoint = midPoint + (x1 - x0);
    }

}

void draw3(double x0 , double x1 , double y0 , double y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(x1-x0);
    double colorConstantG = (color1.g-color0.g)/(x1-x0);
    double colorConstantB = (color1.b-color0.b)/(x1-x0);

    Color c = color0;

    double y = y0;
    double midPoint = (x0 + 1.0) * (y0 - y1) + (y0 - 0.5)* (x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = x0 ; i <= x1 ; i++){

        draw(i,y,c);

        c.r += colorConstantR;
        c.b += colorConstantB;
        c.g += colorConstantG;
        
        //above the line
        if(midPoint > 0.0){
            //go SE
            y--;
            midPoint = midPoint + (y0 - y1) - (x1 - x0);
        }

        else
        //below the line or on the line
            midPoint = midPoint + (y0 - y1);
    }

}

void draw4(double x0 , double x1 , double y0 , double y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(y1-y0);
    double colorConstantG = (color1.g-color0.g)/(y1-y0);
    double colorConstantB = (color1.b-color0.b)/(y1-y0);

    Color c = color0;

    double x = x0;
    double midPoint = (x0 + 0.5)*(y0 - y1) + (y0 - 1) *(x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = y0 ; i >= y1 ; i--){

        draw(x,i,c);

        c.r -= colorConstantR;
        c.b -= colorConstantB;
        c.g -= colorConstantG;
        
        //below the line
        if(midPoint < 0.0){
            //go NE
            x++;
            midPoint = midPoint + (y0 - y1) - (x1 - x0);
        }

        //on the line or above the line
        else
            //go E
            midPoint = midPoint - (x1 - x0);
    }

}

void draw5(double x0 , double x1 , double y0 , double y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(x1-x0);
    double colorConstantG = (color1.g-color0.g)/(x1-x0);
    double colorConstantB = (color1.b-color0.b)/(x1-x0);

    Color c = color0;
    
    for(int i = x0 ; i <=x1 ; i++){

        draw(i,y0,c);
        c.r += colorConstantR;
        c.b += colorConstantB;
        c.g += colorConstantG;
    }

}

void draw6(int x0 , int x1 , int y0 , int y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(x1-x0);
    double colorConstantG = (color1.g-color0.g)/(x1-x0);
    double colorConstantB = (color1.b-color0.b)/(x1-x0);

    Color c = color0;

    double y = y0;
    double midPoint = (x0 - 1.0)*(y0 - y1) + (y0 - 0.5) *(x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = x0 ; i >= x1 ; i--){

        draw(i,y,c);

        c.r -= colorConstantR;
        c.b -= colorConstantB;
        c.g -= colorConstantG;
        
        //below the line
        if(midPoint < 0.0){
            y--;
            midPoint = midPoint - (y0 - y1) - (x1 - x0);
        }

        //on the line or above the line
        else
            //go w
            midPoint = midPoint - (y0 - y1);
    }
}


void draw7(int x0 , int x1 , int y0 , int y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(y1-y0);
    double colorConstantG = (color1.g-color0.g)/(y1-y0);
    double colorConstantB = (color1.b-color0.b)/(y1-y0);

    Color c = color0;

    double x = x0;
    double midPoint = (x0 - 0.5)*(y0 - y1) + (y0 - 1.0) *(x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = y0 ; i >= y1 ; i--){

        draw(x,i,c);

        c.r -= colorConstantR;
        c.b -= colorConstantB;
        c.g -= colorConstantG;
        
        //below the line
        if(midPoint > 0.0){
            x--;
            midPoint = midPoint - (y0 - y1) - (x1 - x0);
        }

        //on the line or above the line
        else
            //go w
            midPoint = midPoint - (x1 - x0);
    }
}


void draw8(int x0 , int x1 , int y0 , int y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(x1-x0);
    double colorConstantG = (color1.g-color0.g)/(x1-x0);
    double colorConstantB = (color1.b-color0.b)/(x1-x0);

    Color c = color0;

    double y = y0;
    double midPoint = (x0 - 1.0)*(y0 - y1) + (y0 + 0.5) *(x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = x0 ; i >= x1 ; i--){

        draw(i,y,c);

        c.r -= colorConstantR;
        c.b -= colorConstantB;
        c.g -= colorConstantG;
        
        //below the line
        if(midPoint > 0.0){
            y++;
            midPoint = midPoint - (y0 - y1) + (x1 - x0);
        }

        //on the line or above the line
        else
            //go w
            midPoint = midPoint - (y0 - y1);
    }
}

void draw9(int x0 , int x1 , int y0 , int y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(y1-y0);
    double colorConstantG = (color1.g-color0.g)/(y1-y0);
    double colorConstantB = (color1.b-color0.b)/(y1-y0);

    Color c = color0;

    double x = x0;
    double midPoint = (x0 - 0.5)*(y0 - y1) + (y0 + 1) *(x1 - x0) + (x0*y1) - (x1*y0);

    for(int i = y0 ; i <= y1 ; i++){

        draw(x,i,c);

        c.r += colorConstantR;
        c.b += colorConstantB;
        c.g += colorConstantG;
        
        //below the line
        if(midPoint < 0.0){
            x--;
            midPoint = midPoint - (y0 - y1) + (x1 - x0);
        }

        //on the line or above the line
        else
            //go w
            midPoint = midPoint + (x1 - x0);
    }

}

void draw10(int x0 , int x1 , int y0 , int y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(x1-x0);
    double colorConstantG = (color1.g-color0.g)/(x1-x0);
    double colorConstantB = (color1.b-color0.b)/(x1-x0);

    Color c = color0;
    
    for(int i = x0 ; i >=x1 ; i--){

        draw(i,y0,c);
        c.r -= colorConstantR;
        c.b -= colorConstantB;
        c.g -= colorConstantG;
    }

}


void draw11(int x0 , int x1 , int y0 , int y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(y1-y0);
    double colorConstantG = (color1.g-color0.g)/(y1-y0);
    double colorConstantB = (color1.b-color0.b)/(y1-y0);

    Color c = color0;
    
    for(int i = y0 ; i <= y1 ; i++){

        draw(x0,i,c);
        c.r += colorConstantR;
        c.b += colorConstantB;
        c.g += colorConstantG;
    }

}

void draw12(int x0 , int x1 , int y0 , int y1, Color color0, Color color1){

    double colorConstantR = (color1.r-color0.r)/(y1-y0);
    double colorConstantG = (color1.g-color0.g)/(y1-y0);
    double colorConstantB = (color1.b-color0.b)/(y1-y0);

    Color c = color0;
    
    for(int i = y0 ; i >= y1 ; i--){

        draw(x0,i,c);
        c.r -= colorConstantR;
        c.b -= colorConstantB;
        c.g -= colorConstantG;
    }

}


int findMax(int a, int b, int c){
    if (a>b){
        if(a>c)
            return a;
        else
            return c;
    }
    else{
        if(b>c)
            return b;
        else
            return c;
    }

}


int findMin(int a, int b, int c){
    if (a<b){
        if(a<c)
            return a;
        else
            return c;
    }
    else{
        if(b<c)
            return b;
        else
            return c;
    }

}

double findArea(int x0,int x1,int x2,int y0,int y1,int y2){

    double areaOfT = (x0 * (y1-y2)) + (x1 * (y2-y0)) + (x2 * (y0-y1));
    areaOfT = areaOfT / 2;
    //areaOfT = abs(areaOfT);
    return areaOfT;

}

void paintForTriangle(int x0,int x1,int x2,int y0,int y1,int y2, Color c0, Color c1, Color c2){

    double alpha; 
    double beta; 
    double gamma;
    
    int xMax = findMax(x0,x1,x2);
    int xMin = findMin(x0,x1,x2);
    
    int yMax = findMax(y0,y1,y2);
    int yMin = findMin(y0,y1,y2);

    
    //Big triangle area
    double areaOfT = findArea(x0,x1,x2,y0,y1,y2);


    for(int y= yMin ; y <= yMax ; y++){
        
        for(int x = xMin ; x <= xMax; x++){

            double alpha = findArea(x,x1,x2,y,y1,y2);
            double beta = findArea(x0,x,x2,y0,y,y2);
            double gamma = findArea(x0,x1,x,y0,y1,y);

            alpha = alpha/areaOfT;
            beta = beta/areaOfT;
            gamma = gamma/areaOfT;

            Color tempColor;

            if(alpha >= 0 && gamma >= 0 && beta >= 0){
                tempColor.r = (int) ((alpha * c0.r) + (beta * c1.r) + (gamma * c2.r)) ;
                tempColor.g = (int) ((alpha * c0.g) + (beta * c1.g) + (gamma * c2.g)) ;
                tempColor.b = (int) ((alpha * c0.b) + (beta * c1.b) + (gamma * c2.b)) ;

                draw(x,y,tempColor); 
            }
         
        }
    }

}

void drawForTriangle(double x0,double x1,double y0,double y1, Color color0, Color color1){
    
    //// DRAW ALWAYS FROM LEFT TO RIGHT
    
    if(x1< x0){
        
        double tempx = x1;
        double tempy = y1;
        Color tempc = color1;
        x1 = x0;
        y1 = y0;
        color1 = color0;
        y0 = tempy;
        x0 = tempx;
        color0 = tempc;
    }
    
    double slope = (y1-y0) / (x1-x0); 


    // draw to north-east with the slope in range 0,1
    if (slope > 0.0 && slope <= 1.0)
        draw1(x0,x1,y0,y1,color0,color1);
    

    // draw to north-east with the slope in range 1,infinty
    else if (slope > 1.0)
        draw2(x0,x1,y0,y1,color0,color1);
    

    // draw to south-east with slope in range 0,-1
    else if (slope < 0.0 && slope > -1.0)
        draw3(x0,x1,y0,y1,color0,color1);
    

    // draw to south-east with slope in range -1, -infinity
    else if (slope < -1.0)
        draw4(x0,x1,y0,y1,color0,color1);
    

    // draw to right with slope = 0
    else
        draw5(x0,x1,y0,y1,color0,color1);

}

Vec3 findNormal(double r_v1[4],double r_v2[4],double r_v3[4]){

    double x0 = r_v1[0]; 
    double y0 = r_v1[1];
    double z0 = r_v1[2];

    double x1 = r_v2[0]; 
    double y1 = r_v2[1];
    double z1 = r_v2[2];

    double x2 = r_v3[0]; 
    double y2 = r_v3[1];
    double z2 = r_v3[2];

    Vec3 v1;
    v1.x = x1-x0;
    v1.y = y1-y0;
    v1.z = z1-z0;

    Vec3 v2;
    v2.x = x2-x0;
    v2.y = y2-y0;
    v2.z = z2-z0;

    // TO DO : control the order of crossProduct 
    
    Vec3 normal = crossProductVec3(normalizeVec3(v1),normalizeVec3(v2));
    return normal;

}

bool isCulled(double r_c[4], double r_v1[4],double r_v2[4],double r_v3[4], double Mcam[4][4], Vec3 camPos){

    Vec3 center;
    center.x = r_c[0];
    center.y = r_c[1];
    center.z = r_c[2];

    double camPosMatrix[4];
    camPosMatrix[0] = camPos.x;
    camPosMatrix[1] = camPos.y;
    camPosMatrix[2] = camPos.z;
    camPosMatrix[3] = 1.0;
    
    double McamPos[4];
    multiplyMatrixWithVec4d(McamPos,Mcam,camPosMatrix);
    perspectiveDivide(McamPos);

    // TO DO : control the P which position of Camera 
    Vec3 p;
    p.x = McamPos[0];
    p.y = McamPos[1];
    p.z = McamPos[2];

    Vec3 v;
    v.x = -(center.x - p.x);
    v.y = -(center.y - p.y);
    v.z = -(center.z - p.z);

    Vec3 normal = findNormal(r_v1, r_v2, r_v3);

    if (dotProductVec3(v,normal) >= 0) 
        return true;
    else
        return false;

}

void renderForEachmodel(Model model,std::vector<Mytriangle> &mytriangles,double Mvp[3][4],double Mper[4][4],double Mcam[4][4], Vec3 camPos){


    double UnitMatrix[4][4];
    makeIdentityMatrix(UnitMatrix);

    double Mmodel[4][4];
    makeIdentityMatrix(Mmodel);


    for (int i = 0; i < model.numberOfTransformations; i++){
        
        int id = model.transformationIDs[i];
        char type = model.transformationTypes[i];  
        
        if(type == 'r'){
            
            double R[4][4];
            double result[4][4];
            createRotationMatrix(rotations[id], R);

            multiplyMatrixWithMatrix(result, R, Mmodel);

            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    Mmodel[j][k] = result[j][k];
                     
        }

        else if(type == 't'){

            double T[4][4];
            double result[4][4];
            createTranslationMatrix(translations[id], T);

            multiplyMatrixWithMatrix(result, T, Mmodel);

            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    Mmodel[j][k] = result[j][k];
             
        }

        else{

            double S[4][4];
            double result[4][4];
            createScaleMatrix(scalings[id], S);

            multiplyMatrixWithMatrix(result, S, Mmodel);

            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    Mmodel[j][k] = result[j][k];
        }
    
    }

    double X[4][4];
    multiplyMatrixWithMatrix(X, Mmodel, UnitMatrix);

    double McamMmodel[4][4];
    multiplyMatrixWithMatrix(McamMmodel, Mcam, X);


    double MperMcamMmodel[4][4];
    multiplyMatrixWithMatrix(MperMcamMmodel, Mper, McamMmodel);


    for (int i = 0; i < mytriangles.size(); i++){
  
        double mper_r_c[4];
        double mper_r_v1[4];
        double mper_r_v2[4];
        double mper_r_v3[4];

        multiplyMatrixWithVec4d(mper_r_c, MperMcamMmodel, mytriangles[i].center);
        multiplyMatrixWithVec4d(mper_r_v1, MperMcamMmodel, mytriangles[i].vertex1);
        multiplyMatrixWithVec4d(mper_r_v2, MperMcamMmodel, mytriangles[i].vertex2);
        multiplyMatrixWithVec4d(mper_r_v3, MperMcamMmodel, mytriangles[i].vertex3);
        
        perspectiveDivide(mper_r_c);
        perspectiveDivide(mper_r_v1);
        perspectiveDivide(mper_r_v2);
        perspectiveDivide(mper_r_v3);

        double r_c[3];
        double r_v1[3];
        double r_v2[3];
        double r_v3[3];


        multiply_3_4_MatrixWithVec4d(r_c, Mvp, mper_r_c);
        multiply_3_4_MatrixWithVec4d(r_v1, Mvp, mper_r_v1);
        multiply_3_4_MatrixWithVec4d(r_v2, Mvp, mper_r_v2);
        multiply_3_4_MatrixWithVec4d(r_v3, Mvp, mper_r_v3);

        int v_id_0 = mytriangles[i].vertexIds[0];
        int v_id_1 = mytriangles[i].vertexIds[1];
        int v_id_2 = mytriangles[i].vertexIds[2];

        int v_colorId_0 = vertices[v_id_0].colorId;
        int v_colorId_1 = vertices[v_id_1].colorId;
        int v_colorId_2 = vertices[v_id_2].colorId;

        Color v0_color = colors[v_colorId_0];
        Color v1_color = colors[v_colorId_1];
        Color v2_color = colors[v_colorId_2];

        double x0 = (int) r_v1[0];
        double y0 = (int) r_v1[1];
          
        double x1 = (int) r_v2[0];
        double y1 = (int) r_v2[1];

        double x2 = (int) r_v3[0];
        double y2 = (int) r_v3[1];

        bool isTriangleCulled = false;

        if(backfaceCullingSetting){

            isTriangleCulled = isCulled(r_c,r_v1,r_v2,r_v3, Mcam, camPos);
            
            // If culling enable then look for isTriangleCulled if so then do nothing : continue; 
            if(isTriangleCulled)
                continue;
        }
        if (!model.type){

            drawForTriangle(x0,x1,y0,y1,v0_color,v1_color);
            drawForTriangle(x0,x2,y0,y2,v0_color,v2_color);
            drawForTriangle(x1,x2,y1,y2,v1_color,v2_color);

        }
        
          
        if(model.type){
            // Do painting for solids
            paintForTriangle(x0,x1,x2,y0,y1,y2,v0_color,v1_color,v2_color);
        }

    }
    
}

void initializeImage(Camera cam) {
    int i, j;

    for (i = 0; i < cam.sizeX; i++)
        for (j = 0; j < cam.sizeY; j++) {
            image[i][j].r = backgroundColor.r;
            image[i][j].g = backgroundColor.g;
            image[i][j].b = backgroundColor.b;

        }

}

void forwardRenderingPipeline(Camera cam){

    double Mper[4][4];
    double Mcam[4][4];
    double Mvp[3][4];

    createMper(cam, Mper);
    createMcam(cam, Mcam);
    createMvp(cam, Mvp);

    for(int i = 0 ; i < numberOfModels ; i++){

        std::vector<Mytriangle> mytriangles;
        for (int j = 0; j < models[i].numberOfTriangles; j++){

            Mytriangle t = findMytriangle(models[i].triangles[j]);
            mytriangles.push_back(t);
        }

        renderForEachmodel(models[i], mytriangles, Mvp, Mper, Mcam, cam.pos);
    }

}

int main(int argc, char **argv) {
    int i, j;

    if (argc < 2) {
        std::cout << "Usage: ./rasterizer <scene file> <camera file>" << std::endl;
        return 1;
    }

    // read camera and scene files
    readSceneFile(argv[1]);
    readCameraFile(argv[2]);

    image = 0;

    for (i = 0; i < numberOfCameras; i++) {

        // allocate memory for image
        if (image) {
            for (j = 0; j < cameras[i].sizeX; j++) {
                delete image[j];
            }

            delete[] image;
        }

        image = new Color*[cameras[i].sizeX];

        if (image == NULL) {
            std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
            exit(1);
        }

        for (j = 0; j < cameras[i].sizeX; j++) {
            image[j] = new Color[cameras[i].sizeY];
            if (image[j] == NULL) {
                std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
                exit(1);
            }
        }


        // initialize image with basic values
        initializeImage(cameras[i]);

        // do forward rendering pipeline operations
        forwardRenderingPipeline(cameras[i]);

        // generate PPM file
        writeImageToPPMFile(cameras[i]);

        // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
        // Notice that os_type is not given as 1 (Ubuntu) or 2 (Windows), below call doesn't do conversion.
        // Change os_type to 1 or 2, after being sure that you have ImageMagick installed.
        convertPPMToPNG(cameras[i].outputFileName, 99);
    }

    return 0;

}