#version 410
layout(location = 0) in vec3 position;

// Data from CPU 
uniform mat4 MVP; // ModelViewProjection Matrix
uniform mat4 MV; // ModelView idMVPMatrix
uniform vec4 cameraPosition;
uniform float heightFactor;

// Texture-related data
uniform sampler2D rgbTexture;
uniform int widthTexture;
uniform int heightTexture;

uniform float wFloat;
uniform float hFloat;


// Output to Fragment Shader
out vec2 textureCoordinate; // For texture-color
out vec3 vertexNormal; // For Lighting computation
out vec3 ToLightVector; // Vector from Vertex to Light;
out vec3 ToCameraVector; // Vector from Vertex to Camera;


void findVertex(in float px, in float pz, out vec3 vertex){
    
    vec2 txtCoor;
    txtCoor.x = -px/wFloat;
    txtCoor.y = -pz/hFloat;

    vec4 txtCol = texture(rgbTexture, txtCoor);
    float yValue = heightFactor*((0.2126*txtCol.x) + (0.7152*txtCol.y) + (0.0722*txtCol.z));
    vertex = vec3(px,yValue,pz);

}


void findNormalofTriangle(in vec3 v1 , in vec3 v2, in vec3 v3 , out vec3 n){
    
    vec3 u = v2 - v1;
    vec3 w = v3 - v1;

    n = cross(u,w);
    n = normalize(n);

}


void findLightVectors(){


      
  // Light position calculations 
  
  vec3 Ipos = vec3(wFloat/2.0, wFloat+hFloat , hFloat/2.0);
  vec3 cpos = vec3(cameraPosition.x, cameraPosition.y, cameraPosition.z);
  vec3 mypos;
  findVertex(position.x,position.z, mypos);
  

  ToLightVector = normalize(Ipos - mypos);
  ToCameraVector = normalize(cpos - mypos);

  

}



void main(){

  
  findLightVectors();


  // H calculations

  textureCoordinate.x = -position.x/wFloat;
  textureCoordinate.y = -position.z/hFloat;

  vec4 textureColor = texture(rgbTexture, textureCoordinate);


  float yValue = heightFactor*((0.2126*textureColor.x) + (0.7152*textureColor.y) + (0.0722*textureColor.z));

  vec4 model = vec4(position.x,yValue,position.z,1.0);
  model = MV * model ;
  model = MVP * model;


  gl_Position =  vec4(model.x, model.y, model.z, model.w);

  
  // Normal Calculations

  vertexNormal = vec3(0.0,0.0,0.0);

  if(position.z == 0.0){

    if(position.x == 0.0){
      // Find normal of triangle a

      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x,position.z+1.0, v2);
      vec3 v3;
      findVertex(position.x+1.0, position.z, v3);

      vec3 normal;
      findNormalofTriangle(v1,v2,v3,normal);

      normal = normalize(normal);
      vertexNormal = normal;

    }


    else if(position.x == wFloat){
      // Find normals of triangles x,y

      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x-1.0f , position.z, v2);
      vec3 v3;
      findVertex(position.x-1.0f, position.z+1.0f, v3);
      vec3 v4;
      findVertex(position.x, position.z+1.0f , v4);

      vec3 normal1;
      findNormalofTriangle(v1,v2,v3,normal1);
      vec3 normal2;
      findNormalofTriangle(v1,v3,v4,normal2);

      vec3 sumOfNormals;
      sumOfNormals = normal1+normal2;
      sumOfNormals = normalize(sumOfNormals);

      vertexNormal = sumOfNormals;

    }


    else{
      // Find normals of triangles (Example : c,d,e for vertice 2)

      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x-1.0f , position.z, v2);
      vec3 v3;
      findVertex(position.x-1.0f, position.z+1.0f, v3);
      vec3 v4;
      findVertex(position.x, position.z+1.0f , v4);
      vec3 v5;
      findVertex(position.x+1.0f, position.z , v5);

      vec3 normal1;
      findNormalofTriangle(v1,v2,v3,normal1);
      vec3 normal2;
      findNormalofTriangle(v1,v3,v4,normal2);
      vec3 normal3;
      findNormalofTriangle(v1,v4,v5,normal3);

      vec3 sumOfNormals;
      sumOfNormals = normal1+normal2+normal3;
      sumOfNormals = normalize(sumOfNormals);

      vertexNormal = sumOfNormals;


    }
  }


  else if(position.z == hFloat){
    

    if(position.x == 0.0){

      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x+1.0f , position.z, v2);
      vec3 v3;
      findVertex(position.x+1.0f, position.z-1.0f, v3);
      vec3 v4;
      findVertex(position.x, position.z - 1.0f , v4);

      vec3 normal1;
      findNormalofTriangle(v1,v2,v3,normal1);
      vec3 normal2;
      findNormalofTriangle(v1,v3,v4,normal2);

      vec3 sumOfNormals;
      sumOfNormals = normal1+normal2;
      sumOfNormals = normalize(sumOfNormals);

      vertexNormal = sumOfNormals;

      

    }
    else if(position.x == wFloat){
      
      // Find normal of triangle sigma

      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x,position.z-1.0, v2);
      vec3 v3;
      findVertex(position.x-1.0, position.z, v3);

      vec3 normal;
      findNormalofTriangle(v1,v2,v3,normal);

      normal = normalize(normal);
      vertexNormal = normal;



    }

    else{
      // Find normals of triangles (Example : beta,q,p for t)


      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x+1.0f , position.z, v2);
      vec3 v3;
      findVertex(position.x+1.0f, position.z-1.0f, v3);
      vec3 v4;
      findVertex(position.x, position.z - 1.0f , v4);
      vec3 v5;
      findVertex(position.x - 1.0f, position.z , v5);

      vec3 normal1;
      findNormalofTriangle(v1,v2,v3,normal1);
      vec3 normal2;
      findNormalofTriangle(v1,v3,v4,normal2);
      vec3 normal3;
      findNormalofTriangle(v1,v4,v5,normal3);

      vec3 sumOfNormals;
      sumOfNormals = normal1+normal2+normal3;
      sumOfNormals = normalize(sumOfNormals);

      vertexNormal = sumOfNormals;

    }
  }

  else if(position.x == 0.0){
    // Find normals of (Example : a,b,g for vertice 1)


      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x , position.z + 1.0f, v2);
      vec3 v3;
      findVertex(position.x+1.0f, position.z, v3);
      vec3 v4;
      findVertex(position.x +1.0f, position.z - 1.0f , v4);
      vec3 v5;
      findVertex(position.x, position.z -1.0f , v5);

      vec3 normal1;
      findNormalofTriangle(v1,v2,v3,normal1);
      vec3 normal2;
      findNormalofTriangle(v1,v3,v4,normal2);
      vec3 normal3;
      findNormalofTriangle(v1,v4,v5,normal3);

      vec3 sumOfNormals;
      sumOfNormals = normal1+normal2+normal3;
      sumOfNormals = normalize(sumOfNormals);

      vertexNormal = sumOfNormals;


  }

  else if(position.x == wFloat){
    //Find normals of (Example: y,w,z for s)

      vec3 v1;
      findVertex(position.x, position.z, v1);
      vec3 v2;
      findVertex(position.x , position.z - 1.0f, v2);
      vec3 v3;
      findVertex(position.x - 1.0f, position.z, v3);
      vec3 v4;
      findVertex(position.x - 1.0f, position.z + 1.0f , v4);
      vec3 v5;
      findVertex(position.x, position.z + 1.0f , v5);

      vec3 normal1;
      findNormalofTriangle(v1,v2,v3,normal1);
      vec3 normal2;
      findNormalofTriangle(v1,v3,v4,normal2);
      vec3 normal3;
      findNormalofTriangle(v1,v4,v5,normal3);

      vec3 sumOfNormals;
      sumOfNormals = normal1+normal2+normal3;
      sumOfNormals = normalize(sumOfNormals);

      vertexNormal = sumOfNormals;

  }

  else{
    
    // Find normals of I,II,III,IV,V,VI



      // IV
      vec3 v1_1;
      findVertex(position.x, position.z, v1_1);
      vec3 v1_2;
      findVertex(position.x,position.z+1.0, v1_2);
      vec3 v1_3;
      findVertex(position.x+1.0, position.z, v1_3);

      vec3 normal1;
      findNormalofTriangle(v1_1,v1_2,v1_3,normal1);

      
      //III
      vec3 v2_1;
      findVertex(position.x, position.z, v2_1);
      vec3 v2_2;
      findVertex(position.x+1.0,position.z, v2_2);
      vec3 v2_3;
      findVertex(position.x+1.0, position.z-1.0, v2_3);

      vec3 normal2;
      findNormalofTriangle(v2_1,v2_2,v2_3,normal2);


      //II
      vec3 v3_1;
      findVertex(position.x, position.z, v3_1);
      vec3 v3_2;
      findVertex(position.x+1.0,position.z-1.0, v3_2);
      vec3 v3_3;
      findVertex(position.x, position.z-1.0, v3_3);

      vec3 normal3;
      findNormalofTriangle(v3_1,v3_2,v3_3,normal3);


      //I
      vec3 v4_1;
      findVertex(position.x, position.z, v4_1);
      vec3 v4_2;
      findVertex(position.x,position.z-1.0, v4_2);
      vec3 v4_3;
      findVertex(position.x-1.0, position.z, v4_3);

      vec3 normal4;
      findNormalofTriangle(v4_1,v4_2,v4_3,normal4);


      //VI
      vec3 v5_1;
      findVertex(position.x, position.z, v5_1);
      vec3 v5_2;
      findVertex(position.x-1.0,position.z, v5_2);
      vec3 v5_3;
      findVertex(position.x-1.0, position.z+1.0, v5_3);

      vec3 normal5;
      findNormalofTriangle(v5_1,v5_2,v5_3,normal5);


      vec3 v6_1;
      findVertex(position.x, position.z, v6_1);
      vec3 v6_2;
      findVertex(position.x-1.0,position.z+1.0, v6_2);
      vec3 v6_3;
      findVertex(position.x, position.z+1.0, v6_3);

      vec3 normal6;
      findNormalofTriangle(v6_1,v6_2,v6_3,normal6);

      vec3 sumOfNormals = normal1+normal2+normal3+normal4+normal5+normal6;
      sumOfNormals = normalize(sumOfNormals);

      vertexNormal = sumOfNormals;





  }
    
}
