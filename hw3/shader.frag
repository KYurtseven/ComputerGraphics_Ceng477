#version 410

// Output Color
out vec4 color;

uniform mat4 MVP; // ModelViewProjection Matrix
uniform mat4 MV; // ModelView idMVPMatrix
uniform vec4 cameraPosition;

// Texture-related data;
uniform sampler2D rgbTexture;
uniform int widthTexture;
uniform int heightTexture;

// Data from Vertex Shader
in vec2 textureCoordinate;
in vec3 vertexNormal; // For Lighting computation
in vec3 ToLightVector; // Vector from Vertex to Light;
in vec3 ToCameraVector; // Vector from Vertex to Camera;

void main() {

  // Assignment Constants below
  // get the texture color
  


  vec4 textureColor = texture(rgbTexture,textureCoordinate);

  // apply Phong shading by using the following parameters
  vec4 ka = vec4(0.25f,0.25f,0.25f,1.0f); // reflectance coeff. for ambient
  vec4 Ia = vec4(0.3f,0.3f,0.3f,1.0f); // light color for ambient
  vec4 Id = vec4(1.0f, 1.0f, 1.0f, 1.0f); // light color for diffuse
  vec4 kd = vec4(1.0f, 1.0f, 1.0f, 1.0f); // reflectance coeff. for diffuse
  vec4 Is = vec4(1.0f, 1.0f, 1.0f, 1.0f); // light color for specular
  vec4 ks = vec4(1.0f, 1.0f, 1.0f, 1.0f); // reflectance coeff. for specular
  int specExp = 100; // specular exponent

  // compute ambient component
  vec4 ambient = ka*Ia;
  
  
  // compute diffuse component
  float cos;
  cos = dot(ToLightVector,vertexNormal);
  cos = max(0.0, cos);
  
  vec4 diffuse = (kd*cos*Id);


  // compute specular component
  vec3 h = ToLightVector + ToCameraVector;
  h = normalize(h);
  
  float cos2 = dot(vertexNormal, h);
  cos2 = max(cos2,0.0);
  cos2 = pow(cos2,specExp);
  
  vec4 specular = (ks*cos2*Is);

  color = vec4(clamp( textureColor.xyz * vec3(ambient + diffuse + specular), 0.0, 1.0), 1.0);

}
