// IMPORTANT COMMENT
// In order to compile the code please add glm folder to
// where source codes exist. 
// glm folder is same as shared glm lib codes by teaching asistant

#include "helper.h"
#include <math.h>
#include <vector>
#include "./glm/glm.hpp" // GL Math library header
#include "./glm/gtc/matrix_transform.hpp"
#include "./glm/gtc/type_ptr.hpp"
#include <GL/glu.h>
#include <GL/glut.h>
static GLFWwindow * win = NULL;

int widthTexture, heightTexture;

// Constant Key settings
int keys[1024];

// Constant vertex and index vectors
std::vector<int> indicesVector;
std::vector<glm::vec3> verticesVector;

// Increase or Decrease with O-L
float hFactor = 10.0;

// change with F
bool fullscreen = false;

// Shaders
GLuint idProgramShader;
GLuint idFragmentShader;
GLuint idVertexShader;
GLuint idJpegTexture;
GLuint idMVPMatrix;

struct Camera{
  float fov;
  float nearDistance;
  float farDistance;
  float speed;
  float rotateAngle;
  glm::vec3 cameraPosition;
  glm::vec3 cameraUp;
  glm::vec3 cameraGaze;
  glm::vec3 front;
  glm::vec3 left;
  float verticalAngle;
  float horizantalAngle;
};

Camera cam;


static void errorCallback(int error,
  const char * description) {
  fprintf(stderr, "Error: %s\n", description);

}


void createVertexData(){

    for (int row = 0; row < heightTexture; row++){
      for (int col = 0; col < widthTexture; col++){


      float xpos = col;
      float ypos = 0.0f;
      float zpos = row;


      verticesVector.push_back(glm::vec3(xpos, ypos, zpos));
    }
  }



}


void createIndexData(){

    for (int row = 0; row < heightTexture-1; row++){

      for (int col = 0; col < widthTexture-1; col++){
      
      int topLeftIndexNum = (row*widthTexture + col );
      int topRightIndexNum = (row*widthTexture + col  + 1);
      int bottomLeftIndexNum = ((row+1) * widthTexture  + col);
      int bottomRightIndexNum = ((row+1)*widthTexture + col + 1);

      indicesVector.push_back(topLeftIndexNum);
      indicesVector.push_back(bottomLeftIndexNum);
      indicesVector.push_back(topRightIndexNum);

      indicesVector.push_back(topRightIndexNum);
      indicesVector.push_back(bottomLeftIndexNum);
      indicesVector.push_back(bottomRightIndexNum);



    }
  }



}

void setUniforms(glm::mat4 projection, glm::mat4 view, glm::vec3 camPos){

  float w = (float) widthTexture;
  float h = (float) heightTexture;

  glm::vec4 cameraPosition = glm::vec4(camPos, 1.0f);

  glUniformMatrix4fv(glGetUniformLocation(idProgramShader,"MVP"),1,GL_FALSE, glm::value_ptr(projection));
  glUniformMatrix4fv(glGetUniformLocation(idProgramShader,"MV"),1,GL_FALSE, glm::value_ptr(view));
  glUniform4fv(glGetUniformLocation(idProgramShader,"cameraPosition"),1, glm::value_ptr(cameraPosition));



  glUniform1f(glGetUniformLocation(idProgramShader,"heightFactor"), hFactor);
  glUniform1f(glGetUniformLocation(idProgramShader,"wFloat"), w);
  glUniform1f(glGetUniformLocation(idProgramShader,"hFloat"), h);
  
}

void initializeCamera(Camera &cam){
  
  cam.fov = 45.0f;
  cam.nearDistance = 0.1f;
  cam.farDistance = 1000.0f;
  cam.speed = 0.0f;
  cam.rotateAngle = 1.0f;
  cam.cameraPosition = glm::vec3(widthTexture/2.0f, widthTexture/10.0f, -widthTexture/4.0f);
  cam.cameraUp = glm::vec3(0.0f,1.0f,0.0f);
  cam.cameraGaze = glm::vec3(0.0f,0.0f,1.0f);
  cam.front = glm::vec3(widthTexture/2.0f, widthTexture/10.0f, -widthTexture/4.0f+0.1f);
  cam.left=glm::normalize(glm::cross(cam.cameraUp,cam.cameraGaze));
  //glViewport(0,0,widthTexture,heightTexture);
}

void setCamera(Camera &cam){

  //glViewport(0,0,widthTexture,heightTexture);  
  cam.cameraPosition = cam.cameraPosition + (cam.speed*cam.cameraGaze);
  cam.front = cam.cameraPosition + (cam.nearDistance * cam.cameraGaze);

  glm::mat4 projection = glm::perspective(glm::radians(cam.fov), 1.0f, cam.nearDistance, cam.farDistance);
  glm::mat4 view = glm::lookAt(cam.cameraPosition,cam.front,cam.cameraUp);

  setUniforms(projection, view, cam.cameraPosition);

}

void bindVertexIndexData(){

  unsigned int VBO, VAO, EBO;

  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glGenBuffers(1, &EBO);

  glBindVertexArray(VAO);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, verticesVector.size()*sizeof(glm::vec3), &verticesVector.front(), GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);


  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indicesVector.size()*sizeof(int), &indicesVector.front(), GL_STATIC_DRAW);


  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);

  glBindVertexArray(VAO); 

}

void drawData(){

    glDrawElements(GL_TRIANGLES, indicesVector.size(), GL_UNSIGNED_INT, 0);

}

void goRight(){

  glm::mat4 R(1);
  R = glm::rotate(R,glm::radians(-cam.rotateAngle),cam.cameraUp);
  cam.cameraGaze = glm::normalize(glm::vec3(R*glm::vec4(cam.cameraGaze,1.0f)));
  cam.left = glm::normalize(glm::cross(cam.cameraUp,cam.cameraGaze));


}

void goLeft(){

  glm::mat4 R(1);
  R = glm::rotate(R,glm::radians(cam.rotateAngle),cam.cameraUp);
  cam.cameraGaze = glm::normalize(glm::vec3(R*glm::vec4(cam.cameraGaze,1.0f)));
  cam.left = glm::normalize(glm::cross(cam.cameraUp,cam.cameraGaze));



}

void goForward(){

  glm::mat4 R(1);
  R = glm::rotate(R,glm::radians(-cam.rotateAngle),cam.left);
  cam.cameraGaze = glm::normalize(glm::vec3(R*glm::vec4(cam.cameraGaze,1.0f)));
  cam.cameraUp = glm::normalize(glm::cross(cam.cameraGaze,cam.left));


}

void goBackward(){

  glm::mat4 R(1);
  R = glm::rotate(R,glm::radians(cam.rotateAngle),cam.left);
  cam.cameraGaze = glm::normalize(glm::vec3(R*glm::vec4(cam.cameraGaze,1.0f)));
  cam.cameraUp = glm::normalize(glm::cross(cam.cameraGaze,cam.left));

    
}

void setFullScreen(GLFWwindow* window){

  glfwSetWindowMonitor(window, glfwGetPrimaryMonitor(), 0, 0, glfwGetVideoMode(glfwGetPrimaryMonitor())->width, glfwGetVideoMode(glfwGetPrimaryMonitor())->height, 1);

}


void setWindowScreen(GLFWwindow* window){

  glfwSetWindowMonitor(window, NULL, 100, 100, 600, 600, 1);

}

void screenChange(){
    
  if(fullscreen){
    fullscreen = false;
    setWindowScreen(win);
  }
  
  else{
    fullscreen = true;
    setFullScreen(win);
  }


}

void checkEvent(Camera &cam){

  if(keys[GLFW_KEY_O]){
        
        keys[GLFW_KEY_O] = 0;
        hFactor += 0.5;
  }

  else if(keys[GLFW_KEY_L]){

        keys[GLFW_KEY_L] = 0;
        hFactor -= 0.5;
  }

  else if (keys[GLFW_KEY_W]){

        keys[GLFW_KEY_W] = 0;
        goForward();
        
  }

  else if (keys[GLFW_KEY_S]){
        keys[GLFW_KEY_S] = 0;
        goBackward();
  }

  else if (keys[GLFW_KEY_A]){

        keys[GLFW_KEY_A] = 0;
        goLeft();
  }

  else if (keys[GLFW_KEY_D]){

        keys[GLFW_KEY_D] = 0;
        goRight();
  }

  else if (keys[GLFW_KEY_U]){

        keys[GLFW_KEY_U] = 0;
        cam.speed += 0.1;
  }

  else if (keys[GLFW_KEY_J]){

        keys[GLFW_KEY_J] = 0;
        cam.speed -= 0.1;
  }

  else if (keys[GLFW_KEY_F]){

        keys[GLFW_KEY_F] = 0;
        screenChange();

  }



}



void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode){
  
  
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
     glfwSetWindowShouldClose(window, GL_TRUE);
  
  if (key >= 0 && key < 1024){
    
    if (action == GLFW_PRESS)
      keys[key] = 1;
    
    else if (action == GLFW_RELEASE)
      keys[key] = 0;
  }



}


void framebuffer_size_callback(GLFWwindow* window, int width, int height){

    glViewport(0, 0, width, height);
}


int main(int argc, char * argv[]) {

  if (argc != 2) {
    printf("Only one texture image expected!\n");
    exit(-1);
  }

  glfwSetErrorCallback(errorCallback);

  if (!glfwInit()) {
    exit(-1);
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

  win = glfwCreateWindow(600, 600, "CENG477 - HW3", NULL, NULL);

  if (!win) {
    glfwTerminate();
    exit(-1);
  }
  glfwMakeContextCurrent(win);
  glfwSetFramebufferSizeCallback(win, framebuffer_size_callback);
  
  GLenum err = glewInit();
  if (err != GLEW_OK) {
    fprintf(stderr, "Error: %s\n", glewGetErrorString(err));

    glfwTerminate();
    exit(-1);
  }


  initShaders();
  glUseProgram(idProgramShader);
  initTexture(argv[1], & widthTexture, & heightTexture);

  createVertexData();
  createIndexData();
  bindVertexIndexData();
  glfwSetKeyCallback(win, key_callback);

  glViewport(0,0,600,600);

  //glfwSetInputMode(win, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
  //glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

  initializeCamera(cam);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glFrontFace(GL_CCW);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  while (!glfwWindowShouldClose(win)) {
    
    //Check if any event ocuured or not (Press Key etc.)
    checkEvent(cam);
    glfwPollEvents();

    //Clear Screen
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    //Set Camera
    setCamera(cam);
  
    //Render Data
    drawData();

    //Swap Buffer
    glfwSwapBuffers(win);

}

  glfwDestroyWindow(win);
  glfwTerminate();

  return 0;

}



