#include "getopt.h"

#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <chrono>
#include <thread>
#include <regex>
#include <atomic>
#include <future>

#include <cstdio>
#include <csignal>

#include <chrono>
#include <thread>
#include <iostream>

#include "TelGL.hh"
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "linenoise.h"
#include "myrapidjson.h"

static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -data     [PATH]   path to data file (input)
  -geometry [PATH]   path to geometry file (input)
)";

int main(int argc, char **argv){
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] =
    {
     { "help",       no_argument,       &do_help,      1  },
     { "verbose",    no_argument,       &do_verbose,   1  },
     { "data",       required_argument, NULL,           'd' },
     { "geometry",   required_argument, NULL,           'g' },
     { 0, 0, 0, 0 }};

  std::string data_path;
  std::string geometry_path;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
    switch (c) {
    case 'd':
      data_path = optarg;
      break;
    case 'g':
      geometry_path = optarg;
      break;
      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      std::fprintf(stderr,"case 1\n");
      exit(1);
      break;
    case ':':
      std::fprintf(stderr,"case :\n");
      exit(1);
      break;
    case '?':
      std::fprintf(stderr,"case ?\n");
      exit(1);
      break;
    default:
      std::fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }

  if(do_help){
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    exit(0);
  }

  if(geometry_path.empty() || data_path.empty()){
    std::fprintf(stdout, "Error: input option.\n");
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    std::exit(1);
  }


  //////////////////////////////////////////
  glfwSetErrorCallback([](int error, const char* description){
                         fprintf(stderr, "Error: %s\n", description);});
  if (!glfwInit()){
    std::exit(1);
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_SAMPLES, 4);
  float sWinWidth = 1000;
  float sWinHeight = 480;
  GLFWwindow* window = glfwCreateWindow((int)sWinWidth, (int)sWinHeight, "OpenGL telescope example", NULL, NULL);
  if (!window){
    glfwTerminate();
    std::exit(0);
  }
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int width, int height){
                                           glViewport(0, 0, width, height);});
  ////^^^/////////////////////////////////

  std::string jsstr_geo = altel::TelGL::readFile(geometry_path);
  JsonDocument jsd_geo = altel::TelGL::createJsonDocument(jsstr_geo);
  altel::TelGL telgl(jsd_geo);

  float deltaTime = 0.0f; // time between current frame and last frame
  float lastFrame = 0.0f;
  float lastWinRatio = sWinWidth/sWinHeight;
  JsonDocument jsd_trans = altel::TelGL::createTransformExample();
  glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f,  -1000.0f);
  glm::vec3 worldCenter = glm::vec3(0.0f, 0.0f,  0.0f);
  glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);


  while (!glfwWindowShouldClose(window)){
    int width;
    int height;
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;
    glfwGetFramebufferSize(window, &width, &height);
    float currentWinRatio = width / (float) height;
    {
      if(glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        continue;
      if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
      if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        cameraPos -= (cameraPos - worldCenter) * (deltaTime * 1.0f);
      if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        cameraPos += (cameraPos - worldCenter) * (deltaTime * 1.0f);
      if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::cross((cameraPos - worldCenter), glm::vec3(0.0f, 1.0f,  0.0f)) * (deltaTime * 2.0f);
      if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::cross((cameraPos - worldCenter), glm::vec3(0.0f, 1.0f,  0.0f)) * (deltaTime * 2.0f);
      if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos -= glm::cross((cameraPos - worldCenter), glm::vec3(1.0f, 0.0f,  0.0f)) * (deltaTime * 2.0f);
      if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos += glm::cross((cameraPos - worldCenter), glm::vec3(1.0f, 0.0f,  0.0f)) * (deltaTime * 2.0f);

      jsd_trans["trans"]["lookat"]["eye"]["x"]= cameraPos[0];
      jsd_trans["trans"]["lookat"]["eye"]["y"]= cameraPos[1];
      jsd_trans["trans"]["lookat"]["eye"]["z"]= cameraPos[2];
      jsd_trans["trans"]["lookat"]["center"]["x"]= worldCenter[0];
      jsd_trans["trans"]["lookat"]["center"]["y"]= worldCenter[1];
      jsd_trans["trans"]["lookat"]["center"]["z"]= worldCenter[2];
      jsd_trans["trans"]["lookat"]["up"]["x"]= cameraUp[0];
      jsd_trans["trans"]["lookat"]["up"]["y"]= cameraUp[1];
      jsd_trans["trans"]["lookat"]["up"]["z"]= cameraUp[2];
      jsd_trans["trans"]["persp"]["fov"]= 60;
      jsd_trans["trans"]["persp"]["ratio"]= currentWinRatio;
      jsd_trans["trans"]["persp"]["near"]= 0.1;
      jsd_trans["trans"]["persp"]["far"]= 2000;
      telgl.updateTransform(jsd_trans);
    }

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    {
      std::string jsstr_data = altel::TelGL::readFile(data_path);
      JsonDocument jsd_data = altel::TelGL::createJsonDocument(jsstr_data);

      if(jsd_data.HasMember("tracks")){
        telgl.drawTracks(jsd_data);
      }
      if(jsd_data.HasMember("hits")){
        telgl.drawHits(jsd_data);
      }
      if(jsd_data.HasMember("detectors")){
        telgl.drawDetectors(jsd_data);
      }
      else{
        telgl.drawDetectors();
      }
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
    lastWinRatio=currentWinRatio;
  }

  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
