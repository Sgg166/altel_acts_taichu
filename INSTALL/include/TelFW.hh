#pragma once

#include <cstdio>
#include <csignal>

#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <chrono>
#include <thread>
#include <regex>
#include <iostream>
#include <future>

#include <unistd.h>

#include "TelGL.hh"

class GLFWwindow;
class TelFW{
 public:
  TelFW(float sWinWidth, float sWinHeight, const std::string& title);
  ~TelFW();
  GLFWwindow* get();
  int stopAsync();

  template<typename T>
  bool startAsync(T* t,
                  std::function<int(T*, GLFWwindow*)> beginHook,
                  std::function<int(T*, GLFWwindow*)> clearHook,
                  std::function<int(T*, GLFWwindow*)> drawHook
                  ){
    if(!m_fut.valid()){
      m_fut = std::async(std::launch::async, &TelFW::loopFun<T>,
                         this, t,  beginHook, clearHook, drawHook);
    }
    return true;
  }
  template<typename T>
  int loopFun(T* t,
              std::function<int(T*, GLFWwindow*)> beginHook,
              std::function<int(T*, GLFWwindow*)> clearHook,
              std::function<int(T*, GLFWwindow*)> drawHook);
  static void initializeGLFW();
  static void terminateGLFW();
  static GLFWwindow* createWindow(float sWinWidth, float sWinHeight, const std::string& title);
  static void deleteWindow(GLFWwindow* window);
  // static getFramebufferSize(GLFWwindow* window, int* width,  int* height);
  // static getKey(GLFWwindow* window, int);
private:
  bool checkifclose();
  void clear();
  void flush();
private:
  int m_width;
  int m_height;
  std::string m_title;
  bool m_runflag;
  std::future<int> m_fut;
  GLFWwindow* m_win;
  static size_t s_count;
};

template<typename T>
int TelFW::loopFun(T* t,
                   std::function<int(T*, GLFWwindow*)> beginHook,
                   std::function<int(T*, GLFWwindow*)> clearHook,
                   std::function<int(T*, GLFWwindow*)> drawHook){
  m_runflag = true;
  if(!s_count){
    initializeGLFW();
  }
  m_win = createWindow(m_width, m_height, m_title.c_str());
  s_count++;

  beginHook(t, m_win);
  while (!checkifclose()){
    // std::cout<< ".";
    int ch = clearHook(t, m_win);
    if(ch==0) continue;
    if(ch>0){
    clear();
    }
    int dh = drawHook(t, m_win);
    flush();// todo: split flush and poll
  }

  deleteWindow(m_win);
  s_count--;
  m_runflag = false;
  std::cout<< "returned loop"<<std::endl;
  return 0;
}
