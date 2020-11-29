#version 430 core

layout(points) in;
layout(line_strip, max_vertices = 16) out;

layout (std140) struct GeoDataLayer{
  ivec4 id; // id, rev , rev, rev
  vec4  pos; // vec3, pad
  vec4  size;     // size_x, size_y, size_z, pad
  vec4  color;    // vec3, pad
  vec4  pitch;    // pitch_x pitch_y pitch_z/thick_z
  uvec4 npixel;   // pixel_x pixel_y pixel_z/always1
  mat4  trans;    // model
};

layout (std140) uniform GeoData{
  GeoDataLayer  ly[20]; //
} geoData;

layout (std140) uniform TransformData{
  mat4 model;
  mat4 view;
  mat4 proj;
} transformData;


in int nl[];
out vec3 fColor;

void main(){
  vec3   pos;
  vec3   thick;
  vec3   color;
  vec3   pitch;
  uvec3  npixel;
  mat4   trans;
  int    n;

  n = -1;
  for (int k = 0; k < 20; ++k){
    if(nl[0] == geoData.ly[k].id.x ){
      n = k;
      break;
    }
  }
  if(n == -1){
    return ;
  }
  for (int k = 0; k < 20; ++k){
    if(n ==k){
      pos    = geoData.ly[k].pos.xyz;
      thick  = geoData.ly[k].size.xyz;
      color  = geoData.ly[k].color.xyz;
      pitch  = geoData.ly[k].pitch.xyz;
      npixel = geoData.ly[k].npixel.xyz;
      trans  = geoData.ly[k].trans;
      break;
    }
  }
  fColor = color;

  mat4 pvmMatrix = transformData.proj * transformData.view * transformData.model * trans;

  vec4 gp = pvmMatrix * vec4(pos, 1.0);

  vec3 halfThick = thick/2.0;
  //                                     N/S       U/D       E/W
  //                                     X         Y         Z
  vec4 NEU = pvmMatrix * vec4(  halfThick.x,  halfThick.y,  halfThick.z, 0.0);
  vec4 NED = pvmMatrix * vec4(  halfThick.x, -halfThick.y,  halfThick.z, 0.0);
  vec4 NWU = pvmMatrix * vec4(  halfThick.x,  halfThick.y, -halfThick.z, 0.0);
  vec4 NWD = pvmMatrix * vec4(  halfThick.x, -halfThick.y, -halfThick.z, 0.0);
  vec4 SEU = pvmMatrix * vec4( -halfThick.x,  halfThick.y,  halfThick.z, 0.0);
  vec4 SED = pvmMatrix * vec4( -halfThick.x, -halfThick.y,  halfThick.z, 0.0);
  vec4 SWU = pvmMatrix * vec4( -halfThick.x,  halfThick.y, -halfThick.z, 0.0);
  vec4 SWD = pvmMatrix * vec4( -halfThick.x, -halfThick.y, -halfThick.z, 0.0);
  gl_Position = gp + NED;
  EmitVertex();
  gl_Position = gp + NWD;
  EmitVertex();
  gl_Position = gp + SWD;
  EmitVertex();
  gl_Position = gp + SED;
  EmitVertex();
  gl_Position = gp + SEU;
  EmitVertex();
  gl_Position = gp + SWU;
  EmitVertex();
  gl_Position = gp + NWU;
  EmitVertex();
  gl_Position = gp + NEU;
  EmitVertex();
  gl_Position = gp + NED;
  EmitVertex();
  gl_Position = gp + SED;
  EmitVertex();
  gl_Position = gp + SEU;
  EmitVertex();
  gl_Position = gp + NEU;
  EmitVertex();
  gl_Position = gp + NWU;
  EmitVertex();
  gl_Position = gp + NWD;
  EmitVertex();
  gl_Position = gp + SWD;
  EmitVertex();
  gl_Position = gp + SWU;
  EmitVertex();
  EndPrimitive();
}
