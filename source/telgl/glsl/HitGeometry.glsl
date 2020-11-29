#version 430 core

layout(points) in;
in vec3 vColor[];

layout(line_strip, max_vertices = 2) out;
out vec3 fColor;

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

void main(){

  vec3 pos_hit = gl_in[0].gl_Position.xyz;
  mat4 pvmMatrix = transformData.proj * transformData.view * transformData.model;

  vec4 gp = pvmMatrix * vec4(pos_hit, 1.0);

  vec4 offset = pvmMatrix * vec4( 0.0,  0.0,  8.0, 0.0);
  //fColor = vColor[0].rgb;
  fColor = vec3(1.0, 0.0, 0.0);
  // fColor = vec3(vColor[0].b , vColor[0].r, vColor[0].g); //shift color
  gl_Position = gp + offset;
  EmitVertex();
  gl_Position = gp - offset;
  EmitVertex();
  EndPrimitive();

  // fColor = vec3(1.0, 1.0, 1.0);

  // // fColor = vec3(vColor[0].b , vColor[0].r, vColor[0].g); //shift color
  // vec4 halfThick(0.2, 0.2, 1.0, 0.0);
  // //                                     N/S       U/D       E/W
  // //                                     X         Y         Z
  // vec4 NEU = pvmMatrix * vec4(  halfThick.x,  halfThick.y,  halfThick.z, 0.0);
  // vec4 NED = pvmMatrix * vec4(  halfThick.x, -halfThick.y,  halfThick.z, 0.0);
  // vec4 NWU = pvmMatrix * vec4(  halfThick.x,  halfThick.y, -halfThick.z, 0.0);
  // vec4 NWD = pvmMatrix * vec4(  halfThick.x, -halfThick.y, -halfThick.z, 0.0);
  // vec4 SEU = pvmMatrix * vec4( -halfThick.x,  halfThick.y,  halfThick.z, 0.0);
  // vec4 SED = pvmMatrix * vec4( -halfThick.x, -halfThick.y,  halfThick.z, 0.0);
  // vec4 SWU = pvmMatrix * vec4( -halfThick.x,  halfThick.y, -halfThick.z, 0.0);
  // vec4 SWD = pvmMatrix * vec4( -halfThick.x, -halfThick.y, -halfThick.z, 0.0);
  // gl_Position = gp + NED;
  // EmitVertex();
  // gl_Position = gp + NWD;
  // EmitVertex();
  // gl_Position = gp + SWD;
  // EmitVertex();
  // gl_Position = gp + SED;
  // EmitVertex();
  // gl_Position = gp + SEU;
  // EmitVertex();
  // gl_Position = gp + SWU;
  // EmitVertex();
  // gl_Position = gp + NWU;
  // EmitVertex();
  // gl_Position = gp + NEU;
  // EmitVertex();
  // gl_Position = gp + NED;
  // EmitVertex();
  // gl_Position = gp + SED;
  // EmitVertex();
  // gl_Position = gp + SEU;
  // EmitVertex();
  // gl_Position = gp + NEU;
  // EmitVertex();
  // gl_Position = gp + NWU;
  // EmitVertex();
  // gl_Position = gp + NWD;
  // EmitVertex();
  // gl_Position = gp + SWD;
  // EmitVertex();
  // gl_Position = gp + SWU;
  // EmitVertex();
  // EndPrimitive();

}
