#version 430 core
in vec4 pos;
// pos.xyz position
// pos.w  mode: 0 global-mm, 1 layerlocal-mm-center, 2 layerlocal-mm-corner, 3 layerlocal-pixel-corner

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
  vec4 pos_global;
  // pos.w  mode: 0 global-mm, 1 layerlocal-mm-center, 2 layerlocal-mm-corner, 3 layerlocal-pixel-corner
  if(abs(pos.w) < 0.1){ //global-mm
    pos_global = vec4(pos.xyz, 1.0);
    fColor = vec3(0.0, 0.0, 0.0);
  }
  else{
    GeoDataLayer geo;
    bool found = false;
    for (int k = 0; k < 20; ++k){
      if( abs(pos.z - geoData.ly[k].id.x) < 0.1){
        geo = geoData.ly[k];
        found = true;
        break;
      }
    }
    if(!found){
      return ;
    }
    if(abs(pos.w-1.0) < 0.1){ // local-mm-center-zero
      pos_global = geo.trans * vec4(geo.pos.xy + pos.xy, geo.pos.z, 1.0);
    }
    else if(abs(pos.w-2.0) < 0.1){ // local-mm-corner-zero
      pos_global = geo.trans * vec4(geo.pos.xy - geo.size.xy/2 + pos.xy, geo.pos.z, 1.0);
    }
    else if(abs(pos.w-3.0) < 0.1){ // local-pixel-corner-zero
      pos_global = geo.trans * vec4(geo.pos.xy - geo.pitch.xy*geo.npixel.xy/2.0 + geo.pitch.xy*pos.xy, geo.pos.z, 1.0);
    }
    else{  // unknown, do nothing
      return;
    }
    fColor = geo.color.rgb;
  }

  mat4 pvmMatrix = transformData.proj * transformData.view * transformData.model;
  gl_Position = pvmMatrix * vec4(pos_global.xyz, 1.0);
}
