// For visualizing T6+ FQ1 and T6- FQ2
#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

#declare cam=CAM;
camera {
	orthographic
	location <0,0,100>
	right -CAM*x*image_width/image_height
	up CAM*z
	sky <0,1,0>
	look_at <0,0,0>
}
background{rgb 1}

light_source{<-3000,2500,7000> color rgb <0.81,0.79,0.79>}
light_source{<5000,-1500,4500> color rgb <0.41,0.43,0.43>}

#declare f0=finish{phong 0.6 specular 0.2 ambient 0.15}
#declare t0=texture{pigment{rgb <1.,0.75,0.80>} finish{f0}} // PINK, lethal
#declare t1=texture{pigment{rgb <0, 0.2, 0.8>} finish{f0}} // BLUE, nonlethal
#declare t2=texture{pigment{rgb <0.95,0.22,0.18>} finish{f0}} // RED, nonlethal
#declare t3=texture{pigment{rgb <0, 0.2, 0.8>} finish{f0}} // BLUE, nonlethal
#declare t4=texture{pigment{rgb <0.75,0.5,0.9>} finish{f0}}
#declare t5=texture{pigment{rgb <0.63,0.75,0.72>} finish{f0}}
#declare t6=texture{pigment{rgb <0.75,0.5,0.9>} finish{f0}}
#declare t7=texture{pigment{rgb <0.63,0.75,0.72>} finish{f0}}
#declare r=1.0;

#include "bac.pov"
