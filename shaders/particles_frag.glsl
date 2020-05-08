#version 150 core

layout (std140) uniform objectBlock
{
	vec4 fillColor;    // set by bov_points_set_color()
	vec4 outlineColor; // set by bov_points_set_outline_color()
	vec2 localPos;     // set by bov_points_set_pos()
	vec2 localScale;   // set by bov_points_scale()
	float width;       // set by bov_points_set_width()
	float marker;      // set by bov_poitns_set_markers()
	float outlineWidth;// set by bov_points_set_outline_width()
	int space_type;    // set by bov_points_set_space_type()
                       // 0: normal sizes, 1: unzoomable, 2: unmodifable pixel size
};

layout (std140) uniform worldBlock
{
	vec2 resolution;
	vec2 translate;
	float zoom;
};

in vec2 posGeom;
flat in vec2 speedGeom;
flat in vec4 dataGeom;
flat in float pixelSize; // = 2.0 / (min(resolution.x, resolution.y) * zoom)

uniform sampler2D textu;

out vec4 outColor;

void main() {
	float sdf = length(posGeom) - width;
	vec2 alpha = smoothstep(-pixelSize, pixelSize, -sdf - vec2(0.0f, outlineWidth));
	float distToCenter = 1-(length(abs(localPos-posGeom))/width);
	vec4 m = mix(outlineColor, dataGeom, alpha.y);
	vec3 col = vec3(m.r, m.g, m.b);
	col = col/length(col);
	col = col/length(col);
	outColor.r = col.x;
	outColor.g = col.y;
	outColor.b = col.z;
	if(marker == 3){//Light particules
		if(posGeom.x < 0 && posGeom.y < 0){
			if(distToCenter > 0.99){
				outColor.a = 1;
			} else
				outColor.a = (distToCenter)*0.8 + 0.2;//*0.2 + 0.8;
		} else {
			outColor.a = 0;
		}
		outColor.g = 0;
	} else if(speedGeom.x == -1000) {//Fluid
		outColor.a = 1;
		//outColor.a = (distToCenter+0.5)*0.2 + 0.6;
		if(distToCenter < 0.3)
			outColor.a = 0;
		outColor.g = 0;
		/*
		if(distToCenter > 0.5){
			outColor.a = 0.8;
		}else if(distToCenter < 0.5){
			outColor.a = 0;
		}*/
		outColor.a *= alpha.x;
	} else if(speedGeom.x == -2000){//Shadow particles
		outColor = m;
		outColor.a = 1;
		outColor.a *= alpha.x;		
	} else if(speedGeom.x >= 4000){//Continious field
		outColor.r = distToCenter*(speedGeom.x - 4000);
		outColor.g = 1;//1 - speedGeom.x + 4000;///////////Mettre le bleu evidemment
		outColor.b = distToCenter;//1 - speedGeom.x + 4000;
		outColor.a = 0.2;//((distToCenter*distToCenter));/////////////ca ne va pas
		outColor.a *= alpha.x;
	} else {
		outColor.g = 0;
		outColor.a = 1;
		outColor.a *= alpha.x;
	}
}