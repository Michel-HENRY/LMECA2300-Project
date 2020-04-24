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

out vec4 outColor;

void main()
{
	float sdf = length(posGeom) - width;
	vec2 alpha = smoothstep(-pixelSize, pixelSize, -sdf - vec2(0.0f, outlineWidth));
	float distToCenter = 1-(length(abs(localPos-posGeom))/width);
	//outColor = vec4(distToCenter, 0, 0, 1);
	//outColor = dataGeom*distToCenter*2;
	//outColor = distToCenter*alpha.y;
	vec4 m = mix(outlineColor, dataGeom, alpha.y);
	//outColor = m;
	vec3 col = vec3(m.r, m.g, m.b);
	col = col/length(col);
	col = col/length(col);
	outColor.r = col.x;
	outColor.g = col.y;
	outColor.b = col.z;
	if(marker == 3){//Light particules
		if(posGeom.x > 0 && posGeom.y > 0){
			if(distToCenter > 0.95){
				outColor.a = 1;
			} else
				outColor.a = (distToCenter)*0.6 + 0.4;//*0.2 + 0.8;
		} else {
			outColor.a = 0;
		}
	//Faudra encore ajouter un cas pour quand on plot la pression etc
	} else if(speedGeom.x == -1000) {//Il ne l'a pas encore je ne sais pas pq
		outColor.a = 1;
		outColor.a = (distToCenter+0.5)*0.2 + 0.6;/////////Utiliser pour le fluide
		if(distToCenter > 0.5){
			outColor.a = 0.8;
		}else if(distToCenter < 0.3){
			outColor.a = 0;
		}
		//outColor.a = 0.6;
		outColor.a *= alpha.x;
	} else {
		/*
		float distToLight = 1 - length(posGeom + width/2)/width;//- length(abs(localPos + width/2)/2) + 0.5;
		if(gl_Color.a == 0){
			outColor.a = 0.5;
		} else {
			outColor.a = 1;
		}*/
		outColor.a = 1;
		outColor.a *= alpha.x;
	}
}