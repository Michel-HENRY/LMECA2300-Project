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

uniform sampler2D framebuffer;

out vec4 outColor;

float min(float a, float b){
	if(a < b){
		return a;
	} else {
		return b;
	}
}

float abs(float a){
	if(a < 0){
		return -a;
	} else {
		return a;
	}
}

void main() {
	//outColor = vec4((gl_FragCoord.xy/resolution).x, 0, 0, 1);//texture(framebuffer, gl_FragCoord.xy/resolution).rgba; // directly show frambuffer
	vec4 col = texture(framebuffer, gl_FragCoord.xy/resolution).rgba;
	vec3 fcol = texture(framebuffer, gl_FragCoord.xy/resolution).rgb;
	float intensity = fcol.r/fcol.b;// + 0.5;
	vec2 pos = gl_FragCoord.xy/resolution;
	if(col.g != 0){//We are into a pression
		float v1 = 3.5*(intensity-0.7);
		float v2 = 1.25*intensity;
		float v3 = min(0.5,intensity)*2.0;

		col.r = 1.5 - 4.0 * abs(intensity - 0.75);
		col.g = 1.5 - 4.0 * abs(intensity - 0.5 );
		col.b = 1.5 - 4.0 * abs(intensity - 0.25);
		col.a = 1;
	} else if(col.b != 0){
		col.r = 0; col.b = 0; col.g = 0;
	}
	outColor = col;
}

