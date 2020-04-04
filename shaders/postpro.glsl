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

void main()
{
	outColor = texture(framebuffer, gl_FragCoord.xy/resolution).rgba; // directly show frambuffer
}

