 /*************************************************************************
  * BOV 0.1
  * A wrapper around OpenGL and GLFW (www.glfw.org) to draw simple 2D
  * graphics.
  *------------------------------------------------------------------------
  * Copyright (c) 2019-2020 Célestin Marot <marotcelestin@gmail.com>
  *
  * This software is provided 'as-is', without any express or implied
  * warranty. In no event will the authors be held liable for any damages
  * arising from the use of this software.
  *
  * Permission is granted to anyone to use this software for any purpose,
  * including commercial applications, and to alter it and redistribute it
  * freely, subject to the following restrictions:
  *
  * 1. The origin of this software must not be misrepresented; you must not
  *    claim that you wrote the original software. If you use this software
  *    in a product, an acknowledgment in the product documentation would
  *    be appreciated but is not required.
  *
  * 2. Altered source versions must be plainly marked as such, and must not
  *    be misrepresented as being the original software.
  *
  * 3. This notice may not be removed or altered from any source
  *    distribution.
  *
  *************************************************************************/

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

in vec2 posGeom;
flat in vec2 speedGeom;
flat in vec4 dataGeom;
flat in float pixelSize;

out vec4 outColor;

void main()
{
	float sdf = length(posGeom) - width;
	vec2 alpha = smoothstep(-pixelSize, pixelSize, -sdf - vec2(0.0f, outlineWidth));
	outColor = mix(outlineColor, dataGeom, alpha.y);
	outColor.a *= alpha.x;
}