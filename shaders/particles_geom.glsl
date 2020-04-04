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
	vec4 fillColor;
	vec4 outlineColor;
	vec2 localPos;
	vec2 localScale;
	float width;
	float marker;
	float outlineWidth;
	// float rotation;
	int space_type; // 0: normal sizes, 1: unzoomable, 2: unmodifable pixel size
};

layout (std140) uniform worldBlock
{
	vec2 resolution;
	vec2 translate;
	float zoom;
};

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

in vec2 speedVert[1];
in vec4 dataVert[1];

out vec2 posGeom;
flat out vec2 speedGeom;
flat out vec4 dataGeom;
flat out float pixelSize;

void main()
{
	float minRes = min(resolution.x, resolution.y);
	vec2 resRatio = minRes / resolution;

	// for particles, space_type has no effect
	vec2 scaling = resRatio * zoom;
	vec2 translation = resRatio * zoom * (localPos + translate);
	pixelSize = 2.0 / (minRes * zoom);

	float w = width + pixelSize;

	vec2 p = gl_in[0].gl_Position.xy * localScale;

	vec2 center = p * scaling + translation;
	vec2 upRight = center + w * scaling;
	vec2 downLeft = center - w * scaling;

	// a bit of culling
	if(any(lessThanEqual(upRight, -vec2(1.0))) ||
	   any(greaterThanEqual(downLeft, vec2(1.0))))
		return;

	vec2 upLeft = vec2(downLeft.x, upRight.y);
	vec2 downRight = vec2(upRight.x, downLeft.y);

	speedGeom = speedVert[0];
	dataGeom = dataVert[0];

	posGeom = vec2(-w, w);
	gl_Position = vec4(upLeft, 0.0, 1.0);
	EmitVertex();
	posGeom = vec2(-w);
	gl_Position = vec4(downLeft, 0.0, 1.0);
	EmitVertex();
	posGeom = vec2(w);
	gl_Position = vec4(upRight, 0.0, 1.0);
	EmitVertex();
	posGeom = vec2(w, -w);
	gl_Position = vec4(downRight, 0.0, 1.0);
	EmitVertex();
	EndPrimitive();
}