#version 150 core
in vec2 pos;
in vec2 speed;
in vec4 data;

out vec2 speedVert;
out vec4 dataVert;

void main()
{
	speedVert = speed;
	dataVert = data;
	gl_Position = vec4(pos, 0.0, 1.0);
}