static GLchar fullscreen_vert[]={"#version 150 core\n"
"const vec2 coords[3] = vec2[]\n"
"(\n"
"  vec2(-1.0, -1.0),\n"
"  vec2( 3.0, -1.0),\n"
"  vec2(-1.0,  3.0)\n"
");\n"
"void main() { gl_Position = vec4(coords[gl_VertexID], 0.0, 1.0); }\n"
};
