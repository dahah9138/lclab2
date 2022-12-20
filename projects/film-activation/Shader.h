#ifndef FILM_SHADER_H
#define FILM_SHADER_H

#include <lclab2.h>

struct Shader {
	Shader() {
		const char* vertexShaderSource = "#version 330 core\n"
			"layout (location = 0) in vec3 aPos;\n"
			"uniform mat4 MVP;\n"
			"void main()\n"
			"{\n"
			"   gl_Position = MVP * vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
			"}\0";
		const char* fragmentShaderSource = "#version 330 core\n"
			"out vec4 FragColor;\n"
			"uniform vec3 color;\n"
			"void main()\n"
			"{\n"
			"   FragColor = vec4(color, 1.0f);\n"
			"}\n\0";

		// vertex shader
		int vertexShader = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
		glCompileShader(vertexShader);
		// check for shader compile errors

		// fragment shader
		int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
		glCompileShader(fragmentShader);
		// check for shader compile errors

		// link shaders
		program = glCreateProgram();
		glAttachShader(program, vertexShader);
		glAttachShader(program, fragmentShader);
		glLinkProgram(program);
		// check for linking errors

		glDeleteShader(vertexShader);
		glDeleteShader(fragmentShader);
	}
	~Shader() {
		glDeleteProgram(program);
	}
	int program;
};


#endif