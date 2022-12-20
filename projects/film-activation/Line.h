#ifndef FILM_LINE_H
#define FILM_LINE_H

#include <lclab2.h>
using namespace Magnum;

class Line {
    int *shaderProgram = NULL;
    unsigned int VBO, VAO;
    bool owner = true;
    std::vector<float> vertices;
    Vector3 startPoint;
    Vector3 endPoint;
public:
    Line() = default;
    Line(Vector3 start, Vector3 end, int &shaderProgram_) {
        Create(start, end, shaderProgram_);
    }

    void Create(Vector3 start, Vector3 end, int& shaderProgram_) {
        startPoint = start;
        endPoint = end;
        shaderProgram = &shaderProgram_;

        vertices = {
             start.x(), start.y(), start.z(),
             end.x(), end.y(), end.z(),

        };

        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    void Move(Line&& line) {
        shaderProgram = line.shaderProgram;
        VAO = line.VAO;
        VBO = line.VBO;
        vertices = line.vertices;
        startPoint = line.startPoint;
        endPoint = line.endPoint;

        line.owner = false;
    }

    int draw(Matrix4 &mvp, Vector3 &color) {
        if (!shaderProgram) return 0;
        glUseProgram(*shaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(*shaderProgram, "MVP"), 1, GL_FALSE, &mvp[0][0]);
        glUniform3fv(glGetUniformLocation(*shaderProgram, "color"), 1, &color[0]);

        glBindVertexArray(VAO);
        glDrawArrays(GL_LINES, 0, 2);
        return 1;
    }

    ~Line() {
        if (!owner) return;
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
    }
};

#endif