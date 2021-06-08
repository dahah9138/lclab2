#ifndef TORUS_H
#define TORUS_H

#include "Sheet.h"


namespace LC {

struct LC_API Torus: public Sheet {


    void Init() override;
    void Draw(const Magnum::Containers::Optional<Magnum::ArcBall>& arcball, const Magnum::Matrix4& projection) override;

};

}


#endif
