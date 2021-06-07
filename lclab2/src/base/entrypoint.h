#pragma once


#ifdef LC_PLATFORM_WIN32

namespace LC
{ 
    extern application *createApplication(int argc, char** argv);
}

int main(int argc, char** argv) {
    
    LC::Logger::init();

    auto app = LC::createApplication(argc, argv);
    int ex = app->exec();
    delete app;
    return ex;
}

#endif