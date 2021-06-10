#pragma once


#ifdef LC_PLATFORM_WIN32 || LC_PLATFORM_UNIX || LC_PLATFORM_MACOS

namespace LC { 
    extern Application *createApplication(int argc, char** argv);
}

int main(int argc, char** argv) {
    
    LC::Logger::init();

    auto app = LC::createApplication(argc, argv);
    int ex = app->exec();
    delete app;
    return ex;
}

#endif
