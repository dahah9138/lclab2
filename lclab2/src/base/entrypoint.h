#pragma once


#if defined(LC_PLATFORM_WIN32) || defined(LC_PLATFORM_UNIX) || defined(LC_PLATFORM_MACOS)

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
