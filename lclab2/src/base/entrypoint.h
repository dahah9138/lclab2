#pragma once


#if defined(LC_PLATFORM_WIN32) || defined(LC_PLATFORM_UNIX) || defined(LC_PLATFORM_MACOS)

namespace LC {
#ifndef LC_CONSOLE_APP
    extern Application *createApplication(int argc, char** argv);
#else
	extern ConsoleApplication *createConsoleApplication(int argc, char** argv);
#endif
}

int main(int argc, char** argv) {
    
    LC::Logger::init();
#ifndef LC_CONSOLE_APP
    auto app = LC::createApplication(argc, argv);
#else
	auto app = LC::createConsoleApplication(argc, argv);
#endif
    int ex = app->exec();
    delete app;
    return ex;
}

#endif
