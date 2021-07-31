#ifndef CONSOLE_APPLICATION_H
#define CONSOLE_APPLICATION_H

#include "core.h"
#include "logger.h"
#include "solver/Solver.h"
#include "Header.h"
#include "Arguments.h"

namespace LC {
	class LC_API ConsoleApplication {
    public:
		
        explicit ConsoleApplication(const Arguments &args);
		virtual ~ConsoleApplication();


		virtual void Run();
		int exec();

		// Solver
    	std::unique_ptr<Solver> _solver;
		// Header
		Header _header;
		Arguments _arguments;

	};
	
}


#endif
