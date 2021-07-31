#include "ConsoleApplication.h"

namespace LC {
	ConsoleApplication::ConsoleApplication(const Arguments &args) : _arguments{args} {
		if (!args.Validate()) args.DisplayHelp();
	}


	ConsoleApplication::~ConsoleApplication() {

	}

	void ConsoleApplication::Run() {
		// Do nothing
	}
	
	int ConsoleApplication::exec() {
		Run();
		return 1;
	}
}
