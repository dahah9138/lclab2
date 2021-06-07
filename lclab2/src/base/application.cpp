#include "application.h"

namespace LC
{
	application::application(const Arguments& arguments): Platform::Application{arguments} {
	

	}
	application::application(const Arguments& arguments, const Configuration& configuration) : Platform::Application{ arguments, configuration } {

	}

	application::~application() {
		Debug{} << "Terminating application.\n";
	}
}