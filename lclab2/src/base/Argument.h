#ifndef ARGUMENT_H
#define ARGUMENT_H

#include "core.h"

namespace LC {
	struct LC_API Argument {
		bool ParseArg(const std::string &arg) {
			// Validate the argument
			if (!Validate(arg)) return false;

			// Valid argument, parse
			flag = arg[1];
			variable = arg.substr(2);
		}
		bool Validate(const std::string& arg) {
			if (arg.size() < 2) return false;
			else if (arg[0] != '-') return false;
			else if (!isalpha(arg[1])) return false;
			
			return true;
		}

		char flag;
		std::string variable;
		std::string help;
		bool optional = false;
	};
}


#endif