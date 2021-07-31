#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include "Argument.h"
#include "logger.h"

namespace LC {
	struct LC_API Arguments {
		enum class Error {
			None = 0,
			Required = BIT(1),
			Parse = BIT(2)
		};

		Arguments(int argc, char** argv) {
			if (argc == 1) return;
			arguments.resize(argc - 1);

			// Fill arguments

			for (int i = 1; i < argc; i++) {
				arguments[i - 1].ParseArg(argv[i]);
			}
		}

		Arguments& AddRequired(char flag, const std::string& help) {
			
			options.insert({ flag, {Required, help} });
			
			// Search for flag
			Argument* arg = Search(flag);
			if (arg) {
				arg->optional = Required;
				arg->help = help;
			}
			else {
				// Error, did not pass required argument
				LC_CORE_WARN("Failed to pass required flag -{0} (Hint: {1})", flag, help.c_str());
				error = static_cast<Error>(static_cast<int>(Error::Required) | static_cast<int>(error));
			}
			return *this;
		}

		Arguments& AddOptional(char flag, const std::string& help) {
			
			options.insert({ flag, {Optional, help} });

			// Search for flag
			Argument* arg = Search(flag);
			if (arg) {
				arg->optional = Optional;
				arg->help = help;
			}
			return *this;
		}

		// Return the pointer in argument array
		// Returns 0 if not found
		Argument* Search(char flag) {
			for (auto& arg : arguments) {
				if (arg.flag == flag) return &arg;
			}
			return 0;
		}

		bool Exists(char flag) const {
			for (auto& arg : arguments) {
				if (arg.flag == flag) return true;
			}
			return false;
		}

		std::string Get(char flag) {
			Argument* arg = Search(flag);
			if (arg) return arg->variable;
			else if (options[flag].first == Required){
				LC_CORE_WARN("Failed to find flag -{0}", flag);
			}
			return std::string();
		}

		// Check that all required options have been found
		bool Validate() const {
			for (const auto& opt : options) {
				if (opt.second.first == Required) {
					bool exists = Exists(opt.first);
					if (!exists) return false;
				}
			}
			return true;
		}

		void DisplayHelp() const {
			LC_CORE_INFO("Help <-flag>: Tip");
			// Parse through all help commands
			for (const auto& opt : options) {
				LC_CORE_INFO("<-{0}>: {1}", opt.first, opt.second.second.c_str());
			}
		}

		std::vector<Argument> arguments;
		std::map<char, std::pair<bool, std::string>> options;
		const bool Optional = true;
		const bool Required = false;
		Error error = Error::None;
	};
}


#endif