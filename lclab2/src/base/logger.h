#pragma once

#include "spdlog/spdlog.h"
#include "core.h"

// Wrapper class for spdlog
namespace LC
{
	class LC_API Logger {

	public:
		static void init();

		inline static std::shared_ptr<spdlog::logger>& getCoreLogger() { return s_CoreLogger; }
		inline static std::shared_ptr<spdlog::logger>& getClientLogger() { return s_ClientLogger; }

	private:
		static std::shared_ptr<spdlog::logger> s_CoreLogger;
		static std::shared_ptr<spdlog::logger> s_ClientLogger;
	};
}
// Some simple error message macros

#define LC_CORE_INFO(...)  LC::Logger::getCoreLogger()->info(__VA_ARGS__)
#define LC_CORE_WARN(...)  LC::Logger::getCoreLogger()->warn(__VA_ARGS__)
#define LC_CORE_ERROR(...) LC::Logger::getCoreLogger()->error(__VA_ARGS__)
#define LC_CORE_CRITICAL(...) LC::Logger::getCoreLogger()->critical(__VA_ARGS__)

#define LC_INFO(...)  LC::Logger::getClientLogger()->info(__VA_ARGS__)
#define LC_WARN(...)  LC::Logger::getClientLogger()->warn(__VA_ARGS__)
#define LC_ERROR(...) LC::Logger::getClientLogger()->error(__VA_ARGS__)
#define LC_CRITICAL(...) LC::Logger::getClientLogger()->critical(__VA_ARGS__)
