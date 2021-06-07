#include "logger.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace LC
{
	std::shared_ptr<spdlog::logger> Logger::s_CoreLogger;
	std::shared_ptr<spdlog::logger> Logger::s_ClientLogger;

	void Logger::init() {
		spdlog::set_pattern("%^[%T] %n: %v%$"); // timestamp logger, (client or engine), message
		s_CoreLogger = spdlog::stdout_color_mt("lclab2");
		s_CoreLogger->set_level(spdlog::level::trace);

		s_ClientLogger = spdlog::stdout_color_mt("client");
		s_ClientLogger->set_level(spdlog::level::trace);
	}
}