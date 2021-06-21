#ifndef HEADER_OBJECT_H
#define HEADER_OBJECT_H

#include "core.h"
#include "logger.h"


namespace LC { 

	struct LC_API Header {
		struct HeaderObject {
			std::string variable;
			std::size_t size;
			std::size_t location;
		};

		enum class Version { V1 = 1 };

		void write(const std::string& file);
		void read(const std::string& file);

		std::vector<HeaderObject> headerObjects;
		Version version;
	};


	
}

#endif