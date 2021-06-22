#ifndef HEADER_OBJECT_H
#define HEADER_OBJECT_H

#include "core.h"
#include "logger.h"


namespace LC { 

	struct LC_API Header {
		struct HeaderObject {
			std::string variable;
			std::size_t size_in_bytes;
			std::size_t location;
		};

		enum class Version { V1 = 1 };

		void write(const std::string& file);
		void read(const std::string& file);

		// Need to store read file in case user specifies 'save'
		// instead of 'save as'
		std::string readFile;

		// HeaderObject and pointer to data pair
		// Initialized to zero from calling write
		std::vector<HeaderObject> headerObjects;
		Version version = Version::V1;
	};


	
}

#endif