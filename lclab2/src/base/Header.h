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
		enum class Option { None = 0, Read = BIT(1), Write = BIT(2) };

		void read(const std::string& file);

		// Reads from readFile
		void read();

		void write(const std::string& file);

		// Writes from writeFile
		void write();
		
		void readBody();
		void writeBody();

		void sortHeaderObjects();
		bool ValidateHeaderObjects();

		void* passData(std::size_t index);
		Header& setData(void *ptr, std::size_t index);

		// Need to store read file in case user specifies 'save'
		// instead of 'save as'
		std::string readFile;
		std::string writeFile;

		std::vector<std::pair<HeaderObject, void*>> headerObjects;
		Version version = Version::V1;
	};


	
}

#endif