#ifndef HEADER_OBJECT_H
#define HEADER_OBJECT_H

#include "core.h"
#include "logger.h"


namespace LC { 

	struct LC_API Header {
		struct HeaderObject {
			HeaderObject() = default;
			HeaderObject(const std::string& var, const std::size_t& sz, const std::size_t& loc) : variable(var), size_in_bytes(sz), location(loc), copy(true) {}
			HeaderObject(const std::string& var, const std::size_t& sz) : variable(var), size_in_bytes(sz), location(0), copy(true) {}
			std::string variable;
			std::size_t size_in_bytes;
			std::size_t location;
			// true -> External data was passed to the HeaderObject
			bool copy;
		};

		enum class Version { V1 = 1 };
		enum class Option { None = 0, Read = BIT(1), Write = BIT(2) };

		~Header();

		void read(const std::string& file);

		// Reads from readFile
		void read();

		void write(const std::string& file);

		// Writes from writeFile
		void write();
		
		void readBody();
		void writeBody();
		void clean();

		void sortHeaderObjects();
		bool ValidateHeaderObjects();

		// Iterates through the data using the starting index passed.
		// The index is incremented by one each time.
		void* passData(std::size_t &index);
		void* passData(const std::string& id);
		// Set location to position in headerObjects vector
		Header& addObject(const std::pair<HeaderObject, void*> &obj);
		// Set location to position in headerObjects vector
		Header& addObject(HeaderObject hobj, void* ptr);
		Header& addObject(HeaderObject hobj, void *ptr, std::size_t& id);
		Header& setData(void *ptr, std::size_t index);

		// Feed object, data pairs to header
		Header& operator << (const std::pair<HeaderObject, void*>& obj);

		// Need to store read file in case user specifies 'save'
		// instead of 'save as'
		std::string readFile;
		std::string writeFile;

		std::vector<std::pair<HeaderObject, void*>> headerObjects;
		Version version = Version::V1;
	};

	typedef std::pair<Header::HeaderObject, void*> HeaderPair;


	
}

#endif