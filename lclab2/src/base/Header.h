#ifndef HEADER_OBJECT_H
#define HEADER_OBJECT_H

#include "core.h"
#include "logger.h"


namespace LC {

	struct LC_API Header {
		struct HeaderObject {
			HeaderObject() = default;
			HeaderObject(const std::string& var, const std::size_t& sz, const std::size_t& loc)
				: variable(var), size_in_bytes(sz), location(loc), copy(true) {
			}
			HeaderObject(const std::string& var, const std::size_t& sz)
				: variable(var), size_in_bytes(sz), location(0), copy(true) {
			}

			std::string variable{};
			std::size_t size_in_bytes = 0;
			std::size_t location = 0;
			bool copy = true;
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

		template <typename T>
		Header& addObject(const std::string& name, T& value) {
			return addObject(HeaderObject{ name, sizeof(T) }, static_cast<void*>(&value));
		}

		template <typename T>
		Header& addObject(const std::string& name, T* ptr, std::size_t count) {
			return addObject(HeaderObject{ name, sizeof(T) * count }, static_cast<void*>(ptr));
		}

		template <typename T>
		bool getObject(const std::string& name, T& outValue, const T& defaultValue) {
			void* raw = passData(name);
			if (raw) {
				std::memcpy(&outValue, raw, sizeof(T));
				delete[] static_cast<char*>(raw);
				return true;
			}
			outValue = defaultValue;
			LC_CORE_INFO("{0} initialized", name.c_str());
			return false;
		}

		template <typename T, std::size_t N>
		bool getArray(const std::string& name, std::array<T, N>& outArray, const T& defaultValue) {
			void* raw = passData(name);
			if (raw) {
				std::memcpy(outArray.data(), raw, sizeof(T) * N);
				delete[] static_cast<char*>(raw);
				return true;
			}
			outArray.fill(defaultValue);
			LC_CORE_INFO("{0} initialized", name.c_str());
			return false;
		}

		template <typename T>
		std::unique_ptr<T[]> getBuffer(const std::string& name, std::size_t count, bool zeroInitIfMissing = true) {
			void* raw = passData(name);
			if (raw) {
				auto buf = std::unique_ptr<T[]>(new T[count]);
				std::memcpy(buf.get(), raw, sizeof(T) * count);
				delete[] static_cast<char*>(raw);
				return buf;
			}
			LC_CORE_WARN("{0} was not found", name.c_str());
			if (zeroInitIfMissing) {
				LC_CORE_INFO("{0} initialized", name.c_str());
				return std::unique_ptr<T[]>(new T[count]());
			}
			return nullptr;
		}

		std::string readFile;
		std::string writeFile;

		std::vector<std::pair<HeaderObject, void*>> headerObjects;
		Version version = Version::V1;

	};
	typedef std::pair<Header::HeaderObject, void*> HeaderPair;
}

#endif