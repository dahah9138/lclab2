#include "Header.h"

namespace LC {

	Header::~Header() {
		for (auto& elem : headerObjects) {
			if (elem.second && !elem.first.copy) {
				delete[] elem.second;
			}
			elem.second = 0;
		}
	}

	void Header::write(const std::string& file) {

		writeFile = file;

		std::size_t size_t_in_bytes = sizeof(std::size_t);
		std::size_t size_of_char_in_bytes = sizeof(char);

		std::ofstream ofile(file.c_str(), std::ios::out | std::ios::binary);

		if (!ofile) {
			LC_CORE_WARN("Failed to open file <{0}>", file.c_str());
			return;
		}

		// Write file version

		ofile.write((char*)&version, sizeof(Version));

		// Write number of header objects
		{
			std::size_t numObj = headerObjects.size();
			ofile.write((char*)&numObj, sizeof(std::size_t));
		}

		// Order header objects based on location specified
		sortHeaderObjects();

		// Validate header objects
		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", readFile.c_str());
		}

		// Write header information

		for (const auto& hObj : headerObjects) {

			std::size_t variable_size = hObj.first.variable.size();

			// Write size of variable name
			ofile.write((char*)&variable_size, size_t_in_bytes);

			// Read variable
			ofile.write((char*)&hObj.first.variable[0], variable_size * size_of_char_in_bytes);

			// Read byte size of obj
			ofile.write((char*)&hObj.first.size_in_bytes, size_t_in_bytes);
			
			// Read location in body file
			ofile.write((char*)&hObj.first.location, size_t_in_bytes);
		}

		ofile.close();

		if (!ofile.good()) {
			LC_CORE_WARN("Write error occurred in file <{0}>", file.c_str());
		}
	}

	void Header::read(const std::string& file) {

		std::size_t size_t_in_bytes = sizeof(std::size_t);
		std::size_t size_of_char_in_bytes = sizeof(char);

		readFile = file;

		std::ifstream ifile(file.c_str(), std::ios::in | std::ios::binary);

		if (!ifile) {
			LC_CORE_WARN("Failed to open file <{0}>", file.c_str());
			return;
		}

		// Read file version
		ifile.read((char*)&version, sizeof(Version));

		// Read number of header objects
		{
			std::size_t numObj;
			ifile.read((char*)&numObj, size_t_in_bytes);

			headerObjects.resize(numObj);
		}

		// Read header information

		for (auto & hObj : headerObjects) {

			// Read size of variable name
			std::size_t variable_size;
			ifile.read((char*)&variable_size, size_t_in_bytes);

			hObj.first.variable.resize(variable_size);

			// Read variable
			ifile.read((char*)&hObj.first.variable[0], variable_size * size_of_char_in_bytes);

			// Read byte size of obj
			ifile.read((char*)&hObj.first.size_in_bytes, size_t_in_bytes);
			
			// Read location in body file
			ifile.read((char*)&hObj.first.location, size_t_in_bytes);
		}
	

		ifile.close();

		if (!ifile.good()) {
			LC_CORE_WARN("Read error occurred in file <{0}>", file.c_str());
		}

		// Order header objects based on location specified
		sortHeaderObjects();

		// Validate header objects
		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", readFile.c_str());
		}
	}

	void Header::readBody() {

		std::size_t size_t_in_bytes = sizeof(std::size_t);
		std::size_t size_of_char_in_bytes = sizeof(char);

		std::ifstream ifile(readFile.c_str(), std::ios::in | std::ios::binary);

		if (!ifile) {
			LC_CORE_WARN("Failed to open file <{0}>", readFile.c_str());
			return;
		}

		// Read file version
		{
			enum class Version v;
			ifile.read((char*)&v, sizeof(Version));
		}

		// Read number of header objects
		std::size_t numObj;
		{
			ifile.read((char*)&numObj, size_t_in_bytes);
		}

		// Read header information

		{
			char buffer[256];

			for (auto& hObj : headerObjects) {

				// Read size of variable name
				std::size_t variable_size;
				ifile.read((char*)&variable_size, size_t_in_bytes);

				// Read variable
				ifile.read(buffer, variable_size * size_of_char_in_bytes);

				// Read byte size of obj
				ifile.read(buffer, size_t_in_bytes);

				// Read location in body file
				ifile.read(buffer, size_t_in_bytes);
			}
		}

		// Order header objects based on location specified
		sortHeaderObjects();

		// Validate header objects
		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", readFile.c_str());
		}

		// Read body information
		for (int i = 0; i < numObj; i++) {

			// Allocate data
			headerObjects[i].second = new char[headerObjects[i].first.size_in_bytes];
			headerObjects[i].first.copy = false;

			// Read data
			ifile.read((char*)headerObjects[i].second, headerObjects[i].first.size_in_bytes);
		}

		ifile.close();

		if (!ifile.good()) {
			LC_CORE_WARN("Read error occurred in file <{0}>", readFile.c_str());
		}

	}

	void Header::writeBody() {

		std::size_t size_t_in_bytes = sizeof(std::size_t);
		std::size_t size_of_char_in_bytes = sizeof(char);

		std::ofstream ofile(writeFile.c_str(), std::ios::out | std::ios::binary | std::ios::app);

		if (!ofile) {
			LC_CORE_WARN("Failed to open file <{0}>", writeFile.c_str());
			return;
		}



		// Order header objects based on location specified
		sortHeaderObjects();

		// Validate header objects
		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", writeFile.c_str());
		}

		// Write body information
		for (int i = 0; i < headerObjects.size(); i++) {

			// Write data
			if (headerObjects[i].second) {
				ofile.write((char*)headerObjects[i].second, headerObjects[i].first.size_in_bytes);
			}
			else {
				LC_CORE_WARN("Abort: Invalid header object <{1}> when writing body to file <{0}>", writeFile.c_str(), 
					headerObjects[i].first.variable.c_str());
				return;
			}
		}

		ofile.close();

		if (!ofile.good()) {
			LC_CORE_WARN("Write error occurred in file <{0}>", writeFile.c_str());
		}

	}

	void Header::sortHeaderObjects() {

		std::sort(headerObjects.begin(), headerObjects.end(),
			[](const std::pair<HeaderObject, void*>& a, const std::pair<HeaderObject, void*>& b){
				return a.first.location < b.first.location; 
			});

	}

	bool Header::ValidateHeaderObjects() {

		for (int i = 1; i < headerObjects.size(); i++) {
			if (headerObjects[i].first.location <= headerObjects[i - 1].first.location) {
				LC_CORE_WARN("Error: Locations {1} and {0}", headerObjects[i].first.location, headerObjects[i - 1].first.location);
				return false;
			}
		}

		return true;
	}

	void* Header::passData(std::size_t &index) {
		// return data at that index and set to null
		if (headerObjects.size() > index) {
			headerObjects[index].first.copy = true;
			void *tmp = headerObjects[index].second;
			headerObjects[index++].second = 0;
			return tmp;
		}
		else return 0;
	}

	Header& Header::setData(void* ptr, std::size_t index) {
		if (index < headerObjects.size()) {
			headerObjects[index].second = ptr;
		}
		return *this;
	}

	Header& Header::addObject(const HeaderPair& obj) {
		std::pair<HeaderObject, void*> o = obj;
		o.first.location = headerObjects.size();
		headerObjects.emplace_back(o);
		return *this;
	}

	Header& Header::addObject(HeaderObject hobj, void* ptr, std::size_t& loc) {

		hobj.location = loc++;

		headerObjects.emplace_back(HeaderPair(hobj, ptr));
		return *this;
	}

	Header& Header::addObject(HeaderObject hobj, void* ptr) {

		hobj.location = headerObjects.size();

		headerObjects.emplace_back(HeaderPair(hobj, ptr));
		return *this;
	}

	Header& Header::operator << (const HeaderPair& obj) {
		return addObject(obj);
	}

	void Header::clean() {
		Header tmp{};
		// Need to go through list and remove any entries that
		// have non-null ptrs and were not passed to an external ptr
		for (auto& obj : headerObjects) {
			if (obj.second && !obj.first.copy) {
				LC_CORE_INFO("Deleting unused loaded variable [{0}]", obj.first.variable.c_str());
				delete[] obj.second;
			}
			obj.second = 0;
		}
		headerObjects.swap(tmp.headerObjects);
	}

	void Header::read() {
		read(readFile);
	}

	void Header::write() {
		write(writeFile);
	}


}