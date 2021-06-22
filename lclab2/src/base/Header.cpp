#include "Header.h"

namespace LC {

	void Header::write(const std::string& file) {

		std::size_t size_t_in_bytes = sizeof(std::size_t);
		std::size_t size_of_char_in_bytes = sizeof(char);

		std::ofstream ofile(file.c_str(), std::ios::out | std::ios::binary);

		if (!ofile) {
			LC_WARN("Failed to open file <{0}>", file.c_str());
			return;
		}

		// Write file version

		ofile.write((char*)&version, sizeof(Version));

		// Write number of header objects
		{
			std::size_t numObj = headerObjects.size();
			ofile.write((char*)&numObj, sizeof(std::size_t));
		}

		// Write header information

		for (const auto& hObj : headerObjects) {

			std::size_t variable_size = hObj.variable.size();

			// Write size of variable name
			ofile.write((char*)&variable_size, size_t_in_bytes);

			// Read variable
			ofile.write((char*)&hObj.variable[0], variable_size * size_of_char_in_bytes);

			// Read byte size of obj
			ofile.write((char*)&hObj.size_in_bytes, size_t_in_bytes);
			
			// Read location in body file
			ofile.write((char*)&hObj.location, size_t_in_bytes);
		}

		ofile.close();

		if (!ofile.good()) {
			LC_WARN("Write error occurred in file <{0}>", file.c_str());
		}
	}

	void Header::read(const std::string& file) {

		std::size_t size_t_in_bytes = sizeof(std::size_t);
		std::size_t size_of_char_in_bytes = sizeof(char);

		readFile = file;

		std::ifstream ifile(file.c_str(), std::ios::in | std::ios::binary);

		if (!ifile) {
			LC_WARN("Failed to open file <{0}>", file.c_str());
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

			hObj.variable.resize(variable_size);

			// Read variable
			ifile.read((char*)&hObj.variable[0], variable_size * size_of_char_in_bytes);

			// Read byte size of obj
			ifile.read((char*)&hObj.size_in_bytes, size_t_in_bytes);
			
			// Read location in body file
			ifile.read((char*)&hObj.location, size_t_in_bytes);
		}
		
		LC_CORE_INFO("Here");

		ifile.close();

		if (!ifile.good()) {
			LC_WARN("Read error occurred in file <{0}>", file.c_str());
		}
	}
}