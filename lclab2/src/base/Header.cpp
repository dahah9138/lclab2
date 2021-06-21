#include "Header.h"

namespace LC {

	void Header::write(const std::string& file) {

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

		for (const auto& obj : headerObjects) {
			ofile.write((char*)&obj, sizeof(HeaderObject));
		}

		ofile.close();

		if (!ofile.good()) {
			LC_WARN("Write error occurred in file <{0}>", file.c_str());
		}
	}

	void Header::read(const std::string& file) {

		std::ifstream ifile(file.c_str(), std::ios::out | std::ios::binary);

		if (!ifile) {
			LC_WARN("Failed to open file <{0}>", file.c_str());
			return;
		}

		// Read file version
		ifile.read((char*)&version, sizeof(Version));

		// Read number of header objects
		{
			std::size_t numObj;
			ifile.read((char*)&numObj, sizeof(std::size_t));
			headerObjects.resize(numObj);
		}

		// Read header information

		{
			std::size_t headerObjectsSize = sizeof(HeaderObject) * headerObjects.size();
			ifile.read((char*)&headerObjects[0], headerObjectsSize);
		}
		
		ifile.close();

		if (!ifile.good()) {
			LC_WARN("Read error occurred in file <{0}>", file.c_str());
		}
	}
}