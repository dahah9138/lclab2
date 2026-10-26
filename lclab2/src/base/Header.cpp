#include "Header.h"

namespace LC {

	Header::~Header() {
		clean();
	}

	void Header::write(const std::string& file) {
		writeFile = file;

		std::ofstream ofile(file, std::ios::out | std::ios::binary);
		if (!ofile) {
			LC_CORE_WARN("Failed to open file <{0}>", file.c_str());
			return;
		}

		sortHeaderObjects();

		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", file.c_str());
			return;
		}

		ofile.write(reinterpret_cast<const char*>(&version), sizeof(Version));

		std::size_t numObj = headerObjects.size();
		ofile.write(reinterpret_cast<const char*>(&numObj), sizeof(std::size_t));

		for (const auto& hObj : headerObjects) {
			std::size_t variable_size = hObj.first.variable.size();

			ofile.write(reinterpret_cast<const char*>(&variable_size), sizeof(std::size_t));
			if (variable_size > 0) {
				ofile.write(hObj.first.variable.data(), variable_size * sizeof(char));
			}
			ofile.write(reinterpret_cast<const char*>(&hObj.first.size_in_bytes), sizeof(std::size_t));
			ofile.write(reinterpret_cast<const char*>(&hObj.first.location), sizeof(std::size_t));
		}

		ofile.close();

		if (!ofile.good()) {
			LC_CORE_WARN("Write error occurred in file <{0}>", file.c_str());
		}
	}

	void Header::read(const std::string& file) {
		readFile = file;

		std::ifstream ifile(file, std::ios::in | std::ios::binary);
		if (!ifile) {
			LC_CORE_WARN("Failed to open file <{0}>", file.c_str());
			return;
		}

		clean();

		ifile.read(reinterpret_cast<char*>(&version), sizeof(Version));

		std::size_t numObj = 0;
		ifile.read(reinterpret_cast<char*>(&numObj), sizeof(std::size_t));
		if (!ifile.good()) {
			LC_CORE_WARN("Failed to read header count in file <{0}>", file.c_str());
			return;
		}

		headerObjects.resize(numObj);

		for (auto& hObj : headerObjects) {
			std::size_t variable_size = 0;
			ifile.read(reinterpret_cast<char*>(&variable_size), sizeof(std::size_t));
			if (!ifile.good()) break;

			hObj.first.variable.resize(variable_size);
			if (variable_size > 0) {
				ifile.read(&hObj.first.variable[0], variable_size * sizeof(char));
			}

			ifile.read(reinterpret_cast<char*>(&hObj.first.size_in_bytes), sizeof(std::size_t));
			ifile.read(reinterpret_cast<char*>(&hObj.first.location), sizeof(std::size_t));
			hObj.second = nullptr;
		}

		ifile.close();

		if (!ifile.good()) {
			LC_CORE_WARN("Read error occurred in file <{0}>", file.c_str());
			return;
		}

		sortHeaderObjects();

		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", readFile.c_str());
			return;
		}
	}

	void Header::readBody() {
		std::ifstream ifile(readFile, std::ios::in | std::ios::binary);
		if (!ifile) {
			LC_CORE_WARN("Failed to open file <{0}>", readFile.c_str());
			return;
		}

		ifile.seekg(sizeof(Version), std::ios::beg);

		std::size_t numObj = 0;
		ifile.read(reinterpret_cast<char*>(&numObj), sizeof(std::size_t));

		if (!ifile.good() || numObj != headerObjects.size()) {
			LC_CORE_WARN("Abort: Header object count mismatch in file <{0}>", readFile.c_str());
			return;
		}

		for (std::size_t i = 0; i < numObj; ++i) {
			std::size_t variable_size = 0;
			ifile.read(reinterpret_cast<char*>(&variable_size), sizeof(std::size_t));
			ifile.seekg(variable_size * sizeof(char) + 2 * sizeof(std::size_t), std::ios::cur);
		}

		sortHeaderObjects();

		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", readFile.c_str());
			return;
		}

		for (std::size_t i = 0; i < numObj; ++i) {
			if (headerObjects[i].second && !headerObjects[i].first.copy) {
				delete[] static_cast<char*>(headerObjects[i].second);
			}

			headerObjects[i].second = new char[headerObjects[i].first.size_in_bytes];
			headerObjects[i].first.copy = false;

			ifile.read(static_cast<char*>(headerObjects[i].second), headerObjects[i].first.size_in_bytes);
		}

		ifile.close();

		if (!ifile.good()) {
			LC_CORE_WARN("Read error occurred in file <{0}>", readFile.c_str());
		}
	}

	void Header::writeBody() {
		std::ofstream ofile(writeFile, std::ios::out | std::ios::binary | std::ios::app);
		if (!ofile) {
			LC_CORE_WARN("Failed to open file <{0}>", writeFile.c_str());
			return;
		}

		sortHeaderObjects();

		if (!ValidateHeaderObjects()) {
			LC_CORE_WARN("Abort: Invalid locations specified in file <{0}>", writeFile.c_str());
			return;
		}

		for (std::size_t i = 0; i < headerObjects.size(); ++i) {
			if (headerObjects[i].second) {
				ofile.write(static_cast<const char*>(headerObjects[i].second), headerObjects[i].first.size_in_bytes);
			}
			else {
				LC_CORE_WARN("Abort: Invalid header object <{1}> when writing body to file <{0}>",
					writeFile.c_str(), headerObjects[i].first.variable.c_str());
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
			[](const std::pair<HeaderObject, void*>& a, const std::pair<HeaderObject, void*>& b) {
				return a.first.location < b.first.location;
			});
	}

	bool Header::ValidateHeaderObjects() {
		for (std::size_t i = 1; i < headerObjects.size(); ++i) {
			if (headerObjects[i].first.location <= headerObjects[i - 1].first.location) {
				LC_CORE_WARN("Error: Locations {1} and {0}", headerObjects[i].first.location, headerObjects[i - 1].first.location);
				return false;
			}
		}
		return true;
	}

	void* Header::passData(std::size_t& index) {
		if (index < headerObjects.size()) {
			headerObjects[index].first.copy = true;
			void* tmp = headerObjects[index].second;
			headerObjects[index++].second = nullptr;
			return tmp;
		}
		return nullptr;
	}

	void* Header::passData(const std::string& id) {
		for (auto hObj = headerObjects.begin(); hObj != headerObjects.end(); ++hObj) {
			if (hObj->first.variable == id) {
				hObj->first.copy = true;
				void* tmp = hObj->second;
				hObj->second = nullptr;
				headerObjects.erase(hObj);
				return tmp;
			}
		}

		LC_CORE_WARN("ID [{0}] not found", id.c_str());
		return nullptr;
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
		for (auto& obj : headerObjects) {
			if (obj.second && !obj.first.copy) {
				LC_CORE_INFO("Deleting unused loaded variable [{0}]", obj.first.variable.c_str());
				delete[] static_cast<char*>(obj.second);
			}
			obj.second = nullptr;
		}
		headerObjects.clear();
	}

	void Header::read() {
		read(readFile);
	}

	void Header::write() {
		write(writeFile);
	}

}