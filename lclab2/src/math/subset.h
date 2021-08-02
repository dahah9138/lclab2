#ifndef SUBSET_H
#define SUBSET_H

#include <string>
#include <vector>

namespace LC { namespace Math {

	template <typename T>
	struct subset {
		enum class subset_type { empty, boundary, interior };
		subset(std::string subsetname = "emptyset") : name(subsetname) {}
		
		subset& operator = (const subset& sub) {
			copy(sub);
		}

		void copy(const subset& sub) {
			indices = sub.indices;
			type = sub.type;
		}

		T operator ()(const std::size_t& i) {
			return indices[i];
		}

		std::vector<T> indices;
		std::string name;
		subset_type type = subset_type::empty;
	};

}}


#endif