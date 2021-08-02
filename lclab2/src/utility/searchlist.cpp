#include "searchlist.h"

namespace LC { namespace Utility {
	template <> uisearchlist create_search_list(std::vector<uirange_pair>&);
	template <> isearchlist create_search_list(std::vector<irange_pair>&);

	template <> uisearchlist create_search_list(const uirange_pair&);
	template <> isearchlist create_search_list(const irange_pair&);
}}