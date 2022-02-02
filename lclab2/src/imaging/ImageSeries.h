#ifndef IMAGE_SERIES_H
#define IMAGE_SERIES_H

#include "BMP.h"
#include <string>

namespace LC { namespace Imaging {
	
	struct ImageSeries {
		struct COLOR {
			uint8_t R,G,B,A;
		};
		ImageSeries() = default;
		ImageSeries(std::int32_t w, std::int32_t h, std::string file = "SERIES");

		void operator = (const ImageSeries& rhs);
		
		// order 0 means (x,y) maps to x + X * y in data
		// order 1 means (x,y) maps to y + Y * x in data
		void GenerateFrame(COLOR *data, std::int32_t w, std::int32_t h, bool order = 0);
		void WriteFrame();
		void GenerateAndWriteFrame(COLOR *data, std::int32_t w, std::int32_t h, bool order = 0);
		
		std::int32_t width{ 0 };
		std::int32_t height{ 0 };
		BMP bmp;
		std::string write_file {"SERIES"};
		std::int32_t count{ 0 };
	};
	
}}


#endif