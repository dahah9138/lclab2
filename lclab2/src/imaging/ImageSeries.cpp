#include "ImageSeries.h"

namespace LC { namespace Imaging {
	
	ImageSeries::ImageSeries(std::int32_t w, std::int32_t h, std::string file) : width(w), height(h), write_file(file), count(0) {
		bmp = BMP(w, h);
	}
	
	void ImageSeries::GenerateFrame(COLOR *data, std::int32_t w, std::int32_t h, bool order) {
		
		if (w != width || h != height) {
			width = w;
			height = h;
			bmp = BMP(w, h);
		}
		
		// Create the frame
		for (std::int32_t x = 0; x < width; x++) {
			for (std::int32_t y = 0; y < height; y++) {
				int index = order ? y + height * x : x + width * y;
				bmp.set_pixel(x, y, data[index].B, data[index].G, data[index].R, data[index].A);
			}
		}
		
	}
	
	void ImageSeries::WriteFrame() {
		std::string name = write_file + "_" + std::to_string(++count) + ".bmp";
		bmp.write(name.c_str());
	}
	
	void ImageSeries::GenerateAndWriteFrame(COLOR *data, std::int32_t w, std::int32_t h, bool order) {
		GenerateFrame(data, w, h, order);
		WriteFrame();
	}

	void ImageSeries::operator = (const ImageSeries& rhs) {
		width = rhs.width;
		height = rhs.height;
		bmp = rhs.bmp;
		write_file = rhs.write_file;
		count = rhs.count;
	}
	
}}