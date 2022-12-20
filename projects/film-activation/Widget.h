#ifndef sandbox_WIDGET_H
#define sandbox_WIDGET_H

#include <lclab2.h>
	
using namespace Magnum;
using namespace Math::Literals;
using Key = Magnum::Platform::Sdl2Application::KeyEvent::Key;

struct StrainRecorder {
	StrainRecorder() {
		int numPoints = (end_value - start_value) / change_in_UdS;
		strain_values.reserve(numPoints);
		UdS_values.reserve(numPoints);
	}

	StrainRecorder(LC::scalar start, LC::scalar end, LC::scalar changeUdS) : start_value(start), end_value(end), change_in_UdS(changeUdS) {
		int numPoints = (end_value - start_value) / change_in_UdS;
		if (!numPoints) return;
		strain_values.reserve(numPoints);
		UdS_values.reserve(numPoints);
	}

	void Add(LC::scalar strain) {
		if (start_value == end_value) return;
		strain_values.push_back(strain);
		if (UdS_values.empty()) UdS_values.push_back(start_value);
		else UdS_values.push_back(UdS_values.back() + change_in_UdS);
	}

	void Write() {

		if (strain_values.empty() || UdS_values.empty()) return;

		std::ofstream ofile(save_file, std::ios::out | std::ios::binary);

		if (!ofile.is_open()) return;

		int size_of_scalar = LC::SIZE_OF_SCALAR;
		int size_of_int = sizeof(int);
		int N = strain_values.size();

		// Record number of entries
		ofile.write((char*)&N, size_of_int);
		ofile.write((char*)&film_thickness, size_of_scalar);
		// Record UdS values
		ofile.write((char*)&UdS_values[0], N * size_of_scalar);
		// Record strain values
		ofile.write((char*)&strain_values[0], N * size_of_scalar);

		ofile.close();
	}

	bool Finished() {
		if (start_value == end_value) return true;
		if (UdS_values.back() >= end_value) return true;
		else return false;
	}

	LC::scalar film_thickness = 0.1; // um
	LC::scalar change_in_UdS = 10.; // kPa
	LC::scalar start_value = -570.; // kPa
	LC::scalar end_value = 570.; // kPa
	LC::scalar tolerance = 1e-5;
	std::vector<LC::scalar> strain_values, UdS_values;
	std::string save_file = "D:\\dev\\lclab2\\data\\test\\film_benchmark\\benchmark1.bin";
};

struct SurfaceSelector {
	struct Axis {
		char x;
		std::array<bool, 2> draw;
	};

	void DrawGUI() {
		ImGui::TextColored({0.f, 1.f, 0.f, 1.f}, "Film surface selection");
		for (auto& ax : axes) {

			std::string ax_str = std::string(1, ax.x);

			std::string label = std::string("+") + ax_str;
			ImGui::Checkbox(label.c_str(), &ax.draw[0]);

			ImGui::SameLine();

			label = std::string("-") + ax_str;
			ImGui::Checkbox(label.c_str(), &ax.draw[1]);

		}
	}

	std::array<Axis, 3> axes = { Axis{ 'x', { true, true } }, Axis{ 'y', { true, true } }, Axis{ 'z', { true, true } } };
};

struct Widget {
	

	// Show main window
	bool showSettings = true;
	bool updateFilm = false;
	bool globalStrain = false;
	float dT;
	int updateCycle = 10;
	bool continuousUpdate = false;
	bool drawFilm = false;
	bool drawZProfile = true;
	bool updateGraphics = true;
	float pitch = 3.f; // um

	int radioZProfileAxis = 0;
	int radioZProfileIndex = -1;

	StrainRecorder strain_recorder;
	SurfaceSelector surface_selector;
};


#endif