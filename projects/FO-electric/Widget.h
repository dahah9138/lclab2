#ifndef FOELECTRIC_WIDGET_H
#define FOELECTRIC_WIDGET_H

#include <lclab2.h>
#include <list>

using namespace Magnum;
using namespace Math::Literals;
using Key = Magnum::Platform::Sdl2Application::KeyEvent::Key;

struct Widget {

	// From example
	bool showDemoWindow = false;
	bool showAnotherWindow = false;

	bool showSettings = true;

	// For simulation
	bool relax = false;
	bool loadedFromFile = false;
	bool updateImageFromLoad = false;
	bool autoRelax = false;
	// State dependent on relax and autoRelax
	bool continuousRelax = false;

	int axis = 2;

	// Num cycles before next draw call
	int cycle = 10;

	std::array<int, 3> iPlane;
	bool midplane = true;

	bool updateImage = false;
	bool nonlinear = false;
	float nonlinTheta = 0.0f;
	float voltage = 0.0f;
	int voltage_iterations = 500;

	// Parameters that can be modified in GUI
	int topological_charge = 1;
	int npp = 30;
	std::array<float, 3> celldims = { 3.0f, 3.0f, 3.0f };
	std::array<int, 3> boundaries = { 1, 1, 0 };

	int interpolate = 1;
	std::array<int, 3> shrink_interval_begin = { 0, 0, 0 };
	std::array<int, 3> shrink_interval_end = { 1, 1, 1 };

	std::list<LC::scalar> energy_series;
	std::vector<LC::scalar> energy_series_vec;
	std::list<LC::scalar> series_x_axis;
	std::vector<LC::scalar> series_x_axis_vec;


	bool POM = false;
	float alpha = 0.5f;

	bool GPU = false;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 5.0, "um" };
};


#endif