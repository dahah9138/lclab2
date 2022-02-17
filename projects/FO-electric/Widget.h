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
	bool showLCINFO = false;

	Magnum::Color4 clearColor{0.5f, 0.5f, 0.5f, 0.5f};

	bool showSettings = true;
	bool showPOMSettings = false;
	bool showPreimageSettings = false;
	bool showModificationWindow = false;
	bool showNonlinearSettings = false;

	int ptheta = 0.0f;
	int pphi = 0.0f;
	float preimage_alpha = 1.0f;
	float isoLevel = 0.0625f;
	Magnum::Vector3 preimage_translate = Magnum::Vector3{ 0.0f, 0.0f, 0.0f };

	bool drawSurfaces = true;

	// For simulation
	bool relax = false;
	bool loadedFromFile = false;
	bool updateImageFromLoad = false;
	bool autoRelax = false;
	// State dependent on relax and autoRelax
	bool continuousRelax = false;

	float energyErrorThreshold = 1e-8f;


	int axis = 2;

	// Num cycles before next draw call
	int cycle = 10;

	std::array<int, 3> iPlane;
	bool midplane = true;

	bool updateImage = false;
	bool global_preimage_alpha = false;
	bool nonlinear = false;
	bool nonlinCircular = false;
	float nonlinTheta = 0.0f;
	float voltage = 0.0f;
	int voltage_iterations = 500;

	bool savePOM = false;
	bool saveNonlinear = false;
	bool sampleCharge = false;
	int radio_save_nonlin_xsection = 0;
	std::string savePOM_loc;
	std::string saveNonlin_loc;

	float computed_topological_charge = 0.0f;


	// Total en. Radio button
	int radioEn = 1;

	// Parameters that can be modified in GUI
	int topological_charge = 1;
	int npp = 20;
	std::array<float, 3> celldims = { 3.0f, 3.0f, 3.0f };
	std::array<int, 3> boundaries = { 1, 1, 0 };

	int interpolate = 1;
	std::array<int, 3> shrink_interval_begin = { 1, 1, 1 };
	std::array<int, 3> shrink_interval_end = { 2, 2, 2 };

	std::list<LC::scalar> energy_series;
	std::vector<LC::scalar> energy_series_vec;
	std::list<LC::scalar> series_x_axis;
	std::vector<LC::scalar> series_x_axis_vec;


	bool POM = false;
	float alpha = 0.75f;

	bool GPU = true;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 5.0, "um" };
};


#endif