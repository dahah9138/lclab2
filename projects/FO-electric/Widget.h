#ifndef FOELECTRIC_WIDGET_H
#define FOELECTRIC_WIDGET_H

#include <lclab2.h>
#include <filesystem>
#include <list>

using namespace Magnum;
using namespace Math::Literals;
using Key = Magnum::Platform::Sdl2Application::KeyEvent::Key;


struct Widget {
	struct KnotInteraction {
		/*
			Position selector to choose positions where to spawn solitons
		*/
	struct PositionGui {
		std::vector<Vector3> position_array;
		Vector3 pot_pos;
		void SubGUI() {
			ImGui::InputFloat3("Position##PositionGUI", &pot_pos[0]);

			if (ImGui::Button("Add position##PositionGUI")) {
				position_array.emplace_back(pot_pos);
			}
			ImGui::Separator();
			// Display initialized positions
			int ctr = 0;
			bool triggered = false;
			int trigger_index = -1;
			for (auto& p : position_array) {
				++ctr;
				ImGui::Text("[%d]: (%f, %f, %f)", ctr, p[0], p[1], p[2]);
				ImGui::SameLine();
				if (ImGui::Button(std::string("Delete##PositionGUI" + std::to_string(ctr)).c_str()) && !triggered) {
					triggered = true;
					trigger_index = ctr - 1;
				}
			}

			if (triggered)
				position_array.erase(position_array.begin() + trigger_index);

		}

	};

	struct EFieldSwitch {
		// Value to switch the Efield to
		float efield = 0.f;
		// number of frames to apply this efield
		int dframes = 0;
	};

	struct VectorPreimage {
		int theta = 0;
		int phi = 0;
	};

		bool showWindow = false;
		PositionGui pos_gui;
		// Knot color [default = red]
		Magnum::Color4 knotColor{ 1.f, 0.f, 0.f, 1.f };
		// Isovalue for knot surface
		float isoValue = 0.58f;
		float gradSval = 0.15f;
		float tolerance = 0.0f;
		int node_density = 35;
		// Upsampling for background grid to compute vortex lines
		float upsampling_multiplier = 2.0f;
		// Upper (magnitude) band threshold
		float sb_pos_iso_ratio_upper = 100.f;
		float sb_neg_iso_ratio_upper = 100.f;
		// Lower (magnitude) band threshold
		float sb_pos_iso_ratio_lower = 0.1f;
		float sb_neg_iso_ratio_lower = 0.1f;

		std::vector<EFieldSwitch> switching_data;
		std::vector< VectorPreimage > preimages = { {0,0}, {180,0} };
		float preimage_isovalue = 0.062;
		// Continously alternate between
		bool cycleEfield = 0;

		int saturation = 2;
		int sb_saturation = 2;
		float point_density = 1.f;
		
		// Start the interaction
		bool start = false;

		// Start from initial conditions
		int useInitialConditions = 0;

		// Initial interaction conditions
		float theta0 = 45.f; // pi/4
		float phi0 = 180.f;
		float seperation = 1.8f;
		std::array<float, 3> CELL = { 6.f, 6.f, 5.f };

		// File location
		std::string file_loc = LCLAB2_ROOT_PATH + std::string("/data/knot/");
		// File name
		std::string file, subdirectory;

		// Use global path
		bool global_path = true;

		const static size_t fname_buffer_size = 128;
		char fname_buffer[fname_buffer_size];

		const static size_t fsubdir_buffer_size = 128;
		char fsubdir_buffer[fsubdir_buffer_size];
		
		// Number of [frames] to save
		int nFrames = 1;
		// frame offset if more frames are desired
		int nFrameOffset = 0;
		// Number of iterations to prerelax the simulation volume
		int nPrerelax = 0;
		// Relaxation iterations per frame
		int nRelax_per_frame = 0;
		// Backup rate (number of frames between backing up data)
		int nBackup_rate = 25;
		// Relaxation rate (negative is underrelaxation)
		float relaxRate = -0.2f;
		// Color convention
		// 0 - nematic
		// 1 - chi field
		int color_convention = 0;
		bool saveStringKnots = 1;

		bool GUI() {
			if (!showWindow)
				return showWindow;

			ImGui::Begin("Knot interaction##Knot-interaction", &showWindow);

			ImGui::PushItemWidth(200.f);
			ImGui::TextColored(ImVec4(0.0f,1.0f,0.0f,1.0f), "Path: %s", PathToFile().c_str());
			ImGui::InputText("Subdirectory##Knot-interaction", fsubdir_buffer, fsubdir_buffer_size);
			ImGui::PopItemWidth();
			ImGui::SameLine();
			ImGui::Checkbox("Global path", &global_path);

			// Copy buffer to string
			subdirectory = std::string(fsubdir_buffer);
			ImGui::PushItemWidth(200.f);
			ImGui::InputText("File name##Knot-interaction", fname_buffer, fname_buffer_size);
			ImGui::PopItemWidth();
			ImGui::SameLine();
			ImGui::Checkbox("Save string knots", &saveStringKnots);

			// Copy buffer to string
			file = std::string(fname_buffer);

			ImGui::PushItemWidth(100.f);
			ImGui::InputInt("Node density##Knot-interaction", &node_density);
			ImGui::SameLine();
			ImGui::InputFloat("Upsampling##Knot-interaction", &upsampling_multiplier);
			ImGui::InputInt("Backup rate##Knot-interaction", &nBackup_rate);
			ImGui::PopItemWidth();
			ImGui::PushItemWidth(50.f);
			//ImGui::InputFloat("Knot isovalue##Knot-interaction", &isoValue);
			//ImGui::InputFloat("grad(S) isovalue##Knot-interaction", &gradSval);
			//ImGui::InputFloat("SB iso upper bd (p)##Knot-interaction", &sb_pos_iso_ratio_upper);
			//ImGui::InputFloat("SB iso lower bd (p)##Knot-interaction", &sb_pos_iso_ratio_lower);
			//ImGui::InputFloat("SB iso upper bd (n)##Knot-interaction", &sb_neg_iso_ratio_upper);
			//ImGui::InputFloat("SB iso lower bd (n)##Knot-interaction", &sb_neg_iso_ratio_lower);
			//ImGui::InputFloat("Chi field tolerance##Knot-interaction", &tolerance);
			ImGui::PushItemWidth(70.f);
			ImGui::InputFloat("Theta (deg)##Knot-interation", &theta0);
			ImGui::SameLine();
			ImGui::InputFloat("Phi (deg)##Knot-interation", &phi0);
			ImGui::PopItemWidth();
			ImGui::SameLine();
			ImGui::InputFloat("Separation##Knot-interaction", &seperation);
			ImGui::InputFloat("Knot point density##Knot-interaction", &point_density);

			ImGui::Separator();
			ImGui::Text("Construct Efield sequence");
			ImGui::SameLine();
			ImGui::Checkbox("Cycle##KNot-interaction-efield", &cycleEfield);
			if (ImGui::Button("Add Switch##Knot-interaction-field")) {
				// Check that the duration is not zero
				switching_data.push_back({ 0.f, 1 });
			}

			std::vector<int> remove_switch;

			// Display current data
			int total_frames = 0;
			for (int i = 0; i < switching_data.size(); i++) {

				if (switching_data[i].dframes <= 0)
					switching_data[i].dframes = 1;

				std::string str_id = "S" + std::to_string(i);
				std::string str_efield = std::string("Efield##Knot-interaction-field") + std::to_string(i);
				std::string str_duration = std::string("Duration##Knot-interaction-field") + std::to_string(i);
				std::string str_delete = std::string("X##Knot-interaction-field") + std::to_string(i);
				ImGui::Text(str_id.c_str());
				ImGui::SameLine();
				ImGui::PushItemWidth(100.f);
				ImGui::InputFloat(str_efield.c_str(), &switching_data[i].efield);
				ImGui::SameLine();
				ImGui::InputInt(str_duration.c_str(), &switching_data[i].dframes);
				ImGui::PopItemWidth();
				ImGui::SameLine();
				
				if (ImGui::Button(str_delete.c_str())) {
					remove_switch.push_back(i);
				}
				else { // Add to the total frames
					total_frames += switching_data[i].dframes;
				}
			}

			// parse through switches to remove one at a time
			while (!remove_switch.empty()) {
				switching_data.erase(switching_data.begin() + remove_switch.back());
				remove_switch.pop_back();
			}

			if (ImGui::Button("Update frame count##Knot-interaction")) {
				nFrames = total_frames;
			}

			// Calculate frames
			ImGui::Separator();

			ImGui::InputFloat("Relax rate##Knot-interaction", &relaxRate);
			ImGui::PopItemWidth();
			ImGui::PushItemWidth(100.f);
			ImGui::InputInt("Saturation##Knot-interaction", &saturation);
			ImGui::InputInt("Saturation (SB)##knot-interaction", &sb_saturation);
			ImGui::InputInt("Iterations (prerelax)##Knot-interaction", &nPrerelax);
			ImGui::InputInt("Iterations/frame##Knot-interaction", &nRelax_per_frame);
			ImGui::InputInt("Frames##Knot-interaction", &nFrames);
			ImGui::InputInt("Frame Offset##Knot-interaction", &nFrameOffset);
			ImGui::Text("Coloring");
			ImGui::SameLine();
			ImGui::RadioButton("Nematic##Knot-interaction", &color_convention, 0);
			ImGui::SameLine();
			ImGui::RadioButton("Chi##Knot-interaction", &color_convention, 1);
			ImGui::SameLine();
			ImGui::RadioButton("Handedness##Knot-interaction", &color_convention, 2);
			ImGui::PopItemWidth();

			// Check whether to use initial conditions or not
			// True -> generate using thet0 and phi0
			// False -> use what is already in the simulation scope
			ImGui::Text("Initialization");
			ImGui::Separator();
			ImGui::RadioButton("Preset##Knot-interaction", &useInitialConditions, 0);
			ImGui::SameLine();
			ImGui::RadioButton("Dimer##Knot-interaction", &useInitialConditions, 1);
			ImGui::SameLine();
			ImGui::RadioButton("Generalized##Knot-interaction", &useInitialConditions, 2);
			ImGui::SameLine();
			ImGui::RadioButton("Dimerizer##Knot-interaction", &useInitialConditions, 3);

			if (useInitialConditions == 2 || useInitialConditions == 3) {
				ImGui::InputFloat3("Cell size##Knot-interaction", &CELL[0]);
				pos_gui.SubGUI();
			}

			if (ImGui::Button("Generate##Knot-interaction"))
				start = true;

			ImGui::End();
			// Return the status of showWindow
			return showWindow;
		}

		bool Dispatch() {
			if (!start) return false;
			else {
				// Set start to false and begin dispatch
				start = false;
				return true;
			}
		}

		std::array<LC::scalar, 3> GetCell() const {
			return { CELL[0], CELL[1], CELL[2] };
		}

		std::string FileName() {
			std::string result;
			if (global_path) { // Use subdirectory as the global path
				result = subdirectory + "/" + file;
			}
			else { // Use local path from build
				result = file_loc + subdirectory + "/" + file;
			}
			return result;
		}

		std::string PathToFile() {
			std::string result;
			if (global_path) { // Use subdirectory as the global path
				result = subdirectory;
			}
			else { // Use local path from build
				result = file_loc + subdirectory;
			}
			return result;
		}

		std::string InfoFile(std::string inf = "info.txt") {
			std::string result;
			if (global_path) { // Use subdirectory as the global path
				result = subdirectory + "/" + inf;
			}
			else { // Use local path from build
				result = file_loc + subdirectory + "/" + inf;
			}
			return result;
		}

	};


	// From example
	bool showDemoWindow = false;
	bool showLCINFO = false;

	Magnum::Color4 clearColor{0.5f, 0.5f, 0.5f, 0.5f};

	float specular = 0.f;
	float diffuse = 0.1f;
	float ambient = 1.f;

	LC::NematicArray::DrawType nematicDrawType = LC::NematicArray::DrawType::Cone;

	bool showSettings = true;
	bool showPOMSettings = false;
	bool showPreimageSettings = false;
	bool showModificationWindow = false;
	bool showInteractionEnergyWindow = false;
	bool showNonlinearSettings = false;
	bool showVortexKnotSettings = false;
	bool couplePlaneAndNematic = true;
	int chiColorScheme = true;
	int S2colors = 1;
	std::array<float, 3> pionComponents = { 0.f, 0.f, 1.0f };

	bool showZProfileWindow = false;
	bool drawProfile = false;
	bool generateProfile = false;
	float thermalRandomization = 0.0f;
	float expansionCoeff = 0.f;
	float contractionCoeff = 1.0f;
	int radioZProfileAxis = 0;
	int radioZProfileIndex = -1;
	float temperature = 293.0f;
	float omegabar = 0.5f;
	float inversion_temp = 314.0f;
	int chain_units = 5;

	bool generateKnots = false;
	int maxVortexComponents = 1;
	int minComponentSize = 5;
	float knotCompletionDist = 1.5f;
	int knotRefinementIterations = 0;
	int knotInitialCutoffIterations = 15;

	float lehmanArclength = 0.628f; // (2 * pi / 10)
	float lehman_dz = 0.1f;
	float total_lehman_z_dist = 1.f; // 1 pitch
	float total_lehman_forward_dist = 0.3f; // 0.3 pitch
	int lehman_upsample = 3;

	int tilt_angle = 0;
	int tilt_direction = 0;
	int translationNumber = 1;
	bool helicalTranslation = true;
	float helicalLayerOffset = 0.0f;
	float separationDistance = 2.f;
	int interactionThetaPoints = 1;
	int interactionPhiPoints = 30;
	int interactionIterations = 200;
	bool singleInteractionHeliknoton = false;
	float interactionOmegaOffset = 0.f;
	int interactionNPP = 20;
	bool forceInteractionCellSize = true;
	std::array<float, 3> interactionCellSize = { 6.f,6.f,6.f };
	bool interactionSymmetry = 1;

	int ptheta = 0.0f;
	int pphi = 0.0f;
	int smoothingIterations = 20;
	// Laplacian smoothing parameter
	float smoothingAlpha = 0.25f;
	// Taubin smoothing parameters
	float smoothingLambda = 0.35f;// 0.33f;
	float smoothingMu = -0.34f;
	int smoothingType = 2;
	float preimage_alpha = 1.0f;
	float isoLevel = 0.0625f;
	Magnum::Vector3 preimage_translate = Magnum::Vector3{ 0.0f, 0.0f, 0.0f };

	bool drawSurfaces = true;
	bool showCellBd = false;

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
	bool showVisualProperties = false;
	bool lehmann_tool_window = false;
	bool lehmann_reset_helical_bg = true;
	bool saveNonlinear = false;
	bool sampleCharge = false;
	int radio_save_nonlin_xsection = 0;
	std::string savePOM_loc;
	std::string saveNonlin_loc;

	float computed_topological_charge = 0.0f;


	// Total en. Radio button
	// 0 - total energy
	// 1 - functional derivative en
	int radioEn = 0;

	// Parameters that can be modified in GUI
	int topological_charge = 1;
	int npp = 20;
	int chirality = 1;
	std::array<float, 3> celldims = { 3.0f, 3.0f, 3.0f };
	std::array<int, 3> boundaries = { 1, 1, 0 };
	int max_graph_points = 10000;

	// 0 = heliknoton, 1 = hopfion
	int hopfion_type = 0;

	int interpolate = 1;
	std::array<int, 3> shrink_interval_begin = { 1, 1, 1 };
	std::array<int, 3> heliknotonInsertLocation = { 1, 1, 1 };
	bool hideHeliknotonInsertIndicator = true;
	std::array<int, 3> shrink_interval_end = { 2, 2, 2 };

	std::list<LC::precision_scalar> energy_series;
	std::vector<LC::precision_scalar> energy_series_vec;
	std::list<LC::precision_scalar> series_x_axis;
	std::vector<LC::precision_scalar> series_x_axis_vec;


	bool POM = false;
	float alpha = 1.0f;

	// Default pitch in micrometers
	LC::SIscalar pitch = { 4.5, "um" };

	KnotInteraction knot_interaction_handle;

	bool multiplane_window = false;
};


#endif