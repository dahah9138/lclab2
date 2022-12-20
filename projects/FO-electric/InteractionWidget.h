#ifndef INTERACTION_WIDGET_H
#define INTERACTION_WIDGET_H

#include <lclab2.h>

struct InteractionWidget {
	InteractionWidget() {
		if (interaction_fname.capacity() < interaction_fname_buffer_size)
			interaction_fname.reserve(interaction_fname_buffer_size);
	}
	// Display Interaction GUI
	void Display() {
		ImGui::PushItemWidth(100.f);
		ImGui::InputFloat("Separation Distance", &separationDistance);
		ImGui::InputInt("Theta points", &interactionThetaPoints);
		ImGui::SameLine();
		ImGui::InputInt("Phi points", &interactionPhiPoints);
		ImGui::SameLine();
		ImGui::InputInt("Node density (per pitch)", &interactionNPP);
		ImGui::Checkbox("Force cell size", &forceInteractionCellSize);
		ImGui::SameLine();
		ImGui::InputFloat("Cell size", &interactionCellSize);
		ImGui::InputInt("Relaxation iterations", &interactionIterations);
		ImGui::PopItemWidth();
	}
	// Heliknoton separation distance (relative to center of mass)
	float separationDistance = 2.5f;
	// Number of polar points
	int interactionThetaPoints = 1;
	// Number of azimuthal points
	int interactionPhiPoints = 10;
	// Node density of interaction volume
	int interactionNPP = 20;
	// (1) force cell size, (0) automatically determine cell size
	bool forceInteractionCellSize = 0;
	// Size of interaction cell (if forced)
	float interactionCellSize = 6.f;
	// How much to relax by
	int interactionIterations = 200;
	std::string interaction_fname = "interaction.bin";
	const size_t interaction_fname_buffer_size = 250;
};




#endif