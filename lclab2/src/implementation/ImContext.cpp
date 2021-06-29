#include "ImContext.h"

#include <stdio.h>

namespace LC { namespace ImPlotIntegration {
	
	ImPlotContext* CreateContext() {
		return ImPlot::CreateContext();
	}
	
}}