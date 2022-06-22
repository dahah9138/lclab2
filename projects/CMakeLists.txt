# Personal projects to add


if (CONSOLE_ONLY)
	add_subdirectory(c-elastic)
	add_subdirectory(c-rbf-fd)
else()
	add_subdirectory(rbf-fd-e)
	add_subdirectory(rbf-fd)
	add_subdirectory(FO-electric)
	add_subdirectory(film-activation)
endif()
	