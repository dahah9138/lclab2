Building:
<Dependencies>
- magnum <https://github.com/mosra/magnum>
- magnum-plugins <https://github.com/mosra/magnum-plugins>
- magnum-integration <https://github.com/mosra/magnum-integration>
- corrade <https://github.com/mosra/corrade>
- imgui <https://github.com/ocornut/imgui>
- implot <https://github.com/epezent/implot>
- portable-file-dialogs <https://github.com/samhocevar/portable-file-dialogs>
- spdlog <https://github.com/gabime/spdlog>
- Eigen <https://github.com/libigl/eigen>
- SDL2 <https://github.com/libsdl-org/SDL>
- hemi <https://github.com/harrism/hemi>
- CUDA
- cmake
- vtk

<Implemented code>
- kNN <https://github.com/vincentfpgarcia/kNN-CUDA>


* For macOS and linux
- If you do not have SDL installed follow the instructions below
	# 	sudo pacman -S sdl2             # on ArchLinux
	# 	sudo apt install libsdl2-dev    # on Ubuntu / Debian
	# 	brew install sdl2               # on macOS (via Homebrew)
	
* On Windows
- I recommend installing MSVC 2022 and 2026

<Directions>
1. 'git clone' this project into where you want it
2. 'cd' into the root directory
3. 'git submodule update --init --recursive'
4. Create a build folder in the root directory
5. 'cd' to build
6. 'cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..'
[Unix]
7. 'make'
[Windows]
7. 'cmake --build . --config Release'
Congratulations!


==== 6/7/2026 ====
Attempting to implement vtk into foel.exe

1. Created folder 'C:\vtk-worksace'
2. cd to vtk-workspace
3. Run the command git clone --branch v9.3.0 https://gitlab.kitware.com/vtk/vtk.git source
4. Build with cuda flags using the command:

cmake ../source ^
  -DCMAKE_INSTALL_PREFIX="C:/local/VTK" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_CXX_STANDARD=17 ^
  -DVTKm_ENABLE_CUDA=ON ^
  -DVTK_USE_CUDA=ON ^
  -DVTK_MODULE_ENABLE_VTK_ParallelDIY=NO
  
5. Compile with the command (after -j is number of cores, adjust as needed):
	cmake --build . --config Release -j 8
6. Install the build with:
	cmake --install . --config Release
7. Back in the lclab2 build directory, run:
	cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DVTK_DIR="C:/local/VTK/lib/cmake/vtk-9.3"