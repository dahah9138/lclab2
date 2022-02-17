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

<Implemented code>
- kNN <https://github.com/vincentfpgarcia/kNN-CUDA>


* For macOS and linux
- If you do not have SDL installed follow the instructions below
	# 	sudo pacman -S sdl2             # on ArchLinux
	# 	sudo apt install libsdl2-dev    # on Ubuntu / Debian
	# 	brew install sdl2               # on macOS (via Homebrew)

<Directions>
1. 'git clone' this project into where you want it
2. 'cd' into the root directory
3. 'git submodule update --init --recursive'
4. Create a build folder in the root directory
5. 'cd' to build
6. 'cmake ..'
[Unix]
7. 'make'
[Windows]
7. 'cmake --build .'

Congratulations!