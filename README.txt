Building:
<Dependencies>
- magnum <https://github.com/mosra/magnum>
- magnum-plugins <https://github.com/mosra/magnum-plugins>
- corrade <https://github.com/mosra/corrade>
- Eigen <https://github.com/libigl/eigen>
- SDL

* For macOS and linux
- If you do not have SDL installed follow the instructions below
	# 	sudo pacman -S sdl2             # on ArchLinux
	# 	sudo apt install libsdl2-dev    # on Ubuntu / Debian
	# 	brew install sdl2               # on macOS (via Homebrew)

<Directions>
1. 'git clone' this project into where you want it
2. 'cd' into the root directory
3. 'git submodule update --init'
4. Create a build folder in the root directory
5. 'cd' to build
6. 'cmake ..'
[Unix]
7. 'make'
[Windows]
7. 'cmake --build .'

Congratulations!