Building:
<Dependencies>
- magnum <https://github.com/mosra/magnum.git>
- magnum-plugins <https://github.com/mosra/magnum-plugins.git>
- corrade <https://github.com/mosra/corrade.git>
- Eigen <https://github.com/libigl/eigen.git>
- SDL

* For macOS and linux
- If you do not have SDL installed follow the instructions below
	# 	sudo pacman -S sdl2             # on ArchLinux
	# 	sudo apt install libsdl2-dev    # on Ubuntu / Debian
	# 	brew install sdl2               # on macOS (via Homebrew)

<Directions>
1. git clone the latest versions of magnum and corrade into external
2. git clone the latest version of eigen into external
3. git clone the latest version of magnum-plugins
4. Create a build folder in the root directory
5. cd to build
6. cmake ..
[Unix]
7. make
[Windows]
7. cmake --build .

Congratulations!
