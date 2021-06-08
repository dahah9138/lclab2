Building:
<Dependencies>
- magnum <https://github.com/mosra/magnum.git>
- magnum-plugins <https://github.com/mosra/magnum-plugins.git>
- corrade <https://github.com/mosra/corrade.git>
- SDL

* For macOS and linux
- If you do not have SDL installed follow the instructions below
	1. 	sudo pacman -S sdl2             # on ArchLinux
	2. 	sudo apt install libsdl2-dev    # on Ubuntu / Debian
	3. 	brew install sdl2               # on macOS (via Homebrew)

<Directions>
1. git clone the latest versions of magnum and corrade into external
2. git clone the latest version of magnum-plugins
3. Create a build folder in the root directory
4. cd to build
5. cmake ..
[Unix]
6. make
[Windows]
6. cmake --build .

Congratulations!
