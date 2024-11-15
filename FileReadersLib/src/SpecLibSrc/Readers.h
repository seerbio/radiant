//
// Created by andrewnichols on 11/15/24.
//

#ifndef READERS_H
#define READERS_H

#include "FileReadersLib_Exports.h"

#include <string>
#include <vector>

class FILEREADERSLIB_EXPORTS Readers {

public:

	template <class F, class T>
  	static void readVector(F &input, std::vector<T> &vec) {
		int size = 0; input.read(reinterpret_cast<char*>(&size), sizeof(int));
		if (size) {
			vec.resize(size);
			input.read(reinterpret_cast<char*>(&(vec[0])), size * sizeof(T));
		}
	}

	template<class F> void
    static readString(F &input, std::string &string) {
		int size = 0; input.read(reinterpret_cast<char*>(&size), sizeof(int));
		if (size) {
			string.resize(size);
			input.read((char*)&(string.front()), size);
		}
	}

	template <class F, class T> void
    static readArray(F &input, std::vector<T> &array) {
		int size = 0; input.read(reinterpret_cast<char*>(&size), sizeof(int));
		if (size) {
			array.resize(size);
			for (int i = 0; i < size; i++) {
				array[i].read(input);
			}
		}
	}

	template <class F, class T>
    static void readArray(F &input, std::vector<T> &array, int version) {
		int size = 0; input.read(reinterpret_cast<char*>(&size), sizeof(int));
		if (size) {
			array.resize(size);
			for (int i = 0; i < size; i++) {
				array[i].read(input, version);
			}
		}
	}

	template<class F>
    static void readStrings(F &input, std::vector<std::string> &strings) {
		int size = 0; input.read(reinterpret_cast<char*>(&size), sizeof(int));
		if (size) {
			strings.resize(size);
			for (int i = 0; i < size; i++) {
				Readers::readString(input, strings[i]);
			}
		}
	}

};


#endif //READERS_H
