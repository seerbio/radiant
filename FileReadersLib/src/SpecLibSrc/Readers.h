//
// Created by andrewnichols on 11/15/24.
//

#ifndef READERS_H
#define READERS_H

#include "FileReadersLib_Exports.h"

#include <string>
#include <vector>
#include <ios>
#include <QFile>
#include <algorithm>
#include <cstring>
#include <cstdio>

/**
 * Sequential view over a memory mapped file.  It deliberately exposes the
 * small subset of std::istream used by the speclib format readers, so the
 * format parser does not need a second implementation.
 */
class FILEREADERSLIB_EXPORTS MappedFileInput {
public:
    explicit MappedFileInput(const QString &filePath) : m_file(filePath) {
        if (!m_file.open(QIODevice::ReadOnly)) {
            return;
        }
        m_size = m_file.size();
        if (m_size > 0) {
            m_data = m_file.map(0, m_size);
            if (m_data == nullptr) {
                m_file.close();
            }
        }
    }

    ~MappedFileInput() {
        if (m_data != nullptr) {
            m_file.unmap(m_data);
        }
    }

    MappedFileInput(const MappedFileInput &) = delete;
    MappedFileInput &operator=(const MappedFileInput &) = delete;

    bool isOpen() const { return m_file.isOpen(); }
    std::streamsize read(char *destination, std::streamsize count) {
        const std::streamsize available = static_cast<std::streamsize>(m_size - m_offset);
        const std::streamsize amount = std::min(count, std::max<std::streamsize>(0, available));
        if (amount > 0) {
            std::memcpy(destination, m_data + m_offset, static_cast<size_t>(amount));
            m_offset += static_cast<qint64>(amount);
        }
        return amount;
    }
    int peek() const {
        return m_offset < m_size ? static_cast<unsigned char>(m_data[m_offset]) : EOF;
    }

private:
    QFile m_file;
    uchar *m_data = nullptr;
    qint64 m_size = 0;
    qint64 m_offset = 0;
};

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

	template<class F>
    static void readString(F &input, std::string &string) {
		int size = 0; input.read(reinterpret_cast<char*>(&size), sizeof(int));
		if (size) {
			string.resize(size);
			input.read((char*)&(string.front()), size);
		}
	}

	template <class F, class T>
    static void readArray(F &input, std::vector<T> &array) {
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
