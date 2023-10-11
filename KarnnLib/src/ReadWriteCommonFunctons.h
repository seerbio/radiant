//
// Created by anichols on 8/23/23.
//

#ifndef PYTHIADIACPP_READWRITECOMMONFUNCTONS_H
#define PYTHIADIACPP_READWRITECOMMONFUNCTONS_H

#include "KarnnLib_Exports.h"

#include <string>
#include <vector>

class KARNNLIB_EXPORTS ReadWriteCommonFunctons {

public:

    template <class F, class T>
    static void write_vector(F &out, std::vector<T> &v) {
        int size = v.size();
        out.write((char*)&size, sizeof(int));
        if (size) out.write((char*)&(v[0]), size * sizeof(T));
    }

    template <class F, class T>
    static void read_vector(F &in, std::vector<T> &v) {
        int size = 0; in.read((char*)&size, sizeof(int));
        if (size) {
            v.resize(size);
            in.read((char*)&(v[0]), size * sizeof(T));
        }
    }

    template<class F>
    static void write_string(F &out, std::string &s) {
        int size = s.size();
        out.write((char*)&size, sizeof(int));
        out.write((char*)&(s[0]), size);
    }

    template<class F>
    static void read_string(F &in, std::string &s) {
        int size = 0; in.read((char*)&size, sizeof(int));
        if (size) {
            s.resize(size);
            in.read((char*)&(s[0]), size);
        }
    }

    template <class F, class T>
    static void write_array(F &out, std::vector<T> &a) {
        int size = a.size();
        out.write((char*)&size, sizeof(int));
        for (int i = 0; i < size; i++) a[i].write(out);
    }

    template <class F, class T>
    static void read_array(F &in, std::vector<T> &a) {
        int size = 0; in.read((char*)&size, sizeof(int));
        if (size) {
            a.resize(size);
            for (int i = 0; i < size; i++) a[i].read(in);
        }
    }

    template<class F>
    static void write_strings(F &out, std::vector<std::string> &strs) {
        int size = strs.size();
        out.write((char*)&size, sizeof(int));
        for (int i = 0; i < size; i++) write_string(out, strs[i]);
    }

    template<class F>
    static void read_strings(F &in, std::vector<std::string> &strs) {
        int size = 0; in.read((char*)&size, sizeof(int));
        if (size) {
            strs.resize(size);
            for (int i = 0; i < size; i++) read_string(in, strs[i]);
        }
    }


};


#endif //PYTHIADIACPP_READWRITECOMMONFUNCTONS_H
