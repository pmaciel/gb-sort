/*
 * Copyright (c) 2022 Pedro Maciel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cxxopts.hpp"


template <class T>
void print(const T& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<typename T::value_type>(std::cout, " "));
    std::cout << std::endl;
}


struct midpoint_t {
    midpoint_t(double _point, int _meta = 0) : point(_point), meta(_meta) {}
    double point;
    int meta;

    bool operator<(const midpoint_t& m) const { return point < m.point || (point == m.point && meta < m.meta); }
    friend std::ostream& operator<<(std::ostream& out, const midpoint_t& m) {
        out << m.point << '/' << m.meta;
        return out;
    }
};


#if 0
struct Area : std::array<double, 4> {
    Area(double _N, double _W, double _S, double _E) : std::array<double, 4>{_N, _W, _S, _E} {
        assert(-90. <= S() && S() <= N() && N() <= 90.);
        assert(W() <= E() && E() <= W() + 360.);
    }
    Area() : Area(90., 0., -90., 360.) {}
    double N() const { return operator[](0); }
    double W() const { return operator[](1); }
    double S() const { return operator[](2); }
    double E() const { return operator[](3); }
};
#endif


struct Grid {
    size_t Nj() const { return Nj_; }
    size_t Ni(size_t j) const { return Ni_[j]; }

    static Grid* build(const std::string& name);

protected:
    Grid(std::vector<size_t>&& Ni) : Nj_(Ni.size()), Ni_(Ni) {}

    std::vector<size_t> Ni_;
    size_t Nj_;
};


struct ReducedGrid : Grid {
    ReducedGrid(std::vector<size_t>&& Ni) : Grid(std::move(Ni)) {}
};


struct RegularGrid : Grid {
    RegularGrid(size_t Ni, size_t Nj) : Grid(std::vector<size_t>(Nj, Ni)) {}
};


Grid* Grid::build(const std::string& name) {
    std::smatch match;  // Note: first sub_match is the whole string

    const static std::regex octahedral("[Oo]([1-9][0-9]*)");
    const static std::regex regular_gg("[Ff]([1-9][0-9]*)");
    const static std::regex regular_ll("[L]([1-9][0-9]*)x([1-9][0-9]*)");

    if (std::regex_match(name, match, octahedral)) {
        assert(match.size() == 2);
        auto N = static_cast<size_t>(std::stol(match[1].str()));
        assert(N > 0);

        std::vector<size_t> Ni(2 * N);
        auto a = Ni.begin();
        auto b = Ni.rbegin();
        for (size_t i = 0; i < N; ++i) {
            *(a++) = *(b++) = 20 + i * 4;
        }

        return new ReducedGrid(std::move(Ni));
    }

    if (std::regex_match(name, match, regular_gg)) {
        assert(match.size() == 2);
        auto N = static_cast<size_t>(std::stol(match[1].str()));
        assert(N > 0);

        return new RegularGrid(N * 4, N * 2);
    }

    if (std::regex_match(name, match, regular_ll)) {
        assert(match.size() == 3);
        auto Ni = static_cast<size_t>(std::stol(match[1].str()));
        auto Nj = static_cast<size_t>(std::stol(match[2].str()));
        assert(Ni > 0 && Nj > 0);

        return new RegularGrid(Ni, Nj);
    }

    throw std::runtime_error("Unrecognized grid '" + name + "'");
}


int main(int argc, const char* argv[]) {
    try {
        // options
        std::unique_ptr<cxxopts::Options> parser(
            new cxxopts::Options(argv[0], " - grid-box intersections interpolation method"));

        parser->add_options()("h,help", "Print help");
        parser->add_options()("i,input", "Input grid", cxxopts::value<std::string>()->default_value("O12"));
        parser->add_options()("o,output", "Output grid", cxxopts::value<std::string>()->default_value("O6"));
        parser->add_options()("I,input-area", "Input grid area",
                              cxxopts::value<std::string>()->default_value("90/0/-90/360"));
        parser->add_options()("O,output-area", "Output grid area",
                              cxxopts::value<std::string>()->default_value("90/0/-90/360"));

        parser->parse_positional({"input", "output"});

        auto options = parser->parse(argc, argv);

        if (options.count("help")) {
            std::cout << parser->help() << std::endl;
            return 0;
        }

        // build input and output grids

        std::unique_ptr<Grid> Gi(Grid::build(options["input"].as<std::string>()));
        std::unique_ptr<Grid> Go(Grid::build(options["output"].as<std::string>()));

        // build test data to sort

        std::vector<midpoint_t> m;
        for (double i : {1.1, 1.2}) {
            for (int j : {3, 1}) {
                m.emplace_back(midpoint_t{i, j});
            }
        }

        print(m);
        std::sort(m.begin(), m.end());
        print(m);

        std::vector<double> v = {1, 3, 5, 2, 4, 6};
        print(v);

        std::inplace_merge(v.begin(), v.begin() + 3, v.end());
        print(v);
    }
    catch (const cxxopts::exceptions::exception& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
