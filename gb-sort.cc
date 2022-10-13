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


using spacing_t = std::vector<double>;

template <typename T>
using iterator_t = decltype(std::begin(std::declval<T&>()));


template <class ForwardIt>
void linear_spacing_n(ForwardIt first, size_t count, double _a, double _b, bool endpoint) {
    assert(1 < count && _a != _b);
    const auto dx = (_b - _a) / static_cast<double>(count - (endpoint ? 1 : 0));

    for (size_t i = 0; i < count; ++i, ++first) {
        *first = _a + i * dx;
    }
};


struct Area : protected std::array<double, 4> {
    Area(double _N, double _W, double _S, double _E) : std::array<double, 4>{_N, _W, _S, _E} { check(); }

    Area(const std::string& NWSE) {
        const static std::string x = "([+-]?([0-9]+([.][0-9]*)?|[.][0-9]+))";
        const static std::regex NWSEx(x + "/" + x + "/" + x + "/" + x);

        std::smatch match;  // Note: first sub_match is the whole string
        assert(std::regex_match(NWSE, match, NWSEx));
        assert(match.size() == 13);

        operator[](0) = std::stod(match[1].str());
        operator[](1) = std::stod(match[4].str());
        operator[](2) = std::stod(match[7].str());
        operator[](3) = std::stod(match[10].str());

        check();
    }

    void check() const {
        assert(-90. <= S() && S() <= N() && N() <= 90.);
        assert(W() <= E() && E() <= W() + 360.);
    }

    double N() const { return operator[](0); }
    double W() const { return operator[](1); }
    double S() const { return operator[](2); }
    double E() const { return operator[](3); }
};


struct Grid {
    static Grid* build(const std::string& name, const Area&);

    size_t Nj() const { return Nn_.size(); }
    size_t Ni(size_t j) const { return Nn_[j]; }
    const Area& area() const { return area_; }

protected:
    Grid(const Area& area) : area_(area) {}

    const Area area_;
    std::vector<size_t> Nn_;
    std::vector<double> Xj_;
};


struct OGrid : Grid {
    OGrid(size_t N, const Area& area) : Grid(area) {
        assert(N > 0);

        Nn_.resize(2 * N);
        auto Na = Nn_.begin();
        auto Nb = Nn_.rbegin();
        for (size_t i = 0; i < N; ++i, ++Na, ++Nb) {
            *Na = *Nb = 20 + i * 4;
        }

        Xj_.resize(2 * N);
        auto Xa = Xj_.begin();
        auto Xb = Xj_.rbegin();
        auto dx = 90. / static_cast<double>(N);
        for (size_t i = 0; i < N; ++i, ++Xa, ++Xb) {
            *Xa = 90. - dx * (i + 0.5);  // just an approximation
            *Xb = -*Xa;
        }
    }
};


struct FGrid : Grid {
    FGrid(size_t N, const Area& area) : Grid(area) {
        assert(N > 0);
        Nn_.assign(2 * N, 4 * N);

        Xj_.resize(2 * N);
        auto Xa = Xj_.begin();
        auto Xb = Xj_.rbegin();
        auto dx = 90. / static_cast<double>(N);
        for (size_t i = 0; i < N; ++i, ++Xa, ++Xb) {
            *Xa = 90. - dx * (i + 0.5);  // just an approximation
            *Xb = -*Xa;
        }
    }
};


struct LLGrid : Grid {
    LLGrid(size_t Ni, size_t Nj, const Area& area) : Grid(area) {
        assert(Ni > 0 && Nj > 0);
        Nn_.assign(Ni, Nj);

        Xj_.resize(Nj);
        linear_spacing_n(Xj_.begin(), Nj, area.N(), area.S(), true);
    }
};


Grid* Grid::build(const std::string& name, const Area& area) {
    const static std::regex octahedral("[Oo]([1-9][0-9]*)");
    const static std::regex regular_gg("[Ff]([1-9][0-9]*)");
    const static std::regex regular_ll("LL([1-9][0-9]*)x([1-9][0-9]*)");

    std::smatch match;  // Note: first sub_match is the whole string
    if (std::regex_match(name, match, octahedral)) {
        assert(match.size() == 2);
        return new OGrid(static_cast<size_t>(std::stol(match[1].str())), area);
    }

    if (std::regex_match(name, match, regular_gg)) {
        assert(match.size() == 2);
        return new FGrid(static_cast<size_t>(std::stol(match[1].str())), area);
    }

    if (std::regex_match(name, match, regular_ll)) {
        assert(match.size() == 3);
        auto Ni = static_cast<size_t>(std::stol(match[1].str()));
        auto Nj = static_cast<size_t>(std::stol(match[2].str()));
        return new LLGrid(Ni, Nj, area);
    }

    throw std::runtime_error("Unrecognized grid '" + name + "'");
}


int main(int argc, const char* argv[]) {
    try {
        // options
        std::unique_ptr<cxxopts::Options> parser(
            new cxxopts::Options(argv[0], " - grid-box intersections interpolation method"));


        const std::string globe = "90/0/-90/360";
        parser->add_options()("h,help", "Print help");
        parser->add_options()("input-grid", "Input grid", cxxopts::value<std::string>()->default_value("LL3x3"));
        parser->add_options()("input-area", "Input area", cxxopts::value<std::string>()->default_value(globe));
        parser->add_options()("output-grid", "Output grid", cxxopts::value<std::string>()->default_value("O6"));
        parser->add_options()("output-area", "Output area", cxxopts::value<std::string>()->default_value(globe));

        parser->parse_positional({"input-grid", "output-grid"});

        auto options = parser->parse(argc, argv);

        if (options.count("help")) {
            std::cout << parser->help() << std::endl;
            return 0;
        }

        // build input and output grids

        auto Ai = options["input-area"].as<std::string>();
        auto Ao = options["output-area"].as<std::string>();

        std::unique_ptr<Grid> Gi(Grid::build(options["input-grid"].as<std::string>(), Area(Ai)));
        std::unique_ptr<Grid> Go(Grid::build(options["output-grid"].as<std::string>(), Area(Ao)));

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
