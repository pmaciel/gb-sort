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
    midpoint_t(double _x = 0., int _i = 0) : x(_x), i(_i) {}
    double x;
    int i;

    friend std::ostream& operator<<(std::ostream& out, const midpoint_t& m) {
        out << m.x << '/' << m.i;
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


void fill_midpoints(iterator_t<std::vector<midpoint_t>>& first, const std::vector<double>& points, double a, double b,
                    int label) {
    // Note: advances first by points.size() + 1
    // Note: assumes constant increment between point coordinates
    assert(1 < points.size());

    const auto it = first;
    const auto dx = (points[1] - points.front()) / 2.;
    for (const auto& x : points) {
        *first++ = {x - dx, label};
    }

    *it      = {a, label};
    *first++ = {b, label};
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

    bool operator==(const Area& other) const {
        return N() == other.N() && W() == other.W() && S() == other.S() && E() == other.E();
    }

    bool includesNorthPole() const { return N() == 90.; }
    bool includesSouthPole() const { return S() == -90.; }
    bool isPeriodicWestEast() const { return E() == W() + 360.; }
    bool isGlobal() const { return includesNorthPole() && includesSouthPole() && isPeriodicWestEast(); }

    double N() const { return operator[](0); }
    double W() const { return operator[](1); }
    double S() const { return operator[](2); }
    double E() const { return operator[](3); }
};


static const std::string GLOBE_STR = "90/0/-90/360";
static const Area GLOBE(GLOBE_STR);


struct Grid {
    virtual ~Grid() = default;

    static Grid* build(const std::string& name, const Area&);

    size_t Nj() const { return N_.size(); }
    size_t Ni(size_t j) const { return N_[j]; }
    const Area& area() const { return area_; }

    virtual std::vector<double> Xj() const = 0;

    std::vector<double> Xi(size_t j) const {
        auto N = Ni(j);
        std::vector<double> x(N);
        linear_spacing_n(x.begin(), N, area_.W(), area_.E(), !area_.isPeriodicWestEast());
        return x;
    }

protected:
    Grid(const Area& area) : area_(area) {}

    const Area area_;
    std::vector<size_t> N_;
};


struct GaussianGrid : Grid {
protected:
    using Grid::Grid;

    std::vector<double> Xj() const override {
        auto N = Nj() / 2;
        assert(0 < N && N * 2 == Nj());

        std::vector<double> x(2 * N);
        auto a  = x.begin();
        auto b  = x.rbegin();
        auto dx = 90. / static_cast<double>(N);
        for (size_t i = 0; i < N; ++i, ++a, ++b) {
            *a = 90. - dx * (i + 0.5);  // just an approximation
            *b = -*a;
        }

        return x;
    }
};


struct OGrid : GaussianGrid {
    OGrid(size_t N, const Area& area) : GaussianGrid(area) {
        assert(area == GLOBE);  // Note: for simplicity
        assert(N > 0);

        N_.resize(2 * N);
        auto a = N_.begin();
        auto b = N_.rbegin();
        for (size_t i = 0; i < N; ++i, ++a, ++b) {
            *a = *b = 20 + i * 4;
        }
    }
};


struct FGrid : GaussianGrid {
    FGrid(size_t N, const Area& area) : GaussianGrid(area) {
        assert(area == GLOBE);  // Note: for simplicity
        assert(N > 0);
        N_.assign(2 * N, 4 * N);
    }
};


struct LLGrid : Grid {
    LLGrid(size_t Ni, size_t Nj, const Area& area) : Grid(area) {
        assert(Ni > 0 && Nj > 0);
        N_.assign(Ni, Nj);
    }

    std::vector<double> Xj() const override {
        std::vector<double> x(Nj());
        linear_spacing_n(x.begin(), Nj(), area_.N(), area_.S(), true);
        return x;
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


        parser->add_options()("h,help", "Print help");
        parser->add_options()("input-grid", "Input grid", cxxopts::value<std::string>()->default_value("O9"));
        parser->add_options()("input-area", "Input area", cxxopts::value<std::string>()->default_value(GLOBE_STR));
        parser->add_options()("output-grid", "Output grid", cxxopts::value<std::string>()->default_value("O45"));
        parser->add_options()("output-area", "Output area", cxxopts::value<std::string>()->default_value(GLOBE_STR));

        parser->parse_positional({"input-grid", "output-grid"});

        auto options = parser->parse(argc, argv);

        if (options.count("help")) {
            std::cout << parser->help() << std::endl;
            return 0;
        }


        // input and output grids

        auto Ain  = options["input-area"].as<std::string>();
        auto Aout = options["output-area"].as<std::string>();

        std::unique_ptr<Grid> Gin(Grid::build(options["input-grid"].as<std::string>(), Area(Ain)));
        std::unique_ptr<Grid> Gout(Grid::build(options["output-grid"].as<std::string>(), Area(Aout)));


        // Grid-box latitude edges (j-direction midpoints)
        std::vector<midpoint_t> Mj(Gin->Nj() + 1 + Gout->Nj() + 1);

        auto it = Mj.begin();
        fill_midpoints(it, Gin->Xj(), Gin->area().N(), Gin->area().S(), 0);
        fill_midpoints(it, Gout->Xj(), Gout->area().N(), Gout->area().S(), 1);
        assert(it == Mj.end());

        print(Mj);
        std::cout << "---" << std::endl;

        std::inplace_merge(
            Mj.begin(), Mj.begin() + Gin->Nj() + 1, Mj.end(),
            [](const midpoint_t& a, const midpoint_t& b) { return a.x > b.x || (a.x == b.x && a.i < b.i); });

        print(Mj);
        std::cout << "---" << std::endl;
    }
    catch (const cxxopts::exceptions::exception& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
