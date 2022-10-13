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

    Area() : Area(90., 0., -90., 360.) {}

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
    size_t Nj() const { return NN_.size(); }
    size_t Ni(size_t j) const { return NN_[j]; }

    static Grid* build(const std::string& name, const Area&);

protected:
    Grid(const Area& area) : area_(area) {}

    Area area_;
    std::vector<size_t> NN_;
};


struct OGrid : Grid {
    OGrid(size_t N, const Area& area) : Grid(area) {
        assert(N > 0);

        NN_.resize(2 * N);
        auto a = NN_.begin();
        auto b = NN_.rbegin();
        for (size_t i = 0; i < N; ++i) {
            *(a++) = *(b++) = 20 + i * 4;
        }
    }
};


struct FGrid : Grid {
    FGrid(size_t N, const Area& area) : Grid(area) {
        assert(N > 0);
        NN_.assign(2 * N, 4 * N);
    }
};


struct LLGrid : Grid {
    LLGrid(size_t Ni, size_t Nj, const Area& area) : Grid(area) {
        assert(Ni > 0 && Nj > 0);
        NN_.assign(Ni, Nj);
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
        parser->add_options()("input-grid", "Input grid", cxxopts::value<std::string>()->default_value("O12"));
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
