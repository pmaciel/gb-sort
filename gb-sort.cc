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
#include <iostream>
#include <iterator>
#include <memory>
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


int main(int argc, const char* argv[]) {
    try {
        // options
        std::unique_ptr<cxxopts::Options> parser(
            new cxxopts::Options(argv[0], " - grid-box intersections interpolation method"));

        parser->add_options()("help", "Print help");
        parser->add_options()("i,input", "Input grid", cxxopts::value<std::string>()->default_value("O12"));
        parser->add_options()("o,output", "Output grid", cxxopts::value<std::string>()->default_value("O6"));

        parser->parse_positional({"input", "output"});

        auto options = parser->parse(argc, argv);

        if (options.count("help")) {
            std::cout << parser->help() << std::endl;
            return 0;
        }

        auto input  = options["input"].as<std::string>();
        auto output = options["output"].as<std::string>();


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
