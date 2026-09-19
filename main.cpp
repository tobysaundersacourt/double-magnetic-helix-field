#include <cmath>
#include <numbers>
#include <iostream>
#include <list>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <fstream>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include "utils.hpp"

using scalar = double;
using vector = std::vector<scalar>;
using list = std::list<scalar>;

static const scalar pi = std::numbers::pi_v<scalar>;
static const scalar tau = 2 * pi;

struct Config {
    scalar I = 100000.0;
    scalar a = 0.0015;
    scalar b = 0.00075;
    scalar wavevector = 448.80;
    scalar outer_r = 0.003;
    scalar delta = 1e-7;
    int N = 4;
    int index_cap = 500;
    int amb_thetas = 100;
};

static std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static Config read_config(const std::string& filename) {
    Config cfg;
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Could not open config file: " + filename);
    }

    std::unordered_map<std::string, std::string> values;
    std::string line;
    while (std::getline(in, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = trim(line);
        if (line.empty()) continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        values[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }

    auto get_scalar = [&](const char* key, scalar current) {
        auto it = values.find(key);
        return it == values.end() ? current : static_cast<scalar>(std::stod(it->second));
    };
    auto get_int = [&](const char* key, int current) {
        auto it = values.find(key);
        return it == values.end() ? current : std::stoi(it->second);
    };

    cfg.I = get_scalar("I", cfg.I);
    cfg.a = get_scalar("a", cfg.a);
    cfg.b = get_scalar("b", cfg.b);
    cfg.wavevector = get_scalar("wavevector", cfg.wavevector);
    cfg.outer_r = get_scalar("outer_r", cfg.outer_r);
    cfg.delta = get_scalar("delta", cfg.delta);
    cfg.N = get_int("N", cfg.N);
    cfg.index_cap = get_int("index_cap", cfg.index_cap);
    cfg.amb_thetas = get_int("amb_thetas", cfg.amb_thetas);

    if (!(cfg.I > 0)) throw std::runtime_error("I must be positive");
    if (!(cfg.a > 0)) throw std::runtime_error("a must be positive");
    if (!(cfg.b > 0 && cfg.b < cfg.a)) throw std::runtime_error("b must satisfy 0 < b < a");
    if (std::abs(cfg.wavevector) < 1e-15) throw std::runtime_error("wavevector must be nonzero");
    if (!(cfg.outer_r > cfg.a + cfg.b)) throw std::runtime_error("outer_r must exceed a + b");
    if (!(cfg.delta > 0 && cfg.delta < cfg.b)) throw std::runtime_error("delta must satisfy 0 < delta < b");
    if (cfg.N <= 0) throw std::runtime_error("N must be positive");
    if (cfg.index_cap < cfg.N) throw std::runtime_error("index_cap must be at least N");
    if (cfg.amb_thetas < 4) throw std::runtime_error("amb_thetas must be at least 4");

    const int halfsize = static_cast<int>(cfg.b / cfg.delta);
    if (halfsize < 2) throw std::runtime_error("b/delta must be at least 2");
    return cfg;
}

scalar c_function(const scalar& r, const scalar& a, const scalar& b) {
    return utils::acos((r * r + a * a - b * b) / (2 * a * r), pi);
}

int main(int argc, char** argv) {
    try {
        const std::string config_file = argc > 1 ? argv[1] : "config.txt";
        const std::string output_file = argc > 2 ? argv[2] : "data/out.csv";
        const Config cfg = read_config(config_file);

        const scalar I = cfg.I;
        const scalar a = cfg.a;
        const scalar b = cfg.b;
        const scalar wavevector = cfg.wavevector;
        const scalar outer_r = cfg.outer_r;
        const scalar delta = cfg.delta;
        const int N = cfg.N;
        const scalar J = I / N / (pi * b * b);
        const int index_cap = cfg.index_cap;
        const int amb_thetas = cfg.amb_thetas;

        const scalar amb = a - b;
        const scalar apb = a + b;
        const int halfsize = static_cast<int>(b / delta);
        const int size = 2 * halfsize;

        vector c_cache(size, 0);
        vector r_cache(size, 0);
        // +1 prevents the original index+1 endpoint write from going out of bounds.
        vector c_integral_cache(size + 1, 0);

        list settings = {wavevector, a, b, J, static_cast<scalar>(N)};
        list coordinates;
        list thetas_out;
        list c_integrals_t;
        list field;

        const int mode_count = 1 + (index_cap - N) / N;
        int mode_index = 0;

        for (int n = N; n <= index_cap; n += N) {
            ++mode_index;
            const scalar n3 = std::pow(static_cast<scalar>(n), 3);
            const bool build = (n == N);
            vector q_cache(size, 0);
            vector w_cache(size, 0);
            scalar i_prime_integral = 0;
            scalar k_prime_integral = 0;
            vector i_prime_cache(size, 0);
            vector k_prime_cache(size, 0);
            vector xss_cache(size, 0);

            int i_prime_integral_index = 0;
            int k_prime_integral_index = size - 1;

            const scalar nk = n * wavevector;
            const scalar nk3d = std::pow(nk, 3) * delta;

            if (build) {
                r_cache[0] = amb;
                r_cache[size - 1] = apb;
            }

            vector nkr_cache(size, 0);

            while (i_prime_integral_index < halfsize) {
                if (build) {
                    c_cache[i_prime_integral_index] = c_function(r_cache[i_prime_integral_index], a, b);
                    c_cache[k_prime_integral_index] = c_function(r_cache[k_prime_integral_index], a, b);
                    c_integral_cache[i_prime_integral_index + 1] =
                        c_integral_cache[i_prime_integral_index]
                        + delta * r_cache[i_prime_integral_index] * c_cache[i_prime_integral_index];
                }

                xss_cache[i_prime_integral_index] = nk3d * std::pow(r_cache[i_prime_integral_index], 2)
                    * std::sin(n * c_cache[i_prime_integral_index]);
                xss_cache[k_prime_integral_index] = nk3d * std::pow(r_cache[k_prime_integral_index], 2)
                    * std::sin(n * c_cache[k_prime_integral_index]);

                nkr_cache[i_prime_integral_index] = nk * r_cache[i_prime_integral_index];
                nkr_cache[k_prime_integral_index] = nk * r_cache[k_prime_integral_index];

                i_prime_cache[i_prime_integral_index] = boost::math::cyl_bessel_i_prime(n, nkr_cache[i_prime_integral_index]);
                k_prime_cache[k_prime_integral_index] = boost::math::cyl_bessel_k_prime(n, nkr_cache[k_prime_integral_index]);

                i_prime_integral += xss_cache[i_prime_integral_index] * i_prime_cache[i_prime_integral_index];
                k_prime_integral += xss_cache[k_prime_integral_index] * k_prime_cache[k_prime_integral_index];

                k_prime_cache[i_prime_integral_index] = boost::math::cyl_bessel_k_prime(n, nkr_cache[i_prime_integral_index]);
                i_prime_cache[k_prime_integral_index] = boost::math::cyl_bessel_i_prime(n, nkr_cache[k_prime_integral_index]);

                q_cache[i_prime_integral_index] += k_prime_cache[i_prime_integral_index] * i_prime_integral;
                q_cache[k_prime_integral_index] += i_prime_cache[k_prime_integral_index] * k_prime_integral;

                w_cache[i_prime_integral_index] += boost::math::cyl_bessel_k(n, nkr_cache[i_prime_integral_index]) * i_prime_integral;
                w_cache[k_prime_integral_index] += boost::math::cyl_bessel_i(n, nkr_cache[k_prime_integral_index]) * k_prime_integral;

                if (build) {
                    r_cache[i_prime_integral_index + 1] = r_cache[i_prime_integral_index] + delta;
                    r_cache[k_prime_integral_index - 1] = r_cache[k_prime_integral_index] - delta;
                }

                ++i_prime_integral_index;
                --k_prime_integral_index;
            }

            while (i_prime_integral_index < size) {
                i_prime_integral += xss_cache[i_prime_integral_index] * i_prime_cache[i_prime_integral_index];
                k_prime_integral += xss_cache[k_prime_integral_index] * k_prime_cache[k_prime_integral_index];

                if (build) {
                    c_integral_cache[i_prime_integral_index + 1] =
                        c_integral_cache[i_prime_integral_index]
                        + delta * r_cache[i_prime_integral_index] * c_cache[i_prime_integral_index];
                }

                q_cache[i_prime_integral_index] += k_prime_cache[i_prime_integral_index] * i_prime_integral;
                q_cache[k_prime_integral_index] += i_prime_cache[k_prime_integral_index] * k_prime_integral;

                w_cache[i_prime_integral_index] += boost::math::cyl_bessel_k(n, nkr_cache[i_prime_integral_index]) * i_prime_integral;
                w_cache[k_prime_integral_index] += boost::math::cyl_bessel_i(n, nkr_cache[k_prime_integral_index]) * k_prime_integral;

                ++i_prime_integral_index;
                --k_prime_integral_index;
            }

            list::iterator iterator;
            if (!build) iterator = field.begin();

            scalar r = amb;
            while (r > 0) {
                const scalar nkr = n * wavevector * r;
                const scalar i_prime = boost::math::cyl_bessel_i_prime(n, nkr);
                const scalar i = boost::math::cyl_bessel_i(n, nkr);
                const scalar q = i_prime * k_prime_integral;
                const scalar w = i * k_prime_integral;
                int thetas = static_cast<int>(amb_thetas * (r / amb));
                if (thetas == 0) thetas = 1;
                const scalar delta_theta = tau / thetas;
                for (scalar theta = 0; theta < tau; theta += delta_theta) {
                    if (build) {
                        coordinates.push_back(r); coordinates.push_back(theta);
                        field.push_back(std::sin(n * theta) * q / n3);
                        field.push_back(std::cos(n * theta) * w / n3);
                        c_integrals_t.push_back(0);
                        thetas_out.push_back(static_cast<scalar>(thetas));
                    } else {
                        *iterator += std::sin(n * theta) * q / n3; ++iterator;
                        *iterator += std::cos(n * theta) * w / n3; ++iterator;
                    }
                }
                r -= delta;
            }

            r = amb;
            for (int index = 0; index < size; ++index) {
                const scalar nkr = n * wavevector * r;
                const scalar q = q_cache[index];
                const scalar w = w_cache[index];
                const int thetas = std::max(1, static_cast<int>(amb_thetas * (r / amb)));
                const scalar delta_theta = tau / thetas;
                for (scalar theta = 0; theta < tau; theta += delta_theta) {
                    if (build) {
                        coordinates.push_back(r); coordinates.push_back(theta);
                        field.push_back(std::sin(n * theta) * q / n3);
                        field.push_back(std::cos(n * theta) * w / n3);
                        c_integrals_t.push_back(c_integral_cache[index]);
                        thetas_out.push_back(static_cast<scalar>(thetas));
                    } else {
                        *iterator += std::sin(n * theta) * q / n3; ++iterator;
                        *iterator += std::cos(n * theta) * w / n3; ++iterator;
                    }
                }
                r += delta;
            }

            r = apb;
            while (r < outer_r) {
                const scalar nkr = n * wavevector * r;
                const scalar k_prime = boost::math::cyl_bessel_k_prime(n, nkr);
                const scalar k = boost::math::cyl_bessel_k(n, nkr);
                const scalar q = k_prime * i_prime_integral;
                const scalar w = k * i_prime_integral;
                const int thetas = std::max(1, static_cast<int>(amb_thetas * (r / amb)));
                const scalar delta_theta = tau / thetas;
                for (scalar theta = 0; theta < tau; theta += delta_theta) {
                    if (build) {
                        coordinates.push_back(r); coordinates.push_back(theta);
                        field.push_back(std::sin(n * theta) * q / n3);
                        field.push_back(std::cos(n * theta) * w / n3);
                        c_integrals_t.push_back(c_integral_cache[size - 1]);
                        thetas_out.push_back(static_cast<scalar>(thetas));
                    } else {
                        *iterator += std::sin(n * theta) * q / n3; ++iterator;
                        *iterator += std::cos(n * theta) * w / n3; ++iterator;
                    }
                }
                r += delta;
            }

            const int percent = static_cast<int>(100.0 * mode_index / mode_count);
            std::cout << "PROGRESS " << percent << "\n" << std::flush;
        }

        const std::filesystem::path out_path(output_file);
        if (out_path.has_parent_path()) std::filesystem::create_directories(out_path.parent_path());
        utils::clear(output_file);
        utils::append(output_file, settings);
        utils::append(output_file, coordinates);
        utils::append(output_file, c_integrals_t);
        utils::append(output_file, field);
        utils::append(output_file, thetas_out);
        std::cout << "OUTPUT " << output_file << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR " << e.what() << "\n";
        return 1;
    }
}
