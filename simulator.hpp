#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <array>
#include <random>
#include "./toml11/toml.hpp"
#include "./particle.hpp"

namespace generalized_langevin {
    class Simulator {
        public:
            Simulator(const std::string& input_setup_file_path);
            void run() noexcept;
        private:
            //係数と物理定数
            double friction_coefficient;
            double coupling_coefficient;
            double K_b;

            Particle bath;//熱浴のダミー粒子
            Particle particle;
            std::ofstream out_coordinate;//粒子の座標の出力先
            std::ofstream out_energy;//エネルギーの出力先
            std::mt19937 random_engine;
            std::size_t step_num;
            std::size_t save_step_num;
            double delta_t;
            double temperature;
            std::normal_distribution<> xi_engine;
            std::array<double, 3> xi_t;
            std::array<double, 3> xi_tph;

            void step() noexcept;
            //粒子の座標と速度を求める関数
            std::array<double, 3> calculate_coordinate(Particle p) noexcept;
            std::array<double, 3> calculate_velocity(Particle p, Particle new_p) noexcept;
            //熱浴の座標と速度をランジュバン方程式に従って求める関数
            std::array<double, 3> langevin_coordinate(Particle b) noexcept;
            std::array<double, 3> langevin_velocity(Particle b) noexcept;
            void write_output() noexcept;
    };//Simulator

    

}//generalized_langevin

#endif