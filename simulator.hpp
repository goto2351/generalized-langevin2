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

    Simulator::Simulator(const std::string& input_setup_file_path) {
        const auto input_setup_file = toml::parse(input_setup_file_path);

        //定数の読み込み
        friction_coefficient = toml::find<double>(input_setup_file, "constants", "friction_coefficient");
        coupling_coefficient = toml::find<double>(input_setup_file, "constants", "coupling_coefficient");
        K_b = toml::find<double>(input_setup_file, "constants", "K_b");

        //熱浴と粒子の初期化
        bath.x = toml::find<double>(input_setup_file, "bath", "x");
        bath.y = toml::find<double>(input_setup_file, "bath", "y");
        bath.z = toml::find<double>(input_setup_file, "bath", "z");
        bath.vx = 0.0;
        bath.vy = 0.0;
        bath.vz = 0.0;
        bath.mass = toml::find<double>(input_setup_file, "bath", "mass");
        particle.x = toml::find<double>(input_setup_file, "particle", "x");
        particle.y = toml::find<double>(input_setup_file, "particle", "y");
        particle.z = toml::find<double>(input_setup_file, "particle", "z");
        particle.vx = 0.0;
        particle.vy = 0.0;
        particle.vz = 0.0;
        particle.mass = toml::find<double>(input_setup_file, "particle", "mass");

        //アウトプットファイルを開く
        const auto project_name = toml::find<std::string>(input_setup_file, "meta_data", "project_name");
        const auto working_path = toml::find<std::string>(input_setup_file, "meta_data", "working_path");
        out_coordinate.open(working_path + "/" + project_name + ".xyz");
        if(!out_coordinate) {
            std::cerr << "cannot open:" << working_path + "/" + project_name + ".xyz" << std::endl;
            std::exit(1);
        }
        out_energy.open(working_path + "/" + project_name + "_energy.txt");
        if(!out_energy) {
            std::cerr << "cannot open:" << working_path + "/" + project_name + "_energy.txt" << std::endl;
            std::exit(1);
        }

        const auto random_seed = toml::find<std::size_t>(input_setup_file, "meta_data", "random_seed");
        random_engine.seed(random_seed);

        step_num = toml::find<std::size_t>(input_setup_file, "meta_data", "step_num");
        save_step_num = toml::find<std::size_t>(input_setup_file, "meta_data", "save_step_num");
        temperature = toml::find<double>(input_setup_file, "meta_data", "temperature");
        delta_t = toml::find<double>(input_setup_file, "meta_data", "delta_t");

        std::normal_distribution<> init_xi_engine(0.0, std::sqrt((2.0*friction_coefficient*K_b*temperature*delta_t)/bath.mass));
        xi_engine = init_xi_engine;

        xi_t = {
            xi_engine(random_engine),
            xi_engine(random_engine),
            xi_engine(random_engine)
        };
        xi_tph = {
            xi_engine(random_engine),
            xi_engine(random_engine),
            xi_engine(random_engine)
        };
    }//constructor

}//generalized_langevin

#endif