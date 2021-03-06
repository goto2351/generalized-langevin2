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
            double equilibrium_length;

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
            std::array<double, 3> calculate_coordinate(Particle p, Particle b) noexcept;
            std::array<double, 3> calculate_velocity(Particle p, Particle new_p, Particle b) noexcept;
            //熱浴の座標と速度をランジュバン方程式に従って求める関数
            std::array<double, 3> langevin_coordinate(Particle b) noexcept;
            std::array<double, 3> langevin_velocity(Particle b) noexcept;
            void write_output(std::size_t step_index) noexcept;
    };//Simulator

    Simulator::Simulator(const std::string& input_setup_file_path) {
        const auto input_setup_file = toml::parse(input_setup_file_path);

        //定数の読み込み
        friction_coefficient = toml::find<double>(input_setup_file, "constants", "friction_coefficient");
        coupling_coefficient = toml::find<double>(input_setup_file, "constants", "coupling_coefficient");
        K_b = toml::find<double>(input_setup_file, "constants", "K_b");
        equilibrium_length = toml::find<double>(input_setup_file, "constants", "equilibrium_length");

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
        out_energy.open(working_path + "/" + project_name + "_energy.csv");
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

    void Simulator::run() noexcept {
        write_output(0);
        for (std::size_t step_index = 1; step_index <= step_num; ++step_index) {
            step();
            if (step_index%save_step_num == 0) {
                write_output(step_index);
            }
        }
    }

    void Simulator::step() noexcept {
        Particle new_bath = bath;
        Particle new_particle = particle;

        const auto [new_bath_x, new_bath_y, new_bath_z] = langevin_coordinate(bath);
        new_bath.x = new_bath_x;
        new_bath.y = new_bath_y;
        new_bath.z = new_bath_z;
        const auto [new_bath_vx, new_bath_vy, new_bath_vz] = langevin_velocity(bath);
        new_bath.vx = new_bath_vx;
        new_bath.vy = new_bath_vy;
        new_bath.vz = new_bath_vz;

        const auto [new_particle_x, new_particle_y, new_particle_z] = calculate_coordinate(particle, bath);
        new_particle.x = new_particle_x;
        new_particle.y = new_particle_y;
        new_particle.z = new_particle_z;
        const auto [new_particle_vx, new_particle_vy, new_particle_vz] = calculate_velocity(particle, new_particle, bath);
        new_particle.vx = new_particle_vx;
        new_particle.vy = new_particle_vy;
        new_particle.vz = new_particle_vz;

        bath = new_bath;
        particle = new_particle;

        xi_t = xi_tph;
        xi_tph = {
            xi_engine(random_engine),
            xi_engine(random_engine),
            xi_engine(random_engine)
        };
    }

    std::array<double, 3> Simulator::langevin_coordinate(Particle b) noexcept {
        //速度Verlet法で熱浴の次の時刻の座標を求める
        const double next_x = b.x + b.vx*delta_t*(1.0-(friction_coefficient*delta_t)/2.0) + ((delta_t*delta_t)/2.0)*xi_t[0];
        const double next_y = b.y + b.vy*delta_t*(1.0-(friction_coefficient*delta_t)/2.0) + ((delta_t*delta_t)/2.0)*xi_t[1];
        const double next_z = b.z + b.vz*delta_t*(1.0-(friction_coefficient*delta_t)/2.0) + ((delta_t*delta_t)/2.0)*xi_t[2];

        return {next_x, next_y, next_z};
    }

    std::array<double, 3> Simulator::langevin_velocity(Particle b) noexcept {
        const double term1 = 1.0 - (friction_coefficient*delta_t)/2.0;
        const double term2 = 1.0 - (friction_coefficient*delta_t)/2.0 + ((friction_coefficient*delta_t)/2.0)*((friction_coefficient*delta_t)/2.0);

        const double next_vx = b.vx*term1*term2 + (delta_t/2.0)*term2*(xi_t[0] + xi_tph[0]);
        const double next_vy = b.vx*term1*term2 + (delta_t/2.0)*term2*(xi_t[1] + xi_tph[1]);
        const double next_vz = b.vx*term1*term2 + (delta_t/2.0)*term2*(xi_t[2] + xi_tph[2]);

        return {next_vx, next_vy, next_vz};
    }

    std::array<double, 3> Simulator::calculate_coordinate(Particle p, Particle b) noexcept {
        std::array<double, 3> vec;
        vec[0] = p.x - b.x;
        vec[1] = p.y - b.y;
        vec[2] = p.x - b.z;
        const double distance = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
        
        const double term1 = (-1.0*coupling_coefficient*(distance - equilibrium_length)) /distance;
        std::array<double,3> f;
        f[0] = term1*vec[0];
        f[1] = term1*vec[1];
        f[2] = term1*vec[2];

        const double next_x = p.x + p.vx*delta_t + (f[0]/p.mass)*delta_t*delta_t/2.0;
        const double next_y = p.y + p.vy*delta_t + (f[1]/p.mass)*delta_t*delta_t/2.0;
        const double next_z = p.z + p.vz*delta_t + (f[2]/p.mass)*delta_t*delta_t/2.0;

        return {next_x, next_y, next_z};
    }

    std::array<double, 3> Simulator::calculate_velocity(Particle p, Particle new_p, Particle b) noexcept {
        std::array<double, 3> vec;
        vec[0] = p.x - b.x;
        vec[1] = p.y - b.y;
        vec[2] = p.x - b.z;
        const double distance = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);  
        const double term1 = (-1.0*coupling_coefficient*(distance - equilibrium_length)) /distance;
        std::array<double,3> f;
        f[0] = term1*vec[0];
        f[1] = term1*vec[1];
        f[2] = term1*vec[2];

        std::array<double, 3> next_vec;
        next_vec[0] = new_p.x - b.x;
        next_vec[1] = new_p.y - b.y;
        next_vec[2] = new_p.z - b.z;
        const double next_distance = std::sqrt(next_vec[0]*next_vec[0] + next_vec[1]*next_vec[1] + next_vec[2]*next_vec[2]);
        const double next_term1 = (-1.0*coupling_coefficient*(next_distance - equilibrium_length)) /next_distance;
        std::array<double, 3> next_f;
        next_f[0] = next_term1*next_vec[0];
        next_f[1] = next_term1*next_vec[1];
        next_f[2] = next_term1*next_vec[2];

        const double next_vx = p.vx + (delta_t/2.0)*(next_f[0] + f[0])/p.mass;
        const double next_vy = p.vy + (delta_t/2.0)*(next_f[1] + f[1])/p.mass;
        const double next_vz = p.vz + (delta_t/2.0)*(next_f[2] + f[2])/p.mass;

        return {next_vx, next_vy, next_vz};
    }

    void  Simulator::write_output(std::size_t step_index) noexcept {
        //座標の書き出し
        out_coordinate << "2" << std::endl;
        out_coordinate << std::endl;
        out_coordinate << "H " << bath.x << " " << bath.y << " " << bath.z << std::endl;
        out_coordinate << "C " << particle.x << " " << particle.y << " " << particle.z << std::endl;

        //運動エネルギーの書き出し
        double kinetic_energy = particle.mass*(std::pow(particle.vx, 2.0) + std::pow(particle.vy, 2.0) + std::pow(particle.vz, 2.0))/2.0;
        out_energy << step_index << "," << kinetic_energy << std::endl;
    }

}//generalized_langevin

#endif