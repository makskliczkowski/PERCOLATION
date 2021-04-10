#pragma once
#ifndef GENERAL_H
#define GENERAL_H

#include <time.h>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <string>
#include <omp.h>
#include <numeric>
#include <filesystem>
#include <future>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include "RandomAlghorithms.h"

// for different systems
static const char* kPathSeparator =
#ifdef _WIN32
"\\";
#else
"/";
#endif

namespace general {

	typedef std::vector<std::vector<std::vector<std::tuple<int, int>>>> neighborsStructure;

	/* FUNCTIONS */
	template <typename T>
	std::string to_string_prec(const T a_value, const int n = 3) {
		std::ostringstream out;
		out.precision(n);
		out << std::fixed << a_value;
		return out.str();
	}
	std::vector<std::string> split_str(std::string s, std::string delimiter);



	/// <summary>
	/// Virutal lattice class, it will provide us with a strict geometry that we need
	/// LATER - THINK ABOUT HOW TO MAKE IT MORE GENERAL AND PROVIDE VISUALIZATION HOW THE LATTICE ELEMENTS ARE NUMBERED
	/// </summary>
	class lattice2D {
	public:
		enum lattice_types {
			square,
			triangle,
			hexagonal
		};
	protected:
		/* lattice parameters */
		int Ns;
		lattice2D::lattice_types type;
		neighborsStructure nearest_neighbours;
		neighborsStructure next_nearest_neighbours;
	public:
		virtual ~lattice2D() = default;
		virtual std::unique_ptr<general::lattice2D> clone() const = 0; // pure virtual clone
		virtual std::unique_ptr<general::lattice2D> move_clone() = 0;
		/* Virtual Getters */
		virtual int get_Lx() = 0;
		virtual int get_Ly() = 0;
		/* Getters */
		int get_Ns();

		lattice2D::lattice_types get_type();
		std::string getString_type();
		std::tuple<int,int> get_nn(int x,int y, int nei_num);
		std::tuple<int, int> get_nnn(int x,int y, int nei_num);
		int get_nn_number(int x,int y);
		int get_nnn_number(int x,int y);
		/* Calculators */
		virtual void calculate_nn_pbc() = 0;
		virtual void calculate_nnn_pbc() = 0;
		virtual void calculate_nn() = 0;
		virtual void calculate_nnn() = 0;
	};
	/// <summary>
	/// Square lattice inherits from lattice
	/// </summary>
	class square_lattice :public virtual lattice2D {
	private:
		int Lx;
		int Ly;
	public:
		square_lattice();
		square_lattice(int Lx, int Ly);
		~square_lattice() final;
		square_lattice(const square_lattice& A);
		square_lattice(square_lattice&& A) noexcept;
		/* Getters */
		int get_Lx() override;
		int get_Ly() override;
		/* Calculators */
		void calculate_nn_pbc() override;
		void calculate_nnn_pbc() override;
		void calculate_nn() override;
		void calculate_nnn() override;

		/* Clone functions */
		virtual std::unique_ptr<general::lattice2D> clone() const override {
			return std::unique_ptr<general::lattice2D>(new general::square_lattice(*this));
		}
		virtual std::unique_ptr <general::lattice2D > move_clone() override {
			return std::unique_ptr<general::lattice2D>(new general::square_lattice(std::move(*this)));
		}
	};
	/// <summary>
	/// Triangular lattice
	/// Dimensions are again given the same as in square
	/// </summary>
	class triangle_lattice :public virtual lattice2D {
	private:
		int Lx;
		int Ly;
	public:
		triangle_lattice();
		triangle_lattice(int Lx, int Ly);
		~triangle_lattice() final;
		triangle_lattice(const triangle_lattice& A);
		triangle_lattice(triangle_lattice&& A) noexcept;
		/* Getters */
		int get_Lx() override;
		int get_Ly() override;
		/* Calculators */
		void calculate_nn_pbc() override;
		void calculate_nnn_pbc() override;
		void calculate_nn() override;
		void calculate_nnn() override;

		/* Clone functions */
		virtual std::unique_ptr<general::lattice2D> clone() const override {
			return std::unique_ptr<general::lattice2D>(new general::triangle_lattice(*this));
		}
		virtual std::unique_ptr <general::lattice2D > move_clone() override {
			return std::unique_ptr<general::lattice2D>(new general::triangle_lattice(std::move(*this)));
		}
	};


	/// <summary>
	/// Map for lattices
	/// </summary>
	std::unordered_map < std::string, general::lattice2D::lattice_types > const lattice_parse_table = {
		{"square",general::lattice2D::square},
		{"triangle",general::lattice2D::triangle},
		{"hexagonal",general::lattice2D::hexagonal},
	};
}

#endif