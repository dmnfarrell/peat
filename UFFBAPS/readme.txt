/*
#
# UFFBAPS - Unified Force Field for Binding And Protein Stability
# Copyright (C) 2010 Jens Erik Nielsen & Chresten Soendergaard
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information: 
# Email: Jens.Nielsen_at_gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
*/
How to install and run:

- Install c++ boost (http://www.boost.org/)

- Changed /usr/include/boost/numeric/ublas/lu.hpp (or whatever path
  you installed to) from:

#ifndef BOOST_UBLAS_LU_H
#define BOOST_UBLAS_LU_H

to:

#ifndef BOOST_UBLAS_LU_H
#define BOOST_UBLAS_LU_H

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>

I know this sucks but I couldn't get to work otherwise...

- Run make in the trunk dir - you might have to install the executable
  manually if permissions prevents make from doing it.

- Copy the run_dir (with sub dirs) to the location you want run the
  program in

- Edit run.cfg in your run dir to tell uffbaps which structures to
  use. There are some examples to get you started.

- Run the executable as 'Uffbaps run.cfg'

- Results are printed to standard out and octave files 'stability.m'
  and 'binding.m' for mutant stability changes and protein-ligand
  binding affinity respectively.

