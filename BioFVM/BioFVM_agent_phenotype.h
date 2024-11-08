/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#ifndef __BioFVM_agent_phenotype_h__
#define __BioFVM_agent_phenotype_h__

#include <string>

#include "BioFVM_microenvironment.h"

namespace BioFVM
{

class Agent_Phenotype;

class Secretion
{
 private:
 public:
	Microenvironment* pMicroenvironment; 
	
	std::vector<double> secretion_rates; 
	std::vector<double> uptake_rates; 
	std::vector<double> saturation_densities;
	std::vector<double> net_export_rates; 
	
	// in the default constructor, we'll size to the default microenvironment, if 
	// specified. (This ties to BioFVM.) 
	Secretion(); // done 

	// use this to properly size the secretion parameters to the microenvironment in 
	// pMicroenvironment
	void sync_to_current_microenvironment( void ); // done 
	
	void advance( Basic_Agent* pCell, Agent_Phenotype& phenotype , double dt ); 
	
	// use this to properly size the secretion parameters to the microenvironment 
	void sync_to_microenvironment( Microenvironment* pNew_Microenvironment ); // done 
	
	void set_all_secretion_to_zero( void ); // NEW
	void set_all_uptake_to_zero( void ); // NEW
	void scale_all_secretion_by_factor( double factor ); // NEW
	void scale_all_uptake_by_factor( double factor ); // NEW

	// ease of access
	double& secretion_rate( std::string name ); 
	double& uptake_rate( std::string name ); 
	double& saturation_density( std::string name ); 
	double& net_export_rate( std::string name );  	
};

class Agent_Phenotype
{
public:
    Secretion secretion;

	virtual ~Agent_Phenotype() = default;
};

}

#endif
