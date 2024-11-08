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

class Molecular
{
	private:
	public: 
		Microenvironment* pMicroenvironment; 
	
		// model much of this from Secretion 
		Molecular(); 
 	
		// we'll set this to replace BioFVM's version		
		std::vector<double> internalized_total_substrates; 

		// for each substrate, a fraction 0 <= f <= 1 of the 
		// total internalized substrate is released back inot
		// the environment at death 
		std::vector<double> fraction_released_at_death; 

		// for each substrate, a fraction 0 <= f <= 1 of the 
		// total internalized substrate is transferred to the  
		// predatory cell when ingested 
		std::vector<double> fraction_transferred_when_ingested; 
		
		/* prototyping / beta in 1.5.0 */ 
		// Boolean, Integer, and Double parameters
/*		
		std::vector<bool> bools; 
		std::unordered_map<std::string,int> bool_name_map; 
		std::string& bool_name( int i ); 
		std::vector<std::string> bool_units; 
		void resize_bools( int n ); 
		int add_bool( std::string name , std::string units , bool value ); 
		bool& access_bool( std::string name ); 
		
		std::vector<int> ints; 
		std::unordered_map<std::string,int> int_name_map; 
		std::string& int_name( int i ); 
		std::vector<std::string> int_units; 
		int& access_int( std::string name ); 
		
		std::vector<int> doubles; 
		std::unordered_map<std::string,int> double_name_map; 
		std::string& double_name( int i ); 
		std::vector<std::string> double_units; 
		double& access_double( std::string name ); 
*/
	
		// use this to properly size the secretion parameters to the 
		// microenvironment in molecular.pMicroenvironment. 
		void sync_to_current_microenvironment( void ); // done 
		
//		void advance( Basic_Agent* pCell, Phenotype& phenotype , double dt ); 
		
		// use this to properly size the secretion parameters to the microenvironment in 
		// pMicroenvironment
		void sync_to_microenvironment( Microenvironment* pNew_Microenvironment ); // done 
		
		// ease of access 
		double&  internalized_total_substrate( std::string name ); 
		
};

class Agent_Phenotype
{
public:
    Secretion secretion;
	Molecular molecular;

	virtual ~Agent_Phenotype() = default;
};

}

#endif
