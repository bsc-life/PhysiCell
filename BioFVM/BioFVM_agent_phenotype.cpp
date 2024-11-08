/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2024, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "BioFVM_agent_phenotype.h"
#include "BioFVM_basic_agent.h"

namespace BioFVM
{


Secretion::Secretion()
{
	pMicroenvironment = get_default_microenvironment(); 
	
	sync_to_current_microenvironment(); 
	return; 
}

void Secretion::sync_to_current_microenvironment( void )
{
	if( pMicroenvironment )
	{
		sync_to_microenvironment( pMicroenvironment ); 
	}
	else
	{
		secretion_rates.resize( 0 , 0.0 ); 
		uptake_rates.resize( 0 , 0.0 ); 
		saturation_densities.resize( 0 , 0.0 ); 
		net_export_rates.resize( 0, 0.0 ); 
	}
	return; 
}
	
void Secretion::sync_to_microenvironment( Microenvironment* pNew_Microenvironment )
{
	pMicroenvironment = pNew_Microenvironment;
	
	secretion_rates.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
	uptake_rates.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
	saturation_densities.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
	net_export_rates.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
	
	return; 
}

void Secretion::advance( Basic_Agent* pCell, Agent_Phenotype& phenotype , double dt )
{
	// if this phenotype is not associated with a cell, exit 
	if( pCell == NULL )
	{ return; }

	// if there is no microenvironment, attempt to sync. 
	if( pMicroenvironment == NULL )
	{
		// first, try the cell's microenvironment
		if( pCell->get_microenvironment() )
		{
			sync_to_microenvironment( pCell->get_microenvironment() ); 
		}
		// otherwise, try the default microenvironment
		else
		{
			sync_to_microenvironment( get_default_microenvironment() ); 
		}

		// if we've still failed, return. 
		if( pMicroenvironment == NULL ) 
		{
			return; 
		}
	}

	// now, call the BioFVM secretion/uptake function 
	
	pCell->simulate_secretion_and_uptake( pMicroenvironment , dt ); 
	
	return; 
}

void Secretion::set_all_secretion_to_zero( void )
{
	for( int i=0; i < secretion_rates.size(); i++ )
	{
		secretion_rates[i] = 0.0; 
		net_export_rates[i] = 0.0; 
	}
	return; 
}

void Secretion::set_all_uptake_to_zero( void )
{
	for( int i=0; i < uptake_rates.size(); i++ )
	{ uptake_rates[i] = 0.0; }
	return; 
}

void Secretion::scale_all_secretion_by_factor( double factor )
{
	for( int i=0; i < secretion_rates.size(); i++ )
	{
		secretion_rates[i] *= factor; 
		net_export_rates[i] *= factor; 
	}
	return; 
}

void Secretion::scale_all_uptake_by_factor( double factor )
{
	for( int i=0; i < uptake_rates.size(); i++ )
	{ uptake_rates[i] *= factor; }
	return; 
}

// ease of access
double& Secretion::secretion_rate( std::string name )
{
	int index = microenvironment.find_density_index(name); 
	return secretion_rates[index]; 
}

double& Secretion::uptake_rate( std::string name ) 
{
	int index = microenvironment.find_density_index(name); 
	return uptake_rates[index]; 
}

double& Secretion::saturation_density( std::string name ) 
{
	int index = microenvironment.find_density_index(name); 
	return saturation_densities[index]; 
}

double& Secretion::net_export_rate( std::string name )  
{
	int index = microenvironment.find_density_index(name); 
	return net_export_rates[index]; 
}

}
