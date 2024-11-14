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

Molecular::Molecular()
{
	pMicroenvironment = get_default_microenvironment(); 
	sync_to_current_microenvironment(); 

	return; 
}

void Molecular::sync_to_current_microenvironment( void )
{
	if( pMicroenvironment )
	{
		sync_to_microenvironment( pMicroenvironment ); 
	}
	else
	{
		internalized_total_substrates.resize( 0 , 0.0 ); 
		fraction_released_at_death.resize( 0 , 0.0 ); 
		fraction_transferred_when_ingested.resize( 0, 0.0 ); 
	}
	return; 
}
	
void Molecular::sync_to_microenvironment( Microenvironment* pNew_Microenvironment )
{
	pMicroenvironment = pNew_Microenvironment;
	
	int number_of_densities = pMicroenvironment->number_of_densities() ; 

	internalized_total_substrates.resize( number_of_densities , 0.0 ); 
	fraction_released_at_death.resize( number_of_densities , 0.0 ); 
	fraction_transferred_when_ingested.resize( number_of_densities , 0.0 ); 
	
	return; 
}

// ease of access 
double&  Molecular::internalized_total_substrate( std::string name )
{
	int index = microenvironment.find_density_index(name); 
	return internalized_total_substrates[index]; 
}

/*
void Molecular::advance( Basic_Agent* pCell, Phenotype& phenotype , double dt )
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

	// make sure the associated cell has the correct rate vectors 
	if( pCell->internalized_substrates != &internalized_substrates )
	{
		// copy the data over 
		internalized_substrates = *(pCell->internalized_substrates);
		// remove the BioFVM copy 
		delete pCell->internalized_substrates; 
		// point BioFVM to this one  
		pCell->internalized_substrates = &internalized_substrates; 
	}

	// now, call the functions 
//	if( pCell->functions.internal_substrate_function )
//	{ pCell->functions.internal_substrate_function( pCell,phenotype,dt);  }
//	if( pCell->functions.molecular_model_function )
//	{ pCell->functions.molecular_model_function( pCell,phenotype,dt);  }


	return; 
}
*/

Motility::Motility()
{
	is_motile = false; 
	
	persistence_time = 1.0;
	migration_speed = 1.0;
	
	migration_bias_direction.resize( 3 , 0.0 ); 
	migration_bias = 0.0; 
		
	restrict_to_2D = false; 
	
	// update_migration_bias_direction = NULL; 
	
	motility_vector.resize( 3 , 0.0 ); 
	
	chemotaxis_index = 0; 
	chemotaxis_direction = 1; 
	
	sync_to_current_microenvironment(); 
	
	return; 
}

void Motility::sync_to_current_microenvironment( void )
{
	Microenvironment* pMicroenvironment = get_default_microenvironment(); 
	if( pMicroenvironment )
	{ sync_to_microenvironment( pMicroenvironment ); } 
	else
	{ chemotactic_sensitivities.resize( 1 , 0.0 ); }

	return; 
}

void Motility::sync_to_microenvironment( Microenvironment* pNew_Microenvironment )
{
	chemotactic_sensitivities.resize( pNew_Microenvironment->number_of_densities() , 0.0 ); 
	return; 
}

double& Motility::chemotactic_sensitivity( std::string name )
{
	int n = microenvironment.find_density_index(name); 
	return chemotactic_sensitivities[n]; 
}

Volume::Volume()
{
	// reference parameter values for MCF-7, in cubic microns 
	fluid_fraction = 0.75;  

	total = 2494; 
	fluid = fluid_fraction * total; 
	solid = total-fluid; 
	
	nuclear = 540.0;

	
	nuclear_fluid = fluid_fraction * nuclear; 
	nuclear_solid = nuclear - nuclear_fluid;

	cytoplasmic = total - nuclear;
	cytoplasmic_fluid = fluid_fraction*cytoplasmic; 
	cytoplasmic_solid = cytoplasmic - cytoplasmic_fluid; 
	
	// rates are in units of 1/min 
	cytoplasmic_biomass_change_rate = 0.27 / 60.0; 
	nuclear_biomass_change_rate = 0.33 / 60.0; 
	fluid_change_rate = 3.0 / 60.0;

	calcified_fraction = 0.0;
	
	calcification_rate = 0.0; 
	
	target_solid_cytoplasmic = cytoplasmic_solid;
	target_solid_nuclear = nuclear_solid;
	target_fluid_fraction = fluid_fraction;
	
	cytoplasmic_to_nuclear_ratio = cytoplasmic / ( 1e-16 + nuclear);
	target_cytoplasmic_to_nuclear_ratio = cytoplasmic_to_nuclear_ratio; 
	
	// the cell bursts at these volumes 
	relative_rupture_volume = 2.0; 
		// as fraction of volume at entry to the current phase
	rupture_volume = relative_rupture_volume * total; // in volume units 
	
	return; 
};


void Volume::multiply_by_ratio( double ratio )
{
	total *= ratio;
	solid *= ratio;
	fluid *= ratio;
	
	nuclear *= ratio;
	nuclear_fluid *= ratio;
	nuclear_solid *= ratio;
	
	cytoplasmic *= ratio;
	cytoplasmic_fluid *= ratio;
	cytoplasmic_solid *= ratio;
	
	rupture_volume *= ratio; 
	
	target_solid_nuclear *= ratio;
	target_solid_cytoplasmic *= ratio; 	
	
	return; 
}

void Volume::divide( void )
{
	multiply_by_ratio( 0.5 ); 
	return; 
}


Geometry::Geometry()
{
	// reference values for MCF-7, based on 
	// volume = 2494 cubic microns
	// nuclear volume = 540 cubic microns 
	radius = 8.412710547954228; 
	nuclear_radius = 5.051670902881889; 
	surface_area = 889.3685284131693; 
	
	polarity = 0.0; 
	return; 
}

void Geometry::update_radius( Basic_Agent* pCell, Agent_Phenotype& phenotype, double dt )
{
	static double four_thirds_pi =  4.188790204786391;
	radius = phenotype.volume.total; 
	radius /= four_thirds_pi; 
	radius = pow( radius , 0.333333333333333333333333333333333333333 ); 
	return; 
}

void Geometry::update_nuclear_radius( Basic_Agent* pCell, Agent_Phenotype& phenotype, double dt )
{
	static double four_thirds_pi = 4.188790204786391;
	nuclear_radius = phenotype.volume.nuclear; 
	nuclear_radius /= four_thirds_pi; 
	nuclear_radius = pow( nuclear_radius , 0.333333333333333333333333333333333333333 ); 
	return; 
}

void Geometry::update_surface_area( Basic_Agent* pCell, Agent_Phenotype& phenotype, double dt )
{
	// 4pi / (4pi/3)^(2/3)
	static double the_constant = 4.835975862049409; 
	surface_area = pow( phenotype.volume.total , 0.666666666666667 );
	surface_area /= the_constant; 
	
	return; 
}

void Geometry::update( Basic_Agent* pCell, Agent_Phenotype& phenotype, double dt )
{
	update_radius(pCell,phenotype,dt); 
	update_nuclear_radius(pCell,phenotype,dt);
	
	// surface area = 4*pi*r^2 = (4/3)*pi*r^3 / (r/3)	
	surface_area = phenotype.volume.total; 
	surface_area /= radius; 
	surface_area *= 3.0; 
	return; 
}

}