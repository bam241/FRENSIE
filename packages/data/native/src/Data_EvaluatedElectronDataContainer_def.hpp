//---------------------------------------------------------------------------//
//!
//! \file   Data_EvaluatedElectronDataContainer_def.hpp
//! \author Luke Kersting
//! \brief  The native eedl container template defs.
//!
//---------------------------------------------------------------------------//

#ifndef DATA_EVALUATED_ELECTRON_DATA_CONTAINER_DEF_HPP
#define DATA_EVALUATED_ELECTRON_DATA_CONTAINER_DEF_HPP

// Boost Includes
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/nvp.hpp>

namespace Data{

// Save the data to an archive
template<typename Archive>
void EvaluatedElectronDataContainer::save( Archive& ar, 
						  const unsigned version) const
{
  ar & boost::serialization::make_nvp( "atomic_number", d_atomic_number );
  ar & boost::serialization::make_nvp( "subshells", d_subshells );
  ar & boost::serialization::make_nvp( "subshell_occupancies", 
				       d_subshell_occupancies );
  ar & boost::serialization::make_nvp( "subshell_binding_energies", 
				       d_subshell_binding_energies );
  ar & boost::serialization::make_nvp( "relaxation_transitions", 
				       d_relaxation_transitions );
  ar & boost::serialization::make_nvp( "relaxation_vacancies", 
				       d_relaxation_vacancies );
  ar & boost::serialization::make_nvp( "relaxation_particle_energies", 
				       d_relaxation_particle_energies );
  ar & boost::serialization::make_nvp( "relaxation_probabilities", 
				       d_relaxation_probabilities );
  ar & boost::serialization::make_nvp( "cutoff_angle", 
                                        d_cutoff_angle );
  ar & boost::serialization::make_nvp( "angular_energy_grid",
				                        d_angular_energy_grid );
/*
  ar & boost::serialization::make_nvp( "number_of_discrete_angles", 
                                        d_number_of_discrete_angles );
  ar & boost::serialization::make_nvp( "hard_elastic_angles", 
                                        d_hard_elastic_angles );
  ar & boost::serialization::make_nvp( "hard_elastic_pdf", 
                                        d_hard_elastic_pdf );
*/
  ar & boost::serialization::make_nvp( "analog_elastic_angles", 
                                        d_analog_elastic_angles );
  ar & boost::serialization::make_nvp( "analog_elastic_pdf", 
                                        d_analog_elastic_pdf );
  ar & boost::serialization::make_nvp( 
            "screened_rutherford_normalization_constant", 
            d_screened_rutherford_normalization_constant );
  ar & boost::serialization::make_nvp( "moliere_screening_constant", 
                                        d_moliere_screening_constant );
  ar & boost::serialization::make_nvp( 
            "moment_preserving_elastic_discrete_angles", 
            d_moment_preserving_elastic_discrete_angles );
  ar & boost::serialization::make_nvp( "moment_preserving_elastic_weights", 
                                        d_moment_preserving_elastic_weights );
  ar & boost::serialization::make_nvp( "electroionization_energy_grid",
				                        d_electroionization_energy_grid );
  ar & boost::serialization::make_nvp( "electroionization_recoil_energy",
				                        d_electroionization_recoil_energy );
  ar & boost::serialization::make_nvp( "electroionization_recoil_pdf",
				                        d_electroionization_recoil_pdf );
  ar & boost::serialization::make_nvp( "bremsstrahlung_energy_grid",
				                        d_bremsstrahlung_energy_grid );
  ar & boost::serialization::make_nvp( "bremsstrahlung_photon_energy",
				                        d_bremsstrahlung_photon_energy );
  ar & boost::serialization::make_nvp( "bremsstrahlung_photon_pdf",
				                        d_bremsstrahlung_photon_pdf );
  ar & boost::serialization::make_nvp( "atomic_excitation_energy_grid",
				                        d_atomic_excitation_energy_grid );
  ar & boost::serialization::make_nvp( "atomic_excitation_energy_loss",
				                        d_atomic_excitation_energy_loss );
  ar & boost::serialization::make_nvp( "electron_energy_grid",
				                        d_electron_energy_grid );
  ar & boost::serialization::make_nvp( "cutoff_elastic_cross_section", 
                                        d_cutoff_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "cutoff_elastic_cross_section_threshold_index", 
            d_cutoff_elastic_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( 
            "screened_rutherford_elastic_cross_section", 
            d_screened_rutherford_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "screened_rutherford_elastic_cross_section_threshold_index", 
            d_screened_rutherford_elastic_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( "total_elastic_cross_section", 
                                        d_total_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "total_elastic_cross_section_threshold_index", 
            d_total_elastic_cross_section_threshold_index );
/*
  ar & boost::serialization::make_nvp( "hard_elastic_cross_section", 
                                        d_hard_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "hard_elastic_cross_section_threshold_index", 
            d_hard_elastic_cross_section_threshold_index );
*/
  ar & boost::serialization::make_nvp( 
            "moment_preserving_elastic_cross_section", 
            d_moment_preserving_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "moment_preserving_elastic_cross_section_threshold_index", 
            d_moment_preserving_elastic_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( 
            "electroionization_subshell_cross_section", 
            d_electroionization_subshell_cross_section );
  ar & boost::serialization::make_nvp( 
            "electroionization_subshell_cross_section_threshold_index", 
            d_electroionization_subshell_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( "bremsstrahlung_cross_section", 
                                        d_bremsstrahlung_cross_section );
  ar & boost::serialization::make_nvp( 
            "bremsstrahlung_cross_section_threshold_index", 
            d_bremsstrahlung_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( "atomic_excitation_cross_section", 
                                        d_atomic_excitation_cross_section );
  ar & boost::serialization::make_nvp( 
            "atomic_excitation_cross_section_threshold_index", 
            d_atomic_excitation_cross_section_threshold_index );
}
  
// Load the data from an archive
template<typename Archive>
void EvaluatedElectronDataContainer::load( Archive& ar, 
						  const unsigned version )
{
  ar & boost::serialization::make_nvp( "atomic_number", d_atomic_number );
  ar & boost::serialization::make_nvp( "subshells", d_subshells );
  ar & boost::serialization::make_nvp( "subshell_occupancies", 
				       d_subshell_occupancies );
  ar & boost::serialization::make_nvp( "subshell_binding_energies", 
				       d_subshell_binding_energies );
  ar & boost::serialization::make_nvp( "relaxation_transitions", 
				       d_relaxation_transitions );
  ar & boost::serialization::make_nvp( "relaxation_vacancies", 
				       d_relaxation_vacancies );
  ar & boost::serialization::make_nvp( "relaxation_particle_energies", 
				       d_relaxation_particle_energies );
  ar & boost::serialization::make_nvp( "relaxation_probabilities", 
				       d_relaxation_probabilities );
  ar & boost::serialization::make_nvp( "cutoff_angle", 
                                        d_cutoff_angle );
  ar & boost::serialization::make_nvp( "angular_energy_grid",
				                        d_angular_energy_grid );
/*
  ar & boost::serialization::make_nvp( "number_of_discrete_angles", 
                                        d_number_of_discrete_angles );
  ar & boost::serialization::make_nvp( "hard_elastic_angles", 
                                        d_hard_elastic_angles );
  ar & boost::serialization::make_nvp( "hard_elastic_pdf", 
                                        d_hard_elastic_pdf );
*/
  ar & boost::serialization::make_nvp( "analog_elastic_angles", 
                                        d_analog_elastic_angles );
  ar & boost::serialization::make_nvp( "analog_elastic_pdf", 
                                        d_analog_elastic_pdf );
  ar & boost::serialization::make_nvp( 
            "screened_rutherford_normalization_constant", 
            d_screened_rutherford_normalization_constant );
  ar & boost::serialization::make_nvp( "moliere_screening_constant", 
                                        d_moliere_screening_constant );
  ar & boost::serialization::make_nvp( 
            "moment_preserving_elastic_discrete_angles", 
            d_moment_preserving_elastic_discrete_angles );
  ar & boost::serialization::make_nvp( "moment_preserving_elastic_weights", 
                                        d_moment_preserving_elastic_weights );
  ar & boost::serialization::make_nvp( "electroionization_energy_grid",
				                        d_electroionization_energy_grid );
  ar & boost::serialization::make_nvp( "electroionization_recoil_energy",
				                        d_electroionization_recoil_energy );
  ar & boost::serialization::make_nvp( "electroionization_recoil_pdf",
				                        d_electroionization_recoil_pdf );
  ar & boost::serialization::make_nvp( "bremsstrahlung_energy_grid",
				                        d_bremsstrahlung_energy_grid );
  ar & boost::serialization::make_nvp( "bremsstrahlung_photon_energy",
				                        d_bremsstrahlung_photon_energy );
  ar & boost::serialization::make_nvp( "bremsstrahlung_photon_pdf",
				                        d_bremsstrahlung_photon_pdf );
  ar & boost::serialization::make_nvp( "atomic_excitation_energy_grid",
				                        d_atomic_excitation_energy_grid );
  ar & boost::serialization::make_nvp( "atomic_excitation_energy_loss",
				                        d_atomic_excitation_energy_loss );
  ar & boost::serialization::make_nvp( "electron_energy_grid",
				                        d_electron_energy_grid );
  ar & boost::serialization::make_nvp( "cutoff_elastic_cross_section", 
                                        d_cutoff_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "cutoff_elastic_cross_section_threshold_index", 
            d_cutoff_elastic_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( 
            "screened_rutherford_elastic_cross_section", 
            d_screened_rutherford_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "screened_rutherford_elastic_cross_section_threshold_index", 
            d_screened_rutherford_elastic_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( "total_elastic_cross_section", 
                                        d_total_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "total_elastic_cross_section_threshold_index", 
            d_total_elastic_cross_section_threshold_index );
/*
  ar & boost::serialization::make_nvp( "hard_elastic_cross_section", 
                                        d_hard_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "hard_elastic_cross_section_threshold_index", 
            d_hard_elastic_cross_section_threshold_index );
*/
  ar & boost::serialization::make_nvp( 
            "moment_preserving_elastic_cross_section", 
            d_moment_preserving_elastic_cross_section );
  ar & boost::serialization::make_nvp( 
            "moment_preserving_elastic_cross_section_threshold_index", 
            d_moment_preserving_elastic_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( 
            "electroionization_subshell_cross_section", 
            d_electroionization_subshell_cross_section );
  ar & boost::serialization::make_nvp( 
            "electroionization_subshell_cross_section_threshold_index", 
            d_electroionization_subshell_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( "bremsstrahlung_cross_section", 
                                        d_bremsstrahlung_cross_section );
  ar & boost::serialization::make_nvp( 
            "bremsstrahlung_cross_section_threshold_index", 
            d_bremsstrahlung_cross_section_threshold_index );
  ar & boost::serialization::make_nvp( "atomic_excitation_cross_section", 
                                        d_atomic_excitation_cross_section );
  ar & boost::serialization::make_nvp( 
            "atomic_excitation_cross_section_threshold_index", 
            d_atomic_excitation_cross_section_threshold_index );
}

} // end Data namespace

#endif // end DATA_EVALUATED_ELECTRON_DATA_CONTAINER_DEF_HPP

//---------------------------------------------------------------------------//
// end Data_EvaluatedElectronDataContainer_def.hpp
//---------------------------------------------------------------------------//
