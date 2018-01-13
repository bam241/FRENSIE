//---------------------------------------------------------------------------//
//!
//! \file   Geometry_Model.hpp
//! \author Alex Robinson
//! \brief  The geometry model base class declaration
//!
//---------------------------------------------------------------------------//

#ifndef GEOMETRY_MODEL_HPP
#define GEOMETRY_MODEL_HPP

// Std Lib Includes
#include <memory>
#include <iostream>

// Boost Includes
#include <boost/serialization/split_member.hpp>
#include <boost/units/systems/cgs/length.hpp>
#include <boost/units/systems/cgs/volume.hpp>

// FRENSIE Includes
#include "Geometry_Navigator.hpp"
#include "Geometry_EstimatorType.hpp"
#include "Geometry_ParticleType.hpp"
#include "Utility_Set.hpp"
#include "Utility_Map.hpp"
#include "Utility_Vector.hpp"
#include "Utility_Tuple.hpp"
#include "Utility_UnitTraits.hpp"
#include "Utility_QuantityTraits.hpp"
#include "Utility_SerializationHelpers.hpp"

namespace Geometry{

//! The model base class
class Model
{

public:

  //! The internal cell handle type
  typedef Navigator::InternalCellHandle InternalCellHandle;
  
  //! The internal surface handle type
  typedef Navigator::InternalSurfaceHandle InternalSurfaceHandle;
  
  //! The internal material handle type
  typedef unsigned long long InternalMaterialHandle;
  
  //! The internal estimator handle type
  typedef unsigned long long InternalEstimatorHandle;

  //! The length unit
  typedef Navigator::LengthUnit LengthUnit;

  //! The length quantity
  typedef Navigator::Length Length;

  //! The volume unit
  typedef boost::units::cgs::volume VolumeUnit;

  //! The volume quantity
  typedef boost::units::quantity<VolumeUnit> Volume;

  /*! The density unit
   *
   * This unit ignores the "mass" component since a density in the model
   * can be defined in terms of mass (negative value) or atoms. The material
   * definition is needed to convert between the two types of densities.
   */
  typedef typename Utility::UnitTraits<VolumeUnit>::InverseUnit DensityUnit;

  //! The density quantity
  typedef boost::units::quantity<DensityUnit> Density;

  //! The material id set type
  typedef std::set<InternalMaterialHandle> MaterialIdSet;

  //! The cell id set type
  typedef Navigator::CellIdSet CellIdSet;

  //! The cell id material id map type
  typedef std::map<InternalCellHandle,InternalMaterialHandle> CellIdMatIdMap;

  //! The cell id density map type
  typedef std::map<InternalCellHandle,Density> CellIdDensityMap;

  //! The cell id array type
  typedef std::vector<InternalCellHandle> CellIdArray;

  //! The cell estimator data type
  typedef std::tuple<EstimatorType,ParticleType,CellIdArray> CellEstimatorData;

  //! The cell estimator id data map type
  typedef std::map<InternalEstimatorHandle,CellEstimatorData> CellEstimatorIdDataMap;

  //! Constructor
  Model()
  { /* ... */ }

  //! Destructor
  virtual ~Model()
  { /* ... */ }

  //! Get the model name
  virtual std::string getName() const = 0;

  //! Check if this is an advanced model
  virtual bool isAdvanced() const;

  //! Check if the model has cell estimator data
  virtual bool hasCellEstimatorData() const = 0;

  //! Get the material ids
  virtual void getMaterialIds( MaterialIdSet& material_ids ) const = 0;

  //! Get the cells
  virtual void getCells( CellIdSet& cells,
                         const bool include_void_cells,
                         const bool include_termination_cells ) const = 0;

  //! Get the cell material ids
  virtual void getCellMaterialIds(
                                CellIdMatIdMap& cell_id_mat_id_map ) const = 0;

  //! Get the cell densities
  virtual void getCellDensities( CellIdDensityMap& cell_density_map ) const = 0;

  //! Get the cell estimator data
  virtual void getCellEstimatorData(
                CellEstimatorIdDataMap& cell_estimator_id_data_map ) const = 0;

  //! Check if a cell exists
  virtual bool doesCellExist( const InternalCellHandle cell ) const = 0;

  //! Check if the cell is a termination cell
  virtual bool isTerminationCell( const InternalCellHandle cell ) const = 0;

  //! Check if a cell is void
  virtual bool isVoidCell( const InternalCellHandle cell ) const = 0;

  //! Get the cell volume
  virtual Volume getCellVolume( const InternalCellHandle cell ) const = 0;

  //! The invalid cell handle
  static InternalCellHandle invalidCellHandle();

  //! The invalid surface handle
  static InternalSurfaceHandle invalidSurfaceHandle();

  //! The invalid material handle
  static InternalMaterialHandle invalidMaterialHandle();

  //! The invalid estimator handle
  static InternalEstimatorHandle invalidEstimatorHandle();

  //! Create a raw, heap-allocated navigator
  virtual Geometry::Navigator* createNavigatorAdvanced() const = 0;
  
  //! Create a navigator
  std::shared_ptr<Geometry::Navigator> createNavigator() const;

private:

  // Save the model to an archive
  template<typename Archive>
  void save( Archive& ar, const unsigned version ) const
  { /* ... */ }

  // Load the model from an archive
  template<typename Archive>
  void load( Archive& ar, const unsigned version )
  { /* ... */ }

  BOOST_SERIALIZATION_SPLIT_MEMBER();

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;
};

// Check if this is an advanced model
inline bool Model::isAdvanced() const
{
  return false;
}

// Create a navigator
inline std::shared_ptr<Geometry::Navigator> Model::createNavigator() const
{
  return std::shared_ptr<Geometry::Navigator>(
                                             this->createNavigatorAdvanced() );
}
  
} // end Geometry namespace

BOOST_SERIALIZATION_ASSUME_ABSTRACT_CLASS( Model, Geometry );
BOOST_SERIALIZATION_CLASS_VERSION( Model, Geometry, 0 );

#endif // end GEOMETRY_MODEL_HPP

//---------------------------------------------------------------------------//
// end Geometry_Model.hpp
//---------------------------------------------------------------------------//
