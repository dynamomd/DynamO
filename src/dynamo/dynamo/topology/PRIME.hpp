/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <dynamo/topology/topology.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <vector>

namespace dynamo {
  /*! \brief A Topology class describing the PRIME potential.

    This class also exposes members containing all interaction
    parameters for the PRIME potentials in DynamO. 
    
    A number of parameters are not available in the original
    publications. These include some of the masses of sidechain sites.

    Sources of data are:
    [1] "α-Helix formation: Discontinuous molecular dynamics on an
    intermediate-resolution protein model", Smith and Hall (2001)
    http://onlinelibrary.wiley.com/doi/10.1002/prot.1100/full
    
    [2] "Solvent effects on the conformational transition of a model
    polyalanine peptide", Nguyen, Marchut and Hall (2004)
    http://onlinelibrary.wiley.com/doi/10.1110/ps.04701304/full
    
    [3] "Spontaneous Formation of Twisted Aβ16-22 Fibrils in Large-Scale
    Molecular-Dynamics Simulations", Cheon, Chang and Hall (2011)
    http://www.cell.com/biophysj/fulltext/S0006-3495%2811%2901018-6
    
    [4] "Influence of temperature on formation of perfect tau fragment
    fibrils using PRIME20/DMD simulations", Cheon, Chang and Hall (2012)
    http://onlinelibrary.wiley.com/doi/10.1002/pro.2141/full
    
    [5] "Extending the PRIME model for protein aggregation to all 20
    amino acids", Cheon, Chang and Hall (2010)
    http://onlinelibrary.wiley.com/doi/10.1002/prot.22817/full
    
    [6] "Assembly of a tetrameric α-helical bundle: Computer simulations
    on an intermediate-resolution protein model", Smith and Hall (2001)
    http://onlinelibrary.wiley.com/doi/10.1002/prot.1103/abstract

    All values listed after _PRIME_bond_tolerance are taken from [1]
    unless otherwise noted.
  */
  class TPRIME: public Topology
  {
  public:
    
    /*! \brief  An enumeration of the bead/site types in the PRIME potential. */
    enum PRIME_site_type {
      NH, CH, CO, A, C, D, E, F, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
      GROUP_COUNT
    };
    
    /*! \brief A mapping of \ref PRIME_site_type enumerations to
        string representations for output. 
    */
    static const std::vector<std::string> PRIME_site_names;

    /*! \brief Masses of each site in the PRIME potential.

      Sourced from: [3] for CH, NH, CO, K, L, V, F, A and E. [4] for
      Q, I and Y. Other values calculated from molecular weights.
     */
    static const double _PRIME_masses[22];
    
    /*! \brief Unbonded interaction well-depths for the PRIME potential.

      Sourced from [5]. Set to zero if the interaction is a
      hard-sphere interaction between sites.
    */
    static const double _PRIME_well_depths[22*22];

    /*! \brief Unbonded interaction well-diameters. 

      Largely sourced from [5].
      
      For SC-BB interactions, bead diameters from _PRIME_diameters multipled by 1.5.
      
      Set to 0.000 if the interaction is hard sphere.
    */
    static const double _PRIME_well_diameters[22*22];

    /*! \brief PRIME bead hard-sphere diameters.

      SC-SC diameters sourced from [5]. BB-BB diameters from [1].
      SC-BB diameters mix SC bead size of 5.0 linearly with the BB-BB
      diameters.  to the alanine residue.
    */
    static const double _PRIME_diameters[22*22];
      
    /*\brief [Pseudo]Bond distances from backbone sites to SC sites.
      
      While the original PRIME authors use differing values for each
      SC site, the values are not publicly available and we find this
      to be a suitable simplification
    */
    static const double _PRIME_SC_BB_bond_lengths[3];

    /*! \brief The allowed fluctuation of the bond distance.
    
      There is an older value of 0.02 given in [1], which is relevant
      for validation purposes; However, the latest value is given in
      [2] as it gives more realistic Ramachandran plots.
    */
    static const double _PRIME_bond_tolerance;

    /*! \brief Bond lengths in the backbone.
    
      These bond lengths only apply to beads seperated by 1 backbone
      bond.
      
      We make this a symmetric tensor to simplify the lookups but the zero
      entries should never be used.
    */
    static const double _PRIME_BB_bond_lengths[9];

    /*! \brief Pesudo-bond lengths in the backbone.  

      These only apply to beads seperated by TWO backbone bonds.  We
      make this a symmetric tensor to simplify the lookups. The zero
      entries should never be used.
    */
    static const double _PRIME_pseudobond_lengths[9];

    /*! \brief Pesudo-bond length for the CH-CH backbone interaction.  
      
      This only applies to CH-CH beads/sites separated by THREE backbone bonds.
      
      This is the special pseudobond length between the CH-CH groups. It
      is the only pseudobond at this distance.
    */
    static const double _PRIME_CH_CH_pseudobond_length;
    
    /*! \brief Scaling of interactions in the backbone at a distance
        of 3 bonds (excl. CH-CH interactions).
    
      This is the scaling factor used on the bead diameters if they
      are closer than four bonds on the same chain.

      This value is taken from the PRIME publication. [1]
    */
    static const double _PRIME_3_bonds_scale_factor;

    /*! \brief Scaling of interactions in the backbone at a distance
        of 4 bonds.
	
	This was found to be necessary for alpha-helix Hbond
	formation. (GIL can you clarify this?)
    */
    static const double _PRIME_4_bonds_scale_factor;

    /*! \brief Backbone-to-backbone hydrogen bonding well diameter.

      This is the maximum distance for the CO and NH sites in a
      Hydrogen bond.  The old value from PRIME [1] of 4.20 has been
      replaced/updated in PRIME20 [4].
    */	
    static const double _PRIME_HB_well_diameter;

    /*! \brief Auxillary pair maximum distance. */	
    static const double _PRIME_HB_aux_min_distances[9];

    /*! \brief An enumeration used to identify where in the chain a
        bead is located. */
    enum BeadLocation { NH_END, MID, CO_END };

    struct BeadData {
      BeadData(PRIME_site_type t, size_t r, BeadLocation e = MID): 
	bead_type(t), residue(r), location(e) {}

      //Define sorting for storing in maps (note location does not affect the sorting order!)
      bool operator<(const BeadData& op) const { return bead_type<op.bead_type || (!(op.bead_type<bead_type) && (residue<op.residue)); }
      bool operator==(const BeadData& op) const { return (residue == op.residue) && (bead_type == op.bead_type); }

      PRIME_site_type bead_type;
      size_t residue;
      BeadLocation location;
    };
    
    typedef boost::bimaps::bimap<boost::bimaps::unordered_set_of<size_t>, boost::bimaps::unordered_set_of<BeadData> > BeadTypeMap;

    TPRIME(const magnet::xml::Node&, dynamo::Simulation*, unsigned int ID);

    virtual ~TPRIME() {}
  
    virtual void operator<<(const magnet::xml::Node&);

    BeadData getBeadInfo(size_t ID) const { 
      const auto it = _types->left.find(ID);
#ifdef DYNAMO_DEBUG
      if (it == _types->left.end())
	M_throw() << "Particle " << ID << " has no bead data for " << getName();
#endif
      return it->second;
    }
    
    size_t getBeadID(BeadData data) const { 
      const auto it = _types->right.find(data);
#ifdef DYNAMO_DEBUG
      if (it == _types->right.end())
	M_throw() << "Particle " << ID << " has no bead data for " << getName();
#endif
      return it->second;
    }

  protected:
    std::shared_ptr<BeadTypeMap> _types;
    std::vector<std::pair<size_t, std::string> > _configData;
    virtual void outputXML(magnet::xml::XmlStream&) const;
  };

  inline size_t hash_value(const dynamo::TPRIME::BeadData& x) {
    std::size_t hash_value = 0;
    boost::hash_combine(hash_value, x.residue);
    boost::hash_combine(hash_value, x.bead_type);
    return  hash_value;
  }
}

namespace std {
  template <> struct hash<dynamo::TPRIME::BeadData>
  {
    size_t operator()(const dynamo::TPRIME::BeadData& x) const {
      return dynamo::hash_value(x);
    }
  };
}
