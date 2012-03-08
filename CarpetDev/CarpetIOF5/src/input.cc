#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <hdf5.h>
#include <F5/F5F.h>
#include <F5/F5R.h>
#include <F5/F5iterate.h>
#include <F5/F5uniform.h>

#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"

#include <bbox.hh>
#include <vect.hh>

#include <carpet.hh>

#include "iof5.hh"



// F5 helper
namespace {
  char const*
  ArrayTypeName(ArrayType const array_type)
  {
    switch (array_type) {
    case UnknownArrayType:            return "Unknown";
    case Contiguous:                  return "Contiguous";
    case SeparatedCompound:           return "SeparatedCompound";
    case Constant:                    return "Constant";
    case FragmentedContiguous:        return "FragmentedContiguous";
    case FragmentedSeparatedCompound: return "FragmentedSeparatedCompound";
    case DirectProduct:               return "DirectProduct";
    case IndexPermutation:            return "IndexPermutation";
    case UniformSampling:             return "UniformSampling";
    case FragmentedUniformSampling:   return "FragmentedUniformSampling";
    default:                          return "(illegal)";
    }
  }
}



namespace CarpetIOF5 {
  
  using namespace std;
  using namespace Carpet;
  
  
  
  // Use a class for reading data, so that we have an easy way to pass
  // variables between the various iterators
  class input_iterator_t {
    cGH const *cctkGH;
    
    vector<bool> const input_var; // whether to input this variable
    
    double time;
    char const *gridname;
    char const *topologyname;
    int index_depth;            // 0=vertex, 1=cell
    int topological_dimension;
    int reflevel;
    char const *fieldname;
    ArrayType fieldtype;   // determines how to interpret the fragment
    int varindex;
    char const *fragmentname;
    int map, component;
    
    // string chartname;
    
  public:
    input_iterator_t (cGH const *const cctkGH_,
                      vector<bool> const& input_var_)
      : cctkGH(cctkGH_),
        input_var(input_var_),
        time(nan),
        gridname(NULL),
        topologyname(NULL), index_depth(-1), topological_dimension(-1),
        reflevel(-1),
        fieldname(NULL), fieldtype(UnknownArrayType),
        varindex(-1),
        fragmentname(NULL),
        map(-1), component(-1)
    {
    }
    
    
    
  private:
    
    void read_timeslice (F5Path *const path)
    {
      indent_t indent;
      cout << indent << "time=" << time << "\n";
      
#if 0
      // TODO: Loop over all charts that exist for the grid, or at
      // least remember how many maps there are. (This is also written
      // as part of the grid structure.)
      // F5iterate_grid_atlas?
      for (int m=0; m<maps; ++m) {
        indent_t indent;
        chartname = generate_chartname(cctkGH, m);
        cout << indent << "chart=\"" << chartname << "\"\n";
        F5iterate_grids
          (path, NULL, grid_iterator, this, NULL, chartname.c_str());
      }
#endif
      
      F5iterate_grids (path, NULL, grid_iterator, this, NULL, NULL);
      
      // Synchronise
      BEGIN_REFLEVEL_LOOP(cctkGH) {
#if 0
        int groups[] = {
          CCTK_GroupIndex("GRID::COORDINATES")
        };
        int const ngroups = sizeof groups / sizeof *groups;
        for (int i=0; i<ngroups; ++i) {
          assert (groups[i]>=0);
        }
        iret = CCTK_SyncGroupsI (cctkGH, ngroups, groups);
        assert (iret >= 0);
#endif
#warning "TODO: don't modify boundaries"
        BEGIN_LOCAL_MAP_LOOP(cctkGH, CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            DECLARE_CCTK_ARGUMENTS;
#pragma omp parallel
            CCTK_LOOP3_BOUNDARIES(boundaries, cctkGH,
                                  i,j,k, 1,1,1,
                                  cctk_lsh[0]-1,cctk_lsh[1]-1,cctk_lsh[2]-1)
            {
              int const ind3d = CCTK_GFINDEX3D(cctkGH, i,j,k);
              x[ind3d] = 0.0;
              y[ind3d] = 0.0;
              z[ind3d] = 0.0;
              r[ind3d] = 0.0;
            } CCTK_ENDLOOP3_BOUNDARIES(boundaries);
          } END_LOCAL_COMPONENT_LOOP;
        } END_LOCAL_MAP_LOOP;
      } END_REFLEVEL_LOOP;
    }
  
    void read_grid (F5Path *const path)
    {
      indent_t indent;
      cout << indent << "grid=\"" << gridname << "\"\n";
      // F5iterate_vertex_fields (path, NULL, field_iterator, this, NULL, NULL);
      F5iterate_topologies (path, NULL, topology_iterator, this);
    }
    
    void read_topology (F5Path *const path)
    {
      indent_t indent;
      
      herr_t herr;
      
      cout << indent
           << "topology=\"" << topologyname << "\""
           << " (" << (index_depth==0 ? "vertex" : "cell") << ")\n"
           << indent
           << "topological dimension=" << topological_dimension << "\n";
      
      
      
      // Ignore topologies that are only an alias for another topology
      H5G_stat_t stat;
      herr = H5Gget_objinfo (path->Grid_hid, topologyname, false, &stat);
      assert (not herr);
      if (stat.type == H5G_LINK) {
        char linkval[100000];
        herr = H5Gget_linkval
          (path->Grid_hid, topologyname, sizeof linkval, linkval);
        assert (not herr);
        indent_t indent;
        cout << indent << "alias for topology \"" << linkval << "\"\n"
             << indent << "ignoring this topology\n";
        return;
      }
      
      
      
      // Determine refinement level from the topology
      // (Could also use topology name instead)
      hsize_t hreffact[FIBER_MAX_RANK];
      int const iret =
        F5LAget_dimensions(path->Topology_hid,
                           FIBER_HDF5_REFINEMENT_INFO, hreffact);
      assert (iret == dim);
      hsize_t hreffact2[FIBER_MAX_RANK];
      void *const pret = F5Tpermute_dimensions(path, dim, hreffact2, hreffact);
      assert (pret);
      ivect const reffact = h2v(hreffact2);
      int rl;
      for (rl=0; rl<reflevels; ++rl) {
        if (all(reffact == Carpet::spacereffacts.AT(rl))) break;
      }
      assert (rl<reflevels);
      reflevel = rl;
      cout << indent << "refinement level " << reflevel << "\n";
      
      
      
      // F5iterate_topology_fields
      //   (path, NULL, field_iterator, this, chartname.c_str(), NULL); 
      F5iterate_topology_fields (path, NULL, field_iterator, this, NULL, NULL); 
    }
    
    void read_field (F5Path *const path)
    {
      indent_t indent;
      cout << indent << "field=\"" << fieldname << "\"\n";
      
      interpret_fieldname(cctkGH, fieldname, varindex);
#warning "TODO: check all variables in the group"
#warning "TODO: loop over all variables in the group"
      if (varindex>=0 and input_var.at(varindex)) {

        int major_version, minor_version, release_version;
        fieldtype =
          F5Fget_field_enum(path, fieldname,
                            &major_version, &minor_version, &release_version);
        cout << indent << "field_type=" << ArrayTypeName(fieldtype) << "\n";
        
        int const iret = F5Fopen(path, fieldname);
        assert(iret);
        
        // Do we need to iterate over fragments?
        int const is_fragmented = F5Fis_fragmented(path, fieldname);
        cout << indent
             << (is_fragmented ? "fragmented" : "not fragmented") << "\n";
        if (is_fragmented) {
          F5iterate_field_fragments (path, NULL, fragment_iterator, this);
        } else {
          read_fragment (path);
        }
        
        F5Fclose(path);
        
      } else {
        indent_t indent;
        cout << indent << "ignoring this field\n";
      }
    }
    
    void read_fragment (F5Path *const path)
    {
      indent_t indent;
      
      int iret;
      void *pret;
      
      if (fragmentname) {
        cout << indent << "fragment=\"" << fragmentname << "\"\n";
      } else {
        cout << indent << "no fragment\n";
      }
      
      int major_version, minor_version, release_version;
      fieldtype =
        F5Fget_field_enum(path, fieldname,
                          &major_version, &minor_version, &release_version);
      cout << indent << "field_type=" << ArrayTypeName(fieldtype) << "\n";
      
      
      
      // Determine map and component from fragment name
      if (fragmentname) {
        interpret_fragmentname (cctkGH, fragmentname, map, component);
      } else {
        map = 0;
        component = CCTK_MyProc(cctkGH);
      }
      cout << indent << "map " << map << " component " << component << "\n";
      
      // This routine has a bug, and does not recognise
      // "SeparatedCompound" as separated
      //int const is_separated = F5Fis_separatedcompound(path, fieldname);
      int const is_separated =
        fieldtype == FragmentedSeparatedCompound or
        fieldtype == SeparatedCompound;
      cout << indent
           << (is_separated ? "separated" : "contiguous") << "\n";
      hid_t const type_id = F5Fget_type(path);
      if (F5Tis_convertible(type_id, H5T_NATIVE_DOUBLE)) {
        cout << indent << "compound type: scalar\n";
        read_variable (path, NULL, varindex);
      } else if (F5Tis_convertible(type_id, F5T_VEC3_DOUBLE)) {
        cout << indent << "compound type: vector\n";
        // This assumes a separated vector; we don't support
        // contiguous vectors
        assert(is_separated);
        read_variable (path, "Dx", varindex+0);
        read_variable (path, "Dy", varindex+1);
        read_variable (path, "Dz", varindex+2);
      } else if (F5Tis_convertible(type_id, F5T_METRIC33_DOUBLE)) {
        cout << indent << "compound type: symmetric tensor\n";
        // Not yet implemented
        assert(0);
      } else {
        // Unknown tensor type
        assert(0);
      }
    }
    
    void read_variable (F5Path *const path, char const *const name,
                        int const var)
    {
      indent_t indent;
      
      herr_t herr;
      int ierr;
      int iret;
      void *pret;
      
      assert (path);
      assert (var >= 0);
      
      cout << indent << "dataset=\"" << (name ? name : "(null)") << "\"\n";
      
      assert (var >= 0);
      {
        char *const fullname = CCTK_FullName(var);
        cout << indent << "variable=\"" << fullname << "\"\n";
        free (fullname);
      }
      
      
      
      // Determine fragment properties
      
      // A dataset can either be the same as a field, or can be in a
      // group containing all fragments, and/or can be in a group
      // containing all compound (vector) elements.
      // The possible group hierarchies are (the rightmost element is
      // the dataset):
      //    .../FIELD
      //    .../FIELD/ELEMENT
      //    .../FIELD/FRAGMENT
      //    .../FIELD/FRAGMENT/ELEMENT
      // The field is the Cactus group/variable name, a fragment
      // describes map and component, the element describes the tensor
      // element.
      // We describe this structure with two booleans:
      //    bool fragment_is_group: the fragment is a group
      //    bool field_is_dataset:  the field is a dataset
      // Note that the F5 library will already have openend the field,
      // either as group or as dataset.
      
      // There is a group for all fragments if (a) we expect a
      // fragment, and (b) there exists a group with this name. (We
      // probably could examine F5's internal state as well, but
      // looking for an HDF5 group is the most robust approach.)
      bool fragment_is_group = false;
      if (fragmentname) {
        H5O_info_t info;
        herr = H5Oget_info_by_name (path->Field_hid, fragmentname,
                                    &info, H5P_DEFAULT);
        assert (not herr);
        fragment_is_group = info.type == H5O_TYPE_GROUP;
      }
      cout << indent << "fragment_is_group=" << fragment_is_group << "\n";
      hid_t fragment;
      if (fragment_is_group) {
        fragment = H5Gopen (path->Field_hid, fragmentname, H5P_DEFAULT);
        assert (fragment >= 0);
      } else {
        fragment = path->Field_hid;
      }
      
      // The field consists of a single dataset if we cannot open the
      // element
      hid_t element = -1;
      char const *datasetname = NULL;
      if (name) {
        datasetname = name;
      } else if (fragmentname and not fragment_is_group) {
        datasetname = fragmentname;
      }
      if (datasetname) {
        H5E_BEGIN_TRY {
          element = H5Dopen(fragment, datasetname, H5P_DEFAULT);
        } H5E_END_TRY;
        assert (element >= 0);
      }
      bool const field_is_dataset = element < 0;
      cout << indent << "field_is_dataset=" << field_is_dataset << "\n";
      if (field_is_dataset) {
        // Can't open the fragment -- we assume the field consists of
        // a single element and has already been opened
        element = fragment;
      }
      assert (element>=0);
      
      // Check index depth
      int index_depth_;
      iret = F5Tget_index_depth(path, &index_depth_);
      assert (iret);
      assert (index_depth_ == index_depth);
      
      // Read the fragment offset. This is stored with the dataset
      // group.
      ivect foff = 0;
      if (fragmentname) {
        hsize_t hoff[FIBER_MAX_RANK];
        iret = F5LAget_dimensions(fragment_is_group ? fragment : element,
                                  FIBER_FRAGMENT_OFFSET_ATTRIBUTE, hoff);
        assert (iret == dim);
        hsize_t hoff2[FIBER_MAX_RANK];
        pret = F5Tpermute_dimensions(path, dim, hoff2, hoff);
        assert (pret);
        foff = h2v(hoff2);
      }
      assert (all(foff>=0));
      
#if 0
      // Read the fragment size. This is stored with the field -- why
      // is this different from the offset?
      hsize_t hlen[FIBER_MAX_RANK];
      iret =
        F5LAget_dimensions(path->Field_hid,
                           FIBER_FIELD_DATASPACE_DIMENSIONS_ATTRIBUTE, hlen);
      assert (iret == dim);
      hsize_t hlen2[FIBER_MAX_RANK];
      pret = F5Tpermute_dimensions(path, dim, hlen2, hlen);
      assert (pret);
      ivect const flen = h2v(hlen2);
      assert (all(flen>=0));
#endif
      hid_t const space = H5Dget_space(element);
      assert (space>=0);
      iret = H5Sget_simple_extent_ndims(space);
      assert (iret == dim);
      hsize_t hlen[dim];
      iret = H5Sget_simple_extent_dims(space, hlen, NULL);
      hsize_t hlen2[dim];
      pret = F5Tpermute_dimensions(path, dim, hlen2, hlen);
      assert (pret);
      ivect const flen = h2v(hlen2);
      assert (all(flen>=0));
      herr = H5Sclose (space);
      assert (not herr);
      
      ibbox const fbox (foff, foff+flen-1, 1);
      {
        indent_t indent;
        cout << indent
             << "dataset bbox is " << foff << ":" << foff+flen-1 << "\n";
      }
      
      
      
      // Examine Cactus variable
      int const group = CCTK_GroupIndexFromVarI(var);
      cGroup groupdata;
      ierr = CCTK_GroupData(group, &groupdata);
      assert (not ierr);
      
      // TODO
      assert (groupdata.grouptype == CCTK_GF);
      assert (groupdata.vartype == CCTK_VARIABLE_REAL);
      assert (groupdata.disttype == CCTK_DISTRIB_DEFAULT);
      assert (groupdata.stagtype == 0);
      assert (groupdata.dim == dim);
      
      // TODO: This does not work; we need to read the data, and later
      // distribute it onto the respective Cactus variables
      ENTER_LEVEL_MODE(cctkGH, reflevel) {
        ENTER_SINGLEMAP_MODE(cctkGH, map, CCTK_GF) {
          ENTER_LOCAL_MODE(cctkGH, component, CCTK_GF) {
            DECLARE_CCTK_ARGUMENTS;
            
            ivect lbnd, lsh, lghosts, ughosts;
            ivect imin, imax, ioff, ilen;
            for (int d=0; d<dim; ++d) {
              lbnd[d]    = cctk_lbnd[d];
              lsh[d]     = cctk_lsh[d];
              // F5 counts the total overlap, which is the sum of the
              // ghost zones on this and on the adjacent component
              lghosts[d] = cctk_bbox[2*d  ] ? 0 : 2*cctk_nghostzones[d];
              ughosts[d] = cctk_bbox[2*d+1] ? 0 : 2*cctk_nghostzones[d];
              imin[d] = 0;
              imax[d] = lsh[d];
              // Do not read in ghost zones
              imin[d] += lghosts[d] / 2;
              imax[d] -= ughosts[d] / 2;
              ioff[d] = lbnd[d] + imin[d];
              ilen[d] = imax[d] - imin[d];
            }
            ibbox const lbox (lbnd, lbnd+lsh-1, 1);
            ibbox const ibox (ioff, ioff+ilen-1, 1);
            ibbox const ovlp = ibox & fbox;
            
            if (ovlp.empty()) {
              indent_t indent;
              cout << indent << "dataset does not intersect component " << component << "; skipping component\n";
              continue;
            }
            
            indent_t indent;
            cout << indent << "dataset intersects component " << component << "\n"
                 << indent << "reading bbox " << ovlp.lower() << ":" << ovlp.upper() << " of " << ibox.lower() << ":" << ibox.upper() << "\n";
            
            #warning "TODO: assert the whole component is read"
            
            int const proc = vhh.AT(Carpet::map)->processor(reflevel, component);
#warning "TODO: send dataset to destination process"
            assert (proc == dist::rank());
            
            int const timelevel = 0;
            void *const data = CCTK_VarDataPtrI(cctkGH, timelevel, var);
            assert (data);
            
            
            
            hvect const hzero(0);
            hvect mlsh;
            pret = F5Tpermute_dimensions
              (path, dim, &mlsh[0], &v2h(lbox.shape())[0]);
            assert (pret);
            hid_t const memspace = H5Screate_simple(dim, &mlsh[0], NULL);
            assert (memspace >= 0);
            hvect mmin;
            pret = F5Tpermute_dimensions
              (path, dim, &mmin[0], &v2h(ovlp.lower()-lbox.lower())[0]);
            assert (pret);
            hvect mlen;
            pret = F5Tpermute_dimensions
              (path, dim, &mlen[0], &v2h(ovlp.shape())[0]);
            assert (pret);
            assert (all(mmin >= hzero));
            assert (all(mmin+mlen <= mlsh));
            herr = H5Sselect_hyperslab
              (memspace, H5S_SELECT_SET, &mmin[0], NULL, &mlen[0], NULL);
            assert (not herr);
            // cout << "mlsh=" << mlsh << " mmin=" << mmin << " mlen=" << mlen << "\n";
            
            hvect flsh;
            pret = F5Tpermute_dimensions
              (path, dim, &flsh[0], &v2h(fbox.shape())[0]);
            assert (pret);
            hid_t const filespace = H5Screate_simple(dim, &flsh[0], NULL);
            assert (filespace >= 0);
            hvect fmin;
            pret = F5Tpermute_dimensions
              (path, dim, &fmin[0], &v2h(ovlp.lower()-fbox.lower())[0]);
            assert (pret);
            hvect const flen = mlen;
            assert (all(fmin >= hzero));
            assert (all(fmin+flen <= flsh));
            herr = H5Sselect_hyperslab
              (filespace, H5S_SELECT_SET, &fmin[0], NULL, &flen[0], NULL);
            // cout << "flsh=" << flsh << " fmin=" << fmin << " flen=" << flen << "\n";
            
            herr = H5Dread (element, H5T_NATIVE_DOUBLE, memspace, filespace,
                            H5P_DEFAULT, data);
            assert (not herr);
            
            herr = H5Sclose(memspace);
            assert (not herr);
            herr = H5Sclose(filespace);
            assert (not herr);
            
          } LEAVE_LOCAL_MODE;
        } LEAVE_SINGLEMAP_MODE;
      } LEAVE_LEVEL_MODE;
      
      
      
      if (not field_is_dataset) {
        herr = H5Dclose(element);
        assert (not herr);
      }
      
      if (fragment_is_group) {
        herr = H5Gclose(fragment);
        assert (not herr);
      }
    }
    
    
    
  public:
    
    void iterate (hid_t const object)
    {
      F5iterate_timeslices (object, NULL, timeslice_iterator, this);
    }
    
    static
    herr_t timeslice_iterator (F5Path *const path, double const time,
                               void *const userdata)
    {
      input_iterator_t *const iterator = (input_iterator_t*)userdata;
      iterator->time = time;
      iterator->read_timeslice (path);
      return 0;
      // finish (successfully) after the first timeslice
      // return 1;
    }
    
    static
    herr_t grid_iterator (F5Path *const path, char const *const gridname,
                          void *const userdata)
    {
      input_iterator_t *const iterator = (input_iterator_t*)userdata;
      iterator->gridname = gridname;
      iterator->read_grid (path);
      return 0;
    }
    
    static
    herr_t topology_iterator (F5Path *const path,
                              char const *const topologyname,
                              int const index_depth,
                              int const topological_dimension,
                              void *const userdata)
    {
      input_iterator_t *const iterator = (input_iterator_t*)userdata;
      iterator->topologyname = topologyname;
      iterator->index_depth = index_depth;
      iterator->topological_dimension = topological_dimension;
      iterator->read_topology (path);
      iterator->topologyname = NULL;
      iterator->index_depth = -1;
      iterator->topological_dimension = -1;
      return 0;
    }
    
    static
    herr_t field_iterator (F5Path *const path, char const *const fieldname,
                           void *const userdata)
    {
      input_iterator_t *const iterator = (input_iterator_t*)userdata;
      iterator->fieldname = fieldname;
      iterator->read_field (path);
      iterator->fieldname = NULL;
      return 0;
    }
    
    static
    herr_t fragment_iterator (F5Path *const path,
                              char const *const fragmentname,
                              void *const userdata)
    {
      input_iterator_t *const iterator = (input_iterator_t*)userdata;
      iterator->fragmentname = fragmentname;
      iterator->read_fragment (path);
      iterator->fragmentname = NULL;
      return 0;
    }
    
  };                            // class input_iterator_t
  
  
  
  void input (cGH const *const cctkGH,
               hid_t const file,
               vector<bool> const& input_var)
  {
    input_iterator_t iterator(cctkGH, input_var);
    iterator.iterate (file);
  }
  
} // end namespace CarpetIOF5
