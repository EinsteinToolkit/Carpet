/*@@
  @file      hdf5_convert_from_sdf.cc
  @date      Tue 24 August 2004
  @author    Thomas Radke
  @date      2021-08-04
  @author    Erik Schnetter
  @desc
             Utility program to convert SDF datafiles into IOHDF5
             datafiles.
  @enddesc
@@*/

// HDF5
#include <hdf5.h>

// SDF
#include <bbhutil.h>
#include <sdf_priv.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  herr_t herr;

  std::cout << "\n"
            << "----------------------------\n"
            << "Carpet SDF-to-HDF5 Converter\n"
            << "----------------------------\n"
            << "\n";

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <varname> <timelevel> <outputfile> <inputfiles...>\n"
              << "   e.g. " << argv[0]
              << " ADMBase::metric 0 gxx.h5 foobar.sdf\n";
    return 1;
  }

  const std::string cactus_varname = argv[1];
  const int timelevel = std::atoi(argv[2]);
  const std::string outputfile = argv[3];
  std::vector<std::string> inputfiles;
  for (int i = 4; i < argc; ++i)
    inputfiles.push_back(argv[i]);

  std::cout << "varname " << cactus_varname << "\n"
            << "timelevel " << timelevel << "\n"
            << "outputfile " << outputfile << "\n"
            << "inputfiles ";
  for (const auto &inputfile : inputfiles)
    std::cout << " " << inputfile;
  std::cout << "\n";

  if (cactus_varname.find("::") == std::string::npos) {
    std::cerr << "varname must contain a \"::\" sequence\n";
    return 1;
  }

  // Global properties

  std::set<double> timesteps;
  double origin[3] = {1.0 / 0.0, 1.0 / 0.0, 1.0 / 0.0};
  double nigiro[3] = {-1.0 / 0.0, -1.0 / 0.0, -1.0 / 0.0};
  double coarse_spacing[3] = {0.0, 0.0, 0.0};
  std::map<std::array<double, 3>, std::set<std::array<int, 6> > > bboxes;

  // First pass: Collect some metadata
  {
    int block = 0;
    for (const auto &inputfile : inputfiles) {
      gft_sdf_file_data *const infile =
          gft_open_sdf_file(const_cast<char *>(inputfile.c_str()));
      if (!infile) {
        std::cerr << "Could not open SDF input file \"" << inputfile << "\"\n";
        return 2;
      }

      double timestep;
      int version;
      int rank;
      int dsize;
      int csize;
      char *varname;
      char *coordname;
      char *tag;
      int *dims;
      double *bbox;
      double *coords;
      double *data;
      while (low_read_sdf_stream(1, infile->fp, &timestep, &version, &rank,
                                 &dsize, &csize, &varname, &coordname, &tag,
                                 &dims, &bbox, &coords, &data)) {

        timesteps.insert(timestep);

        int npoints = 1;
        for (int d = 0; d < rank; ++d)
          npoints *= dims[d];
        assert(npoints == dsize);

        using std::min;
        for (int d = 0; d < rank; ++d)
          origin[d] = min(origin[d], bbox[2 * d]);

        using std::max;
        for (int d = 0; d < rank; ++d)
          nigiro[d] = max(nigiro[d], bbox[2 * d + 1]);

        assert(rank <= 3);
        double spacing[3];
        for (int d = 0; d < rank; ++d)
          spacing[d] = (bbox[2 * d + 1] - bbox[2 * d]) / (dims[d] - 1);

        using std::max;
        for (int d = 0; d < rank; ++d)
          coarse_spacing[d] = max(coarse_spacing[d], spacing[d]);

        std::array<double, 3> spacing_array;
        for (int d = 0; d < 3; ++d)
          spacing_array[d] = d < rank ? spacing[d] : 0.0;
        std::array<int, 6> bbox_array;
        for (int d = 0; d < 6; ++d)
          bbox_array[d] = d < 2 * rank ? bbox[d] : 0;
        bboxes[spacing_array].insert(bbox_array);

        ++block;
      }

      gsfd_close(infile);
    } // for inputfile
  }

  // Find origin_time and delta_time

  const double timestep_epsilon = 1.0e-6;

  // std::cout << "timesteps\n";
  // for (const auto &ts : timesteps)
  //   std::cout << "  " << ts << "\n";

  std::vector<double> timesteps1;
  for (const auto &ts : timesteps)
    if (timesteps1.empty() || ts > timesteps1.back() + timestep_epsilon)
      timesteps1.push_back(ts);
  std::cout << "timesteps\n";
  for (const auto &ts : timesteps1)
    std::cout << "  " << ts << "\n";

  assert(!timesteps1.empty());
  const double origin_time = timesteps1.at(0);
  const double delta_time =
      timesteps1.size() >= 2 ? timesteps1.at(1) - timesteps1.at(0) : 1.0;
  std::cout << "origin_time " << origin_time << "\n"
            << "delta_time " << delta_time << "\n";

  // Check whether the bboxes cover the whole domain

  std::cout << "origin " << origin[0] << " " << origin[1] << " " << origin[2]
            << "\n"
            << "nigiro " << nigiro[0] << " " << nigiro[1] << " " << nigiro[2]
            << "\n";

  int gsh[3] = {1, 1, 1};
  for (int d = 0; d < 3; ++d)
    gsh[d] = std::lrint((nigiro[d] - origin[d]) / coarse_spacing[d]) + 1;
  int di[3];
  di[0] = 1;
  di[1] = di[0] * gsh[0];
  di[2] = di[1] * gsh[1];
  const int np = di[2] * gsh[2];
  std::vector<int> read_count_vec(np, 0);
  const auto read_count = [&](int i, int j, int k) -> int & {
    assert(i >= 0 && i < gsh[0]);
    assert(j >= 0 && j < gsh[1]);
    assert(k >= 0 && k < gsh[2]);
    return read_count_vec.at(i * di[0] + j * di[1] + k * di[2]);
  };

  const hid_t outfile =
      H5Fcreate(outputfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (outfile < 0) {
    std::cerr << "Could not create HDF5 output file \"" << outputfile << "\"\n";
    return 3;
  }

  // add a dummy group so that the HDF5 file is recognized as Carpet format
  const hid_t group = H5Gcreate(outfile, "Parameters and Global Attributes",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  {
    const hid_t scalar_space = H5Screate(H5S_SCALAR);

    {
      const hid_t attr = H5Acreate(group, "nioprocs", H5T_NATIVE_INT,
                                   scalar_space, H5P_DEFAULT, H5P_DEFAULT);
      const int nioprocs = 1;
      herr = H5Awrite(attr, H5T_NATIVE_INT, &nioprocs);
      herr = H5Aclose(attr);
    }

    {
      const hid_t attr = H5Acreate(group, "carpet_reflevels", H5T_NATIVE_INT,
                                   scalar_space, H5P_DEFAULT, H5P_DEFAULT);
      const int carpet_reflevels = bboxes.size();
      herr = H5Awrite(attr, H5T_NATIVE_INT, &carpet_reflevels);
      herr = H5Aclose(attr);
    }

    {
      const hid_t attr = H5Acreate(group, "GH$iteration", H5T_NATIVE_INT,
                                   scalar_space, H5P_DEFAULT, H5P_DEFAULT);
      const int GH_iteration = 0;
      herr = H5Awrite(attr, H5T_NATIVE_INT, &GH_iteration);
      herr = H5Aclose(attr);
    }

    {
      const hid_t attr = H5Acreate(group, "main loop index", H5T_NATIVE_INT,
                                   scalar_space, H5P_DEFAULT, H5P_DEFAULT);
      const int main_loop_index = 0;
      herr = H5Awrite(attr, H5T_NATIVE_INT, &main_loop_index);
      herr = H5Aclose(attr);
    }

    {
      const hid_t attr =
          H5Acreate(group, "carpet_global_time", H5T_NATIVE_DOUBLE,
                    scalar_space, H5P_DEFAULT, H5P_DEFAULT);
      const int carpet_global_time = origin_time;
      herr = H5Awrite(attr, H5T_NATIVE_INT, &carpet_global_time);
      herr = H5Aclose(attr);
    }

    {
      const hid_t attr =
          H5Acreate(group, "carpet_delta_time", H5T_NATIVE_DOUBLE, scalar_space,
                    H5P_DEFAULT, H5P_DEFAULT);
      const int carpet_delta_time = delta_time;
      herr = H5Awrite(attr, H5T_NATIVE_INT, &carpet_delta_time);
      herr = H5Aclose(attr);
    }

    herr = H5Sclose(scalar_space);
  }

  herr = H5Gclose(group);

  std::set<std::string> datasetnames;

  // Second pass: Convert data
  {
    int block = 0;
    for (const auto &inputfile : inputfiles) {
      gft_sdf_file_data *const infile =
          gft_open_sdf_file(const_cast<char *>(inputfile.c_str()));
      if (!infile) {
        std::cerr << "Could not open SDF input file \"" << inputfile << "\"\n";
        return 2;
      }

      double timestep;
      int version;
      int rank;
      int dsize;
      int csize;
      char *varname;
      char *coordname;
      char *tag;
      int *dims;
      double *bbox;
      double *coords;
      double *data;
      while (low_read_sdf_stream(1, infile->fp, &timestep, &version, &rank,
                                 &dsize, &csize, &varname, &coordname, &tag,
                                 &dims, &bbox, &coords, &data)) {

        const int iteration = lrint((timestep - origin_time) / delta_time);

        assert(rank <= 3);
        double spacing[3];
        for (int d = 0; d < rank; ++d)
          spacing[d] = (bbox[2 * d + 1] - bbox[2 * d]) / (dims[d] - 1);

        using std::log2, std::lrint;
        const int reflevel = lrint(log2(coarse_spacing[0] / spacing[0]));
        assert(reflevel >= 0 && reflevel < 100);
        for (int d = 0; d < rank; ++d)
          assert(lrint(log2(coarse_spacing[d] / spacing[d])) == reflevel);

        std::cout << "Processing dataset \"" << varname << "\"\n"
                  << "  timestep " << timestep << "\n"
                  << "  iteration " << iteration << "\n"
                  << "  version " << version << "\n"
                  << "  rank " << rank << "\n"
                  << "  dsize " << dsize << "\n"
                  << "  csize " << csize << "\n"
                  << "  varname \"" << varname << "\"\n"
                  << "  coordname \"" << coordname << "\"\n"
                  << "  tag ";
        if (tag)
          std::cout << "\"" << tag << "\"";
        else
          std::cout << "(null)";
        std::cout << "\n"
                  << "  dims";
        for (int d = 0; d < rank; ++d)
          std::cout << " " << dims[d];
        std::cout << "\n"
                  << "  bbox";
        for (int d = 0; d < 2 * rank; ++d)
          std::cout << " " << bbox[d];
        std::cout << "\n"
                  << "  spacing";
        for (int d = 0; d < rank; ++d)
          std::cout << " " << spacing[d];
        std::cout << "\n"
                  << "  origin";
        for (int d = 0; d < rank; ++d)
          std::cout << " " << origin[d];
        std::cout << "\n"
                  << "  coarse_spacing";
        for (int d = 0; d < rank; ++d)
          std::cout << " " << coarse_spacing[d];
        std::cout << "\n"
                  << "  reflevel " << reflevel << "\n"
                  << std::flush;

        hsize_t h5dims[3];
        for (int d = 0; d < rank; ++d)
          h5dims[rank - 1 - d] = dims[d];
        const hid_t block_space = H5Screate_simple(rank, h5dims, nullptr);

        const std::string varname0 = varname;
        const std::string varname1 = varname0.substr(0, varname0.rfind("_"));
        const std::string component_string =
            varname0.substr(varname0.rfind("_") + 1);
        const int component = std::atoi(component_string.c_str());

        std::ostringstream buf;
        // The variable name needs to contain "::"
        buf << cactus_varname
            // << "::" << varname1
            << " it=" << iteration << " tl=" << timelevel << " rl=" << reflevel
            << " c=" << component;
        const std::string datasetname = buf.str();
        std::cout << "datasetname " << datasetname << "\n";
        assert(!datasetnames.count(datasetname));
        datasetnames.insert(datasetname);

        const hid_t dataset =
            H5Dcreate(outfile, datasetname.c_str(), H5T_NATIVE_DOUBLE,
                      block_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        herr = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, data);
        assert(!herr);
        herr = H5Sclose(block_space);
        assert(!herr);

        {
          const hid_t scalar_space = H5Screate(H5S_SCALAR);

          {
            const hid_t attr =
                H5Acreate(dataset, "carpet_mglevel", H5T_NATIVE_INT,
                          scalar_space, H5P_DEFAULT, H5P_DEFAULT);
            const int carpet_mglevel = 0;
            herr = H5Awrite(attr, H5T_NATIVE_INT, &carpet_mglevel);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr =
                H5Acreate(dataset, "level", H5T_NATIVE_INT, scalar_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            herr = H5Awrite(attr, H5T_NATIVE_INT, &reflevel);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr =
                H5Acreate(dataset, "timestep", H5T_NATIVE_INT, scalar_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            herr = H5Awrite(attr, H5T_NATIVE_INT, &iteration);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr =
                H5Acreate(dataset, "group_timelevel", H5T_NATIVE_INT,
                          scalar_space, H5P_DEFAULT, H5P_DEFAULT);
            const int group_timelevel = 0;
            herr = H5Awrite(attr, H5T_NATIVE_INT, &group_timelevel);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr =
                H5Acreate(dataset, "time", H5T_NATIVE_DOUBLE, scalar_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            herr = H5Awrite(attr, H5T_NATIVE_DOUBLE, &timestep);
            herr = H5Aclose(attr);
          }

          {
            const hid_t datatype = H5Tcopy(H5T_C_S1);

            // const char *const name = varname;
            const char *const name = cactus_varname.c_str();
            herr = H5Tset_size(datatype, std::strlen(name) + 1);
            const hid_t attr =
                H5Acreate(dataset, "name", datatype, scalar_space, H5P_DEFAULT,
                          H5P_DEFAULT);
            herr = H5Awrite(attr, datatype, name);
            herr = H5Aclose(attr);
            herr = H5Tclose(datatype);
          }

          herr = H5Sclose(scalar_space);
        }

        // "cctk_bbox"
        // "cctk_nghostzones"

        {
          const hsize_t size = rank;
          const hid_t vector_space = H5Screate_simple(1, &size, NULL);

          {
            const hid_t attr =
                H5Acreate(dataset, "origin", H5T_NATIVE_DOUBLE, vector_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            double origin[3];
            for (int d = 0; d < rank; ++d)
              origin[d] = bbox[2 * d];
            herr = H5Awrite(attr, H5T_NATIVE_DOUBLE, origin);
            herr = H5Aclose(attr);
          }
          {
            const hid_t attr =
                H5Acreate(dataset, "delta", H5T_NATIVE_DOUBLE, vector_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            double delta[3];
            for (int d = 0; d < rank; ++d)
              delta[d] = (bbox[2 * d + 1] - bbox[2 * d]) / (dims[d] - 1);
            herr = H5Awrite(attr, H5T_NATIVE_DOUBLE, delta);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr =
                H5Acreate(dataset, "iorigin", H5T_NATIVE_INT, vector_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            int iorigin[3];
            using std::lrint;
            for (int d = 0; d < rank; ++d) {
              assert(lrint(1000 * (bbox[2 * d] - origin[d]) /
                           (coarse_spacing[d] / (1 << reflevel))) %
                         1000 ==
                     0);
              iorigin[d] = lrint((bbox[2 * d] - origin[d]) /
                                 (coarse_spacing[d] / (1 << reflevel)));
            }
            herr = H5Awrite(attr, H5T_NATIVE_INT, iorigin);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr =
                H5Acreate(dataset, "ioffset", H5T_NATIVE_INT, vector_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            const int ioffset[3] = {0, 0, 0};
            herr = H5Awrite(attr, H5T_NATIVE_INT, ioffset);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr =
                H5Acreate(dataset, "ioffsetdenom", H5T_NATIVE_INT, vector_space,
                          H5P_DEFAULT, H5P_DEFAULT);
            const int ioffsetdenom[3] = {1, 1, 1};
            herr = H5Awrite(attr, H5T_NATIVE_INT, ioffsetdenom);
            herr = H5Aclose(attr);
          }

          herr = H5Sclose(vector_space);
        }

        if (reflevel == 0) {
          int iorigin[3];
          using std::lrint;
          for (int d = 0; d < rank; ++d)
            iorigin[d] = lrint((bbox[2 * d] - origin[d]) /
                               (coarse_spacing[d] / (1 << reflevel)));
          for (int k = 0; k < dims[2]; ++k) {
            for (int j = 0; j < dims[1]; ++j) {
              for (int i = 0; i < dims[0]; ++i) {
                ++read_count(iorigin[0] + i, iorigin[1] + j, iorigin[2] + k);
              }
            }
          }
        }

        ++block;
      }

      gsfd_close(infile);
    } // for inputfile
  }

  herr = H5Fclose(outfile);

  std::cout << "Unread grid points:\n";
  for (int k = 0; k < gsh[2]; ++k) {
    for (int j = 0; j < gsh[1]; ++j) {
      for (int i = 0; i < gsh[0]; ++i) {
        if (read_count(i, j, k) == 0)
          std::cout << i << " " << j << " " << k << " " << read_count(i, j, k)
                    << "\n";
      }
    }
  }

  return 0;
}
