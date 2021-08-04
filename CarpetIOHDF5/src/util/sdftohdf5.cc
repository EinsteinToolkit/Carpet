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

#include <cassert>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

int main(int argc, char **argv) {
  herr_t herr;

  std::cout << "\n"
            << "----------------------------\n"
            << "Carpet SDF-to-HDF5 Converter\n"
            << "----------------------------\n"
            << "\n";

  if (argc != 7) {
    std::cerr << "Usage: " << argv[0]
              << " <groupname> <iteration> <timelevel> <ntimelevels> "
                 "<inputfile> <outputfile>\n"
              << "   e.g. " << argv[0]
              << " ADMBase::metric 0 0 3 foobar.sdf gxx.h5\n";
    return 1;
  }

  const char *const groupname = argv[1];
  // const int iteration = std::atoi(argv[2]);
  const int timelevel = std::atoi(argv[3]);
  const int ntimelevels = std::atoi(argv[4]);

  gft_sdf_file_data *const infile = gft_open_sdf_file(argv[5]);
  if (!infile) {
    std::cerr << "Could not open SDF input file \"" << argv[5] << "\"\n";
    return 2;
  }

  const hid_t outfile =
      H5Fcreate(argv[6], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (outfile < 0) {
    std::cerr << "Could not create HDF5 output file \"" << argv[6] << "\"\n";
    return 3;
  }

  // add a dummy group so that the HDF5 file is recognized as Carpet format
  const hid_t group = H5Gcreate(outfile, "Parameters and Global Attributes",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  {
    const hid_t dataspace = H5Screate(H5S_SCALAR);
    const hid_t attr = H5Acreate(group, "nioprocs", H5T_NATIVE_INT, dataspace,
                                 H5P_DEFAULT, H5P_DEFAULT);
    const int nioprocs = 1; // TODO
    herr = H5Awrite(attr, H5T_NATIVE_INT, &nioprocs);
    herr = H5Aclose(attr);
    herr = H5Sclose(dataspace);
  }
  herr = H5Gclose(group);

  {
    int component = 0;
    const double time_epsilon = 1.0e-6;
    bool have_origin_time = false;
    double origin_time;
    bool have_delta_time;
    double delta_time;
    int iteration = -1;
    bool have_origin = false;
    double origin[3];
    bool have_coarse_spacing = false;
    double coarse_spacing[3];
    int reflevel = -1;

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
      ++component;

      if (!have_origin_time)
        origin_time = timestep;
      have_origin_time = true;

      std::cout << std::setprecision(17) << "timestep=" << timestep
                << " origin_time=" << origin_time << "\n";
      assert(timestep >= origin_time - time_epsilon);

      if (timestep > origin_time + time_epsilon) {
        if (!have_delta_time)
          delta_time = timestep - origin_time;
        have_delta_time = true;

        using std::lrint;
        assert(lrint(1000000 * (timestep - origin_time) / delta_time) %
                   1000000 ==
               0);
      }

      using std::lrint;
      if (!have_delta_time)
        iteration = 0;
      else
        iteration = lrint((timestep - origin_time) / delta_time);

      if (!have_origin)
        for (int d = 0; d < rank; ++d)
          origin[d] = bbox[2 * d];
      have_origin = true;

      for (int d = 0; d < rank; ++d)
        assert(bbox[2 * d] >= origin[d]);

      assert(rank <= 3);
      double spacing[3];
      for (int d = 0; d < rank; ++d)
        spacing[d] = (bbox[2 * d + 1] - bbox[2 * d]) / (dims[d] - 1);

      if (!have_coarse_spacing)
        for (int d = 0; d < rank; ++d)
          coarse_spacing[d] = spacing[d];
      have_coarse_spacing = true;

      for (int d = 0; d < rank; ++d)
        assert(spacing[d] <= coarse_spacing[d]);

      using std::lrint;
      for (int d = 0; d < rank; ++d)
        assert(lrint(1000000 * coarse_spacing[d] / spacing[d]) % 1000000 == 0);

      assert(rank >= 1);
      using std::log2, std::lrint;
      reflevel = lrint(log2(coarse_spacing[0] / spacing[0]));
      assert(reflevel >= 0 && reflevel < 100);
      for (int d = 0; d < rank; ++d)
        assert(lrint(log2(coarse_spacing[d] / spacing[d])) == reflevel);

      std::cout << "Processing dataset \"" << varname << "\"\n"
                << "  timestep " << timestep << "\n"
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
                << "  reflevel " << reflevel << "\n";

      hsize_t h5dims[3];
      for (int d = 0; d < rank; ++d)
        h5dims[d] = dims[d];
      const hid_t dataspace = H5Screate_simple(rank, h5dims, nullptr);

      std::ostringstream buf;
      {
        const int component = 0;
        buf << varname << " it=" << iteration << " tl=" << timelevel
            << " rl=" << reflevel << " c=" << component;
      }
      const std::string datasetname = buf.str();

      const hid_t dataset =
          H5Dcreate(outfile, datasetname.c_str(), H5T_NATIVE_DOUBLE, dataspace,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      herr = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      data);
      assert(!herr);
      herr = H5Sclose(dataspace);
      assert(!herr);

      ///

      {
        const hid_t dataspace = H5Screate(H5S_SCALAR);

        {
          const hid_t attr = H5Acreate(dataset, "level", H5T_NATIVE_INT,
                                       dataspace, H5P_DEFAULT, H5P_DEFAULT);
          herr = H5Awrite(attr, H5T_NATIVE_INT, &reflevel);
          herr = H5Aclose(attr);
        }

        {
          const hid_t attr = H5Acreate(dataset, "timestep", H5T_NATIVE_INT,
                                       dataspace, H5P_DEFAULT, H5P_DEFAULT);
          herr = H5Awrite(attr, H5T_NATIVE_INT, &iteration);
          herr = H5Aclose(attr);
        }

        {
          const hid_t attr = H5Acreate(dataset, "time", H5T_NATIVE_DOUBLE,
                                       dataspace, H5P_DEFAULT, H5P_DEFAULT);
          herr = H5Awrite(attr, H5T_NATIVE_DOUBLE, &timestep);
          herr = H5Aclose(attr);
        }

        {
          const hid_t datatype = H5Tcopy(H5T_C_S1);
          herr = H5Tset_size(datatype, std::strlen(varname) + 1);
          const hid_t attr = H5Acreate(dataset, "name", datatype, dataspace,
                                       H5P_DEFAULT, H5P_DEFAULT);
          herr = H5Awrite(attr, datatype, varname);
          herr = H5Aclose(attr);
          herr = H5Tclose(datatype);
        }

        herr = H5Sclose(dataspace);

        // "cctk_bbox"
        // "cctk_nghostzones"

        {
          const hsize_t size = rank;
          const hid_t dataspace = H5Screate_simple(1, &size, NULL);

          {
            const hid_t attr = H5Acreate(dataset, "origin", H5T_NATIVE_DOUBLE,
                                         dataspace, H5P_DEFAULT, H5P_DEFAULT);
            double origin[3];
            for (int d = 0; d < rank; ++d)
              origin[d] = bbox[2 * d];
            herr = H5Awrite(attr, H5T_NATIVE_DOUBLE, origin);
            herr = H5Aclose(attr);
          }
          {
            const hid_t attr = H5Acreate(dataset, "delta", H5T_NATIVE_DOUBLE,
                                         dataspace, H5P_DEFAULT, H5P_DEFAULT);
            double delta[3];
            for (int d = 0; d < rank; ++d)
              delta[d] = (bbox[2 * d + 1] - bbox[2 * d]) / (dims[d] - 1);
            herr = H5Awrite(attr, H5T_NATIVE_DOUBLE, delta);
            herr = H5Aclose(attr);
          }

          {
            const hid_t attr = H5Acreate(dataset, "iorigin", H5T_NATIVE_INT,
                                         dataspace, H5P_DEFAULT, H5P_DEFAULT);
            int iorigin[3];
            using std::lrint;
            for (int d = 0; d < rank; ++d)
              iorigin[d] = lrint((bbox[2 * d] - origin[d]) /
                                 (coarse_spacing[d] / (1 << reflevel)));
            herr = H5Awrite(attr, H5T_NATIVE_INT, iorigin);
            herr = H5Aclose(attr);
          }
        }
      }

      ///

      ++component;
    }
  }

  return 0;
}
