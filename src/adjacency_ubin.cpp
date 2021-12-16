/// @file     adjacency_ubin.cpp
/// @brief    Generate a (sparse) adjacency matrix for a point cloud.
///           This program generates a list of all pairs of points
///           that lie within a user-specified cutoff distance.
///           This program uses a simple uniform binning algorithm.
/// @author   Andrew Jewett, Steve Plimpton
/// @license  GPL-2.0
///
/// @note  To change the number of dimensions in which the point cloud lives,
/// (ie, the number of coordinates on each line of the structure file)
/// then compile using following compiler flag "-DG_DIM=4" (4 dimensions).
/// (You can also edit the "global_defs.h" file.)


#include <cassert> //needed for assert()
#include <cstring> //needed for strcmp()
#include <cstdlib> //needed for atol()
#include <map>
#include <vector>
using namespace std;


#include "global_defs.h"
#include "bins.h"
#include "io.h"


// global variables
const char g_program_name[]   = "adjacency_ubin";
const char g_version_string[] = "0.4.0";
const char g_date_string[]    = "<2020-4-19>";




struct Settings
{

  enum CoordFileFormat {COORD_FORMAT_RAW,
                        COORD_FORMAT_XYZ};

  CoordFileFormat coord_file_format;
  Real rcut; //cutoff distance for whether two atoms are "in contact"
  Icrd aUserNumBins[g_dim]; //let the user choose the number of bins

  Settings() {
    coord_file_format = COORD_FORMAT_RAW;
    rcut = 1.5;
    for (int d = 0; d < g_dim; d++)
      aUserNumBins[d] = UNINITIALIZED;
  }

}; //struct Settings



void PrintAdjacency(const map<pair<Iatm,Iatm>,long>& m,
                    ostream& out_file)
{
  for(map<pair<Iatm,Iatm>,long>::const_iterator p=m.begin(); p!=m.end(); p++)
  {
    pair<Iatm,Iatm> ijpair = p->first;
    Iatm ia = ijpair.first;
    Iatm ja = ijpair.second;
    long count = p->second;
    out_file << ia+1 << " " << ja+1 << " " << count << "\n";
  }
}





// THIS FUNCTION IS INEFFICIENT.  I SHOULD REWRITE IT
void 
CalcAdjacency(map<pair<Iatm,Iatm>,long>& m,
              ConstVect *aaX,  // <--> equivalent to "double (*x)[g_dim],"
              Iatm num_atoms,
              Real rcut,
              Bins &bins)
{
  Icrd aIcrd[g_dim];
  Icrd aJcrd[g_dim];
  Real rcutsq = rcut*rcut;

  for (Iatm ia = 0;  ia < num_atoms;  ia++)
  {
    Ibin ib = bins.aBinFromIatm[ia]; //the bin-id of bin containing atom ia
    // Find the array of coordinates (aIcrd) for that bin (ib):
    bins.IcrdFromIbin(aIcrd, ib);

    aJcrd[0] = ICRD_UNINITIALIZED; // <-signal this is the first neighbor bin.
    while (bins.NextNeighborBin(aIcrd, aJcrd))
    {
      // loop over all of the atoms in that bin
      Ibin jb  = bins.IbinFromIcrd(aJcrd);
      // original version:
      //Iatm ja = bins.aHeadFromIbin[jb]; // "ja" is the atom-id of a neighbor
      //                                  // candidate located in bin "jb"
      // memory efficient version:
      Ibin jbb = jb % bins.num_bins;
      Iatm ja  = bins.aHeadFromIbin[jbb];

      while (ja >= 0)
      {
        if ((DistanceSqd(aaX[ia], aaX[ja]) < rcutsq) && (ja > ia))
          m[ pair<Iatm,Iatm>(ia, ja) ] = 1;
        ja = bins.aNextFromIatm[ja];
      }
    }
  } //for (Iatm ia = 0;  ia < num_atoms;  ia++)
} //CalcAdjacency()



# if 0
// COMMENTING OUT
void 
CalcAdjacency(map<pair<Iatm,Iatm>,long>& m,
              ConstVect *aaX,  // <--> equivalent to "double (*x)[g_dim],"
              Iatm num_atoms,
              Real rcut,
              Bins &bins)
{
  Icrd aIcrd[g_dim];
  Icrd aJcrd[g_dim];
  for (Iatm ia = 0;  ia < num_atoms;  ia++)
  {
    Ibin ib = bins.aBinFromIatm[ia]; //the bin-id of bin containing atom ia
    //   Note: This is equivalent to the next two commented lines:
    //  IcrdFromX(aXX[ia], aIcrd);
    //  Ibin ib = Bins.IbinFromX(aax[ia]);

    Ibin jb=UNINITIALIZED; //jb is bin-id of one of the nearby bins
    while ((jb = bins.NextNeighborBin(ib, jb)) >= 0)
    {
      // loop over all of the atoms in that bin
      Iatm ja = bins.aHeadFromIbin[jb]; //ja is the atom-id an atom in bin jb
      while (ja >= 0)
      {
        if (DistanceSqd(aaX[ia], aaX[ja]) < rcut*rcut)
          m[ pair<Iatm,Iatm>(ia, ja) ] = 1;
        ja = bins.aNextFromIatm[ja];
      }
    }
  } //for (Iatm ia = 0;  ia < num_atoms;  ia++)
} //CalcAdjacency()
#endif //#if 0



void 
ParseArgs(int argc,
          char **argv,
          Settings &settings)
{

  stringstream explanation;
  explanation 
    <<"\n"
    "Explanation:\n"
    "\n"
    "  This program reads a coordinate file containing one (or more) structures\n"
    "  (in .raw or .xyz format), and generates a sparse matrix of contacts\n"
    "  i1 j2 count1\n"
    "  i2 j2 count2\n"
    "  i3 j3 count3\n"
    "   :  :   :\n"
    "  where \"i\" and \"j\" refer to monomers in the structure(s), and\n"
    "  \"count\" refers to the number of times those two monomers were found\n"
    "  to be in contact with eachother in this structure(s).\n"
    "    (\"in contact\" <--> to be spatially separated by a distance <= rcut)\n"
    "\n"
    "Typical Usage:\n"
    "\n"
    "\n"
    <<"  "<< g_program_name<<" -r rcut [optional arguments..] < coordinate_file > matrix_file\n"
    "\n"
    "\n"
    "Optional arguments:\n"
    "\n"
    "  -r rcut   <- specify the threshold contact distance (1.5 by default).\n"
    "\n"
    "  -raw         <-The default input file format: 3-column space-\n"
    "                 delimited text file with blank lines separating frames.\n"
    "                 Note: \".RAW\" format is assumed by default.\n"
    "  -xyz         <-Instead, assume the input coordinate file uses \".XYZ\" format.\n"
    "                 (This feature is not well tested.  Hopefully it works.)\n"
    "  -bins N      <-Specify the number of bins (aka voxels) in each direction.\n"
    "                 Alternatley, you can specify 3 numbers (in 3 dimensions)\n"
    "                 if you want it to vary in different directions.\n"
    "                 (This option has not been tested and may not work. 2013-60-6)\n"
    "\n";

  try {
    int which_arg = 1;
    while (which_arg < argc)
    {
      if (strcmp(argv[which_arg], "-raw") == 0)
      {
        settings.coord_file_format = Settings::COORD_FORMAT_RAW;
        which_arg += 1;
      }
      else if (strcmp(argv[which_arg], "-xyz") == 0)
      {
        settings.coord_file_format = Settings::COORD_FORMAT_XYZ;
        which_arg += 1;
      }
      else if ((strcmp(argv[which_arg], "-r") == 0) ||
	       (strcmp(argv[which_arg], "-d") == 0) ||
	       (strcmp(argv[which_arg], "-rcut") == 0))
      {
        if (which_arg+1 == argc) { 
          stringstream errmsg;
          errmsg << "\nError: expected a number following \""<<argv[which_arg]
                 << "\"\n";
          throw ArgParseErr(errmsg.str());
        }
        settings.rcut = atof(argv[which_arg+1]);
        which_arg += 2;
      }
      else if (strcmp(argv[which_arg], "-bins") == 0)
      {
        if (which_arg+1 == argc) { 
          stringstream errmsg;
          errmsg << "\nError: expected at least one number following \""<<argv[which_arg]
                 << "\"\n";
          throw ArgParseErr(errmsg.str());
        }

        vector<Icrd> vUserNumBins; //<-Store numbers here while reading arg list
        int j = which_arg+1;
        while ((j < argc) && 
               (strlen(argv[j]) != 0) &&
               isdigit(argv[j][0]))
        {
          vUserNumBins.push_back(atoi(argv[j]));
          j++;
        }
        // Then copy the numbers into the settings.aUserNumBins[] array:
        if (vUserNumBins.size() == 1)
          // That means use the same number of bins in all 3 directions
          for (int d = 0; d < g_dim; ++d) 
            settings.aUserNumBins[d] = vUserNumBins[0];
        else if (vUserNumBins.size() == g_dim) {
          for (int d = 0; d < g_dim; ++d) 
            settings.aUserNumBins[d] = vUserNumBins[d];
        }
        else {
          stringstream errmsg;
          errmsg << "Error: Wrong number of integers ("
                 << vUserNumBins.size()
                 << ")\n"
                 << "       following \""<<argv[which_arg]<<"\" argument.\n"
                 << "       Expected either 1 or "<<g_dim<<" integers.\n";
          throw ArgParseErr(errmsg.str());
        }
        which_arg = j;
      }
      else
      {
        stringstream errmsg;
        errmsg << "\nError: Unrecognized argument: \""<<argv[which_arg]<<"\"\n";
        throw ArgParseErr(errmsg.str());
      }
    } //while (which_arg < argc)
  }
  catch (ArgParseErr& e)
  {
    cerr << "--------------------------------------------------\n";
    cerr << "       Error in argument list: \n" << e.what() << endl;
    cerr << "--------------------------------------------------\n";
    cerr << "\n";
    //cerr << explanation.str() << endl;
    exit(-1);
  }

  //---- finished parsing the argument list ----

} //ParseArgs()






int
main(int argc, char **argv)
{

  //---- Load the coordinate file into a large buffer ----

  try 
  {
    long long line_count = 1; //counts all lines including blanks and comments

    map<pair<Iatm,Iatm>,long> m_tot;

    cerr << g_program_name   << ", v"
         << g_version_string << " "
         << g_date_string    << "\n";

    // Process the argument list
    Settings settings;
    ParseArgs(argc, argv, settings);

    cerr << "  reading input file...\n";

    long long frame_counter = 0;

    // Read one structure from the file and calculate its adjacency matrix.
    while(cin)
    {
      Vect *aaX = NULL;

      // Read in the next frame
      long num_atoms = UNINITIALIZED; //number of atoms in each snapshot/frame

      if (settings.coord_file_format == Settings::COORD_FORMAT_RAW)
        num_atoms = ReadCoordsRAW(cin, &aaX, line_count);
      else if (settings.coord_file_format == Settings::COORD_FORMAT_XYZ)
        num_atoms = ReadCoordsXYZ(cin, &aaX, line_count);
      else
        assert(0);

      if (num_atoms > 0)
        cerr << "    finished reading frame " << frame_counter+1 << endl;

      // We could be at a point in the file with trailing whitespace.  Check
      // to make sure that we did actually read something before we continue.

      if (num_atoms > 0) {

        Bins bins(aaX, 
                  num_atoms, 
                  settings.rcut,
                  settings.aUserNumBins);

        map<pair<Iatm,Iatm>,long> m; // m = adjacency matrix (sparse)

        //    Note: If you are unfamilliar with C++ "maps" or "pairs", 
        //          they are similar to python dictionaries and tuples.  See:
        //    http://www.cplusplus.com/reference/map/map/
        //    http://www.cplusplus.com/reference/utility/pair/

        CalcAdjacency(m,
                      aaX,
                      num_atoms,
                      settings.rcut,
                      bins);


	// --- COMMENTING OUT: making a separate matrix for each structure 
	// ---              is slow and uses a lot of disk space.          
        //Now print the adjacency matrix
        //stringstream fname;
        //fname << "out_adjacency_" << (frame_counter+1) << ".dat";
        //fstream out_file(fname.str().c_str(), ofstream::out);
        //PrintAdjacency(m, out_file);
        //out_file.close();
        // ---


        // Add this adjacency matrix to the total adjacency matrix.
        for(map<pair<Iatm,Iatm>,long>::const_iterator p = m.begin();
            p != m.end();
            p++)
        {
          pair<Iatm,Iatm> ijpair = p->first;
          long count = p->second;
          if (m_tot.find(ijpair) != m_tot.end())
            m_tot[ijpair] += count;
          else
            m_tot[ijpair] = count;
        }

        frame_counter++;
        //cerr << " done" << endl;
        cerr << "    finished processing frame " << frame_counter << endl;
      }
      else
        assert(! cin);

      delete [] aaX;

    } //while(cin)


    if (frame_counter > 0) {
      // Now print the total adjacency matrix
      cerr << "  Printing adjacency matrix" << endl;
      PrintAdjacency(m_tot, cout);
    }

  } //try
  catch (InputErr& e)
  {
    cerr << "       Error in input file(s): \n" << e.what() << endl;
    exit(-1);
  }

} //main()

