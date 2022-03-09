#ifndef _ANTENNA_POSITIONS__
#define _ANTENNA_POSITIONS__

#include <vector>
#include <string>

using namespace std;

class CBgFits;

class InputMapping
{
public :
   int input;
   int antenna;
   string szAntName;
   char pol;
   int delta;
   int flag;  
   
   double x;
   double y;
   double z;

   InputMapping() 
   : input(-1), antenna(-1), pol('U'), delta(0), flag(0), x(0), y(0), z(0)
   {};

//   static int read_mapping_file( std::vector<InputMapping>& inputs , const char* filename="instr_config.txt" );
};

class CAntennaPositions : public std::vector<InputMapping> 
{
public :
   CAntennaPositions();
   ~CAntennaPositions();
   
   int ReadAntennaPositions(const char* filename, bool bConvertToXYZ=false);
   int CalculateUVW( CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, bool bSaveFits=false ); // in meters for now 
};

#endif 