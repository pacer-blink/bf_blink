#include "antenna_positions.h"
#include <myparser.h>
#include <myfile.h>
#include <mystrtable.h>
#include <bg_fits.h>

CAntennaPositions::CAntennaPositions()
{}

CAntennaPositions::~CAntennaPositions()
{}

int CAntennaPositions::CalculateUVW( CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, bool bSaveFits /*=false*/ )
{
   if( size() <= 0 ){
      return 0;
   }

   // initalise :
   // set NaN to be consistent with CASA dump (was 0 but CASA dump script sets NaN for flagged antennas and upper half when option --conjugate is not used):
   fits_vis_u.SetNaN();
   fits_vis_v.SetNaN();
   fits_vis_w.SetNaN();

   int n_baselines = 0;
   for(int i=0;i<size();i++){
      InputMapping& ant1 = (*this)[i];
      
      for(int j=(i+1);j<size();j++){
         InputMapping& ant2 = (*this)[j];         
         
         double u = (ant1.x - ant2.x);
         double v = (ant1.y - ant2.y);
         double w = (ant1.z - ant2.z);
         
         // j,i (instead of i,j) to be consistent with CASA UVW array:
         fits_vis_u.setXY( j, i, u );
         fits_vis_v.setXY( j, i, v );
         fits_vis_w.setXY( j, i, w );
         
         n_baselines++;
      }
   }
   
   if( bSaveFits ){
      if( fits_vis_u.WriteFits( "u.fits" ) ){
        printf("ERROR : could not write output file u.fits\n");        
      }
      if( fits_vis_v.WriteFits( "v.fits" ) ){
        printf("ERROR : could not write output file v.fits\n");        
      }
      if( fits_vis_w.WriteFits( "w.fits" ) ){
        printf("ERROR : could not write output file w.fits\n");        
      }
   }
   
   return n_baselines;
}

int CAntennaPositions::ReadAntennaPositions(const char* filename, bool bConvertToXYZ /*=false*/ )
{
   if( filename && strlen(filename) ){
      if( !MyFile::DoesFileExist( filename ) ){
         printf("ERROR: filename %s does not exist\n",filename);
         return -1;
      }
   }else{
      printf("ERROR : empty filename provided -> cannot continue\n");
      return -1;      
   }

   clear();
   int n_ant_index = 0;
   MyFile file(filename);
   const char* pLine;
   if(!file.IsOpened()){
      file.Open( filename );
   }
   while(pLine = file.GetLine(TRUE)){
      if(strlen(pLine) && pLine[0]!='#'){
         MyParser pars=pLine;
         CMyStrTable items;
         pars.GetItems( items );
         
         if( strcmp(items[0].c_str(),"#") && items.size()>=4 ){
            InputMapping input;            
            input.szAntName = items[0].c_str();
            input.x = atof( items[1].c_str() );
            input.y = atof( items[2].c_str() );
            input.z = atof( items[3].c_str() );
            
            push_back( input );
         }
      }
   }
   file.Close();        

   return size();
}


