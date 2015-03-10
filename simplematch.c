#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "superalign.h"

int main(int argc, char *argv[])
{
  StarCatRec StarCatBaseProp,StarCatProp;
  SuperAlignParmRec SuperAlignParmProp;
  float ShiftX,ShiftY,RotAngle;
  char FN1[100],FN2[100],TransFN[100];
  FILE *fout;
  float CosF,SinF;

  sscanf(argv[1],"%s",FN1);
  sscanf(argv[2],"%s",FN2);
  sscanf(argv[3],"%g",&ShiftX);
  sscanf(argv[4],"%g",&ShiftY);
  sscanf(argv[5],"%g",&RotAngle);
  sscanf(argv[6],"%s",TransFN);
  
  ReadMatchinFiles(FN1,&StarCatBaseProp,0.0,0.0,0.0);
  ReadMatchinFiles(FN2,&StarCatProp,ShiftX,ShiftY,RotAngle);
  SuperAlignParmProp.TypDist = 0.5/20;
  SuperAlignParmProp.MaxDist = 2.0/20;
  SuperAlignParmProp.TolPair = 8.0/20;
  SuperAlignParmProp.MaxShiftErr = 200.0/20;
  SuperAlignParmProp.MaxRotErr = 3.0;

  TryAlignment(&StarCatBaseProp,&StarCatProp,&SuperAlignParmProp,
	       &ShiftX,&ShiftY,&RotAngle);

  fout = fopen(TransFN,"w");
  fprintf(fout,"linear\n");
  CosF = cos(M_PI*RotAngle/180);
  SinF = sin(M_PI*RotAngle/180);
  fprintf(fout,"a %.3f\n",ShiftX);
  fprintf(fout,"b %.6f\n",CosF);
  fprintf(fout,"c %.6f\n",SinF);
  fprintf(fout,"d %.3f\n",ShiftY);
  fprintf(fout,"e %.6f\n",-SinF);
  fprintf(fout,"f %.6f\n",CosF);
  fclose(fout);
}

