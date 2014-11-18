/* 
Package Name: Superalign 1.0
Author: Rychard J. Bouwens
Written: February 2003
Modified: April 2005

Brief Description: 

Usage: 

*/

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "arrays.h"
#include "superalign.h"
#include "chdecomp.c"
#include "chsvx.c"

/* 

Here is the required format for the EXPOSURE file:
[# of exposures]
[exposure #1]  [x-shift guess #1] [y-shift guess #1] [rot-angle guess #1]
[exposure #2]  [x-shift guess #2] [y-shift guess #2] [rot-angle guess #2]
[exposure #3]  [x-shift guess #3] [y-shift guess #3] [rot-angle guess #3]
....


*/

extern void CopyStarCatProp(StarCatRec *PStarCatPropD,
			    StarCatRec *PStarCatPropS)
/* This routine copies the contents of *PStarCatPropS to *PStarCatPropD */
{
  int I;

  strcpy(PStarCatPropD->FN,PStarCatPropS->FN);
  PStarCatPropD->NObj = PStarCatPropS->NObj;
  PStarCatPropD->HDR_DX = PStarCatPropS->HDR_DX;
  PStarCatPropD->HDR_DY = PStarCatPropS->HDR_DY;
  PStarCatPropD->HDR_ROT = PStarCatPropS->HDR_ROT;
  PStarCatPropD->PX = (float *)calloc(PStarCatPropD->NObj,sizeof(float));
  PStarCatPropD->PY = (float *)calloc(PStarCatPropD->NObj,sizeof(float));
  PStarCatPropD->PM = (float *)calloc(PStarCatPropD->NObj,sizeof(float));
  for (I=0;I<PStarCatPropD->NObj;I++) {
    PStarCatPropD->PX[I] = PStarCatPropS->PX[I];
    PStarCatPropD->PY[I] = PStarCatPropS->PY[I];
    PStarCatPropD->PM[I] = PStarCatPropS->PM[I];
  }
}

extern void DestructStarCatProp(StarCatRec *PStarCatProp)
/* This routine deallocates the memory used by *PStarCatProp */
{
  free(PStarCatProp->PX);
  free(PStarCatProp->PY);
  free(PStarCatProp->PM);
}

extern void GroupExposures(int *PNPoint, PointingRec **PPPointingProp,
			   AllExposureRec *PAllExposureProp,
			   float GroupSize, int FirstRefImage)
/* This routine takes all the exposures in *PAllExposureProp and breaks
   them into distinct groups.  Each group is placed in a separate
   PointingRec structure.  *PNPoint contains the number of groups which
   were determined.
   GroupSize gives the maximum separation that two pointings can have
      to be considered part of the same group.
   FirstRefImage indicates whether the first exposure is the reference image and
      should be place in its own group.  */
{
  int I,NExp,Count,Count2,MaxCount2;
  int *NNeigh,**NeighList;
  float *PX,*PY,*PR,Dist;
  int *PCheck;  /* Indicates whether an exposure has been assigned a
                   group yet, 0 == not assigned, 1 = assigned */
  int Wh,J,K,Wh2,Loop;

  fprintf(stdout,"Grouping All Exposure into Sets of Pointings...\n");
  NExp = PAllExposureProp->NExp;
  PX = (float *)calloc(NExp,sizeof(float));
  PY = (float *)calloc(NExp,sizeof(float));
  PR = (float *)calloc(NExp,sizeof(float));
  PCheck = (int *)calloc(NExp,sizeof(int));
  
  
  for (I=0;I<NExp;I++) {
    ExposureRec *PExposureProp;

    PExposureProp = &PAllExposureProp->PExposureProp[I];
    PX[I] = PExposureProp->hdr_dx;
    PY[I] = PExposureProp->hdr_dy;
    PR[I] = PExposureProp->hdr_rot;
  }
  NNeigh = (int *)calloc(NExp,sizeof(int));
  NeighList = (int **)calloc(NExp,sizeof(int *));
  NeighborList(PX,PY,NExp,PX,PY,NExp,GroupSize*2,
	       NNeigh,NeighList);
  *PNPoint = NExp;
  *PPPointingProp = (PointingRec *)calloc(NExp,sizeof(PointingRec));
  *PNPoint = 0;

  if (FirstRefImage) {
    PointingRec *PPointingProp;
    
    PPointingProp = &(*PPPointingProp)[*PNPoint];
    PPointingProp->NExp = 1;
    PPointingProp->PStarCatProp = (StarCatRec *)
      calloc(PPointingProp->NExp,sizeof(StarCatRec));
    CopyStarCatProp(&PPointingProp->PStarCatProp[0],
		    &PAllExposureProp->PExposureProp[0].StarCatProp);
    PCheck[0] = 1;
    (*PNPoint)++;
  }

  for (Loop=0;Loop<3;Loop++)
    for (I=0;I<NExp;I++) {
      if (!PCheck[I]) {
        /* Exposure not assigned yet, so go through a see how many 
           neighbors it has */
	Count = 0;
	MaxCount2 = 0;
	for (J=0;J<NNeigh[I];J++) {
	  Wh = NeighList[I][J];
	  if (!PCheck[Wh]) {	  
	    Dist = sqrt((PX[I] - PX[Wh])*(PX[I] - PX[Wh]) +
			(PY[I] - PY[Wh])*(PY[I] - PY[Wh]));
	    if (fabs(Dist) < GroupSize)
	      Count++;
	    
            /* Count the number of neighbors each neighbor has */
	    Count2 = 0;  
	    for (K=0;K<NNeigh[Wh];K++) {
	      Wh2 = NeighList[Wh][K];
	      
	      if (!PCheck[Wh2]) {	  
		Dist = sqrt((PX[Wh] - PX[Wh2])*(PX[Wh] - PX[Wh2]) +
			    (PY[Wh] - PY[Wh2])*(PY[Wh] - PY[Wh2]));
		if (fabs(Dist) < GroupSize)
		  Count2++;
	      }
	    }
	    if (Count2 > MaxCount2)
	      MaxCount2 = Count2;
	  }
	}
        /* If the current object has a greater or equal number of
           neighbors to any of its neighbors, use it to define a
           group/pointing */
	if (Count >= MaxCount2) {
	  PointingRec *PPointingProp;
	  
	  PPointingProp = &(*PPPointingProp)[*PNPoint];
	  PPointingProp->NExp = Count;
	  PPointingProp->PStarCatProp = (StarCatRec *)
	    calloc(PPointingProp->NExp,sizeof(StarCatRec));
	  Count = 0;
	  for (J=0;J<NNeigh[I];J++) {
	    Wh = NeighList[I][J];
	    
	    if (!PCheck[Wh]) {	  
	      Dist = sqrt((PX[I] - PX[Wh])*(PX[I] - PX[Wh]) +
			  (PY[I] - PY[Wh])*(PY[I] - PY[Wh]));
	      if (fabs(Dist) < GroupSize) {
		CopyStarCatProp(&PPointingProp->PStarCatProp[Count],
				&PAllExposureProp->PExposureProp[Wh].StarCatProp);
		PCheck[Wh] = 1;
		Count++;
	      }
	    }
	  }
	  (*PNPoint)++;
	}
      }
    }
  fprintf(stdout,"Found %i Pointings...\n",*PNPoint);
  
  if (*PNPoint > NExp)
    *PPPointingProp = (PointingRec *)
      realloc(*PPPointingProp,NExp*sizeof(PointingRec));
  free(PX);
  free(PY);
  free(PR);
  free(PCheck);
  for (I=0;I<NExp;I++)
    free(NeighList[I]);
  free(NeighList);
  free(NNeigh);
}

typedef struct {
  int *NNeigh;      /* Number of Neighbors for specific objects in
                       a catalog */
  int **NeighList;  /* Indices to the Neighbors for specific objects 
                       in a catalog */
} NeighRec;

extern void MakeCoordListProp(StarCatRec *PStarCatProp,
			      CoordListRec *PCoordListProp)
/* Transform the object positions from their positions on the original
   frames (*PStarCatProp) to their position in the universal frame 
   (*PCoordListProp) */
{
  int I;
  float cosf,sinf;

  PCoordListProp->NObj = PStarCatProp->NObj;
  PCoordListProp->PX = (float *)calloc(PCoordListProp->NObj,sizeof(float));
  PCoordListProp->PY = (float *)calloc(PCoordListProp->NObj,sizeof(float));
  PCoordListProp->PM = (float *)calloc(PCoordListProp->NObj,sizeof(float));
  cosf = cos(M_PI*PStarCatProp->HDR_ROT/180);
  sinf = sin(M_PI*PStarCatProp->HDR_ROT/180);
  for (I=0;I<PCoordListProp->NObj;I++) {
    float TX,TY;

    TX = PStarCatProp->PX[I];
    TY = PStarCatProp->PY[I];
    
    PCoordListProp->PX[I] = cosf*TX + sinf*TY + 
      + PStarCatProp->HDR_DX;
    PCoordListProp->PY[I] = -sinf*TX + cosf*TY + 
      + PStarCatProp->HDR_DY;
    PCoordListProp->PM[I] = PStarCatProp->PM[I];
  }
}

extern void DestructCoordListProp(CoordListRec *PCoordListProp)
/* Deallocate memory in *PCoordListProp */
{
  free(PCoordListProp->PX);
  free(PCoordListProp->PY);
  free(PCoordListProp->PM);
}

extern void OutputCoordListProp(CoordListRec *PCoordListProp, char *FN)
/* Output positions in *PCoordListProp to file FN */
{
  FILE *fout;
  int I;

  fout = fopen(FN,"w");
  for (I=0;I<PCoordListProp->NObj;I++) {
    fprintf(fout,"%i %g %g %g\n",I,
	    PCoordListProp->PX[I],PCoordListProp->PY[I],
	    PCoordListProp->PM[I]);
  }
  fclose(fout);
}

extern void OutputStarCatProp(StarCatRec *PStarCatProp, char *FN)
/* Output positions in *PStarCatProp to file FN */
{
  FILE *fout;
  int I;

  fout = fopen(FN,"w");
  for (I=0;I<PStarCatProp->NObj;I++) {
    fprintf(fout,"%i %g %g %g\n",I,
	    PStarCatProp->PX[I],PStarCatProp->PY[I],
	    PStarCatProp->PM[I]);
  }
  fclose(fout);
}

extern int CountMatches(CoordListRec *PCoordListProp0, int J,
			CoordListRec *PCoordListProp, int NExp,
			NeighRec *PNeighProp0,
			int MakeMatchList, int **PPWhMatch,
			int NMatch)
/* Count the number of matches to object J from PCoordListProp0 
   in all NExp of PCoordListProp.

   PNeighProp0 contains all neighbor lists relative to all NExp 

   if MakeMatchList == 1, record changes in **PWhMatch 
*/
{
  int K,L,Count;
  
  Count = 0;
  for (K=0;K<NExp;K++) {
    CoordListRec *PCoordListProp1;
    int Good,Wh,WhG;
    float DiffX,DiffY,Diff;
    NeighRec *PNeighProp;
    
    PCoordListProp1 = &PCoordListProp[K];
    PNeighProp = &PNeighProp0[K];
    Good = 0;
    for (L=0;L<PNeighProp->NNeigh[J];L++) {
      Wh = PNeighProp->NeighList[J][L];
      DiffX = PCoordListProp0->PX[J]-PCoordListProp1->PX[Wh];
      DiffY = PCoordListProp0->PY[J]-PCoordListProp1->PY[Wh];
      Diff = DiffX*DiffX + DiffY*DiffY;
      if (Diff < 1.5) {
	WhG = Wh;
	Good++;
      }
    }
    if (Good==1) {
      if (MakeMatchList) {
	PPWhMatch[K][NMatch] = WhG;
      }
      Count++;
    }
  }
  return Count;
}

extern int FlagMatches(int J, CoordListRec *PCoordListProp, int NExp,
		       NeighRec *PNeighProp, int **PPCheck)
{
  int K,L,Wh;
  
  for (K=0;K<NExp;K++) {
    NeighRec *PNeighProp0;

    PNeighProp0 = &PNeighProp[K];
    for (L=0;L<PNeighProp0->NNeigh[J];L++) {
      Wh = PNeighProp0->NeighList[J][L];
      PPCheck[K][Wh] = 1;
    }
  }
}

extern void GStackStarCatProp(int NExp, StarCatRec *PStarCatProp,
			      int Necessary, SuperAlignParmRec *PSAPP,
			      int *PNMatch, int **PPWhMatch)
/* This routine groups all of the object lists from each exposure
   into a composite list.  Necessary specifies the number of occurrences 
   an object must have to make it into the final list.   If Necessary
   == 1, then the operation is effectively a logical OR-ing of all
   object lists.  If Necessary == NExp, then the operation is
   effectively a logical AND-ing of all object lists.

   The number of objects in the final list is *PNMatch.   

   */
{
  int N1Max,I,MaxObj,J,K,L;
  int **PPCheck;
  float TX,TY;
  NeighRec **PPNeighProp;
  CoordListRec *PCoordListProp;

  N1Max = NExp;

  /* if an object must occur at least 4 times, then require that 
     use first 6 exposures as a basis for the final object list.  
     As a result, any object which is to appear in the final 
     "stacked" list must occur in the first six exposures. */
  if (Necessary > 3)
    N1Max = 6;
  
  /* First convert PStarCatProp into the universal frame 
     and initialize all of the objects in each list as unassigned */

  PPCheck = (int **)calloc(NExp,sizeof(int *));
  PCoordListProp = (CoordListRec *)
    calloc(NExp,sizeof(CoordListRec));
  MaxObj = 0;
  for (I=0;I<NExp;I++) {
    PPCheck[I] = (int *)calloc(PStarCatProp[I].NObj,sizeof(int));
    MakeCoordListProp(&PStarCatProp[I],&PCoordListProp[I]);
    MaxObj += PStarCatProp[I].NObj;
  } 

  /* Go through and make neighbor lists for all the objects in
     the first N1Max exposures */

  PPNeighProp = (NeighRec **)calloc(N1Max,sizeof(NeighRec *));
  for (I=0;I<N1Max;I++) {    
    StarCatRec *PStarCatProp0;
    CoordListRec *PCoordListProp0;

    PStarCatProp0 = &PStarCatProp[I];
    PCoordListProp0 = &PCoordListProp[I];
    PPNeighProp[I] = (NeighRec *)calloc(NExp,sizeof(NeighRec));
    for (J=0;J<NExp;J++) {
      NeighRec *PNeighProp;
      StarCatRec *PStarCatProp1;
      CoordListRec *PCoordListProp1;

      PStarCatProp1 = &PStarCatProp[J];
      PCoordListProp1 = &PCoordListProp[J];
      PNeighProp = &PPNeighProp[I][J];
      PNeighProp->NNeigh = (int *)calloc(PStarCatProp0->NObj,sizeof(int));
      PNeighProp->NeighList = 
	(int **)calloc(PStarCatProp0->NObj,sizeof(int *));
      NeighborList(PCoordListProp0->PX,PCoordListProp0->PY,
		   PCoordListProp0->NObj,
		   PCoordListProp1->PX,PCoordListProp1->PY,
		   PCoordListProp1->NObj,PSAPP->MaxDist*0.7,
		   PNeighProp->NNeigh,PNeighProp->NeighList);
    }
  }

  /* Go through all the catalogs and associate all objects with
     no object (i.e., = -1) in the final stacking.  */

  for (I=0;I<NExp;I++) {
    PPWhMatch[I] = (int *)calloc(MaxObj,sizeof(int));
    for (J=0;J<MaxObj;J++)
      PPWhMatch[I][J] = -1;
  }
  *PNMatch = 0;

  /* Go through first N1Max Exposures and see if the objects on
     these exposures match up with more than (or equal) 
     Necessary.  If so, then identify this object as a match
     and flag it in PPCheck to remove this object (and those
     it matches) from the total pool of objects.  */

  for (I=0;I<N1Max;I++) {    
    StarCatRec *PStarCatProp0;
    CoordListRec *PCoordListProp0;

    PStarCatProp0 = &PStarCatProp[I];
    PCoordListProp0 = &PCoordListProp[I];
    for (J=0;J<PStarCatProp0->NObj;J++) {
      if (!PPCheck[I][J]) {
	int MaxCount,Count2,Count;
	
	MaxCount = 0;
	Count = CountMatches(PCoordListProp0,J,PCoordListProp,NExp,
			     PPNeighProp[I],0,NULL,0);
	if (Count>=Necessary) {
	  for (K=0;K<N1Max;K++) {
	    StarCatRec *PStarCatProp1;
	    CoordListRec *PCoordListProp1;
	    NeighRec *PNeighProp;
	    int Good;
	    
	    PStarCatProp1 = &PStarCatProp[K];
	    PCoordListProp1 = &PCoordListProp[K];
	    PNeighProp = &PPNeighProp[I][K];
	    Good = 0;
	    for (L=0;L<PNeighProp->NNeigh[J];L++) {
	      int Wh;

	      Wh = PNeighProp->NeighList[J][L];
	      
	      if (!PPCheck[K][Wh]) {
		Count2 = CountMatches(&PCoordListProp[K],Wh,
				      PCoordListProp,NExp,
				      PPNeighProp[K],0,NULL,0);
		if (Count2 > MaxCount)
		  MaxCount = Count2;
	      }
	    }
	  }
	  if (Count >= MaxCount) {
	    CountMatches(PCoordListProp0,J,
			 PCoordListProp,NExp,
			 PPNeighProp[I],1,PPWhMatch,
			 *PNMatch);
	    (*PNMatch)++;
	    FlagMatches(J,PCoordListProp,NExp,PPNeighProp[I],
			PPCheck);
	  }
	}
      }
    }
  }
  for (I=0;I<NExp;I++)
    DestructCoordListProp(&PCoordListProp[I]);
  free(PCoordListProp);
  for (I=0;I<NExp;I++)
    free(PPCheck[I]);
  free(PPCheck);
  for (I=0;I<N1Max;I++) {
    StarCatRec *PStarCatProp0;

    PStarCatProp0 = &PStarCatProp[I];
    for (J=0;J<NExp;J++) {
      NeighRec *PNeighProp;

      PNeighProp = &PPNeighProp[I][J];
      for (K=0;K<PStarCatProp0->NObj;K++)
	free(PNeighProp->NeighList[K]);
      free(PNeighProp->NeighList);
      free(PNeighProp->NNeigh);      
    }
    free(PPNeighProp[I]);
  }
  free(PPNeighProp);
  if (*PNMatch>0) {
    for (I=0;I<NExp;I++) {
      if (NExp < MaxObj)
	PPWhMatch[I] = (int *)
	  realloc(PPWhMatch[I],(*PNMatch)*sizeof(int));
    }
  }
  else {
    for (I=0;I<NExp;I++) {
      free(PPWhMatch[I]);
    }
    free(PPWhMatch);
  }
}

extern void GStackStarCatProp2(int NExp, StarCatRec *PStarCatProp,
			       StarCatRec *PNStarCatProp,
			       int Necessary, SuperAlignParmRec *PSAPP,
			       char *FN, int OutputStars)
/* OutputStars generates an output file from each of the inputs
   which contains the objects which are not rejected in the final stack.
   The output name is the same are the input name for the star list, but
   with ".stars" appended.
   The average positions and magnitudes of specific objects in the 
   stack are placed in *PNStarCatProp. */
{
  int **PPWhMatch,I,J,NMatch;
  CoordListRec *PCoordListProp;
  
  strcpy(PNStarCatProp->FN,FN);
  PNStarCatProp->HDR_DX = 0.0;
  PNStarCatProp->HDR_DY = 0.0;
  PNStarCatProp->HDR_ROT = 0.0;

  PPWhMatch = (int **)calloc(NExp,sizeof(int *));
  GStackStarCatProp(NExp,PStarCatProp,Necessary,PSAPP,
		    &NMatch,PPWhMatch);
  PCoordListProp = (CoordListRec *)
    calloc(NExp,sizeof(CoordListRec));
  for (I=0;I<NExp;I++) {
    MakeCoordListProp(&PStarCatProp[I],&PCoordListProp[I]);
  }

  if (OutputStars) {   
    for (I=0;I<NExp;I++) {
      FILE *fout;
      char FN[200];
      int Count;      
      
      Count = 0;
      sprintf(FN,"%s.stars",PStarCatProp[I].FN);      
      fout = fopen(FN,"w");
      for (J=0;J<NMatch;J++)
	if (PPWhMatch[I][J] >= 0) {	  
	  fprintf(fout,"%i %g %g %g\n",
		  Count,PStarCatProp[I].PX[PPWhMatch[I][J]],
		  PStarCatProp[I].PY[PPWhMatch[I][J]],
		  PStarCatProp[I].PM[PPWhMatch[I][J]]);
	  Count++;	 
	}
      fclose(fout);
    }
  }
  

  PNStarCatProp->NObj = NMatch;
  PNStarCatProp->PX = (float *)calloc(NMatch,sizeof(float));
  PNStarCatProp->PY = (float *)calloc(NMatch,sizeof(float));
  PNStarCatProp->PM = (float *)calloc(NMatch,sizeof(float));
  for (I=0;I<NMatch;I++) {
    float SumX,SumY,SumM;
    int Count;

    SumX = 0;
    SumY = 0;
    SumM = 0;
    Count = 0;
    for (J=0;J<NExp;J++) {
      if (PPWhMatch[J][I] >= 0) {
	SumX += PCoordListProp[J].PX[PPWhMatch[J][I]];
	SumY += PCoordListProp[J].PY[PPWhMatch[J][I]];
	SumM += PCoordListProp[J].PM[PPWhMatch[J][I]];
	Count++;
      }
    }
    SumX /= Count;
    SumY /= Count;
    SumM /= Count;
    PNStarCatProp->PX[I] = SumX;
    PNStarCatProp->PY[I] = SumY;
    PNStarCatProp->PM[I] = SumM;
  }
  for (I=0;I<NExp;I++)
    DestructCoordListProp(&PCoordListProp[I]);
  free(PCoordListProp);
  if (NMatch > 0) {
    for (I=0;I<NExp;I++) {
      free(PPWhMatch[I]);
    }
    free(PPWhMatch);
  }
}

extern void MergeStarCatProp(int NPoint, StarCatRec *PStarCatProp,
			     StarCatRec *PStarCatPropOut, char *FN,
			     SuperAlignParmRec *PSAPP)
/* Take NPoint different star lists and output the merged / stacked
   list in *PStarCatPropOut */
{
  GStackStarCatProp2(NPoint,PStarCatProp,
		     PStarCatPropOut,1,PSAPP,FN,0);
}

extern void MakeStackStarCatProp(int NExp, StarCatRec *PStarCatProp,
				 StarCatRec *PNStarCatProp, char *FN,
				 int OutputStars, SuperAlignParmRec *PSAPP)
{
  int Necessary;

  /* Necessary specifies the number of occurrences a given object must
     have in the final stack to make it qualify as a "real" object and
     not as a likely cosmic ray. */

  if (NExp == 1) Necessary = 1;
  else if (NExp == 2) Necessary = 2;
  else if (NExp == 3) Necessary = 2;
  else if (NExp == 4) Necessary = 2;
  else if (NExp == 5) Necessary = 2;
  else if (NExp == 6) Necessary = 2;
  else if (NExp == 7) Necessary = 3;
  else if (NExp == 8) Necessary = 3;
  else if (NExp == 9) Necessary = 3;
  else if (NExp == 10) Necessary = 3;
  else if (NExp == 11) Necessary = 3;
  else if (NExp == 12) Necessary = 3;
  else if (NExp == 13) Necessary = 4;
  else if (NExp == 14) Necessary = 4;
  else if (NExp == 15) Necessary = 4;
  else if (NExp == 16) Necessary = 4;
  else if (NExp == 17) Necessary = 4;
  else if (NExp == 18) Necessary = 4;
  else if (NExp == 19) Necessary = 4;
  else if (NExp == 20) Necessary = 4;
  else if (NExp > 20) 
    Necessary = NExp/6+1;
  else if (NExp > 31) 
    Necessary = NExp/8+2;

  fprintf(stdout,"Stacking Images to generate star list.\n");
  fprintf(stdout,"GOOD STAR = SEEN in %i out of %i images.\n",
	  Necessary,NExp);
  GStackStarCatProp2(NExp,PStarCatProp,PNStarCatProp,Necessary,
		     PSAPP,FN,OutputStars);
  fprintf(stdout,"%i stars in list.\n",PNStarCatProp->NObj);
}

extern int TransOK(StarCatRec *PStarCatProp,
		   StarCatRec *PTStarCatProp, int Stand,
		   SuperAlignParmRec *PSAPP)
/* Compares the already existing stacked catalog *PStarCatProp
   with the catalog to be combined with it *PTStarCatProp (with 
   an implicitly new HDR_DX,HDR_DY,HDR_ROT).  Routine returns
   1 if successful and 0 if not.  The variable Stand indicates the
   standard that should be upheld in this asssessment.  
   Use Stand = 1 for a high standard, Stand = 0 for a comparably
   lower standard, Stand = -1 for a very low standard */
{
  StarCatRec *PStarCatProp0,StarCatProp;
  int OK,MinObj;

  if (PStarCatProp->NObj < PTStarCatProp->NObj)
    MinObj = PStarCatProp->NObj;
  else
    MinObj = PTStarCatProp->NObj;
    
  PStarCatProp0 = (StarCatRec *)calloc(2,sizeof(StarCatRec));
  CopyStarCatProp(&PStarCatProp0[0],PStarCatProp);
  CopyStarCatProp(&PStarCatProp0[1],PTStarCatProp);
  GStackStarCatProp2(2,PStarCatProp0,&StarCatProp,2,PSAPP,"TEST",0);
  OK = 0;
  if (Stand==1) {
    if (StarCatProp.NObj > 10) OK = 1;
    if ((StarCatProp.NObj > MinObj / 3) || (StarCatProp.NObj > 4.5)) OK = 1;
  }
  else if (Stand == 0) {
    if (StarCatProp.NObj > 2.9) OK = 1;
  }
  else if (Stand == -1) {
    if (StarCatProp.NObj > 1.9) OK = 1;
  }
  DestructStarCatProp(&PStarCatProp0[0]);
  DestructStarCatProp(&PStarCatProp0[1]);
  DestructStarCatProp(&StarCatProp);
  free(PStarCatProp0);

  return OK;
}

extern void AttemptAlign(StarCatRec *PStarCatProp1, 
			 StarCatRec *PStarCatProp2, int *POK,
			 SuperAlignParmRec *PSAPP, int Stand)
/* Attempts to align *PStarCatProp2 relative to *PStarCatProp1,
   Return *POK = 1 if the alignment was successful and
   *POK = 0 if the alignment failed.
   Stand gives the standard level this alignment will be required to meet.
*/
{
  float cosf,sinf,BestdX,BestdY,BestRotAngle,NX,NY;
  StarCatRec StarCatProp;

  TryAlignment(PStarCatProp1,PStarCatProp2,PSAPP,
	       &BestdX,&BestdY,&BestRotAngle);
  CopyStarCatProp(&StarCatProp,PStarCatProp2);
  /* Since BestdX, BestdY, and BestRotAngle give the best-fit relative 
     alignment parameters of *PStarCatProp2 to *PStarCatProp1, it is 
     necessary to convert these quantities from the *PStarCatProp1 
     frame to the universal frame */
  StarCatProp.HDR_ROT = PStarCatProp1->HDR_ROT + BestRotAngle;
  if (StarCatProp.HDR_ROT > 180)
    StarCatProp.HDR_ROT -= 360.0;
  else if (StarCatProp.HDR_ROT < -180)
    StarCatProp.HDR_ROT += 360.0;
  cosf = cos(PStarCatProp1->HDR_ROT*M_PI/180);
  sinf = sin(PStarCatProp1->HDR_ROT*M_PI/180);
  NX = BestdX*cosf + BestdY*sinf;
  NY = -BestdX*sinf + BestdY*cosf;
  StarCatProp.HDR_DX = NX + PStarCatProp1->HDR_DX;
  StarCatProp.HDR_DY = NY + PStarCatProp1->HDR_DY;

  /* Check to see if there are a sufficient number of matches to
     "guarantee" that this is a good alignment. */
  *POK = TransOK(PStarCatProp1,&StarCatProp,Stand,PSAPP);
  if (*POK) {
    DestructStarCatProp(PStarCatProp2);
    CopyStarCatProp(PStarCatProp2,&StarCatProp);
    fprintf(stdout,"Alignment Successful\n");
    fprintf(stdout,"Solution: dx = %.2f, dy = %.2f, theta = %.2f\n",
	    StarCatProp.HDR_DX,StarCatProp.HDR_DY,
	    StarCatProp.HDR_ROT);
  }
  else {
    fprintf(stdout,"Alignment Failed\n");
    fprintf(stdout,"Attempted Solution: dx = %.2f, dy = %.2f, theta = %.2f\n",
	    StarCatProp.HDR_DX,StarCatProp.HDR_DY,
	    StarCatProp.HDR_ROT);
  }
  DestructStarCatProp(&StarCatProp);
}

extern float CalcCurrentError(int NPoint, StarCatRec *PStarCatProp,
			      int When, SuperAlignParmRec *PSAPP)
/* Calculates the RMS error on aligned objects for NPoint number of
   different exposures.  Undistorted object positions are given in
   the PStarCatProp array.  The variable "When" only has an influence
   of what is sent to stdout.  For When == 0, it is "Start:" and
   for When == 1, it is "End:" */
{
  int I,**PPWhMatch,NMatch,WhJ,WhK,J,K,Count;
  CoordListRec *PCoordListProp;
  float dX,dY,Diff,Err;

  PPWhMatch = (int **)calloc(NPoint,sizeof(int *));
  GStackStarCatProp(NPoint,PStarCatProp,1,PSAPP,&NMatch,PPWhMatch);

  PCoordListProp = (CoordListRec *)
    calloc(NPoint,sizeof(CoordListRec));
  for (I=0;I<NPoint;I++) {
    MakeCoordListProp(&PStarCatProp[I],&PCoordListProp[I]);
  }
  Count = 0;
  Diff = 0.0;
  for (I=0;I<NMatch;I++) {
    for (J=0;J<NPoint;J++)
      for (K=J+1;K<NPoint;K++) {
	if ((PPWhMatch[J][I] >= 0)&&(PPWhMatch[K][I] >= 0)) {
	  WhJ = PPWhMatch[J][I];
	  WhK = PPWhMatch[K][I];

	  dX = (PCoordListProp[J].PX[WhJ] - PCoordListProp[K].PX[WhK]);
	  dY = (PCoordListProp[J].PY[WhJ] - PCoordListProp[K].PY[WhK]);
	  Diff += dX*dX+dY*dY;
	  Count++;
	}
      }
  }
  if (NMatch > 0) {
    for (I=0;I<NPoint;I++) {
      free(PPWhMatch[I]);
    }
    free(PPWhMatch);
  }
  for (I=0;I<NPoint;I++)
    DestructCoordListProp(&PCoordListProp[I]);
  free(PCoordListProp);
  Err = sqrt(Diff/Count);
  if (When==0)
    fprintf(stdout,"Start: RMS Error = %.3f pixels\n",Err);
  else 
    fprintf(stdout,"End: RMS Error = %.3f pixels\n",Err);

  return Err;
}

extern void ImproveSolution(int NPoint, StarCatRec *PStarCatProp,
			    int *perr, SuperAlignParmRec *PSAPP)
/* Starting with a number (NPoint) of exposures (stored in the array
   PStarCatProp) which are already aligned with each other.  This
   routine simultaneously varies the x and y-shifts and rotation
   angles for (NPoint-1) of these exposures.  Since one of these 
   x,y shifts and rotation angles is arbitrary, cannot vary all of
   them.  Since the solution should already be close, this routine 
   works by linearizing the computation and then using that to 
   compute the best-fit parameters.
  
   *perr == 1 if some error occurs in the computation. */
{
  int NMatch,**PPWhMatch,I,J,K,Dim;
  float **A,*B,**M,Sum,*p,*X,OErr,FErr;
  CoordListRec *PCoordListProp;

  PCoordListProp = (CoordListRec *)
    calloc(NPoint,sizeof(CoordListRec));
  for (I=0;I<NPoint;I++) {
    MakeCoordListProp(&PStarCatProp[I],&PCoordListProp[I]);
  }

  fprintf(stdout,"Improving overall solution.\n");

  PPWhMatch = (int **)calloc(NPoint,sizeof(int *));
  GStackStarCatProp(NPoint,PStarCatProp,2,PSAPP,&NMatch,PPWhMatch);

  OErr = CalcCurrentError(NPoint,PStarCatProp,0,PSAPP);

  Dim = (NPoint-1)*3;
  allocFloatArray(&A,Dim,Dim);
  B = (float *)calloc(Dim,sizeof(float));
  for (I=0;I<NMatch;I++) {
    for (J=0;J<NPoint;J++)
      for (K=J+1;K<NPoint;K++) {
	if ((PPWhMatch[J][I] >= 0)&&(PPWhMatch[K][I] >= 0)) {
	  int WhJ,WhK;
	  float RelRotJ,RelRotK,cosJ,sinJ,cosK,sinK;
	  float A0,B0,C0,D0,E0,F0;
	  int J0,K0;

	  WhJ = PPWhMatch[J][I];
	  WhK = PPWhMatch[K][I];
	  RelRotJ = PStarCatProp[J].HDR_ROT;
	  RelRotK = PStarCatProp[K].HDR_ROT;
	  cosJ = cos(RelRotJ*M_PI/180);
	  sinJ = sin(RelRotJ*M_PI/180);
	  cosK = cos(RelRotK*M_PI/180);
	  sinK = sin(RelRotK*M_PI/180);
	  A0 = -PStarCatProp[J].PX[WhJ]*sinJ + PStarCatProp[J].PY[WhJ]*cosJ;
	  B0 = -PStarCatProp[K].PX[WhK]*sinK + PStarCatProp[K].PY[WhK]*cosK;
	  C0 = PStarCatProp[J].HDR_DX-PStarCatProp[K].HDR_DX +
	    (PStarCatProp[J].PX[WhJ]*cosJ + PStarCatProp[J].PY[WhJ]*sinJ) -
	    (PStarCatProp[K].PX[WhK]*cosK + PStarCatProp[K].PY[WhK]*sinK);
	  D0 = -PStarCatProp[J].PX[WhJ]*cosJ -PStarCatProp[J].PY[WhJ]*sinJ;
	  E0 = -PStarCatProp[K].PX[WhK]*cosK -PStarCatProp[K].PY[WhK]*sinK;
	  F0 = PStarCatProp[J].HDR_DY-PStarCatProp[K].HDR_DY +
	    (-PStarCatProp[J].PX[WhJ]*sinJ + PStarCatProp[J].PY[WhJ]*cosJ) -
	    (-PStarCatProp[K].PX[WhK]*sinK + PStarCatProp[K].PY[WhK]*cosK);
	  J0 = J-1;
	  K0 = K-1;
	  if (J0>=0) {
	    A[J0*3][J0*3] += (A0*A0 + D0*D0)/1e4;
	    A[J0*3][J0*3+1] += A0/1e2;
	    A[J0*3][J0*3+2] += D0/1e2;
	    A[J0*3][K0*3] += (-A0*B0 - D0*E0)/1e4;
	    A[J0*3][K0*3+1] += -A0/1e2;
	    A[J0*3][K0*3+2] += -D0/1e2;
	    B[J0*3] += (-A0*C0-D0*F0)/1e2;
	    A[J0*3+1][J0*3] += A0/1e2;
	    A[J0*3+1][J0*3+1] += 1;
	    A[J0*3+1][K0*3] += -B0/1e2;
	    A[J0*3+1][K0*3+1] += -1.0;
	    B[J0*3+1] += -C0;
	    A[J0*3+2][J0*3] += D0/1e2;
	    A[J0*3+2][J0*3+2] += 1;
	    A[J0*3+2][K0*3] += -E0/1e2;
	    A[J0*3+2][K0*3+2] += -1.0;
	    B[J0*3+2] += - F0;
	  }
	  if (J0>=0) {
	    A[K0*3][J0*3] += (-A0*B0 - D0*E0)/1e4;
	    A[K0*3][J0*3+1] += -B0/1e2;
	    A[K0*3][J0*3+2] += -E0/1e2;
	  }
	  A[K0*3][K0*3] += (B0*B0 + E0*E0)/1e4;
	  A[K0*3][K0*3+1] += B0/1e2;
	  A[K0*3][K0*3+2] += E0/1e2;
	  B[K0*3] += (B0*C0+E0*F0)/1e2;
	  if (J0>=0) {
	    A[K0*3+1][J0*3] += -A0/1e2;
	    A[K0*3+1][J0*3+1] += -1;
	  }
	  A[K0*3+1][K0*3] += B0/1e2;
	  A[K0*3+1][K0*3+1] += 1.0;
	  B[K0*3+1] += C0;
	  if (J0>=0) {
	    A[K0*3+2][J0*3] += -D0/1e2;
	    A[K0*3+2][J0*3+2] += -1;
	  }
	  A[K0*3+2][K0*3] += E0/1e2;
	  A[K0*3+2][K0*3+2] += 1.0;
	  B[K0*3+2] += F0;
	}
      }
  }
  if (NMatch > 0) {
    for (I=0;I<NPoint;I++) {
      free(PPWhMatch[I]);
    }
    free(PPWhMatch);
  }
  p = (float *)calloc(Dim,sizeof(float));
  X = (float *)calloc(Dim,sizeof(float));
  chdecomp(A,Dim,p,perr);
  if (*perr)
    return;
  chsvx(A,Dim,p,B,X);
  for (J=1;J<NPoint;J++) {
    PStarCatProp[J].HDR_ROT += (180/M_PI)*X[(J-1)*3]/1e2;
    PStarCatProp[J].HDR_DX += X[(J-1)*3+1];
    PStarCatProp[J].HDR_DY += X[(J-1)*3+2];
  }
  free(p);
  free(X);
  free(B);
  freeFloatArray(A,Dim,Dim);
  FErr = CalcCurrentError(NPoint,PStarCatProp,1,PSAPP);
  if (FErr < 0.99 * OErr)
    *perr = 1;
  else
    *perr = 0;
}

extern void ImproveSolutionMP(int NPoint, StarCatRec *PStarCatProp,
			      SuperAlignParmRec *PSAPP)
{
  int err;
  
  err = 1;
  while (err) {
    ImproveSolution(NPoint,PStarCatProp,&err,PSAPP);
  }
}

extern void ImproveSolutionAA(int NPoint, StarCatRec *PStarCatProp,
			      StarCatRec *PStackStarCatProp,
			      SuperAlignParmRec *PSAPP)
{
  int I,OK;
  float OErr,NErr;

  OErr = CalcCurrentError(NPoint,PStarCatProp,0,PSAPP);
  for (I=0;I<NPoint;I++) {
    AttemptAlign(PStackStarCatProp,&PStarCatProp[I],&OK,PSAPP,1);
  }
  fprintf(stdout,"Start: RMS Error = %.3f pixels\n",OErr);
  NErr = CalcCurrentError(NPoint,PStarCatProp,1,PSAPP);
}

extern void AlignEntirePointing(PointingRec *PPointingProp, int Point,
				SuperAlignParmRec *PSAPP, int *POK, int Iter)
/* 
   The variable Point does not affect operation except appearing in
   the output.  */
{
  StarCatRec *PNStarCatProp;
  char FN[150];
  int NAlign,Wh1,Wh2,*PUsed,I,Try,TotalIter,Stand;
  long Global;

  Global = -5-Iter;
  PNStarCatProp = (StarCatRec *)
    calloc(PPointingProp->NExp,sizeof(StarCatRec));
  PUsed = (int *)calloc(PPointingProp->NExp,sizeof(int));
  NAlign = 0;

  Try = 0;
  TotalIter = 0;
  while (NAlign < PPointingProp->NExp) {
    StarCatRec StackStarCatProp;
    int Found,OK;

    if (NAlign > 0) {
      sprintf(FN,"CRREJ-STARSTACK-%i",Point);
      MakeStackStarCatProp(NAlign,PNStarCatProp,
			   &StackStarCatProp,FN,0,PSAPP);
    }
    OK = 0;
    while (!OK) {
      if (NAlign == 0) {
	float TR;
	
	TR = ran3(&Global);
	Wh1 = (int)(PPointingProp->NExp*TR+Try)%PPointingProp->NExp;
	CopyStarCatProp(&StackStarCatProp,
			&PPointingProp->PStarCatProp[Wh1]);
      }
      Found = 0;
      while (!Found) {
	Try++;
        Wh2 = (int)(PPointingProp->NExp*ran3(&Global)+Try)
	  %PPointingProp->NExp;
        if ((!PUsed[Wh2])&&(Wh2!=Wh1))
	  Found = 1;
      }
      Stand = 1;
      if (TotalIter > 50)
	Stand = 0;
      if (TotalIter > 80)
	Stand = -1;
      AttemptAlign(&StackStarCatProp,
		   &PPointingProp->PStarCatProp[Wh2],&OK,PSAPP,Stand);
      TotalIter++;
      if (TotalIter > PPointingProp->NExp*3+150) {
	*POK = 0;
	return;
      }      
    }
    if (OK) {
      if (NAlign == 0) {
	CopyStarCatProp(&PNStarCatProp[NAlign],
			&PPointingProp->PStarCatProp[Wh1]);
	PUsed[Wh1] = 1;
	NAlign++;
      }

      CopyStarCatProp(&PNStarCatProp[NAlign],
		      &PPointingProp->PStarCatProp[Wh2]);
      PUsed[Wh2] = 1;
      NAlign++;
    }
    DestructStarCatProp(&StackStarCatProp);
  }
  /*  ImproveSolutionMP(NAlign,PNStarCatProp);*/
  for (I=0;I<PPointingProp->NExp;I++) {
    DestructStarCatProp(&PPointingProp->PStarCatProp[I]);
    CopyStarCatProp(&PPointingProp->PStarCatProp[I],
		    &PNStarCatProp[I]);
    DestructStarCatProp(&PNStarCatProp[I]);
  }
  free(PNStarCatProp);
  free(PUsed);
  *POK = 1;
}

extern void GetCOM(StarCatRec *PStarCatProp, float *PCX, float *PCY)
/* Get "center-of-mass" for single exposure */
{
  CoordListRec CoordListProp;
  int I;
  float Sum;

  MakeCoordListProp(PStarCatProp,&CoordListProp);
  *PCX = 0;
  *PCY = 0;
  Sum = 0;
  for (I=0;I<CoordListProp.NObj;I++) {
    *PCX += CoordListProp.PX[I];
    *PCY += CoordListProp.PY[I];
    Sum++;
  }
  *PCX /= (Sum+1e-5);
  *PCY /= (Sum+1e-5);
}

extern void FindCOM(int NPointTBD, StarCatRec *PStarCatPropTBD,
		    float *PCX, float *PCY)
/* Find center-of-mass for group of exposures.  */
{
  float SumX,SumY;
  int I;

  SumX = 0;
  SumY = 0;
  for (I=0;I<NPointTBD;I++) {
    StarCatRec *PStarCatProp;
    float TX,TY;
    
    PStarCatProp = &PStarCatPropTBD[I];
    GetCOM(PStarCatProp,&TX,&TY);
    SumX += TX;
    SumY += TY;
  }
  *PCX = SumX / NPointTBD;
  *PCY = SumY / NPointTBD;
}

extern void FindClosestPoint(float CX, float CY,
			     int NPointTBD, StarCatRec *PStarCatPropTBD,
			     int *PWh, int Attempts)
{
  float DiffX,DiffY,Diff,MinDiff;
  int I,Best;
  long Global;

  Global = -2-Attempts;
  for (I=0;I<NPointTBD;I++) {
    StarCatRec *PStarCatProp;
    float TX,TY;
    
    PStarCatProp = &PStarCatPropTBD[I];
    GetCOM(PStarCatProp,&TX,&TY);
    DiffX = TX - CX;
    DiffY = TY - CY;
    Diff = sqrt(DiffX*DiffX + DiffY*DiffY);    
    if (Attempts>12) Diff += ran3(&Global)*500;
    else if (Attempts>8) Diff += ran3(&Global)*200;
    else if (Attempts>5) Diff += ran3(&Global)*100;
    else if (Attempts>0) Diff += ran3(&Global)*50;

    if ((I==0)||(Diff < MinDiff)) {
      MinDiff = Diff;
      Best = I;
    }
  }
  *PWh = Best;
}

extern void RemoveStarCatProp(StarCatRec *PStarCatProp, int Wh,
			      int *PNum)
{
  int I;

  DestructStarCatProp(&PStarCatProp[Wh]);
  for (I=Wh;I<*PNum-1;I++)
    CopyStarCatProp(&PStarCatProp[I],&PStarCatProp[I+1]);
  (*PNum)--;
}
			
extern void FixPointings(int NPoint, PointingRec *PPointingProp0,
			 StarCatRec *PStarCatPropN)
/* PPointingProp0 is the array of all the Pointings/Groups
   StarCatPropN is the array of all stacked star lists incorporated
   into the final alignment solution */
{
  int I,J;
  
  /* Loop over pointings/groups and apply x and y-shifts, rotations
     determined for the group to individual exposures within the
     group */

  for (I=0;I<NPoint;I++) {
    PointingRec *PPointingProp;
    int Wh;
    char SStr[100];
    
    sprintf(SStr,"POINT-%i",I);
    for (J=0;J<NPoint;J++)
      if (!strcmp(PStarCatPropN[J].FN,SStr))
	Wh = J;
    
    PPointingProp = &PPointingProp0[I];
    for (J=0;J<PPointingProp->NExp;J++) {
      StarCatRec *PStarCatProp;
      float NewRot,NewdX,NewdY,cosf,sinf,NX,NY;      
      
      PStarCatProp = &PPointingProp->PStarCatProp[J];
      /* Since PStarCatProp is only given relative to its own frame
         (PStarCatPropN[Wh]), it needs to be converted to the universal
         coordinate system.  */         
      NewRot = PStarCatPropN[Wh].HDR_ROT + PStarCatProp->HDR_ROT;
      cosf = cos(PStarCatPropN[Wh].HDR_ROT*M_PI/180);
      sinf = sin(PStarCatPropN[Wh].HDR_ROT*M_PI/180);
      NX = PStarCatProp->HDR_DX*cosf + PStarCatProp->HDR_DY*sinf;
      NY = -PStarCatProp->HDR_DX*sinf + PStarCatProp->HDR_DY*cosf;
      PStarCatProp->HDR_DX = NX + PStarCatPropN[Wh].HDR_DX;
      PStarCatProp->HDR_DY = NY + PStarCatPropN[Wh].HDR_DY;
      PStarCatProp->HDR_ROT = NewRot;
    }
  }
}

extern void OutputPointings(int NPoint, PointingRec *PPointingProp0,
			    StarCatRec *PStarCatPropN2)
/* Put all the StarCatProp from all PPointingProp0 into a single
   array *PStarCatPropN2 */
{
  int I,J,Tot;
  
  Tot = 0;
  for (I=0;I<NPoint;I++) {
    PointingRec *PPointingProp;
    
    PPointingProp = &PPointingProp0[I];
    for (J=0;J<PPointingProp->NExp;J++) {
      StarCatRec *PStarCatProp,*PStarCatPropN;

      PStarCatProp = &PPointingProp->PStarCatProp[J];
      PStarCatPropN = &PStarCatPropN2[Tot];
      CopyStarCatProp(PStarCatPropN,PStarCatProp);
      Tot++;
    }
  }
  fprintf(stderr,"%i Exposures Total.\n",Tot);
}

extern void ZeroPointings(StarCatRec *PStarCatProp,
			  AllExposureRec *PAllExposureProp,
			  StarCatRec *PStarCatPropGroup, int Num)
/* Shift and Rotate pointings such that the reference image,
   e.g., the first image in the exposure list, has shifts of 0.0,0.0 and
   a rotation angle of 0. */
{
  int I,Wh;
  StarCatRec StarCatProp;
  float ShiftX,ShiftY,RotateAng;

  for (I=0;I<PAllExposureProp->NExp;I++) {
    if (!strcmp(PStarCatProp[I].FN,PAllExposureProp->PExposureProp[0].FN)) {
      Wh = I;
    }
  }
  CopyStarCatProp(&StarCatProp,&PStarCatProp[Wh]);
  for (I=0;I<PAllExposureProp->NExp;I++) {
    DetermineGuess(&StarCatProp,&PStarCatProp[I],&ShiftX,&ShiftY,&RotateAng);
    PStarCatProp[I].HDR_DX = ShiftX;
    PStarCatProp[I].HDR_DY = ShiftY;
    PStarCatProp[I].HDR_ROT = RotateAng;
  }
  for (I=0;I<Num;I++) {
    DetermineGuess(&StarCatProp,&PStarCatPropGroup[I],
		   &ShiftX,&ShiftY,&RotateAng);
    PStarCatPropGroup[I].HDR_DX = ShiftX;
    PStarCatPropGroup[I].HDR_DY = ShiftY;
    PStarCatPropGroup[I].HDR_ROT = RotateAng;
  }
  DestructStarCatProp(&StarCatProp);
}

extern void OutputOffsets(int NExp, StarCatRec *PStarCatPropN,
			  char *FixedPointFN)
/* Output all offsets from PStarCatPropN files to file with the
   FixedPointFN. */
{
  int I;
  FILE *fout;
  
  fout = fopen(FixedPointFN,"w");
  
  for (I=0;I<NExp;I++) {
    StarCatRec *PStarCatProp;
    
    PStarCatProp = &PStarCatPropN[I];
    fprintf(fout,"%s %.2f %.2f %.3f\n",PStarCatProp->FN,
	    PStarCatProp->HDR_DX,PStarCatProp->HDR_DY,
	    PStarCatProp->HDR_ROT);
  }
  fclose(fout);
}

extern void CentralCommand(AllExposureRec *PAllExposureProp,
			   char *TotStarListFN,
			   char *FixedPointFN,
			   SuperAlignParmRec *PSAPP,
			   float GroupSize, int FirstRefImage)
{
  int NPoint,I,*PUsed,J;
  PointingRec *PPointingProp0;
  StarCatRec *PStarCatPropO;
  int NumO,NumN,Wh,err;
  StarCatRec *PStarCatPropN,*PStarCatPropN2;
  StarCatRec StackStarCatProp;
  float CX,CY;
  int OK,Abort;
  
  /* Break up all the original exposures into separate pointings/groups */

  GroupExposures(&NPoint,&PPointingProp0,PAllExposureProp,GroupSize,FirstRefImage);

  PStarCatPropO = (StarCatRec *)calloc(NPoint,sizeof(StarCatRec));
  PStarCatPropN = (StarCatRec *)calloc(NPoint,sizeof(StarCatRec));
 
  /* Loop over the different pointings and work out the alignments of
     all the original exposures */

  for (I=0;I<NPoint;I++) {
    PointingRec *PPointingProp;
    char FN[150];

    PPointingProp = &PPointingProp0[I];
    fprintf(stdout,"Aligning Pointing Set #%i of %i...\n",I,
	    NPoint);
    fprintf(stdout,"%i exposures within Group.\n",PPointingProp->NExp);
    for (J=0;J<PPointingProp->NExp;J++) {
      fprintf(stdout,"%s\n",PPointingProp->PStarCatProp[J].FN);
    }
    if (PPointingProp->NExp>1) {
      int Iter;

      Iter = 0;
      OK = 0;
      while (!OK) {
	AlignEntirePointing(PPointingProp,I,PSAPP,&OK,Iter);
	Iter++;
	if ((!OK)&&(Iter > 0)) {
	  fprintf(stderr,"Failed to align exposures within group!");
	  exit(1);	
	}
      }
    }
    sprintf(FN,"POINT-%i",I);
    MakeStackStarCatProp(PPointingProp->NExp,
			 PPointingProp->PStarCatProp,
			 &PStarCatPropO[I],FN,1,PSAPP);
  }
  fprintf(stdout,"Done Aligning Images at Individual Pointings...\n");
  fprintf(stdout,"Now to align them relative to each other...\n\n");

  /* Find center-of-mass for all the pointings */
  FindCOM(NPoint,PStarCatPropO,&CX,&CY);

  /* Find pointing which is closest to the center-of-mass */
  FindClosestPoint(CX,CY,NPoint,PStarCatPropO,&Wh,0);

  if (FirstRefImage)
    Wh = 0;

  fprintf(stdout,"Centering on %s.\n",PStarCatPropO[Wh].FN);
  NumO = NPoint;
  CopyStarCatProp(&PStarCatPropN[0],&PStarCatPropO[Wh]);
  RemoveStarCatProp(PStarCatPropO,Wh,&NumO);
  NumN = 1;

  /* Continue to accrete pointings into mosaic until no pointings left */
  Abort = 0;
  while (NumO > 0) {
    StarCatRec StackStarCatProp;
    int OK,Attempts;

    /* Find New Center-of-Mass */
    FindCOM(NumN,PStarCatPropN,&CX,&CY);
    /* Stack the star lists from all the pointings which have been
       incorporated into the mosaic */
    MergeStarCatProp(NumN,PStarCatPropN,&StackStarCatProp,
		     "ALLSTARS",PSAPP);
    Attempts = 0;

    OK = 0;
    while (!OK) {
      int Stand;
      /* Find Pointing/Group which is closest to this Center-of-Mass */
      FindClosestPoint(CX,CY,NumO,PStarCatPropO,&Wh,Attempts);
      fprintf(stdout,"Next to %s.\n",PStarCatPropO[Wh].FN);

      /* Attempt to align new pointing relative to stacked star list. */
      Stand = 0;
      if (Attempts > 25)
	Stand = -1;
      AttemptAlign(&StackStarCatProp,&PStarCatPropO[Wh],&OK,PSAPP,Stand);
      if (!OK) {
	fprintf(stderr,"Alignment onto central mosaic failed!\n");
	Attempts++;
	if (Attempts > 50) {
	  fprintf(stderr,"Alignment failed.\n");
	  Abort = 1;
	  OK = 1;
	}
      }
    }

    if (!Abort) {
      /* Add star list from current pointing to list of accreted star lists */
      CopyStarCatProp(&PStarCatPropN[NumN],&PStarCatPropO[Wh]);
      NumN++;
      RemoveStarCatProp(PStarCatPropO,Wh,&NumO);
      /* Attempt to iteratively improve the alignment among all pointings */
      ImproveSolutionMP(NumN,PStarCatPropN,PSAPP);
      DestructStarCatProp(&StackStarCatProp);
    }
    else {
      fprintf(stderr,"Adding Remaining Pointings to the Mosaic (without alignment)!\n");
      while (NumO > 0) {
	FindClosestPoint(CX,CY,NumO,PStarCatPropO,&Wh,Attempts);
	CopyStarCatProp(&PStarCatPropN[NumN],&PStarCatPropO[Wh]);
	NumN++;
	RemoveStarCatProp(PStarCatPropO,Wh,&NumO);
      }
    }
  }
  FixPointings(NPoint,PPointingProp0,PStarCatPropN);

  /* Put all the object lists in a single array */
  PStarCatPropN2 = (StarCatRec *)calloc(PAllExposureProp->NExp,
					sizeof(StarCatRec));
  OutputPointings(NPoint,PPointingProp0,PStarCatPropN2);
  /* Adjust the shifts/rotations of all exposures such that the reference
     image is at 0,0 and no rotation. */
  ZeroPointings(PStarCatPropN2,PAllExposureProp,
		PStarCatPropN,NumN); 
  /* Stack all the exposures to derive a object list for the whole
     mosaic */
  MergeStarCatProp(NumN,PStarCatPropN,&StackStarCatProp,"ALLSTARS",PSAPP);
  OutputOffsets(PAllExposureProp->NExp,PStarCatPropN2,FixedPointFN);

  /* Output a stacked version of all previous exposures */
  OutputStarCatProp(&StackStarCatProp,TotStarListFN);  
}

extern void InvertTrans(float dx, float dy, float dth,
			float *pdx, float *pdy, float *pdth)
{
  float sinf,cosf;

  sinf = sin(dth*M_PI/180);
  cosf = cos(dth*M_PI/180);
  *pdx = dy*sinf - dx*cosf;
  *pdy = -dx*sinf - dy*cosf;
  *pdth = -dth;
}

extern void NoInvertTrans(float dx, float dy, float dth,
			  float *pdx, float *pdy, float *pdth)
{
  float sinf,cosf;

  *pdx = dx;
  *pdy = dy;
  *pdth = dth;
}

extern void ReadMatchinList(char *FN, AllExposureRec *PAllExposureProp,
			    int *PFirstRefImage)
/* Read in Exposure List and other specifications from FN */
{
  FILE *fin;
  float dx,dy,rot;
  char FirstRefImage[50];
  int I;

  fin = fopen(FN,"r");
  fscanf(fin,"%i %i%[^\n]",&PAllExposureProp->NExp,PFirstRefImage);
  fgetc(fin);
  PAllExposureProp->PExposureProp = (ExposureRec *)
    calloc(PAllExposureProp->NExp,sizeof(ExposureRec));
  for (I=0;I<PAllExposureProp->NExp;I++) {
    ExposureRec *PExposureProp;

    PExposureProp = &PAllExposureProp->PExposureProp[I];
    fscanf(fin,"%s %g %g %g[^\n]",PExposureProp->FN,&dx,&dy,&rot);
    NoInvertTrans(dx,dy,rot,&PExposureProp->hdr_dx,&PExposureProp->hdr_dy,
		  &PExposureProp->hdr_rot);
    fgetc(fin);
  }
  fclose(fin);
  for (I=0;I<PAllExposureProp->NExp;I++) {
    ExposureRec *PExposureProp;

    PExposureProp = &PAllExposureProp->PExposureProp[I];
    ReadMatchinFiles(PExposureProp->FN,&PExposureProp->StarCatProp,
		     PExposureProp->hdr_dx,PExposureProp->hdr_dy,
		     PExposureProp->hdr_rot);
  }
}

int main(int argc, char *argv[])
{
  AllExposureRec AllExposureProp;
  SuperAlignParmRec SuperAlignParmProp;
  int FirstRefImage;

  /* Parameters described in superalign.h */
  SuperAlignParmProp.TypDist = 0.5/20;
  SuperAlignParmProp.MaxDist = 4.0/20;
  SuperAlignParmProp.TolPair = 8.0/20;
  SuperAlignParmProp.MaxShiftErr = 200.0/20;
  SuperAlignParmProp.MaxRotErr = 2.0;
  SuperAlignParmProp.GroupSize = 200.0/20;

  ReadMatchinList(argv[1],&AllExposureProp,&FirstRefImage);
  CentralCommand(&AllExposureProp,argv[2],argv[3],
		 &SuperAlignParmProp,SuperAlignParmProp.GroupSize,
		 FirstRefImage);
}
