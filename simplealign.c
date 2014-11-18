
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "superalign.h"

void hpsort(long n, float *ra, int *ind);

typedef enum {
  TakeBound, Warning
} RespType;

extern int HighSigPowerOf2(int N)
/* Determine smallest power of 2 which is just greater or equal to N */
{
  int N2;
  
  if (N > 1) {
    N2 = 1;
    while (N2 < N) {
      N2 *= 2;
    }
    N2 /= 2;
  }
  else N2 = 1;

  return N2;
}

extern double BinSearchF(float *PArr, double X, int Steps, 
			 RespType Response)
{
  /*Response is TakeBound if just take boundary
            is warning if stop
   PArr should go from 0..Steps*/
  int C2, Ind;
  int D, Increasing;

  if ((X >= PArr[0] && X <= PArr[Steps]) ||
      (X <= PArr[0] && X >= PArr[Steps])) {
    C2 = HighSigPowerOf2(Steps+1);
    Ind = C2 / 2;
    D = 0;
    Increasing = (PArr[Steps] > PArr[0]);
    while (!D && ((Increasing && ((PArr[C2] < X) || (PArr[C2 - 1] > X))) ||
		  (!Increasing && ((PArr[C2] > X) || (PArr[C2 - 1] < X))))) {
      if (((PArr[C2] > X) && Increasing) || ((PArr[C2] < X) && !Increasing))
	C2 -= Ind;
      else
	C2 += Ind;
      if (C2 > Steps)
	C2 -= Ind;
      Ind /= 2;
      if (Ind == 0)
	D = 1;
    }
    if (C2 == 0)
      return 0.0;
    else if (PArr[C2] != PArr[C2-1])
      return (C2 - 1 + (X - PArr[C2 - 1]) / (PArr[C2] - PArr[C2 - 1]));
    else 
      return C2 - 1 + 0.5;
  } else if (Response == Warning) {
    fprintf(stdout,"Outside boundary in BinSearch\n");
	fprintf(stdout,"%f %f %f\n",X,PArr[0],PArr[Steps]);
    exit(1);
  } else if (((X < PArr[0])&&(PArr[1]>PArr[0]))||
	     ((X > PArr[0])&&(PArr[1]<PArr[0])))
    return 0.0;
  else
    return Steps;
  return 0;
}

extern void NeighborList(float *X1, float *Y1, int N1,
			 float *X2, float *Y2, int N2,
			 float Tol, int *NNeigh,
			 int **NeighList);

extern float CalcSigma(float dX, float dY, float Orient,
		       float *X1, float *Y1, float *M1, int NL1,
		       float *X2, float *Y2, float *M2, int NL2,
		       int *NNeigh, int **NeighList, 
		       SuperAlignParmRec *PSAPP,
		       int Flag, int *POK1, int *POK2, float *PNum)
/* This routine provides a score to a possible alignment between
   two object lists.  dX and dY denote the x/y shifts which should
   be applied to the object positions and Orient denotes the rotation
   to be applied.

   The x,y positions of the objects in list 1 are given by X1,Y1 and
   those in list 2 are given by X2,Y2.  The scalefree fluxes of objects
   in these lists are given by M1 and M2 respectively.  N1 and N2 denote
   the total number of objects in those lists.  

   The scoring system for matches is different here from the clipped
   chi-squared system, which gives huge weight to outliers.  This can be
   problematic since some objects are not real or have centroids which
   are not easy to define (because of different blending with neighbors).
   Here is the scoring system is a positive one which gives huge
   rewards for each positive match.  The scoring for each object is 
   a sharply increasing function of decreasing distance until the 
   distance match becomes typical for the quality of data or instrument.

   The weighting is:
                            1                         1
              ----------------------------- x  -----------------------------
              TypDist * TypDist + Dist*Dist    (1/flux**0.33 + 1/flux2**0.33)

   Set Flag == 1, if you want objects which appear as matches 
   to be checked in the POK1 and POK2 arrays.  

   NNeigh and NeighList are the neighbor lists for individual objects in
   the first list.

   *PNum is an approximate estimate on the number of objects which 
   are important for the current score.  

   This routine returns a score rating the alignment.  */
{
  int I,J;
  float TDX,TDY,Dist,ScoreA,Score,ScoreSq,CosF,SinF;

  CosF = cos(M_PI*Orient/180);
  SinF = sin(M_PI*Orient/180);
  Score = 0;
  ScoreSq = 0;
  if (Flag) {
    for (I=0;I<NL1;I++) 
      POK1[I] = 0;
    for (I=0;I<NL2;I++) 
      POK2[I] = 0;
  }
  for (I=0;I<NL1;I++) {    
    for (J=0;J<NNeigh[I];J++) {
      float NX2,NY2,OX2,OY2;

      OX2 = X2[NeighList[I][J]];
      OY2 = Y2[NeighList[I][J]];
      NX2 = dX+OX2*CosF+OY2*SinF;
      NY2 = dY-OX2*SinF+OY2*CosF;
      TDX = X1[I] - NX2;
      TDY = Y1[I] - NY2;
      Dist = sqrt(TDX*TDX+TDY*TDY);
      ScoreA = 1 / (PSAPP->TypDist*PSAPP->TypDist+Dist*Dist) / 
	(1/pow(M1[I]+1,0.3) + 1/pow(M2[NeighList[I][J]]+1,0.3));
      Score += ScoreA;
      ScoreSq += ScoreA*ScoreA;
      if (Dist < PSAPP->MaxDist) {
	if (Flag) {
	  POK1[I] = 1;
	  POK2[NeighList[I][J]] = 1;
	}
      }
    }
  }
  *PNum = Score * Score / (ScoreSq+1e-15);
  if (*PNum < 3)
    Score *= exp(-(3-*PNum)*2);
  return Score;
}

#define NFPRINTF(w,x) {fprintf(w, "\33[1M> %s\n\33[1A",x);}

extern void FindSolution(float *dX, float *dY, float *Orient,
			 float *X1, float *Y1, float *M1, int NL1,
			 float *X2, float *Y2, float *M2, int NL2,
			 int *NNeigh, int **NeighList,
			 SuperAlignParmRec *PSAPP)
/* Tune up the best fit *dX, *dY, and *Orient that matches
   the (X1,Y1,M1) array with the (X2,Y2,M2) array.  Since the solution
   should already be very close, no need to compute the number of neighbors
   or neighbor lists for each object -- contained in NNeigh and NNeighList */
{
  float Score,dX1,dY1,Orient1,D,BestScore,NumC,BestNumC;
  int I;
  long Global2;

  Score = CalcSigma(*dX,*dY,*Orient,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
		    NNeigh,NeighList,PSAPP,0,NULL,NULL,&NumC);
  BestScore = Score;
  BestNumC = NumC;
  Global2 = -1;
  /* Use a simulated annealing type procedure to determine the
     solution with the highest score */
  for (I=0;I<1000;I++) {
    D = PSAPP->MaxDist*exp(-(float)I/200)/4;
    Orient1 = *Orient + (0.5-ran3(&Global2))*D/100;
    dX1 = *dX + (0.5-ran3(&Global2))*D;
    dY1 = *dY + (0.5-ran3(&Global2))*D;
    Score = CalcSigma(dX1,dY1,Orient1,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
		      NNeigh,NeighList,PSAPP,0,NULL,NULL,&NumC);
    if (Score > BestScore) {
      BestScore = Score;
      *dX = dX1;
      *dY = dY1;
      *Orient = Orient1;
    }
    else {
      Score = BestScore;
      NumC = BestNumC;
    }
  }
}

extern void RemoveClosePairs(float *X1, float *Y1, int NL1, 
			     int Where,
			     SuperAlignParmRec *PSAPP)
/* Removes objects which are closer than "SAPP->TolDist" pixels from each
   other.  In principle, this avoids having blending type issues 
   affect the alignment. */
{
  int I,J;
  int *NNeigh,**NeighList,*POK1;

  NNeigh = (int *)calloc(NL1,sizeof(int));
  NeighList = (int **)calloc(NL1,sizeof(int *));
  
  NeighborList(X1,Y1,NL1,X1,Y1,NL1,PSAPP->MaxDist,
	       NNeigh,NeighList);
  POK1 = (int *)calloc(NL1,sizeof(int));
  for (I=0;I<NL1;I++)
    POK1[I] = 0;
  for (I=0;I<NL1;I++)
    for (J=0;J<NNeigh[I];J++) 
      if (I!=NeighList[I][J]) {
	float dX,dY,Dist;
	
	dX = X1[I]-X1[NeighList[I][J]];
	dY = Y1[I]-Y1[NeighList[I][J]];
	Dist = dX*dX+dY*dY;
	if (Dist < PSAPP->TolPair) {
	  POK1[I] = 1;
	  POK1[NeighList[I][J]] = 1;
	}
      }
  for (I=0;I<NL1;I++)
    if (POK1[I]) {
      X1[I] = -1e5*Where;
      Y1[I] = -1e5*Where;
    }
  for (I=0;I<NL1;I++)
    free(NeighList[I]);
  free(NNeigh);
  free(NeighList);
  free(POK1);
}

extern void OutputXYM(float *X1, float *Y1, float *M1, int N1, char *FN)
{
  FILE *fout;
  int I;

  fout = fopen(FN,"w");
  for (I=0;I<N1;I++)
    fprintf(fout,"%g %g %g\n",X1[I],Y1[I],M1[I]);
  fclose(fout);
}

extern void DetermineGuess(StarCatRec *PStarCatBaseProp,
			   StarCatRec *PStarCatProp,
			   float *PShiftX, float *PShiftY,
			   float *PRotateAng)
/* Since we will be aligning *PStarCatProp relative to
   *PStarCatBaseProp, determine the shifts and rotation angles
   necessary to transform the objects in *PStarCatProp to
   match those in *PStarCatBaseProp */
{
  float ShiftX,ShiftY,RotateAng,TX,TY,cosf,sinf,NX,NY;

  *PRotateAng = PStarCatProp->HDR_ROT - PStarCatBaseProp->HDR_ROT;
  cosf = cos(PStarCatBaseProp->HDR_ROT*M_PI/180);
  sinf = sin(PStarCatBaseProp->HDR_ROT*M_PI/180);
  TX = PStarCatProp->HDR_DX - PStarCatBaseProp->HDR_DX;
  TY = PStarCatProp->HDR_DY - PStarCatBaseProp->HDR_DY;
  NX = TX*cosf - TY*sinf;
  NY = TX*sinf + TY*cosf;
  *PShiftX = NX;
  *PShiftY = NY;
}

extern void ShowGuess(StarCatRec *PStarCatBaseProp,
		      StarCatRec *PStarCatProp,
		      float ShiftdX, float ShiftdY, float RotateAng)
/* This routine converts the transformation variables from
   being relative to the *PStarCatBaseProp transformation to 
   being relative to the universal frame. */
{
  float cosf,sinf,NX,NY,HDR_DX,HDR_DY,HDR_ROT;

  HDR_ROT = RotateAng + PStarCatBaseProp->HDR_ROT;
  cosf = cos(PStarCatBaseProp->HDR_ROT*M_PI/180);
  sinf = sin(PStarCatBaseProp->HDR_ROT*M_PI/180);
  NX = ShiftdX*cosf + ShiftdY*sinf;
  NY = -ShiftdX*sinf + ShiftdY*cosf;
  HDR_DX = NX + PStarCatBaseProp->HDR_DX;
  HDR_DY = NY + PStarCatBaseProp->HDR_DY;
  if (HDR_ROT > 180) HDR_ROT -= 360.0;
  else if (HDR_ROT < -180) HDR_ROT += 360.0;

  fprintf(stdout,"Guess: dx = %.2f, dy = %.2f, theta = %.2f (Rel)\n",
	  HDR_DX,HDR_DY,HDR_ROT);
}

extern int TryAlignment(StarCatRec *PStarCatBaseProp,
			StarCatRec *PStarCatProp,
			SuperAlignParmRec *PSAPP,
			float *PBestdX, float *PBestdY,
			float *PBestRotAngle)
/* BestX > 0 ... means PStarCatProp is to the right of PStarCatBaseProp */
/* BestX < 0 ... means PStarCatProp is to the left of PStarCatBaseProp */
/* BestRotAngle > 0 means PStarCatProp is CCW of PStarCatBaseProp */
/* BestRotAngle < 0 means PStarCatProp is CW of PStarCatBaseProp */
{
  FILE *fout;
  int NL1,NL2,I,J;
  float CosF,SinF,ShiftX,ShiftY,RotateAng;
  float *X1,*Y1,*M1,*X2,*Y2,*M2,BestdX,BestdY,BestOrient,BestScore,Score;
  float *TX2,*TY2,NumC,YDiff,XDiff,TempX2,TempY2,BestNumC;
  int *POK1,*POK2,*ind,*ind2;
  int *NNeigh,**NeighList;
  float *PAng1,*PAng2,*PAngT,*PDist1,*PDist2,*PDistT,Ang;
  int K,K1,LL,HL,L,Tot,NTot;

  DetermineGuess(PStarCatBaseProp,PStarCatProp,
		 &ShiftX,&ShiftY,&RotateAng);
  Verbose = 1;
  
  fprintf(stdout,"Aligning %s (%i objs) to %s (%i objs)\n",
	  PStarCatProp->FN,PStarCatProp->NObj,
	  PStarCatBaseProp->FN,PStarCatBaseProp->NObj);
  ShowGuess(PStarCatBaseProp,PStarCatProp,ShiftX,ShiftY,RotateAng);
  NL1 = PStarCatBaseProp->NObj;
  X1 = (float *)calloc(NL1,sizeof(float));
  Y1 = (float *)calloc(NL1,sizeof(float));
  M1 = (float *)calloc(NL1,sizeof(float));
  POK1 = (int *)calloc(NL1,sizeof(int));
  for (I=0;I<NL1;I++) {
    X1[I] = PStarCatBaseProp->PX[I];
    Y1[I] = PStarCatBaseProp->PY[I];
    M1[I] = pow(10.0,-0.4*PStarCatBaseProp->PM[I]);
  }
  CosF = cos(M_PI*RotateAng/180);
  SinF = sin(M_PI*RotateAng/180);
  NL2 = PStarCatProp->NObj;
  TX2 = (float *)calloc(NL2,sizeof(float));
  TY2 = (float *)calloc(NL2,sizeof(float));
  X2 = (float *)calloc(NL2,sizeof(float));
  Y2 = (float *)calloc(NL2,sizeof(float));
  M2 = (float *)calloc(NL2,sizeof(float));
  POK2 = (int *)calloc(NL2,sizeof(int));
  for (I=0;I<NL2;I++) {
    TX2[I] = ShiftX+PStarCatProp->PX[I]*CosF+
      PStarCatProp->PY[I]*SinF;
    TY2[I] = ShiftY-PStarCatProp->PX[I]*SinF+
      PStarCatProp->PY[I]*CosF;
    X2[I] = PStarCatProp->PX[I];
    Y2[I] = PStarCatProp->PY[I];
    M2[I] = pow(10.0,-0.4*PStarCatProp->PM[I]);
  }

  OutputXYM(X1,Y1,M1,NL1,"A");
  OutputXYM(TX2,TY2,M2,NL2,"B");
  OutputXYM(X2,Y2,M2,NL2,"C");
  RemoveClosePairs(X1,Y1,NL1,1,PSAPP);
  RemoveClosePairs(TX2,TY2,NL2,-1,PSAPP);

  NNeigh = (int *)calloc(NL1,sizeof(int));
  NeighList = (int **)calloc(NL1,sizeof(int *));
  NeighborList(X1,Y1,NL1,TX2,TY2,NL2,PSAPP->MaxShiftErr,
	       NNeigh,NeighList);
  BestScore = 0;
  ind2 = (int *)calloc(NL2,sizeof(int));
  PAngT = (float *)calloc(NL2,sizeof(float));
  PDistT = (float *)calloc(NL2,sizeof(float));
  PAng2 = (float *)calloc(NL2,sizeof(float));
  PDist2 = (float *)calloc(NL2,sizeof(float));
  PAng1 = (float *)calloc(NL1,sizeof(float));
  PDist1 = (float *)calloc(NL1,sizeof(float));

  ind = (int *)calloc(NL1,sizeof(int));
  hpsort(NL1,M1,ind);
  BestdX = 0.0;
  BestdY = 0.0;
  BestOrient = 0.0;
  for (I=NL1-1;I>=0;I--) {
    int Wh;

    Wh = ind[I];
    if (!strcmp(PStarCatBaseProp->FN,
		PStarCatProp->FN)) {
      NNeigh[Wh] = 0;
      I = 0;
    }
    
    for (J=0;J<NNeigh[Wh];J++) {
      float Score,dX,dY;
      
      for (K=0;K<NL2;K++) {
	XDiff = X2[K] - X2[NeighList[Wh][J]];
	YDiff = Y2[K] - Y2[NeighList[Wh][J]];
	PAngT[K] = (180/M_PI)*atan2(YDiff,XDiff);
	PDistT[K] = sqrt(YDiff*YDiff+XDiff*XDiff);
      }
      hpsort(NL2,PAngT,ind2);
      for (K=0;K<NL2;K++) {
	PAng2[K] = PAngT[ind2[K]];
	PDist2[K] = PDistT[ind2[K]];
      }
      for (K=0;K<NL1;K++) {
	XDiff = X1[K]-X1[Wh];
	YDiff = Y1[K]-Y1[Wh];
	PAng1[K] = (180/M_PI)*atan2(YDiff,XDiff);
	PDist1[K] = sqrt(YDiff*YDiff+XDiff*XDiff);
	
	for (K1=-1;K1<2;K1++) {
	  LL = BinSearchF(PAng2,PAng1[K]+RotateAng-PSAPP->MaxRotErr+K1*360,
			  NL2-1,TakeBound);
	  HL = BinSearchF(PAng2,PAng1[K]+RotateAng+PSAPP->MaxRotErr+K1*360,
			  NL2-1,TakeBound)+1;
	  if (LL < 0) {
	    LL = BinSearchF(PAng2,PAng1[K]+RotateAng-PSAPP->MaxRotErr+K1*360,
			    NL2-1,TakeBound);
	    HL = BinSearchF(PAng2,PAng1[K]+RotateAng+PSAPP->MaxRotErr+K1*360,
			    NL2-1,TakeBound)+1;	    
	  }

	  for (L=LL;L<HL;L++) {
	    if ((fabs(PAng1[K]+RotateAng+K1*360-PAng2[L])<PSAPP->MaxRotErr)&&
		(fabs(PDist1[K] - PDist2[L]) < PSAPP->MaxDist)) {
	      Ang = PAng2[L]-PAng1[K];
	      CosF = cos(M_PI*Ang/180);
	      SinF = sin(M_PI*Ang/180);
	      TempX2 = X2[NeighList[Wh][J]]*CosF+Y2[NeighList[Wh][J]]*SinF;
	      TempY2 = -X2[NeighList[Wh][J]]*SinF+Y2[NeighList[Wh][J]]*CosF;
	      dX = X1[Wh]-TempX2;
	      dY = Y1[Wh]-TempY2;
	      Score = CalcSigma(dX,dY,Ang,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
				NNeigh,NeighList,PSAPP,
				0,NULL,NULL,&NumC);
	      if (Score > BestScore) {
		BestScore = CalcSigma(dX,dY,Ang,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
				      NNeigh,NeighList,PSAPP,
				      0,NULL,NULL,&NumC);
		BestdX = dX;
		BestdY = dY;
		BestOrient = Ang;
		BestNumC = NumC;
		if (NumC > 18) {
                  /* if more than 30 matches are obtained, this is good
                     enough and we will simply abort. */
		  I = 0;
		  J = NNeigh[Wh];
		  L = HL;
		  K1 = 2;
		  K = NL1;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  /* Attempt one last iterative improvement to best transformation
     obtained above */
  FindSolution(&BestdX,&BestdY,&BestOrient,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
	       NNeigh,NeighList,PSAPP);

  Score = CalcSigma(BestdX,BestdY,BestOrient,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
		    NNeigh,NeighList,PSAPP,1,POK1,POK2,&NumC);
  free(X1);
  free(Y1);
  free(M1);
  free(POK1);
  free(TX2);
  free(TY2);
  free(X2);
  free(Y2);
  free(M2);
  free(POK2);
  free(ind);
  free(ind2);
  free(PAngT);
  free(PDistT);
  free(PAng2);
  free(PDist2);
  free(PAng1);
  free(PDist1);
  free(NNeigh);
  for (I=0;I<NL1;I++)
    free(NeighList[I]);
  free(NeighList);

  *PBestdX = BestdX;
  *PBestdY = BestdY;
  *PBestRotAngle = BestOrient;
  fprintf(stdout,"Converged on %g good stars.\n",NumC);
  /*  fprintf(stdout,"Solution: dx = %.2f, dy = %.2f, theta = %.2f\n",
      BestdX,BestdY,BestOrient);*/
}

extern int ImproveAlignment(StarCatRec *PStarCatBaseProp,
			    StarCatRec *PStarCatProp,
			    SuperAlignParmRec *PSAPP,
			    float *PBestdX, float *PBestdY,
			    float *PBestRotAngle)
/* BestX > 0 ... means PStarCatProp is to the right of PStarCatBaseProp */
/* BestX < 0 ... means PStarCatProp is to the left of PStarCatBaseProp */
/* BestRotAngle > 0 means PStarCatProp is CCW of PStarCatBaseProp */
/* BestRotAngle < 0 means PStarCatProp is CW of PStarCatBaseProp */
{
  FILE *fout;
  int NL1,NL2,I,J;
  float CosF,SinF,ShiftX,ShiftY,RotateAng;
  float *X1,*Y1,*M1,*X2,*Y2,*M2,BestdX,BestdY,BestOrient,BestScore,Score;
  float *TX2,*TY2,NumC,YDiff,XDiff,TempX2,TempY2,BestNumC;
  int *POK1,*POK2,*ind,*ind2;
  int *NNeigh,**NeighList;
  float *PAng1,*PAng2,*PAngT,*PDist1,*PDist2,*PDistT,Ang;

  DetermineGuess(PStarCatBaseProp,PStarCatProp,
		 &ShiftX,&ShiftY,&RotateAng);
  Verbose = 1;
  
  fprintf(stdout,"Aligning %s (%i objs) to %s (%i objs)\n",
	  PStarCatProp->FN,PStarCatProp->NObj,
	  PStarCatBaseProp->FN,PStarCatBaseProp->NObj);
  ShowGuess(PStarCatBaseProp,PStarCatProp,ShiftX,ShiftY,RotateAng);
  NL1 = PStarCatBaseProp->NObj;
  X1 = (float *)calloc(NL1,sizeof(float));
  Y1 = (float *)calloc(NL1,sizeof(float));
  M1 = (float *)calloc(NL1,sizeof(float));
  POK1 = (int *)calloc(NL1,sizeof(int));
  for (I=0;I<NL1;I++) {
    X1[I] = PStarCatBaseProp->PX[I];
    Y1[I] = PStarCatBaseProp->PY[I];
    M1[I] = pow(10.0,-0.4*PStarCatBaseProp->PM[I]);
  }
  CosF = cos(M_PI*RotateAng/180);
  SinF = sin(M_PI*RotateAng/180);
  NL2 = PStarCatProp->NObj;
  TX2 = (float *)calloc(NL2,sizeof(float));
  TY2 = (float *)calloc(NL2,sizeof(float));
  X2 = (float *)calloc(NL2,sizeof(float));
  Y2 = (float *)calloc(NL2,sizeof(float));
  M2 = (float *)calloc(NL2,sizeof(float));
  POK2 = (int *)calloc(NL2,sizeof(int));
  for (I=0;I<NL2;I++) {
    TX2[I] = ShiftX+PStarCatProp->PX[I]*CosF+
      PStarCatProp->PY[I]*SinF;
    TY2[I] = ShiftY-PStarCatProp->PX[I]*SinF+
      PStarCatProp->PY[I]*CosF;
    X2[I] = PStarCatProp->PX[I];
    Y2[I] = PStarCatProp->PY[I];
    M2[I] = pow(10.0,-0.4*PStarCatProp->PM[I]);
  }

  RemoveClosePairs(X1,Y1,NL1,1,PSAPP);
  RemoveClosePairs(TX2,TY2,NL2,-1,PSAPP);

  NNeigh = (int *)calloc(NL1,sizeof(int));
  NeighList = (int **)calloc(NL1,sizeof(int *));
  NeighborList(X1,Y1,NL1,TX2,TY2,NL2,PSAPP->MaxDist,
	       NNeigh,NeighList);

  BestdX = ShiftX;
  BestdY = ShiftY;
  BestOrient = RotateAng;

  FindSolution(&BestdX,&BestdY,&BestOrient,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
	       NNeigh,NeighList,PSAPP);

  Score = CalcSigma(BestdX,BestdY,BestOrient,X1,Y1,M1,NL1,X2,Y2,M2,NL2,
		    NNeigh,NeighList,PSAPP,1,POK1,POK2,&NumC);
  free(X1);
  free(Y1);
  free(M1);
  free(POK1);
  free(TX2);
  free(TY2);
  free(X2);
  free(Y2);
  free(M2);
  free(POK2);
  *PBestdX = BestdX;
  *PBestdY = BestdY;
  *PBestRotAngle = BestOrient;
  fprintf(stdout,"Converged on %g good stars.\n",NumC);
  /*  fprintf(stdout,"Solution: dx = %.2f, dy = %.2f, theta = %.2f\n",
      BestdX,BestdY,BestOrient);*/
}

extern void ReadMatchinFiles(char *FN, StarCatRec *PStarCatProp,
			     float HDR_DX, float HDR_DY, float HDR_ROT)
{
  float **PX;
  int NL,I;

  strcpy(PStarCatProp->FN,FN);
  PStarCatProp->HDR_DX = HDR_DX;
  PStarCatProp->HDR_DY = HDR_DY;
  PStarCatProp->HDR_ROT = HDR_ROT;
  ReadASCIITable(&NL,4,&PX,FN);
  PStarCatProp->NObj = NL;
  PStarCatProp->PX = (float *)calloc(NL,sizeof(float));
  PStarCatProp->PY = (float *)calloc(NL,sizeof(float));
  PStarCatProp->PM = (float *)calloc(NL,sizeof(float));
  for (I=0;I<NL;I++) {
    PStarCatProp->PX[I] = PX[I][1];
    PStarCatProp->PY[I] = PX[I][2];
    PStarCatProp->PM[I] = PX[I][3];
  }
  freeFloatArray(PX,4,NL);
}

extern char *strsub(char *STR1, char *Str, int Index, int Length)
{
  int I;
  for (I=0; I<Length; I++) 
    STR1[I] = Str[Index-1+I];
  STR1[Length] = '\0';
  return STR1;
}

