#include <stdio.h>
#include <stdlib.h>
#include "arrays.h"
#include <math.h>

typedef struct {
  int **NBin,***Index;
  int DimX,DimY;
} BinType;

extern int BinDim(int MaxN)
{
  return (int)(sqrt((float)MaxN/4)+1);
}

extern void SetUpBins(BinType *PB, float DimX, float DimY,
		      float *X, float *Y, float *Tol, int BDimX,
		      int BDimY, int N)
{
  int J,K,L;
  float dX,dY;
  
  PB->DimX = BDimX;
  PB->DimY = BDimY;
  dX = DimX / PB->DimX;
  dY = DimY / PB->DimY;
  
  allocIntArray(&PB->NBin,PB->DimY,PB->DimX);
  for (J=0;J<2;J++) {
    for (K=0;K<PB->DimX;K++)
      for (L=0;L<PB->DimY;L++)
	(PB->NBin)[K][L] = 0;
    for (K=0;K<N;K++) {
      int BIX,EIX,BIY,EIY,IX,IY;
      
      if (K==149) {
	int II;
	II = 1;
      }
      BIX = (int)((X[K] - Tol[K]) / dX);
      EIX = (int)((X[K] + Tol[K]) / dX) + 1;
      BIY = (int)((Y[K] - Tol[K]) / dY);
      EIY = (int)((Y[K] + Tol[K]) / dY) + 1;
      for (IX=BIX;IX<EIX;IX++)
	for (IY=BIY;IY<EIY;IY++) 
	  if ((IX >= 0) && (IX < PB->DimX) &&
	      (IY >= 0) && (IY < PB->DimY)) {
	    if (J==1)
	      (PB->Index)[IX][IY][(PB->NBin)[IX][IY]] = K;
	    (PB->NBin)[IX][IY] += 1;
	  }
    }
    if (J==0) {
      PB->Index = (int ***)calloc(PB->DimX,sizeof(int **));
      for (K=0;K<PB->DimX;K++) {
	(PB->Index)[K] = (int **)calloc(PB->DimY,sizeof(int *));
	for (L=0;L<PB->DimY;L++)
	  (PB->Index)[K][L] = (int *)
	    calloc((PB->NBin)[K][L],sizeof(int));
      }
    }
  }
}

extern void CloseBins(BinType *PB)
{
  int J,K;

  for (J=0;J<PB->DimX;J++) {
    for (K=0;K<PB->DimY;K++) 
      free((PB->Index)[J][K]);
    free((PB->Index)[J]);
  }
  free(PB->Index);
  freeIntArray(PB->NBin,PB->DimY,PB->DimX);
}

extern void XYRange(float *X1, float *Y1, int N1,
		    float *MinX, float *MinY, float *MaxX, float *MaxY)
{
  int I;

  for (I=0;I<N1;I++) {
    if ((I==0) || X1[I] < *MinX)
      *MinX = X1[I];
    if ((I==0) || X1[I] > *MaxX)
      *MaxX = X1[I];
    if ((I==0) || Y1[I] < *MinY)
      *MinY = Y1[I];
    if ((I==0) || Y1[I] > *MaxY)
      *MaxY = Y1[I];
  }
}

extern void SensibleShift(float *X1, float *Y1, int N1,
			  float *MinX, float *MinY, float *DimX, float *DimY)
{
  int I;
  float MaxX,MaxY;

  for (I=0;I<N1;I++) {
    if ((I==0) || X1[I] < *MinX)
      *MinX = X1[I];
    if ((I==0) || X1[I] > MaxX)
      MaxX = X1[I];
    if ((I==0) || Y1[I] < *MinY)
      *MinY = Y1[I];
    if ((I==0) || Y1[I] > MaxY)
      MaxY = Y1[I];
  }
  for (I=0;I<N1;I++) {
    X1[I] -= *MinX;
    Y1[I] -= *MinY;
  }
  if (N1 > 0) {
    *DimX = (MaxX - *MinX)*1.01;
    *DimY = (MaxY - *MinY)*1.01;
    if (*DimX < 0.001) 
      *DimX = 0.001;
    if (*DimY < 0.001) 
      *DimY = 0.001;
  }    
}

void SensibleShift2(float *X1, float *Y1, int N1,
		    float *X2, float *Y2, int N2,
		    float *MinX, float *MinY,
		    float *DimX, float *DimY)
{
  int I;
  float MaxX,MaxY;

  MaxX = 0;
  MaxY = 0;
  for (I=0;I<N1;I++) {
    if ((I==0) || X1[I] < *MinX)
      *MinX = X1[I];
    if ((I==0) || X1[I] > MaxX)
      MaxX = X1[I];
    if ((I==0) || Y1[I] < *MinY)
      *MinY = Y1[I];
    if ((I==0) || Y1[I] > MaxY)
      MaxY = Y1[I];
  }
  for (I=0;I<N2;I++) {
    if ((I+N1==0) || X2[I] < *MinX)
      *MinX = X2[I];
    if ((I+N1==0) || X2[I] > MaxX)
      MaxX = X2[I];
    if ((I+N1==0) || Y2[I] < *MinY)
      *MinY = Y2[I];
    if ((I+N1==0) || Y2[I] > MaxY)
      MaxY = Y2[I];
  }
  for (I=0;I<N1;I++) {
    X1[I] -= *MinX;
    Y1[I] -= *MinY;
  }
  for (I=0;I<N2;I++) {
    X2[I] -= *MinX;
    Y2[I] -= *MinY;
  }
  if (N1+N2 > 1) {
    *DimX = (MaxX - *MinX)*1.01;
    *DimY = (MaxY - *MinY)*1.01;
    if (*DimX < 0.001) 
      *DimX = 0.001;
    if (*DimY < 0.001) 
      *DimY = 0.001;
  }
  else if (N1+N2 >= 0) {
    if (MaxX > 0) *DimX = MaxX;
    else *DimX = 1.0;
    if (MaxY > 0) *DimY = MaxY;
    else *DimY = 1.0;
  }
}

extern void NeighborList(float *X1, float *Y1, int N1,
			 float *X2, float *Y2, int N2,
			 float Tol, int *NNeigh,
			 int **NeighList)
{
  int TN,N[2];
  float *X[2],*Y[2],*TolA;
  int I,J,K,L;
  BinType B[2];
  float MinX,MinY,DimX,DimY;
  int *Count;
  float *CX1,*CY1,*CX2,*CY2;

  CX1 = (float *)calloc(N1,sizeof(float));
  CY1 = (float *)calloc(N1,sizeof(float));
  for (I=0;I<N1;I++) {
    CX1[I] = X1[I];
    CY1[I] = Y1[I];
  }
  CX2 = (float *)calloc(N2,sizeof(float));
  CY2 = (float *)calloc(N2,sizeof(float));
  for (I=0;I<N2;I++) {
    CX2[I] = X2[I];
    CY2[I] = Y2[I];
  }

  if (N1 > N2) TN = N1;
  else TN = N2;
    
  SensibleShift2(CX1,CY1,N1,CX2,CY2,N2,&MinX,&MinY,&DimX,&DimY);

  if (DimX < Tol) DimX = Tol;
  if (DimY < Tol) DimY = Tol;
  N[0] = N1;
  N[1] = N2;
  X[0] = CX1;
  X[1] = CX2;
  Y[0] = CY1;
  Y[1] = CY2;

  TolA = (float *)calloc(TN,sizeof(float));
  SetUpBins(&B[0],DimX,DimY,X[0],Y[0],TolA,BinDim(TN),BinDim(TN),N[0]);
  for (I=0;I<TN;I++)
    TolA[I] = Tol;
  SetUpBins(&B[1],DimX,DimY,X[1],Y[1],TolA,BinDim(TN),BinDim(TN),N[1]);
  free(TolA);

  for (L=0;L<2;L++) {
    for (I=0;I<B[0].DimX;I++)
      for (J=0;J<B[0].DimY;J++)
	for (K=0;K<B[0].NBin[I][J];K++) {
	  int VK;

	  for (VK=0;VK<B[1].NBin[I][J];VK++) {
	    float DiffX,DiffY,Dist;
	    
	    DiffX = X[0][B[0].Index[I][J][K]] -
	      X[1][B[1].Index[I][J][VK]];
	    DiffY = Y[0][B[0].Index[I][J][K]] -
	      Y[1][B[1].Index[I][J][VK]];
	    Dist = sqrt(DiffX*DiffX+DiffY*DiffY);
	    if (Dist < Tol) {
	      if (L==0) NNeigh[B[0].Index[I][J][K]]++;
	      else {
		NeighList[B[0].Index[I][J][K]][Count[B[0].Index[I][J][K]]]
		  = B[1].Index[I][J][VK];
		Count[B[0].Index[I][J][K]]++;
	      }
	    }
	  }
	}
    if (L==0) {
      Count = (int *)calloc(N1,sizeof(int));
      for (I=0;I<N1;I++)
	NeighList[I] = (int *)calloc(NNeigh[I],sizeof(int));
    }
    else free(Count);
  }
  for (I=0;I<2;I++)
    for (J=0;J<N[I];J++) {
      X[I][J] += MinX;
      Y[I][J] += MinY;
    }
    
  for (I=0;I<2;I++)
    CloseBins(&B[I]);
  free(CX1);
  free(CY1);
  free(CX2);
  free(CY2);
}

extern void MatchUp(float *X1, float *Y1, int N1,
		    float *X2, float *Y2, int N2,
		    float Tol, int *PIdent)
{
  float DiffX,DiffY,Dist,BestD;
  int *NNeigh,**NeighList;
  int I,J;

  NNeigh = (int *)calloc(N1,sizeof(int));
  NeighList = (int **)calloc(N1,sizeof(int *));
  NeighborList(X1,Y1,N1,X2,Y2,N2,Tol,NNeigh,NeighList); 
  for (I=0;I<N1;I++) {
    PIdent[I] = -1;
    for (J=0;J<NNeigh[I];J++) {
      DiffX = X1[I] - X2[NeighList[I][J]];
      DiffY = Y1[I] - Y2[NeighList[I][J]];
      Dist = sqrt(DiffX*DiffX+DiffY*DiffY);
      if ((J==0)||(Dist<BestD)) {
	BestD = Dist;
	PIdent[I] = NeighList[I][J];
      }
    }
  }
  free(NNeigh);
  for (I=0;I<N1;I++)
    free(NeighList[I]);
  free(NeighList);
}

extern void InvertMatchUp(int N1, int N2, int *PIdent1, int *PIdent2)
{
  int I;

  for (I=0;I<N1;I++)
    PIdent1[I] = -1;
  for (I=0;I<N2;I++)
    if (PIdent2[I]>-1)
      PIdent1[PIdent2[I]] = I;
}

extern FILE *fopenv(char *PathName, char *Access)
{
  FILE *TEMP;

  TEMP = fopen(PathName,Access);
  if (TEMP == NULL) {
    if (Access[0] == 'r')
      fprintf(stderr,"**ERROR**: Can't read file %s\n",PathName);
    else if (Access[0] == 'w')
      fprintf(stderr,"**ERROR**: Can't write file %s\n",PathName);
    perror(PathName);
    exit(1);
  }
  return TEMP;
}

extern void GetLines(int *PNL, char *FStr)
{
  char Str[600];
  FILE *fin;

  *PNL = 0;
  fin = fopenv(FStr,"r");
  while (EOF != fscanf(fin,"%[^\n]",Str)) {
    fgetc(fin);
    (*PNL)++;
  }
  fclose(fin);
}

extern void ReadASCIITable(int *PNL, int NumEntry, float ***PX, char *FStr)
{
  FILE *fin;
  char Str[600],TStr[600];
  float **X;

  *PNL = 0;
  fin = fopenv(FStr,"r");
  while (EOF != fscanf(fin,"%[^\n]",Str)) {
    int i,bi,ei,j,tooearly;
    fgetc(fin);
    
    i = 0;
    tooearly = 0;
    for (j=0;j<NumEntry;j++) {      
      while ((Str[i] == ' ')||(Str[i] == ':')||(Str[i]==9)) i++;
      bi = i;
    
      while ((Str[i] != ' ')&&(Str[i] != ':')&&(Str[i]!=9)&&(Str[i] != 0)) i++;
      ei = i;
      if ((Str[i] ==0)&&(j<NumEntry-1)) {
	j = NumEntry;
	tooearly = 1;
      }
    }
    if (tooearly == 0) (*PNL)++;
  }
  fclose(fin);
  allocFloatArray(PX,NumEntry,*PNL);
  X = *PX;

  *PNL = 0;
  fin = fopenv(FStr,"r");
  while (EOF != fscanf(fin,"%[^\n]",Str)) {
    int i,bi,ei,j,tooearly;
    fgetc(fin);
    
    i = 0;
    tooearly = 0;
    for (j=0;j<NumEntry;j++) {      
      while ((Str[i] == ' ')||(Str[i] == ':')||(Str[i]==9)) i++;
      bi = i;
    
      while ((Str[i] != ' ')&&(Str[i] != ':')&&(Str[i]!=9)&&(Str[i] != 0)) i++;
      ei = i;
      if ((Str[i] ==0)&&(j<NumEntry-1)) {
	j = NumEntry;
	tooearly = 1;
      }
    }
    if (tooearly == 0) {
      i = 0;
      for (j=0;j<NumEntry;j++) {      
	while ((Str[i] == ' ')||(Str[i] == ':')||(Str[i]==9)) i++;
	bi = i;
	
	while ((Str[i] != ' ')&&(Str[i] != ':')&&(Str[i]!=9)&&(Str[i] != 0)) 
	  i++;
	ei = i;

	strsub(TStr,Str,bi+1,ei-bi);
	sscanf(TStr,"%g",&X[*PNL][j]);
      }
      
      (*PNL)++;
    }
  }
  fclose(fin);
}

