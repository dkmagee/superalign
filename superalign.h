int Verbose;

/* CoordListRec are the structures which contain the object positions after transformation to the universal frame */

typedef struct {
  int NObj;
  float *PX;  /* x position of object (transformed frame) */
  float *PY;  /* y position of object (transformed frame) */
  float *PM;  /* mag of object -- not currently used */
} CoordListRec;

/* StarCatRec is the structure which contains the (undistorted) object positions in the original frame -- except for HDR_DX,HDR_DY,HDR_ROT, it contents should not change during operation */

typedef struct {
  char FN[150];  /* Name of file where all these object positions, etc. are given */
  int NObj;  /* Number of compact objects */
  float *PX; /* x position of object */
  float *PY; /* y position of object */
  float *PM; /* mag of object -- not used */
  float HDR_DX; /* current x-shift that needs to be applied to objects in this frame to align it with the universal solution */
  float HDR_DY; /* current y-shift that needs to be applied to objects in this frame to align it with the universal solution */
  float HDR_ROT;  /* current rotation angle that needs to be applied to objects in this frame to align it with the universal solution */
} StarCatRec;

/* PointingRec contains all the exposures which pertain to a specific grouping */

typedef struct {
  int NExp;
  StarCatRec *PStarCatProp;
} PointingRec;

/* ExposureRec is the structure which contains the (undistorted) object positions in the original frame */

typedef struct {
  char FN[150];
  StarCatRec StarCatProp;
  float hdr_dx; /* x-shift that needs to be applied to the contents of StarCatProp */
  float hdr_dy; /* y-shift that needs to be applied to the contents of StarCatProp */
  float hdr_rot;  /* rotation angle that needs to be applied to the contents of StarCatProp */
} ExposureRec;

/* AllExposureRec contains all the exposures that need to be aligned with respect to each other */

typedef struct {
  int NExp;   /* Number of Exposures */
  ExposureRec *PExposureProp;  
} AllExposureRec;

typedef struct {
  float TypDist;   /* Typical RMS error for objects which are well aligned */
  float MaxDist;   /* MaxDist is equal to the distance two objects can be 
                      apart to even register as matches.  Even if objects
                      register as matches -- this does not imply they will
                      have a huge effect on the final score. */
  float TolPair;   /* Objects with Neighbors Closer than TolPair pixels 
                      are thrown out */
  float MaxShiftErr;  /* Maximum Pixel Error in X/Y Shifts */
  float MaxRotErr; /* Maximum Error in Rotation Angle (Degrees) */
  float GroupSize; /* Maximum Distance Pointings can have to be grouped 
                      together */
} SuperAlignParmRec;

extern void ReadASCIITable(int *PNL, int NumEntry, float ***PX, char *FStr);

extern int TryAlignment(StarCatRec *PStarCatBaseProp,
			StarCatRec *PStarCatProp,
			SuperAlignParmRec *PSAPP,
			float *PBestdX, float *PBestdY,
			float *PBestRotAngle);

extern int ImproveAlignment(StarCatRec *PStarCatBaseProp,
			    StarCatRec *PStarCatProp,
			    SuperAlignParmRec *PSAPP,
			    float *PBestdX, float *PBestdY,
			    float *PBestRotAngle);

extern void NeighborList(float *X1, float *Y1, int N1,
			 float *X2, float *Y2, int N2,
			 float Tol, int *NNeigh,
			 int **NeighList);

extern void DetermineGuess(StarCatRec *PStarCatBaseProp,
			   StarCatRec *PStarCatProp,
			   float *PShiftX, float *PShiftY,
			   float *PRotateAng);


extern void ReadMatchinFiles(char *FN, StarCatRec *PStarCatProp,
			     float HDR_DX, float HDR_DY, float HDR_ROT);

extern float ran3(long *idum);

extern void chdecomp(float **a, int n, float p[], int *err);

extern void chsvx(float **a, int n, float p[], float b[], float x[]);
