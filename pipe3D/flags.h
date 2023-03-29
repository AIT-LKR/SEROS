#ifndef FLAGS_H
#define FLAGS_H

const int flagConstraintSolid =  0;   // cells which may not ever erode (border and buffer)
const int flagSolid           =  5;   // BB cells not in contact with fluid cells
const int flagBB              =  85;  // BB cells directly in contact with fluid cells
const int flagSedimenting     =  95;  // wet or solid cells which change into BB
const int flagWet             =  128; // fluid cells directly in contact with BB cells
const int flagEroding         =  140; // BB cells which change into fluid
const int flagBulk            =  170; // fluid cells not in contact with BB cells
const int flagConstraintBB    =  220;  // 'frozen' areas within seros region
const int flagBuffer          =  255; // fluid cells which may not ever sediment
const int flagInlet           =  flagBuffer;
const int flagOutlet          =  flagBuffer;


#endif
