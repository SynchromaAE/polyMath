#define GROUPS 256
#define SLOTS 128
#define MAXSEQ 2048
#define VARIATIONS 6

/* polyMath~ - sequential polyrhythm generator
 * Edward Kelly 2017 - 2021
 * 
 * TO DO:1) seqInSlot for writing sequences 1 element at a time, and associated initSeqSlot, with addition of seq.filled copying in scramble
 * TO DO:2) moveInSlot and associated sequence re-writing function
 * TO DO:3) it would be nice to be able to toggle between tuple and linear phase sequence by performing some sort of quantize of a linear seq
 * TO DO:4) Scramble groups, and also scramble an individual group or internal sequences of groups
 * TO DO:5) Insert group
 *       6) Nested subgroup e.g. make a 3 inside a 5:
 * 20 30 30 30 20 20   |      |    |    |    |      |      |
 * 16ths               |        |        |        |        |
 * ...easy to work out (7.5th or 2x15th) but the ability to command "split two 5ths into 3" is...just bump the sequence down one, change one entry and create a new entry after it (into the "duplication")...or just write the whole sequence (and groups) from the change to the end!
 */

/* data types
 * full sequences:
 * dType --- dataOut
 * 0     --- eOff: {phaseOff[0], phaseOff[1], phaseOff[2] ... phaseOff[len-1]}
 * 1     --- eSize: {sizePhase[0], sizePhase[1], sizePhase[2] ... sizePhase[len-1]}
 * 2     --- groupStep: {gStep[0], gStep[1], gStep[2] ... gStep[len-1]}
 * 3     --- groupNum:  {gNum[0], gNum[1], gNum[2] ... gNum[len-1]}
 * 4     --- eJoin (is the event joined to the next event) - not fully implemented yet
 * 5     --- jSize (length of the join - this could be eliminated if eJoin was jSize)
 *
 * 11    --- pAcc1
 * 12    --- eAcc1
 * 13    --- pAcc2
 * 14    --- eAcc2
 * 15    --- pAcc3
 * 16    --- eAcc3
 * 17    --- pAcc4
 * 18    --- eAcc4
 * 19    --- pAcc5
 * 20    --- eAcc5
 * 21    --- pAcc6
 * 22    --- eAcc6
 * 23    --- pAcc7
 * 24    --- eAcc7
 * 25    --- pAcc8
 * 26    --- eAcc8
 * group values:
 * 91    --- {gType, nGroups, cycles, len, remains, 0}
 * 92    --- group numerators
 * 93    --- group denominators
 * 94    --- group offsets (phase)
 * 95    --- group size in phase
 * 96    --- group start (events)
 * 97    --- 
 *       --- 
 * 98    --- {slot, length}
 * 99    --- {slot, variation}
 * since any phase value has a reciprocal, gType 0 could have only whole denom (e.g 4/16) 
 * and gType 1 fractional (e.g. 1/5.33333). gType 1 would mostly be 1/something e.g.
 * {1/3.25, 1/2.18, 1/4.28096} adds up to 1
 * sequence lengths that are not equal to an integer number of phase cycles will
 * be dealt with in the grand (external to this object) structure using the isobars object. 
 */


#ifdef __APPLE__
#include <sys/types.h>
#include <sys/time.h>
#endif

#include <stdlib.h>
#include <time.h>

#include "m_pd.h" 

#define EVENTLIST 27

static t_class *polyMath_tilde_class;

typedef struct drand48_data *randomize;

typedef struct {
  long tv_sec;
  long tv_usec;
} timeval;

typedef struct _groups
{
  int gType[SLOTS];                  // NEW (and as yet undefined at 1 Dec 2018) tuples (gType = 0) or
                                     // linear phase (gType = 1)
  int nGroups[SLOTS];                // number of groups in this sequence
  int gStart[GROUPS * SLOTS];        // where in the sequence does each group start?
  t_atom n[GROUPS * SLOTS];          // numerator of the time sig (fraction)
  t_atom d[GROUPS * SLOTS];          // denominator of the time sig
  int cycles[SLOTS];                 // number of cycles of phasor~
  t_atom offset[GROUPS * SLOTS];     // group offset in phase
  t_atom size[GROUPS * SLOTS];       // group size in phase
  t_atom sizeInv[GROUPS * SLOTS];    // 1 / size - so that look-up table can be used instead of math at runtime
  //t_atom rLength[GROUPS * SLOTS];  // real length of sequence calculated to float precision - precision is now set globally
  t_atom remains[GROUPS * SLOTS];    // how much of the phase is left (in case of a 0 entry - skip one and if still 0 use remainder)
                                     // also, if an incomplete list is issued, the last value will be based on this!
                                     // THIS IS NOT WELL DEFINED IN THE CODE AS OF 17th December 2017
  int fillGroup[SLOTS];              // if a fill group is used for an incomplete setGroups call, the index will be stored here
  int isUnFilled[SLOTS];               // when the slots are uninitialized it is a 1
} t_groups;                         

typedef struct _vars
{
  int gType[SLOTS * VARIATIONS];
  int nGroups[SLOTS * VARIATIONS];                // number of groups in this sequence
  int gStart[GROUPS * SLOTS * VARIATIONS];        // where in the sequence does each group start?
  t_atom n[GROUPS * SLOTS * VARIATIONS];          // numerator of the time sig (fraction)
  t_atom d[GROUPS * SLOTS * VARIATIONS];          // denominator of the time sig
  int cycles[SLOTS * VARIATIONS];                 // number of cycles of phasor~
  t_atom offset[GROUPS * SLOTS * VARIATIONS];     // group offset in phase
  t_atom size[GROUPS * SLOTS * VARIATIONS];       // group size in phase
  t_atom sizeInv[GROUPS * SLOTS * VARIATIONS];    // 1 / size - so that look-up table can be used instead of math at runtime
  //t_atom rLength[GROUPS * SLOTS * VARIATIONS];    // real length of sequence calculated to float precision
  t_atom remains[GROUPS * SLOTS * VARIATIONS];    // how much of the phase * cycles is left (in case of a 0 entry - skip one and if still 0 use remainder)
                                     // also, if an incomplete list is issued, the last value will be based on this!
  int fillGroup[SLOTS];              // if a fill group is used for an incomplete setGroups call, the index will be stored here
  int swaps[MAXSEQ];
  int swapsRef[MAXSEQ * 2];
  int swapped[MAXSEQ];
  int groupSwaps[GROUPS];

} t_vars;                

typedef struct _variations
{
  int len[SLOTS * VARIATIONS];                   // how many events are in the current sequence
  int variations[SLOTS * VARIATIONS];
  int nGroups[SLOTS * VARIATIONS];
  int excludes[SLOTS * VARIATIONS * MAXSEQ]; // where a join has been implemented, these should not be available to scramble
  t_atom varStep[MAXSEQ * SLOTS * VARIATIONS];
  t_atom allStep[MAXSEQ * SLOTS * VARIATIONS];
  //NEW FOR seqInSlot, Jan 2019
  t_atom filled[MAXSEQ * SLOTS * VARIATIONS];   // for seqInSlot - this will determine which elements are filled, autofilled, and unfilled
  t_atom groupStep[MAXSEQ * SLOTS * VARIATIONS]; // which element of the current group are we in?
  t_atom groupNum[MAXSEQ * SLOTS * VARIATIONS];  // which sequential group are we in?
  //  t_atom oGStep[MAXSEQ * SLOTS * VARIATIONS]; // which element of the original group are we in?
  //  t_atom oGNum[MAXSEQ * SLOTS * VARIATIONS];  // which original sequential group are we in?
  t_atom eSize[MAXSEQ * SLOTS * VARIATIONS];     // individual event size in phase
  t_atom eOff[MAXSEQ * SLOTS * VARIATIONS];      // individual event offset in phase
  t_atom eJoin[MAXSEQ * SLOTS * VARIATIONS];     // events joined together (e.g. 3/16ths as one note) - affects eSizeInv (below). Normally 1.
  t_atom jSize[MAXSEQ * SLOTS * VARIATIONS];     // joined event(s) size in phase
  t_atom eAcc1[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom eAcc2[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom eAcc3[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom eAcc4[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc1[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc2[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc3[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc4[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom eSizeInv[MAXSEQ * SLOTS * VARIATIONS];  // 1 / eSize for chunk phase output. If eJoin > 1 then there will be e.g. 1 / (3/16) followed by two '0'
  t_atom denom[MAXSEQ * SLOTS * VARIATIONS];
  t_atom varOff[MAXSEQ * SLOTS * VARIATIONS];
  t_atom grpOff[MAXSEQ * SLOTS * VARIATIONS];

  t_atom debugList[17];

  //new for 2019
  t_atom altOff[MAXSEQ * SLOTS * VARIATIONS];      // alternative event offset in phase
  
  t_atom eAcc5[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom eAcc6[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom eAcc7[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom eAcc8[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc5[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc6[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc7[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event
  t_atom pAcc8[MAXSEQ * SLOTS * VARIATIONS];     // accent or parameter storage for event

} t_variations;

typedef struct _sequences
{//allStep filled groupStep groupNum eSize eOff eJoin jSize eAcc1-8 pAcc1-8 eSizeInv denom altOff
  int len[SLOTS];                   // how many events are in the current sequence
  t_atom allStep[MAXSEQ * SLOTS];   // which event of the total sequence is this?
  //NEW FOR seqInSlot, Jan 2019
  t_atom filled[MAXSEQ * SLOTS];   // for seqInSlot - this will determine which elements are filled, autofilled, and unfilled
  t_atom groupStep[MAXSEQ * SLOTS]; // which element of the current group are we in?
  t_atom groupNum[MAXSEQ * SLOTS];  // which sequential group are we in?
  t_atom eSize[MAXSEQ * SLOTS];     // individual event size in phase
  t_atom eOff[MAXSEQ * SLOTS];      // individual event offset in phase
  t_atom eJoin[MAXSEQ * SLOTS];     // events joined together (e.g. 3/16ths as one note) - affects eSizeInv (below). Normally 1.
  t_atom jSize[MAXSEQ * SLOTS];     // joined event(s) size in phase
  t_atom eAcc1[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom eAcc2[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom eAcc3[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom eAcc4[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc1[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc2[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc3[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc4[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom eSizeInv[MAXSEQ * SLOTS];  // 1 / eSize for chunk phase output. If eJoin > 1 then there will be e.g. 1 / (3/16) followed by two '0'
                                    // entries. '0' entries will then cause the algorithm to default to the 1 / (3/16) or more like 3 * (1/16)
  // but do they? more work needed...
  t_atom pList1[2];
  t_atom pList2[2];
  t_atom pList3[2];
  t_atom pList4[2];
  t_atom wrapCycles1[MAXSEQ];
  t_atom wrapCycles2[MAXSEQ];

  t_atom denom[MAXSEQ * SLOTS];

  //new for 2019
  t_atom altOff[MAXSEQ * SLOTS];      // alternative event offset in phase

  t_atom eAcc5[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom eAcc6[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom eAcc7[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom eAcc8[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc5[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc6[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc7[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pAcc8[MAXSEQ * SLOTS];     // accent or parameter storage for event
  t_atom pList5[2];
  t_atom pList6[2];
  t_atom pList7[2];
  t_atom pList8[2];

} t_sequences;                      

typedef struct _polyMath_tilde
{
  t_object x_obj;
  t_float f;
  t_groups grp;
  t_sequences seq;
  t_variations var;
  t_vars vGrp;

  int SEQSIZE, VARSIZE, GROUPSIZE, VGROUPSIZE;
  
  int firstStart;

  t_atom outList[MAXSEQ];
  t_atom eventList[EVENTLIST];
  t_atom dList[2]; // next duration / phase
  
  unsigned short int seed16v[3];
  unsigned long int timeSeed;

  //rounder
  int iRound;
  t_float fRound, rDiff;
  //writeGroup
  t_float Wsize, WsizeRem, WESize, groupOffset, Woff, WSInv, WOffAcc;
  int d, Gstart, IsizeRem;
  //int groupWrite;
  //setGroups
  t_float Goff, Gcycle, Eoff, Grem;
  int c, e, a, b, h, g, i;
  int Pslot, Pac, PLStep, PStep;
  t_float PGcyc;
  //getVariables
  t_float clockOut, E_Acc1, E_Acc2, E_Acc3, E_Acc4, Pacc1, Pacc2, Pacc3, Pacc4, E_Acc5, E_Acc6, E_Acc7, E_Acc8, Pacc5, Pacc6, Pacc7, Pacc8, Pthis, PJoin;
  int Gnm;
  t_float Gstep, ESize, ESInv, Gn, Gd, GSize, GSInv;
  int cycles;
  t_float InVal, PreVal, TotVal, PStepOff;
  //leftovers, perform and joins:
  int barNew, join, Joined, myBug, pJoin, PJoined, JoinVal, JoinTot, j, k, RW, l, m;
  int slot, Wstep, Icycle, Gstp, GroupStart, PSlot, JGstt, JGnm;
  int jFlag, jFirst, joinSuccess, sortFlag; // FLAGS
  int JSlot, JGrp, JLoc, JLen, JBuf, JGst, initSlot;
  // care must be taken to reset jFlag and jFirst when entering new slot, or writing groups! 30/10/2017
  t_float PESize, PEOff, PESInv, JPESI, group, JGSize, JESize, JJoin, JGn, JGd, JGt;
  int maxGrp, fGroup;
  //aternate signal outs for event seg~
  int altOut, altNum, eChanged, preChange, eMult, altEarly;
  t_float percentVal, eOut, eVal;

  //groupThisSlot / jumpTo / jumpNext
  int thisSlot, changeSlot, changeVar, NStep;
  t_float wrapSubVal, nextShotVal;
  
  int JlastOffset, JnextOffset;
  int JlastLen, JnextLen, JiWrap, JnextFlag, JlocateFlag;
  t_float JlastCycle, JnextCycle;
  t_float JsizeNext, JoffNext, JwrapCycle;

  //scramble - REFACTOR FROM HERE - which are eventually used?
  int o, p, q, r, s, t, u;
  int scramMeth, seqLen, halfSeq, swapsNum, doSwaps, swapNdx1, swapNdx2, swapFlag, iFSwapsNum, offsetVar, grpOffset, noRepeats;

  int variation, thisVar, scramSlot, varTest, varPerf, scrambling;
  t_float copyVal, swapVal, seqProb, getD;

  t_float fSeqLen, fHalfSeq, fSwapsNum;
  int copyWell, scramWell, swapWell, groupWell, varWrite;
  double randNum1, randNum2, randNum3;
  //int extraSwapFlag; // deprecated: using rounder() function instead

  int VGnm, VGCount, VEJoin;
  t_float Vd, VESize, VGSize, VGSizeInv, VEOff, VGOff, VJSize, VJoin, VVStep, VVLast, VLastD, VONext, instant, gInstant, eVVal, varOff, VOffG, VPESI, VOff;
  int swapVal1, swapVal2, newVar, setInstant;

  int autoThreshold;
  t_float cycleDiff, sizeThreshold, halfSize, sizeFrac;
  
  //int wrapNextState, wrapThisLen, wrapNextCycle, nextSlotFlag, unScrambleNext, nextScramble, scrambleNext, SStep
  //t_float cycleWrap, cycleMod, nextSlotWrap, wrapOff, wrapOffSize, floatNextCycle, nextSlotVal, nextSlotWrap, slotSwitchVal
  int lastLen, nextLen, nextSlot, nextSlotVal, lastSlot, nextVar, validJumpState;
  int zeroNextPhase, zeroNextVar, zeroNextSlot, lastVar, wrapLen, swapState; // 8th Jan 2018 - preparation for change-slot/variation-on-next-event
  t_float offNext, thisInVal; // 8th Jan 2018 - preparation for change-slot/variation-on-next-event
  int jumpSlotAtEnd, jumpVarAtEnd;

  //polyMath_tilde_getSeq
  int getSlot, getVar, getPar, v, getVarNum;
  int grpOff, seqOff, seqGrpOff, lenSeq, lenGrp;
  t_float getSeqVal, getGrpVal;

  //polyMath_tilde_seqInSlot
  int seqGrpOffset;
  int seqSlotOffset;
  int w, sType, seqPos, seqNum, seqDen, prevSNum, prevSDen;
  t_float seqPhase, seqPOff, prevSPhase;

  //polyMath_tilde_seqInSlot and swap variable
  int swapSlot, swapVar, swapLoc, swapShift, x, swapLength;
  t_float swapP, swapE;
  int isSwapList;
  int isShuffled;

  int pageNum;
  t_float durBeat, barBeat, dur1, dur2, BPM, dPhase;
  int altLen, y, z;
  //int pageFlag, pageNum, lastPage;
  //float fPageNum;

  //groupScramble
  int GSMode, GSSlot, GSVar, GSDestVar, GSISwap;
  t_float GSScramRand, GSFSwap;
  
  //slotLen
  int isLength, getSlotLen;
  
  t_clock *fOut, *early, *pageTurner;
  t_outlet *clock, *subclock; // from v1
  t_outlet *cycle, *newgroup, *newbar, *p1, *p2, *p3, *p4, *p5, *p6, *p7, *p8, *groupnum, *num, *denom; // from v1
  t_outlet *eventLengthPhase, *eventLengthNum, *alt, *eChange, *eAlt, *page;
  t_outlet *dataOut, *dType, *durFirst, *durAlt; // list outlets to communicate with app
  // new 30th Oct 2017, for joined-clock events in phase (e.g. 0.125) and num (e.g. 2x16ths) and alt versions (for priming playback subpatches)
  // new 28th November 2018, sequences output from rightmost outlet...I have yet to write any code for this (15:21PM, 28th Nov 2018)
} t_polyMath_tilde;

int rounder(t_polyMath_tilde *x, t_float f, int limit) // limiting round function
{
  x->fRound = f;
  x->iRound = (int)f;
  x->rDiff = (t_float)x->iRound - x->fRound;
  if(x->rDiff >= 0.5) x->iRound++;
  if(x->iRound > limit) x->iRound = limit;
  return(x->iRound);
}

static void getVariables(t_polyMath_tilde *x)
{
  x->clockOut = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.allStep);
  x->E_Acc1 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc1);
  x->E_Acc2 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc2);
  x->E_Acc3 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc3);
  x->E_Acc4 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc4);
  x->Pacc1 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc1);
  x->Pacc2 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc2);
  x->Pacc3 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc3);
  x->Pacc4 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc4);
  x->E_Acc5 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc5);
  x->E_Acc6 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc6);
  x->E_Acc7 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc7);
  x->E_Acc8 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eAcc8);
  x->Pacc5 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc5);
  x->Pacc6 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc6);
  x->Pacc7 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc7);
  x->Pacc8 = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.pAcc8);
  SETFLOAT(&x->seq.pList1[0], x->Pacc1); SETFLOAT(&x->seq.pList1[1], x->E_Acc1);
  SETFLOAT(&x->seq.pList2[0], x->Pacc2); SETFLOAT(&x->seq.pList2[1], x->E_Acc2);
  SETFLOAT(&x->seq.pList3[0], x->Pacc3); SETFLOAT(&x->seq.pList3[1], x->E_Acc3);
  SETFLOAT(&x->seq.pList4[0], x->Pacc4); SETFLOAT(&x->seq.pList4[1], x->E_Acc4);
  SETFLOAT(&x->seq.pList5[0], x->Pacc5); SETFLOAT(&x->seq.pList5[1], x->E_Acc5);
  SETFLOAT(&x->seq.pList6[0], x->Pacc6); SETFLOAT(&x->seq.pList6[1], x->E_Acc6);
  SETFLOAT(&x->seq.pList7[0], x->Pacc7); SETFLOAT(&x->seq.pList7[1], x->E_Acc7);
  SETFLOAT(&x->seq.pList8[0], x->Pacc8); SETFLOAT(&x->seq.pList8[1], x->E_Acc8);
  if(x->myBug == 4) post("P2 = %f, E2 = %f, Location = %d",atom_getfloatarg(0,2,x->seq.pList2),atom_getfloatarg(1,2,x->seq.pList2),x->slot * MAXSEQ + x->PStep);
  //x->Pthis = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eSize);
  x->PJoin = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eJoin);
  //  if(x->Pthis == 0 && x->PJoin > 1) x->Pthis = x->PJoin;

  //assignment of PJoined happens here, and then the value is manipulated in perform. See "FLAGS"
  //  if(x->PJoin > 1) x->PJoined = x->PJoin; // see below
  x->Gnm = (int)atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.groupNum);
  x->Gstep = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.groupStep);
  // trying this in perform, since it now inhabits a signal outlet:
  //x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eOff);
  x->PESize = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eSize);
  x->PESInv = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eSizeInv);
  // 2017 30th October:
  // FLAGS::::::In the clockGen_tilde_perform function
  /* When a value of an element of x->seq.eJoin is encountered that is greater than 1, its value is passed to x->PJoined.
   * The JPESI value is set to 1 / the length of PESiz * x->PJoined.
   * A flag, firstJoin is used to indicate that this was just set - if(x->firstJoin == 0) x->firstJoin = 1;. Another flag, jFlag, is set to 1.
   * Next step, the firstJoin flag is revoked - if(x->firstJoin == 1) x->firstJoin = 0; 
   * but if x->jFlag remains 1 and x->PJoined > 0, x->JPESI is not reset.
   * if x->PJoined == 0 (i.e. the else statement after if(x->PJoined > 0) then x->JPESI is reset to x->PESiz
   *
   * This function, getVariables, is kept as minimal as possible, although cout will include a new variable (and outlet) or two
   * eventLengthNum - rightmost outlet, so that other aspects can be kept at bay, or allowed to proceed if necessary via spigots in Pd
   */
  x->Gn = atom_getfloatarg(x->slot * GROUPS + x->Gnm, x->GROUPSIZE, x->grp.n);
  x->Gd = atom_getfloatarg(x->slot * GROUPS + x->Gnm, x->GROUPSIZE, x->grp.d);
  x->GSize = atom_getfloatarg(x->slot * GROUPS + x->Gnm, x->GROUPSIZE, x->grp.size);
  x->GSInv = atom_getfloatarg(x->slot * GROUPS + x->Gnm, x->GROUPSIZE, x->grp.sizeInv);
  x->Goff = atom_getfloatarg(x->slot * GROUPS + x->Gnm, x->GROUPSIZE, x->grp.offset);
  x->Grem = atom_getfloatarg(x->slot, SLOTS, x->grp.remains);
  //if(x->GSize > 0) x->GSInv = 1 / x->GSize;
  x->cycles = x->grp.cycles[x->slot];
  x->eChanged = 0;
}

static void getVariations(t_polyMath_tilde *x)
{
  x->clockOut = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.allStep);
  x->E_Acc1 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc1);
  x->E_Acc2 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc2);
  x->E_Acc3 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc3);
  x->E_Acc4 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc4);
  x->Pacc1 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc1);
  x->Pacc2 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc2);
  x->Pacc3 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc3);
  x->Pacc4 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc4);
  x->E_Acc5 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc5);
  x->E_Acc6 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc6);
  x->E_Acc7 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc7);
  x->E_Acc8 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc8);
  x->Pacc5 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc5);
  x->Pacc6 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc6);
  x->Pacc7 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc7);
  x->Pacc8 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc8);
  SETFLOAT(&x->seq.pList1[0], x->Pacc1); SETFLOAT(&x->seq.pList1[1], x->E_Acc1);
  SETFLOAT(&x->seq.pList2[0], x->Pacc2); SETFLOAT(&x->seq.pList2[1], x->E_Acc2);
  SETFLOAT(&x->seq.pList3[0], x->Pacc3); SETFLOAT(&x->seq.pList3[1], x->E_Acc3);
  SETFLOAT(&x->seq.pList4[0], x->Pacc4); SETFLOAT(&x->seq.pList4[1], x->E_Acc4);
  SETFLOAT(&x->seq.pList5[0], x->Pacc5); SETFLOAT(&x->seq.pList5[1], x->E_Acc5);
  SETFLOAT(&x->seq.pList6[0], x->Pacc6); SETFLOAT(&x->seq.pList6[1], x->E_Acc6);
  SETFLOAT(&x->seq.pList7[0], x->Pacc7); SETFLOAT(&x->seq.pList7[1], x->E_Acc7);
  SETFLOAT(&x->seq.pList8[0], x->Pacc8); SETFLOAT(&x->seq.pList8[1], x->E_Acc8);
  if(x->myBug == 4) post("P2 = %f, E2 = %f, Location = %d",atom_getfloatarg(0,2,x->seq.pList2),atom_getfloatarg(1,2,x->seq.pList2),x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep);
  x->PJoin = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eJoin);
  //assignment of PJoined happens here, and then the value is manipulated in perform. See "FLAGS"
  //  if(x->PJoin > 1) x->PJoined = x->PJoin; // see below
  x->Gnm = (int)atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.groupNum);
  //  if(x->Gnm != x->PrevG) x->PStepOff = 0;// see below
  x->Gstep = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.groupStep);
  // trying this in perform, since it now inhabits a signal outlet:
  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eOff);
  x->PESize = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eSize);
  x->PESInv = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eSizeInv);
  x->Gn = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE + x->Gnm, x->VGROUPSIZE, x->vGrp.n);
  x->Gd = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE + x->Gnm, x->VGROUPSIZE, x->vGrp.d);
  x->GSize = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE + x->Gnm, x->VGROUPSIZE, x->vGrp.size);
  x->GSInv = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE + x->Gnm, x->VGROUPSIZE, x->vGrp.sizeInv);
  x->Goff = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE + x->Gnm, x->VGROUPSIZE, x->vGrp.offset);
  x->Grem = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE, x->VGROUPSIZE, x->vGrp.remains);
  //if(x->GSize > 0) x->GSInv = 1 / x->GSize;
  x->cycles = x->vGrp.cycles[x->slot + x->varPerf * SLOTS];
  x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.varOff);
  x->VOffG = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.grpOff);
  x->eChanged = 0;
}

void polyMath_tilde_pageTurn(t_polyMath_tilde *x)
{
  outlet_float(x->page, (t_float)x->pageNum);
}

void polyMath_tilde_eChange(t_polyMath_tilde *x)
{
  int varSlot;
  outlet_float(x->eAlt, (t_float)!x->altNum);
  if(x->changeSlot)
    {
      if(x->PStep >= x->seq.len[x->slot] - 1)
	{
	  x->nextSlotVal = x->nextSlot * MAXSEQ;
	  varSlot = 0;
	}
      else
	{
	  x->nextSlotVal - x->nextSlot * MAXSEQ + x->PStep + 1;
	  varSlot = 0;
	}
    }
  else if(x->changeVar)
    {
      if(x->PStep >= x->var.len[x->slot + x->thisVar * SLOTS])
	{
	  x->nextSlotVal = x->nextSlot * MAXSEQ + x->varPerf * x->SEQSIZE;
	  //x->eOut = atom_getfloatarg(x->nextSlot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.eAcc1);
	  varSlot = 1;
	}
      else
	{
	  x->nextSlotVal = x->nextSlot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep;
	  //x->eOut = atom_getfloatarg(x->nextSlot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc1);
	  varSlot = 1;
	}
    } 
  else if(!x->scrambling)
    {
      if(x->PStep >= x->seq.len[x->slot] - 1)
	{
	  x->nextSlotVal = x->slot * MAXSEQ;
	  //x->eOut = atom_getfloatarg(x->slot * MAXSEQ, x->SEQSIZE, x->seq.eAcc1);
	  varSlot = 0;
	}
      else
	{
	  x->nextSlotVal = x->slot * MAXSEQ + x->PStep + 1;
	  //x->eOut = atom_getfloatarg(x->slot * MAXSEQ + x->PStep + 1, x->SEQSIZE, x->seq.eAcc1);
	  varSlot = 0;
	}
    }
  else
    {
      if(x->PStep >= x->var.len[x->slot + x->thisVar * SLOTS])
	{
	  x->nextSlotVal = x->slot * MAXSEQ + x->varPerf * x->SEQSIZE;
	  //x->eOut = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.eAcc1);
	  varSlot = 1;
	}
      else
	{
	  x->nextSlotVal = x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep;
	  //x->eOut = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eAcc1);
	  varSlot = 1;
	}
    }
  //x->PESize = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eSize);
  //x->PESize = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eSize);
  if(varSlot == 0)
    {
      x->dPhase = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eSize);
      if(x->altNum == 0)
	{
	  x->dur2 = x->barBeat * x->dPhase;
	  SETFLOAT(&x->dList[0], x->dur2);
	  SETFLOAT(&x->dList[1], x->dPhase);
	  outlet_list(x->durAlt, gensym("list"), 2, x->dList);
	}
      else
	{
	  x->dur1 = x->barBeat * x->dPhase;
	  SETFLOAT(&x->dList[0], x->dur1);
	  SETFLOAT(&x->dList[1], x->dPhase);
	  outlet_list(x->durFirst, gensym("list"), 2, x->dList);
	}
      for(x->y = 0; x->y < x->altLen; x->y++)
	{
	  x->z = x->y * 2;
	  //switch(x->y)
	  //{
	  if(x->y == 0)
	    {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc1);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc1);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	    }
	    else if(x->y == 1)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc2);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc2);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 2)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc3);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc3);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 3)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc4);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc4);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 4)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc5);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc5);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 5)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc6);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc6);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 6)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc7);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc7);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 7)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.pAcc8);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->SEQSIZE, x->seq.eAcc8);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	      //default:
	      //break;
	      //}
	}
    }
  else
    {
      x->dPhase = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eSize);
      if(x->altNum == 0)
	{
	  x->dur2 = x->barBeat * x->dPhase;
	  SETFLOAT(&x->dList[0], x->dur2);
	  SETFLOAT(&x->dList[1], x->dPhase);
	  outlet_list(x->durAlt, gensym("list"), 2, x->dList);
	}
      else
	{
	  x->dur1 = x->barBeat * x->dPhase;
	  SETFLOAT(&x->dList[0], x->dur1);
	  SETFLOAT(&x->dList[1], x->dPhase);
	  outlet_list(x->durFirst, gensym("list"), 2, x->dList);
	}
      for(x->y = 0; x->y < x->altLen; x->y++)
	{
	  x->z = x->y * 2;
	  //switch(x->y)
	  //{
	    if(x->y == 0)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc1);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc1);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 1)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc2);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc2);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 2)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc3);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc3);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 3)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc4);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc4);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 4)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc5);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc5);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 5)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc6);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc6);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 6)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc7);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc7);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	    else if(x->y == 7)
	      {
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.pAcc8);
	      SETFLOAT(&x->outList[x->z], x->eOut);
	      x->eOut = atom_getfloatarg(x->nextSlotVal, x->VARSIZE, x->var.eAcc8);
	      SETFLOAT(&x->outList[x->z+1], x->eOut);
	      //break;
	      }
	      //default:
	      //break;
	      //}
	}
    }
  //outlet_float(x->eChange, x->eOut);
  outlet_list(x->eChange, gensym("list"), x->altLen * 2, x->outList);
}

void polyMath_tilde_preOut(t_polyMath_tilde *x, t_floatarg f)
{
  x->percentVal = f < 1 ? 0.1 : f > 99 ? 0.99 : f * 0.01;
}

void polyMath_tilde_altEarly(t_polyMath_tilde *x, t_floatarg f)
{
  x->altEarly = f > 0 ? 1 : 0;
}

void polyMath_tilde_cout(t_polyMath_tilde *x)
{
  SETFLOAT(&x->outList[0], (t_float)x->slot);
  SETFLOAT(&x->outList[1], (t_float)x->variation);
  outlet_float(x->dType, 99);
  outlet_list(x->dataOut, gensym("list"), 2, x->outList);
  /* if(x->jFirst > 0) outlet_float(x->eventLengthNum, (t_float)x->JoinVal);
  else if(x->PJoined <= 0) outlet_float(x->eventLengthNum, 1);
  if(x->jFirst > 0) outlet_float(x->eventLengthPhase, x->JESize);
  else if(x->PJoined <= 0) outlet_float(x->eventLengthPhase, x->PESize); */
  outlet_float(x->alt, (t_float)x->altNum);
  if(x->PJoined > 0 && x->jFirst == 1) outlet_float(x->eventLengthNum, (t_float)x->JoinVal);
  else if(x->jFlag == 0) outlet_float(x->eventLengthNum, 1.0f);
  if(x->PJoined > 0 && x->jFirst == 1) outlet_float(x->eventLengthPhase, x->JESize);
  else if(x->jFlag == 0) outlet_float(x->eventLengthPhase, x->PESize);
  outlet_float(x->denom, x->Gd);
  outlet_float(x->num, x->Gn);
  outlet_float(x->groupnum, (t_float)x->Gnm);
  outlet_list(x->p8, gensym("list"), 2, x->seq.pList8);
  outlet_list(x->p7, gensym("list"), 2, x->seq.pList7);
  outlet_list(x->p6, gensym("list"), 2, x->seq.pList6);
  outlet_list(x->p5, gensym("list"), 2, x->seq.pList5);
  outlet_list(x->p4, gensym("list"), 2, x->seq.pList4);
  outlet_list(x->p3, gensym("list"), 2, x->seq.pList3);
  outlet_list(x->p2, gensym("list"), 2, x->seq.pList2);
  outlet_list(x->p1, gensym("list"), 2, x->seq.pList1);
  if(x->barNew > 0)
    {
      outlet_bang(x->newbar);
      x->barNew = 0;
    }
  if(x->Gstep == 0)
    {
      outlet_bang(x->newgroup);
    }
  outlet_float(x->cycle, (t_float)x->cycles);
  outlet_float(x->subclock, x->Gstep);
  outlet_float(x->clock, (t_float)x->PStep + x->PStepOff);
}

int writeGroup(t_polyMath_tilde *x, int group)
{
  //moving this back to writeGroup
  x->Wsize = x->Gn / x->Gd;
  if(x->myBug > 0) post("Wsize = %f",x->Wsize);
  SETFLOAT(&x->grp.n[x->slot * GROUPS + group],x->Gn);
  SETFLOAT(&x->grp.d[x->slot * GROUPS + group],x->Gd);
  SETFLOAT(&x->grp.size[x->slot * GROUPS + group],x->Wsize);
  SETFLOAT(&x->grp.sizeInv[x->slot * GROUPS + group],1/x->Wsize);
  SETFLOAT(&x->grp.offset[x->slot * GROUPS + group],x->groupOffset);
  //end trying this in setGroups...
  x->WESize = x->Wsize;// / x->Gn; //whu? 
  x->getD = atom_getfloatarg(x->slot * GROUPS + group, x->GROUPSIZE, x->grp.d);
  //x->grp.gStart[group] = x->Gstart; // trying in setGroups instead
  if(atom_getfloatarg(x->slot * GROUPS + group, x->GROUPSIZE, x->grp.n) == 0 ||
     atom_getfloatarg(x->slot * GROUPS + group, x->GROUPSIZE, x->grp.d) == 0 ||
     atom_getfloatarg(x->slot * GROUPS + group, x->GROUPSIZE, x->grp.size) == 0)
    {
      //return(0);
      if(x->myBug == 3)
	{
	  post("n=%d, d=%d, size=%f", atom_getfloatarg(x->slot * GROUPS + group, x->GROUPSIZE, x->grp.n),
	       atom_getfloatarg(x->slot * GROUPS + group, x->GROUPSIZE, x->grp.d),
	       atom_getfloatarg(x->slot * GROUPS + group, x->GROUPSIZE, x->grp.size));
	}
      post("Exiting due to invalid entries = n, d or size == 0");
      return(3);
    }
  else
    {
      if(x->WESize <= 0)
	{
	  post("Exiting due to size <= 0");
	  return(2);
	}
      else
	{
	  x->WOffAcc = atom_getfloatarg(x->slot * MAXSEQ + x->Gstart - 1, x->SEQSIZE, x->seq.eOff) + atom_getfloatarg(x->slot * MAXSEQ + x->Gstart - 1, x->SEQSIZE, x->seq.eSize);
	  for(x->d = 0; x->d < (int)x->Gn; x->d++)
	    {
	      x->Wstep = x->seq.len[x->slot] + x->d;
	      if(x->myBug == 101) post("x->Wstep = %d",x->Wstep);
	      x->Woff = x->groupOffset + (x->WESize * (t_float)x->d) + x->WOffAcc;
	      x->WSInv = 1 / x->Wsize;
	      SETFLOAT(&x->seq.eSize[x->slot * MAXSEQ + x->Wstep], x->WESize);
	      SETFLOAT(&x->seq.eJoin[x->slot * MAXSEQ + x->Wstep], 0); // joins are set separately
	      SETFLOAT(&x->seq.jSize[x->slot * MAXSEQ + x->Wstep], x->WESize);
	      SETFLOAT(&x->seq.eSizeInv[x->slot * MAXSEQ + x->Wstep], 1 / x->WESize);
	      SETFLOAT(&x->seq.eOff[x->slot * MAXSEQ + x->Wstep], x->Woff);
	      SETFLOAT(&x->seq.allStep[x->slot * MAXSEQ + x->Wstep], (t_float)x->d + x->Gstart);
	      SETFLOAT(&x->seq.groupStep[x->slot * MAXSEQ + x->Wstep], (t_float)x->d);
	      SETFLOAT(&x->seq.groupNum[x->slot * MAXSEQ + x->Wstep], (t_float)group);
	      SETFLOAT(&x->seq.denom[x->slot * MAXSEQ + x->Wstep], x->getD);
	      if(x->myBug == 1) post("Step = %d, GStep = %d, WESize = %f, Woff = %f, Write: %d",x->seq.len[x->slot] + x->d, x->d, x->WESize, x->Woff, x->slot * MAXSEQ + x->Wstep);
	      if(x->myBug == 101) post("step: %d, gStep: %d, seq.eSize: %f, seq.eOff: %f",x->Wstep, (t_int)atom_getfloatarg(x->slot * MAXSEQ + x->Wstep, x->SEQSIZE, x->seq.groupStep), atom_getfloatarg(x->slot * MAXSEQ + x->Wstep, x->SEQSIZE, x->seq.eSize), atom_getfloatarg(x->slot * MAXSEQ + x->Wstep, x->SEQSIZE, x->seq.eOff));
	    }
	  x->seq.len[x->slot] += (int)x->Gn;
	  x->groupOffset += x->WESize * x->Gn;
	  return(1);
	}
    }
}

int reWriteSeq(t_polyMath_tilde *x)
{
  for(x->j = x->Gstp + (int)x->JJoin; x->j < x->seq.len[x->JSlot]; x->j++)
    {
      x->k = x->j - (int)x->JJoin;
      SETFLOAT(&x->seq.eSize[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eSize));
      SETFLOAT(&x->seq.eJoin[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eJoin));
      SETFLOAT(&x->seq.jSize[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.jSize));
      SETFLOAT(&x->seq.eSizeInv[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eSizeInv));
      SETFLOAT(&x->seq.eOff[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eOff));
      SETFLOAT(&x->seq.allStep[x->JSlot * MAXSEQ + x->k], (t_float)x->k);
      SETFLOAT(&x->seq.groupNum[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.groupNum));
      // next line should have JGt offset
      SETFLOAT(&x->seq.groupStep[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.groupStep) - x->JGt);
      SETFLOAT(&x->seq.eAcc1[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc1));
      SETFLOAT(&x->seq.eAcc2[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc2));
      SETFLOAT(&x->seq.eAcc3[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc3));
      SETFLOAT(&x->seq.eAcc4[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc4));
      SETFLOAT(&x->seq.pAcc1[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc1));
      SETFLOAT(&x->seq.pAcc2[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc2));
      SETFLOAT(&x->seq.pAcc3[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc3));
      SETFLOAT(&x->seq.pAcc4[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc4));
      SETFLOAT(&x->seq.eAcc5[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc5));
      SETFLOAT(&x->seq.eAcc6[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc6));
      SETFLOAT(&x->seq.eAcc7[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc7));
      SETFLOAT(&x->seq.eAcc8[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.eAcc8));
      SETFLOAT(&x->seq.pAcc5[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc5));
      SETFLOAT(&x->seq.pAcc6[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc6));
      SETFLOAT(&x->seq.pAcc7[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc7));
      SETFLOAT(&x->seq.pAcc8[x->JSlot * MAXSEQ + x->k], atom_getfloatarg(x->JSlot * MAXSEQ + x->j, x->SEQSIZE, x->seq.pAcc8));
    }
  return(1);
}

void polyMath_tilde_addGroup(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
 if(argc == 3)
   {
     x->slot = atom_getfloat(argv);
     x->slot = x->slot < 0 ? 0 : x->slot > 127 ? 127 : x->slot;
     x->maxGrp = x->grp.nGroups[x->slot];
     x->Gstart = x->grp.gStart[x->slot * GROUPS + x->maxGrp - 1];
     if(x->myBug == 1) post("x->maxGrp = %d, x->Gstart = %d",x->maxGrp,x->Gstart);
     x->Gstart += (int)atom_getfloatarg(x->slot * GROUPS + x->maxGrp - 1,x->GROUPSIZE, x->grp.n);
     if(x->myBug == 101) post("x->Gstart + x->grp.n = %d",x->Gstart);
     x->Gn = atom_getfloat(argv + 1);
     x->Gd = atom_getfloat(argv + 2);
     x->fGroup = x->grp.fillGroup[x->slot];
     if(x->fGroup > 0) // what does this do?
       {
	 x->c = x->fGroup;
       }
     else
       {
	 if(x->grp.isUnFilled[x->slot] == 1)
	   {
	     x->c = 0; // do we need to store this variable in an array?
	 // it strikes me that you can only add a group after a setGroups routine!
	   }
	 else x->c = x->maxGrp;
	 if(x->myBug == 101) post("x->c == %d, x->maxGrp = %d",x->c,x->maxGrp);
       }
     x->groupOffset = atom_getfloatarg(x->slot * GROUPS + x->maxGrp - 1, x->GROUPSIZE, x->grp.offset) + atom_getfloatarg(x->slot * GROUPS + x->maxGrp - 1, x->GROUPSIZE, x->grp.size);
     if(x->Gn <= 0 || x->Gd <= 0)
       {
	 post("Error: numerator and denominator must both be > 0");
	 SETFLOAT(&x->grp.remains[x->c + x->slot * GROUPS], x->Grem);
       }
     else 
       {
	 if(x->grp.isUnFilled[x->slot] == 1)
	   {
	     x->grp.isUnFilled[x->slot] = 0;
	   }
	 /*x->Wsize = x->Gn / x->Gd;
	 x->grp.nGroups[x->slot] = x->c + 1;
	 SETFLOAT(&x->grp.n[x->slot * GROUPS + x->c],x->Gn);
	 SETFLOAT(&x->grp.d[x->slot * GROUPS + x->c],x->Gd);
	 SETFLOAT(&x->grp.size[x->slot * GROUPS + x->c],x->Wsize);
	 SETFLOAT(&x->grp.sizeInv[x->slot * GROUPS + x->c],1/x->Wsize);
	 SETFLOAT(&x->grp.offset[x->slot * GROUPS + x->c],x->groupOffset); */
	      //end of trying this here - it didn't work
	 x->h = writeGroup(x,x->c);
	 if(x->h > 0)
	   {
	     x->grp.gStart[x->slot * GROUPS + x->c] = x->Gstart;
	     if(x->myBug == 1) post("group write exit code: %d",x->h);
	     x->c++;
	     if(x->myBug == 1) post("x->Gn = %d, x->Gstart = %d",(int)x->Gn,(int)x->Gstart);
	     x->Gstart += x->Gn;
	     x->Gcycle += (x->Gn / x->Gd);
	     x->grp.nGroups[x->slot]++;
	     if(x->myBug == 101) post("x->c = %d, x->grp.n = %f, x->grp.d = %f",x->c, atom_getfloatarg(x->slot * GROUPS + x->c, x->GROUPSIZE, x->grp.n), atom_getfloatarg(x->slot * GROUPS + x->c, x->GROUPSIZE, x->grp.d));
		  //x->seq.len[x->slot] += (int)x->Gn;
	   }
	 else post("group write unsuccessful");
       }
   }
 else post("Incorrect arguments to addGroup");
}

void polyMath_tilde_autoThreshold(t_polyMath_tilde *x, t_floatarg f)
{
  x->autoThreshold = f != 0 ? 1 : 0;
}

void polyMath_tilde_sizeThreshold(t_polyMath_tilde *x, t_floatarg f)
{
  x->sizeThreshold = f > 0.000001 ? f : 0.000001;
}

void polyMath_tilde_sizeFrac(t_polyMath_tilde *x, t_floatarg f)
{
  x->sizeFrac = f > 0.001 ? f : 0.5;
}

// the next two functions should be moved to the start of the code, and initSeqSlot should be
// called after the sequence is initialized thus: eventSeq_tilde_initSeqSlot(x, x->slot, 1);
// The variation writing functions should have the values for seq.filled copied into var.filled
int swapEventList(t_polyMath_tilde *x, int location, int slot, int var)
{
  x->isSwapList = 0;
  int varSeq = var - 1;
  if(var > 0)
    {
      SETFLOAT(&x->eventList[0], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.allStep));
      SETFLOAT(&x->eventList[1], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.filled));
      SETFLOAT(&x->eventList[2], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.groupStep));
      SETFLOAT(&x->eventList[3], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.groupNum));
      SETFLOAT(&x->eventList[4], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eSize));
      SETFLOAT(&x->eventList[5], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eOff));
      SETFLOAT(&x->eventList[6], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eJoin));
      SETFLOAT(&x->eventList[7], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.jSize));
      SETFLOAT(&x->eventList[8], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc1));
      SETFLOAT(&x->eventList[9], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc2));
      SETFLOAT(&x->eventList[10], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc3));
      SETFLOAT(&x->eventList[11], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc4));
      SETFLOAT(&x->eventList[12], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc5));
      SETFLOAT(&x->eventList[13], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc6));
      SETFLOAT(&x->eventList[14], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc7));
      SETFLOAT(&x->eventList[15], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eAcc8));
      SETFLOAT(&x->eventList[16], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc1));
      SETFLOAT(&x->eventList[17], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc2));
      SETFLOAT(&x->eventList[18], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc3));
      SETFLOAT(&x->eventList[19], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc4));
      SETFLOAT(&x->eventList[20], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc5));
      SETFLOAT(&x->eventList[21], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc6));
      SETFLOAT(&x->eventList[22], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc7));
      SETFLOAT(&x->eventList[23], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.pAcc8));
      SETFLOAT(&x->eventList[24], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.eSizeInv));
      SETFLOAT(&x->eventList[25], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.denom));
      SETFLOAT(&x->eventList[26], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.altOff));
      x->isSwapList = 1;
    }
  else
    {
      SETFLOAT(&x->eventList[0], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.allStep));
      SETFLOAT(&x->eventList[1], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.filled));
      SETFLOAT(&x->eventList[2], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.groupStep));
      SETFLOAT(&x->eventList[3], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.groupNum));
      SETFLOAT(&x->eventList[4], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eSize));
      SETFLOAT(&x->eventList[5], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eOff));
      SETFLOAT(&x->eventList[6], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eJoin));
      SETFLOAT(&x->eventList[7], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.jSize));
      SETFLOAT(&x->eventList[8], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc1));
      SETFLOAT(&x->eventList[9], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc2));
      SETFLOAT(&x->eventList[10], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc3));
      SETFLOAT(&x->eventList[11], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc4));
      SETFLOAT(&x->eventList[12], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc5));
      SETFLOAT(&x->eventList[13], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc6));
      SETFLOAT(&x->eventList[14], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc7));
      SETFLOAT(&x->eventList[15], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eAcc8));
      SETFLOAT(&x->eventList[16], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc1));
      SETFLOAT(&x->eventList[17], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc2));
      SETFLOAT(&x->eventList[18], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc3));
      SETFLOAT(&x->eventList[19], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc4));
      SETFLOAT(&x->eventList[20], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc5));
      SETFLOAT(&x->eventList[21], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc6));
      SETFLOAT(&x->eventList[22], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc7));
      SETFLOAT(&x->eventList[23], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.pAcc8));
      SETFLOAT(&x->eventList[24], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.eSizeInv));
      SETFLOAT(&x->eventList[25], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.denom));
      SETFLOAT(&x->eventList[26], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.altOff));
      x->isSwapList = 1;
    }
  return(x->isSwapList);
}

int addEventList(t_polyMath_tilde *x, int location, int len, int slot, int var)
{
  int isAdded = 0;
  int varSeq = var - 1;
  if(var > 0)
    {
      SETFLOAT(&x->var.allStep[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(0, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.filled[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(1, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.groupStep[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(2, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.groupNum[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(3, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eSize[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(4, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eOff[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(5, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eJoin[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(6, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.jSize[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(7, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc1[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(8, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc2[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(9, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc3[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(10, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc4[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(11, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc5[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(12, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc6[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(13, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc7[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(14, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eAcc8[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(15, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc1[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(16, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc2[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(17, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc3[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(18, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc4[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(19, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc5[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(20, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc6[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(21, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc7[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(22, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.pAcc8[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(23, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.eSizeInv[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(24, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.denom[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(25, EVENTLIST, x->eventList));
      SETFLOAT(&x->var.altOff[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(26, EVENTLIST, x->eventList));
      isAdded = 1;
    }
  else
    {
      SETFLOAT(&x->seq.allStep[slot * MAXSEQ + location], atom_getfloatarg(0, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.filled[slot * MAXSEQ + location], atom_getfloatarg(1, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.groupStep[slot * MAXSEQ + location], atom_getfloatarg(2, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.groupNum[slot * MAXSEQ + location], atom_getfloatarg(3, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eSize[slot * MAXSEQ + location], atom_getfloatarg(4, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eOff[slot * MAXSEQ + location], atom_getfloatarg(5, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eJoin[slot * MAXSEQ + location], atom_getfloatarg(6, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.jSize[slot * MAXSEQ + location], atom_getfloatarg(7, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc1[slot * MAXSEQ + location], atom_getfloatarg(8, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc2[slot * MAXSEQ + location], atom_getfloatarg(9, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc3[slot * MAXSEQ + location], atom_getfloatarg(10, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc4[slot * MAXSEQ + location], atom_getfloatarg(11, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc5[slot * MAXSEQ + location], atom_getfloatarg(12, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc6[slot * MAXSEQ + location], atom_getfloatarg(13, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc7[slot * MAXSEQ + location], atom_getfloatarg(14, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eAcc8[slot * MAXSEQ + location], atom_getfloatarg(15, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc1[slot * MAXSEQ + location], atom_getfloatarg(16, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc2[slot * MAXSEQ + location], atom_getfloatarg(17, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc3[slot * MAXSEQ + location], atom_getfloatarg(18, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc4[slot * MAXSEQ + location], atom_getfloatarg(19, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc5[slot * MAXSEQ + location], atom_getfloatarg(20, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc6[slot * MAXSEQ + location], atom_getfloatarg(21, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc7[slot * MAXSEQ + location], atom_getfloatarg(22, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.pAcc8[slot * MAXSEQ + location], atom_getfloatarg(23, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.eSizeInv[slot * MAXSEQ + location], atom_getfloatarg(24, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.denom[slot * MAXSEQ + location], atom_getfloatarg(25, EVENTLIST, x->eventList));
      SETFLOAT(&x->seq.altOff[slot * MAXSEQ + location], atom_getfloatarg(26, EVENTLIST, x->eventList));
      isAdded = 1;
    }
  return(isAdded);
}

int oneToTheRightOrLeft(t_polyMath_tilde *x, int location, int len, int slot, int var, int dir)
{ //tests are carried out in the parent function to make sure this is a valid swap
  int swapShuffle = 0;
  int location2 = location - 1;
  if(dir == -1) location2 = location + 1;
  int varSeq  = var - 1;
  if(var > 0)
    {
      if(location >= len)
	{
	  //case for moving to the end - or is that impossible?
	}
      else
	{
	  SETFLOAT(&x->var.allStep[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.allStep));
	  SETFLOAT(&x->var.filled[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.filled));
	  SETFLOAT(&x->var.groupStep[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.groupStep));
	  SETFLOAT(&x->var.groupNum[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.groupNum));
	  SETFLOAT(&x->var.eSize[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eSize));
	  SETFLOAT(&x->var.eOff[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eOff));
	  SETFLOAT(&x->var.eJoin[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eJoin));
	  SETFLOAT(&x->var.jSize[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.jSize));
	  SETFLOAT(&x->var.eAcc1[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc1));
	  SETFLOAT(&x->var.eAcc2[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc2));
	  SETFLOAT(&x->var.eAcc3[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc3));
	  SETFLOAT(&x->var.eAcc4[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc4));
	  SETFLOAT(&x->var.eAcc5[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc5));
	  SETFLOAT(&x->var.eAcc6[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc6));
	  SETFLOAT(&x->var.eAcc7[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc7));
	  SETFLOAT(&x->var.eAcc8[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc8));
	  SETFLOAT(&x->var.pAcc1[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc1));
	  SETFLOAT(&x->var.pAcc2[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc2));
	  SETFLOAT(&x->var.pAcc3[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc3));
	  SETFLOAT(&x->var.pAcc4[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc4));
	  SETFLOAT(&x->var.pAcc5[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc5));
	  SETFLOAT(&x->var.pAcc6[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc6));
	  SETFLOAT(&x->var.pAcc7[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc7));
	  SETFLOAT(&x->var.pAcc8[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc8));
	  SETFLOAT(&x->var.eSizeInv[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eSizeInv));
	  SETFLOAT(&x->var.denom[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.denom));
	  SETFLOAT(&x->var.altOff[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.altOff));
	  swapShuffle++;
	}
    }
  else
    {
      if(location == 0)
	{
	  //case for moving to the end, or is that impossible? 
	}
      else
	{
	  SETFLOAT(&x->seq.allStep[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.allStep));
	  SETFLOAT(&x->seq.filled[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.filled));
	  SETFLOAT(&x->seq.groupStep[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.groupStep));
	  SETFLOAT(&x->seq.groupNum[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.groupNum));
	  SETFLOAT(&x->seq.eSize[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eSize));
	  SETFLOAT(&x->seq.eOff[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eOff));
	  SETFLOAT(&x->seq.eJoin[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eJoin));
	  SETFLOAT(&x->seq.jSize[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.jSize));
	  SETFLOAT(&x->seq.eAcc1[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc1));
	  SETFLOAT(&x->seq.eAcc2[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc2));
	  SETFLOAT(&x->seq.eAcc3[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc3));
	  SETFLOAT(&x->seq.eAcc4[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc4));
	  SETFLOAT(&x->seq.eAcc5[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc5));
	  SETFLOAT(&x->seq.eAcc6[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc6));
	  SETFLOAT(&x->seq.eAcc7[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc7));
	  SETFLOAT(&x->seq.eAcc8[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc8));
	  SETFLOAT(&x->seq.pAcc1[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc1));
	  SETFLOAT(&x->seq.pAcc2[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc2));
	  SETFLOAT(&x->seq.pAcc3[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc3));
	  SETFLOAT(&x->seq.pAcc4[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc4));
	  SETFLOAT(&x->seq.pAcc5[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc5));
	  SETFLOAT(&x->seq.pAcc6[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc6));
	  SETFLOAT(&x->seq.pAcc7[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc7));
	  SETFLOAT(&x->seq.pAcc8[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc8));
	  SETFLOAT(&x->seq.eSizeInv[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eSizeInv));
	  SETFLOAT(&x->seq.denom[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.denom));
	  SETFLOAT(&x->seq.altOff[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.altOff));
	  swapShuffle++;
	}
    }
  return(swapShuffle);
}

int onePToTheLeftOrRight(t_polyMath_tilde *x, int location, int len, int slot, int var, int P, int direction)
{ //tests are carried out in the parent function to make sure this is a valid swap
  int swapShuffle = 0;
  int location2 = 0;
  if(direction == -1) location2 = location - 1;
  else location2 = location + 1;
  int varSeq  = var - 1;
  if(var > 0)
    {
      if(location >= len)
	{
	  //case for moving to the end - or is that impossible?
	}
      else
	{
	  switch(P)
	    {
	    case(1):
	      SETFLOAT(&x->var.pAcc1[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc1));
	      break;
	    case(2):
	      SETFLOAT(&x->var.pAcc2[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc2));
	      break;
	    case(3):
	      SETFLOAT(&x->var.pAcc3[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc3));
	      break;
	    case(4):
	      SETFLOAT(&x->var.pAcc4[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc4));
	      break;
	    case(5):
	      SETFLOAT(&x->var.pAcc5[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc5));
	      break;
	    case(6):
	      SETFLOAT(&x->var.pAcc6[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc6));
	      break;
	    case(7):
	      SETFLOAT(&x->var.pAcc7[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc7));
	      break;
	    case(8):
	      SETFLOAT(&x->var.pAcc8[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.pAcc8));
	      break;
	    default:
	      break;
	    }
	  swapShuffle++;
	}
    }
  else
    {
      if(location == 0)
	{
	  //case for moving to the end, or is that impossible? 
	}
      else
	{
	  switch(P)
	    {
	    case(1):
	      SETFLOAT(&x->seq.pAcc1[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc1));
	      break;
	    case(2):
	      SETFLOAT(&x->seq.pAcc2[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc2));
	      break;
	    case(3):
	      SETFLOAT(&x->seq.pAcc3[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc3));
	      break;
	    case(4):
	      SETFLOAT(&x->seq.pAcc4[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc4));
	      break;
	    case(5):
	      SETFLOAT(&x->seq.pAcc5[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc5));
	      break;
	    case(6):
	      SETFLOAT(&x->seq.pAcc6[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc6));
	      break;
	    case(7):
	      SETFLOAT(&x->seq.pAcc7[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc7));
	      break;
	    case(8):
	      SETFLOAT(&x->seq.pAcc8[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.pAcc8));
	      break;
	    default:
	      break;
	    }
	  swapShuffle++;
	}
    }
  return(swapShuffle);
}

int oneEToTheLeftOrRight(t_polyMath_tilde *x, int location, int len, int slot, int var, int P, int direction)
{ //tests are carried out in the parent function to make sure this is a valid swap
  int swapShuffle = 0;
  int location2 = 0;
  if(direction == -1) location2 = location - 1;
  else location2 = location + 1;
  int varSeq  = var - 1;
  if(var > 0)
    {
      if(location >= len)
	{
	  //case for moving to the end - or is that impossible?
	}
      else
	{
	  switch(P)
	    {
	    case(1):
	      SETFLOAT(&x->var.eAcc1[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc1));
	      break;
	    case(2):
	      SETFLOAT(&x->var.eAcc2[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc2));
	      break;
	    case(3):
	      SETFLOAT(&x->var.eAcc3[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc3));
	      break;
	    case(4):
	      SETFLOAT(&x->var.eAcc4[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc4));
	      break;
	    case(5):
	      SETFLOAT(&x->var.eAcc5[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc5));
	      break;
	    case(6):
	      SETFLOAT(&x->var.eAcc6[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc6));
	      break;
	    case(7):
	      SETFLOAT(&x->var.eAcc7[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc7));
	      break;
	    case(8):
	      SETFLOAT(&x->var.eAcc8[slot * MAXSEQ + varSeq * x->SEQSIZE + location], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location2, x->VARSIZE, x->var.eAcc8));
	      break;
	    default:
	      break;
	    }
	  swapShuffle++;
	}
    }
  else
    {
      if(location == 0)
	{
	  //case for moving to the end, or is that impossible? 
	}
      else
	{
	  switch(P)
	    {
	    case(1):
	      SETFLOAT(&x->seq.eAcc1[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc1));
	      break;
	    case(2):
	      SETFLOAT(&x->seq.eAcc2[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc2));
	      break;
	    case(3):
	      SETFLOAT(&x->seq.eAcc3[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc3));
	      break;
	    case(4):
	      SETFLOAT(&x->seq.eAcc4[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc4));
	      break;
	    case(5):
	      SETFLOAT(&x->seq.eAcc5[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc5));
	      break;
	    case(6):
	      SETFLOAT(&x->seq.eAcc6[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc6));
	      break;
	    case(7):
	      SETFLOAT(&x->seq.eAcc7[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc7));
	      break;
	    case(8):
	      SETFLOAT(&x->seq.eAcc8[slot * MAXSEQ + location], atom_getfloatarg(slot * MAXSEQ + location2, x->SEQSIZE, x->seq.eAcc8));
	      break;
	    default:
	      break;
	    }
	  swapShuffle++;
	}
    }
  return(swapShuffle);
}

void polyMath_tilde_swapElement(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  int process = 0;
  int shuffled = 0;
  if(argc == 2)
    {
      x->swapSlot = x->slot;
      x->swapVar = x->varPerf;
      x->swapLoc = (int)atom_getfloat(argv);
      x->swapShift = (int)atom_getfloat(argv+1);
      if(x->varTest > 0)
	{
	  x->swapLength = x->var.len[x->swapSlot + x->swapVar * SLOTS];
	  if(x->swapLength < x->swapLoc + x->swapShift)
	    {
	      post("You cannot shift an element beyond the end of the sequence! Length = %d, newLoc = %d",x->var.len[x->swapSlot + x->swapVar * SLOTS], x->swapLoc + x->swapShift);
	    }
//allStep filled groupStep groupNum eSize eOff eJoin jSize eAcc1-8 pAcc1-8 eSizeInv denom altOff
	  else
	    if(x->swapLoc >= 0 && x->swapLoc < x->swapLength)
	      {
		x->isSwapList = swapEventList(x, x->swapLoc, x->swapSlot, x->varTest);
		if(x->swapShift > 0)
		  {
		    for(x->x=x->swapLoc; x->x < x->swapLoc + x->swapShift; x->x++)
		      {
			//jiggle events
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, x->varTest, -1);
		      }
		    //int addEventList(t_polyMath_tilde *x, int location, int len, int slot, int var)
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, x->varTest);
		  }
		else if(x->swapShift < 0)
		  {
		    for(x->x = x->swapLoc; x->x > x->swapLoc + x->swapShift; x->x--)
		      {
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, x->varTest, 0);
		      }
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, x->varTest);
		  }
		if(!process)
		  {
		    post("Element was not moved successfully!");
		  }
		//if(shuffled == fabs(x->swapShift)
	      }
	    else
	      {
		post("You cannot move an element that doesn't exist, i.e. is beyond the sequence! Length = %d, Location = %d",x->var.len[x->swapSlot + x->swapVar * SLOTS], x->swapLoc);
	      }
	}
      else
	{
	  x->swapLength = x->seq.len[x->swapSlot];
	  if(x->swapLength < x->swapLoc + x->swapShift)
	    {
	      post("You cannot shift an element beyond the end of the sequence! Length = %d, newLoc = %d",x->var.len[x->swapSlot + x->swapVar * SLOTS], x->swapLoc + x->swapShift);
	    }
//allStep filled groupStep groupNum eSize eOff eJoin jSize eAcc1-8 pAcc1-8 eSizeInv denom altOff
	  else
	    if(x->swapLoc >= 0 && x->swapLoc < x->swapLength)
	      {
		x->isSwapList = swapEventList(x, x->swapLoc, x->swapSlot, 0);
		if(x->swapShift > 0)
		  {
		    for(x->x=x->swapLoc; x->x < x->swapLoc + x->swapShift; x->x++)
		      {
			//jiggle events
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, 0, -1);
		      }
		    //int addEventList(t_polyMath_tilde *x, int location, int len, int slot, int var)
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, x->varTest);
		  }
		else if(x->swapShift < 0)
		  {
		    for(x->x = x->swapLoc; x->x > x->swapLoc + x->swapShift; x->x--)
		      {
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, 0, 0);
		      }
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, 0);
		  }
		if(!process)
		  {
		    post("Element was not moved successfully!");
		  }
		//if(shuffled == fabs(x->swapShift)
	      }
	  //then deal with argc == 4 and argc == 6
	}
    }
  else if(argc == 4)
    {
      x->swapSlot = (int)atom_getfloat(argv);
      x->swapVar = (int)atom_getfloat(argv+1);
      x->swapLoc = (int)atom_getfloat(argv+2);
      x->swapShift = (int)atom_getfloat(argv+3);
      if(x->swapSlot < 0 || x->swapSlot > SLOTS)
	{
	  post("slot is out of range: %d",x->swapSlot);
	}
      else if(x->swapVar < 0 || x->swapVar > VARIATIONS)
	{
	  post("variation is out of range: %d",x->swapVar);
	}
      else if(x->swapVar > 0)
	{
	  x->swapVar--;
	  x->swapLength = x->var.len[x->swapSlot + x->swapVar * SLOTS];
	  if(x->swapLength < x->swapLoc + x->swapShift)
	    {
	      post("You cannot shift an element beyond the end of the sequence! Length = %d, newLoc = %d",x->var.len[x->swapSlot + x->swapVar * SLOTS], x->swapLoc + x->swapShift);
	    }
//allStep filled groupStep groupNum eSize eOff eJoin jSize eAcc1-8 pAcc1-8 eSizeInv denom altOff
	  else
	    if(x->swapLoc >= 0 && x->swapLoc < x->swapLength)
	      {
		x->isSwapList = swapEventList(x, x->swapLoc, x->swapSlot, x->varTest);
		if(x->swapShift > 0)
		  {
		    for(x->x=x->swapLoc; x->x < x->swapLoc + x->swapShift; x->x++)
		      {
			//jiggle events
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, x->varTest, -1);
		      }
		    //int addEventList(t_polyMath_tilde *x, int location, int len, int slot, int var)
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, x->varTest);
		  }
		else if(x->swapShift < 0)
		  {
		    for(x->x = x->swapLoc; x->x > x->swapLoc + x->swapShift; x->x--)
		      {
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, x->varTest, 0);
		      }
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, x->varTest);
		  }
		if(!process)
		  {
		    post("Element was not moved successfully!");
		  }
		//if(shuffled == fabs(x->swapShift)
	      }
	    else
	      {
		post("You cannot move an element that doesn't exist, i.e. is beyond the sequence! Length = %d, Location = %d",x->var.len[x->swapSlot + x->swapVar * SLOTS], x->swapLoc);
	      }
	}
      else
	{
	  x->swapLength = x->seq.len[x->swapSlot];
	  if(x->swapLength < x->swapLoc + x->swapShift)
	    {
	      post("You cannot shift an element beyond the end of the sequence! Length = %d, newLoc = %d",x->var.len[x->swapSlot + x->swapVar * SLOTS], x->swapLoc + x->swapShift);
	    }
//allStep filled groupStep groupNum eSize eOff eJoin jSize eAcc1-8 pAcc1-8 eSizeInv denom altOff
	  else
	    if(x->swapLoc >= 0 && x->swapLoc < x->swapLength)
	      {
		x->isSwapList = swapEventList(x, x->swapLoc, x->swapSlot, 0);
		if(x->swapShift > 0)
		  {
		    for(x->x=x->swapLoc; x->x < x->swapLoc + x->swapShift; x->x++)
		      {
			//jiggle events
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, 0, -1);
		      }
		    //int addEventList(t_polyMath_tilde *x, int location, int len, int slot, int var)
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, x->varTest);
		  }
		else if(x->swapShift < 0)
		  {
		    for(x->x = x->swapLoc; x->x > x->swapLoc + x->swapShift; x->x--)
		      {
			shuffled += oneToTheRightOrLeft(x, x->x, x->swapLength, x->swapSlot, 0, 0);
		      }
		    process = addEventList(x, x->swapLoc + x->swapShift, x->swapLength, x->swapSlot, 0);
		  }
		if(!process)
		  {
		    post("Element was not moved successfully!");
		  }
		//if(shuffled == fabs(x->swapShift)
	      }
	  //then deal with argc == 4 and argc == 6
	}
    }
  else if(argc == 6)
    {
      //here we will allow the swapping of only the pAcc and eAcc elements
      x->swapSlot = (int)atom_getfloat(argv);
      x->swapVar = (int)atom_getfloat(argv+1);
      x->swapLoc = (int)atom_getfloat(argv+2);
      x->swapShift = (int)atom_getfloat(argv+3);
      x->swapP = atom_getfloat(argv+4);
      x->swapE = atom_getfloat(argv+5);

    }
  //this is a placeholder for the function that will re-organise sequences
}

void polyMath_tilde_initSeqSlot(t_polyMath_tilde *x, t_floatarg newSeqSlot, t_floatarg isSeq)
{
  int seqSlot = (int)newSeqSlot;
  int seqIs = (int)isSeq;
  if(newSeqSlot < SLOTS && newSeqSlot >= 0)
    {
      if(isSeq == 0)
	{
	  for(x->w=0;x->w < MAXSEQ;x->w++)
	    {
	      SETFLOAT(&x->seq.filled[seqSlot * MAXSEQ + x->w],0);
	    }
	}
      else if(isSeq == 1)
	{
	  for(x->w = 0; x->w < x->seq.len[seqSlot]; x->w++)
	    {
	      SETFLOAT(&x->seq.filled[seqSlot * MAXSEQ + x->w],1);
	    }
	  for(x->w = x->seq.len[seqSlot]; x->w < MAXSEQ; x->w++)
	    {
	      SETFLOAT(&x->seq.filled[seqSlot * MAXSEQ + x->w],0);
	    }
	}
    }
  else
    {
      post("Sequence slots must be from 0 to %d",SLOTS);
    }
}

/* New group type - sequence type group. New samples can be loaded into a sequence buffer, or a single sample can be used as the buffer.
 * Events are positioned in the sequence as a timeline, with gaps where line output is static (!alt) or reset to 0 (alt) between.
 * Event lengths are still measured in phase-fractions, but in this case they may be any value. 
 * Wouldn't this be easier if we built another object? But that would require transfer of data between this one if converting
 * from this to the other...
 * ...this would have to happen within this object anyway, so let's do it...but it requires an extra if statement to be evaluated
 * sample-by-sample. Another object would be cheaper. So how to put two perform functions into one object?
 * 2018 - 11 - 13 - I'll work it out!
 */
void polyMath_tilde_seqInSlot(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  x->seqGrpOffset = x->slot * GROUPS;
  x->seqSlotOffset = x->slot * MAXSEQ;
  // incomplete...
  /* type 0 sequences are tuple-based (numeral, denominator) and do not have any non-integer-ratio gaps or elements
   * type 1 sequences are linear time, and elements do not have to conform to integer-ratios. Gaps or sections of audio between  
   * one event and the next are interpreted (ratios of 1/p where p = phase are configured)
   * we need to be able to insert samples into a sequence as well as segmenting long audio files.
   * 
   */
  if(argc == 4)
    {
      x->sType = (int)atom_getfloat(argv+3);
      if(x->sType == 0)
	{
	  //post("Each event needs 3 elements: position, numeral, tuple!");
	  //OPTION / METHOD idea 1:
	  /* Normally a sort of step entry inserts tuples one-by-one, with the position incrementing each time.
	   * If an element is inserted later on in the sequence, the gap is divided up evenly into tuples...
	   * These may be set differently later on, with say a gap of 3/6:
	   * We set the first element to be 1/12 and the remaining two elements become 5/24 long each.
	   * An internal state is needed to record the edit status of each element
	   */
	  x->seqPos = (int)atom_getfloat(argv);
	  x->seqNum = (int)atom_getfloat(argv+1);
	  x->seqDen = (int)atom_getfloat(argv+2);
	  if(x->seqNum < 1 || x->seqDen < 1)
	    {
	      post("Numerator and denominator must be >= 1!");
	    }
	  else if(x->seqPos < 0 || x->seqPos >= MAXSEQ)
	    {
	      post("Position must be an integer from 0 thru %d!",MAXSEQ - 1);
	    }
	  else
	    {
	      if(x->seqPos > x->seq.len[x->slot])
		{
		  //the decision is made thus: if the position is after the end of the sequence, an
		  //event (or more than one) is generated from the previous event's settings,
		  //but with a 'silent' flag (one of the eAcc/pAcc variables)
		  //x->prevSNum = x->grp.n[x->grp.nGroups[x->slot] - 1];
		  x->prevSDen = (int)atom_getfloatarg(x->grp.nGroups[x->slot] - 1, SLOTS, x->grp.d);
		  int seqDiff = x->seqPos - x->seq.len[x->slot];
		  x->seqPhase = (t_float)x->seqNum / (t_float)x->seqDen;
		  //for(x=0;x<seqDiff;x++)
		  //{
		  //  write values of previous seq (silent)
		  //}
		}
	      else if(x->seqPos == x->seq.len[x->slot])
		{
		  x->seqPOff = atom_getfloatarg(x->seqPos - 1 + x->slot * MAXSEQ, x->SEQSIZE, x->seq.eOff) + atom_getfloatarg(x->seqPos - 1 + x->slot * MAXSEQ, x->SEQSIZE, x->seq.eSize);
		  x->seqPhase = (t_float)x->seqNum / (t_float)x->seqDen;
		}
	      else
		{

		  x->seqPhase = (t_float)x->seqNum / (t_float)x->seqDen;
		}
	    }
	}
      else if(x->sType == 1)
	{
	  x->seqPos = (int)atom_getfloat(argv);
	  x->seqNum = (int)atom_getfloat(argv+1);
	  x->seqDen = (int)atom_getfloat(argv+2);
	  
	}
    }
}

/* When moving, the graphics need to respond separately, so that you see the resulting sequence before 
 * you commit (e.g. when you move an element 1 to the right, the element that was there moves to the left)
 */

/*typedef struct _vars
{
  int nGroups[SLOTS * VARIATIONS];                // number of groups in this sequence
  int gStart[GROUPS * SLOTS * VARIATIONS];        // where in the sequence does each group start?
  t_atom n[GROUPS * SLOTS * VARIATIONS];          // numerator of the time sig (fraction)
  t_atom d[GROUPS * SLOTS * VARIATIONS];          // denominator of the time sig
  int cycles[SLOTS * VARIATIONS];                 // number of cycles of phasor~
  t_atom offset[GROUPS * SLOTS * VARIATIONS];     // group offset in phase
  t_atom size[GROUPS * SLOTS * VARIATIONS];       // group size in phase
  t_atom sizeInv[GROUPS * SLOTS * VARIATIONS];    // 1 / size - so that look-up table can be used instead of math at runtime
  //t_atom rLength[GROUPS * SLOTS * VARIATIONS];    // real length of sequence calculated to foat precision
  t_atom remains[GROUPS * SLOTS * VARIATIONS];    // how much of the phase * cycles is left (in case of a 0 entry - skip one and if still 0 use rem
                                     // also, if an incomplete list is issued, the last value will be based on this!
  int fillGroup[SLOTS];              // if a fill group is used for an incomplete setGroups call, the index will be stored here
  int swaps[MAXSEQ];
  int swapsRef[MAXSEQ * 2];
  int swapped[MAXSEQ];*/
void polyMath_tilde_groupScramble(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  int nGroups;
  t_float fNGroups;
  //args: slot, var, destVar, scramRand, (mode CHANGED BEHAVIOUR)
  //if mode == 0, just do that percentage of scrambles
  //if mode == 1, ensure exactly that number of unique displacements is carried out (NOT IMPLEMENTED)
  //if mode == 2, scramble all of them, ensuring that none end up in the same place (order) as before
  //defaults:
  x->GSScramRand = 0.5;

  if(argc == 4)
    {
      x->GSSlot = (t_int)atom_getfloat(argv+0);
      x->GSVar = (t_int)atom_getfloat(argv+1);
      x->GSDestVar = (t_int)atom_getfloat(argv+2);
      x->GSScramRand = atom_getfloat(argv+3);
      x->GSMode = 0;
    }
  else if(argc == 3)
    {
      x->GSSlot = (t_int)atom_getfloat(argv+0);
      x->GSVar = (t_int)atom_getfloat(argv+1);
      x->GSDestVar = (t_int)atom_getfloat(argv+2);
      //x->GSScramRand = atom_getfloat(argv+3);
      x->GSMode = 2;
    }
      //  else if(argc == 4)
      // {
      // }
  if(x->GSSlot >= SLOTS || x->GSSlot < 0) post("slot must be a whole number from 0 to %d",SLOTS - 1);
  else if(x->GSVar < 0 || x->GSVar > VARIATIONS) post("var must bo 0 (no-var) or a whole number from 1 to %d",VARIATIONS);
  else if(x->GSDestVar < 0 || x->GSDestVar > VARIATIONS) post("dest var must bo 0 (no-var) or a whole number from 1 to %d",VARIATIONS);
  else if(x->GSScramRand < 0 || x->GSScramRand > 1) post("randomness must be a floating point number from 0 to 1");
  else
    {
      if(x->GSMode == 0)
	{
	  /*x->randNum1 = drand48();
	    x->fSwapsNum = x->fSeqLen * prob;
	    x->swapsNum = rounder(x,x->fSwapsNum,x->seqLen - 1);
	    x->doSwaps = x->swapsNum;*/
	  if(x->GSVar == 0)
	    {
	      
	    }
	  post("groupScramble Mode 0: NOT YET IMPLEMENTED!");
	}
      else if(x->GSMode == 1)
	{
	  post("groupScramble Mode 1: NOT YET IMPLEMENTED!");
	}
      else if(x->GSMode == 2)
	{
	  if(x->GSVar > 0) x->swapsNum = x->vGrp.nGroups[x->GSSlot + x->GSVar * SLOTS];
	  else x->swapsNum = x->grp.nGroups[x->GSSlot];
	  x->fSwapsNum = (t_float)x->swapsNum;
	  for(x->o = 0; x->o < x->swapsNum; x->o++)
	    {
	      x->vGrp.groupSwaps[x->o] = -1;
	    }
	  while(x->swapsNum)
	    {
	      x->o = 0;
	      x->randNum1 = drand48();
	      x->randNum2 = drand48();
	      x->GSFSwap = x->randNum1 * x->fSwapsNum;
	      x->GSISwap = (t_int)x->GSFSwap;
	      if(x->o != x->GSISwap)
		{
		  x->vGrp.groupSwaps[x->GSISwap] = x->o;
		  x->o++;
		  x->swapsNum--;
		}
	      else if(x->swapsNum == 1 && x->o == x->GSISwap)
		{
		  //kludge? recursion? ignore? ignore!
		  x->vGrp.groupSwaps[x->GSISwap] = x->o;
		  //x->o++;
		  x->swapsNum--;		  
		}
	    }
	  //t_int scramSuccess = scramGroup(x, x->GSSlot, x->
	}
    }
  //else post("Insufficient arguments to groupScramble!");  
}

/*  int varSeq = var - 1;
  if(var > 0)
    {
    SETFLOAT(&x->eventList[0], atom_getfloatarg(slot * MAXSEQ + varSeq * x->SEQSIZE + location, x->VARSIZE, x->var.allStep));
  else
    {
      SETFLOAT(&x->eventList[0], atom_getfloatarg(slot * MAXSEQ + location, x->SEQSIZE, x->seq.allStep));

  x->Gn = atom_getfloatarg(x->slot * GROUPS + x->Gnm, x->GROUPSIZE, x->grp.n);

  x->Gn = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE + x->Gnm, x->VGROUPSIZE, x->vGrp.n);

*/

t_int scramGroup(t_polyMath_tilde *x, t_int slot, t_int sVar, t_int dVar, t_int groupsNum)
{
  t_int isVar = sVar > 0 ? 1 : 0;
  t_int sourceOffset = isVar ? slot * GROUPS + sVar * x->GROUPSIZE : slot * GROUPS;
  t_int destinOffset = slot * GROUPS + dVar * x->GROUPSIZE;
  t_int seqOffset = isVar ? slot * MAXSEQ + sVar * x->SEQSIZE : slot * MAXSEQ; //x->VARSIZE
  t_int varOffset = slot * MAXSEQ + dVar * x->SEQSIZE; //x->SEQSIZE
}


void polyMath_tilde_groupInSlot(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  x->grp.gType[x->slot] = 0;
  int groupOffset = x->slot * GROUPS;
  int slotOffset = x->slot * MAXSEQ;
  int i, mark;
  x->s = 0;
  x->Woff = 0;
  mark = 0;
  x->grp.nGroups[x->slot] = 0;
  x->seq.len[x->slot] = 0;
  for(x->c = 0; x->c < argc; x->c += 2)
    {
      x->Gn = atom_getfloat(argv + x->c);
      x->Gd = atom_getfloat(argv + x->c + 1);
      //post("Group in slot: %f / %f",x->Gn,x->Gd);
      //post("slot = %d",x->slot);
      if(x->Gd > 0.0f && x->Gn > 0.0f)
	{
	  x->r = x->c / 2;
	  SETFLOAT(&x->grp.n[groupOffset + x->r],x->Gn);
	  SETFLOAT(&x->grp.d[groupOffset + x->r],x->Gd);
	  x->halfSize = (1 / x->Gd) * x->sizeFrac;
	  if(x->autoThreshold && x->halfSize < x->sizeThreshold) x->sizeThreshold = x->halfSize; 
	}
      else post("values not greater than 0!");
      x->seq.len[x->slot] += (int)x->Gn;
      x->grp.gStart[groupOffset + x->r] = mark;
      x->t = 0;
      x->GSize = 0;
      x->WESize = 1 / x->Gd;
      x->WSInv = 1 / x->WESize;
      //post("WESize = %f, WSI = %f",x->WESize,x->WSInv);
      //for(x->s; x->s < x->s + (int)x->Gn; x->s++)
      for(x->s = 0;x->s < (int)x->Gn; x->s++)
	{
	  if(x->Gd > 0.0f && x->Gd > 0.0f)
	    {
	      SETFLOAT(&x->seq.eSize[slotOffset + mark + x->s], x->WESize);
	      SETFLOAT(&x->seq.eOff[slotOffset + mark + x->s], x->Woff);
	      SETFLOAT(&x->seq.eSizeInv[slotOffset + mark + x->s], x->WSInv);
	      SETFLOAT(&x->seq.denom[slotOffset + mark + x->s], x->Gd);

	      SETFLOAT(&x->seq.allStep[slotOffset + mark + x->s], (t_float)mark + (t_float)x->s);
	      SETFLOAT(&x->seq.groupStep[slotOffset + mark + x->s], (t_float)x->t);
	      SETFLOAT(&x->seq.groupNum[slotOffset + mark + x->s], (t_float)x->r);
	      SETFLOAT(&x->seq.eJoin[slotOffset + mark + x->s], 1);
	      SETFLOAT(&x->seq.jSize[slotOffset + mark + x->s], x->WESize);
	      x->GSize += x->WESize;
	      x->Woff += x->WESize;
	      x->t++;
		  
	    }
	  else post("You can't have size == 0.000000");
	}
      mark += x->s;
      SETFLOAT(&x->grp.offset[groupOffset + x->r],x->Goff);
      SETFLOAT(&x->grp.size[groupOffset + x->r],x->GSize);
      SETFLOAT(&x->grp.sizeInv[groupOffset + x->r],1 / x->GSize);
	  
      x->Goff = x->Woff;
      x->grp.nGroups[x->slot]++;
      x->s++;
    }
  x->Icycle = (int)x->Goff;
  x->Gcycle = x->Goff;
  x->cycleDiff = x->Gcycle - (t_float)x->Icycle;
  if(x->cycleDiff > x->sizeThreshold && 1 - x->cycleDiff > x->sizeThreshold)
    {
      if(x->myBug == 14) post("x->Gcycle - (t_float)x->Icycle = %f",x->Gcycle - (t_float)x->Icycle);
      if(x->myBug == 10 || x->myBug == 14) post("Gcycle = %f, (t_float)Icycle = %f",x->Gcycle,(t_float)x->Icycle);
      x->Icycle++;
      x->Grem = (t_float)x->Icycle - x->Gcycle;
      if(x->myBug == 10 || x->myBug == 14)
	{
	  post("START");
	  post("x->Goff = %f",x->Goff);
	  post("x->Grem = %f",x->Grem);
	  post("");
	}
      x->Gd = 1 / x->Grem;
      x->Gn = 1;
      SETFLOAT(&x->grp.n[groupOffset + x->r + 1],x->Gn);
      SETFLOAT(&x->grp.d[groupOffset + x->r + 1],x->Gd);
      SETFLOAT(&x->grp.size[groupOffset + x->r + 1],x->Grem);
      SETFLOAT(&x->grp.sizeInv[groupOffset + x->r + 1],x->Gd);
      SETFLOAT(&x->grp.offset[groupOffset + x->r + 1],x->Goff);
      x->grp.nGroups[x->slot]++;
      x->grp.gStart[groupOffset + x->r + 1] = mark;
      
      SETFLOAT(&x->seq.eSize[slotOffset + mark], x->Grem);
      SETFLOAT(&x->seq.eOff[slotOffset + mark], x->Goff);
      SETFLOAT(&x->seq.eSizeInv[slotOffset + mark], x->Gd);
      SETFLOAT(&x->seq.denom[slotOffset + mark], x->Gd);

      SETFLOAT(&x->seq.allStep[slotOffset + mark], (t_float)mark);
      SETFLOAT(&x->seq.groupStep[slotOffset + mark], 0);
      SETFLOAT(&x->seq.groupNum[slotOffset + mark], (t_float)x->r + 1);
      SETFLOAT(&x->seq.eJoin[slotOffset + mark], 1);
      SETFLOAT(&x->seq.jSize[slotOffset + mark], x->Grem);
      mark++;
    }
  else if(x->myBug == 14) post("sizeThreshold = %f, difference = either %f or %f",x->sizeThreshold,x->cycleDiff, 1 - x->cycleDiff);
  x->grp.cycles[x->slot] = x->Icycle;
  
  if(x->myBug == 10)
    {
      for(x->q = 0; x->q < mark; x->q++)
	{
	  post("eSize %f, eOff %f, eSI %f, den %f, as %f, gs %f, gn %f, eJ %f, jS %f",atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eSize),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eOff),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eSizeInv),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.denom),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.allStep),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.groupStep),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.groupNum),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eJoin),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.jSize));
	}
      post("");
      post("x->grp.cycles[%d] = %d",x->slot,x->grp.cycles[x->slot]);
      post("");
      post("nGroups[%d] = %d",x->slot,x->grp.nGroups[x->slot]);
      post("");
      for(x->q = 0; x->q < x->grp.nGroups[x->slot]; x->q++)
	{
	  post("start = %d, Gn = %f, Gd = %f, offset = %f, size = %f, sizeInv = %f",x->grp.gStart[groupOffset + x->q],atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.n),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.d),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.offset),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.size),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.sizeInv));
	}
    }
}

void polyMath_tilde_thisSlot(t_polyMath_tilde *x, t_floatarg f)
{
  x->thisSlot = f < 0 ? 0 : f >= SLOTS ? SLOTS - 1 : (int)f;
}

void polyMath_tilde_groupThisSlot(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  int groupOffset = x->thisSlot * GROUPS;
  int slotOffset = x->thisSlot * MAXSEQ;
  int i, mark;
  x->s = 0;
  x->Woff = 0;
  mark = 0;
  x->seq.len[x->thisSlot] = 0;
  x->grp.nGroups[x->thisSlot] = 0;
  for(x->c = 0; x->c < argc; x->c += 2)
    {
      x->Gn = atom_getfloat(argv + x->c);
      x->Gd = atom_getfloat(argv + x->c + 1);
      //post("Group in slot: %f / %f",x->Gn,x->Gd);
      //post("slot = %d",x->thisSlot);
      if(x->Gd > 0.0f && x->Gn > 0.0f)
	{
	  x->r = x->c / 2;
	  SETFLOAT(&x->grp.n[groupOffset + x->r],x->Gn);
	  SETFLOAT(&x->grp.d[groupOffset + x->r],x->Gd);
	  x->halfSize = (1 / x->Gd) * x->sizeFrac;
	  if(x->autoThreshold && x->halfSize < x->sizeThreshold) x->sizeThreshold = x->halfSize; 
	  x->seq.len[x->thisSlot] += (int)x->Gn;
	  if(x->myBug == 11) post("seq.len[x->slot] = %d",x->seq.len[x->slot]);
	  x->grp.gStart[groupOffset + x->r] = mark;
	  x->t = 0;
	  x->GSize = 0;
	  x->WESize = 1 / x->Gd;
	  x->WSInv = 1 / x->WESize;
	  //post("WESize = %f, WSI = %f",x->WESize,x->WSInv);
	  //for(x->s; x->s < x->s + (int)x->Gn; x->s++)
	  for(x->s = 0;x->s < (int)x->Gn; x->s++)
	    {
	      if(x->Gd > 0.0f && x->Gd > 0.0f)
	      {
		  SETFLOAT(&x->seq.eSize[slotOffset + mark + x->s], x->WESize);
		  SETFLOAT(&x->seq.eOff[slotOffset + mark + x->s], x->Woff);
		  SETFLOAT(&x->seq.eSizeInv[slotOffset + mark + x->s], x->WSInv);
		  SETFLOAT(&x->seq.denom[slotOffset + mark + x->s], x->Gd);

		  SETFLOAT(&x->seq.allStep[slotOffset + mark + x->s], (t_float)mark + (t_float)x->s);
		  SETFLOAT(&x->seq.groupStep[slotOffset + mark + x->s], (t_float)x->t);
		  SETFLOAT(&x->seq.groupNum[slotOffset + mark + x->s], (t_float)x->r);
		  SETFLOAT(&x->seq.eJoin[slotOffset + mark + x->s], 1);
		  SETFLOAT(&x->seq.jSize[slotOffset + mark + x->s], x->WESize);
		  x->GSize += x->WESize;		  
		  x->Woff += x->WESize;
		  x->t++;
		}
	      else post("You can't have size <= 0 - Gn = %d, Gd = %d",(int)x->Gn,(int)x->Gd);
	    }
	  mark += x->s;
	  SETFLOAT(&x->grp.offset[groupOffset + x->r],x->Goff);
	  SETFLOAT(&x->grp.size[groupOffset + x->r],x->GSize);
	  SETFLOAT(&x->grp.sizeInv[groupOffset + x->r],1 / x->GSize);
	  
	  x->Goff = x->Woff;
	  x->grp.nGroups[x->thisSlot]++;
	  x->s++;
	  x->grp.isUnFilled[x->thisSlot] = 0;
	}
      else post("values not greater than 0!");
    }

  x->Icycle = (int)x->Goff;
  x->Gcycle = x->Goff;
  x->cycleDiff = x->Gcycle - (t_float)x->Icycle;
  if(x->cycleDiff > x->sizeThreshold && 1 - x->cycleDiff > x->sizeThreshold)
    {
      if(x->myBug == 10 || x->myBug == 14) post("Gcycle = %f, (t_float)Icycle = %f",x->Gcycle,(t_float)x->Icycle);
      if(x->myBug == 14) post("Gcycle - (t_float)Icycle = %f",x->Gcycle - (t_float)x->Icycle);
      x->Icycle++;
      x->Grem = (t_float)x->Icycle - x->Gcycle;
      if(x->myBug == 10 || x->myBug == 14)
	{
	  post("START");
	  post("x->Goff = %f",x->Goff);
	  post("x->Grem = %f",x->Grem);
	  post("");
	}
      x->Gd = 1 / x->Grem;
      x->Gn = 1;
      SETFLOAT(&x->grp.n[groupOffset + x->r + 1],x->Gn);
      SETFLOAT(&x->grp.d[groupOffset + x->r + 1],x->Gd);
      SETFLOAT(&x->grp.size[groupOffset + x->r + 1],x->Grem);
      SETFLOAT(&x->grp.sizeInv[groupOffset + x->r + 1],x->Gd);
      SETFLOAT(&x->grp.offset[groupOffset + x->r + 1],x->Goff);
      x->grp.nGroups[x->thisSlot]++;
      x->grp.gStart[groupOffset + x->r + 1] = mark;
      
      SETFLOAT(&x->seq.eSize[slotOffset + mark], x->Grem);
      SETFLOAT(&x->seq.eOff[slotOffset + mark], x->Goff);
      SETFLOAT(&x->seq.eSizeInv[slotOffset + mark], x->Gd);
      SETFLOAT(&x->seq.denom[slotOffset + mark], x->Gd);

      SETFLOAT(&x->seq.allStep[slotOffset + mark], (t_float)mark);
      SETFLOAT(&x->seq.groupStep[slotOffset + mark], 0);
      SETFLOAT(&x->seq.groupNum[slotOffset + mark], (t_float)x->r + 1);
      SETFLOAT(&x->seq.eJoin[slotOffset + mark], 1);
      SETFLOAT(&x->seq.jSize[slotOffset + mark], x->Grem);
      mark++;
    }
  else if(x->myBug == 14) post("sizeThreshold = %f, difference = either %f or %f",x->sizeThreshold,x->cycleDiff, 1 - x->cycleDiff);
  x->grp.cycles[x->thisSlot] = x->Icycle;
  
  if(x->myBug == 10)
    {
      for(x->q = 0; x->q < mark; x->q++)
	{
	  post("eSize %f, eOff %f, eSI %f, den %f, as %f, gs %f, gn %f, eJ %f, jS %f",atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eSize),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eOff),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eSizeInv),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.denom),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.allStep),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.groupStep),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.groupNum),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.eJoin),atom_getfloatarg(slotOffset + x->q, x->SEQSIZE, x->seq.jSize));
	}
      post("");
      post("x->grp.cycles[%d] = %d",x->thisSlot,x->grp.cycles[x->thisSlot]);
      post("");
      post("nGroups[%d] = %d",x->thisSlot,x->grp.nGroups[x->thisSlot]);
      post("");
      for(x->q = 0; x->q < x->grp.nGroups[x->thisSlot]; x->q++)
	{
	  post("start = %d, Gn = %f, Gd = %f, offset = %f, size = %f, sizeInv = %f",x->grp.gStart[groupOffset + x->q],atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.n),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.d),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.offset),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.size),atom_getfloatarg(groupOffset + x->q, x->GROUPSIZE, x->grp.sizeInv));
	}
    }
}

void polyMath_tilde_setGroups(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
// 31/10 We need an init 1 slot method so that the slot can be wiped before setting groups
{
  x->Goff = 0;
  x->Gcycle = 0;
  x->Icycle = 0;
  x->Eoff = 0;
  if(x->Grem <= 0) x->Grem = 1;
  x->groupOffset = 0;
  //x->seq.len[x->slot] = 0;
  x->c = 0;
  x->Gstart = 0;
  if(argc >= 3)
    {
      x->slot = atom_getfloat(argv);
      x->slot = x->slot < 0 ? 0 : x->slot > SLOTS - 1 ? SLOTS - 1 : x->slot;
      //      x->Eoffset = 0;
      for(x->c = 0; x->c < (argc - 1) / 2; x->c++)
	//      for(x->c = 0; x->c < (argc - 1) / 2; x->c++)
	{
	  x->Gn = atom_getfloat(argv + 1 + (x->c * 2));
	  x->Gd = atom_getfloat(argv + 2 + (x->c * 2));
	  if(x->myBug == 12) post("x->Gn = %f, x->Gd = %f, offset = %d",x->Gn,x->Gd,1 + (x->c * 2));
	  if(x->Gn <= 0 || x->Gd <= 0)
	    {
	      post("Error: numerator and denominator must both be > 0");
	      SETFLOAT(&x->grp.remains[x->c + x->slot * GROUPS], x->Grem);
	      SETFLOAT(&x->grp.n[x->c + x->slot * GROUPS], 1);
	      SETFLOAT(&x->grp.d[x->c + x->slot * GROUPS], 1 / x->Grem);
	    }
	  else
	    {
	      x->grp.isUnFilled[x->slot] = 0;
	      //reverting this to writeGroup
	      /* x->Wsize = x->Gn / x->Gd;
	      if(x->myBug > 0) post("Wsize = %f",x->Wsize);
	      SETFLOAT(&x->grp.n[x->slot * GROUPS + x->c],x->Gn);
	      SETFLOAT(&x->grp.d[x->slot * GROUPS + x->c],x->Gd);
	      SETFLOAT(&x->grp.size[x->slot * GROUPS + x->c],x->Wsize);
	      SETFLOAT(&x->grp.sizeInv[x->slot * GROUPS + x->c],1/x->Wsize);
	      SETFLOAT(&x->grp.offset[x->slot * GROUPS + x->c],x->groupOffset);
	      //end of trying this here
	      x->h = writeGroup(x,x->c); */
	      if(x->Gn > x->Gd)
		{
		  x->WsizeRem = x->Gn / x->Gd;
		  x->IsizeRem = (int)x->WsizeRem;
		  x->Grem = 1 - (x->WsizeRem - (t_float)x->IsizeRem);
		}
	      if(x->h > 0)
		{
		  x->grp.gStart[x->slot * GROUPS + x->c] = x->Gstart;
		  if(x->myBug > 0) post("group write exit code: %d",x->h);
		  x->c++;
		  if(x->myBug == 1) post("x->Gn = %d, x->Gstart = %d",(int)x->Gn,(int)x->Gstart);
		  x->Gstart += x->Gn;
		  x->Gcycle += (x->Gn / x->Gd);
		  x->seq.len[x->slot] += (int)x->Gn;
		}
	      else post("group write unsuccessful");
	    }
	}
      //remember - the entry after the final group (remainder entry if incomplete, or the last complete entry)
      // last must be all zeroes?
      //
      if(x->myBug > 0) post("x->c = %d",x->c);
      x->Icycle = (int)x->Gcycle;
      if(x->Gcycle - (t_float)x->Icycle > 0.001)
	{
	  x->Grem = x->Gcycle - (t_float)x->Icycle;
	  //x->c++;
	  SETFLOAT(&x->grp.n[x->c],1);
	  SETFLOAT(&x->grp.d[x->c],1 / x->Grem);

	  x->Gn = 1;
	  x->Gd = 1 / x->Grem;
	  x->h = writeGroup(x,x->c);
	  x->grp.fillGroup[x->slot] = x->c;
	  x->Gcycle = (t_float)x->Icycle + 1;
	  if(x->myBug == 5) post("Icycle = %d, Gcycle = %f, x->Grem = %f",x->Icycle,x->Gcycle,x->Grem);
	  x->grp.nGroups[x->slot] = (argc - 1) / 2 + 1;
	  x->seq.len[x->slot]++;
	}
      else
	{
	  x->grp.nGroups[x->slot] = (argc - 1) / 2;
	  x->grp.fillGroup[x->slot] = 0;
	}
      /*if(x->Gcycle > (t_float)((int)x->Gcycle))
	{
	  x->grp.cycles[x->slot] = (int)x->Gcycle + 1;
	}
	else*/ x->grp.cycles[x->slot] = (int)x->Gcycle;
      post("Gcycle = %f",x->Gcycle);
      SETFLOAT(&x->grp.remains[x->slot],x->Grem); // if there is a gap at the end, it is this long
      //x->seq.len[x->slot] = x->e + 1; // length of the sequence
    }
  //  getVariables(x); // perhaps we might not do this here!
}

void polyMath_tilde_setP(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  if(argc == 5) // Pslot, step, (p1/2/3/4), pNum, pVal
    {
      x->PSlot = (int)atom_getfloat(argv);
      x->Pac = (int)atom_getfloat(argv+2);
      switch(x->Pac)
	{
	case(1):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  if(x->myBug > 0)
	    {
	      post("x->Location = %d",x->PLStep + x->PSlot * MAXSEQ);
	      post("x->PSlot = %d",x->PSlot);
	    }
	  SETFLOAT(&x->seq.pAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(2):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc2[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc2[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(3):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc3[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc3[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(4):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc4[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc4[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	case(5):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc5[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc5[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	case(6):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc6[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc6[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	case(7):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc7[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc7[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	case(8):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc8[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc8[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	default:
	  break;
	}
    }
  else
    {
      post("pSet takes a 5 element list: [slotNum, stepNum, p1/2/3/4/5/6/7/8, pNum, pVal");
    }
}

void polyMath_tilde_setPOnly(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  if(argc == 4) // Pslot, step, (p1/2/3/4), pNum
    {
      x->PSlot = (int)atom_getfloat(argv);
      x->Pac = (int)atom_getfloat(argv+2);
      switch(x->Pac)
	{
	case(1):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(2):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  if(x->myBug > 0)
	    {
	      post("x->Location = %d",x->PLStep + x->PSlot * MAXSEQ);
	      post("x->PSlot = %d",x->PSlot);
	    }
	  SETFLOAT(&x->seq.pAcc2[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc2[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(3):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc3[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc3[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(4):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc4[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc4[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(5):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc5[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(6):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc6[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(7):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc7[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	case(8):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  SETFLOAT(&x->seq.pAcc8[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  //SETFLOAT(&x->seq.eAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+4));
	  break;
	default:
	  break;
	}
    }
  else
    {
      post("pSetOnly takes a 4 element list: [slotNum, stepNum, p1/2/3/4/5/6/7/8, pNum");
    }
}

void polyMath_tilde_setVOnly(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  if(argc == 4) // Pslot, step, (p1/2/3/4), vNum
    {
      x->PSlot = (int)atom_getfloat(argv);
      x->Pac = (int)atom_getfloat(argv+2);
      switch(x->Pac)
	{
	case(1):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  //SETFLOAT(&x->seq.pAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	case(2):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  if(x->myBug > 0)
	    {
	      post("x->Location = %d",x->PLStep + x->PSlot * MAXSEQ);
	      post("x->PSlot = %d",x->PSlot);
	    }
	  //SETFLOAT(&x->seq.pAcc2[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc2[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	case(3):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  //SETFLOAT(&x->seq.pAcc3[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc3[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	case(4):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  //SETFLOAT(&x->seq.pAcc4[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc4[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	case(5):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  //SETFLOAT(&x->seq.pAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc5[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	case(6):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  //SETFLOAT(&x->seq.pAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc6[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	case(7):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  //SETFLOAT(&x->seq.pAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc7[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	case(8):
	  x->PLStep = (int)atom_getfloat(argv+1);
	  x->PLStep = x->PLStep >= MAXSEQ ? MAXSEQ - 1 : x->PLStep < 0 ? 0 : x->PLStep;
	  //SETFLOAT(&x->seq.pAcc1[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  SETFLOAT(&x->seq.eAcc8[x->PLStep + x->PSlot * MAXSEQ],atom_getfloat(argv+3));
	  break;
	default:
	  break;
	}
    }
  else
    {
      post("vSetOnly takes a 4 element list: [slotNum, stepNum, p1/2/3/4/5/6/7/8, pVal");
    }
}

/* data types:
 * dType --- data
 * 99    --- {slot, variation}
 * 98    --- {slot, length}
 */

void polyMath_tilde_slot(t_polyMath_tilde *x, t_floatarg f)
{
  x->slot = f < 0 ? 0 : f >= SLOTS ? SLOTS - 1 : (int)f;
  //SETFLOAT(&x->seq.wrapCycles1[x->s],x->JoffNext);
  //x->outList[0] = (t_float)x->slot;
  //x->outList[1] = (t_float)x->variation;
  SETFLOAT(&x->outList[0], (t_float)x->slot);
  SETFLOAT(&x->outList[1], (t_float)x->variation);
  outlet_float(x->dType, 99);
  outlet_list(x->dataOut, gensym("list"), 2, x->outList);
}

void polyMath_tilde_slotLen(t_polyMath_tilde *x, t_floatarg f)
{
  x->getSlotLen = f < 0 ? 0 : f >= SLOTS ? SLOTS - 1 : (int)f;
  //int isLength, int getSlot;
  x->isLength = x->seq.len[x->getSlotLen];
  SETFLOAT(&x->outList[0], (t_float)x->getSlotLen);
  SETFLOAT(&x->outList[1], (t_float)x->isLength);
  outlet_float(x->dType, 98);
  outlet_list(x->dataOut, gensym("list"), 2, x->outList);  
}

//useless: lastLen
void polyMath_tilde_jumpNext(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  //int lastOffset, nextOffset;
  x->validJumpState = 0;
  //lastOffset = x->slot * MAXSEQ; nextOffset = x->slot * MAXSEQ; // just initialize in case no args presented
  if(argc == 2)
    {
      x->lastVar = x->varPerf;
      x->lastSlot = x->slot;
      x->nextSlot = (int)atom_getfloat(argv);
      x->nextSlot = x->nextSlot < 0 ? 0 : x->nextSlot >= SLOTS ? SLOTS - 1 : x->nextSlot; 
      x->nextVar = (int)atom_getfloat(argv+1);
      x->nextVar = x->nextVar < 0 ? 0 : x->nextVar >= VARIATIONS ? VARIATIONS : x->nextVar;
      if(x->scrambling == 1) x->JlastOffset = x->lastSlot * MAXSEQ + x->lastVar * x->SEQSIZE;
      else x->JlastOffset = x->lastSlot * MAXSEQ;
      if(x->nextVar > 0)
	{
	  x->JnextOffset = x->nextSlot * MAXSEQ + x->varPerf * x->SEQSIZE;
	  if(x->myBug == 10) post("x->JnextOffset = %d, x->nextSlot = %d, x->varPerf = %d",x->JnextOffset,x->nextSlot,x->varPerf);
	}
      else x->JnextOffset = x->nextSlot * MAXSEQ;
      if(x->nextVar == 0)
	{
	  if(x->seq.len[x->nextSlot] == 0)
	    {
	      post("Invalid jump state - sequence not defined {:-(");
	      x->validJumpState = 0;
	    }
	  else
	    {
	      post("Jumping to slot %d",x->nextSlot);
	      x->validJumpState = 1;
	    }
	}
      else
	{
	  if(x->var.len[x->slot + (x->nextVar - 1) * SLOTS] == 0)
	    {
	      post("Invalid jump state - variation not defined [;-(");
	      x->validJumpState = 0;
	    }
	  else
	    {
	      post("Jumping to variation %d of slot %d",x->nextVar,x->nextSlot);
	      x->validJumpState = 1;
	    }
	}    
    }
  else if(argc == 1)
    {
      x->nextVar = 0;
      x->lastVar = x->varPerf;
      x->lastSlot = x->slot;
      x->nextSlot = (int)atom_getfloat(argv);
      x->nextSlot = x->nextSlot < 0 ? 0 : x->nextSlot >= SLOTS ? SLOTS - 1 : x->nextSlot;
      if(x->myBug == 10) post("nextSlot = %d",x->nextSlot);
      x->nextVar = 0;
      //x->varPerf = 0;
      if(x->scrambling == 1) x->JlastOffset = x->lastSlot * MAXSEQ + x->lastVar * x->SEQSIZE;
      else x->JlastOffset = x->lastSlot * MAXSEQ;
      x->JnextOffset = x->nextSlot * MAXSEQ;
      if(x->seq.len[x->nextSlot] == 0)
	{
	  post("Invalid jump state - sequence not defined {:-(");
	  x->validJumpState = 0;
	}
      else
	{
	  post("Jumping to slot %d",x->nextSlot);
	  x->validJumpState = 1;
	  x->varPerf = 0;
	}
    }
  if(x->validJumpState)
    {
      x->varTest = x->nextVar;
      x->zeroNextSlot = 0;
      x->zeroNextVar = 0;
      x->changeSlot = 0;
      x->changeVar = 0;
      if(x->nextVar == 0)
	{
	  x->jumpSlotAtEnd = 1;	
	  x->jumpVarAtEnd = 0;
	}
      else
	{
	  x->jumpVarAtEnd = 1;
	  x->jumpSlotAtEnd = 0;
	}
    }
}

void polyMath_tilde_jumpTo(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv) // Jan 17th 2018
{
  //int lastOffset, nextOffset;
  //int lastLen, nextLen, iWrap, nextFlag, locateFlag;
  //t_float lastCycle, nextCycle;
  //t_float sizeLast, offLast,
  //  sizeNext, offNext, wrapCycle;
  x->validJumpState = 0;
  //lastOffset = x->slot * MAXSEQ; x->JnextOffset = x->slot * MAXSEQ; // just initialize in case no args presented
  if(argc == 2)
    {
      x->lastVar = x->varPerf;
      x->lastSlot = x->slot;
      x->nextSlot = (int)atom_getfloat(argv);
      x->nextSlot = x->nextSlot < 0 ? 0 : x->nextSlot >= SLOTS ? SLOTS - 1 : x->nextSlot; 
      x->nextVar = (int)atom_getfloat(argv+1);
      x->nextVar = x->nextVar < 0 ? 0 : x->nextVar >= VARIATIONS ? VARIATIONS : x->nextVar;
      if(x->scrambling == 1) x->JlastOffset = x->lastSlot * MAXSEQ + x->lastVar * x->SEQSIZE;
      else x->JlastOffset = x->lastSlot * MAXSEQ;
      if(x->nextVar > 0)
	{//should this be x->varPerf or x->nextVar?
	  x->varPerf = x->nextVar - 1;
	  x->varTest = x->nextVar;
	  x->JnextOffset = x->nextSlot * MAXSEQ + x->varPerf * x->SEQSIZE;
	  //x->JnextOffset = x->nextSlot * MAXSEQ + (x->nextVar - 1) * x->SEQSIZE;
	  if(x->myBug == 10) post("x->JnextOffset = %d, x->nextSlot = %d, x->varPerf = %d",x->JnextOffset,x->nextSlot,x->varPerf);
	}
      else x->JnextOffset = x->nextSlot * MAXSEQ;
      if(x->nextVar == 0)
	{
	  if(x->seq.len[x->nextSlot] == 0)
	    {
	      post("Invalid jump state - sequence not defined {:-(");
	      x->validJumpState = 0;
	    }
	  else
	    {
	      post("Jumping to slot %d",x->nextSlot);
	      x->varTest = x->nextVar;
	      x->validJumpState = 1;
	    }
	}
      else
	{
	  if(x->var.len[x->slot + (x->nextVar - 1) * SLOTS] == 0)
	    {
	      post("Invalid jump state - variation not defined [;-(");
	      x->validJumpState = 0;
	    }
	  else
	    {
	      x->varPerf = x->nextVar > 0 ? x->nextVar - 1 : 0;
	      post("Jumping to variation %d of slot %d",x->nextVar,x->nextSlot);
	      x->validJumpState = 1;
	      x->varTest = x->nextVar;
	    }
	}    
    }
  else if(argc == 1)
    {
      x->nextVar = 0;
      x->lastVar = x->varPerf;
      x->lastSlot = x->slot;
      x->nextSlot = (int)atom_getfloat(argv);
      x->nextSlot = x->nextSlot < 0 ? 0 : x->nextSlot >= SLOTS ? SLOTS - 1 : x->nextSlot;
      if(x->myBug == 10) post("nextSlot = %d",x->nextSlot);
      x->nextVar = 0;
      //x->varPerf = 0;
      if(x->scrambling == 1) x->JlastOffset = x->lastSlot * MAXSEQ + x->lastVar * x->SEQSIZE;
      else x->JlastOffset = x->lastSlot * MAXSEQ;
      x->JnextOffset = x->nextSlot * MAXSEQ;
      if(x->seq.len[x->nextSlot] == 0)
	{
	  post("Invalid jump state - sequence not defined {:-(");
	  x->validJumpState = 0;
	}
      else
	{
	  post("Jumping to slot %d",x->nextSlot);
	  x->validJumpState = 1;
	  x->varPerf = 0;
	  x->varTest = 0;        ;
	}
    }
  if(argc > 0 && x->validJumpState)
    {
      x->thisInVal = x->InVal;
      if(x->scrambling)
	{
	  x->JlastLen = x->var.len[x->lastSlot + x->lastVar * SLOTS];
	  x->JlastCycle = (t_float)x->vGrp.cycles[x->lastSlot + x->lastVar * SLOTS];
	}
      else
	{
	  x->JlastLen = x->seq.len[x->lastSlot];
	  x->JlastCycle = (t_float)x->grp.cycles[x->lastSlot];
	}
      if(x->nextVar == 0)
	{
	  x->JnextLen = x->seq.len[x->nextSlot];
	  x->JnextCycle = (t_float)x->grp.cycles[x->nextSlot];
	}
      else
	{
	  x->JnextLen = x->var.len[x->nextSlot + x->varPerf * SLOTS];
	  x->JnextCycle = (t_float)x->vGrp.cycles[x->nextSlot + x->nextVar * SLOTS];
	}
      x->JlocateFlag = 1;
      if(!x->scrambling)
	{
	  if(x->PStep == x->seq.len[x->lastSlot] - 1)
	    {
	      if(x->nextVar > 0)
		{
		  x->zeroNextSlot = 0;
		  x->zeroNextVar = 1;
		  x->JlocateFlag = 0;
		  if(x->myBug == 10) post("!x->scrambling && zeroNextVar");
		}
	      else
		{
		  x->zeroNextSlot = 1;
		  x->zeroNextVar = 0;
		  x->JlocateFlag = 0;
		  if(x->myBug == 10) post("!x->scrambling && zeroNextSlot");
		}
	    }
	}
      else
	{
	  if(x->PStep == x->var.len[x->lastSlot + x->lastVar * SLOTS] - 1)
	    {
	      if(x->nextVar > 0)
		{
		  x->zeroNextSlot = 0;
		  x->zeroNextVar = 1;
		  x->JlocateFlag = 0;
		  if(x->myBug == 10) post("x->scrambling && zeroNextVar");
		}
	      else
		{
		  x->zeroNextSlot = 1;
		  x->zeroNextVar = 0;
		  x->JlocateFlag = 0;
		  if(x->myBug == 10) post("x->scrambling && zeroNextSlot");
		}
	    }
	}
      if(x->JlocateFlag)
	{
	  if(x->myBug == 10) post("Locating transition point");
	  if(x->JlastCycle > x->JnextCycle)
	    {
	      if(x->nextVar > 0)
		{
		  if(x->JnextLen > 0)
		    {
		      x->Woff = 0;
		      for(x->s = 0; x->Woff < x->JlastCycle; x->s++)
			{
			  x->JsizeNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->VARSIZE, x->var.eSize);
			  x->JoffNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->VARSIZE, x->var.varOff);
			  x->JiWrap = (x->s / x->JnextLen) * (int)x->JnextCycle;
			  x->Woff += x->JsizeNext;
			  if(x->myBug == 10) post("x->s = %d, x->JsizeNext = %f, x->JoffNext = %f",x->s,x->JsizeNext,x->JoffNext);
			  SETFLOAT(&x->seq.wrapCycles1[x->s],x->JoffNext);
			  SETFLOAT(&x->seq.wrapCycles2[x->s],(t_float)x->JiWrap);
			}
		    }
		  else post("can't do - x->JnextLen = %d",x->JnextLen);
		}
	      else
		{
		  if(x->JnextLen > 0)
		    {
		      x->Woff = 0;
		      for(x->s = 0; x->Woff < x->JlastCycle; x->s++)
			{
			  x->JsizeNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->SEQSIZE, x->seq.eSize);
			  x->JoffNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->SEQSIZE, x->seq.eOff);
			  if(x->myBug == 10) post("x->s = %d, x->JsizeNext = %f, x->JoffNext = %f",x->s,x->JsizeNext,x->JoffNext);
			  x->JiWrap = (x->s / x->JnextLen) * (int)x->JnextCycle;
			  x->Woff += x->JsizeNext;
			  SETFLOAT(&x->seq.wrapCycles1[x->s],x->JoffNext);
			  SETFLOAT(&x->seq.wrapCycles2[x->s],(t_float)x->JiWrap);
			}
		    }
		  else post("can't do - x->JnextLen = %d",x->JnextLen);		  
		}
	    }
  /* If the lastCycle > nextCycle and x->InVal > nextCycle, then the value of x->seq.wrapCycles2[n] will need to be subtracted
   * from x->PGcyc when the transition is made to the next slot / var.
   * A new function should be added to scramble - another value for the input list indicating to scramble e.g 0.5 of the sequence
   * or to scramble the sequence 3 times, etc.
   */
	  else
	    {
	      if(x->nextVar > 0)
		{
		  if(x->JnextLen > 0)
		    {
		      x->Woff = 0;
		      for(x->s = 0; x->Woff < x->JnextCycle; x->s++)
			{
			  x->JsizeNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->VARSIZE, x->var.eSize);
			  x->JoffNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->VARSIZE, x->var.varOff);
			  if(x->myBug == 10) post("x->s = %d, x->JsizeNext = %f, x->JoffNext = %f",x->s,x->JsizeNext,x->JoffNext);
			  x->JiWrap = (x->s / x->JnextLen) * (int)x->JnextCycle;
			  x->Woff += x->JsizeNext;
			  SETFLOAT(&x->seq.wrapCycles1[x->s],x->JoffNext);
			  SETFLOAT(&x->seq.wrapCycles2[x->s],(t_float)x->JiWrap);
			}
		    }
		  else post("Can't do! x->JnextLen = %d",x->JnextLen);
		}
	      else
		{
		  if(x->JnextLen > 0)
		    {
		      x->Woff = 0;
		      for(x->s = 0; x->Woff < x->JnextCycle; x->s++)
			{
			  x->JsizeNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->SEQSIZE, x->seq.eSize);
			  x->JoffNext = atom_getfloatarg(x->JnextOffset + (x->s % x->JnextLen), x->SEQSIZE, x->seq.eOff);
			  if(x->myBug == 10) post("x->s = %d, x->JsizeNext = %f, x->JoffNext = %f",x->s,x->JsizeNext,x->JoffNext);
			  x->JiWrap = (x->s / x->JnextLen) * (int)x->JnextCycle;
			  x->Woff += x->JsizeNext;
			  SETFLOAT(&x->seq.wrapCycles1[x->s],x->JoffNext);
			  SETFLOAT(&x->seq.wrapCycles2[x->s],(t_float)x->JiWrap);
			}
		    }
		  else post("can't do - x->JnextLen = %d",x->JnextLen);
		}
	    }
	  if(x->JnextLen > 0)
	    {
	      x->wrapLen = x->s + 1;
	      x->JnextFlag = 0;
	      for(x->s = 0; x->s < x->wrapLen; x->s++)
		{
		  if(!x->JnextFlag)
		    {
		      x->JoffNext = atom_getfloatarg(x->s,MAXSEQ,x->seq.wrapCycles1);
		      x->JwrapCycle = atom_getfloatarg(x->s,MAXSEQ,x->seq.wrapCycles2);
		      if(x->myBug == 10) post("x->JoffNext + x->JwrapCycle = %f, thisInVal + PGcyc = %f",x->JoffNext + x->JwrapCycle, x->thisInVal + x->PGcyc);
		      if(x->JoffNext + x->JwrapCycle > x->thisInVal + x->PGcyc)
			{
			  x->nextShotVal = x->JoffNext + x->JwrapCycle;
			  x->wrapSubVal = x->JwrapCycle;
			  x->NStep = x->s;
			  if(x->nextVar > 0)
			    {
			      x->changeSlot = 0;
			      x->changeVar = 1;
			    }
			  else
			    {
			      x->changeSlot = 1;
			      x->changeVar = 0;
			    }
			  x->JnextFlag = 1;
			}
		    }
		}
	    }
	  else post("Terminating procedure when length = %d",x->JnextLen);
	}
    }
}

/* ALL SCRAMBLE CODE GOES HERE */

void polyMath_tilde_variation(t_polyMath_tilde *x, t_floatarg f) // OBSOLETE when starting from !0 - use _jumpTo instead
{
  //THIS WILL NEED TO BE REVISITED after the completion of the re-factoring of polyMath_tilde_nextSlot in the perf function!
  x->varTest = f < 0 ? 0 : f > 8 ? 8 : (int)f;
  if(x->varTest > 0)
    {
      x->lastVar = x->varPerf;
      x->varPerf = x->varTest - 1;
      if(x->var.len[x->slot + x->varPerf * SLOTS] > 0) x->scrambling = 1;
    }
}

int copySeq(t_polyMath_tilde *x, int slot, int varOffset)
{
  x->copyWell = 1;
  for(x->o = 0; x->o < x->seq.len[slot]; x->o++)
    {
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.allStep);
      SETFLOAT(&x->var.varStep[varOffset + x->o],x->copyVal);
      SETFLOAT(&x->var.allStep[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.groupStep);
      SETFLOAT(&x->var.groupStep[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.groupNum);
      SETFLOAT(&x->var.groupNum[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eSize);
      if(x->myBug == 8) post("x->copyVal = %f",x->copyVal);
      SETFLOAT(&x->var.eSize[varOffset + x->o],x->copyVal);
      if(x->myBug == 15)
	{
	  post("CHECK COPY VALUES:");
	  post("varOffset = %d",varOffset);
	  post("eSize = %f",x->copyVal);
	  //  for(x->o = 0; x->o < x->seq.len[slot]; x->o++)
	  //x->offsetVar = x->scramSlot * MAXSEQ + x->thisVar * x->SEQSIZE;
	  post("scramSlot = %d, thisVar = %d",(varOffset - (x->thisVar * x->SEQSIZE)) / MAXSEQ, (varOffset - (x->scramSlot * MAXSEQ)) / x->SEQSIZE);
	  post("x->o = %d",x->o);
	}
      if(atom_getfloatarg(varOffset + x->o, x->VARSIZE, x->var.eSize) == 0)
	{
	  post("Event size must be greater than 0 :-( copyVal = %f",x->copyVal);
	  x->copyWell = 0;
	}
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eOff);
      SETFLOAT(&x->var.eOff[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eJoin);
      SETFLOAT(&x->var.eJoin[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.jSize);
      SETFLOAT(&x->var.jSize[varOffset + x->o],x->copyVal);
      if(atom_getfloatarg(varOffset + x->o, x->VARSIZE, x->var.jSize) == 0)
	{
	  post("Join size must be greater than 0 :-(");
	  x->copyWell = 0;
	}
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eSizeInv);
      SETFLOAT(&x->var.eSizeInv[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc1);
      SETFLOAT(&x->var.eAcc1[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc1);
      SETFLOAT(&x->var.pAcc1[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc2);
      SETFLOAT(&x->var.eAcc2[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc2);
      SETFLOAT(&x->var.pAcc2[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc3);
      SETFLOAT(&x->var.eAcc3[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc3);
      SETFLOAT(&x->var.pAcc3[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc4);
      SETFLOAT(&x->var.eAcc4[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc4);
      SETFLOAT(&x->var.pAcc4[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc5);
      SETFLOAT(&x->var.eAcc5[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc5);
      SETFLOAT(&x->var.pAcc5[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc6);
      SETFLOAT(&x->var.eAcc6[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc6);
      SETFLOAT(&x->var.pAcc6[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc7);
      SETFLOAT(&x->var.eAcc7[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc7);
      SETFLOAT(&x->var.pAcc7[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc8);
      SETFLOAT(&x->var.eAcc8[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc8);
      SETFLOAT(&x->var.pAcc8[varOffset + x->o],x->copyVal);
      x->copyVal = atom_getfloatarg(slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.denom);
      SETFLOAT(&x->var.denom[varOffset + x->o],x->copyVal);
      //if(x->myBug == 9) post("varOffset + x->o = %d",varOffset + x->o);//post("d = %f", x->copyVal);
    }
  return(x->copyWell);
}

int scrambleSeq(t_polyMath_tilde *x, int varOffset, int slot)
{
  x->scramWell = 1;
  for(x->q = 0; x->q < x->doSwaps; x->q++) 
    {
      //int scrambleSeq(t_polyMath_tilde *x, int slot, int varOffset, int var, int len, t_float prob)
      //here is where we rewrite VARIATION SEQUENCES
      x->swapVal1 = x->vGrp.swapsRef[x->q];
      x->swapVal2 = x->vGrp.swapsRef[x->q + MAXSEQ];
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.allStep);
      SETFLOAT(&x->var.varStep[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.groupStep);
      SETFLOAT(&x->var.groupStep[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.groupNum);
      SETFLOAT(&x->var.groupNum[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eSize);
      if(x->myBug == 15)
	{
	  post("CHECK VALUES:");
	  post("varOffset = %d",varOffset);
	  post("eSize = %f",x->swapVal);
	  //x->offsetVar = x->scramSlot * MAXSEQ + x->thisVar * x->SEQSIZE;
	  post("scramSlot = %d, thisVar = %d",(varOffset - (x->thisVar * x->SEQSIZE)) / MAXSEQ, (varOffset - (x->scramSlot * MAXSEQ)) / x->SEQSIZE);
	  post("x->q = %d",x->q);
	}
      SETFLOAT(&x->var.eSize[varOffset + x->swapVal2],x->swapVal);
      if(atom_getfloatarg(varOffset + x->swapVal1, x->VARSIZE, x->var.eSize) == 0)
	{
	  post("Event size must be greater than 0 '-(");
	  x->scramWell = 0;
	}
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eOff);
      SETFLOAT(&x->var.eOff[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eJoin);
      SETFLOAT(&x->var.eJoin[varOffset + x->swapVal2],x->swapVal);
      if(x->swapVal > 1)
	{
          for(x->r = 0; x->r < (int)x->swapVal; x->r++)
	    {
	      SETFLOAT(&x->var.eJoin[varOffset + x->swapVal1 + x->r],1);
	      //excludes?
	      if(x->r > 0) SETFLOAT(&x->var.eJoin[varOffset + x->swapVal2 + x->r],1);
	    }
	}
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.jSize);
      SETFLOAT(&x->var.jSize[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eSizeInv);
      SETFLOAT(&x->var.eSizeInv[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc1);
      SETFLOAT(&x->var.eAcc1[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc1);
      SETFLOAT(&x->var.pAcc1[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc2);
      SETFLOAT(&x->var.eAcc2[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc2);
      SETFLOAT(&x->var.pAcc2[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc3);
      SETFLOAT(&x->var.eAcc3[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc3);
      SETFLOAT(&x->var.pAcc3[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc4);
      SETFLOAT(&x->var.eAcc4[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc4);
      SETFLOAT(&x->var.pAcc4[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc5);
      SETFLOAT(&x->var.eAcc5[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc5);
      SETFLOAT(&x->var.pAcc5[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc6);
      SETFLOAT(&x->var.eAcc6[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc6);
      SETFLOAT(&x->var.pAcc6[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc7);
      SETFLOAT(&x->var.eAcc7[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc7);
      SETFLOAT(&x->var.pAcc7[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eAcc8);
      SETFLOAT(&x->var.eAcc8[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.pAcc8);
      SETFLOAT(&x->var.pAcc8[varOffset + x->swapVal2],x->swapVal);
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.denom);
      SETFLOAT(&x->var.denom[varOffset + x->swapVal2],x->swapVal);
      if(x->swapVal == 0) x->scramWell = 0;
      if(x->myBug == 8) post("denom = %d",(int)x->swapVal);
    }
  return(x->scramWell);
}

int regroup(t_polyMath_tilde *x, int grpOffset, int varOffset, int slot, int var, int len)
{
  x->groupWell = 1;
  x->VGnm = 0;
  x->VGCount = 0;
  x->VGSize = 0;
  for(x->p = 0; x->p < len; x->p++)
    {
      //here is where we rewrite GROUPS
      //int VGnm, VGCount, VEJoin;
      //t_float Vd, VESize, VGSize, VEOff, VGOff, VJSize, VJoin, VVStep, VVLast, VONext;
      x->VESize = atom_getfloatarg(varOffset + x->p, x->VARSIZE, x->var.eSize);
      //x->VJSize = atom_getfloatarg(varOffset + x->p, x->VARSIZE, x->var.jSize);
      x->VVStep = atom_getfloatarg(varOffset + x->p, x->VARSIZE, x->var.varStep);
      x->Vd = atom_getfloatarg(varOffset + x->p, x->VARSIZE, x->var.denom);
      x->VEOff = atom_getfloatarg(varOffset + x->p, x->VARSIZE, x->var.eOff);      
      x->vGrp.cycles[slot + var * SLOTS] = x->grp.cycles[slot];
      if(x->p == 0)
	{
	  x->VGSize += x->VESize;
	  x->VGOff = x->VEOff;
	  SETFLOAT(&x->vGrp.offset[grpOffset],x->VGOff);
	  SETFLOAT(&x->vGrp.size[grpOffset],x->VGSize);
	  if(x->VGSize <= 0) x->scramWell = 0;
	  else
	    {
	      x->VGSizeInv = 1 / x->VGSize;
	      SETFLOAT(&x->vGrp.sizeInv[grpOffset],x->VGSizeInv);
	    }
	  SETFLOAT(&x->vGrp.n[grpOffset],(t_float)x->VGCount + 1);
	  SETFLOAT(&x->vGrp.d[grpOffset],x->Vd);
	  x->vGrp.gStart[grpOffset] = (int)x->VVStep;
	  x->vGrp.nGroups[slot + var * SLOTS] = 1; //x->VGnm + 1
	  x->VONext = x->VEOff + x->VESize;
	  SETFLOAT(&x->var.groupStep[varOffset],0);
	  SETFLOAT(&x->var.groupNum[varOffset],0);
	  x->VVLast = x->VVStep;
	  x->VLastD = x->Vd;
	}
      else
	{
	  if(x->VEOff != x->VONext)
	    {
	      x->VGnm++;
	      x->VGCount = 0;
	      x->VGSize = x->VESize;
	      if(x->VGSize == 0) x->scramWell = 0;
	      else x->VGSizeInv = 1 / x->VGSize;
	      x->VGOff = x->VEOff;
	      SETFLOAT(&x->vGrp.n[grpOffset + x->VGnm],(t_float)x->VGCount + 1);
	      SETFLOAT(&x->vGrp.d[grpOffset + x->VGnm],x->Vd);
	      x->vGrp.gStart[grpOffset + x->VGnm] = (int)x->VVStep;
	      x->vGrp.nGroups[slot + var * SLOTS] = x->VGnm + 1;
	    }
	  else if(x->VVStep != x->VVLast + 1)
	    {
	      x->VGnm++;
	      x->VGCount = 0;
	      x->VGSize = x->VESize;
	      if(x->VGSize == 0) x->scramWell = 0;
	      else x->VGSizeInv = 1 / x->VGSize;
	      x->VGOff = x->VEOff;
	      SETFLOAT(&x->vGrp.n[grpOffset + x->VGnm],(t_float)x->VGCount + 1);
	      SETFLOAT(&x->vGrp.d[grpOffset + x->VGnm],x->Vd);
	      x->vGrp.gStart[grpOffset + x->VGnm] = (int)x->VVStep;
	      x->vGrp.nGroups[slot + var * SLOTS] = x->VGnm + 1;
	    }
	  else if(x->VLastD != x->Vd)
	    {
	      x->VGnm++;
	      x->VGCount = 0;
	      x->VGSize = x->VESize;
	      if(x->VGSize == 0) x->scramWell = 0;
	      else x->VGSizeInv = 1 / x->VGSize;
	      x->VGOff = x->VEOff;
	      SETFLOAT(&x->vGrp.n[grpOffset + x->VGnm],(t_float)x->VGCount + 1);
	      SETFLOAT(&x->vGrp.d[grpOffset + x->VGnm],x->Vd);
	      x->vGrp.gStart[grpOffset + x->VGnm] = (int)x->VVStep;
	      x->vGrp.nGroups[slot + var * SLOTS] = x->VGnm + 1;
	    }
	  else 
	    {
	      x->VGCount++;
	      x->VGSize += x->VESize;
	      if(x->VGSize == 0) x->scramWell = 0;
	      else x->VGSizeInv = 1 / x->VGSize;
	    }
	  SETFLOAT(&x->vGrp.offset[grpOffset + x->VGnm],x->VGOff);
	  SETFLOAT(&x->vGrp.size[grpOffset + x->VGnm],x->VGSize);
	  SETFLOAT(&x->vGrp.sizeInv[grpOffset + x->VGnm],x->VGSizeInv);
	  SETFLOAT(&x->vGrp.n[grpOffset + x->VGnm],(t_float)x->VGCount + 1);
	  SETFLOAT(&x->vGrp.d[grpOffset + x->VGnm],x->Vd);
	  x->VONext = x->VEOff + x->VESize;
	  SETFLOAT(&x->var.groupStep[varOffset + x->p],(t_float)x->VGCount);
	  SETFLOAT(&x->var.groupNum[varOffset + x->p],x->VGnm);
	  x->VVLast = x->VVStep;
	  x->VLastD = x->Vd;
	}
    }
  return(x->groupWell);
}

int varOffsets(t_polyMath_tilde *x, int varOffset, int len)
{
  x->varWrite = 1;
  x->varOff = 0;
  x->swapVal = 0;
  for(x->p = 0; x->p < len; x->p++)
    {
      SETFLOAT(&x->var.varOff[varOffset + x->p], x->varOff);
      x->swapVal1 = (int)atom_getfloatarg(varOffset + x->p, x->VARSIZE, x->var.groupStep);
      if(x->swapVal1 == 0) x->swapVal = x->varOff;
      SETFLOAT(&x->var.grpOff[varOffset + x->p], x->swapVal);
      x->varOff += atom_getfloatarg(varOffset + x->p, x->VARSIZE, x->var.eSize);
      if(x->myBug == 8) post("varOff = %f, grpOff = %f");
    }
  return(x->varWrite);
}

//int scramLen, scramMeth // scramMeth: 0 = no repeats, 1 = allow repeats
//reverting to clockSeg~ working scramble routine as the new one isn't working
int scrambleSwaps(t_polyMath_tilde *x, int slot, int varOffset, int var, int len, t_float prob, int grpOffset)
{
  x->swapWell = 1;
  x->seqLen = x->var.len[slot + var * SLOTS];
  x->fSeqLen = (t_float)x->seqLen;
  x->fHalfSeq = x->fSeqLen * 0.5;
  x->halfSeq = x->seqLen / 2;
  //if(x->scramMeth == 0) x->swapsNum = 
  for(x->q = 0; x->q < x->seqLen; x->q++)
    {
      x->vGrp.swapped[x->q] = 0;
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->q, x->SEQSIZE, x->seq.eJoin);
      if(x->swapVal > 1)
	{
	  x->vGrp.swaps[x->q] = x->q * -1;
	  while(x->swapVal > 1)
	    {
	      x->q++;
	      x->vGrp.swaps[x->q] = MAXSEQ * -1;
	      x->swapVal--;
	    }
	}
      else
	{
	  x->vGrp.swaps[x->q] = x->q;
	}
    }
  if(x->scramMeth == 0)
    {
      x->randNum1 = drand48();
      x->fSwapsNum = x->fSeqLen * prob;
      x->swapsNum = rounder(x,x->fSwapsNum,x->seqLen - 1);
      x->doSwaps = x->swapsNum;
      //x->doSwaps = x->swapsNum = (int)x->fSwapsNum;
      //if(x->fSwapsNum > (t_float)x->swapsNum && x->fSeqLen > x->fSwapsNum) x->extraSwapFlag = 1;
      //else x->extraSwapFlag = 0;
      x->r = 0;
      if(x->scramMeth == 0) // 3rd January 2018 - have yet to take account of joins in this!
	{
	  //while(x->swapsNum)
	  for(x->u = 0; x->u < x->swapsNum; x->u)
	    {
	      x->randNum2 = drand48();
	      x->randNum3 = drand48();
	      x->swapVal1 = (int)((t_float)x->randNum2 * x->fSeqLen);
	      x->swapVal2 = (int)((t_float)x->randNum3 * x->fSeqLen);
	      // May 30th 2021 - what is the > (MAXSEQ * -1) for? Can we get rid of it? does x->scramMeth refer to repeats?
	      if(x->swapVal1 != x->swapVal2 && x->vGrp.swapped[x->swapVal1] == 0 && x->vGrp.swapped[x->swapVal2] == 0 && x->vGrp.swaps[x->swapVal1])// > (MAXSEQ * -1) && x->vGrp.swaps[x->swapVal2] > (MAXSEQ * -1))
		{
		  // how do we deal with joins here?
		  x->vGrp.swapped[x->swapVal1] = 1;
		  x->vGrp.swapped[x->swapVal2] = 1;
		  x->vGrp.swaps[x->swapVal2] = x->swapVal1;
		  x->vGrp.swaps[x->swapVal1] = x->swapVal2;
		  x->vGrp.swapsRef[x->r] = x->swapVal1;
		  x->vGrp.swapsRef[x->r + MAXSEQ] = x->swapVal2;
		  x->r++;
		  //x->swapsNum--;
		  x->u++;
		  if(x->swapVal1 < 0)
		    {
		      x->s = x->swapVal1 * -1 - 1;
		      x->t = 1;
		      while(x->s)
			{
			  x->vGrp.swaps[x->swapVal2 + x->t] = MAXSEQ * -1;
			  x->vGrp.swaps[x->swapVal1 + x->t] = x->swapVal1 + x->t;
			  x->s--;
			}
		    }
		  else if(x->swapVal2 < 0)
		    {
		      x->s = x->swapVal2 * -1 - 1;
		      x->t = 1;
		      while(x->s)
			{
			  x->vGrp.swaps[x->swapVal1 + x->t] = MAXSEQ * -1;
			  x->vGrp.swaps[x->swapVal2 + x->t] = x->swapVal2 + x->t;
			  x->s--;
			}
		    }
		}
	    }
	}
    }
  //January 3rd 2018 - have yet to write x->scramMeth == 1
  return(x->swapWell);
}
//here is the newer version:
/*int scrambleSwaps(t_polyMath_tilde *x, int slot, int varOffset, int var, int len, t_float prob, int grpOffset)
{
  x->swapWell = 1;
  x->seqLen = x->var.len[slot + var * SLOTS];
  x->fSeqLen = (t_float)x->seqLen;
  x->fHalfSeq = x->fSeqLen * 0.5;
  x->halfSeq = x->seqLen / 2;
  //if(x->scramMeth == 0) x->swapsNum = 
  for(x->q = 0; x->q < x->seqLen; x->q++)
    {
      x->vGrp.swapped[x->q] = 0;
      x->swapVal = atom_getfloatarg(slot * MAXSEQ + x->q, x->SEQSIZE, x->seq.eJoin);
      if(x->swapVal > 1)
	{
	  x->vGrp.swaps[x->q] = x->q * -1;
	  while(x->swapVal > 1)
	    {
	      x->q++;
	      x->vGrp.swaps[x->q] = MAXSEQ * -1;
	      x->swapVal--;
	    }
	}
      else
	{
	  x->vGrp.swaps[x->q] = x->q;
	}
    }
  x->randNum1 = drand48();
  x->fSwapsNum = x->fSeqLen * prob;
  x->swapsNum = rounder(x,x->fSwapsNum,x->seqLen - 1);
  x->doSwaps = x->swapsNum;
  //x->doSwaps = x->swapsNum = (int)x->fSwapsNum;
  //if(x->fSwapsNum > (t_float)x->swapsNum && x->fSeqLen > x->fSwapsNum) x->extraSwapFlag = 1;
  //else x->extraSwapFlag = 0;
  x->r = 0;
  // 3rd January 2018 - have yet to take account of joins in this!
  //while(x->swapsNum)
  for(x->u = 0; x->u < x->swapsNum; x->u++)
    {
      x->randNum2 = drand48();
      x->randNum3 = drand48();
      x->swapVal1 = (int)((t_float)x->randNum2 * x->fSeqLen);
      x->swapVal2 = (int)((t_float)x->randNum3 * x->fSeqLen);
      if(x->scramMeth == 0)
	{
	  if(x->swapVal1 != x->swapVal2 && x->vGrp.swapped[x->swapVal1] == 0 && x->vGrp.swapped[x->swapVal2] == 0 && x->vGrp.swaps[x->swapVal1] > (MAXSEQ * -1) && x->vGrp.swaps[x->swapVal2] > (MAXSEQ * -1))
	    {
	      // how do we deal with joins here?
	      x->vGrp.swapped[x->swapVal1] = 1;
	      x->vGrp.swapped[x->swapVal2] = 1;
	      x->vGrp.swaps[x->swapVal2] = x->swapVal1;
	      x->vGrp.swaps[x->swapVal1] = x->swapVal2;
	      x->vGrp.swapsRef[x->r] = x->swapVal1;
	      x->vGrp.swapsRef[x->r + MAXSEQ] = x->swapVal2;
	      x->r++;
	      //x->swapsNum--;
	      x->u++;
	      if(x->swapVal1 < 0) // is this the join routine?
		{
		  x->s = x->swapVal1 * -1 - 1;
		  x->t = 1;
		  while(x->s)
		    {
		      x->vGrp.swaps[x->swapVal2 + x->t] = MAXSEQ * -1;
		      x->vGrp.swaps[x->swapVal1 + x->t] = x->swapVal1 + x->t;
		      x->s--;
		    }
		}
	      else if(x->swapVal2 < 0)
		{
		  x->s = x->swapVal2 * -1 - 1;
		  x->t = 1;
		  while(x->s)
		    {
		      x->vGrp.swaps[x->swapVal1 + x->t] = MAXSEQ * -1;
		      x->vGrp.swaps[x->swapVal2 + x->t] = x->swapVal2 + x->t;
		      x->s--;
		    }
		}
	    }
	}
      else if(x->scramMeth == 1)
	{
	  if(x->swapVal1 != x->swapVal2 && x->vGrp.swaps[x->swapVal1] > (MAXSEQ * -1) && x->vGrp.swaps[x->swapVal2] > (MAXSEQ * -1))
	    {
	      // how do we deal with joins here?
	      x->vGrp.swapped[x->swapVal1] = 1;
	      x->vGrp.swapped[x->swapVal2] = 1;
	      x->vGrp.swaps[x->swapVal2] = x->swapVal1;
	      x->vGrp.swaps[x->swapVal1] = x->swapVal2;
	      x->vGrp.swapsRef[x->r] = x->swapVal1;
	      x->vGrp.swapsRef[x->r + MAXSEQ] = x->swapVal2;
	      x->r++;
	      //x->swapsNum--;
	      x->u++;
	      if(x->swapVal1 < 0) // is this the join routine?
		{
		  x->s = x->swapVal1 * -1 - 1;
		  x->t = 1;
		  while(x->s)
		    {
		      x->vGrp.swaps[x->swapVal2 + x->t] = MAXSEQ * -1;
		      x->vGrp.swaps[x->swapVal1 + x->t] = x->swapVal1 + x->t;
		      x->s--;
		    }
		}
	      else if(x->swapVal2 < 0)
		{
		  x->s = x->swapVal2 * -1 - 1;
		  x->t = 1;
		  while(x->s)
		    {
		      x->vGrp.swaps[x->swapVal1 + x->t] = MAXSEQ * -1;
		      x->vGrp.swaps[x->swapVal2 + x->t] = x->swapVal2 + x->t;
		      x->s--;
		    }
		}
	    }
	}
    }
  //January 3rd 2018 - have yet to write x->scramMeth == 1
  return(x->swapWell);
  }*/

void polyMath_tilde_noRepeats(t_polyMath_tilde *x, t_floatarg f)
{
  x->noRepeats = f !=0 ? 1 : 0;
}

void polyMath_tilde_scramble(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  if(argc < 1 && argc > 3)
    {
      post("Incorrect arguments to scramble!");
    }
  else
    {
      if(argc == 3)
	{
	  x->scramSlot = (int)atom_getfloat(argv);
	  x->scramSlot = x->scramSlot < 0 ? 0 : x->scramSlot >= SLOTS ?  SLOTS - 1 : x->scramSlot;
	  x->variation = (int)atom_getfloat(argv+1);
	  x->variation = x->variation < 0 ? 0 : x->variation > VARIATIONS ? VARIATIONS : x->variation;
	  x->seqProb = atom_getfloat(argv+2);
	  x->seqProb = x->seqProb < 0 ? 0 : x->seqProb > 1 ? 1 : x->seqProb;
	  x->seqProb *= 0.5;
	  if(x->seqProb > 0.5)
	    {
	      x->seqProb = 0.5;
	      post("Probability must be within the bounds 0, 1!");
	    }
	  if(x->myBug == 9) post("x->scramSlot = %d, x->variation = %d, x->seqProb = %f",x->scramSlot, x->variation, x->seqProb);
	}
      else if(argc == 2)
	{
	  x->scramSlot = x->slot;
	  x->variation = (int)atom_getfloat(argv);
	  x->variation = x->variation < 0 ? 0 : x->variation > VARIATIONS ? VARIATIONS : x->variation;
	  x->seqProb = atom_getfloat(argv+1);
	  x->seqProb = x->seqProb < 0 ? 0 : x->seqProb > 1 ? 1 : x->seqProb;
	  x->seqProb *= 0.5;
	  /*	  if(x->seqProb > 0.5)
	    {
	      x->seqProb = 0.5;
	      post("Probability must be within the bounds 0, 1!");
	      }*/
	}
      else if(argc == 1)
	{
	  x->scramSlot = x->slot;
	  //x->variation = x->variation;
	  // variation is whatever it was last set to - for an instantaneous scramble
	  x->seqProb = atom_getfloat(argv);
	  x->seqProb = x->seqProb < 0 ? 0 : x->seqProb > 1 ? 1 : x->seqProb;
	  x->seqProb *= 0.5;
	  /*	  if(x->seqProb > 0.5)
	    {
	      x->seqProb = 0.5;
	      post("Probability must be within the bounds 0, 1!");
	      }*/
	}
      //	  post("length = %d",length);
      /* We're going to have to watch it here. The value of x->variation will stay the same if the slot is changed */
      if(x->variation == 0) post("You cannot scramble the original sequence (i.e. variation 0)");
      // joined elements should stay joined, so that the value of seqLen will take this into account
      /* int offsetVar, scramSlot, seqLen, remainSeq, o, p, copyWell, scramWell, joinElement, oldLoc, newLoc, scramMeth
       * t_float swapVal, copyVal, seqProbVal, seqRemainRounding
       * double randNum
       *
       */
      else
	{
	  x->seqLen = x->seq.len[x->scramSlot];
	  if(x->myBug == 9) post("x->seqLen = %d", x->seqLen);
	  x->thisVar = x->variation - 1;
	  x->var.len[x->scramSlot + x->thisVar * SLOTS] = x->seqLen;
	  x->offsetVar = x->scramSlot * MAXSEQ + x->thisVar * x->SEQSIZE;
	  if(x->myBug == 15)
	    {
	      post("SCRAMBLE VALUES:");
	      post("x->scramSlot = %d",x->scramSlot);
	      post("x->thisVar = %d",x->thisVar);
	    }
	  x->grpOffset = x->scramSlot * GROUPS + x->GROUPSIZE * x->thisVar;
	  //x->offsetVar = x->scramSlot * MAXSEQ + SLOTS * MAXSEQ * x->thisVar;
	  if(x->myBug == 9) post("x->offsetVar = %d, x->seqLen = %d",x->offsetVar,x->seqLen);
	  if(copySeq(x,x->scramSlot,x->offsetVar) == 1)
	    {
	      post("Sequence copied successfully!");
	      if(scrambleSwaps(x,x->scramSlot,x->offsetVar,x->thisVar,x->seqLen,x->seqProb,x->grpOffset))
		{
		  post("Swaplists compiled!");
		  if(scrambleSeq(x,x->offsetVar,x->scramSlot))
		    {
		      post("Scrambling successful!");
		      if(regroup(x,x->grpOffset,x->offsetVar,x->scramSlot,x->thisVar,x->seqLen))
			{
			  post("Re-grouping successful!");
			  if(varOffsets(x,x->offsetVar,x->seqLen))
			    {
			      post("Var offsets for instant written successfully");
			      x->var.variations[x->slot + x->thisVar * SLOTS] = 1;
			    }
			  else post("Var offset writing unsuccessful ;-(");
			}
		      else post("Re-grouping unsuccessful ;-(");
		    }
		  else post("Scrambling unsuccessful ;-(");
		}
	      else post("Swaplist compilation unsuccessful ;-(");
	    }
	  else post("copy sequence unsuccessful ;-(");
	}
    }
}

void polyMath_tilde_scramMeth(t_polyMath_tilde *x, t_floatarg f)
{
  x->scramMeth = f != 0 ? 1 : 0;
}

void polyMath_tilde_eMult(t_polyMath_tilde *x, t_floatarg f)
{
  x->eMult = f > 0 ? 1 : 0;
}

void polyMath_tilde_preChange(t_polyMath_tilde *x, t_floatarg f) // note that this will need some alteration / intervention in the timeline when nextSlot is implemented in the perf function! - 12th Jan 2018
{
  x->percentVal = f <= 0.1 ? 0.1 : f > 99.9 ? 99.9 : f;
}

// JOINED ELEMENTS
void polyMath_tilde_makeJoin(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  x->joinSuccess = 0;
  // argv: slot, group, location, length - not yet: (, location2, length2 (, location3, length3 etc))
  //could be dicey...if we make a 3 then a 2, the 2 will be ignored
  // make sure we've got the right kind of list:
  if(argc < 4)
    {
      post("makeJoin: you need at least slot, group, location, length");
      x->joinSuccess = 0;
    }
  else if((argc - 4) % 2 == 1)
    {
      post("makeJoin: each join requires at least location, length");
      x->joinSuccess = 0;
    }
  else
    {
      //// MISSING grp.gStart[x->JSlot * GROUPS + x->JGrp] value!!!!
      /// in the following
      // ...this could be key to making it work!
      /* 00:51am 31/10/2017 */
      
      //x->JSlot is an INT
      x->JSlot = (int)atom_getfloat(argv);
      // here we switch to affermative consequences of if statements:
      if(x->JSlot >= 0 && x->JSlot < SLOTS)
	{
	  //JGrp is an INT
	  x->JGrp = (int)atom_getfloat(argv+1);
	  //back to negative consequences:
	  if(x->JGrp >= x->grp.nGroups[x->JSlot] || x->JGrp < 0) post("makeJoin: group must exist in sequence");
	  else
	    {
	      x->JLoc = (int)atom_getfloat(argv+2);
	      x->JLen = (int)atom_getfloat(argv+3);
	      x->JGnm = (int)atom_getfloatarg(x->JGrp + GROUPS * x->JSlot, GROUPS * SLOTS, x->grp.n);
	      if(x->myBug == 7) post("JSlot = %d, JGrp = %d, JLoc = %d, JLen = %d, JGnm = %d",x->JSlot,x->JGrp,x->JLoc,x->JLen,x->JGnm);
	      if(x->myBug == 6) post("JLoc = %d, JLen = %d, JGnm = %d",x->JLoc,x->JLen,x->JGnm);
	      if(x->JLoc >= x->JGnm)
		{
		  post("makeJoin: location is after the group");
		  x->joinSuccess = 0;
		}
	      else if(x->JLoc + x->JLen > x->JGnm)
		{
		  post("makeJoin: location + length is greater than the numerator");
		  x->joinSuccess = 0;
		}
	      else if(x->JLoc < 0 || x->JLen <= 0 || x->JGnm <= 0)
		{
		  post("makeJoin: x->JLoc must be >= 0: %d",x->JLoc);
		  post("makeJoin: x->JLen must be > 0: %d",x->JLen);
		  post("makeJoin: x->JGnm must be > 0: %d",x->JGnm);
		  x->joinSuccess = 0;
		}
	      else // here we go - write the value
		{ // after, we might want to point out overlapping groups (if any)
		  // and possibly have a "safe mode" where errors are sequentially removed e.g. 3 3 1 4 3 1 1 rewritten 3 1 1 4 1 1 1
		  // actually that should be the default!
		  x->JGst = x->grp.gStart[x->JSlot * GROUPS + x->JGrp];
		  if(x->myBug == 7) post("x->JGst = %d",x->JGst);
		  SETFLOAT(&x->seq.eJoin[x->JLoc + x->JGst + x->JSlot * MAXSEQ],(t_float)x->JLen);
		  for(x->j = x->JLoc + 1;x->j < x->JLen; x->j++)
		    {
		      /*if(x->j == 0)
			{
			  SETFLOAT(&x->seq.eJoin[x->j + x->JGst + x->JSlot * MAXSEQ],(t_float)x->JLen);
			  if(x->myBug == 7) post("Location = %d, x->JLen = %d",x->j + x->JGst,x->JLen);
			}
			else*/
		      SETFLOAT(&x->seq.eJoin[x->j + x->JGst + x->JSlot * MAXSEQ],1);
		    }
		  x-> joinSuccess = 1;
		  // The above method should work if the joins are created in one direction then cleaned in another
		  // ...but the consequences to which direction to make and which direction to clean are different.
		  // If errors are found, best to create backwards and clean forwards.
		  // Of course if every statement is musically and functionally complete and correct, it'll all be fine
		  // ...but human behaviour is rarely complete, nor correct!
		}
	    }
	}
      else
	{
	  post("makeJoin: slot must be between 0 and %d",SLOTS - 1);
	  x->joinSuccess = 0;
	}
    }
  if(x->joinSuccess == 1)
    {
      // clean the group!
      for(x->k = 0; x->k < x->JGn; x->k++)
	{
	  x->JLen = (int)atom_getfloatarg(x->k + x->JSlot * MAXSEQ, x->SEQSIZE, x->seq.eJoin);
	  if(x->JLen > 1)
	    {
	      for(x->j = 1; x->j < x->JLen; x->j++) SETFLOAT(&x->seq.eJoin[x->j + x->k + x->JSlot * MAXSEQ],1);
	    }
	  //when an N > 1 value is encountered sequentially, the subsequent N-1 values must be 1 
	  SETFLOAT(&x->seq.jSize[x->k + x->JSlot * MAXSEQ],1/(t_float)x->JLen); // do we need this?
	}
      x->joinSuccess = 0;
    }
}

//This one doesn't work yet
void polyMath_tilde_joinSeq(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  if(argc > 2)
    {
      x->JSlot = (int)atom_getfloat(argv);
      if(x->JSlot < 0 || x->JSlot > 127) post("ERROR: slot must be an integer from 0 to 127!");
      else
	{
	  x->JGrp = (int)atom_getfloat(argv+1);
	  if(x->grp.nGroups[x->JSlot] <= x->JGrp) post("ERROR: group must already exist in sequence!");
	  else
	    {
	      if(x->myBug == 6) post("Into the main routine");
	      x->GroupStart = x->grp.gStart[x->JGrp];
	      x->JGd = atom_getfloatarg(x->JSlot * GROUPS + x->GroupStart, x->GROUPSIZE, x->grp.d);
	      x->JGn = atom_getfloatarg(x->JSlot * GROUPS + x->GroupStart, x->GROUPSIZE, x->grp.n);
	      x->JGnm = (int)x->JGn;
	      x->JESize = x->JGn / x->JGd;
	      x->JoinTot = 0;
	      for(x->g = 0; x->g < argc - 2; x->g++)
		{
		  x->JoinTot += (int)atom_getfloat(argv + 2 + x->g);
		}
	      if(x->myBug == 6) post("JoinTot = %d, JGn = %d, GroupStart = %d",x->JoinTot,(int)x->JGn,x->GroupStart);
	      x->JGn = atom_getfloatarg(x->JSlot * GROUPS + x->GroupStart, x->GROUPSIZE, x->grp.n);
	      if(x->JoinTot != (int)x->JGn)
		{
		  post("ERROR: Joins total is not equal to numerator!");
		  if(x->myBug == 6)
		    {
		      for(x->k = 0; x->k < argc - 2;x->k++)
			{
			  post("%f",atom_getfloat(argv + 2 + x->k));
			}
		    }
		}
	      else
		{
		  x->JGt = 0;
		  x->g = 0;
		  for(x->k = 0; x->k < argc - 2;x->k++)
		    {
  		      /*while(x->JoinTot > 0)
			{ */
		      if(x->myBug == 6) post("Into the while() statements!");
		      x->JJoin = atom_getfloat(argv + 2 + x->g);
		      if(x->JJoin > 1)
			{
			  x->Gstp++;
			  SETFLOAT(&x->seq.eJoin[x->JSlot * MAXSEQ + x->GroupStart + x->g],x->JJoin);
			  SETFLOAT(&x->seq.jSize[x->JSlot * MAXSEQ + x->GroupStart + x->g],x->JJoin * x->JESize);
			  //SETFLOAT(&x->seq.eSize[x->JSlot * MAXSEQ + x->GroupStart + x->g],x->JJoin * x->JESize);
			  //SETFLOAT(&x->seq.eSizeInv[x->JSlot * MAXSEQ + x->GroupStart + x->g],1/(x->JJoin * x->JESize));
			  // at this point it will be necessary to re-write the rest of the sequence
			  x->JoinTot -= (int)x->JJoin;
			  x->JGt += x->JJoin - 1;
			  if(x->JJoin > 1)
			    x->RW = reWriteSeq(x); // get x->JJoin and x->Gstp in reWriteGroups(x)
			}
		      else
			post("ERROR: Join must be >= 1");
		    }
		  for(x->k = x->JGrp + 1; x->k < x->grp.nGroups[x->JSlot]; x->k++)
		    {
		      // reWrite the array in t_groups for gStart
		      x->JGstt = x->grp.gStart[x->JSlot * GROUPS + x->k] - (int)x->JGt;
		      x->grp.gStart[x->JSlot * GROUPS + x->k] = x->JGstt;
		    }
		}
	    }
	}		    
    }
  else
    post("ERROR: not enough arguments to joinSeq");
}

void polyMath_tilde_initSlot(t_polyMath_tilde *x, t_floatarg f)
{
  x->initSlot = f < 0 ? 0 : f > 127 ? 127 : (int)f;  
  x->grp.nGroups[x->initSlot] = 0;
  x->grp.cycles[x->initSlot] = 1;
  for(x->l = 0; x->l < GROUPS; x->l++)
    {
      x->grp.gStart[x->initSlot * GROUPS + x->l] = 0;
      SETFLOAT(&x->grp.n[x->initSlot * GROUPS + x->l],1); // set all to 1/1 to avoid divide-by-zero errors
      SETFLOAT(&x->grp.d[x->initSlot * GROUPS + x->l],1);
      SETFLOAT(&x->grp.offset[x->initSlot * GROUPS + x->l],0);
      SETFLOAT(&x->grp.size[x->initSlot * GROUPS + x->l],1);
      SETFLOAT(&x->grp.sizeInv[x->initSlot * GROUPS + x->l],1);
      //SETFLOAT(&x->grp.rLength[x->initSlot * GROUPS + x->l],1);
      SETFLOAT(&x->grp.remains[x->initSlot * GROUPS + x->l],0);
    }
  x->seq.len[x->initSlot] = 0;
  for(x->l = 0; x->l < MAXSEQ; x->l++)
    {
      SETFLOAT(&x->seq.allStep[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.groupStep[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.groupNum[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eSize[x->initSlot * MAXSEQ + x->l],1);
      SETFLOAT(&x->seq.eSizeInv[x->initSlot * MAXSEQ + x->l],1);
      SETFLOAT(&x->seq.jSize[x->initSlot * MAXSEQ + x->l],1);
      SETFLOAT(&x->seq.eAcc1[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc1[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eAcc2[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc2[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eAcc3[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc3[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eAcc4[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc4[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eAcc5[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc5[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eAcc6[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc6[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eAcc7[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc7[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.eAcc8[x->initSlot * MAXSEQ + x->l],0);
      SETFLOAT(&x->seq.pAcc8[x->initSlot * MAXSEQ + x->l],0);
    }
}

/* perform function should be able to connect joins together
 * two signal outlets: 1) group phase (0-1)
                       2) event phase - may be also joined events phase so 3 events joined together results in a 3-event length ramp
 */

void polyMath_tilde_altOut(t_polyMath_tilde *x, t_floatarg f)
{
  x->altOut = f != 0 ? 1 : 0;
  if(f > 0 && f < 9)
    {
      x->altLen = (int)f;
    }
}

//This is very old! it needs to be revisited last (after unused variables are weeded out of the main struct
void polyMath_tilde_init(t_polyMath_tilde *x, t_symbol *s)
{
  // counters init:
  x->a = x->b = x->c = x->d = x->e = x->g = x->h = 0;//x->i = x->j = x->k = x->l = x->m = x->n = x->o = x->p = x->q = x->v = x->w = x->x = x->y = x->z = 0;
  x->Gstep = 0;
  x->Gn = 16;
  x->Gd = 16;
  x->ESize = 0.0625;
  x->Eoff = 0;
  x->ESInv = 16;
  x->GSize = 1;
  x->GSInv = 1;
  x->Goff = 0;
  x->cycles = 1;
  x->InVal = x->PreVal = 0;
  x->PGcyc = 0;
  // PGcyc, cycles, PEOff, PESize, PJoined, jFlag, jFirst, PJoin, JESize, JPESI = 1 / JESize
  x->PEOff = 0;
  x->JESize = x->PESize = 0.0625;
  x->JPESI = 16;
  x->PJoined = x->jFlag = x->jFirst = 0;
  
  x->slot = x->group = x->Gstp = 0;
  x->Gnm = 0;
  x->Joined = 0;
  x->join = 1;

  x->altOut = 0;
  x->altNum = 0;
  
  for(x->a = 0; x->a < SLOTS; x->a++)
    {
      x->grp.nGroups[x->a] = 1;
      x->grp.cycles[x->a] = 1;
      for(x->b = 0; x->b < MAXSEQ; x->b++)
	if(x->b < GROUPS)
	  {
	    if(x->b == 0)
	      {
		SETFLOAT(&x->grp.n[x->a * x->b],16);
		SETFLOAT(&x->grp.d[x->a * x->b],16);
		SETFLOAT(&x->grp.offset[x->a * x->b],0);
		SETFLOAT(&x->grp.size[x->a * x->b],1);
		SETFLOAT(&x->grp.sizeInv[x->a * x->b],1);
		//SETFLOAT(&x->grp.cycles[x->a * x->b],1);
		//SETFLOAT(&x->grp.rLength[x->a * x->b],1);
		//SETFLOAT(&x->grp.remains[x->a * x->b],0);
	      }
	    else
	      {
		SETFLOAT(&x->grp.n[x->a * x->b],0);
		SETFLOAT(&x->grp.d[x->a * x->b],0);
		SETFLOAT(&x->grp.offset[x->a * x->b],0);
		SETFLOAT(&x->grp.size[x->a * x->b],0);
		SETFLOAT(&x->grp.sizeInv[x->a * x->b],0);
		//SETFLOAT(&x->grp.cycles[x->a * x->b],0);
		//SETFLOAT(&x->grp.rLength[x->a * x->b],0);
		//SETFLOAT(&x->grp.remains[x->a * x->b],0);
	      }
	  }
      if(x->b < 16)
	{
	  SETFLOAT(&x->seq.allStep[x->a * x->b],(t_float)x->b);
	  SETFLOAT(&x->seq.groupStep[x->a * x->b],(t_float)x->b);
	  SETFLOAT(&x->seq.allStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.eSize[x->a * x->b],0.0625);
	  SETFLOAT(&x->seq.eOff[x->a * x->b],0.0625 * (t_float)x->b);
	  SETFLOAT(&x->seq.eJoin[x->a * x->b],1);
	  SETFLOAT(&x->seq.eAcc1[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc2[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc3[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc4[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc5[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc6[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc7[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc8[x->a * x->b],0);
	  SETFLOAT(&x->seq.eSizeInv[x->a * x->b],16);
	}
      else
	{
	  SETFLOAT(&x->seq.allStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.groupStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.allStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.eSize[x->a * x->b],0);
	  SETFLOAT(&x->seq.eOff[x->a * x->b],0);
	  SETFLOAT(&x->seq.eJoin[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc1[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc2[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc3[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc4[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc5[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc6[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc7[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc8[x->a * x->b],0);
	  SETFLOAT(&x->seq.eSizeInv[x->a * x->b],0);
	}
    }
  x->Pthis = 0.0625;
  x->PJoin = 1;
  x->PStep = 0;
  getVariables(x);

  x->barNew = 0;
}

void polyMath_tilde_getSeq(t_polyMath_tilde *x, t_symbol *s, int argc, t_atom *argv)
{
  if(argc == 3)
    {
      // slot, variation, dType
      /* dType slots / vars:
       * When interacting with this object in the app, the data will be modelled both in the app and in here
       * Text-based saving of these internal arrays, and restoration from textfiles needs to be implemented,
       * but it will still be necessary to output data from the object when variations are made, since this
       * happens here and not in Pd or in OF.
       * 0 = eOff, 1 = eSize, 2 = gStep, 3 = gNum (get denom data from grp / vGrp data outputs)
       * 4 = eJoin, 5 = jSize (need to define this better; is it 1s for eJoin and 2+ for jSize? 
       * is it 0,1,2 or 1,2,3 for jSize?)
       * 11 = pAcc1, 12 = eAcc 1, 13 = pAcc2...18 = eAcc4 (we need to add the other 4 channels from eventSeq~, 
       * so keeping 19-26 reserved (in fact, let's skip up to 91)
       * 91 = {gType, nGroups, cycles, length, remains}
       * note that len(gth) comes from sequences / variations 
       * 92 = numerals, 93 = denominators, 
       * 94 = gOffsets, 95 = gSizes, 96 = gStarts
       */
      x->getSlot = atom_getfloat(argv);
      if(x->getSlot >= 0 && x->getSlot < SLOTS)
	{
	  x->getVar = atom_getfloat(argv+1);
	  if(x->getVar >= 0 && x->getVar < VARIATIONS + 1)
	    {
	      if(x->getVar == 0)
		{
		  x->seqOff = x->getSlot * MAXSEQ;
		  x->seqGrpOff = x->getSlot * GROUPS;
		  x->lenSeq = x->seq.len[x->getSlot];
		  x->lenGrp = x->grp.nGroups[x->getSlot];
		  x->getVarNum = -1;
		  x->getPar = atom_getfloat(argv+2);
		  if(x->getPar < 91) // 91+ are group settings, not sequence settings!
		    {
		      if(x->lenSeq < 1)
			{
			  post("This slot has not been filled yet!");
			  outlet_float(x->dType,-1);
			}
		      else
			{
			  switch(x->getPar)
			    {
			    case(0):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eOff));
				}
			      outlet_float(x->dType,0);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(1):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eSize));
				}
			      outlet_float(x->dType,1);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(2):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.groupStep));
				}
			      outlet_float(x->dType,2);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(3):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.groupNum));
				}
			      outlet_float(x->dType,3);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(4):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eJoin));
				}
			      outlet_float(x->dType,4);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(5):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.jSize));
				}
			      outlet_float(x->dType,5);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;

			    case(11):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc1));
				}
			      outlet_float(x->dType,11);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(12):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc1));
				}
			      outlet_float(x->dType,12);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(13):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc2));
				}
			      outlet_float(x->dType,13);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(14):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc2));
				}
			      outlet_float(x->dType,14);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(15):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc3));
				}
			      outlet_float(x->dType,15);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(16):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc3));
				}
			      outlet_float(x->dType,16);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(17):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc4));
				}
			      outlet_float(x->dType,17);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(18):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc4));
				}
			      outlet_float(x->dType,18);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(19):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc5));
				}
			      outlet_float(x->dType,19);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(20):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc5));
				}
			      outlet_float(x->dType,20);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(21):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc6));
				}
			      outlet_float(x->dType,21);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(22):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc6));
				}
			      outlet_float(x->dType,22);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(23):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc7));
				}
			      outlet_float(x->dType,23);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(24):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc7));
				}
			      outlet_float(x->dType,24);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(25):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc8));
				}
			      outlet_float(x->dType,25);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(26):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc8));
				}
			      outlet_float(x->dType,26);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    default:
			      post("That sequence output is undefined (yet!");
			      outlet_float(x->dType,-2);
			      break;
			    }
			}
		    }
		  else if(x->getPar > 90)
		    {
		      if(x->lenGrp < 1)
			{
			  post("This slot has not been filled yet!");
			  outlet_float(x->dType,-1);
			}
		      else
			{
			  switch(x->getPar)
			    {
			    case(91):
			      x->getGrpVal = (t_float)x->grp.gType[x->getSlot];
			      SETFLOAT(&x->outList[0],x->getGrpVal);
			      x->getGrpVal = (t_float)x->grp.nGroups[x->getSlot];
			      SETFLOAT(&x->outList[1],x->getGrpVal);
			      x->getGrpVal = (t_float)x->grp.cycles[x->getSlot];
			      SETFLOAT(&x->outList[2],x->getGrpVal);
			      SETFLOAT(&x->outList[3],x->lenSeq);
			      x->getGrpVal = atom_getfloatarg(x->getSlot,SLOTS,x->grp.remains);
			      SETFLOAT(&x->outList[4],x->getGrpVal);
			      SETFLOAT(&x->outList[5],0);
			      outlet_float(x->dType,91);
			      outlet_list(x->dataOut, gensym("list"), 6, x->outList);
			      break;
			    case(92):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqGrpOff + x->v, x->GROUPSIZE, x->grp.n));
				}
			      outlet_float(x->dType,92);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(93):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqGrpOff + x->v, x->GROUPSIZE, x->grp.d));
				}
			      outlet_float(x->dType,93);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(94):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v], atom_getfloatarg(x->seqGrpOff + x->v, x->GROUPSIZE, x->grp.offset));
				}
			      outlet_float(x->dType,94);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(95):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqGrpOff + x->v, x->GROUPSIZE, x->grp.size));
				}
			      outlet_float(x->dType,95);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(96):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],(t_float)x->grp.gStart[x->seqGrpOff + x->v]);
				}
			      outlet_float(x->dType,96);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    default:
			      post("That group output is undefined (yet!)");
			      outlet_float(x->dType,-2);
			      break;
			    }
			}
		    }
		}
	      else if(x->getVar < 9)
		{ //slot + var * SLOTS
		  x->getVarNum = x->getVar - 1;
		  x->seqOff = x->getSlot * MAXSEQ + x->getVarNum * x->SEQSIZE;
		  x->seqGrpOff = x->getSlot * GROUPS + x->getVar * x->GROUPSIZE; 
		  x->lenSeq = x->seq.len[x->getSlot * VARIATIONS + x->getVarNum];
		  x->lenGrp = x->grp.nGroups[x->getSlot + x->getVarNum * SLOTS];
		  x->getPar = atom_getfloat(argv+2);
      //  x->Pacc2 = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.pAcc2);
      //  x->Goff = atom_getfloatarg(x->slot * GROUPS + x->varPerf * x->GROUPSIZE + x->Gnm, x->VGROUPSIZE, x->vGrp.offset);
		  if(x->getPar < 91) // 91+ are group settings, not sequence settings!
		    {
		      if(x->lenSeq < 1)
			{
			  post("This slot has not been filled yet!");
			  outlet_float(x->dType,-1);
			}
		      else
			{
			  switch(x->getPar)
			    {
			    case(0):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.eOff));
				}
			      outlet_float(x->dType,0);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(1):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.eSize));
				}
			      outlet_float(x->dType,1);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(2):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.groupStep));
				}
			      outlet_float(x->dType,2);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(3):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.groupNum));
				}
			      outlet_float(x->dType,3);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(4):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.eJoin));
				}
			      outlet_float(x->dType,4);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(5):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.jSize));
				}
			      outlet_float(x->dType,5);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;

			    case(11):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.pAcc1));
				}
			      outlet_float(x->dType,11);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(12):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.eAcc1));
				}
			      outlet_float(x->dType,12);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(13):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.pAcc2));
				}
			      outlet_float(x->dType,13);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(14):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.eAcc2));
				}
			      outlet_float(x->dType,14);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(15):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.pAcc3));
				}
			      outlet_float(x->dType,15);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(16):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.eAcc3));
				}
			      outlet_float(x->dType,16);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(17):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.pAcc4));
				}
			      outlet_float(x->dType,17);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(18):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->VARSIZE, x->var.eAcc4));
				}
			      outlet_float(x->dType,18);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(19):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc5));
				}
			      outlet_float(x->dType,19);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(20):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc5));
				}
			      outlet_float(x->dType,20);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(21):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc6));
				}
			      outlet_float(x->dType,21);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(22):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc6));
				}
			      outlet_float(x->dType,22);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(23):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc7));
				}
			      outlet_float(x->dType,23);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(24):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc7));
				}
			      outlet_float(x->dType,24);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(25):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.pAcc8));
				}
			      outlet_float(x->dType,25);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    case(26):
			      for(x->v = 0; x->v < x->lenSeq; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqOff + x->v, x->SEQSIZE, x->seq.eAcc8));
				}
			      outlet_float(x->dType,26);
			      outlet_list(x->dataOut, gensym("list"), x->lenSeq, x->outList);
			      break;
			    default:
			      post("That sequence output is undefined (yet!");
			      outlet_float(x->dType,-2);
			      break;
			    }
			}
		    }
		  else if(x->getPar > 90)
		    {
		      if(x->lenGrp < 1)
			{
			  post("This slot has not been filled yet!");
			  outlet_float(x->dType,-1);
			}
		      else
			{
			  switch(x->getPar)
			    {
			    case(91):
			      x->getGrpVal = (t_float)x->grp.gType[x->getSlot];
			      SETFLOAT(&x->outList[0],x->getGrpVal);
			      x->getGrpVal = (t_float)x->grp.nGroups[x->getSlot];
			      SETFLOAT(&x->outList[1],x->getGrpVal);
			      x->getGrpVal = (t_float)x->grp.cycles[x->getSlot];
			      SETFLOAT(&x->outList[2],x->getGrpVal);
			      SETFLOAT(&x->outList[3],x->lenSeq);
			      x->getGrpVal = atom_getfloatarg(x->getSlot,SLOTS*VARIATIONS,x->grp.remains);
			      SETFLOAT(&x->outList[4],x->getGrpVal);
			      SETFLOAT(&x->outList[5],(t_float)x->getVar);
			      outlet_float(x->dType,91);
			      outlet_list(x->dataOut, gensym("list"), 6, x->outList);
			      break;
			    case(92):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqGrpOff + x->v, x->VGROUPSIZE, x->grp.n));
				}
			      outlet_float(x->dType,92);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(93):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqGrpOff + x->v, x->VGROUPSIZE, x->grp.d));
				}
			      outlet_float(x->dType,93);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(94):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqGrpOff + x->v, x->VGROUPSIZE, x->grp.offset));
				}
			      outlet_float(x->dType,94);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(95):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],atom_getfloatarg(x->seqGrpOff + x->v, x->VGROUPSIZE, x->grp.size));
				}
			      outlet_float(x->dType,95);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    case(96):
			      for(x->v = 0; x->v < x->lenGrp; x->v++)
				{
				  SETFLOAT(&x->outList[x->v],(t_float)x->grp.gStart[x->seqGrpOff + x->v]);
				}
			      outlet_float(x->dType,96);
			      outlet_list(x->dataOut, gensym("list"), x->lenGrp, x->outList);
			      break;
			    default:
			      post("That group output is undefined (yet!)");
			      outlet_float(x->dType,-2);
			      break;
			    }
			}
		    }
		}
	    }
	  else
	    {
	      post("Variation must be 0, or 1 to %d",VARIATIONS);
	    }
	}
      else
	{
	  post("Slot must be from 0 to %d",SLOTS - 1);
	}
    }
  else if(argc == 1)
    {
      
    }
}

// DEBUG CODE
void polyMath_tilde_debug(t_polyMath_tilde *x, t_floatarg myBug)
{
  int bug = (int)myBug;
  t_float val = 0; t_float val2 = 0;
  x->myBug = bug;
  if(bug == 1)
    {
      post("Sequence Offsets:");
      for(x->h = 0;x->h < 16;x->h++)
	{
	  val = atom_getfloatarg(x->slot * MAXSEQ + x->h, x->SEQSIZE, x->seq.eOff);
	  post("%f",val);
	}
    }
  else if(bug == 2)
    {
      post("Event Sizes:");
      for(x->h = 0;x->h < 16;x->h++)
	{
	  val = atom_getfloatarg(x->slot * MAXSEQ + x->h, x->SEQSIZE, x->seq.eSize);
	  post("%f",val);
	}
    }
  else if(bug == 3)
    {
      post("Event allStep:");
      for(x->h = 0;x->h < 16;x->h++)
	{
	  val = atom_getfloatarg(x->slot * MAXSEQ + x->h, x->SEQSIZE, x->seq.allStep);
	  post("%f",val);
	}
    }
  else if(bug == 4)
    {
      post("Event groupStep:");
      for(x->h = 0;x->h < 16;x->h++)
	{
	  val = atom_getfloatarg(x->slot * MAXSEQ + x->h, x->SEQSIZE, x->seq.groupStep);
	  post("%f",val);
	}
    }
  else if(bug == 5)
    {
      x->myBug = 1;
      post("GROUPS_____________________________________________________");
      for(x->h = 0;x->h < 6;x->h++)
	{
	  val = atom_getfloatarg(x->slot * GROUPS + x->h, x->GROUPSIZE, x->grp.n);
	  post("Num: %f",val);
	  val = atom_getfloatarg(x->slot * GROUPS + x->h, x->GROUPSIZE, x->grp.d);
	  post("Den: %f",val);	  
	  val = atom_getfloatarg(x->slot * GROUPS + x->h, x->GROUPSIZE, x->grp.offset);
	  post("Offset: %f",val);	  
	  val = atom_getfloatarg(x->slot * GROUPS + x->h, x->GROUPSIZE, x->grp.size);
	  post("Size: %f",val);
	  val = atom_getfloatarg(x->slot * GROUPS + x->h, x->GROUPSIZE, x->grp.sizeInv);
	  post("SizeInv: %f",val);	  
	  /*val = atom_getfloatarg(x->slot * GROUPS + x->h, x->GROUPSIZE, x->grp.rLength);
	  post("rLength: %f",val);
	  val = atom_getfloatarg(x->slot * GROUPS + x->h, x->GROUPSIZE, x->grp.remains);
	  post("Remains: %f",val);*/
	}
    }
  else if(bug == 6) x->myBug = 6;
  else if(bug == 7)
    {
      x->myBug = 7;
      for(x->h = 0; x->h < x->seq.len[x->slot]; x->h++)
	{
	  val = atom_getfloatarg(x->slot * MAXSEQ + x->h, x->SEQSIZE, x->seq.eJoin);
	  val2 = atom_getfloatarg(x->slot * MAXSEQ + x->h, x->SEQSIZE, x->seq.groupStep);
	  post("join at %d, groupStep = %d, length = %d",x->h, (int)val2, (int)val);
	}
      post("seq.len[slot] = %d",x->seq.len[x->slot]);
      post("grp.nGroups[slot] = %d",x->grp.nGroups[x->slot]);
      for(x->h = 0; x->h < x->grp.nGroups[x->slot]; x->h++) post("group %d, start %d",x->h,x->grp.gStart[x->h + x->slot * GROUPS]);
    }
  else if(bug == 8)
    {
      x->myBug = 8;
      for(x->q = 0; x->q < x->var.len[x->scramSlot + x->thisVar * SLOTS]; x->q++)
	{
	  post("eSize = %f, eOff = %f, eJoin = %d, groupStep = %d, denom = %d, varOff = %f",atom_getfloatarg(x->offsetVar + x->q,x->VARSIZE,x->var.eSize),atom_getfloatarg(x->offsetVar + x->q,x->VARSIZE,x->var.eOff),(int)atom_getfloatarg(x->offsetVar + x->q,x->VARSIZE,x->var.eJoin),(int)atom_getfloatarg(x->offsetVar + x->q,x->VARSIZE,x->var.groupStep),(int)atom_getfloatarg(x->offsetVar + x->q,x->VARSIZE,x->var.denom),atom_getfloatarg(x->offsetVar + x->q,x->VARSIZE,x->var.varOff));
	}
      for(x->p = 0; x->p < x->seq.len[x->scramSlot]; x->p++)
	{
	  post("Denominator = %d",(int)atom_getfloatarg(x->scramSlot * MAXSEQ, x->VARSIZE, x->var.denom));
	}
    }
  else if(bug == 9)
    {
      x->myBug = 9;
      post("x->slot = %d, x->varPerf = %d, ENTRIES:",x->slot,x->varPerf);
      for(x->q = 0; x->q < x->var.len[x->slot + x->varPerf * SLOTS]; x->q++)
	{
	  post("eSize = %f, eOff = %f, eJoin = %d, groupNum = %d, groupStep = %d, denom = %d",atom_getfloatarg(x->varPerf * x->SEQSIZE + x->slot * MAXSEQ + x->q,x->VARSIZE,x->var.eSize),atom_getfloatarg(x->varPerf * x->SEQSIZE + x->slot * MAXSEQ + x->q,x->VARSIZE,x->var.eOff),(int)atom_getfloatarg(x->varPerf * x->SEQSIZE + x->slot * MAXSEQ + x->q,x->VARSIZE,x->var.eJoin),(int)atom_getfloatarg(x->varPerf * x->SEQSIZE + x->slot * MAXSEQ + x->q,x->VARSIZE,x->var.groupNum),(int)atom_getfloatarg(x->varPerf * x->SEQSIZE + x->slot * MAXSEQ + x->q,x->VARSIZE,x->var.groupStep),(int)atom_getfloatarg(x->varPerf * x->SEQSIZE + x->slot * MAXSEQ + x->q,x->VARSIZE,x->var.denom));
	}
    }
  else if(bug == 10)
    {
      x->myBug = 10;
      for(x->q = 0; x->q < x->wrapLen; x->q++)
	{
	  if(x->q < 10) post("swap at _%d: %f",x->q, atom_getfloatarg(x->q,MAXSEQ,x->seq.wrapCycles2));
	  else post("swap at %d: %f",x->q, atom_getfloatarg(x->q,MAXSEQ,x->seq.wrapCycles2));
	}
    }
  else if(bug == 11)
    {
      for(x->q = 0; x->q < 30; x->q++)
	{
	  post("x->seq.len[%d] = %d",x->q,x->seq.len[x->q]);
	}
    }
  else if(bug == 12)
    {
      post("nextSlot = %d, length = %d",x->nextSlot,x->seq.len[x-> nextSlot]);
      for(x->q = 0; x->q < x->seq.len[x->nextSlot]; x->q++)
	{
	  post("x->q = %d, eOff = %f",x->q,atom_getfloatarg(x->nextSlot * MAXSEQ + x->q,x->SEQSIZE,x->seq.eOff));
	}
    }
  else if(bug == 13)
    {
      x->myBug = 13;
      post("len = %d", x->var.len[x->scramSlot + x->thisVar * SLOTS]);
      for(x->q = 0; x->q < x->var.len[x->scramSlot + x->thisVar * SLOTS]; x->q++)
	{
	  post("swapsRef = %d %d",x->vGrp.swapsRef[x->q],x->vGrp.swapsRef[x->q + MAXSEQ]);
	}
    }
  else if(bug == 14) x->myBug = 14;
  else if(bug == 15)
    {
      post("Full swap sequence: ------------------------------------");
      x->myBug = 15;
      for(x->o = 0; x->o < x->seq.len[x->scramSlot]; x->o++)
	{
	  //x->offsetVar = x->scramSlot * MAXSEQ + x->thisVar * x->SEQSIZE;
	  x->swapVal1 = x->vGrp.swapsRef[x->o];
	  x->swapVal = atom_getfloatarg(x->scramSlot * MAXSEQ + x->swapVal1, x->SEQSIZE, x->seq.eSize);
	  post("x->swapVal1 = x->vGrp.swapsRef[x->o] = %d",x->swapVal1);
	  post("x->swapVal = %f",x->swapVal);
	  
	}
    }
  else if(bug == 16)
    {
      for(x->o = 0; x->o < x->seq.len[x->slot]; x->o++)
	{
	  if(x->o < 10) post("index:  %d| p1 %d | e1 %d | p2 %d | e2 %d | p3 %d | e3 %d | p4 %d | e4 %d",x->o,(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc1),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc1),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc2),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc2),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc3),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc3),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc4),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc4));
	}
    }
  else if(bug == 17)
    {
      for(x->o = 0; x->o < x->seq.len[x->slot]; x->o++)
	{
	  if(x->o < 10) post("index:  %d| p5 %d | e5 %d | p6 %d | e6 %d | p7 %d | e7 %d | p8 %d | e8 %d",x->o,(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc5),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc5),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc6),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc6),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc7),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc7),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.pAcc8),(int)atom_getfloatarg(x->slot * MAXSEQ + x->o, x->SEQSIZE, x->seq.eAcc8));
	}
    }
}

// JOIN ROUTINES
static void checkJoinsOut(t_polyMath_tilde *x)
{
  if(x->PJoined > 1) // 2nd November - what's the deal with 1 instead of 0? 14th Nov 2019: This is an error, but we can live with it for now...
    {
      x->PJoined--; // 3, 2, 1
      x->jFirst = 0;
    }
  else
    {
      if(x->PJoin == 1)
	{
	  if(x->jFlag > 0)
	    {
	      x->jFlag = 0;
	      x->JESize = x->PESize;
	      x->JPESI = x->PESInv;
	      if(x->JESize > x->sizeThreshold) clock_delay(x->fOut,0L);
	    }
	  else
	    { // normal step
	      x->JESize = x->PESize;
	      x->JPESI = x->PESInv;
	      if(x->JESize > x->sizeThreshold) clock_delay(x->fOut,0L);
	    }
	}
      else
	{
	  x->jFirst = 1;
	  x->JoinVal = (int)x->PJoin;
	  x->PJoined = x->JoinVal;
	  x->JESize = x->PESize * x->PJoin;
	  x->JPESI = x->PESInv / x->PJoin;
	  if(x->JESize > x->sizeThreshold) clock_delay(x->fOut,0L);
	}
    }
}

static void checkJoinsVarOut(t_polyMath_tilde *x)
{
  if(x->PJoined > 1)
    {
      x->PJoined--; // 3, 2, 1
      x->jFirst = 0;
    }
  else
    {
      if(x->PJoin == 1)
	{
	  x->eChanged = 0;
	  if(x->jFlag > 0)
	    {
	      x->jFlag = 0;
	      x->JESize = x->PESize;
	      x->VPESI = 1 / x->PESize;
	      if(x->JESize > x->sizeThreshold) clock_delay(x->fOut,0L);
	    }
	  else
	    { // normal step
	      x->JESize = x->PESize;
	      x->VPESI = 1 /x->PESize;
	      if(x->JESize > x->sizeThreshold) clock_delay(x->fOut,0L);
	    }
	}
      else // Joined event
	{
	  x->eChanged = 0;
	  x->jFirst = 1;
	  x->JoinVal = (int)x->PJoin;
	  x->PJoined = x->JoinVal;
	  x->JESize = x->PESize * x->PJoin;
	  x->VPESI = 1 / x->JESize; // perhaps it's time to use the jSize array in var?
	  x->jFlag = 1;
	  if(x->JESize > x->sizeThreshold) clock_delay(x->fOut,0L);
	}
    }
}

void polyMath_tilde_setBpm(t_polyMath_tilde *x, t_floatarg f)
{
  if(f > 0)
    {
      x->BPM = f;
      x->durBeat = 60000/f;
      x->barBeat = x->durBeat * 4;
    }
  else post("bpm must be a positive number!");
}

// PERFORM ROUTINE
//13th Jan 2018 - a better idea: while(n--) { if(x->changeFlag) {...} else { (the normal slot or var routines) }
t_int *polyMath_tilde_perform(t_int *w)
{
  t_polyMath_tilde *x = (t_polyMath_tilde *)(w[1]);
  t_float *in         = (t_float *)(w[2]);
  t_float *out3       = (t_float *)(w[3]);
  t_float *alt3       = (t_float *)(w[4]);
  t_float *offset     = (t_float *)(w[5]);
  int n               = (int)(w[6]);
  if(x->firstStart == 1)
    {
      getVariables(x);
      clock_delay(x->fOut,0L);
      x->firstStart = 0;
    }
  if(x->seq.len[x->slot] == 0)
    {
      getVariables(x);
      while(n--)
	{
	  *out3++ = *in++;
	  *alt3++ = 0;
	  *offset++ = 0;
	}
    }
  else while(n--)
    {
      x->InVal = *in++;
      if(!x->scrambling)
	{
	  if(x->InVal < x->PreVal)
	    {
	      if(x->zeroNextPhase)
		{
		  //x->eChanged = 0;
		  //x->slot = x->nextSlot;
		  x->PStep = 0;
		  x->PGcyc = 0;
		  x->zeroNextPhase = 0;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ, x->SEQSIZE, x->seq.eOff);
		  //x->pageNum = 0;
		  //x->pageFlag = 1;
		  getVariables(x);
		  x->scrambling = 0;
		  checkJoinsOut(x);
		  //add join code in here!
		}
	      else if(x->zeroNextSlot)
		{
		  //x->eChanged = 0;
		  x->slot = x->nextSlot;
		  x->PStep = 0;
		  x->PGcyc = 0;
		  x->zeroNextSlot = 0;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ, x->SEQSIZE, x->seq.eOff);
		  //x->pageNum = 0;
		  //x->pageFlag = 1;
		  getVariables(x);
		  x->scrambling = 0;
		  checkJoinsOut(x);
		  //add join code in here!
		}
	      else if(x->zeroNextVar)
		{
		  x->slot = x->nextSlot;
		  x->PStep = 0;
		  x->PGcyc = 0;
		  x->zeroNextVar = 0;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.eOff);
		  x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.varOff);
		  x->scrambling = 1;
		  //x->pageNum = 0;
		  //x->pageFlag = 1;
		  getVariations(x);
		  checkJoinsOut(x);
		}
	      else if((int)x->PGcyc < x->grp.cycles[x->slot] - 1)
		{
		  x->PGcyc++;
		  if(x->changeSlot == 1)
		    {
		      if(x->InVal + x->PGcyc >= x->nextShotVal)
			{
			  x->slot = x->nextSlot;
			  x->PGcyc -= x->wrapSubVal;
			  x->PStep = x->NStep;
			  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eOff);
			  getVariables(x);
			  checkJoinsOut(x);
			  x->changeSlot = 0;
			  x->scrambling = 0;
			}
		    }
		  else if(x->changeVar == 1)
		    {
		      if(x->InVal + x->PGcyc >= x->nextShotVal)
			{
			  x->slot = x->nextSlot;
			  x->PGcyc -= x->wrapSubVal;
			  x->PStep = x->NStep;
			  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eOff);
			  x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.varOff);
			  getVariations(x);
			  checkJoinsVarOut(x);
			  x->changeVar = 0;
			  x->scrambling = 1;
			}
		    }
		  else if(x->InVal + x->PGcyc >= x->PEOff + x->PESize)
		    {
		      x->eChanged = 0;
		      if(x->altOut) x->altNum = !x->altNum;
		      //START November 2nd version 
		      x->PStep++;
		      x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eOff);
		      if(x->myBug == 7) post("PEOff = %f",x->PEOff);
		      getVariables(x);
		      checkJoinsOut(x);
		      //} // END November 2nd version
		    }
		}
	      else
		{ //x->instant = 0; // added Jan 6th 2018
		  x->eChanged = 0;
		  x->PGcyc = 0;
		  x->PEOff = 0;
		  x->PStep = 0;
		  x->PJoined = 0;
		  x->JoinVal = 1; // just in case...but beware the potential source of a bug!
		  x->jFlag = 0;
		  x->jFirst = 0;
                  if(x->jumpSlotAtEnd)
		    {
		      x->slot = x->nextSlot;
		      x->jumpSlotAtEnd = 0;
		      x->barNew = 1;
		      if(x->altOut) x->altNum = !x->altNum;
		      x->scrambling = 0;
		      getVariables(x);
		      checkJoinsOut(x);
		    }
		  else if(x->jumpVarAtEnd)
		    {
		      x->slot = x->nextSlot;
		      x->varPerf = x->nextVar - 1;
		      x->barNew = 1;
		      x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE, x->VARSIZE, x->var.eOff);
		      x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE, x->VARSIZE, x->var.varOff);
		      if(x->altOut) x->altNum = !x->altNum;
		      getVariations(x);
		      x->scrambling = 1;
		      checkJoinsVarOut(x);
		    }
		  else
		    {
		      x->barNew = 1;
		      if(x->altOut) x->altNum = !x->altNum;
		      getVariables(x);
		      checkJoinsOut(x);
		    }
		}
	    }
	  else
	    if(x->changeSlot == 1)
	      {
		if(x->InVal + x->PGcyc >= x->nextShotVal)
		  {
		    x->slot = x->nextSlot;
		    x->PGcyc -= x->wrapSubVal;
		    x->PStep = x->NStep;
		    x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eOff);
		    getVariables(x);
		    checkJoinsOut(x);
		    x->changeSlot = 0;
		    x->scrambling = 0;
		  }
	      }
	    else if(x->changeVar == 1)
	      {
		if(x->InVal + x->PGcyc >= x->nextShotVal)
		  {
		    x->slot = x->nextSlot;
		    x->PGcyc -= x->wrapSubVal;
		    x->PStep = x->NStep;
		    x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eOff);
		    x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.varOff);
		    getVariations(x);
		    checkJoinsVarOut(x);
		    x->changeVar = 0;
		    x->scrambling = 1;
		  }
	      }
	    else if(x->InVal + x->PGcyc >= x->PEOff + x->PESize)
	      { //START November 2nd version
		    //x->instant = x->InVal; // added Jan 6th 2018
		x->PStep++;
		x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eOff);
		if(x->myBug == 7) post("PEOff = %f",x->PEOff);
		getVariables(x);
		if(x->altOut) x->altNum = !x->altNum;
		checkJoinsOut(x);
	      }
	}
      else // PERF_SCRAMBLING -- 13th Jan 2018: we need to add the nextScramble etc code from _nextSlot to this also
	{
	  if(x->InVal < x->PreVal)
	    {
	      if(x->zeroNextPhase)
		{
		  //x->eChanged = 0;
		  //x->slot = x->nextSlot;
		  x->PStep = 0;
		  x->PGcyc = 0;
		  x->zeroNextPhase = 0;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.eOff);
		  x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.varOff);
		  x->scrambling = 1;
		  getVariations(x);
		  checkJoinsVarOut(x);
		  //add join code in here!
		}
	      else if(x->zeroNextSlot)
		{
		  x->slot = x->nextSlot;
		  x->PStep = 0;
		  x->PGcyc = 0;
		  x->zeroNextSlot = 0;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ, x->SEQSIZE, x->seq.eOff);
		  getVariables(x);
		  checkJoinsOut(x);
		  x->scrambling = 0;
		  //add join code in here!
		}
	      else if(x->zeroNextVar)
		{
		  x->slot = x->nextSlot;
		  x->PStep = 0;
		  x->PGcyc = 0;
		  x->zeroNextVar = 0;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.eOff);
		  x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE, x->VARSIZE, x->var.varOff);
		  x->scrambling = 1;
		  getVariations(x);
		  checkJoinsVarOut(x);
		}
	      else if(x->PGcyc + 1 == x->vGrp.cycles[x->slot + x->varPerf * SLOTS])
		{ // reset to 0
		  x->PStep = 0;
		  x->PGcyc = 0;
		  x->PJoined = 0;
		  x->JoinVal = 1; // just in case...but beware the potential source of a bug!
		  x->jFlag = 0;
		  x->jFirst = 0;
                  if(x->jumpSlotAtEnd)
		    {
		      x->slot = x->nextSlot;
		      x->jumpSlotAtEnd = 0;
		      x->barNew = 1;
		      if(x->altOut) x->altNum = !x->altNum;
		      x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE, x->VARSIZE, x->var.eOff);
		      x->scrambling = 0;
		      getVariables(x);
		      checkJoinsOut(x);
		    }
		  else if(x->jumpVarAtEnd)
		    {
		      x->slot = x->nextSlot;
		      x->varPerf = x->nextVar - 1;
		      x->barNew = 1;
		      x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE, x->VARSIZE, x->var.eOff);
		      x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE, x->VARSIZE, x->var.varOff);
		      if(x->altOut) x->altNum = !x->altNum;
		      getVariations(x);
		      x->scrambling = 1;
		      checkJoinsVarOut(x);
		    }
		  else
		    {
		      x->barNew = 1;
		      if(x->altOut) x->altNum = !x->altNum;
		      x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE, x->VARSIZE, x->var.eOff);
		      x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE, x->VARSIZE, x->var.varOff);
		      getVariations(x);
		      checkJoinsVarOut(x);
		    }
		}
	      else
		{
		  x->PGcyc++;
		  if(x->changeSlot == 1)
		    {
		      if(x->InVal + x->PGcyc >= x->nextShotVal)
			{
			  x->slot = x->nextSlot;
			  x->PGcyc -= x->wrapSubVal;
			  x->PStep = x->NStep;
			  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eOff);
			  getVariables(x);
			  checkJoinsOut(x);
			  x->changeSlot = 0;
			  x->scrambling = 0;
			}
		    }
		  else if(x->changeVar == 1)
		    {
		      if(x->InVal + x->PGcyc >= x->nextShotVal)
			{
			  x->slot = x->nextSlot;
			  x->PGcyc -= x->wrapSubVal;
			  x->PStep = x->NStep;
			  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eOff);
			  x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.varOff);
			  getVariations(x);
			  checkJoinsVarOut(x);
			  x->changeVar = 0;
			  x->scrambling = 1;
			}
		    }
		  else if(x->InVal + x->PGcyc >= x->VOff + x->PESize)
		    {
		      x->PStep++;
		      x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eOff);
		      x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE + x->PStep, x->VARSIZE, x->var.varOff);
		      getVariations(x);
		      if(x->altOut) x->altNum = !x->altNum;
		      checkJoinsVarOut(x);
		    }
		}
	    }
	  else if(x->changeSlot == 1)
	    {
	      if(x->InVal + x->PGcyc >= x->nextShotVal)
		{
		  x->slot = x->nextSlot;
		  x->PGcyc -= x->wrapSubVal;
		  x->PStep = x->NStep;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->PStep, x->SEQSIZE, x->seq.eOff);
		  getVariables(x);
		  checkJoinsOut(x);
		  x->changeSlot = 0;
		  x->scrambling = 0;
		}
	    }
	  else if(x->changeVar == 1)
	    {
	      if(x->InVal + x->PGcyc >= x->nextShotVal)
		{
		  x->slot = x->nextSlot;
		  x->PGcyc -= x->wrapSubVal;
		  x->PStep = x->NStep;
		  x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eOff);
		  x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf * x->SEQSIZE + x->PStep, x->VARSIZE, x->var.varOff);
		  getVariations(x);
		  checkJoinsVarOut(x);
		  x->changeVar = 0;
		  x->scrambling = 1;
		}
	    }
	  else if(x->InVal + x->PGcyc >= x->VOff + x->PESize)
	    {
	      x->PStep++;
	      x->PEOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE + x->PStep, x->VARSIZE, x->var.eOff);
	      x->VOff = atom_getfloatarg(x->slot * MAXSEQ + x->varPerf + x->SEQSIZE + x->PStep, x->VARSIZE, x->var.varOff);
	      getVariations(x);
	      if(x->altOut) x->altNum = !x->altNum;
	      checkJoinsVarOut(x);
	    }
	}
      x->TotVal = x->InVal + x->PGcyc;
      if(!x->scrambling)
	{
	  x->eVal = x->TotVal - x->PEOff;
	  x->eVVal = x->eVal * x->JPESI;
	  if(!x->eChanged) // we'll reset x->eChanged in getVariations and getVariables so that it always happens on event reset
	    {
	      if(x->eVVal > x->percentVal)
		{
		  x->eChanged = 1;
		  clock_delay(x->early,0L);
		}
	    }
	  if(x->altOut)
	    {
	      if(!x->altNum)
		{
		  if(x->eMult) *out3++ = x->eVVal;
		  else *out3++ = x->eVal;
		  *alt3++ = 0;
		}
	      else
		{
		  *out3++ = 0;
		  if(x->eMult) *alt3++ = x->eVVal;
		  else *alt3++ = x->eVal;
		}
	    }
	  else
	    {
	      if(x->eMult) *out3++ = x->eVVal;
	      else *out3++ = x->eVal;
	      *alt3++ = 0;
	    }
	}
      else
	{
	  x->eVal = x->TotVal - x->VOff;
	  x->eVVal = x->eVal * x->VPESI;
	  if(!x->eChanged)
	    {
	      if(x->eVVal > x->percentVal)
		{
		  x->eChanged = 1;
		  clock_delay(x->early,0L);
		}
	    }
	  if(x->altOut)
	    {
	      if(!x->altNum)
		{
		  if(x->eMult) *out3++ = x->eVVal;
		  else *out3++ = x->eVal;
		  *alt3++ = 0;
		}
	      else
		{
		  *out3++ = 0;
		  if(x->eMult) *alt3++ = x->eVVal;
		  else *alt3++ = x->eVal;
		}
	    }
	  else
	    {
	      if(x->eMult) *out3++ = x->eVVal;
	      else *out3++ = x->eVal;
	      *alt3++ = 0;
	    }
	}
      *offset++ = x->PEOff;
      x->PreVal = x->InVal;
      if((int)x->TotVal != x->pageNum)
	{
	  x->pageNum = (int)x->TotVal;
	  clock_delay(x->pageTurner, 0L);
	}
    }
  return(w+7);
}

void polyMath_tilde_dsp(t_polyMath_tilde *x, t_signal **sp)
{
  dsp_add(polyMath_tilde_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

static void *polyMath_tilde_new()
{
  t_polyMath_tilde *x = (t_polyMath_tilde *)pd_new(polyMath_tilde_class);
  
  outlet_new(&x->x_obj, gensym("signal"));
  outlet_new(&x->x_obj, gensym("signal"));
  outlet_new(&x->x_obj, gensym("signal"));

  x->clock = outlet_new(&x->x_obj, &s_float);
  x->subclock = outlet_new(&x->x_obj, &s_float);
  x->cycle = outlet_new(&x->x_obj, &s_float);
  x->newgroup = outlet_new(&x->x_obj, &s_bang);
  x->newbar = outlet_new(&x->x_obj, &s_bang);
  x->p1 = outlet_new(&x->x_obj, &s_list);
  x->p2 = outlet_new(&x->x_obj, &s_list);
  x->p3 = outlet_new(&x->x_obj, &s_list);
  x->p4 = outlet_new(&x->x_obj, &s_list);
  x->p5 = outlet_new(&x->x_obj, &s_list);
  x->p6 = outlet_new(&x->x_obj, &s_list);
  x->p7 = outlet_new(&x->x_obj, &s_list);
  x->p8 = outlet_new(&x->x_obj, &s_list);
  x->groupnum = outlet_new(&x->x_obj, &s_float);
  x->num = outlet_new(&x->x_obj, &s_float);
  x->denom = outlet_new(&x->x_obj, &s_float);
  x->eventLengthPhase = outlet_new(&x->x_obj, &s_float);
  x->eventLengthNum = outlet_new(&x->x_obj, &s_float);
  x->alt = outlet_new(&x->x_obj, &s_float);
  x->eChange = outlet_new(&x->x_obj, &s_list);
  x->eAlt = outlet_new(&x->x_obj, &s_float);
  x->dataOut = outlet_new(&x->x_obj, &s_list);
  x->dType = outlet_new(&x->x_obj, &s_float);
  x->page = outlet_new(&x->x_obj, &s_float);
  //*durFirst, *durAlt;
  x->durFirst = outlet_new(&x->x_obj, &s_list);
  x->durAlt = outlet_new(&x->x_obj, &s_list);
  // counters init:
  x->a = x->b = x->c = x->d = x->e = x->g = x->h = x->i = x->j = x->k = 0;// x->l = x->m = x->n = x->o = x->p = x->q = x->v = x->w = x->x = x->y = x->z = 0;

  //31/10/2017 - need to re-evaluate this list. Re-write from scratch so that only used variables are initialized
  // then weed out the unused variables to save memory
  // then re-write polyMath_tilde_init
  x->Gstep = 0;
  x->Gn = 16;
  x->Gd = 16;
  x->ESize = 0.0625;
  x->Eoff = 0;
  x->ESInv = 16;
  x->GSize = 1;
  x->GSInv = 1;
  x->Goff = 0;
  x->cycles = 1;

  x->firstStart = 1;
  
  x->InVal = x->PreVal = 0;
  x->PGcyc = 0;
  // PGcyc, cycles, PEOff, PESize, PJoined, jFlag, jFirst, PJoin, JESize, JPESI = 1 / JESize
  x->PEOff = 0;
  x->JESize = x->PESize = 0.0625;
  x->JPESI = 16;
  x->PJoined = x->jFlag = x->jFirst = 0;
  x->JoinVal = 1;
  
  x->slot = x->group = x->Gstp = 0;
  x->Gnm = 16;
  x->Joined = 0;
  x->join = 1;

  x->PSlot = 0;
  x->PJoin = 1;
  x->Pthis = 0.0625;
  x->PStep = 0;
  
  x->altOut = 0;
  x->altNum = 0;
  x->altLen = 1;

  x->percentVal = 0.8;
  //  x->altEarly = 0;
  x->eMult = 1;

  x->SEQSIZE = SLOTS * MAXSEQ;
  x->GROUPSIZE = SLOTS * GROUPS;
  x->VARSIZE = VARIATIONS * x->SEQSIZE;
  x->VGROUPSIZE = VARIATIONS * x->GROUPSIZE;
  
  unsigned short int seed1 = 12345;
  unsigned short int seed2 = 28374;
  x->timeSeed = time(NULL);
  x->timeSeed = x->timeSeed % 65536;
  x->seed16v[2] = (unsigned short int)x->timeSeed;
  x->seed16v[1] = seed2;
  x->seed16v[0] = seed1;
  seed48(x->seed16v);

  x->scramMeth = 0;
  x->scramSlot = 0;
  x->copyWell = 0;
  x->swapWell = 0;
  x->scramWell = 0;
  x->groupWell = 0;
  x->varOff = 0;
  x->instant = 0;
  x->gInstant = 0;
  x->newVar = 0;
  x->varWrite = 0;
  x->lastVar = 0;

  x->zeroNextPhase = 0;
  x->zeroNextVar = 0;
  x->zeroNextSlot = 0;
  x->changeVar = 0;
  x->changeSlot = 0;
  x->NStep = 0;
  
  x->nextSlot = 0;
  x->lastVar = 0;
  x->noRepeats = 0;
  
  x->validJumpState = 0;
  x->jumpSlotAtEnd = 0;	
  x->jumpVarAtEnd = 0;

  x->cycleDiff = 0;
  x->sizeThreshold = 0.00001;
  x->autoThreshold = 0;
  x->sizeFrac = 0.5;
  
  x->thisSlot = x->slot = 0;
  x->pageNum = -1;

  x->BPM = 60;
  x->durBeat = 1000;
  x->barBeat = 4000;
  x->dur1 = x->dur2 = 250;
  
  for(x->a = 0; x->a < SLOTS; x->a++)
    {
      x->grp.isUnFilled[x->a] = 1;
      x->grp.nGroups[x->a] = 1;
      x->grp.cycles[x->a] = 1;
      for(x->b = 0; x->b < MAXSEQ; x->b++)
	if(x->b < GROUPS)
	  {
	    if(x->a == 0)
	      {
		SETFLOAT(&x->grp.n[x->a * x->b],16);
		SETFLOAT(&x->grp.d[x->a * x->b],16);
		SETFLOAT(&x->grp.offset[x->a * x->b],0);
		SETFLOAT(&x->grp.size[x->a * x->b],1);
		SETFLOAT(&x->grp.sizeInv[x->a * x->b],1);
		x->grp.gStart[x->a * x->b] = 0;
		//SETFLOAT(&x->grp.cycles[x->a * x->b],1);
		//SETFLOAT(&x->grp.rLength[x->a * x->b],1);
		//SETFLOAT(&x->grp.remains[x->a * x->b],0);
	      }
	    else
	      {
		SETFLOAT(&x->grp.n[x->a * x->b],1);
		SETFLOAT(&x->grp.d[x->a * x->b],1);
		SETFLOAT(&x->grp.offset[x->a * x->b],0);
		SETFLOAT(&x->grp.size[x->a * x->b],1);
		SETFLOAT(&x->grp.sizeInv[x->a * x->b],1);
		x->grp.gStart[x->a * x->b] = 0;
		//SETFLOAT(&x->grp.cycles[x->a * x->b],0);
		//SETFLOAT(&x->grp.rLength[x->a * x->b],0);
		//SETFLOAT(&x->grp.remains[x->a * x->b],0);
	      }
	  }
      if(x->b < 16)
	{
	  SETFLOAT(&x->seq.allStep[x->a * x->b],x->b);
	  SETFLOAT(&x->seq.groupStep[x->a * x->b],x->b);
	  SETFLOAT(&x->seq.allStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.eSize[x->a * x->b],0.0625);
	  SETFLOAT(&x->seq.eOff[x->a * x->b],0.0625 * (t_float)x->b);
	  SETFLOAT(&x->seq.eJoin[x->a * x->b],1);
	  SETFLOAT(&x->seq.eAcc1[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc2[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc3[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc4[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc5[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc6[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc7[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc8[x->a * x->b],0);
	  SETFLOAT(&x->seq.pAcc1[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc2[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc3[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc4[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc5[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc6[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc7[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc8[x->a * x->b],-1);
	  SETFLOAT(&x->seq.eSizeInv[x->a * x->b],16);
	}
      else
	{
	  SETFLOAT(&x->seq.allStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.groupStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.allStep[x->a * x->b],0);
	  SETFLOAT(&x->seq.eSize[x->a * x->b],0);
	  SETFLOAT(&x->seq.eOff[x->a * x->b],0);
	  SETFLOAT(&x->seq.eJoin[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc1[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc2[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc3[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc4[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc5[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc6[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc7[x->a * x->b],0);
	  SETFLOAT(&x->seq.eAcc8[x->a * x->b],0);
	  SETFLOAT(&x->seq.pAcc1[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc2[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc3[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc4[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc5[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc6[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc7[x->a * x->b],-1);
	  SETFLOAT(&x->seq.pAcc8[x->a * x->b],-1);
	  SETFLOAT(&x->seq.eSizeInv[x->a * x->b],0);
	}
    }
  //int swaps[MAXSEQ];
  //int swapsRef[MAXSEQ * 2];
  //int swapped[MAXSEQ];
  for(x->s = 0; x->s < MAXSEQ; x->s++)
    {
      SETFLOAT(&x->outList[x->s],0);
      x->vGrp.swaps[x->s] = MAXSEQ * -1;
      x->vGrp.swapsRef[x->s] = -1;
      x->vGrp.swapsRef[x->s + MAXSEQ] = -2;
      x->vGrp.swapped[x->s] = 0;
    }
  for(x->s = 0; x->s < x->VARSIZE; x->s++)
    {
      SETFLOAT(&x->var.allStep[x->s],0);
      SETFLOAT(&x->var.groupStep[x->s],0);
      SETFLOAT(&x->var.allStep[x->s],0);
      SETFLOAT(&x->var.eSize[x->s],0);
      SETFLOAT(&x->var.eOff[x->s],0);
      SETFLOAT(&x->var.eJoin[x->s],1);
      SETFLOAT(&x->var.eAcc1[x->s],0);
      SETFLOAT(&x->var.eAcc2[x->s],0);
      SETFLOAT(&x->var.eAcc3[x->s],0);
      SETFLOAT(&x->var.eAcc4[x->s],0);
      SETFLOAT(&x->var.pAcc1[x->s],-1);
      SETFLOAT(&x->var.pAcc2[x->s],-1);
      SETFLOAT(&x->var.pAcc3[x->s],-1);
      SETFLOAT(&x->var.pAcc4[x->s],-1);
      SETFLOAT(&x->var.eAcc5[x->s],0);
      SETFLOAT(&x->var.eAcc6[x->s],0);
      SETFLOAT(&x->var.eAcc7[x->s],0);
      SETFLOAT(&x->var.eAcc8[x->s],0);
      SETFLOAT(&x->var.pAcc5[x->s],-1);
      SETFLOAT(&x->var.pAcc6[x->s],-1);
      SETFLOAT(&x->var.pAcc7[x->s],-1);
      SETFLOAT(&x->var.pAcc8[x->s],-1);
      SETFLOAT(&x->var.eSizeInv[x->s],16);
      SETFLOAT(&x->var.varOff[x->s],0.0625 * (t_float)x->s);
    }
  for(x->t = 0; x->t < x->VGROUPSIZE; x->t++)
    {
      SETFLOAT(&x->vGrp.n[x->t],0);
      SETFLOAT(&x->vGrp.d[x->t],0);
      SETFLOAT(&x->vGrp.offset[x->t],0);
      SETFLOAT(&x->vGrp.size[x->t],0);
      SETFLOAT(&x->vGrp.sizeInv[x->t],0);
      x->var.variations[x->t] = 0;
    }
  x->fOut = clock_new(x, (t_method)polyMath_tilde_cout);
  x->early = clock_new(x, (t_method)polyMath_tilde_eChange);
  x->pageTurner = clock_new(x, (t_method)polyMath_tilde_pageTurn);
  x->barNew = 0;
  getVariables(x);
  return (x);
}

void polyMath_tilde_setup(void)
{
  polyMath_tilde_class = class_new(gensym("polyMath~"), (t_newmethod)polyMath_tilde_new, 
			      0, sizeof(t_polyMath_tilde), 0, A_GIMME, 0);
    CLASS_MAINSIGNALIN(polyMath_tilde_class, t_polyMath_tilde, f);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_dsp, gensym("dsp"), 0);

    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_init, gensym("init"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_slot, gensym("slot"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_setGroups, gensym("setGroups"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_addGroup, gensym("addGroup"), A_GIMME, 0);

    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_groupInSlot, gensym("groupInSlot"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_groupThisSlot, gensym("groupThisSlot"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_thisSlot, gensym("thisSlot"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_jumpTo, gensym("jumpTo"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_jumpNext, gensym("jumpNext"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_autoThreshold, gensym("autoThreshold"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_sizeThreshold, gensym("sizeThreshold"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_sizeFrac, gensym("sizeFrac"), A_DEFFLOAT, 0);

    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_joinSeq, gensym("setJoins"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_makeJoin, gensym("makeJoin"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_setP, gensym("pSet"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_setPOnly, gensym("pSetOnly"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_setVOnly, gensym("vSetOnly"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_altOut, gensym("altOut"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_preOut, gensym("precent"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_eMult, gensym("eMult"), A_DEFFLOAT, 0);

    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_variation, gensym("variation"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_scramble, gensym("scramble"), A_GIMME, 0);

    class_addmethod(polyMath_tilde_class, (t_method)getVariables, gensym("getVariables"), A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_debug, gensym("debug"), A_DEFFLOAT, 0);

    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_swapElement, gensym("swap"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_seqInSlot, gensym("seqUnit"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_initSeqSlot, gensym("initSeqSlot"), A_DEFFLOAT, A_DEFFLOAT, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_groupScramble, gensym("scrambleGroups"), A_GIMME, 0); //not finished!
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_noRepeats, gensym("noRepeats"), A_DEFFLOAT, 0); //not finished!

    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_getSeq, gensym("getSequence"), A_GIMME, 0);
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_slotLen, gensym("slotLength"), A_GIMME, 0);
    
    class_addmethod(polyMath_tilde_class, (t_method)polyMath_tilde_setBpm, gensym("bpm"), A_DEFFLOAT, 0);
}
