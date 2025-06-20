/* eidors_var.id: calculate the sha1 hash of variable in memory
 *   The goal is to be able to create a hash index of
 *   files and a quick way to determine whether files are
 *   identical
 *
 *   $Id: eidors_var_id.cpp 6587 2023-04-18 18:03:10Z aadler $

 * Documentation 
 * http://www.mathworks.com/support/tech-notes/1600/1605.html
 * http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_external/f11333.html

 * New version uses xxHash
 * Download from:
https://raw.githubusercontent.com/Cyan4973/xxHash/dev/xxhash.h

 * This code was created to solve the issue that Matlab doesn't save
 *   variables to disk consistently
 * OLD TEST CODE:
   %   for i=1:20; t1=mk_stim_patterns(16, 16, '{ad}','{ad}',{}, 10); save(sprintf('t%02d.mat',i),'t1');end
   % for i in t*.mat ; do dd if=$i  skip=100 ibs=1 | xxd > $i.xxd ; done
   % for i in *.xxd ; do sha1sum $i ; done
   % 6 different saved formats exist.
 * NEW TEST CODE:
   %   for i=1:20; t1=mk_stim_patterns(16, 16, '{ad}','{ad}',{}, 10); disp(eidors_var_id(t1)); end

 */

#include <stdio.h>
#include <string.h>
#include <mex.h>

/*
 * Defines to alow stat
 */
#include <sys/stat.h>
#include <sys/types.h>
/*
 * defines for SHA1.c code
 */

// #define USE_SHA1_HASH

#if (defined WIN32)
  #define LITTLE_ENDIAN_DEF 5
#elif (__BYTE_ORDER == __LITTLE_ENDIAN ) 
  #define LITTLE_ENDIAN_DEF 5
#elif (__BYTE_ORDER == __BIG_ENDIAN ) 
  #define BIG_ENDIAN_DEF 5
#endif

#define IF_NULL_ERR(a) if (!a) { \
    mexErrMsgTxt("Memory allocation problem"); } 
#define IF_BADSTATUS_ERR(a) if (0) {} else if (a) { \
    mexErrMsgTxt("syscall returned bad status"); }
// This would be nice but doesn't work for the default matlab compiler
//  mexErrMsgTxt(__FILE__  __LINE__  "syscall returned bad status",

#ifdef _MSC_VER
typedef unsigned short uint16_t; //sizeof(unsigned short)=2 
typedef unsigned long  uint32_t; //sizeof(unsigned short)=2 
#else
#include <stdint.h>
#endif


#define unsigned_int32 uint32_t

#if 0 /* Modern Matlab doesn't require this */
#include <limits.h>
#if ULONG_MAX == 0xffffffff
/* OK 32-bit */
#elif ULONG_MAX == 0xffffffffffffffff
   #warning Make sure that you are compiling this with: \
                    mex -largeArrayDims flag. \
            otherwise matlab gives a stupid runtime error \
            (this message is always printed)
#else
   #warning "System does not appear to be 32 or 64 bit"
#endif
#endif

// Code from Matlab's explore.c function
#if defined(_LP64) || defined (_WIN64)
#ifdef MX_COMPAT_32
   #warning MEX-files compiled on a 64-bit platform that use sparse array functions need to be compiled using -largeArrayDims.
#endif
#endif




#ifdef USE_SHA1_HASH

typedef struct {
    unsigned_int32 state[5];
    unsigned_int32 count[2];
    unsigned char buffer[64];
} hash_context;
#define HW 5
static void
hash_initial( hash_context * c );

static void
hash_process( hash_context * c, unsigned char * data, unsigned len );

static void
hash_final( hash_context * c, unsigned_int32[HW] );

#else
// Here we use xxHash instead

#define XXH_STATIC_LINKING_ONLY
#define XXH_IMPLEMENTATION   /* access definitions */
#include "xxhash.h"
#define hash_context XXH3_state_t

static void
hash_initial( hash_context * c ){
  XXH3_128bits_reset_withSeed(c, (XXH64_hash_t)0);
}

static void
hash_process( hash_context * c, unsigned char * data, unsigned len ){
   XXH3_128bits_update(c, data, len);
}

static void
hash_final( XXH3_state_t* c, unsigned_int32 *digest )
{
   XXH128_hash_t u128s;
  // Code from xxhash 6248
   u128s = XXH3_128bits_digest(c);
   digest[0]= (unsigned long)(u128s.high64 >> 32);
   digest[1]= (unsigned long)(u128s.high64 & 0xffffffff);
   digest[2]= (unsigned long)(u128s.low64 >> 32);
   digest[3]= (unsigned long)(u128s.low64  & 0xffffffff);

  // Get 4 more bytes, put all the
  // state back into the hash. Since there
  // are > 256 bits in the state, we 
  // can get at them here
  // But the state is not stable beyond the 256 bits, it seems
// printf("[%d]\n",sizeof(XXH3_state_t));
// {int i;for(i=0;i<sizeof(XXH3_state_t);i++) printf("%02X",c[i]);}
// {int i;for(i=0;i<10;i++) printf("%02X",c[i]);}
// XXH3_128bits_update(c, c, sizeof(XXH3_state_t)/4);
   XXH3_128bits_update(c, c, 32);
   u128s = XXH3_128bits_digest(c);
   digest[4]= (unsigned long)(u128s.high64 & 0xffffffff);
}
#endif
static void
recurse_hash( hash_context *c, const mxArray *var );


// This is the core function. Iterate over the variable
//   to calculate the in-memory SHA1 hash
// Processing - test if:
//   1. Numeric (full) -> Add to SHA1 hash (including complex part)
//   2. Numeric (Sparse) -> Add to SHA1 hash (including complex part)
//   3. String (or other)
//         - hash string content
//         - check if string points to a function hash that
//   4. Empty -> Ignore (Assume not relevant) 
//   5. Cell -> recursively call for each element
//   6. Struct -> recursively call ( In SORTED order )
//   7. Function pointer -> Convert to string and add
//
// NOTE: this function does not hash the array size and
//    orientation. That means that a vectorized version of
//    the same array will have the same ID

#define sDBL sizeof(double)
#define sSGL sizeof(float)
#define sINT sizeof(int)
#define sSZT sizeof(size_t)
   #undef VERBOSE 
// #define VERBOSE 1
		
// check to see if a given string points to an function
//   on disk *.m file.  If it does -> get the file modification
//   time
void lookupfiletime( hash_context *c, const mxArray *var ) {

  mxArray * lhs[1];
  // STEP 1: see if function is an M-file

  { // STEP 1A: Call exist
    mxArray * rhs[1];
    rhs[0] = (mxArray *) var; // BULLS**T from Matlab. Why throw away const?
//  rhs[1] = mxCreateString("file");
//  This should save time, but breaks badly. More BS from Matlab?
    
   {int retval = mexCallMATLAB(1,lhs, 1, rhs, "exist");
//  mxDestroyArray( rhs[1] );
    if ( retval != 0 ||
         mxGetNumberOfElements(lhs[0])!=1 ||
         !mxIsNumeric( lhs[0] ) ) {
      // var doesn't point to a function -> leave
      mxDestroyArray( lhs[0] );
      return;
   } } // for C
  }

  #ifdef VERBOSE
  { char buf[300]; //it'll do for debugging
    mxGetString( var, buf, mxGetNumberOfElements(var)+1 );
    mexPrintf("exist('%s') => %3.1f\n", buf, (mxGetPr(lhs[0]))[0] );
  }
  #endif
  // STEP 1A: Check if exist( var ) == 2 (m-file)
  { double * 
    pr  = mxGetPr( lhs[0] );
    if ( !pr ||
          pr[0] != 2.0 ) {
        mxDestroyArray( lhs[0] );
        return;
    }
  }

  // STEP 2: Get the path to the file
  mexCallMATLAB(1,lhs, 1, (mxArray **) &var, "which");
  if ( mxGetNumberOfElements(lhs[0])==0 ) {
    // var doesn't point to a function -> leave
    mxDestroyArray( lhs[0] );
    return;
  }

 {int len= mxGetNumberOfElements( *lhs ) + 1;
 {char * fname= (char *) mxMalloc(len* sizeof(char));
  IF_NULL_ERR( fname );

  IF_BADSTATUS_ERR(
     mxGetString( *lhs, fname, len) );
  mxDestroyArray( lhs[0] );

 {struct stat buffer;
  IF_BADSTATUS_ERR(
     stat(fname, &buffer) );

  #ifdef VERBOSE
  mexPrintf("Got string=%s mtime=%d\n", fname, buffer.st_mtime);
  #endif
  hash_process( c, (unsigned char *) &buffer.st_mtime, 
                   sizeof( buffer.st_mtime ) );

  mxFree( fname );
  } } } // for C
}

void hash_struct( hash_context *c, const mxArray *var )
{
  mxArray * sortord;
  { // STEP 1: Sort field names
    mxArray *lhs1[1], *lhs2[2];
    int retval;
    retval = mexCallMATLAB(1,lhs1, 1, (mxArray **) &var, "fieldnames");
    if ( retval != 0 ) {
      return;
    }

    retval = mexCallMATLAB(2,lhs2, 1, lhs1, "sort");
    if ( retval != 0 ) {
      mxDestroyArray( lhs1[0] );
      return;
    }

    mxDestroyArray( lhs1[0] );
    mxDestroyArray( lhs2[0] );
    sortord = lhs2[1];
  }

  { // Step 2: Recurse through the fields
    int i,j;
    double * s_idx;
    int NF, NE;

    s_idx= mxGetPr( sortord );
    if (!s_idx) {
        mxDestroyArray( sortord );
        return;
    }
    
    NF= mxGetNumberOfFields( var );
    NE= mxGetNumberOfElements( var );

    for (i= 0; i< NF; i++) {
      int k= (int) s_idx[i] -1; // -1 because Matlab uses 1 indexing
      const char * fdname = mxGetFieldNameByNumber(var, k);
      if (fdname == NULL) continue;
      #ifdef VERBOSE
        mexPrintf("processing field ( %s ) [%d->%d]:\n", fdname, i, k);
      #endif
      for (j= 0; j< NE; j++) {
        mxArray * fd = mxGetFieldByNumber( var, j, k );
        if (fd == NULL ) {
      #ifdef VERBOSE
          mexPrintf("empty field(%s,%d):",
                    mxGetFieldNameByNumber(var,k), j+1);
      #endif
        } else if (strcmp("name",fdname)==0) {
      #ifdef VERBOSE
           mexPrintf("field => IGNORE ( %s ) [%d->%d]:\n", fdname, i, k);
      #endif
        } else {
          hash_process( c, (unsigned char *) fdname, strlen( fdname ) );
          recurse_hash(c, fd);
        }
      }
    }
  }
  mxDestroyArray( sortord );
}

static void
recurse_hash( hash_context *c, const mxArray *var ) {

  #ifdef VERBOSE    
      mexPrintf("processing var of ClassID ( %d ):", mxGetClassID( var ) );
  #endif
  if ( var == NULL ) {
    #ifdef VERBOSE
       mexPrintf("ignoring element ( NULL ):");
    #endif
  } else
  if ( mxIsSparse(var) ) {
    // sparse variable. We need to hash the numeric data,
    // as well as the row and col index pointers
    int  zero= 0;
    unsigned int nnz, cols, i; 
//  unsigned char *p_jcs, *p_irs;
    mwIndex *p_jcs, *p_irs;

    p_irs = mxGetIr( var );
    p_jcs = mxGetJc( var );
    cols= mxGetN( var );
    nnz = *(mxGetJc(var) + cols ); /* after last element of jcs */

    #ifdef VERBOSE
{
       mxClassID category = mxGetClassID(var);

      printf("st= %ld nnz=%d\n",sSZT,nnz);

      printf("doing sparse on %ld byte (size_t = %ld): Sparse CLASS=",sizeof(mwIndex), sizeof(size_t));
       switch (category) {
          case mxLOGICAL_CLASS:  printf("mxLOGICAL_CLASS\n");   break;
          case mxCHAR_CLASS:     printf("mxCHAR_CLASS\n");      break;
          case mxSTRUCT_CLASS:   printf("mxSTRUCT_CLASS\n");    break;
          case mxCELL_CLASS:     printf("mxCELL_CLASS\n");      break;
          case mxUNKNOWN_CLASS:  printf("mxUNKNOWN_CLASS\n");   break;
          default:               printf("mxdefault\n");         break;
       }
}
    #endif

    hash_process( c, (unsigned char *) p_jcs, sizeof(mwIndex) * cols );
    hash_process( c, (unsigned char *) p_irs, sizeof(mwIndex) * nnz );
// if matlab introduces more sparse classes, then this will fail
    if( mxGetClassID(var) == mxLOGICAL_CLASS ) { // Logical sparse
          mxLogical *pr; // point to imaginary is impossible in logcal
          pr  = (mxLogical *) mxGetPr( var );
    #ifdef VERBOSE
          printf("sparse logical => first elem [i=%d,j=%d,v=%d]", p_jcs[0],p_irs[0],pr[0]);
    #endif
          hash_process( c, (unsigned char *) pr,  sizeof(mxLogical) * nnz );
    } else { // Double sparse
          double *pr,*pi;
          pr  = mxGetPr( var );
          pi  = mxGetPi( var );
    #ifdef VERBOSE
          printf("sparse double  => first elem [i=%d,j=%d,v=%f]", p_jcs[0],p_irs[0],pr[0]);
    #endif
          hash_process( c, (unsigned char *) pr,  sDBL * nnz );
          if ( pi != NULL ) {
             hash_process( c, (unsigned char *) pi, sDBL * nnz );
          }
    }
  } else
  if ( mxIsLogical(var) ) {
    int len= mxGetNumberOfElements( var );
        #ifdef VERBOSE    
          mexPrintf("Logical ClassID ( %d [%d]):", mxGetClassID( var ), len );
        #endif
    mxLogical *pl =  mxGetLogicals( var );
    hash_process( c, (unsigned char *) pl, len );
  } else
  if ( mxIsNumeric(var) ) {
    // full numeric variable. ) We need to hash the numeric data.
    double *pr,*pi;
    int len= mxGetNumberOfElements( var );
        #ifdef VERBOSE    
          mexPrintf("Numeric ClassID ( %d [%d] ):", mxGetClassID( var ), len );
        #endif
    if( mxIsDouble(var) ) {
        #ifdef VERBOSE    
          if (len>0)
          mexPrintf("DBL len=%d, first=%5.3g\n:", len, *mxGetPr( var ) );
        #endif
        len= sDBL * len;
    } else 
    if( mxIsSingle(var) ) {
        #ifdef VERBOSE    
          if (len>0)
          mexPrintf("SGL len=%d, first=%5.3g\n:", len, *mxGetPr( var ) );
        #endif
        len= sSGL * len;
    } else 
    if( mxIsInt32(var) || mxIsUint32(var) ) {
        len= 4 * len;
    } else 
    if( mxIsInt16(var) || mxIsUint16(var) ) {
        len= 2 * len;
    } else 
    if( mxIsInt8(var) || mxIsUint8(var) ) {
        len= 1 * len;
    } else {
        mexErrMsgTxt("Unrecognized numeric type");
    }

    pr = mxGetPr( var );
    pi = mxGetPi( var );

    hash_process( c, (unsigned char *) pr, len );
    if ( pi != NULL ) {
       hash_process( c, (unsigned char *) pi, len );
    }
  } else
  if ( mxIsChar(var) ) {
    // string variable. Each char is packed into 2 bytes in Matlab
#ifdef OCTINTERP_API /* Are we using octave */
#define CHARBYTES 1
#else
#define CHARBYTES 2
#endif
        #ifdef VERBOSE    
          mexPrintf("Char ClassID ( %d ) CHARBYTES=%d nel=%d:", mxGetClassID( var ), CHARBYTES,
                  mxGetNumberOfElements( var ) ); // gives correct len of UTF8 encoded
        #endif
    double * pr = mxGetPr( var );
    hash_process( c, (unsigned char *) pr,
                  CHARBYTES*mxGetNumberOfElements( var ) );

    // If var is a *.m file, add it's modification time
    lookupfiletime( c, var);

  } else
  if ( mxIsCell(var) ) {
    // cell variable. Iterate through elements and recurse
    unsigned int i;
    #ifdef VERBOSE
      mexPrintf("processing cell ( %d ):\n", mxGetNumberOfElements(var));
    #endif
    for (i= 0;
         i< mxGetNumberOfElements( var );
         i++) {
      mxArray * fd = mxGetCell( var, i );
      recurse_hash(c, fd);
    }
  } else
  if ( mxIsStruct(var) ) {
    #ifdef VERBOSE
      mexPrintf("processing struct\n");
    #endif
    hash_struct( c, var);
  } else
  if ( mxIsFunctionHandle(var) ) {
    // function_handle. Get string of fcn name
    /* I can't find any documentation on getting fcn string in mex */
    mxArray *lhs[1];
    // Complete BULLS*** from Matlab here. You need to call Matlab
    // using a non-const pointer, but you get a const one
    mexCallMATLAB(1,lhs, 1, (mxArray **) &var, "func2str");
    if ( mxIsChar(lhs[0]) ) {
      // string variable. Each char is packed into 2 bytes
      double * pr = mxGetPr( *lhs );
      int len= 2* mxGetNumberOfElements( *lhs );
      hash_process( c, (unsigned char *) pr, len );
      // and file changed time of file - if it is a file on disk
      lookupfiletime( c, lhs[0] );
    } else {
      mexErrMsgTxt("eidors_var_id: weird output for function_handle");
    }
    mxDestroyArray( lhs[0] );
  } else
  {
    #ifdef VERBOSE
      mexPrintf("ignoring var of ClassID ( %d ):", mxGetClassID( var ) );
    #endif
  }


}
// Calculate a variable_id, which represents the // content. This can be used to test whether a 
// particular calculation has been done before.
//
// usage: eidors_var_id( var1, var2)
// output: id_87AB ... [40 hex digits]
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  hash_context c;
  
// the xxHash library recomends 
// mallocing the state like this, but
// static allocs are cleaner and don't
// require checking errors. Timing
// test show it exactly the same
//XXH3_state_t* cp= XXH3_createState();
//XXH3_state_t c = *cp;
//XXH3_freeState(cp);
  
  unsigned_int32 digest[5];

  if ( !nrhs )  {
    mexErrMsgTxt("eidors_var_id: requires at least one input");
  }

  hash_initial( &c);

  { int i;
    for( i=0; i< nrhs; i++ ) 
      recurse_hash( &c, prhs[i] );
  }

  hash_final( &c, digest);

#define SHA1BUFSZ (3 + 8*5 + 1)
  { char *
  sha1buf = (char *) mxMalloc(SHA1BUFSZ * sizeof(char));
  snprintf(sha1buf, SHA1BUFSZ, "id_%08X%08X%08X%08X%08X", 
          digest[0], digest[1], digest[2], digest[3], digest[4] );
  plhs[0] = mxCreateString(sha1buf);
  mxFree( sha1buf );
  }  
  return;
}

#ifdef USE_SHA1_HASH
/*
SHA1 Calculation Code

This code is available from:
   http://ds.dial.pipex.com/george.barwood/v8/pegwit.htm
SHA-1 in C
By Steve Reid <steve@edmweb.com>
100% Public Domain

Test Vectors (from FIPS PUB 180-1)
"abc"
  A9993E36 4706816A BA3E2571 7850C26C 9CD0D89D
"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"
  84983E44 1C3BD26E BAAE4AA1 F95129E5 E54670F1
A million repetitions of "a"
  34AA973C D4C4DAA4 F61EEB2B DBAD2731 6534016F
*/


#if !defined(LITTLE_ENDIAN_DEF) && !defined(BIG_ENDIAN_DEF)
  #if defined(_M_IX86) || defined(_M_I86) || defined(__alpha)
    #define LITTLE_ENDIAN_DEF
  #else
    #error "LITTLE_ENDIAN_DEF or BIG_ENDIAN_DEF must be defined"
  #endif
#endif

#define SHA1HANDSOFF /* Copies data before messing with it. */


#define rol(value, bits) (((value) << (bits)) | ((value) >> (32 - (bits))))

/* blk0() and blk() perform the initial expand. */
/* I got the idea of expanding during the round function from SSLeay */
#ifdef LITTLE_ENDIAN_DEF
#define blk0(i) (block->l[i] = (rol(block->l[i],24)&0xFF00FF00) \
    |(rol(block->l[i],8)&0x00FF00FF))
#else
#define blk0(i) block->l[i]
#endif
#define blk(i) (block->l[i&15] = rol(block->l[(i+13)&15]^block->l[(i+8)&15] \
    ^block->l[(i+2)&15]^block->l[i&15],1))

/* (R0+R1), R2, R3, R4 are the different operations used in SHA1 */
#define R0(v,w,x,y,z,i) z+=((w&(x^y))^y)+blk0(i)+0x5A827999+rol(v,5);w=rol(w,30);
#define R1(v,w,x,y,z,i) z+=((w&(x^y))^y)+blk(i)+0x5A827999+rol(v,5);w=rol(w,30);
#define R2(v,w,x,y,z,i) z+=(w^x^y)+blk(i)+0x6ED9EBA1+rol(v,5);w=rol(w,30);
#define R3(v,w,x,y,z,i) z+=(((w|x)&y)|(w&x))+blk(i)+0x8F1BBCDC+rol(v,5);w=rol(w,30);
#define R4(v,w,x,y,z,i) z+=(w^x^y)+blk(i)+0xCA62C1D6+rol(v,5);w=rol(w,30);


/* Hash a single 512-bit block. This is the core of the algorithm. */

static
void SHA1Transform(unsigned_int32 state[5], unsigned char buffer[64])
{
unsigned_int32 a, b, c, d, e;
typedef union {
    unsigned char c[64];
    unsigned_int32 l[16];
} CHAR64LONG16;
CHAR64LONG16* block;
#ifdef SHA1HANDSOFF
static unsigned char workspace[64];
    block = (CHAR64LONG16*)workspace;
    memcpy(block, buffer, 64);
#else
    block = (CHAR64LONG16*)buffer;
#endif
    /* Copy context->state[] to working vars */
    a = state[0];
    b = state[1];
    c = state[2];
    d = state[3];
    e = state[4];
    /* 4 rounds of 20 operations each. Loop unrolled. */
    R0(a,b,c,d,e, 0); R0(e,a,b,c,d, 1); R0(d,e,a,b,c, 2); R0(c,d,e,a,b, 3);
    R0(b,c,d,e,a, 4); R0(a,b,c,d,e, 5); R0(e,a,b,c,d, 6); R0(d,e,a,b,c, 7);
    R0(c,d,e,a,b, 8); R0(b,c,d,e,a, 9); R0(a,b,c,d,e,10); R0(e,a,b,c,d,11);
    R0(d,e,a,b,c,12); R0(c,d,e,a,b,13); R0(b,c,d,e,a,14); R0(a,b,c,d,e,15);
    R1(e,a,b,c,d,16); R1(d,e,a,b,c,17); R1(c,d,e,a,b,18); R1(b,c,d,e,a,19);
    R2(a,b,c,d,e,20); R2(e,a,b,c,d,21); R2(d,e,a,b,c,22); R2(c,d,e,a,b,23);
    R2(b,c,d,e,a,24); R2(a,b,c,d,e,25); R2(e,a,b,c,d,26); R2(d,e,a,b,c,27);
    R2(c,d,e,a,b,28); R2(b,c,d,e,a,29); R2(a,b,c,d,e,30); R2(e,a,b,c,d,31);
    R2(d,e,a,b,c,32); R2(c,d,e,a,b,33); R2(b,c,d,e,a,34); R2(a,b,c,d,e,35);
    R2(e,a,b,c,d,36); R2(d,e,a,b,c,37); R2(c,d,e,a,b,38); R2(b,c,d,e,a,39);
    R3(a,b,c,d,e,40); R3(e,a,b,c,d,41); R3(d,e,a,b,c,42); R3(c,d,e,a,b,43);
    R3(b,c,d,e,a,44); R3(a,b,c,d,e,45); R3(e,a,b,c,d,46); R3(d,e,a,b,c,47);
    R3(c,d,e,a,b,48); R3(b,c,d,e,a,49); R3(a,b,c,d,e,50); R3(e,a,b,c,d,51);
    R3(d,e,a,b,c,52); R3(c,d,e,a,b,53); R3(b,c,d,e,a,54); R3(a,b,c,d,e,55);
    R3(e,a,b,c,d,56); R3(d,e,a,b,c,57); R3(c,d,e,a,b,58); R3(b,c,d,e,a,59);
    R4(a,b,c,d,e,60); R4(e,a,b,c,d,61); R4(d,e,a,b,c,62); R4(c,d,e,a,b,63);
    R4(b,c,d,e,a,64); R4(a,b,c,d,e,65); R4(e,a,b,c,d,66); R4(d,e,a,b,c,67);
    R4(c,d,e,a,b,68); R4(b,c,d,e,a,69); R4(a,b,c,d,e,70); R4(e,a,b,c,d,71);
    R4(d,e,a,b,c,72); R4(c,d,e,a,b,73); R4(b,c,d,e,a,74); R4(a,b,c,d,e,75);
    R4(e,a,b,c,d,76); R4(d,e,a,b,c,77); R4(c,d,e,a,b,78); R4(b,c,d,e,a,79);
    /* Add the working vars back into context.state[] */
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    /* Wipe variables */
    a = b = c = d = e = 0;
}


/* Initialize new context */

static
void hash_initial(hash_context* context)
{
    /* SHA1 initialization constants */
    context->state[0] = 0x67452301;
    context->state[1] = 0xEFCDAB89;
    context->state[2] = 0x98BADCFE;
    context->state[3] = 0x10325476;
    context->state[4] = 0xC3D2E1F0;
    context->count[0] = context->count[1] = 0;
}


/* Run your data through this.
   We do nothing for empty vectors.
   Thus [] is the same as not having the variable
 */
static
void hash_process( hash_context * context, unsigned char * data, unsigned len )
{
unsigned int i, j;
unsigned_int32 blen = ((unsigned_int32)len)<<3;
#if VERBOSE
  {int i;
   printf("%d:[",len);
   for(i=0;i<len; i++) printf("%02X",data[i]);
   printf("]  ");}
#endif
    if(len==0) return;

    j = (context->count[0] >> 3) & 63;
    if ((context->count[0] += blen) < blen ) context->count[1]++;
    context->count[1] += (len >> 29);
    if ((j + len) > 63) {
        memcpy(&context->buffer[j], data, (i = 64-j));
        SHA1Transform(context->state, context->buffer);
        for ( ; i + 63 < len; i += 64) {
            SHA1Transform(context->state, &data[i]);
        }
        j = 0;
    }
    else i = 0;
    memcpy(&context->buffer[j], &data[i], len - i);
}


/* Add padding and return the message digest. */

static
void hash_final( hash_context* context, unsigned_int32 digest[5] )
{
unsigned_int32 i, j;
unsigned char finalcount[8];

    for (i = 0; i < 8; i++) {
        finalcount[i] = (unsigned char)((context->count[(i >= 4 ? 0 : 1)]
         >> ((3-(i & 3)) * 8) ) & 255);  /* Endian independent */
    }
    hash_process(context, (unsigned char *)"\200", 1);
    while ((context->count[0] & 504) != 448) {
        hash_process(context, (unsigned char *)"\0", 1);
    }
    hash_process(context, finalcount, 8);  /* Should cause a SHA1Transform() */
    for (i = 0; i < 5; i++) {
        digest[i] = context->state[i];
    }
    /* Wipe variables */
    i = j = 0;
    memset(context->buffer, 0, 64);
    memset(context->state, 0, 20);
    memset(context->count, 0, 8);
    memset(&finalcount, 0, 8);
#ifdef SHA1HANDSOFF  /* make SHA1Transform overwrite it's own static vars */
    SHA1Transform(context->state, context->buffer);
#endif
}

#endif // USE_SHA1_HASH
