/*
 * SuffixArray2.cpp
 *
 *  Created on: May 16, 2013
 *      Author: jinzhang
 */

#include "SuffixArray2.h"
#include "SuffixArray.h"

////////////// min, max etc. //////////////////////////////////////

#ifndef Max
#define Max(x,y) ((x)>=(y)?(x):(y))
#endif

#ifndef Min
#define Min(x,y) ((x)<=(y)?(x):(y))
#endif

#ifndef Abs
#define Abs(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef PI
#define PI 3.1415927
#endif

// is this the right definition of limit?
inline double limit(double x, double bound)
{
  if      (x >  bound) { return  bound; }
  else if (x < -bound) { return -bound; }
  else                   return x;
}

inline bool leq(int a1, int a2,   int b1, int b2) { // lexic. order for pairs
  return(a1 < b1 || a1 == b1 && a2 <= b2);
}                                                   // and triples
inline bool leq(int a1, int a2, int a3,   int b1, int b2, int b3) {
  return(a1 < b1 || a1 == b1 && leq(a2,a3, b2,b3));
}
// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(int* a, int* b, int* r, uint32_t n, uint32_t K)
{ // count occurrences
  int* c = new int[K + 1];                          // counter array
  for (uint32_t i = 0;  i <= K;  i++) c[i] = 0;         // reset counters
  for (uint32_t i = 0;  i < n;  i++) c[r[a[i]]]++;    // count occurences
  for (uint32_t i = 0, sum = 0;  i <= K;  i++) { // exclusive prefix sums
     int t = c[i];  c[i] = sum;  sum += t;
  }
  for (uint32_t i = 0;  i < n;  i++) b[c[r[a[i]]]++] = a[i];      // sort
  delete [] c;
}

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void suffixArray(int* s, int* SA, uint32_t n, uint32_t K) {
  uint32_t n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2;
  int* s12  = new int[n02 + 3];  s12[n02]= s12[n02+1]= s12[n02+2]=0;
  int* SA12 = new int[n02 + 3]; SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
  int* s0   = new int[n0];
  int* SA0  = new int[n0];

  // generate positions of mod 1 and mod  2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (uint32_t i=0, j=0;  i < n+(n0-n1);  i++) if (i%3 != 0) s12[j++] = i;

  // lsb radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, s+2, n02, K);
  radixPass(SA12, s12 , s+1, n02, K);
  radixPass(s12 , SA12, s  , n02, K);

  // find lexicographic names of triples
  int name = 0, c0 = -1, c1 = -1, c2 = -1;
  for (uint32_t i = 0;  i < n02;  i++) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) {
      name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }

  // recurse if names are not yet unique
  if (name < n02) {
    suffixArray(s12, SA12, n02, name);
    // store unique names in s12 using the suffix array
    for (uint32_t i = 0;  i < n02;  i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (uint32_t i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i;

  // stably sort the mod 0 suffixes from SA12 by their first character
  for (uint32_t i=0, j=0;  i < n02;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, s, n0, K);

  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (uint32_t p=0,  t=n0-n1,  k=0;  k < n;  k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    int i = GetI(); // pos of current offset 12 suffix
    int j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ?
        leq(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
        leq(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
    { // suffix from SA12 is smaller
      SA[k] = i;  t++;
      if (t == n02) { // done --- only SA0 suffixes left
        for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
      }
    } else {
      SA[k] = j;  p++;
      if (p == n0)  { // done --- only SA12 suffixes left
        for (k++;  t < n02;  t++, k++) SA[k] = GetI();
      }
    }
  }
  delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
}




SuffixArray2::SuffixArray2() {
	// TODO Auto-generated constructor stub
	array=NULL;
}

int SuffixArray2::builtArray(char* tmp, uint32_t length) {
/*
    int* s = new int[length+3];
    for(uint32_t i=0;i<length;i++)
    {
    	switch(tmp[i])
    	{
    		case 'A':
    			s[i]=2;
    			break;
    		case 'C':
    			s[i]=3;
    			break;
    		case 'G':
    			s[i]=4;
    			break;
    		case 'T':
    			s[i]=6;
    			break;
    		case 'N':
    			s[i]=5;
    			break;
    		case '$':
    			s[i]=1;
    			break;
    		default:
    			cout<<"there is a letter not in ACGTN$";
    			exit(1);
    			break;
    	}
    }
    array = new int[length+3];
    s[length] = s[length+1] = s[length+2] = array[length] = array[length+1] = array[length+2] = 0;
	suffixArray(s,array,length,20);

	delete [] s;
*/
	//SuffixArray	
	    string input(tmp);
            SuffixArray<uint32_t> sa(input);
	    array = new int[length];
	    for(int i=0;i<length;i++)
		{
			array[i]=sa[i];
		}		
	return 0;
}

int SuffixArray2::getSA(int i) {
	return array[i];
}

SuffixArray2::~SuffixArray2() {
	// TODO Auto-generated destructor stub
	if(array!=NULL)
		delete [] array;
}

