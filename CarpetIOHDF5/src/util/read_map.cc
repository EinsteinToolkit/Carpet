#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cassert>

using namespace std;

int main(void)
{
  vector<char> content;
  int sz = -1;

  if (1) {

    FILE *fh = stdin;

    if (fh) {
      fseek(fh, 0, SEEK_END);
      sz = ftell(fh);
      fseek(fh, 0, SEEK_SET);

      content.resize(sz);

      fread(&content[0], sz, 1, fh);

      fclose(fh);
    }
  }

  if (sz > 0) {
    assert(sz % sizeof(int) == 0);
    int nelems = sz / sizeof(int);
    const int* data = (const int*)(&content[0]);
    for(int i = 0 ; i < nelems ; ) {
      int filenum;
      int lenvname, lenpatchname;

      if (data[i++] != 0x4444) {
        printf("Invalid data found at %d: %d\n", 4*(i-1), data[i-1]);
        exit(1);
      }
      filenum = data[i++];

      int map = data[i++];
      int mglevel = data[i++];
      int reflevel = data[i++];
      int timestep = data[i++];
      int timelevel = data[i++];
      int component = data[i++];
      int rank = data[i++];
      int iorigin[3], ioffset[3], ioffsetdenom[3];
      iorigin[0] = data[i++];
      iorigin[1] = data[i++];
      iorigin[2] = data[i++];
      ioffset[0] = data[i++];
      ioffset[1] = data[i++];
      ioffset[2] = data[i++];
      ioffsetdenom[0] = data[i++];
      ioffsetdenom[1] = data[i++];
      ioffsetdenom[2] = data[i++];
      vector<int> shape;
      shape.resize(rank);
      shape.assign(data+i, data+i+rank);
      i+= 3; // file always write three elements
      // skip: double origin[3], delta[3], time
      // skip: int nghosts[3] bbox[6]
      i += 3*2 + 3*2 + 2 + 3 + 6;
      lenvname = data[i++]; // strlen(objectname)
      lenpatchname = data[i++]; // strlen(varname)

      // cache found indices in a map rather than having many many VarIndex
      // lookups
      const char *vname = (const char*)&data[i];
      i += ((lenvname + 1 + sizeof(int)-1) & ~(sizeof(int)-1)) / sizeof(int);

      const char *patchname = (const char*)&data[i];
      i += ((lenpatchname + 1 + sizeof(int)-1) & ~(sizeof(int)-1)) / sizeof(int);

      if(1) {
        printf("still alive... just read patch, at pos %d of %d\n", i, nelems);
        printf("filenum=%d,"
               "vname=%s,map=%d,mglevel=%d,reflevel=%d,timestep=%d,"
               "timelevel=%d,component=%d,rank=%d,iorigin=(%d,%d,%d),"
               "ioffset=(%d,%d,%d),ioffsetdenom=(%d,%d,%d),shape=(%d,%d,%d),"
               "patchname=%s,\n",
               filenum,
               vname, map, mglevel, reflevel, timestep, timelevel, component,
               rank, iorigin[0], iorigin[1], iorigin[2], ioffset[0], ioffset[1],
               ioffset[2], ioffsetdenom[0], ioffsetdenom[1], ioffsetdenom[2],
               (int)shape[0], (int)shape[1], (int)shape[2], patchname);
      }

      assert(i <= nelems);
    }
  }

  return 0;
}

