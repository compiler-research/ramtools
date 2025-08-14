#include "ttree/RAMRecord.h"

RAMRefs *RAMRecord::fgRnameRefs = 0;
RAMRefs *RAMRecord::fgRnextRefs = 0;
RAMIndex *RAMRecord::fgIndex = 0;

TTree *RAMRecord::GetTree(TFile *file, const char *treeName)
{
   if (!file) {
      ::Error("RAMRecord::GetTree", "file was not opened");
      return 0;
   }
   TTree *t = (TTree *)file->Get(treeName);
   ReadAllRefs();
   ReadIndex();
   return t;
}

void RAMRecord::WriteRefs(const RAMRefs *refs, const char *refname)
{
   if (gFile) {
      if (refs)
         gFile->WriteObjectAny(refs, "RAMRefs", refname);
   } else
      ::Error("RAMRecord::WriteRefs", "no file open");
}

void RAMRecord::ReadRefs(RAMRefs *&refs, const char *refname)
{
   if (gFile) {
      auto r = (RAMRefs *)gFile->Get(refname);
      if (refs)
         delete refs;
      refs = r;
   } else
      ::Error("RAMRecord::ReadRefs", "no file open");
}

void RAMRecord::WriteIndex()
{
   if (gFile) {
      if (fgIndex)
         gFile->WriteObjectAny(fgIndex, "RAMIndex", "Index");
   } else
      ::Error("RAMRecord::WriteIndex", "no file open");
}

void RAMRecord::ReadIndex()
{
   if (gFile) {
      auto index = (RAMIndex *)gFile->Get("Index");
      if (fgIndex)
         delete fgIndex;
      fgIndex = index;
   } else
      ::Error("RAMRecord::ReadIndex", "no file open");
}

RAMRefs::RAMRefs()
{
   fLastId = 0;
   fMaxId = 100;
   fRefVec.reserve(fMaxId);
}

int RAMRefs::GetRefId(const char *rname)
{
   // Convert name to a refid.

   if (rname[0] == '*')
      return -1;

   if (fLastName == rname)
      return fLastId;

   auto it = std::find(fRefVec.begin(), fRefVec.end(), rname);
   if (it != fRefVec.end()) {
      auto index = std::distance(fRefVec.begin(), it);
      fLastId = (int)index;
      fLastName = rname;
      return fLastId;
   }

   if (fLastId + 1 >= fMaxId) {
      fMaxId *= 2;
      fRefVec.reserve(fMaxId);
   }

   fRefVec.push_back(rname);

   fLastId = fRefVec.size() - 1;
   fLastName = rname;

   return fLastId;
}

const char *RAMRefs::GetRefName(int rid)
{
   // Convert refid to name.

   if (rid == -1)
      return "*";

   return fRefVec[rid].c_str();
}

void RAMRefs::Print() const
{
   int size = fRefVec.size();
   printf("RAMRefs vector:\n");
   for (int i = 0; i < size; i++)
      printf("%d: %s\n", i, fRefVec[i].c_str());
}

void RAMIndex::AddItem(int refid, int pos, Long64_t row)
{
   Key_t key = std::make_pair(refid, pos);
   fIndex[key] = row;
}

Long64_t RAMIndex::GetRow(int refid, int pos)
{
   Index_t::iterator low;
   Key_t key = std::make_pair(refid, pos);
   low = fIndex.lower_bound(key);

   if (low == fIndex.end()) {
      return -1; // nothing found
   } else if (low == fIndex.begin()) {
      return low->second;
   } else {
      if ((low->first.first == refid) && (low->first.second == pos)) {
         return low->second;
      } else {
         --low;
         return low->second;
      }
   }
}

void RAMIndex::Print() const
{
   Index_t::const_iterator it = fIndex.begin();
   printf("RAMIndex map:\n");
   while (it != fIndex.end()) {
      printf("%lld: refid=%d, pos=%d\n", it->second, it->first.first, it->first.second);
      ++it;
   }
}

