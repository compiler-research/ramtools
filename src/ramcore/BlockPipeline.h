// The pipeline both converters share: the calling thread reads the input in
// blocks, worker threads fill each block into a fill context of their own, and
// the clusters are committed in block order, so the records keep the order of
// the input.
#ifndef RAMCORE_BLOCKPIPELINE_H
#define RAMCORE_BLOCKPIPELINE_H

#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMNTupleRecord.h"
#include "rntuple/RAMTags.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleFillContext.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleParallelWriter.hxx>

#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <exception>
#include <future>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace ramcore::blocks {

struct BlockOrder {
   RAMCoordinateOrder check;
   bool has_placed = false;
   int32_t first_refid = -1;
   int32_t first_pos = -1;
   uint32_t max_span = 0;
   RAMRefCounts counts;
   RAMFlagCounts flags;

   void Note(int32_t refid, int32_t pos)
   {
      if (refid >= 0 && !has_placed) {
         has_placed = true;
         first_refid = refid;
         first_pos = pos;
      }
      check.Note(refid, pos);
   }
};

struct FileOrder {
   RAMCoordinateOrder check;
   uint32_t max_span = 0;
   RAMRefCounts counts;
   RAMFlagCounts flags;

   void Add(const BlockOrder &b)
   {
      counts.Add(b.counts);
      flags.Add(b.flags);
      if (!b.check.sorted)
         check.sorted = false;
      if (b.has_placed) {
         check.Note(b.first_refid, b.first_pos);
         check.Note(b.check.last_refid, b.check.last_pos);
      }
      if (b.check.seen_unplaced)
         check.seen_unplaced = true;
      max_span = std::max(max_span, b.max_span);
   }
};

struct Block {
   size_t seq = 0;        ///< Position in the input, from 0.
   size_t first_line = 0; ///< 1-based number of the block's first line, for warnings.
   std::vector<char> data;
};

// Bounded FIFO. Not TBB (TTaskGroup, TThreadExecutor): it does not start tasks
// in order, and a task waiting for an earlier block could then deadlock.
class BlockQueue {
   std::mutex m_mutex;
   std::condition_variable m_notEmpty;
   std::condition_variable m_notFull;
   std::deque<Block> m_blocks;
   size_t m_capacity;
   bool m_closed = false;

public:
   explicit BlockQueue(size_t capacity) : m_capacity(std::max<size_t>(capacity, 1)) {}

   void Push(Block &&block)
   {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_notFull.wait(lock, [&] { return m_blocks.size() < m_capacity || m_closed; });
      if (m_closed)
         return;
      m_blocks.push_back(std::move(block));
      m_notEmpty.notify_one();
   }

   /// Returns once a Push would not have to wait, so the reader can read its
   /// next block when there is room for it instead of holding it until then.
   void WaitForRoom()
   {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_notFull.wait(lock, [&] { return m_blocks.size() < m_capacity || m_closed; });
   }

   bool Pop(Block &block)
   {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_notEmpty.wait(lock, [&] { return !m_blocks.empty() || m_closed; });
      if (m_blocks.empty())
         return false;
      block = std::move(m_blocks.front());
      m_blocks.pop_front();
      m_notFull.notify_one();
      return true;
   }

   void Close()
   {
      const std::lock_guard<std::mutex> lock(m_mutex);
      m_closed = true;
      m_notEmpty.notify_all();
      m_notFull.notify_all();
   }

   void Abort()
   {
      const std::lock_guard<std::mutex> lock(m_mutex);
      m_closed = true;
      m_blocks.clear();
      m_notEmpty.notify_all();
      m_notFull.notify_all();
   }
};

struct RefCache {
   std::string name;
   int id = -1;
   bool valid = false;

   int Lookup(RAMNTupleRefs &refs, const std::string &rname)
   {
      if (valid && rname == name)
         return id;
      id = refs.GetRefId(rname);
      name = rname;
      valid = true;
      return id;
   }
};

struct Progress {
   std::mutex mutex;
   std::condition_variable next_turn;
   size_t next_seq = 0; ///< The block whose clusters may be committed next.
   bool failed = false; ///< A worker or the reader gave up; everyone stops.
   FileOrder order;
   size_t records = 0;
   std::vector<uint64_t> block_ends;
   std::vector<std::string> late_headers;

   void Fail()
   {
      const std::lock_guard<std::mutex> lock(mutex);
      failed = true;
      next_turn.notify_all();
   }

   bool Failed()
   {
      const std::lock_guard<std::mutex> lock(mutex);
      return failed;
   }
};

// A worker leaving by exception stops the others; the exception itself
// travels in the worker's future.
class StopOthersOnException {
   Progress &m_progress;
   BlockQueue &m_queue;
   const int m_exceptions = std::uncaught_exceptions();

public:
   StopOthersOnException(Progress &progress, BlockQueue &queue) : m_progress(progress), m_queue(queue) {}
   StopOthersOnException(const StopOthersOnException &) = delete;
   StopOthersOnException &operator=(const StopOthersOnException &) = delete;
   StopOthersOnException(StopOthersOnException &&) = delete;
   StopOthersOnException &operator=(StopOthersOnException &&) = delete;
   ~StopOthersOnException()
   {
      if (std::uncaught_exceptions() > m_exceptions) {
         m_progress.Fail();
         m_queue.Abort();
      }
   }
};

class CloseQueueOnExit {
   BlockQueue &m_queue;

public:
   explicit CloseQueueOnExit(BlockQueue &queue) : m_queue(queue) {}
   CloseQueueOnExit(const CloseQueueOnExit &) = delete;
   CloseQueueOnExit &operator=(const CloseQueueOnExit &) = delete;
   CloseQueueOnExit(CloseQueueOnExit &&) = delete;
   CloseQueueOnExit &operator=(CloseQueueOnExit &&) = delete;
   ~CloseQueueOnExit() { m_queue.Abort(); }
};

// The tables every worker fills; each locks.
struct SharedRefs {
   RAMNTupleRefs rnames;
   RAMNTupleRefs rnexts;
   RAMTagTables tags;
};

/// The model of the records for RNTupleParallelWriter, which needs a bare one.
inline std::unique_ptr<ROOT::RNTupleModel> MakeBareModel()
{
   auto model = ROOT::RNTupleModel::CreateBare();
   model->MakeField<RAMNTupleRecord>("record");
   model->MakeField<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField);
   model->MakeField<std::vector<std::uint8_t>>(RAMNTupleRecord::kNameBlockField);
   return model;
}

/// A worker's loop: fills every block it pops into \p ctx with
/// process(block, out, records, late_headers), which returns the block's
/// BlockOrder, and commits the block's clusters once the blocks before it are
/// committed.
template <typename Process>
void FillBlocks(BlockQueue &queue, Progress &progress, ROOT::RNTupleFillContext &ctx, Process process)
{
   BlockWriter out(ctx.CreateEntry(), ctx.CreateEntry(), [&ctx](ROOT::REntry &e) { ctx.Fill(e); });
   Block block;
   while (queue.Pop(block)) {
      size_t records = 0;
      std::vector<std::string> late_headers;
      const BlockOrder order = process(block, out, records, late_headers);
      // The records are copied into the fill context and `out` now: free the
      // input before compressing and waiting for this block's turn.
      std::vector<char>().swap(block.data);
      // Quality blocks end with the input block, whose clusters are committed as a unit.
      out.Finish();
      const std::vector<uint64_t> block_ends = out.TakeBlockEnds();
      ctx.FlushCluster();

      std::unique_lock<std::mutex> lock(progress.mutex);
      progress.next_turn.wait(lock, [&] { return progress.next_seq == block.seq || progress.failed; });
      if (progress.failed)
         return;
      ctx.CommitStagedClusters();
      progress.order.Add(order);
      for (const uint64_t end : block_ends)
         progress.block_ends.push_back(progress.records + end);
      progress.records += records;
      progress.late_headers.insert(progress.late_headers.end(), late_headers.begin(), late_headers.end());
      progress.next_seq++;
      progress.next_turn.notify_all();
   }
}

/// Starts \p threads workers, each running worker(queue, progress, ctx) on a
/// fill context of \p writer, queues the blocks next(block) reads until it
/// returns false, and returns once every block is committed. An exception in a
/// worker stops the others and comes out of here.
template <typename Worker, typename Next>
void FillInBlocks(ROOT::RNTupleParallelWriter &writer, int threads, Progress &progress, Worker worker, Next next)
{
   // One block waiting is enough: reading a block takes a fraction of the
   // time a worker spends parsing and compressing one, and every block held
   // costs its size in memory.
   BlockQueue queue(/*capacity=*/1);
   // Destruction order on an exception: close the queue, join the futures,
   // then drop the contexts; the caller drops the writer after them.
   std::vector<std::shared_ptr<ROOT::RNTupleFillContext>> contexts;
   std::vector<std::future<void>> workers;
   const CloseQueueOnExit closer(queue);
   for (int i = 0; i < threads; i++) {
      auto ctx = writer.CreateFillContext();
      ctx->EnableStagedClusterCommitting();
      contexts.push_back(ctx);
      workers.push_back(std::async(std::launch::async, [&queue, &progress, &worker, ctx] {
         const StopOthersOnException guard(progress, queue);
         worker(queue, progress, *ctx);
      }));
   }

   size_t seq = 0;
   Block block;
   while (!progress.Failed() && next(block)) {
      block.seq = seq++;
      queue.Push(std::move(block));
      block = Block{};
      queue.WaitForRoom();
   }
   queue.Close();
   for (auto &w : workers)
      w.get();
}

} // namespace ramcore::blocks

#endif // RAMCORE_BLOCKPIPELINE_H
