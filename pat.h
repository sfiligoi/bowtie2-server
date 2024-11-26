/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PAT_H_
#define PAT_H_

#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <condition_variable>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>
#include <cassert>
#include <string>
#include <ctype.h>
#include <vector>
#include <map>
#include <queue>
#include <atomic>
#include <utility>
#include "alphabet.h"
#include "assert_helpers.h"
#include "random_source.h"
#include "threading.h"
#include "qual.h"
#include "search_globals.h"
#include "sstring.h"
#include "ds.h"
#include "read.h"
#include "util.h"
#include "concurrentqueue.h"
#ifdef WITH_ZSTD
#include "zstd_decompress.h"
#endif

#ifdef USE_SRA
#include <ncbi-vdb/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>
#endif

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>

#define SET_BINARY_MODE(fd) setmode(fd, O_BINARY)
#define getc_unlocked _fgetc_nolock
#else
#define SET_BINARY_MODE(fd)
#endif

// full definition in aln_sink.h
// but we only need the pointer here
class AlnSink;
class AlnSinkSam;

/**
 * Classes and routines for reading reads from various input sources.
 */

/**
 * Parameters affecting how reads and read in.
 */
struct PatternParams {
	PatternParams() { }

	PatternParams(
		int format_,
		bool interleaved_,
		bool fileParallel_,
		uint32_t seed_,
		size_t max_buf_,
		bool solexa64_,
		bool phred64_,
		bool intQuals_,
		int trim5_,
		int trim3_,
		pair<short, size_t> trimTo_,
		int sampleLen_,
		int sampleFreq_,
		size_t skip_,
		uint64_t upto_,
		int nthreads_,
		bool fixName_,
		bool preserve_tags_,
		bool align_paired_reads_) :
		format(format_),
		interleaved(interleaved_),
		fileParallel(fileParallel_),
		seed(seed_),
		max_buf(max_buf_),
		solexa64(solexa64_),
		phred64(phred64_),
		intQuals(intQuals_),
		trim5(trim5_),
		trim3(trim3_),
		trimTo(trimTo_),
		sampleLen(sampleLen_),
		sampleFreq(sampleFreq_),
		skip(skip_),
		upto(upto_),
		nthreads(nthreads_),
		fixName(fixName_),
		preserve_tags(preserve_tags_),
		align_paired_reads(align_paired_reads_) { }

	int format;			  // file format
	bool interleaved;	  // some or all of the FASTQ/FASTA reads are interleaved
	bool fileParallel;	  // true -> wrap files with separate PatternComposers
	uint32_t seed;		  // pseudo-random seed
	size_t max_buf;		  // number of reads to buffer in one read
	bool solexa64;		  // true -> qualities are on solexa64 scale
	bool phred64;		  // true -> qualities are on phred64 scale
	bool intQuals;		  // true -> qualities are space-separated numbers
	int trim5;                // amt to hard clip from 5' end
	int trim3;                // amt to hard clip from 3' end
	pair<short, size_t> trimTo;
	int sampleLen;		  // length of sampled reads for FastaContinuous...
	int sampleFreq;		  // frequency of sampled reads for FastaContinuous...
	size_t skip;		  // skip the first 'skip' patterns
	uint64_t upto;		  // max number of queries to read
	int nthreads;		  // number of threads for locking
	bool fixName;		  //
	bool preserve_tags;       // keep existing tags when aligning BAM files
	bool align_paired_reads;
};

/**
 * All per-thread storage for input read data.
 */
struct PerThreadReadBuf {

	PerThreadReadBuf(size_t max_buf) :
		max_buf_(max_buf),
		bufa_(max_buf),
		bufb_(max_buf),
		rdid_()
	{
		bufa_.resize(max_buf);
		bufb_.resize(max_buf);
		reset();
	}

	Read& read_a() { return bufa_[cur_buf_]; }
	Read& read_b() { return bufb_[cur_buf_]; }

	const Read& read_a() const { return bufa_[cur_buf_]; }
	const Read& read_b() const { return bufb_[cur_buf_]; }

	/**
	 * Return read id for read/pair currently in the buffer.
	 */
	TReadId rdid() const {
		assert_neq(rdid_, std::numeric_limits<TReadId>::max());
		return rdid_ + cur_buf_;
	}

	/**
	 * Reset state as though no reads have been read.
	 */
	void reset() {
		cur_buf_ = bufa_.size();
		for(size_t i = 0; i < max_buf_; i++) {
			bufa_[i].reset();
			bufb_[i].reset();
		}
		rdid_ = std::numeric_limits<TReadId>::max();
	}

	/**
	 * Advance cursor to next element
	 */
	void next() {
		assert_lt(cur_buf_, bufa_.size());
		cur_buf_++;
	}

	/**
	 * Return true when there's nothing left for next().
	 */
	bool exhausted() const {
		assert_leq(cur_buf_, bufa_.size());
		return cur_buf_ >= bufa_.size()-1 || bufa_[cur_buf_+1].readOrigBuf.empty();
	}

	/**
	 * Just after a new batch has been loaded, use init to
	 * set cur_buf_ appropriately.
	 */
	void init() {
		cur_buf_ = 0;
	}

	/**
	 * Set read id of first read in buffer.
	 */
	void setReadId(TReadId rdid) {
		rdid_ = rdid;
	}

	const size_t max_buf_; // max # reads to read into buffer at once
	EList<Read> bufa_;	   // Read buffer for mate as
	EList<Read> bufb_;	   // Read buffer for mate bs
	size_t cur_buf_;	   // Read buffer currently active
	TReadId rdid_;		   // index of read at offset 0 of bufa_/bufb_
};

extern void wrongQualityFormat(const BTString& read_name);
extern void tooFewQualities(const BTString& read_name);
extern void tooManyQualities(const BTString& read_name);

/**
 * Encapsulates a synchronized source of patterns; usually a file.
 * Optionally reverses reads and quality strings before returning them,
 * though that is usually more efficiently done by the concrete
 * subclass.  Concrete subclasses should delimit critical sections with
 * calls to lock() and unlock().
 */
class PatternSource {

public:

	PatternSource(const PatternParams& p) :
		pp_(p),
		readCnt_(0),
		mutex() { }

	virtual ~PatternSource() { }

	/**
	 * Implementation to be provided by concrete subclasses.  An
	 * implementation for this member is only relevant for formats
	 * where individual input sources look like single-end-read
	 * sources, e.g., formats where paired-end reads are specified in
	 * parallel read files.
	 */
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		AlnSink* &msink,
		bool batch_a,
		bool lock = false) = 0;

	/**
	 * Finishes parsing a given read.  Happens outside the critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const = 0;

	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 */
	virtual void reset() { readCnt_ = 0; }

	/**
	 * Return a new dynamically allocated PatternSource for the given
	 * format, using the given list of strings as the filenames to read
	 * from or as the sequences themselves (i.e. if -c was used).
	 */
	static PatternSource* patsrcFromStrings(
		const PatternParams& p,
		AlnSink* msink,
		const EList<std::string>& qs);

	/**
	 * Return number of reads light-parsed by this stream so far.
	 */
	TReadId readCount() const { return readCnt_; }

protected:


	// Reference to global input-parsing parameters
	const PatternParams& pp_;

	// The number of reads read by this PatternSource
	volatile TReadId readCnt_;

	// Lock enforcing mutual exclusion for (a) file I/O, (b) writing fields
	// of this or another other shared object.
	MUTEX_T mutex;
};

/**
 * Encapsulates a source of patterns which is an in-memory vector.
 */
class VectorPatternSource : public PatternSource {

public:

	/**
	 * Populate member lists, v_, quals_, names_, etc, with information parsed
	 * from the given list of strings.
	 */
	VectorPatternSource(
		const EList<std::string>& v,
		const PatternParams& p,
		AlnSink* msink);

	virtual ~VectorPatternSource() { }

	/**
	 * Read next batch.  However, batch concept is not very applicable for this
	 * PatternSource where all the info has already been parsed into the fields
	 * in the contsructor.	This essentially modifies the pt as though we read
	 * in some number of patterns.
	 */
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		AlnSink* &msink,
		bool batch_a,
		bool lock = false);

	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 */
	virtual void reset() {
		PatternSource::reset();
		cur_ = skip_;
		paired_ = false;
	}

	/**
	 * Finishes parsing outside the critical section
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

private:

	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		AlnSink* &msink,
		bool batch_a);

	AlnSink* msink_;		// Associated sink
	size_t cur_;			   // index for first read of next batch
	size_t skip_;			   // # reads to skip
	bool paired_;			   // whether reads are paired
	EList<string> tokbuf_;	   // buffer for storing parsed tokens
	EList<Read::TBuf> bufs_;   // per-read buffers
	char nametmp_[20];		   // temp buffer for constructing name
};

/**
 * Parent class for PatternSources that read from a file.
 * Uses unlocked C I/O, on the assumption that all reading
 * from the file will take place in an otherwise-protected
 * critical section.
 */
class CFilePatternSource : public PatternSource {
	enum CompressionType {
		NONE,
		GZIP,
#ifdef WITH_ZSTD
		ZSTD,
#endif
	};
public:
	CFilePatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink) :
		PatternSource(p),
		msink_(msink),
		infiles_(infiles),
		filecur_(0),
		fp_(NULL),
		zfp_(NULL),
#ifdef WITH_ZSTD
		zstdfp_(NULL),
#endif
		is_open_(false),
		skip_(p.skip),
		first_(true),
		compressionType_(CompressionType::NONE)
	{
		assert_gt(infiles.size(), 0);
		errs_.resize(infiles_.size());
		errs_.fill(0, infiles_.size(), false);
		open(); // open first file in the list
		filecur_++;
	}

	/**
	 * Close open file.
	 */
	virtual ~CFilePatternSource() {
		if(is_open_) {
			switch (compressionType_) {
			case CompressionType::NONE:
				assert(fp_ != NULL);
				fclose(fp_);
				break;
			case CompressionType::GZIP:
				assert(zfp_ != NULL);
				gzclose(zfp_);
				break;
#ifdef WITH_ZSTD
			case CompressionType::ZSTD:
				assert(zstdfp_ != NULL);
				zstdClose(zstdfp_);
				break;
#endif
			}
		}
	}

	/**
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.	This function gets called by
	 * all the search threads, so we must handle synchronization.
	 *
	 * Returns pair<bool, int> where bool indicates whether we're
	 * completely done, and int indicates how many reads were read.
	 */
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		AlnSink* &msink,
		bool batch_a,
		bool lock = false);

	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 * Should only be called by the master thread.
	 */
	virtual void reset() {
		PatternSource::reset();
		filecur_ = 0,
		open();
		filecur_++;
	}

protected:

	/**
	 * Light-parse a batch of unpaired reads from current file into the given
	 * buffer.	Called from CFilePatternSource.nextBatch().
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx) = 0;

	/**
	 * Reset state to handle a fresh file
	 */
	virtual void resetForNextFile() { }

	/**
	 * Open the next file in the list of input files.
	 */
	void open();

	int getc_wrapper() {
		int c;

		do {
			if (compressionType_ == CompressionType::GZIP)
				c = gzgetc(zfp_);
#ifdef WITH_ZSTD
			else if (compressionType_ == CompressionType::ZSTD)
				c = zstdGetc(zstdfp_);
#endif
			else
				c = getc_unlocked(fp_);
		} while (c != EOF && c != '\t' && c != '\r' && c != '\n' && !isprint(c));

		return c;
	}

	int ungetc_wrapper(int c) {
		if (compressionType_ == CompressionType::GZIP)
			return gzungetc(c, zfp_);
#ifdef WITH_ZSTD
		else if (compressionType_ == CompressionType::ZSTD)
			return zstdUngetc(c, zstdfp_);
#endif
		else
			return ungetc(c, fp_);
	}

	int zread(voidp buf, unsigned len) {
		int r = gzread(zfp_, buf, len);
		if (r < 0) {
			const char *err = gzerror(zfp_, NULL);
			if (err != NULL) {
				std::cerr << err << std::endl;
			}
		}
		return r;
	}

	bool is_gzipped_file(int fd) {
		if (fd == -1) {
			return false;
		}

		uint8_t byte1, byte2;

		ssize_t r1 = read(fd, &byte1, sizeof(uint8_t));
		ssize_t r2 = read(fd, &byte2, sizeof(uint8_t));

		lseek(fd, 0, SEEK_SET);
                if (r1 == 0 || r2 == 0) {
                        std::cerr << "Unable to read file magic number" << std::endl;
                        return false;
                }

		if (byte1 == 0x1f && byte2 == 0x8b) {
			return true;
		}

		return false;
	}

#ifdef WITH_ZSTD
	bool is_zstd_file(int fd) {
		if (fd == -1)
			return false;

		unsigned magic;

                if (read(fd, &magic, sizeof(unsigned)) != sizeof(unsigned)) {
			std::cerr << "is_zstd_file: unable to read magic number" << std::endl;
			return false;
                }
		lseek(fd, 0, SEEK_SET);

                return magic == 0xfd2fb528;
	}
#endif

	AlnSink* msink_;		// Associated sink
	EList<std::string> infiles_;	 // filenames for read files
	EList<bool> errs_;		 // whether we've already printed an error for each file
	size_t filecur_;		 // index into infiles_ of next file to read
	FILE *fp_;			 // read file currently being read from
	gzFile zfp_;			 // compressed version of fp_
#ifdef WITH_ZSTD
	zstdStrm *zstdfp_;	         // zstd compressed file
#endif
	bool is_open_;			 // whether fp_ is currently open
	TReadId skip_;			 // number of reads to skip
	bool first_;			 // parsing first record in first file?
	char buf_[64*1024];		 // file buffer
	CompressionType compressionType_;

private:

	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		AlnSink* &msink,
		bool batch_a);
};

/**
 * Parent class for PatternSources that read from a file.
 * Uses unlocked C I/O, on the assumption that all reading
 * from a socket will take place in an otherwise-protected
 * critical section.
 */
class SocketPatternSource : public PatternSource {
public:
	SocketPatternSource(
		const PatternParams& p,
		int fd,
		long int max_bytes,
                AlnSink* msink) :
		PatternSource(p),
		read_fd_(dup(fd)),
		max_bytes_(max_bytes),
		bytes_read(0),
		buf_allocated(0),
		msink_(msink),
		fp_(NULL),
		buf_(NULL)
		//zfp_(NULL),
	{
		// TODO: Consider support for compression
		fp_ = fdopen(read_fd_, "rb");
		if (max_bytes<0) {
			// use chunked logic, pre-allocate small-ish buffer
			buf_=(char*) malloc(30000);
			buf_allocated = 30000;
			max_bytes_ = 0; // declare a new buffer need to be read
		}
	}

	/**
	 * Close open file.
	 */
	virtual ~SocketPatternSource() {
		if (buf_!=NULL) {
			free(buf_);
			buf_ = NULL;
		}
		if (fp_!=NULL) {
			fclose(fp_);
			fp_ = NULL;
		}
	}

	/**
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.	This function gets called by
	 * all the search threads, so we must handle synchronization.
	 *
	 * Returns pair<bool, int> where bool indicates whether we're
	 * completely done, and int indicates how many reads were read.
	 */
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		AlnSink* &msink,
		bool batch_a,
		bool lock = false);

protected:

	/**
	 * Light-parse a batch of unpaired reads from current file into the given
	 * buffer.	Called from CFilePatternSource.nextBatch().
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx) = 0;

	// get next element from the buffer
	// assume whole line is in the buffer (no carry overs)
	int getc_buffer_wrapper() {
		int c;

		do {
			if ((max_bytes_>=0) && (bytes_read>=max_bytes_)) return EOF; // already read everythign there was to read
			c = buf_[bytes_read];
			if (c!=EOF) bytes_read++;
		} while (c != EOF && c != '\t' && c != '\r' && c != '\n' && !isprint(c));

		return c;
	}

	// read straight from fp
	int getc_fp_wrapper() {
		int c;

		do {
			if ((max_bytes_>=0) && (bytes_read>=max_bytes_)) return EOF; // already read everythign there was to read
			c = getc_unlocked(fp_);
			if (c!=EOF) bytes_read++;
		} while (c != EOF && c != '\t' && c != '\r' && c != '\n' && !isprint(c));

		return c;
	}

	// assuming chunked transfer, read the hex header len
	// Return -1 on error
	long int read_buf_len() {
		char head[8];
		int hlen=0;
		int c = getc_unlocked(fp_);
		if (!isxdigit(c)) return -1; // the first char must be a hex number
		head[hlen] = c;
		hlen++;
		while (hlen<8) {
			c = getc_unlocked(fp_);
			if (c=='\r') continue; // ignore
			if (c=='\n') {
				// found separator, return the content of the buffer
				head[hlen] = 0; // make it a proper string
				return strtol(head,NULL,16);
			}
			if (!isxdigit(c)) return -1; // else, only hex numbers allowed
			head[hlen] = c;
			hlen++;
		}
		// fprintf(stderr,"WARN: Client sent too big of a chunk\n");
		return -1; // if we got here, the counter was too high
	}

	void next_buffer_chunk() {
		bytes_read = 0;
		// read the length
		long int len = read_buf_len();
		//fprintf(stderr, "Buf len: %i\n",int(len));
		// avoid very large buffers, too
		if ((len<0) || (len>999999)) { // len==0 means end of data, so treat same as error (TODO: improve)
			if (len>=0) fprintf(stderr,"WARN: Chunk too large, aborting request!\n");
			// we know the buffer is at least one byte long
			buf_[0] = EOF;
			max_bytes_ = 1;
			return;
		}
		// make sure we have enough space
		if (len>buf_allocated) {
			buf_ = (char*) realloc(buf_, len);
			buf_allocated = len;
		}
		// read in the buffer content
		for (long int i=0; i<len; i++) {
			int c = getc_unlocked(fp_);
			if (c==EOF) {
				fprintf(stderr,"WARN: Unexpected EOF found!\n");
				// something went wrong, throw away the whole buffer (TODO: improve)
				buf_[0] = EOF;
				max_bytes_ = 1;
				return;
			}
			buf_[i] = c;
		}
		// read and discard the final separator
		{
			int c = getc_unlocked(fp_);
			if (c=='\r') c = getc_unlocked(fp_); // \r is optional
			if (c!='\n') {
				fprintf(stderr,"WARN: Invalid chunked termnation, aborting request\n");
				// something went wrong, throw away the whole buffer (TODO: improve)
				buf_[0] = EOF;
				max_bytes_ = 1;
				return;
			}
		}

		if (len>0) {
			max_bytes_ = len;
		} else {
			// len==0 is special, it means this is the last chunk
			buf_[0] = EOF;
			max_bytes_ = 1;
		}
	}

	int getc_wrapper() {
		if (buf_==NULL) return getc_fp_wrapper(); // not buffered, just read straight from fp
		if (bytes_read<max_bytes_) return getc_buffer_wrapper(); // we have not finished reading from the existing buf, so just use that

		next_buffer_chunk();
		return getc_buffer_wrapper();
	}

	const int read_fd_;		// socker fd
	long int max_bytes_;		// max number of bytes to read from read_fd or buf_
	long int bytes_read;		// how many bytes did we read so far
	long int buf_allocated;		// how many bytes does the buf contain
	AlnSink* msink_;		// Associated sink
	FILE *fp_;			// read file currently being read from
	char *buf_;			// If not NULL, read from the buffer
	//gzFile zfp_;			 // compressed version of fp_
	//CompressionType compressionType_;

private:

	pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		AlnSink* &msink,
		bool batch_a);
};

/**
 * Synchronized concrete pattern source for a list of FASTA files.
 */
class FastaPatternSource : public CFilePatternSource {

public:

	FastaPatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink,
		bool interleaved) :
		CFilePatternSource(infiles, p, msink),
		first_(true),
		interleaved_(interleaved) { }

	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 * Should only be called by the master thread.
	 */
	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize FASTA parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a FASTA batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx);

	/**
	 * Scan to the next FASTA record (starting with >) and return the first
	 * character of the record (which will always be >).
	 */
	static int skipToNextFastaRecord(FileBuf& in) {
		int c;
		while((c = in.get()) != '>') {
			if(in.eof()) return -1;
		}
		return c;
	}

	/**
	 * Reset state to handle a fresh file
	 */
	virtual void resetForNextFile() {
		first_ = true;
	}

	bool first_;
	bool interleaved_;
};

/**
 * Synchronized concrete pattern source for a list of files with tab-
 * delimited name, seq, qual fields (or, for paired-end reads,
 * basename, seq1, qual1, seq2, qual2).
 */
class TabbedPatternSource : public CFilePatternSource {

public:

	TabbedPatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink,
		bool  secondName) :  // whether it's --12/--tab5 or --tab6
		CFilePatternSource(infiles, p, msink),
		secondName_(secondName) { }

	/**
	 * Finalize tabbed parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a batch of tabbed-format reads into given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx);

	bool secondName_;	// true if --tab6, false if --tab5
};

/**
 * Synchronized concrete pattern source for a list of files with tab-
 * delimited name, seq, qual fields (or, for paired-end reads,
 * basename, seq1, qual1, seq2, qual2).
 */
class TabbedSocketPatternSource : public SocketPatternSource {

public:

	TabbedSocketPatternSource(
		const PatternParams& p,
		int fd, 
		long int max_bytes,
                AlnSink* msink) :
		SocketPatternSource(p, fd, max_bytes, msink) {}

	/**
	 * Finalize tabbed parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a batch of tabbed-format reads into given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx);
};

/**
 * Synchronized concrete pattern source for Illumina Qseq files.  In
 * Qseq files, each read appears on a separate line and the tab-
 * delimited fields are:
 *
 * 1. Machine name
 * 2. Run number
 * 3. Lane number
 * 4. Tile number
 * 5. X coordinate of spot
 * 6. Y coordinate of spot
 * 7. Index: "Index sequence or 0. For no indexing, or for a file that
 *	  has not been demultiplexed yet, this field should have a value of
 *	  0."
 * 8. Read number: 1 for unpaired, 1 or 2 for paired
 * 9. Sequence
 * 10. Quality
 * 11. Filter: 1 = passed, 0 = didn't
 */
class QseqPatternSource : public CFilePatternSource {

public:

	QseqPatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink):
		CFilePatternSource(infiles, p, msink) { }

	/**
	 * Finalize qseq parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx);

	EList<std::string> qualToks_;
};

/**
 * Synchronized concrete pattern source for a list of FASTA files where
 * reads need to be extracted from long continuous sequences.
 */
class FastaContinuousPatternSource : public CFilePatternSource {
public:
	FastaContinuousPatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink):
		CFilePatternSource(infiles, p, msink),
		length_(p.sampleLen),
		freq_(p.sampleFreq),
		eat_(length_-1),
		beginning_(true),
		bufCur_(0),
		cur_(0llu),
		last_(0llu)
	{
		assert_gt(freq_, 0);
		resetForNextFile();
		assert_leq(length_, 1024);
	}

	virtual void reset() {
		CFilePatternSource::reset();
		resetForNextFile();
	}

	/**
	 * Finalize FASTA parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx);

	/**
	 * Reset state to be read for the next file.
	 */
	virtual void resetForNextFile() {
		eat_ = length_-1;
		name_prefix_buf_.clear();
		beginning_ = true;
		bufCur_ = 0;
		last_ = cur_;
	}

private:
	const size_t length_; /// length of reads to generate
	const size_t freq_;   /// frequency to sample reads
	size_t eat_;		/// number of characters we need to skip before
						/// we have flushed all of the ambiguous or
						/// non-existent characters out of our read
						/// window
	bool beginning_;	/// skipping over the first read length?
	char buf_[1024];	/// FASTA sequence buffer
	Read::TBuf name_prefix_buf_; /// FASTA sequence name buffer
	char name_int_buf_[20]; /// for composing offsets for names
	size_t bufCur_;		/// buffer cursor; points to where we should
						/// insert the next character
	uint64_t cur_;
	uint64_t last_;     /// number to subtract from readCnt_ to get
						/// the pat id to output (so it resets to 0 for
						/// each new sequence)
};

/**
 * Read a FASTQ-format file.
 * See: http://maq.sourceforge.net/fastq.shtml
 */
class FastqPatternSource : public CFilePatternSource {

public:

	FastqPatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink,
		bool interleaved) :
		CFilePatternSource(infiles, p, msink),
		first_(true),
		interleaved_(interleaved) { }

	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize FASTQ parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx);

	/**
	 * Reset state to be ready for the next file.
	 */
	virtual void resetForNextFile() {
		first_ = true;
	}

	bool first_;		// parsing first read in file
	bool interleaved_;	// fastq reads are interleaved
};

class BAMPatternSource : public CFilePatternSource {
	struct hdr_t {
		uint8_t id1;
		uint8_t id2;
		uint8_t cm;
		uint8_t flg;
		uint32_t mtime;
		uint8_t xfl;
		uint8_t os;
		uint16_t xlen;
	};

	struct ftr_t {
		uint32_t crc32;
		uint32_t isize;
	};

	struct BGZF {
		hdr_t hdr;
		uint8_t cdata[1 << 16];
		ftr_t ftr;
	};

	struct BAMField {
		enum aln_rec_field_name {
			refID,
			pos,
			l_read_name,
			mapq,
			bin,
			n_cigar_op,
			flag,
			l_seq,
			next_refID,
			next_pos,
			tlen,
			read_name,
		};
	};

public:

	BAMPatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink) : 
		CFilePatternSource(infiles, p, msink),
		first_(true),
		alignment_batch(0),
		alignment_offset(0),
		delta_(0),
		pp_(p)
		{
			stream.zalloc = Z_NULL;
			stream.zfree = Z_NULL;
			stream.opaque = Z_NULL;
			alignment_batch.reserve(1 << 16);
		}

	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize BAM parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

	~BAMPatternSource() {
	}


protected:

	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf& pt, AlnSink* &msink, bool batch_a, bool lock = false);

	uint16_t nextBGZFBlockFromFile(BGZF& block);

	/**
	 * Reset state to be ready for the next file.
	 */
	virtual void resetForNextFile() {
		first_ = true;
	}

	bool first_; // parsing first read in file

private:
	virtual std::pair<bool, int> nextBatchFromFile(PerThreadReadBuf&, bool, unsigned) {
		return make_pair(true, 0);
	}

	int decompress_bgzf_block(uint8_t *dst, size_t dst_len, uint8_t *src, size_t src_len);
	std::pair<bool, int> get_alignments(PerThreadReadBuf& pt, bool batch_a, unsigned& readi, bool lock);

        static const int offset[];
	static const uint8_t EOF_MARKER[];

	std::vector<uint8_t> alignment_batch;
	size_t alignment_offset;
	size_t delta_;
	z_stream stream;

	PatternParams pp_;
};

/**
 * Read a Raw-format file (one sequence per line).	No quality strings
 * allowed.  All qualities are assumed to be 'I' (40 on the Phred-33
 * scale).
 */
class RawPatternSource : public CFilePatternSource {

public:

	RawPatternSource(
		const EList<std::string>& infiles,
		const PatternParams& p,
                AlnSink* msink) : 
		CFilePatternSource(infiles, p, msink), first_(true) { }

	virtual void reset() {
		first_ = true;
		CFilePatternSource::reset();
	}

	/**
	 * Finalize raw parsing outside critical section.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

protected:

	/**
	 * Light-parse a batch into the given buffer.
	 */
	virtual std::pair<bool, int> nextBatchFromFile(
		PerThreadReadBuf& pt,
		bool batch_a,
		unsigned read_idx);

	/**
	 * Reset state to be ready for the next file.
	 */
	virtual void resetForNextFile() {
		first_ = true;
	}

private:

	bool first_;
};

/**
 * Abstract parent class for synhconized sources of paired-end reads
 * (and possibly also single-end reads).
 */
class PatternComposer {
public:
	PatternComposer() : mutex_m() { }

	virtual ~PatternComposer() { }

	virtual void reset() = 0;

	/**
	 * Member function override by concrete, format-specific classes.
	 */
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf& pt, AlnSink* &msink) = 0;

	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) = 0;

	/**
	 * Given the values for all of the various arguments used to specify
	 * the read and quality input, create a list of pattern sources to
	 * dispense them.
	 */
	static PatternComposer* setupPatternComposer(
		const EList<std::string>& si,	 // singles, from argv
		const EList<std::string>& m1,	 // mate1's, from -1 arg
		const EList<std::string>& m2,	 // mate2's, from -2 arg
		const EList<std::string>& m12,	 // both mates on each line, from --12
		const EList<std::string>& q,	 // qualities associated with singles
		const EList<std::string>& q1,	 // qualities associated with m1
		const EList<std::string>& q2,	 // qualities associated with m2
#ifdef USE_SRA
		const EList<string>& sra_accs,   // SRA accessions
#endif
		PatternParams& p,		// read-in params
		AlnSink* msink,
		bool verbose);				// be talkative?

protected:

	static void free_EList_pmembers(const EList<PatternSource*>&);

	/// Lock enforcing mutual exclusion for (a) file I/O, (b) writing fields
	/// of this or another other shared object.
	MUTEX_T mutex_m;
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class SoloPatternComposer : public PatternComposer {

public:

	SoloPatternComposer(
		const EList<PatternSource*>* src):
		PatternComposer(),
		cur_(0),
		src_(src)
	{
		assert(src_ != NULL);
		lock_ = false;
		for(size_t i = 0; i < src_->size(); i++) {
			assert((*src_)[i] != NULL);
		}
	}

	virtual ~SoloPatternComposer() {
		free_EList_pmembers(*src_);
		delete src_;
	}

	/**
	 * Reset this object and all the PatternSources under it so that
	 * the next call to nextBatch gets the very first read.
	 */
	virtual void reset() {
		for(size_t i = 0; i < src_->size(); i++) {
			(*src_)[i]->reset();
		}
		cur_ = 0;
	}

	/**
	 * Calls member functions of the individual PatternSource objects
	 * to get more reads.  Since there's no need to keep two separate
	 * files in sync (as there is for DualPatternComposer),
	 * synchronization can be handed by the PatternSource contained
	 * in the src_ field.
	 */
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf& pt, AlnSink* &msink);

	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) {
		return (*src_)[0]->parse(ra, rb, rdid);
	}

protected:
	volatile bool lock_;
	volatile size_t cur_; // current element in parallel srca_, srcb_ vectors
	const EList<PatternSource*>* src_; /// PatternSources for paired-end reads
};

/**
 * Encapsulates a synchronized source of both paired-end reads and
 * unpaired reads, where the paired-end must come from parallel files.
 */
class DualPatternComposer : public PatternComposer {

public:

	DualPatternComposer(
		const EList<PatternSource*>* srca,
		const EList<PatternSource*>* srcb):
		PatternComposer(), cur_(0), srca_(srca), srcb_(srcb)
	{
		assert(srca_ != NULL);
		assert(srcb_ != NULL);
		// srca_ and srcb_ must be parallel
		assert_eq(srca_->size(), srcb_->size());
		lock_ = false;
		// lock_ = false;
		for(size_t i = 0; i < srca_->size(); i++) {
			// Can't have NULL first-mate sources.	Second-mate sources
			// can be NULL, in the case when the corresponding first-
			// mate source is unpaired.
			assert((*srca_)[i] != NULL);
			for(size_t j = 0; j < srcb_->size(); j++) {
				assert_neq((*srca_)[i], (*srcb_)[j]);
			}
		}
	}

	virtual ~DualPatternComposer() {
		free_EList_pmembers(*srca_);
		delete srca_;
		free_EList_pmembers(*srcb_);
		delete srcb_;
	}

	/**
	 * Reset this object and all the PatternSources under it so that
	 * the next call to nextBatch gets the very first read.
	 */
	virtual void reset() {
		for(size_t i = 0; i < srca_->size(); i++) {
			(*srca_)[i]->reset();
			if((*srcb_)[i] != NULL) {
				(*srcb_)[i]->reset();
			}
		}
		cur_ = 0;
	}

	/**
	 * Calls member functions of the individual PatternSource objects
	 * to get more reads.  Since we need to keep the two separate
	 * files in sync, synchronization can be handed at this level, with
	 * one critical section surrounding both calls into the srca_ and
	 * srcb_ member functions.
	 */
	virtual std::pair<bool, int> nextBatch(PerThreadReadBuf& pt, AlnSink* &msink);

	/**
	 * Make appropriate call into the format layer to parse individual read.
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) {
		return (*srca_)[0]->parse(ra, rb, rdid);
	}

protected:

	volatile bool lock_;
	volatile size_t cur_; // current element in parallel srca_, srcb_ vectors
	const EList<PatternSource*>* srca_; // for 1st matesunpaired
	const EList<PatternSource*>* srcb_; // for 2nd mates
};

/**
 * Encapsulates a single thread's interaction with the PatternSource.
 * Most notably, this class holds the buffers into which the
 * PatterSource will write sequences.  This class is *not* threadsafe
 * - it doesn't need to be since there's one per thread.  PatternSource
 * is thread-safe.
 */
class PatternSourcePerThread {

public:

	PatternSourcePerThread(
		PatternComposer& composer,
		const PatternParams& pp) :
		composer_(composer),
		buf_(pp.max_buf),
		pp_(pp),
		msink_(NULL),
		last_batch_(false),
		last_batch_size_(0) { }

	// If it returns true, nextReadPair is non-blocking and should return success
	bool nextReadPairReady() const {return !buf_.exhausted();}

	/**
	 * Use objects in the PatternSource and/or PatternComposer
	 * hierarchies to populate the per-thread buffers.
	 */
	std::pair<bool, bool> nextReadPair();

	AlnSink& msink() {return *msink_; }

	Read& read_a() { return buf_.read_a(); }
	Read& read_b() { return buf_.read_b(); }

	const Read& read_a() const { return buf_.read_a(); }
	const Read& read_b() const { return buf_.read_b(); }

private:

	/**
	 * When we've finished fully parsing and dishing out reads in
	 * the current batch, we go get the next one by calling into
	 * the composition layer.
	 */
	std::pair<bool, int> nextBatch() {
		buf_.reset();
		msink_ = NULL;
		std::pair<bool, int> res = composer_.nextBatch(buf_, msink_);
		buf_.init();
		return res;
	}

	/**
	 * Once name/sequence/qualities have been parsed for an
	 * unpaired read, set all the other key fields of the Read
	 * struct.
	 */
	void finalize(Read& ra);

	/**
	 * Once name/sequence/qualities have been parsed for a
	 * paired-end read, set all the other key fields of the Read
	 * structs.
	 */
	void finalizePair(Read& ra, Read& rb);

	/**
	 * Call into composition layer (which in turn calls into
	 * format layer) to parse the read.
	 */
	bool parse(Read& ra, Read& rb) {
		return composer_.parse(ra, rb, buf_.rdid());
	}

	void trim(Read& r) {
		if (pp_.trimTo.second > 0) {
			switch (pp_.trimTo.first) {
				case 3:
					if (r.patFw.length() > pp_.trimTo.second) {
						r.trimmed5 = r.patFw.length() - pp_.trimTo.second;
						r.patFw.trimEnd(r.trimmed5);
						r.qual.trimEnd(r.trimmed5);
					}
					break;
				case 5:
					if (r.patFw.length() > pp_.trimTo.second) {
						r.trimmed3 = r.patFw.length() - pp_.trimTo.second;
						r.patFw.trimBegin(r.trimmed3);
						r.qual.trimBegin(r.trimmed3);
					}
					break;
			}
		}
	}

	PatternComposer& composer_; // pattern composer
	PerThreadReadBuf buf_;		// read data buffer
	const PatternParams& pp_;	// pattern-related parameters
	AlnSink        *msink_;         // associated sink
	bool last_batch_;			// true if this is final batch
	int last_batch_size_;		// # reads read in previous batch
};

/**
 * Abstract parent factory for PatternSourcePerThreads.
 */
class PatternSourcePerThreadFactory {
public:
	PatternSourcePerThreadFactory(
		PatternComposer& composer,
		const PatternParams& pp) :
		composer_(composer),
		pp_(pp) {}

	/**
	 * Create a new heap-allocated PatternSourcePerThreads.
	 */
	virtual PatternSourcePerThread* create() const {
		return new PatternSourcePerThread(composer_, pp_);
	}

	/**
	 * Create a new heap-allocated vector of heap-allocated
	 * PatternSourcePerThreads.
	 */
	virtual EList<PatternSourcePerThread*>* create(uint32_t n) const {
		EList<PatternSourcePerThread*>* v = new EList<PatternSourcePerThread*>;
		for(size_t i = 0; i < n; i++) {
			v->push_back(new PatternSourcePerThread(composer_, pp_));
			assert(v->back() != NULL);
		}
		return v;
	}

	virtual ~PatternSourcePerThreadFactory() {}

private:
	/// Container for obtaining paired reads from PatternSources
	PatternComposer& composer_;
	const PatternParams& pp_;
};

using namespace moodycamel;

class PatternSourceReadAheadFactory {
public:
	class ReadElement {
	public:
		PatternSourcePerThread *ps;  // not owned, just a pointer
		std::pair<bool, bool>   readResult;
	};

	// Simple wrapper for safely holding the result of PatternSourceReadAheadFactory
	class ReadAhead {
	public:
		ReadAhead(PatternSourceReadAheadFactory& fact) :
			fact_(fact),
			re_(fact.nextReadPair()),
			first_(true)	{}

		~ReadAhead() {
			fact_.returnUnready(re_);
		}

		const std::pair<bool, bool>& nextReadResult() {
			if (first_) {
				// nextReadPair was already called in the psrah constructor
				first_ = false;
			} else {
				re_.readResult = re_.ps->nextReadPair();
			}
			return re_.readResult;
		}

		PatternSourcePerThread* ptr() {return re_.ps;}
	private:
		PatternSourceReadAheadFactory& fact_;
		PatternSourceReadAheadFactory::ReadElement re_;
		bool first_;
	};

	PatternSourceReadAheadFactory(
		PatternComposer& composer,
		const PatternParams& pp, size_t n, bool useCV) :
		psfact_(composer,pp),
		psq_ready_(useCV),
		psq_idle_(useCV,psfact_,n),
		asynct_(readAsync, this) {}

	~PatternSourceReadAheadFactory() {
		returnUnready(NULL); // this will signal asynct_ it is time to quit
		asynct_.join();
	}

	// wait for data, if none in the queue
	ReadElement nextReadPair() {
		return psq_ready_.pop();
	}

	void returnUnready(PatternSourcePerThread* ps) {
		psq_idle_.push(ps);
	}

	void returnUnready(ReadElement& re) {
		returnUnready(re.ps);
	}

private:

	template <typename T>
	class LockedQueueCV {
	public:
		virtual ~LockedQueueCV() {}

		bool empty() {
			bool ret = false;
			{
				std::unique_lock<std::mutex> lk(m_);
				ret = q_.empty();
			}
			return ret;
		}

		void push(T& ps) {
			std::unique_lock<std::mutex> lk(m_);
			q_.push(ps);
			cv_.notify_all();
		}

		// wait for data, if none in the queue
		T pop() {
			T ret;
			{
				std::unique_lock<std::mutex> lk(m_);
				cv_.wait(lk, [this] { return !q_.empty();});
				ret = q_.front();
				q_.pop();
			}
			return ret;
		}

	protected:
		std::mutex m_;
		std::condition_variable cv_;
		std::queue<T> q_;
	};

	class LockedPSQueueCV : public LockedQueueCV<PatternSourcePerThread*> {
	public:
		LockedPSQueueCV(PatternSourcePerThreadFactory& psfact, size_t n) : 
			LockedQueueCV<PatternSourcePerThread*>() {
			for (size_t i=0; i<n; i++) {
				q_.push(psfact.create());
			}
		}

		virtual ~LockedPSQueueCV() {
			while (!q_.empty()) {
				delete q_.front();
				q_.pop();
			}
		}
	};

	class LockedREQueueCV : public LockedQueueCV<ReadElement> {
	public:
		virtual ~LockedREQueueCV() {
			// we actually own ps while in the queue
			while (!q_.empty()) {
				delete q_.front().ps;
				q_.pop();
			}
		}
	};

	template <typename T>
	class LockedQueueMC {
	public:
		virtual ~LockedQueueMC() {}

		bool empty() {
			return q_.size_approx() == 0;
		}

		void push(T& ps) {
			q_.enqueue(ps);
		}

		// wait for data, if none in the queue
		T pop() {
			T ret;

			while (!q_.try_dequeue(ret)) ;
			return ret;
		}

	protected:
		ConcurrentQueue<T> q_;
	};

	class LockedPSQueueMC: public LockedQueueMC<PatternSourcePerThread*> {
	public:
		LockedPSQueueMC(PatternSourcePerThreadFactory& psfact, size_t n) : 
			LockedQueueMC<PatternSourcePerThread*>() {
			for (size_t i=0; i<n; i++) {
				q_.enqueue(psfact.create());
			}
		}

		virtual ~LockedPSQueueMC() {
			PatternSourcePerThread *item;

			while (q_.size_approx() != 0) {
				while (!q_.try_dequeue(item)) ;
				delete item;
			}
		}

	};

	class LockedREQueueMC : public LockedQueueMC<ReadElement> {
	public:
		virtual ~LockedREQueueMC() {
			// we actually own ps while in the queue
			ReadElement item;

			while (q_.size_approx() != 0) {
				while (!q_.try_dequeue(item)) ;
				delete item.ps;
			}

		}
	};

	template <typename T, typename Q1, typename Q2>
	class LockedQueueDuo {
	public:
		LockedQueueDuo(Q1 *pq1, Q2 *pq2) 
		: pq1_(pq1), pq2_(pq2)
		{
			assert( (pq1_!=NULL) || (pq2_!=NULL) );
		}

		virtual ~LockedQueueDuo() {
			if (pq1_!=NULL) delete pq1_;
			if (pq2_!=NULL) delete pq2_;
		}

		bool empty() {
			return (pq1_!=NULL) ? pq1_->empty() : pq2_->empty();
		}

		void push(T& ps) {
			if (pq1_!=NULL) {
				pq1_->push(ps);
			} else {
				pq2_->push(ps);
			}
		}

		T pop() {
			return (pq1_!=NULL) ? pq1_->pop() : pq2_->pop();
		}

	protected:
		// owned, and exactly one should be not NULL
		Q1 *pq1_;
		Q2 *pq2_;
	};

	class LockedPSQueue: public LockedQueueDuo<PatternSourcePerThread*, LockedPSQueueCV, LockedPSQueueMC> {
	public:
		LockedPSQueue(bool useCV, PatternSourcePerThreadFactory& psfact, size_t n) : 
			LockedQueueDuo<PatternSourcePerThread*, LockedPSQueueCV, LockedPSQueueMC>(
				useCV ? new LockedPSQueueCV(psfact,n) : NULL,
				useCV ? NULL : new LockedPSQueueMC(psfact,n)
			) {}
	};

	class LockedREQueue: public LockedQueueDuo<ReadElement, LockedREQueueCV, LockedREQueueMC> {
	public:
		LockedREQueue(bool useCV) : 
			LockedQueueDuo<ReadElement, LockedREQueueCV, LockedREQueueMC>(
				useCV ? new LockedREQueueCV() : NULL,
				useCV ? NULL : new LockedREQueueMC()
			) {}
	};

        static void readAsync(PatternSourceReadAheadFactory *obj) {
		LockedREQueue &psq_ready = obj->psq_ready_;
		LockedPSQueue &psq_idle = obj->psq_idle_;
                while(true) {
			ReadElement re;
			re.ps = psq_idle.pop();
			if (re.ps==NULL) break; // the destructor added this in the queue

			if (re.ps->nextReadPairReady()) {
				// Should never get in here, but just in case
				re.readResult = make_pair(true, false);
			} else {
				re.readResult = re.ps->nextReadPair();
			}
			psq_ready.push(re);
                }
	}

	PatternSourcePerThreadFactory psfact_;
	LockedREQueue psq_ready_;
	LockedPSQueue psq_idle_;
	std::thread asynct_;
};

class PatternSourceServiceFactory {
private:
	class LockedIdleQueueCV;
public:
	class ReadElement {
	public:
		// Nothing read yet
		ReadElement(PatternSourcePerThread *ps_, PatternSourceServiceFactory::LockedIdleQueueCV&   psq_idle) :
			ps(ps_),
			readResult(make_pair(false, false)),
			psq_idle_(psq_idle),
			lastBatch_(false)	{}

		// final, invalid one
		ReadElement(PatternSourceServiceFactory::LockedIdleQueueCV&   psq_idle) :
			ps(NULL),
	       		psq_idle_(psq_idle)	{}

		// we allow the object to be moved around
		ReadElement(const ReadElement& other) = default;
		ReadElement(ReadElement&& other) = default;

		void nextReadPair() {
			readResult = ps->nextReadPair();
			handleLast();
		}

		// Service-internal last flag
		bool isLast() const { return lastBatch_;}

		void returnUnready() {
			if (lastBatch_) {
				// Nothing more to do with this one, cleanup
				if(ps!=NULL) delete ps;
				// inject NULL RE for proper handling in accept
				ps = NULL;
			}
			psq_idle_.push(*this);
		}

		// not owned owned by the object, just pointer
		// cleanup may still be needed on some actions
		PatternSourcePerThread *ps;
		std::pair<bool, bool>   readResult;
	private:
		void handleLast() {
			bool lastBatch = readResult.second;
			if(lastBatch) {
				// Do not expose it to the main loop, it is internal info
				readResult.second = false;
				lastBatch_ = true;
			}
		}

		PatternSourceServiceFactory::LockedIdleQueueCV&   psq_idle_;
		bool lastBatch_;
	};

	// Simple wrapper for safely holding the result of PatternSourceServiceFactory
	class ReadAhead {
	public:
		ReadAhead(PatternSourceServiceFactory& fact) :
			fact_(fact),
			re_(fact.nextReadPair()),
			first_(true)	{}
	
		~ReadAhead() {
			fact_.returnUnready(re_);
		}

		const std::pair<bool, bool>& nextReadResult() {
			if (first_) {
				// nextReadPair was already called in the psrah constructor
				first_ = false;
			} else {
				re_.nextReadPair();
			}
			return re_.readResult;
		}

		PatternSourcePerThread* ptr() {return re_.ps;}
	private:
		PatternSourceServiceFactory& fact_;
		PatternSourceServiceFactory::ReadElement re_;
		bool first_;
	};

	class Config {
	public:
		Config(const char *iname) : 
			index_name(iname),
			seedLen(0),
			maxDpStreak(0),
			seedRounds(0),
			khits(-1) {}

		const char *index_name;
		int seedLen;
		int maxDpStreak;
		int seedRounds;
		int khits; // set to -1 if allHits==true
	};

	PatternSourceServiceFactory(
		int server_port,
		PatternComposer& composer,
		const PatternParams& pp, size_t n_readahead,
		AlnSinkSam &msink,
		const Config& config):
		server_port_(server_port),
		server_backlog_(128),
		config_(config),
		pp_(pp),
		template_msink_(msink),
		base_url_(string("/BT2SRV/")+config.index_name),
		psfact_(composer,pp),
		n_readahead_(n_readahead),
		psq_ready_(),
		listent_(acceptConnections, this) {}

	~PatternSourceServiceFactory() {
		// TODO: Signal listent_ it is time to quit
		listent_.join();
#if 0
		// TODO: proper cleanup
		ReadElement re(NULL); // this will signal readaheadt_ it is time to quit
		returnUnready(re);
#endif
	}

	// wait for data, if none in the queue
	ReadElement nextReadPair() {
		return psq_ready_.pop();
	}

	void returnUnready(ReadElement& re) {
		re.returnUnready();
	}

private:

	template <typename T>
	class LockedQueueCV {
	public:
		virtual ~LockedQueueCV() {}

		bool empty() {
			bool ret = false;
			{
				std::unique_lock<std::mutex> lk(m_);
				ret = q_.empty();
			}
			return ret;
		}

		void push(T& ps) {
			std::unique_lock<std::mutex> lk(m_);
			q_.push(ps);
			cv_.notify_all();
		}

		// wait for data, if none in the queue
		T pop() {
			std::unique_lock<std::mutex> lk(m_);
			cv_.wait(lk, [this] { return !q_.empty();});
			T ret(q_.front());
			q_.pop();
			return ret;
		}

	protected:
		std::mutex m_;
		std::condition_variable cv_;
		std::queue<T> q_;
	};

	class LockedIdleQueueCV : public LockedQueueCV<ReadElement> {
	public:
		virtual ~LockedIdleQueueCV() {
			while (!q_.empty()) {
				delete q_.front().ps;
				q_.pop();
			}
		}
	};

	class LockedReadyQueueCV : public LockedQueueCV<ReadElement> {
	public:
		virtual ~LockedReadyQueueCV() {
			// we actually own ps while in the queue
			while (!q_.empty()) {
				delete q_.front().ps;
				q_.pop();
			}
		}
	};

	// write a string with retries, but silently abort on error
	static void try_write_str(int fd, const char *str) {
		int len = strlen(str);
		while(len>0) {
			int written = ::write(fd,str,len);
			if (written<1) {
				// try once more only
				written = ::write(fd,str,len);
			}

			if (written>0) {
				// get ready for next chunk, if needed
				len-= written;
				str+=written;
			} else {
				break; // just abort
			}
		}
	}

	// write a string with retries, report any errors
	// returns false in case of error
	static bool write_str(int fd, const char *str, const int str_len) {
		int len = str_len;
		while(len>0) {
			int written = ::write(fd,str,len);
			if (written<1) {
				// try once more only
				written = ::write(fd,str,len);
			}

			if (written>0) {
				// get ready for next chunk, if needed
				len-= written;
				str+=written;
			} else {
				return false;
			}
		}
		return true;
	}
	static bool write_str(int fd, const char *str) {
		int len = strlen(str);
		return write_str(fd, str, len);
	}

	// read until \n\n detected
	// can NOT go over \n\n
	// return true if we read right up to t\n\n
        static bool read_header(int fd, char *buf, int& buf_len);

	bool is_legit_align_header(char buf[], int nels);
	bool is_legit_config_header(char buf[], int nels);

	// extract content length from the header string
	// return -1 if cannot find it
	static long int find_content_length(const char str[]);

	// look for the terminator request
	// returns true if it finds one
	static bool find_request_terminator(const char str[]);

	// look for Transfer-Encoding: chunked
	static bool find_chunked_encoding(const char str[]);

	// just return the config
	// if is_header==true, prepend with X-BT2SRV-
	bool reply_config(int fd, bool is_header);

	// this is the real alignment happens
	// We only paritally parsed the header
	// buf contains what we read from fd so far
	// returns false in case of error
	bool align(int fd, long int data_size);

	// read header and pick the right response
        static void serveConnection(PatternSourceServiceFactory *obj, int client_fd);

	// called by the main (listening) thread
	inline void add_client(int client_fd) {
		std::unique_lock<std::mutex> lk(m_); // must be locked since we are updating a shared resource (clients)
		std::thread *client_t = new std::thread(serveConnection,this, client_fd);
		clients_[client_fd] = client_t;
	}

	// should be called by the main thread
	inline void maintain_clients() {
		std::unique_lock<std::mutex> lk(m_); // must be locked since we are updating a shared resource (finalizing)
		for (std::thread* &client_t : finalizing_) {
			client_t->join();
			delete client_t;
			client_t = NULL;
		}
		finalizing_.clear();
	}

	// should be called by the main thread
	inline void final_wait_clients() {
		while (!clients_.empty()) sleep(1); // we need to wait for all the clients to finish, but no need to be efficient in waiting
		maintain_clients(); // final cleanup
	}

	// called by each client thread
	void finalize_client(int client_fd) {
		std::unique_lock<std::mutex> lk(m_); // must be locked since we are updating a shared resource (clients,finalizing)
		auto client_itr = clients_.find(client_fd);
		if (client_itr!=clients_.end()) { // we should always get in here, but just in case
			std::thread *client_t = client_itr->second;
			clients_.erase(client_itr);
			finalizing_.push_back(client_t);
		}
		// close after erase, so it does not get resued by the system
		close_socket(client_fd);
		//fprintf(stderr,"PatternSourceServiceFactory::Client> Closed connection %i\n",client_fd);
	}

	static constexpr int MAX_HEADER_SIZE = 1023;

        static int start_listening(int port, int backlog);
	static void close_socket(int fd);

        static void acceptConnections(PatternSourceServiceFactory *obj);

	const int server_port_;
	const int server_backlog_;
	const Config& config_;
	const PatternParams& pp_;
	AlnSinkSam& template_msink_;

	const std::string base_url_;

	std::mutex m_;
	PatternSourcePerThreadFactory psfact_;
	const unsigned int n_readahead_;
	LockedReadyQueueCV psq_ready_;

	std::map<int,std::thread *> clients_; 	// all active client threads
	std::vector<std::thread *> finalizing_; // inactive threads that still need to joined

	// this must be last, as it gets the obj as parameter during construction
	std::thread listent_; 			// main, listening thread
};

// not a real PatternSource, but it is convenient to keep client and server in the same source file
class PatternSourceWebClient {
private:
	class LockedOrigBufMap;
	class LockedIdleQueueCV;
public:
	// For remembering the exact input text used to define the read pair
	class OrigBuf {
	public:
		char    *readPairOrigBuf;
		uint32_t readPairOrigBuf_capacity;
		uint32_t readaOrigBuf_len; // read_a len
		uint32_t readbOrigBuf_len; // read b len
		// read names are appended behind the origbuf data
		uint32_t readaNameBuf_len; // read_a len
		uint32_t readbNameBuf_len; // read b len

		// Default, empty constructor
		OrigBuf() 
			: readPairOrigBuf(NULL), readPairOrigBuf_capacity(0)
			, readaOrigBuf_len(0), readbOrigBuf_len(0)
			, readaNameBuf_len(0), readbNameBuf_len(0)
		{}

		// copy constructor and assignment
		OrigBuf(const OrigBuf& other)
			: readPairOrigBuf(NULL), readPairOrigBuf_capacity(0)
			, readaOrigBuf_len(other.readaOrigBuf_len), readbOrigBuf_len(other.readbOrigBuf_len)
			, readaNameBuf_len(other.readaNameBuf_len), readbNameBuf_len(other.readbNameBuf_len)
		{
			const uint32_t readbOrigBufPair_len = readOrigBuf_len();
			if ((other.readPairOrigBuf!=NULL) && (readbOrigBufPair_len>0)) {
				origbuf_alloc(readbOrigBufPair_len);
				memcpy(readPairOrigBuf,other.readPairOrigBuf,readbOrigBufPair_len);
			}
		}

		OrigBuf& operator=(const OrigBuf& other) {
			readaOrigBuf_len = other.readaOrigBuf_len;
			readbOrigBuf_len = other.readbOrigBuf_len;
			readaNameBuf_len = other.readaNameBuf_len;
			readbNameBuf_len = other.readbNameBuf_len;
			const uint32_t readbOrigBufPair_len = readOrigBuf_len();
			if ((other.readPairOrigBuf!=NULL) && (readbOrigBufPair_len>0)) {
				origbuf_alloc(readbOrigBufPair_len);
				memcpy(readPairOrigBuf,other.readPairOrigBuf,readbOrigBufPair_len);
			}

			return *this;
		}

		// move constructor and assignment
		OrigBuf(OrigBuf&& other)
			: readPairOrigBuf(std::exchange(other.readPairOrigBuf,(char*)NULL)), readPairOrigBuf_capacity(std::exchange(other.readPairOrigBuf_capacity,0))
			, readaOrigBuf_len(std::exchange(other.readaOrigBuf_len,0)), readbOrigBuf_len(std::exchange(other.readbOrigBuf_len,0))
			, readaNameBuf_len(std::exchange(other.readaNameBuf_len,0)), readbNameBuf_len(std::exchange(other.readbNameBuf_len,0))
		{
		}

		OrigBuf& operator=(OrigBuf&& other) {
			readPairOrigBuf = std::exchange(other.readPairOrigBuf,(char*)NULL);
			readPairOrigBuf_capacity = std::exchange(other.readPairOrigBuf_capacity,0);
			readaOrigBuf_len = std::exchange(other.readaOrigBuf_len,0);
			readbOrigBuf_len = std::exchange(other.readbOrigBuf_len,0);
			readaNameBuf_len = std::exchange(other.readaNameBuf_len,0);
			readbNameBuf_len = std::exchange(other.readbNameBuf_len,0);

			return *this;
		}

		// cleanup 
		void reset() {
			if (readPairOrigBuf!=NULL) delete[] readPairOrigBuf;
			readPairOrigBuf = NULL;
			readPairOrigBuf_capacity = 0;
			readaOrigBuf_len = 0;
			readbOrigBuf_len = 0;
			readaNameBuf_len = 0;
			readbNameBuf_len = 0;
		}

		~OrigBuf() {
			reset();
		}

		void saveOrigBufs(const Read& read_a, const Read& read_b);

		// access
		uint32_t readOrigBuf_len() const {return readaOrigBuf_len+readbOrigBuf_len+readaNameBuf_len+readbNameBuf_len+4;}

		const char *name_buf(bool get_read_b) const {
			if (readPairOrigBuf==NULL) return NULL;
			uint32_t off = readaOrigBuf_len+readbOrigBuf_len+2;
			if (get_read_b) off += readaNameBuf_len+1;
			return readPairOrigBuf+off;
		}

		uint32_t name_len(bool get_read_b) const { return get_read_b? readbNameBuf_len : readaNameBuf_len;}

		const char *orig_buf(bool get_read_b) const {
			if (readPairOrigBuf==NULL) return NULL;
			uint32_t off = 0;
			if (get_read_b) off += readaOrigBuf_len+1;
			return readPairOrigBuf+off;
		}

		uint32_t orig_len(bool get_read_b) const { return get_read_b? readbOrigBuf_len : readaOrigBuf_len;}

	private:
		void origbuf_alloc(uint32_t size);

	};

	class ReadElement {
	private:
		// Buffer containing the tab6-formatted string
		char    *tab6_str;
		uint32_t tab6_capacity;
		uint32_t tab6_len;

		// For remembering the exact input text used to define the read pair
		OrigBuf readPairOrigBuf;
 
	public:
		ReadElement() 
			: tab6_str(NULL), tab6_capacity(0), tab6_len(0)
			, readPairOrigBuf()
		{}

		~ReadElement() {
			reset();
		}

		ReadElement(const ReadElement& other)
			: tab6_str(NULL), tab6_capacity(0), tab6_len(other.tab6_len) 
			, readPairOrigBuf(other.readPairOrigBuf)
		{
			if (tab6_len>0) {
				// alloc only what we need, even if other buf was much larger for historical reasons
				tab6_alloc(tab6_len);
				memcpy(tab6_str, other.tab6_str,tab6_len);
			} else if (other.tab6_str!=NULL) {
				// since tab6_str==NULL is special, make sure mine is not null, either
				tab6_alloc(1);
			}
		}

		ReadElement& operator=(const ReadElement& other) {
			tab6_len = other.tab6_len;
			if (tab6_len>0) {
				// alloc only what we need, even if other buf was much larger for historical reasons
				tab6_alloc(other.tab6_len);
				memcpy(tab6_str, other.tab6_str,tab6_len);
			} else if (other.tab6_str!=NULL) {
                                // since tab6_str==NULL is special, make sure mine is not null, either
                                tab6_alloc(1);
                        }
			readPairOrigBuf = other.readPairOrigBuf;

			return *this;
		}

		ReadElement(ReadElement&& other)
			: tab6_str(std::exchange(other.tab6_str,(char*)NULL)), tab6_capacity(std::exchange(other.tab6_capacity,0)), tab6_len(std::exchange(other.tab6_len,0)) 
			, readPairOrigBuf(std::move(other.readPairOrigBuf))
		{
		}

		ReadElement& operator=(ReadElement&& other) {
			tab6_str = std::exchange(other.tab6_str,(char*)NULL);
			tab6_capacity= std::exchange(other.tab6_capacity,0);
			tab6_len = std::exchange(other.tab6_len,0);
			readPairOrigBuf = std::move(other.readPairOrigBuf);
			return *this;
		}

		// fill this buffer with tab6 data from the two reads
		void readPair2Tab6(PatternSourceWebClient::LockedOrigBufMap& obmap, const Read& read_a, const Read& read_b);

		bool empty() { return tab6_str==NULL; }

		char *buf() { return tab6_str; }
		const char *buf() const { return tab6_str; }
		uint32_t buf_len() const { return tab6_len; }

		const char *name_buf(bool get_read_b) const { return readPairOrigBuf.name_buf(get_read_b); }
		uint32_t name_len(bool get_read_b) const { return readPairOrigBuf.name_len(get_read_b); }

		const char *orig_buf(bool get_read_b) const { return readPairOrigBuf.orig_buf(get_read_b); }
		uint32_t orig_len(bool get_read_b) const { return readPairOrigBuf.orig_len(get_read_b); }

		// =======================
		// mostly for internal use

		void clear_and_alloc(uint32_t size);
		// assumes the buffer is already allocated and large enough
		void append(const char *str, uint32_t str_len);
		// assumes the buffer is already allocated and large enough
		void append(const char chr);

		void reset() {
			if (tab6_str!=NULL) delete[] tab6_str;
			tab6_str = NULL;
			tab6_capacity=0;
			tab6_len=0;
			readPairOrigBuf.reset();
		}
	private:
		void tab6_alloc(uint32_t size);
	};

	class Config {
	public:
		Config(const char *iname) : 
			index_name(iname),
			seedLen(0),
			maxDpStreak(0),
			seedRounds(0),
			khits(-1) {}

		const char *index_name;
		int seedLen;
		int maxDpStreak;
		int seedRounds;
		int khits; // set to -1 if allHits==true
	};

	PatternSourceWebClient(
		const char *server_hostname,
		int server_port,
		OutFileBuf& obuf,
		const Config& config,
		int n_writecache):
		server_hostname_(server_hostname),
		server_port_(server_port),
		config_(config),
		obuf_(obuf),
		obmap_(),
		psq_empty_(n_writecache),
		psq_send_(),
		hasErrors_(false),
		isConnected_(false),
		socket_fd_(fdInit(this)),
		sendt_(sendDataWorker, this),
		receivet_(receiveDataWorker, this) {}

	~PatternSourceWebClient() {
		finalize(false);
		if (sendt_.joinable()) sendt_.join();
		if (receivet_.joinable()) receivet_.join();
		if (isConnected_ && (socket_fd_>=0)) close(socket_fd_);
	}

	// will turn false on error
	bool goodState() const {return isConnected_ && (!hasErrors_);}

	bool isConnected() const {return isConnected_;}

	// will return false on error
	bool addReadPair(const Read& read_a, const Read& read_b) {
		if (!goodState()) return false;
		ReadElement el(psq_empty_.pop());
		if (goodState()) {
			el.readPair2Tab6(obmap_,read_a,read_b);
			psq_send_.push(el);
		}
		return goodState();
	}

	// Wait for all data to be returned
	// Return false if anything went wrong
	bool finalize(bool patient=true) {
		ReadElement el;
		psq_send_.push(el); // this will signal the sendt_ thread to terminate
		if (patient) {
			// once the receiving thread terminates, we know we are done
			if (receivet_.joinable()) receivet_.join();
		}
		return goodState();
	}

	static constexpr int RE_PER_PACKET = 40; // can be small-ish, we just need enough for a couple IP packets, and each line is O(100) bytes

private:
	/*
	 * Use a ordered list to remember the original read names and bufs
	 * since we are sending only the id to the server.
	 *
	 * To minimize lookup costs, we use a direct idx->buf linear mapping,
	 * only resetting the counter when we empty the buffer.
	 *
	 * To avoid waiting (most of the times), we use two internal buffers
	 * so one is typically filled and the other is being emptied.
	 */
	class LockedOrigBufMap {
	public:
		static constexpr uint16_t BUF_CAPACITY = 10000;

		LockedOrigBufMap()
			: buf0_used_idx(0), buf1_used_idx(0)
			, buf0_used_cnt(0), buf1_used_cnt(0)
		{} // use default contructor for the rest

		// save the origbufs and names in the map
		// return index on how to access it
		uint16_t saveOrigBufs(const Read& read_a, const Read& read_b) {
			OrigBuf ob;
			ob.saveOrigBufs(read_a, read_b);
			return take_ownership(std::move(ob));
		}

		// Retrieval API
		// Caller responsible for ensuring idx is (still) valid
		const OrigBuf& lookup(uint16_t idx) const { 
			// no locking needed
			return (idx<BUF_CAPACITY) ? buf0[idx] : buf1[idx-BUF_CAPACITY]; 
		}

		OrigBuf release(uint16_t idx) {
			std::unique_lock<std::mutex> lk(m_);
			OrigBuf out = std::move((idx<BUF_CAPACITY) ? buf0[idx] : buf1[idx-BUF_CAPACITY]);
			if (idx<BUF_CAPACITY) {
				buf0_used_cnt--;
				if (buf0_used_cnt==0) {
					// since we emptied the buffer, we can write again
					buf0_used_idx = 0;
					cv_.notify_all();
				}
			} else {
				buf1_used_cnt--;
				if (buf1_used_cnt==0) {
					// since we emptied the buffer, we can write again
					buf1_used_idx = 0;
					cv_.notify_all();
				}
			}
			return out; // ideally this will be a move after optimization
		}

	private:
		uint16_t take_ownership(OrigBuf&& el) {
			std::unique_lock<std::mutex> lk(m_);
			cv_.wait(lk, [this] { return (buf0_used_idx<BUF_CAPACITY) || (buf1_used_idx<BUF_CAPACITY);});
			uint16_t myidx;
			if (buf0_used_idx<BUF_CAPACITY) {
				myidx = buf0_used_idx;
				buf0[myidx] = std::move(el);
				buf0_used_idx++;
				buf0_used_cnt++;
			} else {
				assert(buf1_used_idx<BUF_CAPACITY);
				myidx = buf1_used_idx;
				buf1[myidx] = std::move(el);
				buf1_used_idx++;
				buf1_used_cnt++;
				myidx+=BUF_CAPACITY;
			}

			return myidx;
		}

		uint16_t buf0_used_idx;
		uint16_t buf1_used_idx;
		uint16_t buf0_used_cnt;
		uint16_t buf1_used_cnt;
		std::array<OrigBuf,BUF_CAPACITY> buf0;
		std::array<OrigBuf,BUF_CAPACITY> buf1;
		std::mutex m_;
		std::condition_variable cv_;
	};

	template <typename T>
	class LockedQueueCV {
	public:
		virtual ~LockedQueueCV() {}

		bool empty() {
			bool ret = false;
			{
				std::unique_lock<std::mutex> lk(m_);
				ret = q_.empty();
			}
			return ret;
		}

		void push(T& ps) {
			std::unique_lock<std::mutex> lk(m_);
			T my = move(ps); // take ownership of the internal buffers
			q_.push(my);
			cv_.notify_all();
		}

		void pushN(int buf_len, T buf[]) {
			std::unique_lock<std::mutex> lk(m_);
			for (int i=0; i<buf_len; i++) {
				T my = move(buf[i]); // take ownership of the internal buffers
				q_.push(my);
			}
			cv_.notify_all();
		}

		// wait for data, if none in the queue
		T pop() {
			std::unique_lock<std::mutex> lk(m_);
			cv_.wait(lk, [this] { return !q_.empty();});
			T ret(q_.front());
			q_.pop();
			return ret;
		}

		// wait for data, if none in the queue
		void popUpToN(int max_read, T buf[], int& buf_len) {
			std::unique_lock<std::mutex> lk(m_);
			cv_.wait(lk, [this] { return !q_.empty();});
			buf_len = 0;
			while (buf_len<max_read) {
				buf[buf_len] = q_.front();
				buf_len++;
				q_.pop();
				if (q_.empty()) break;
			}
		}

	protected:
		std::mutex m_;
		std::condition_variable cv_;
		std::queue<T> q_;
	};

	class LockedSendQueueCV : public LockedQueueCV<ReadElement> {
	public:
		virtual ~LockedSendQueueCV() {
			//while (!q_.empty()) {
			//	ReadElement tmp(q_.pop());
			//}
		}
	};

	class LockedEmptyQueueCV : public LockedQueueCV<ReadElement> {
	public:
		LockedEmptyQueueCV(size_t n) : 
			LockedQueueCV<ReadElement>() {
			for (size_t i=0; i<n; i++) {
				ReadElement tmp;
				q_.push(tmp);
			}
		}

		virtual ~LockedEmptyQueueCV() {
			//while (!q_.empty()) {
			//	ReadElement tmp(q_.pop());
			//}
		}
	};

	// write a string with retries, but silently abort on error
	static void try_write_str(int fd, const char *str) {
		int len = strlen(str);
		while(len>0) {
			int written = ::write(fd,str,len);
			if (written<1) {
				// try once more only
				written = ::write(fd,str,len);
			}

			if (written>0) {
				// get ready for next chunk, if needed
				len-= written;
				str+=written;
			} else {
				break; // just abort
			}
		}
	}

	// write a string with retries, report any errors
	// returns false in case of error
	static bool write_str(int fd, const char *str, const int str_len) {
		int len = str_len;
		while(len>0) {
			int written = ::write(fd,str,len);
			if (written<1) {
				// try once more only
				written = ::write(fd,str,len);
			}

			if (written>0) {
				// get ready for next chunk, if needed
				len-= written;
				str+=written;
			} else {
				return false;
			}
		}
		return true;
	}

	static bool write_str(int fd, const char *str) {
		int len = strlen(str);
		return write_str(fd, str, len);
	}

	// prepend string length when sending
	// also append \r\n at end of str, as per protocol (updates input string in-place)
	// returns false in case of error
	static bool write_chunked_str(int fd, char *str, const int str_len) {
		char hexlen[16];
		str[str_len] = '\r';
		str[str_len+1] = '\n';
		str[str_len+2] = 0;
		sprintf(hexlen,"%x\r\n",str_len);
		bool noerr = write_str(fd,hexlen);
		if (noerr) noerr = write_str(fd,str,str_len+2);
		return noerr;
	}

	// read until \n\n detected, discard content
	// can go over \n\n
        static void finish_header_read(int fd, char *init_buf, int init_len);

	// read until \n\n detected
	// can NOT go over \n\n
	// return true if we read right up to t\n\n
        static bool read_header(int fd, char *buf, int& buf_len);

	// extract content length from the header string
	// return -1 if cannot find it
	static long int find_content_length(const char str[]);

	// look for the terminator request
	// returns true if it finds one
	static bool find_request_terminator(const char str[]);

	static constexpr int MAX_HEADER_SIZE = 1023;

	static void close_socket(int fd);

	static bool socketConnect(int fd, struct addrinfo &res, int port);
	static bool initialHandshake(int fd, const char *hostname, int port, const Config& config);
	static bool parseHeader(int fd, const Config& config);

	// called by contructor, assumes all but fd have been initialized
	static int fdInit(PatternSourceWebClient *obj);

	// thread procedures
	static void sendDataWorker(PatternSourceWebClient *obj);
	static void receiveDataWorker(PatternSourceWebClient *obj);

	const char *server_hostname_;
	const int server_port_;
	const Config& config_;
	OutFileBuf& obuf_;

	LockedOrigBufMap   obmap_;
	LockedEmptyQueueCV psq_empty_;
	LockedSendQueueCV  psq_send_;

	bool hasErrors_;
	bool isConnected_;
	int socket_fd_;

	// these must be last, as it gets the obj as parameter during construction
	std::thread sendt_; 			// thread sending data
	std::thread receivet_; 			// thread receiving the data
};

#ifdef USE_SRA

namespace ngs {
	class ReadCollection;
	class ReadIterator;
}

/**
 * Pattern source for reading directly from the SRA archive.
 */
class SRAPatternSource : public PatternSource {
public:
	SRAPatternSource(
		const EList<string>& sra_accs,
		const PatternParams& p) :
		PatternSource(p),
		sra_accs_(sra_accs),
		sra_acc_cur_(0),
		cur_(0),
		first_(true),
		mutex_m(),
		pp_(p)
	{
		assert_gt(sra_accs_.size(), 0);
		errs_.resize(sra_accs_.size());
		errs_.fill(0, sra_accs_.size(), false);
		open(); // open first file in the list
		sra_acc_cur_++;
	}

	virtual ~SRAPatternSource() {
		delete read_iter_;
	}

	/**
	 * Fill Read with the sequence, quality and name for the next
	 * read in the list of read files.	This function gets called by
	 * all the search threads, so we must handle synchronization.
	 *
	 * Returns pair<bool, int> where bool indicates whether we're
	 * completely done, and int indicates how many reads were read.
	 */
	virtual std::pair<bool, int> nextBatch(
		PerThreadReadBuf& pt,
		bool batch_a,
		bool lock);

	/**
	 * Finishes parsing outside the critical section
	 */
	virtual bool parse(Read& ra, Read& rb, TReadId rdid) const;

	/**
	 * Reset so that next call to nextBatch* gets the first batch.
	 * Should only be called by the master thread.
	 */
	virtual void reset() {
		PatternSource::reset();
		sra_acc_cur_ = 0;
		if (read_iter_)
			delete read_iter_;
		open();
		sra_acc_cur_++;
	}

protected:

	std::pair<bool, int> nextBatchImpl(
		PerThreadReadBuf& pt,
		bool batch_a);

	/**
	 * Open the next file in the list of input files.
	 */
	void open();

	EList<string> sra_accs_; // filenames for read files
	EList<bool> errs_;       // whether we've already printed an error for each file
	size_t sra_acc_cur_;     // index into infiles_ of next file to read
	size_t cur_;             // current read id
	bool first_;

	ngs::ReadIterator* read_iter_;

	/// Lock enforcing mutual exclusion for (a) file I/O, (b) writing fields
	/// of this or another other shared object.
	MUTEX_T mutex_m;

	PatternParams pp_;
};

#endif

#endif /*PAT_H_*/
