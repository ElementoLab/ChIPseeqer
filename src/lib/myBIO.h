/*
 * =====================================================================================
 *
 *       Filename:  myBio.h
 *
 *    Description:  classes related to Bio
 *
 *        Version:  1.0
 *        Created:  05/21/2010 01:44:23 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  David Soong (),
 *        Company:
 *
 * =====================================================================================
 */
#ifndef MYBIO_H
#define MYBIO_H

#include "myio.h"
#include <cstdio>

extern "C"{
	#include "zlib.h"
	#include "third_party/samtools/sam.h"
	#include "third_party/samtools/faidx.h"
}

class gene_str {
	private:
		string repTranscript;	// representative transcript or gene name

	public:
		int start;	// transcription start (5` end whether it is on + or -)
		int end;	// transcription end
		int cdst;	// coding start
		int cden;	// coding end

		string chr;	// chromosome
				// TO DO: maybe change to short int to save space

		char strand;		// strand +/-

		vector <int> exst;	// exon starts	--- NOTE: 0-based in genome browser data
		vector <int> exen;	// exon ends	--- NOTE: 1-based in genome browser data

		gene_str(){
			start= -1;
			end= -1;
			cdst= -1;
			cden= -1;
			strand= 0;
			repTranscript= "";
			strand= 0;
		}

		gene_str(int st, int en, int cds, int cde, char aStrand){
			start= st;
			end= en;
			cdst= cds;
			cden= cde;
			strand= aStrand;
			cerr << "NOT DONE\n";
		}

		gene_str(string st, string en, string cds, string cde, char aStrand){
			start= atoi(st.c_str());
			end= atoi(en.c_str());
			cdst= atoi(cds.c_str());
			cden= atoi(cde.c_str());
			strand= aStrand;
			repTranscript= "";
		}

		void setChr(string c){
			chr= c;
		}

		void setRepTranscript(string tp){
			repTranscript= tp;
		}

		string getRepTranscript(){
			return repTranscript;
		}

		void addExons(const string &numex, const string &st, const string &en);
		bool inExon(const int &pos);
		int getExonId(const int &pos);
		void printExons();

		inline bool isCDRange(const int &pos);
		inline bool isUTR5(const int &pos);
		inline bool isUTR3(const int &pos);

		bool isNULL(){
			return (start == -1 );			// if initialized => start > 0
		}

		/*
		 * TO DO:
		 * 	currently we use the first transcript to represent a gene's structure
		 * 	maybe we can later merge all of a gene's exons to build a more comprehensive model
		 */
		void mergeExons(){
			cerr << "not implemented!! "<< endl;
		}
};

/*
 * between coding start and coding end??
 */
bool gene_str::isCDRange(const int & pos){
	return ( cdst <= pos && pos < cden);
}

bool gene_str::isUTR5(const int &pos){
	if (strand == '+'){
		return (pos >= start && pos < cdst);
	}
	return (pos< end && pos >= cden);
}

bool gene_str::isUTR3(const int &pos){
	if (strand == '-'){
		return (pos >= start && pos < cdst);
	}
	return (pos< end && pos >= cden);
}




/*
 * a data structure storing all RefSeq genes
 */
class genemodel{

	private:
		map <string, map<string, int> > name2refid;	// gene name -> refSeq transcript id
		map <string, gene_str> name2genestr;		// gene name -> gene structure
		map <string, map<string, int> > chr2name;	// chr -> gene names
		map <string, map <int, set <string> > > chrbins2name;	// chromosome bins, each is a small region on the
									// genome and contains a small number of genes (-> speeds up search)
									// 	(instead of using getGenesOnChr, use getGenesInBin())

		// char strand;
		// int  txStart;
		// int  txEnd;
		int BIN_SIZE;
		bool bUseTranscriptName;				// use transcript name or gene name as key
		bool bLoadChr;						// load chromosome?

	public:
		map <string, string> nm2np;				// information from refLink
		map <string, string> np2nm;				// information from refLink

		genemodel(){
			BIN_SIZE= 100000;
			bUseTranscriptName= false;			// by default uses the gene NAME as key
			bLoadChr= false;				// by default doesn't load the chromosome to save memory
		}

		inline int getBinSize(){return BIN_SIZE;}		// bin size in binning the chromosomes
		inline int getBinIdx(int pos){return int(pos/BIN_SIZE);}	// return bin index of pos

		void load(const string & fname);			// load refSeq data file
		void loadRefLink(string fname);				// load NM <-> NP mapping file

		string getTranscripts(const string &genename);		// get transcripts of gene genename

		vector <string> getAllGenes();				// get all genes (or transcript names)
		vector <string> getGenesOnChr(string chr);		// get all genes on chromosome chr

		vector <string> findGenesAt(const string &chr, const int &pos);
		vector <string> findGenesAtPosFromSet(const vector <string> &geneset, int pos);
		string findClosestGene(const string &chr, const int &pos, int &mdist);
		string findFirstGeneAtPosFromSet( const set <string> &geneset, int pos);

		void printGeneStr(const string &genename);		// print gene exon structure

		gene_str getGeneStr(const string &genename);

		map <int, string> getChrGeneStarts(string chr);
		set <string> getGenesInBin(const string chr, int binidx);
		set <string> getGenesInRegion(const string chr, const int &st, const int &en){
		    cerr << "WARNING: getGenesInRegion() NOT IMPLEMENTED!"<< endl;
		}
		vector <string> getChromosomes();

		void setBinsize(int s){
			BIN_SIZE= s;
		}
		void setLoadChr(bool b){
			bLoadChr= b;
		}
		void setUseTranscript(bool key){
			bUseTranscriptName= key;
		}
		friend class gLookup;
};

inline vector <string> genemodel::getChromosomes(){
	return getKeys <string, map<string, int> > (chr2name);
}

inline set <string> genemodel::getGenesInBin(string chr, int binidx){
	return chrbins2name[chr][binidx];
}

/*
 * class to find arbitray objects on the genome
 *
 * objects= restriction enzyme site, TF binding site, enhancer mark, interaction site
 *
 * TO DO:
 * 	continue here
 */
class gLookup: public genemodel{
	public:
		void load(const string &fname);
		void test(){
			cout << "gLookup here\n";
		}
};

/*
 * class to read and store sam data
 * 	focus: alignment level
 */
class sam{
	private:
		bool b_unique;
		bool b_pileup;
		bool b_ntidx;
		bool b_err_only;					// get err rates only (i.e. won't get pileup, read counts, nt indexes)
		map <int, vector <char> > sam_pos2pile;		// store pileup information
		map <int, vector <int> > sam_pos2nt_idx;	// store nt position in read
		map <int, int> sam_pos2cnt;			// assuming we are on the same chromosome
		bool verbose;
		bool bShowAllPos;				// in pileup, show all positions even if they have 0 read count
		vector <int> mismatch_cnt;			// mismatch count
		int num_reads;
		int read_len;

	public:
		sam(){
			b_unique= true;		// by default use unique reads only
			b_pileup= false;	// by default don't generate pileup
			b_ntidx= false;		// by default don't load nt positions
			b_err_only= false;
			verbose= false;		//
			bShowAllPos= false;	//
			num_reads= 0;
			read_len= 0;

		};
		void load(const string samfname, string tgChr="");
		int getReadCount(int sampos);
		map <int, int> getPos2Cnt(){
			return sam_pos2cnt;
		}
		string getPos2PileSeq(const int &sampos){
			return vector2string2 <char> (sam_pos2pile[sampos], "");
		}
		map <string, int> getGeneReadCount(string gene, genemodel &gm);
		void appMappability(string chr_list, string alidir);
		int getReadLength(){return read_len;}
		static int getReadLength(const string & samfname);
		int getNumReads(){return num_reads;}
		void parseCigar(const string & readseq, const string & cigar, map <int, int> &sam_pos2cnt, map <int, vector <char> > & sam_pos2pile, const int &offset);
		void parseMDtag(const string & mdtag);
		void reset();
		void showPileup(const string &chr, const int &start, const int &end, const string hg18file= "");
		vector <int> get_mismatch_counts();
		void show_mismatch_rate();
		void save_mismatch_rate(const string &fname);
		inline void setErrOnly(bool b_err);
		inline void setUnique(bool unq);
		inline bool getUnique(){return b_unique;}
		inline void setPileup(bool pile);
		inline void setNtIdx(bool nt);
		~sam(){
			// cout << "cleaning up SAM\n";
			reset();
		}
		// utility funcitons
		static void findMachineNameX0 (const string &fname, string & machine_code, int & mapcnt_idx);
		static map <string, map<int, int> > loadAllChrSAM(const string &samdir, const string &prefix, const string &suffix);
		static map <string, read_str> loadAlignedReads(const string &fname, map <string, map < int, int > > & counts, map <string, read_str> &subset, bool bUniqMap, int &read_len);
};

inline void sam::setErrOnly(bool b_err){
	b_err_only= b_err;
}

inline void sam::setNtIdx(bool nt){
	b_ntidx= nt;
}

inline void sam::setUnique(bool unq){
	b_unique= unq;
}

inline void sam::setPileup(bool pile){
	b_pileup= pile;
}

/*
 * get read count at position pos from a sam file
 */
inline int sam::getReadCount(int sampos){
	// cerr << "WARNING: not implemented!\n";
	return sam_pos2cnt[sampos];
}



/*
 * class to hold fasta/fastq objects
 */
class sequence{
	private:
		map <string, string> id2seq;
		bool bLoadSeq;
		faidx_t * fidx;
	public:
		set <string> alignedids;
		sequence(){
			fidx= 0;
			bLoadSeq= true;			// by default we load the sequences, but this can
							// also be disabled to load only the sequence ids
		}
		void loadFasta(const string & fname, set <string> subset= set <string>());
		void loadFastq(const string & fname, set <string> subset= set <string>());
		void loadFastqIds(string fname);
		void loadAlignedIds(const string & fname);
		int getNumSeq(){
			return id2seq.size();
		}
		string getSeq(const string & id);

		string getHG18seq(const string &region);
		void loadHG18(string refname="");
		int getLen(string id);
		vector <string> getIds();					// get all ids
		void setLoadSeq(bool loadseq){					// load sequences? or we just load the ids
			bLoadSeq= loadseq;
		}
		void print();
		void reset();
		~sequence(){
			reset();
		}

		// utilities
		static void countACGT(const string &dna, int *& freq);
		static void findPolyA(const string &seq, int &lastpos, int &len);
		static void findPolyT(const string &seq, int &lastpos, int &len);
		static int countN(const string &seq, const char &N, const int pos);
		static string getComplement(string dna);
		static string getReverse(string dna);
		static string getRevComplement(string dna);
};

inline string sequence::getSeq(const string & id){
	if (id2seq.count(id)== 1){
		return id2seq[id];
	}else{
		return "";
	}
}

inline vector <string> sequence::getIds(){
	return getKeys <string, string> (id2seq);
}


/*
 * class for genome browser related stuff
 */
class gbrowser{
	public:
		gbrowser(){};
		static void appPairedEndPlot(string file1, string file2, string outfile, int max_dist, bool bUniqRead);
		static void appendToBedGraph(ofstream &fout, string chr, int chrlen, map <int, int> allpos);
		static void appendToBedGraph(ofstream &fout, string chr, int chrlen, vector <int> allpos);
		static void appendToBedGraph(ofstream &fout, string chr, map <int, int> allpos);
		static string countToBedGraph(string chr, int chrlen, map <int, int> allpos);
		static string countToBedGraph(string chr, int chrlen, vector <int> &allpos);
};

class pileup{
	private:
		map <string, map <int, streamoff> > allpos2idx;		// chr position -> closest file index
		map <string, map <int, streamoff> > allpos2blocksize;		// chr position -> closest file index
		int block_size;						// how many bases (on chr) in a bin
		int maxblocksize;					// size of a file chunk
		char * buffer;						// to store a file chunk
		string	cachedChr;
		string buffer_str;
		string pile_fname;
		map <string, map <int, int> > cachedReads;
		ifstream fh_pile;		// (pile_fname.c_str());

		// global tokens to avoid repeated initializations
		vector <string> gTokens;
		vector <string> gTmpTokens;


	public:
		pileup(){
			// TO DO:
			// 	maybe increase block_size since we have much more RAM
			//
			block_size= 100000;				// unit= bp on chromosomes
			maxblocksize= 10000;				// size of a chunk
		}
		void load(string fname);
		void index(string fname);
		void loadIndex(string fname);
		int  fetchReadCount(string chr, int sampos);
		void getRPKM(string fname, genemodel &gm);

		// TO DO:
		// deconstructor goes here
		~pileup(){
			cout << "deleting pileup now...\n";
			// TO DO:
			// 	close index file
			// 	close pileup file if loaded
			if (fh_pile.is_open()){
				cout << "closing "<< pile_fname << endl;
				fh_pile.close();
			}
			if (buffer){
				free(buffer);
				buffer_str.clear();
			}
		}
};

char codon2aa(const string & codon);
void parseRegion(const string &reg, string &chr, int &st, int &en);

#endif
