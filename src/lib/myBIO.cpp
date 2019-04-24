/*
 * =====================================================================================
 *
 *       Filename:  myBIO.cpp
 *
 *    Description:  Bio-related classes
 *
 *        Version:  1.0
 *        Created:  05/21/2010 01:46:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  David Soong (),
 *        Company:
 *
 * =====================================================================================
 */

#include "myBIO.h"

/*
 * parse a free-text string to retrieve chr, start, end
 */
void parseRegion(const string &reg, string &chr, int &st, int &en){
	vector <string> tokens;
	// TO DO:
	// 	deal with comma-padded numbers
	Tokenize(reg, tokens, ":-");	// chr:111-222
	chr= tokens[0];
	st= atoi(tokens[1].c_str());
	en= atoi(tokens[2].c_str());
}


/*
 * check if pos is in an exon
 */
bool gene_str::inExon(const int& pos){
	for (int i=0; i< (int)exst.size(); i++){
		// note exen (end position) is 1-based
		if (exst[i] <= pos && pos < exen[i])
			return true;
	}
	return false;
}

/*
 * return which exon this positon falls in
 *
 * NOTE:
 * 	exon id here goes from the 5' -> 3' order.
 * 	maybe do it in 3' -> 5' for the reverse strand in ANOTHER FUNCTION!!!
 * 	(don't modify this function, dna->mRNA->protein transcription/translation depends on it)
 *
 */
int gene_str::getExonId(const int& pos){
	for (int i=0; i< (int)exst.size(); i++){
		// note exen (end position) is 1-based
		if (exst[i] <= pos && pos < exen[i])
			return i;
	}
	return -1;
}

/*
 * print the exon structure
 */
void gene_str::printExons(){
	cout << "strand: "<< strand << endl;
	cout << "tx start: "<< start<< endl;
	cout << "cd start: "<< cdst << endl;
	for (int i=0; i< (int)exst.size(); i++){
		printf("Exon %d: %d - %d\n", i, exst[i], exen[i]);
	}
	cout << "cd end: "<< cden << endl;
	cout << "tx end: "<< end << endl;
}

/*
 * add exons from refseq table (source: ucsc genome browser)
 */
void gene_str::addExons(const string &numex, const string &st, const string &en){
	vector <string> vExst;
	vector <string> vExen;
	Tokenize(st, vExst, ",");
	Tokenize(en, vExen, ",");
	int num_exons= atoi(numex.c_str());
	for (int i=0; i < num_exons; i++){
		exst.push_back(atoi(vExst[i].c_str()));
		exen.push_back(atoi(vExen[i].c_str()));
	}
	vExst.clear();
	vExen.clear();
}


/*
 * return the gene structure given a gene name
 * if no match, return an uninitialized gene_str
 * 	=> check for returned gene_str to see if it is NULL
 */
gene_str genemodel::getGeneStr(const string &genename){
	if (name2genestr.count(genename)!= 0){
		return (name2genestr[genename]);
	}
	else{
		gene_str x;
		return x;
	}
}

/*
 * find genes that overlap with position pos on chrosome chr
 */
vector <string> genemodel::findGenesAtPosFromSet(const vector <string> &geneset, int pos){
	vector <string> output= vector <string> ();		// output vector

	for( int i=0; i< (int)geneset.size(); i++){
		// check if positition pos falls in range
		if (name2genestr[geneset[i]].start <= pos && name2genestr[geneset[i]].end > pos){
			output.push_back(geneset[i]);
		}
	}
	return output;
}

/*
 * looks for the first gene that overlaps with this genomic position
 */
string genemodel::findFirstGeneAtPosFromSet(const set <string> &geneset, int pos){
	string output= "";

	set <string> ::iterator it;
	// for( int i=0; i< (int) geneset.size(); i++){
	for ( it= geneset.begin() ; it != geneset.end(); it++ ){
		// check if positition pos falls in range
		// if (name2genestr[geneset[i]].start <= pos && name2genestr[geneset[i]].end > pos){
		//
		if (name2genestr[*it].start <= pos && name2genestr[*it].end > pos){
			// output= geneset[i];
			output= *it;
			return output;
		}
	}
	return output;
}

/*
 * find the gene that's closest to a given position
 */
string genemodel::findClosestGene(const string &chr, const int &pos, int & mdist){
	vector <string> genes= getGenesOnChr(chr);		// genes on chromosome chr
	int mindist= INT_MAX;
	int dist, d1, d2;
	string output= "";

	vector <string>::iterator it;
	// go through the gene list one by one
	for( it = genes.begin(); it != genes.end(); it++ ){
		d1= name2genestr[*it].start - pos;
		d2= name2genestr[*it].end - pos;
		if (d1 <= 0 && d2 >= 0){
			// position in gene body
			output= *it;
			mindist= 0;
			break;
		}else{
			dist= min( abs(d1), abs(d2) );
			if (dist < mindist){
				mindist= dist;
				output= *it;
			}
		}
	}
	mdist= mindist;
	return output;
}

/*
 * find genes that overlap with position pos on chrosome chr
 */
vector <string> genemodel::findGenesAt(const string &chr, const int &pos){
	vector <string> genes= getGenesOnChr(chr);		// genes on chromosome chr
	vector <string> output= vector <string> ();		// output vector

	vector <string>::iterator it;
	for( it = genes.begin(); it != genes.end(); it++ ){
		// check if positition pos falls in range
		// printGeneStr(*it);
		if (name2genestr[*it].start <= pos && name2genestr[*it].end > pos){
			output.push_back(*it);
			// cout << chr << " pos " << pos<< ": "<< *it << " " << name2genestr[*it].start << " - " << name2genestr[*it].end << endl;
		}

	}
	return output;
}


/*
 * print the gene structure of genename
 */
void genemodel::printGeneStr(const string & genename){
	if (name2genestr.count(genename)==0){
		// not found
		cout << genename << " not found\n";
	}else{
		// now print gene structure
		gene_str tmp_gene= name2genestr[genename];
		cout << genename << ": start= " << tmp_gene.start << " end= "<< tmp_gene.end << endl;
	}

}

/*
 * get all gene/transcript names (depending on the flag bUseTranscriptName)
 */
vector <string> genemodel::getAllGenes(){
	vector <string> v= getKeys <string, gene_str> (name2genestr);
	return v;
}

/*
 * get genes on chromosome chr
 *
 * TO DO:
 * 	should sort by chromosomal position to speed up
 */
vector <string> genemodel::getGenesOnChr(string chr){
	// chr2name

	if (chr2name.count(chr) > 0){
		// there are genes on chromosome chr
		vector <string> v= getKeys <string, int> (chr2name[chr]);
		return v;
	}
	else{
		return vector <string>();
	}
}

/*
 * get transcripts belonging to genename
 *
 */
string genemodel::getTranscripts(const string & genename){
	// use chr2name

	cerr << genename << " exists? = " << name2refid.count(genename) << endl;

	if (name2refid.count(genename) > 0){
		// use iterator
		   for( map<string, int>::iterator ii= name2refid[genename].begin(); ii!= name2refid[genename].end(); ++ii)
		   {
		       cerr << (*ii).first << ": " << (*ii).second << endl;
		   }
	}

	return "";

}

/*
 * get the gene start sites (map <int, string>) for all genes on chromosome che
 */
map <int, string> genemodel::getChrGeneStarts(string chr){
	vector <string> allgenes= getGenesOnChr(chr);
	map <int, string> allstarts;
	gene_str tmpgene;
	for (int i=0; i< (int) allgenes.size(); i++){
		tmpgene= getGeneStr(allgenes[i]);
		allstarts[tmpgene.start]= allgenes[i];
	}
	return allstarts;
}

void genemodel::loadRefLink(string fname){
	cerr << "loading refLink from "<< fname << endl;

	string line;
	vector <string> tokens;
	fileSlurp fs(fname);
	while(! fs.eof()){
		line= fs.readLine();
		Tokenize(line, tokens, "\t");
		if (tokens[0].empty() || tokens[1].empty()){
			continue;
		}
		nm2np[tokens[0]]= tokens[1];	// TO DO: here we assume a 1-1 relationship. maybe refine.
		np2nm[tokens[1]]= tokens[0];
	}
	fs.close();
}


/*
 * load a UCSC gene model (refGene...)
 * 	NOTE:
 * 		starting coordinates are 0-based
 * 		ending coordinates are 1-based
 */
void genemodel::load(const string &fname){
	std::cerr << "loading "+ string(fname) << "\n";
	typedef map<string, gene_str> mapType;

	string line;
	int binidx1;
	int binidx2;

	string tmpname;
	string transcript;
	string chr;
	char strand;
	int tx_st;
	int tx_en;

	ifstream myfile(fname.c_str());
	if (myfile.is_open()){
		vector <string> tokens;
		while (! myfile.eof() )
		{
			getline (myfile,line);
			Tokenize(line, tokens, "\t");

			if (tokens.size()!= 0){
			    // add information
			    // TO DO:
			    // 	should be able to load transcript information by adding a switch here
			    tmpname= tokens.at(12);
			    transcript= tokens.at(1);
			    if (bUseTranscriptName){
				    tmpname= tokens.at(1);		// use transcript name as key
				    transcript= tokens.at(12);		// store gene name in the repTranscript slot
			    }
			    chr= tokens.at(2);
			    strand= tokens.at(3)[0];
			    tx_st= atoi(tokens.at(4).c_str());
			    tx_en= atoi(tokens.at(5).c_str());
			    // int cd_st= atoi(tokens.at(6).c_str());
			    // int cd_en= atoi(tokens.at(7).c_str());

			    // gene -> mulriple transcripts
			    // name2refid.insert(pair<string, string>(tmpname, transcript));
			    // name2refid[tmpname][transcript]= 1;

			    // chr -> genes
			    chr2name[chr][tmpname]= 1;
			    // cout << "adding chr2name " << chr << " " << tmpname << " = " << chr2name[chr].size() <<  endl;

			    // bin -> genes
			    // TO DO:
			    // 	check how many genes are in each bin
			    binidx1= int(tx_st/BIN_SIZE);
			    binidx2= int(tx_en/BIN_SIZE);
			    for (int tmpbin= binidx1; tmpbin <= binidx2; tmpbin++){
				if (chrbins2name[chr][tmpbin].size()== 0){
					set <string> tmpgenes;
					tmpgenes.insert(tmpname);
					chrbins2name[chr][tmpbin]= tmpgenes;
					// cout << "initializing " << chr << " bin= " << tmpbin << " with "<< tmpname << endl;
				}else{
					chrbins2name[chr][tmpbin].insert(tmpname);

					// cout << "adding "<< tmpname << " to " << chr << " bin= " << tmpbin << " current bin size= "<<
					//	chrbins2name[chr][tmpbin].size() <<endl;
				}
			    }

			    // retrieve gene structure
			    // 	gene -> start point
			    // 	gene -> end point
			    mapType::iterator iter;
			    iter= name2genestr.find(tmpname);
			    if( iter != name2genestr.end() ) {
				// record located -> check if bounds are different
				// WARNING:
				// 	currently we only load the FIRST TRANSCRIPT if there are multiple isoforms for a gene
				//

				/*
					gene_str tmp_gene= iter->second;
					if (cd_st < tmp_gene.cdst){
						if (DEBUG)
							printf("change in CDS start of %s (%s) old= %d new= %d\n", transcript.c_str(), tmpname.c_str(),
									tmp_gene.cdst, tx_st);
						tmp_gene.cdst= cd_st;
					}
					if (cd_en > tmp_gene.cden){
						if (DEBUG)
							printf("change in CDS end of %s (%s) old= %d new= %d\n", transcript.c_str(), tmpname.c_str(),
									tmp_gene.cden, tx_en);
						tmp_gene.cden= cd_en;
					}
					if (tx_st < tmp_gene.start){
						// cout << "change in start of "<< transcript<< " old= " << tmp_gene.start << " new= "<< tx_st<< endl;
						tmp_gene.start= tx_st;
					}
					if (tx_en > tmp_gene.end){
						// cout << "change in end of "<< transcript<< " old= " << tmp_gene.end << " new= "<< tx_en<< endl;
						tmp_gene.end= tx_en;
					}
					// TO DO:
					// for now we only use the first transcript's exon structure !!
					// 	BUT we update the start and end sites of transcription and coding region
					// tmp_gene.addExons(tokens[8], tokens[9], tokens[10]);
					name2genestr[tmpname]= tmp_gene;
					// cout << "found "+ tmpname + "("+ transcript+") st= "<< tmp_gene.start << " end= " << tmp_gene.end << endl;
				*/
			    }
			    else{
				// cout << "adding new gene_str for "<< tmpname << '\n';
				// cout << "gene " << tmpname << " strand= "<< strand << " rep transcript= "<< transcript<< endl;
				gene_str a_gene(tokens[4], tokens[5], tokens[6], tokens[7], strand);
				a_gene.setRepTranscript(transcript);
				a_gene.addExons(tokens[8], tokens[9], tokens[10]);
				if (bLoadChr)
					a_gene.setChr(chr);
				name2genestr[tmpname]= a_gene;
			    }

			    // cout << "adding "+ tokens.at(12)<< " -> "<< tokens.at(1) << endl;
			}else{
				break;
			}
			// printf("%d elements\n", (int)tokens.size());
			// cout << line << endl;
		}
		myfile.close();
		cerr << "Total: " <<  name2genestr.size() << " genes\n";

	}
	else cerr << "Unable to open file " << fname << endl;
}

void gLookup::load(const string &fname){
	std::cerr << "loading "+ string(fname) << "\n";
	typedef map<string, gene_str> mapType;

	string line;
	int binidx1;
	int binidx2;

	string tmpname;
	string transcript;
	string chr;
	// char strand;
	int tx_st;
	int tx_en;
	int lineidx= 0;

	ifstream myfile(fname.c_str());
	if (myfile.is_open()){
		vector <string> tokens;
		while (! myfile.eof() )
		{
			getline (myfile,line);
			Tokenize(line, tokens, "\t");

			if (tokens.size()!= 0){
			    // add information
			    // TO DO:
			    // 	should be able to load transcript information by adding a switch here
			    tmpname= castToString <int> (lineidx);

			    chr= tokens.at(0);
			    tx_st= atoi(tokens.at(1).c_str());
			    tx_en= atoi(tokens.at(2).c_str());

			    // chr -> genes
			    chr2name[chr][tmpname]= 1;
			    // cout << "adding chr2name " << chr << " " << tmpname << " = " << chr2name[chr].size() <<  endl;

			    // bin -> genes
			    // TO DO:
			    // 	check how many genes are expected in each bin => adjust BIN_SIZE
			    binidx1= int(tx_st/BIN_SIZE);
			    binidx2= int(tx_en/BIN_SIZE);
			    for (int tmpbin= binidx1; tmpbin <= binidx2; tmpbin++){
				if (chrbins2name[chr][tmpbin].size()== 0){
					set <string> tmpgenes;
					tmpgenes.insert(tmpname);
					chrbins2name[chr][tmpbin]= tmpgenes;
					// cout << "initializing " << chr << " bin= " << tmpbin << " with "<< tmpname << endl;
				}else{
					chrbins2name[chr][tmpbin].insert(tmpname);

					// cout << "adding "<< tmpname << " to " << chr << " bin= " << tmpbin << " current bin size= "<<
					//	chrbins2name[chr][tmpbin].size() <<endl;
				}
			    }

			    // retrieve gene structure
			    // 	gene -> start point
			    // 	gene -> end point
			    mapType::iterator iter;
			    iter= name2genestr.find(tmpname);
			    if( iter != name2genestr.end() ) {
				// WARNING:
				// 	currently we only load the FIRST TRANSCRIPT if there are multiple isoforms for a gene
			    }
			    else{
				// cout << "adding new gene_str for "<< tmpname << '\n';
				// cout << "gene " << tmpname << " strand= "<< strand << " rep transcript= "<< transcript<< endl;
				// gene_str a_gene(tokens[4], tokens[5], tokens[6], tokens[7], strand);
				gene_str a_gene(tokens[1], tokens[2], tokens[1], tokens[2], 'x');
				name2genestr[tmpname]= a_gene;
			    }
			    lineidx++;

			    // cout << "adding "+ tokens.at(12)<< " -> "<< tokens.at(1) << endl;
			}else{
				break;
			}
			// printf("%d elements\n", (int)tokens.size());
			// cout << line << endl;
		}
		myfile.close();
		cerr << "Total: " <<  name2genestr.size() << " genes\n";

	}
	else cerr << "Unable to open file " << fname << endl;
}


void sequence::reset(){
	bLoadSeq= true;
	id2seq.clear();
	if (fidx){
		cerr << "sequence::reset\n";
		fai_destroy(fidx);
	}
}

/*
 * count the occurence of nucleotide N in the first or last pos positions
 * 	pos > 0 => first pos positions
 * 	pos < 0 => last pos positions
 * WARNING:
 * 	aborts if pos is out of bound
 */
int sequence::countN(const string &seq, const char & N, const int pos){
	// cout << "seq= " << seq << " (" << seq.length() << ") char= "<< N << " pos= "<< pos << " len+pos= "<< (seq.length()+pos)<<  endl;
	int cnt=0;
	if (pos > 0){
		for (int i=0; i< pos; i++){
			if (seq.at(i) == N)
				cnt++;
		}
	}else{
		// pos < 0  ==> last pos positions
		for (int i= seq.length()-1; i >= int(seq.length()+pos); i--){
			if (seq.at(i) == N)
				cnt++;
		}
	}
	return cnt;
}

/*
 * count occurrences of A, C, G, and T nucleotides in a dna sequence
 */
void sequence::countACGT(const string &dna, int * & freq){
	for (int i=0; i< (int) dna.length(); i++){
		switch (dna.at(i)){
			case 'A':
				freq[0]++; break;
			case 'C':
				freq[1]++; break;
			case 'G':
				freq[2]++; break;
			case 'T':
				freq[3]++; break;
			default:
				cerr << "ERROR: not ACGT sequence\n";
		}
	}
}

/*
 * retrieve a region from the hg18 human genome.
 *
 * NOTE: need to call loadHG18() first
 */
string sequence::getHG18seq(const string & region){
	int len;
	char * seq= fai_fetch(fidx, region.c_str(), &len);
	string outstr(seq);

	// printf("%s\nlength= %d\n", seq, len);
	free(seq);
	return outstr;
}

/*
 * load the hg18 SAM index file
 */
void sequence::loadHG18(string refname){
	if (refname== "")
		refname= expandHome("~/miro/databases/wg/hg18/hg18.fa");

	fidx= fai_load(refname.c_str());

	if (fidx == 0){
		cerr << "Error loading reference genome from "<< refname << "\n";
		exit(1);
	}
}

/*
 * load aligned read ids
 */
void sequence::loadAlignedIds(const string &alnfile){
	string line;
	vector <string> tokens;
	fileSlurp fs(alnfile);
	while(! fs.eof()){
		line= fs.readLine();
		if (line.at(0) == '@')
			continue;
		Tokenize(line, tokens, "\t");
		if (tokens[2] != "*"){
			// this read is aligned; TopHat only reports aligned reads
			alignedids.insert(tokens[0]);
		}
	}
	fs.close();
}

void sequence::print(){
	map <string, string>::iterator it;
	for (it= id2seq.begin(); it != id2seq.end(); it++){
		cout << ">"<< (*it).first << "\t"<< (*it).second<< endl;
	}
}

/*
 * load a FASTA sequence file;
 *
 * TO DO:
 * 	(NOT IMPLEMENTED!!)
 * 	if <subset> is provided, only load the sequences of the specified ids.
 */
void sequence::loadFasta(const string &fname, set <string> subset){
	string id= "", seq= "", line;
	bool bLoadSubset= false;

	if ( subset.size() > 0 ){
		cerr << "WARNING: NOT IMPLEMENTED loadFasta\n";
		cerr << "loading only a subset of ids\n";
		bLoadSubset= true;
	}

	fileSlurp fs(fname);
	while(! fs.eof()){
		line= fs.readLine();
		// cerr << "line: "<< line << "|" << endl;
		if ( startsWith(line, ">") ){	// old: line.at(0) == '>'){
			// cerr << "loading "<< line << endl;
			// if seq exists, store it
			if (seq != ""){
				// TO DO:
				// 	finish loading the subset
				id2seq[id]= seq;
				// cout << id << " -> "<< seq << endl;
			}
			// get id name
			id = line.substr(1);
			seq= "";
			continue;
		}

		// get sequence and store
		seq+= line;
	}
	if (seq != ""){
		id2seq[id]= seq;
	}
	fs.close();
}



/*
 * load a fastq sequence file
 * 	1. use setLoadSeq() to determine if we load only the read ids
 * 	2. when given subset => load a subset of sequences
 * 	3.
 *
 * TO DO:
 * 	load quality scores??
 * 	parse header info?
 * 		i.e. lane,tile_num,x,y,index,pair
 * 	in-memory compression??
 */
void sequence::loadFastq(const string &fname, set <string> subset){
	// vector <string> tokens;
	string id, seq, line;
	bool bLoadSubset= false;

	cerr << "subset size = " << subset.size() << endl;
	if ( subset.size() > 0 ){
		cerr << "loading only a subset of ids\n";
		bLoadSubset= true;
	}

	fileSlurp fs(fname);
	while(! fs.eof()){
		line= fs.readLine();
		if (line.at(0) == '@'){
			// get id name
			id = line.substr(1);

			// get sequence and store
			seq= fs.readLine();

			// skip two lines
			// 	TO DO:
			// 		add filtering out low-quality reads
			fs.readLine();
			fs.readLine();

			if (bLoadSeq){
				// load all sequences (~ 5 GB memory for a sequencing reads file)
				if ( bLoadSubset ){
					// load subset -> compare with subset
					if ( subset.count (id) > 0 ){
						// cout << "adding read id"<< id << endl;
						id2seq[id]= seq;
					}
				}else{
					// load everything
					id2seq[id]= seq;
				}
			}else{
				// only load the ids to save memory
				id2seq[id]= "";
			}
			// cout << ">" << id << "\n" << seq <<endl;
			continue;
		}
	}
	fs.close();
}

/*
 * load only the read ids (but still takes about 3GB memory??)
 */
void sequence::loadFastqIds(string fname){
	setLoadSeq(false);
	loadFastq(fname);
}

/*
 * check if a sequence starts with poly (T)
 * return values:
 * 	lastpos:	last matching position of poly T
 * 	len:		length of poly T
 */
void sequence::findPolyT(const string &seq, int &lastpos, int &len){
	int seqlen= seq.length();
	if (seq.at(0) != 'T' ){
		lastpos= -1;
		len= -1;
	}

	for (int i= 0; i< seqlen; i++){
		// cout << i << "\t"<< seq.at(i) << endl;
		if (seq.at(i) != 'T'){
			// cout << "\tnot T at " << i << endl;;
			lastpos= i- 1;
			break;
		}
	}
	len= (lastpos+1);
}

/*
 * check if a sequence has a poly A tail
 * return values:
 * 	lastpos:	last matching position of poly A
 * 	len:		length of poly A
*/
void sequence::findPolyA(const string &seq, int &lastpos, int &len){
	int seqlen= seq.length();
	if (seq.at(seqlen-1) != 'A' ){
		lastpos= -1;
		len= -1;
	}

	// int lastpos= -1;

	for (int i= seqlen- 1; i>= 0; i--){
		// cout << i << "\t"<< seq.at(i) << endl;
		if (seq.at(i) != 'A'){
			// cout << "\tnot A at " << i << endl;;
			lastpos= i+ 1;
			break;
		}
	}
	// cout << "last pos= "<< lastpos << ", len= "<< (seqlen- lastpos) << endl;
	(len= seqlen- lastpos);
}

string sequence::getRevComplement(string dna){
	dna= sequence::getComplement(dna);
	reverse(dna.begin(), dna.end());
	return dna;
}

string sequence::getReverse(string dna){
	reverse(dna.begin(), dna.end());
	return dna;
}

string sequence::getComplement(string dna){
	for (int i=0; i< (int)dna.length(); i++){
		switch (dna.at(i)){
			case 'A':
			case 'a':
				dna.at(i)= 'T';
				break;
			case 'T':
			case 't':
				dna.at(i)= 'A';
				break;
			case 'C':
			case 'c':
				dna.at(i)= 'G';
				break;
			case 'G':
			case 'g':
				dna.at(i)= 'C';
				break;
			case 'N':
				// do not complement Ns
				break;
			default:
				cerr << "illegal DNA string "<< dna << endl;
		}
	}
	return dna;
}

void sam::reset(){
	// TO DO:
	// reset sam_pos2cnt (map)
	sam_pos2cnt.clear();
	if (b_pileup){
		deleteMapComponents <int, vector <char> > (sam_pos2pile);	// remove map components
		sam_pos2pile.clear();						// clean up
	}
	if (b_ntidx){
		deleteMapComponents <int, vector <int> > (sam_pos2nt_idx);
		sam_pos2nt_idx.clear();
	}
}

/*
 * output the read count information for a specified region in pileup format
 */
void sam::showPileup(const string &chr, const int &start, const int &end, const string hg18file){
	// TO DO:
	// 	add option to load/not load hg18 sequences
	sequence s;
	s.loadHG18( hg18file);			// pre-load hg18 index

	string region= chr+":"+ castToString <int> (start)+ "-"+ castToString <int> (end);

	string wgseq= s.getHG18seq(region);			// NOTE: 0-based coordinates

	for (int sampos= start; sampos< end; sampos++){
		if ( !bShowAllPos && sam_pos2pile.count(sampos) > 0){
			cout << chr << "\t"<< sampos << "\t" << toupper(wgseq.at(sampos-1)) << "\t" << vector2string2 <char> (sam_pos2pile[sampos], "") <<
				"\t" << sam_pos2cnt[sampos];

			if (b_ntidx)
				cout << "\t" << vector2string2 <int> (sam_pos2nt_idx[sampos], ",");

			cout << endl;
		}
	}
}

/*
 * parse the CIGAR string and store the read counts, aligned bases, and nt positions at each genomic position
 *	i.e. if we just want to get the error rates, we can skip this
 *
 * TO DO:
 * 	stream line this function
 */
void sam::parseCigar(const string & readseq, const string & cigar, map <int, int> &sam_pos2cnt, map <int, vector <char> > &sam_pos2pile, const int &offset){
	string run("");
	int startpos= 0;
	int runval;
	int j;
	int sam_pos;
	int last_readpos= 0;

	for (int i=0; i< (int)cigar.length(); i++){
		// is alphabet?
		if (isalpha(cigar.at(i))){
			// cout << cigar.at(i) << " -> test \n";
			runval= atoi(run.c_str());
			switch (cigar.at(i)){
				// TO DO:
				// 	add code to deal with INDELs here
				//
				case 'M':
					// mark on chromosome
					// cout << "\tM run= "<< run << endl;
					// cout << "\t\tmarking "<< startpos << " - " << (startpos+ runval -1)<< endl;
					for (j=startpos ; j < (startpos+ runval); j++){
						sam_pos2cnt[j + offset]++;
					}

					if (b_pileup){
						for (j= 0; j < runval; j++){
							sam_pos= j+startpos+offset;
							// cout << "adding " << readseq.at(last_readpos+j) << " to position "<< (j+ startpos+ offset) << endl;
							sam_pos2pile[sam_pos].push_back(readseq.at(last_readpos+j));
						}
					}
					// TO DO:
					// add read quality and mapping quality

					// position in read, i.e. sam_pos2readidx
					if (b_ntidx){
						for (j= 0; j < runval; j++){
							sam_pos= j+startpos+offset;
							sam_pos2nt_idx[sam_pos].push_back(last_readpos+j);
						}
					}

					// TO DO:
					// 	add option to store the read ids
					startpos+= runval;
					last_readpos+= runval;
					run= "";
					break;
				case 'N':
					// skip
					// cout << "\tN run= "<< run << endl;
					// cout << "\t\tskipping "<< startpos << " - " << (startpos+ runval -1)<< endl;
					startpos= startpos+ runval;
					run= "";
					break;
				default:
					cerr << "Not recognized: "<< cigar.at(i) << endl;
			}
		}else{
			// number
			// cout << cigar.at(i) << endl;
			run+= cigar.at(i);
		}
	}
}

/*
 * parse MD tags, estimate error rates, and reconstruct the reference bases
 * 	TO DO:
 * 		check if mdtag starts wtih MD:Z
 * 		reconstruct the reference bases in order to calculate position specific error rate
 */
void sam::parseMDtag(const string & mdtag){
	// MD tag:
	// 	e.g. MD:Z:76 => all perfectly matched
	// 	e.g. MD:Z:2T72 => 2 matching bases -> T -> 72 matching bases

	string run;
	int runval=0;
	int current_pos=0;

	for (int i=0; i< (int) mdtag.length(); i++){
		// is alphabet? i.e. found a mismatch
		if ( isalpha(mdtag.at(i)) ){
			runval= atoi(run.c_str());
			if ( mdtag.at(i) == 'A' || mdtag.at(i) == 'T' || mdtag.at(i) == 'C' || mdtag.at(i) == 'G' ){
				current_pos+= runval;
				// cout << mdtag << " found mismatch " << mdtag.at(i) << " at position "<< current_pos<<  endl;
				mismatch_cnt[current_pos]++;
				current_pos++;
				run= "";

			}else{
				if ( mdtag.at(i) != 'N' )
					cerr << "MDtag not recognized: "<< mdtag.at(i) << "\t"<< mdtag << endl;
			}
		}else{
			// number
			// cout << cigar.at(i) << endl;
			run+= mdtag.at(i);
		}
	}
}

/*
 * get mismatch counts
 */
vector <int> sam::get_mismatch_counts(){
	vector <int> mmrate;
	for (int i=0; i< read_len; i++){
		mmrate.push_back(mismatch_cnt[i]);
	}
	return mmrate;
}

void sam::show_mismatch_rate(){
	cout << "position\tcount\tpercentage\ttotal_num_reads\n";
	for (int i=0; i< read_len; i++){
		cout << i << "\t"<< mismatch_cnt[i] << "\t" << double(mismatch_cnt[i])/num_reads<< "\t"<< num_reads<< endl;
	}
}

/*
 * store the mismatch rates per read position
 *
 */
void sam::save_mismatch_rate(const string &fname){
	ofstream fout(fname.c_str());
	if (! fout.is_open() ){
		cerr << "Error saving to "<< fname << endl;
	}else{
		fout << "# position\tcount\tpercentage\ttotal_num_reads\n";
		for (int i=0; i< read_len; i++){
			fout << i << "\t"<< mismatch_cnt[i] << "\t" << double(mismatch_cnt[i])/num_reads<< "\t"<< num_reads<< endl;
		}
	}
	fout.close();
}

/*
 * automatically determine the read length by going through the first 200 lines in a SAM file and finding the longest sequence length
 */
int sam::getReadLength(const string & samfname){
	vector <string> tokens;

	int idx;
	int len=0, tmplen=0;
	string line;

	fileSlurp fs(samfname, MB);
	while(! fs.eof()){
		line= fs.readLine();
		if ( startsWith(line, "@") ){
			continue;
		}

		Tokenize_strict(line, tokens, "\t");

		tmplen= tokens[9].length();

		if (tmplen > len )
			len= tmplen;
		idx++;
		if (idx == 200){
			break;
		}
	}
	fs.close();
	return len;
}

/*
 * load sam file
 *
 * TO DO:
 * 	maybe specify which chromosome we want to get read information for
 */
void sam::load(const string samfname, string tgChr){
	// TO DO:
	// 	assuming the SAM file has @... header
	cerr << "opening SAM file "<< samfname << "\nunique= "<< b_unique <<endl;

	reset();	// clean up previously loaded sam data (sam_pos2cnt, sam_pos2pile)

	vector <string> tokens;
	string line;
	string cigar;
	string mdtag;

	vector <int> pos2mark;
	map <int, map <int, bool > > readStartPos;
	int strand;
	int sam_mappos;
	int bitflag;
	bool bSpecifiedChr= false;		// if false -> assume all come from same chromosome
	if (tgChr != "")
		bSpecifiedChr= true;		// if true -> only retain those matching this chromosome
						//

	// TO DO:
	// 	remember not to overwrite the read length every time a different chromosome is loaded
	read_len= getReadLength(samfname);
	mismatch_cnt.assign(read_len, 0);		// initialize mismatch counts
	cerr << "read length= "<< read_len << endl;

	//bool bErrOnly= true;			// get err rates only (i.e. won't get pileup, read counts, nt indexes)
	fileSlurp fs(samfname);
	while(! fs.eof()){
		line= fs.readLine();
		if ( startsWith(line, "@") ){
			continue;
		}

		Tokenize_strict(line, tokens, "\t");

		if (bSpecifiedChr){
			// this might be a whole genome sam file
			// we want to look at only tgChr 
			// TO DO:
			// 	=> WILL EVENTUALLY USE SORTED SAM FILE ONLY
			if (tgChr != tokens[2]){
				// cerr << "not "<< tgChr << endl;
				continue;
			}
		}

		//	get strand and position
		bitflag = atoi(tokens[1].c_str());
		strand= (bitflag & ( 1 << 4) || bitflag & (1 << 5));

		// convert to genomic positions (sam, 1-based) and store
		sam_mappos= atoi(tokens[3].c_str());

		// check if this read is unique (based on strand and position)
		if (b_unique ){
			if (readStartPos.count(sam_mappos) == 0 || readStartPos[sam_mappos].count(strand) == 0){
				readStartPos[sam_mappos][strand]= true;
				// cout << "adding " << tokens[0] << " " << sam_mappos << " strand= "<< strand << endl;
			}else{
				// seen this position already, skip
				// cout << "seen " << tokens[0] << " " << sam_mappos << " strand= "<< strand << " ... skipped ...\n";
				continue;
			}
		}

		// parse CIGAR string
		cigar= tokens[5];

		// get MD tag
		// 	parse md tag to reconstruct the reference bases
		mdtag= tokens[tokens.size()-1];
		if (mdtag.substr(0, 4) != "MD:Z"){
			// it might not always be the last field. double check!
			cerr << "Couldn't find an MD tag!!! \n";
		}else{
			parseMDtag(stringToUpper(mdtag.substr(5)));
		}


		// cout << tokens[2] << "\t" << tokens[3] << "\t"<< cigar << "\t"<< mdtag << endl;

		// TO DO:
		// 	now we assume everything is on the same chromosome...
		// 	modify to generalize later
		// if we only want to get the error rates, skip calculating the read counts, etc.
		if ( !b_err_only )
			parseCigar(tokens[9], cigar, sam_pos2cnt, sam_pos2pile, sam_mappos);

		num_reads++;
	}
	fs.close();

	// here it's 1-based coordinate as reported by *.sam
	// TO DO:
	// 	continue here
	// 	1. calculate RPKM for a gene (geneModel) from SAM file
	// 	2. should be able to clean up sam_pos2cnt when switching to a different chromosome !!
	// cout << "52549466 => " << sam_pos2cnt[52549466] << endl;
	deleteMapComponents <int, map <int, bool> > (readStartPos);	// remove map components
}

/*
 * get read counts for a gene using UCSC information and SAM alignment file (total, 5'UTR, coding, 3'UTR)
 *
 * return value:
 * 	map <string, int>
 *
 * data structure:
 * 	map[utr5]= 5'utr count
 * 	map[utr3]= 3'utr count
 * 	map[coding]= coding count
 * 	map[total]= total count
 *	map[utr5_len]= 5'utr length
 *	map[utr3_len]= 3'utr length
 *	map[coding_len]= coding length
 *	map[total_len]= total length
 *
 *	NOTE:
 *		UCSC: start= 0-based end= 1-based
 *		SAM:  all 1-based
 */
map <string, int> sam::getGeneReadCount(string gene, genemodel &gm){
	// get gene structure
	// 	5'UTR, coding, 3'UTR (based on strand info)
	gene_str tmpgene= gm.getGeneStr(gene);
	if (tmpgene.isNULL()){
		cerr << "WARNING: no such gene "<< gene << endl;
		return map <string, int> ();			// return empty map (TO DO: a better way??)
	}

	// TO DO:
	// 	make sure chromosome matches record

	if (verbose){
		cout << "transcript: "<< tmpgene.getRepTranscript() << endl;
		cout << "strand: "<< tmpgene.strand << endl;
		cout << "tx start: "<< tmpgene.start << endl;
		cout << "cd start: "<< tmpgene.cdst << endl;
	}
	int cnt, utrup= 0, utrdown= 0, body= 0;
	int utr5, utr3;
	int utrup_len=0, utrdown_len=0, total_len= 0;
	int utr5_len, utr3_len;

	// iterate through the exons
	for (int i=0; i< (int)tmpgene.exst.size(); i++){
		if (verbose)
			printf("Exon %d: %d - %d\n", i, tmpgene.exst[i], tmpgene.exen[i]);

		// walk through the exons (in 0-based coordinates)
		for (int j= tmpgene.exst[i] ; j < tmpgene.exen[i]; j++){
			if (sam_pos2cnt.count(j+1)> 0){							// SAM is 1-based
				cnt= sam_pos2cnt[j+1];
				// cout << "sam pos: "<< (j+1) << ": "<< cnt << "\t"<< vector2string2 <char> (sam_pos2pile[j+1], "") << endl;
				if (j < tmpgene.cdst){
					// UTR 1
					utrup+= cnt;
				}else if ( j >= tmpgene.cden ){
					// UTR 2, i.e. j > cden-1
					utrdown+= cnt;
				}else{
					body+= cnt;
				}
			}
		}

		// calculate 5'utr, coding, 3'utr lengths
		if (tmpgene.exst[i] < tmpgene.cdst ){
			//
			utrup_len+= ( min( tmpgene.exen[i]- 1 , tmpgene.cdst ) - tmpgene.exst[i]);	// start= 0-based
		}
		if (tmpgene.exen[i] > tmpgene.cden ){
			utrdown_len+= ( tmpgene.exen[i] - max ( tmpgene.cden, tmpgene.exst[i] + 1) );	// end= 1-based
		}
		total_len+= (tmpgene.exen[i] - tmpgene.exst[i] ); // i.e. (tmpgene.exen[i] -1 - tmpgene.exst[i]  + 1 );
	}

	if (verbose){
		cout << "cd end: "<< tmpgene.cden << endl;
		cout << "tx end: "<< tmpgene.end << endl;
	}

	if (tmpgene.strand== '+'){
		utr5= utrup;
		utr3= utrdown;
		utr5_len= utrup_len;
		utr3_len= utrdown_len;
	}else if (tmpgene.strand == '-'){
		utr5= utrdown;
		utr3= utrup;
		utr5_len= utrdown_len;
		utr3_len= utrup_len;
	}else{
		cerr << "ERROR: unknown strand value "<< endl;
	}

	// prepare output
	if (verbose)
		cout << "5'UTR: " << utr5 << " (len= "<< utr5_len << ") coding: "<< body << " (len= "<< (total_len - utr5_len- utr3_len) <<
			") 3'UTR: "<< utr3 << " (len= "<< utr3_len<< ")"<< endl;
	map <string, int> cnt_data;
		cnt_data["utr5"]= utr5;
		cnt_data["utr3"]= utr3;
		cnt_data["coding"]= body;
		cnt_data["total"]= (utr5+ utr3 + body);
		cnt_data["utr5_len"]= utr5_len;
		cnt_data["utr3_len"]= utr3_len;
		cnt_data["coding_len"]= (total_len - utr5_len - utr3_len);
		cnt_data["total_len"]= total_len;
	return cnt_data;
}

map <string, map<int,int> > sam::loadAllChrSAM(const string &samdir, const string &prefix, const string &suffix){
	map <string, map<int,int> > output;
	HG18 hg18;
	vector <string> chrs= hg18.getChrs();
	for (int i=0; i< (int)chrs.size(); i++){
		if (chrs[i].find("_") != string::npos){
			// if (chrs[i] != "chr1"){
			// found a _random
			continue;
		}
		// cerr << "WARNING: loaded only chr1's sam file\n";

		string readfile= samdir+ "/"+ prefix+ chrs[i]+ suffix;
		cerr << "loading "<< chrs[i]<< " from "<< readfile<< endl;
		sam a_sam;
		a_sam.load(readfile);
		output[chrs[i]]= a_sam.getPos2Cnt();		// TO DO: this is a copy of pos2cnt or a reference to??
	}
	return output;
}

/*
 * get machine name and column index of X0:i
 */
void sam::findMachineNameX0 (const string &fname, string & machine_code, int & mapcnt_idx){
	// int mapcnt_idx= -1;			// which column contains the number of best hits
	string line;
	vector <string> tokens, fields;
	fileSlurp tmpfs(fname, MB);
	while ( !tmpfs.eof() ){
		line= tmpfs.readLine();
		if (startsWith(line, "@"))
			continue;

		Tokenize(line, tokens, "\t");

		if (tokens[2] == "*")
			// nothing matched
			continue;

		Tokenize(tokens[0], fields, ":");

		// update machine name
		machine_code= fields[0];
		// machine_code_len= machine_code.length()+1;			// machine instrument code length +1

		for (int i= 11; i< (int) tokens.size() ; i++){
			if (startsWith(tokens[i], "X0:i:")){
				mapcnt_idx= i;
				// cerr << "mapcnt_idx= "<< mapcnt_idx<< "\t" << tokens[i]<<  endl;
				break;
			}
		}
		if (mapcnt_idx != -1)
			break;
	}
	tmpfs.close();
}

/*
 * load the aligned reads and return a read_str object
 *
 * 	subset: load only reads in this subset (i.e. for paired end reads, load reads that are aligned in both)
 * 	bUniqMap: load only the unambiguously mapped reads
 *
 *  example alignment foramt:
 *  HWI-EAS261_0009:1:1:7:265#0     16      chr10   122042710       25      76M     *       0       0
 */
map <string, read_str> sam::loadAlignedReads(const string &fname, map <string, map < int, int > > & counts, map <string, read_str> & subset, bool bUniqMap, int & read_len){
	cerr << "opening "<< fname << endl;
	// int read_len;
	read_len= getReadLength(fname);
	cerr << "read length= "<< read_len << endl;

	// cerr << "assuming all reads are 76M ... fix later... (now we don't calculate the density)" << endl;
	int linecnt= 0;

	string machine_code;							// Illumina instrument code, e.g. "HWI-EAS261_0009:";
	int machine_code_len;							// Illumina instrument code len, i.e. machine_code.length();

	map <string, read_str> id2readstr;
	string readid;

	double vm, rss;

	int bitflag;
	char strand;

	string line;
	vector <string> tokens;
	vector <string> fields;

	// find X0 tag column number using the first 1000 lines
	//
	int mapcnt_idx= -1;			// which column contains the number of best hits
	fileSlurp tmpfs(fname, MB);
	while ( !tmpfs.eof() ){
		line= tmpfs.readLine();
		if (line.at(0)== '@')
			continue;

		Tokenize(line, tokens, "\t");

		if (tokens[2] == "*")
			// nothing matched
			continue;

		Tokenize(tokens[0], fields, ":");

		// update machine name
		machine_code= fields[0];
		machine_code_len= machine_code.length()+1;			// machine instrument code length +1

		for (int i= 11; i< (int) tokens.size() ; i++){
			if (startsWith(tokens[i], "X0:i:")){
				mapcnt_idx= i;
				cerr << "mapcnt_idx= "<< mapcnt_idx<< "\t" << tokens[i]<<  endl;
				break;
			}
		}
		if (mapcnt_idx != -1)
			break;
	}
	tmpfs.close();
	cerr << "machine_code: "<< machine_code << endl;
	cerr << "machine_code_len: "<< machine_code_len<< endl;

	// TO DO:
	// 	maybe add option to rewind in fileSlurp
	// open read file using 50MB of buffer
	fileSlurp fs(fname, 50* MB);
	while(! fs.eof()){
		linecnt++;
		if (linecnt % 10000 ==0){
			process_mem_usage(vm, rss);
			vm/=1024;
			rss/=1024;
			cerr << "read "<< linecnt << "  lines, VM: "<< vm << "mb; RSS: "<< rss << "mb     \r";
		}

		line= fs.readLine();
		if (line.at(0)== '@'){
			continue;
		}
		// get read ids; here we assume paired-end reads don't have the trailing 1/2 tag
		// 	e.g. HWI-EAS261_0009:1:1:0:636#0
		// TO DO:
		// 	1. simplify the read name to reduce memory usage
		Tokenize(line, tokens, "\t");

		if (tokens[2] == "*"){
			// nothing matched
			continue;
		}else{
			// load subset only?
			readid= tokens[0].substr(machine_code_len);	// TO DO: convert to numeric read structure
			if (subset.size()>0 && subset.count(readid)==0){
				// subset provided. this read not in subset
				continue;
			}

			bitflag = atoi(tokens[1].c_str());
			strand= (bitflag & ( 1 << 4) || bitflag & (1 << 5))?'-':'+';	// 0: forward 1: reverse

			// TO DO:
			// 	filter out multiply aligned reads
			// 	BWA: look for tag X0 (number of best hits)
			// 		e.g. X0:i:20
			// 	Bowtie:
			// 		the XM:I field is set to <int>+1.

			if (bUniqMap && tokens[mapcnt_idx] != "X0:i:1"){
				// read not uniquely aligned
				// cout << readid << " is not uniquely mapped: "<< tokens[mapcnt_idx]<< endl;
				continue;
			}

			read_str rs;
				rs.chr= tokens[2];		// TO DO:	use numeric idx to save space
				rs.strand= strand;
				rs.start= atoi(tokens[3].c_str());
				rs.end= rs.start+ read_len ;		// TO DO:
								// here we assume ALL reads align without gaps to the genome (i.e. 76M)
								// 	-> fix this later (parse CIGAR to get actual span on genome)
			id2readstr[readid]= rs;

			// make peak density plot
			// 	-> moved to appBedToDensity()
			// TO DO:
			// 	seems very memory intensive... skipped for now...
			/*
				for (int i= rs.start; i< rs.end; i++){
					counts[rs.chr][i]++;
				}
			*/
		}
	}
	fs.close();
	return id2readstr;
}



/*
 * TO DO:
 * 	consider cases where two or more genes overlap...
 *
 * TO DO:
 * 	seems very slow...
 * 	... maybe looking up read counts from indexed pileip file is faster
 */
void pileup::getRPKM(string pilefname, genemodel &gm){
	string line;
	// ifstream myfile(pilefname.c_str());

	fileSlurp fh(pilefname);

	string chr= "";
	string oldchr= "";
	int pos;

	// TO DO:
	// 	maybe use fileSlurp() to speed up
	//
	if (fh.is_open()){
		int cnt= 0;
		vector <string> tokens;
		set <string> chr_gene_set;
		gene_str	current_gene;
		string 		current_gene_name;
		string		old_gene_name;
		map <string, int> rpkm;
		int binid;
		int oldbinid;
		int binsize= gm.getBinSize();

		while (! fh.eof()){
			cnt++;
			// getline (myfile,line);
			line= fh.readLine();
			tokens.clear();
			Tokenize(line, tokens, "\t");

			chr= tokens[0];
			pos= atoi(tokens[1].c_str());
			binid= int(pos/binsize);

			// cout << "chr= "<< chr << " -> "<< chr.compare("klasdjf") << endl;;
			if (binid != oldbinid){
				// unseen chromosome -> load gene starts
				// chr_gene_set= gm.getGenesOnChr(chr);
				chr_gene_set= gm.getGenesInBin(chr, binid);

				// printf("loading chromosome genes for %s bin %d\n", chr.c_str(), binid);
			}
			oldchr= chr;
			oldbinid= binid;

			// look for gene that overlaps with this position
			current_gene_name= gm.findFirstGeneAtPosFromSet(chr_gene_set, pos);

			if (current_gene_name == ""){
				// no gene here -> continue
				continue;
			}


			// get gene structure if current_gene_name is different from previous gene
			if (current_gene_name != old_gene_name){
				// cout << "loading gene structure for "<< current_gene_name << endl;
				current_gene= gm.getGeneStr(current_gene_name);
			}

			// TO DO:
			// 	add in3UTR() and in5UTR()
			// 	-> take into consideration the strand
			if (current_gene.inExon(pos)){
				// calculate exon rpkm
				rpkm[current_gene_name]+= atoi(tokens[3].c_str());
				//  cout << chr << " " << pos << " is in exon of "<< current_gene_name << " -> rpkm = "<< rpkm[current_gene_name]<< endl;
			}

			old_gene_name= current_gene_name;

			// load read information here

			// printf("%d elements\n", tokens.size());
			// printf("%s\t%s\t%s\n", tokens.at(0).c_str(), tokens.at(1).c_str(), tokens.at(2).c_str());
			// cout << line << endl;
		}
		fh.close();
		cout << "Total: "<< cnt << " lines\n";

		// continue here !!
		// TO DO:
		// return RPKM as an object
	}
	else cout << "Unable to open file " << pilefname << endl;
}

/*
 * TO DO:
 * retrieves a region of the genome from the pileup file
 * 	1. make sure block_size is consistent
 * 	2. maybe store block_size in cidx file
 *
 * */
int pileup::fetchReadCount(string chr, int sampos){
	// check if it's in cache -> return value if yes
	if (chr != cachedChr){
		// TO DO:
		// 	c++ smart enough to get rid of all content maps???
		cerr << "new chr ("<< chr << ") -> cleaning up cache...\n";
		cachedReads.clear();
	}

	// cout << "looking for "<< chr << " " << sampos << endl;
	if (cachedReads[chr].count(sampos)> 0){
		// cout << "found "<< chr << " " << sampos << " => "<< cachedReads[chr][sampos]<< endl;
		return cachedReads[chr][sampos];
	}

	int idx= int(sampos/block_size)*block_size;
	int ans= -1;

	// cout << "opening "<< pile_fname << endl;
	string line;

	// use global gTokens to speed up
	int tmppos;
	int tmpread;
	string tmpchr;

	// open pileup file
	// 	moved open() to loadIndex() so we don't open and close all the time

	if (allpos2idx[chr].count(idx) > 0){
		fh_pile.seekg( allpos2idx[chr][idx]);
		// continue here
		//	1. load/cache N positions
		//	2. make sure position is converted to 1-based
		//	3. option to turn caching on/off
		//	4. find out where this block ends and SLURP in a huge block to reduce I/O

		// cout << "loading block "<< chr << " - "<< idx << " size= "<< allpos2blocksize[chr][idx]<< endl;
		fh_pile.read(buffer, allpos2blocksize[chr][idx]);
		buffer[fh_pile.gcount()]= '\0';
		buffer_str.assign(buffer);

		// getline (fh_pile,line);
		gTokens.clear();
		Tokenize(buffer_str, gTokens, "\n");
		for (int i=0; i< (int)gTokens.size(); i++){
			gTmpTokens.clear();
			Tokenize(gTokens[i], gTmpTokens, "\t");
			tmpchr= gTmpTokens[0];
			tmppos= atoi(gTmpTokens[1].c_str());
			tmpread= atoi(gTmpTokens[3].c_str());

			// caching goes here
			cachedReads[tmpchr][tmppos]= tmpread;

			// cout << "line "<< i<< " " << gTokens[i] << endl;
		}
		// TO DO:
		// 	determine when to clean up cache
		// 		-> when we change to another chromosome?
		cachedChr= chr;
		if (cachedReads[chr].count(sampos)> 0){
			// cout << "scanning: found "<< chr << " " << sampos << " => "<< cachedReads[chr][sampos]<< endl;
			ans= cachedReads[chr][sampos];
		}

		// cout << "record begins at "<< chr << " " << idx << endl;
		// cout << line << endl;
	}
	return ans;
}


/*
 * TO DO:
 * load indexed pileup file
 */
void pileup::loadIndex(string pilefname){
	// load index file in one go
	vector <string> tokens;
	string line;
	int pos, oldpos, tmpblock;
	streamoff fidx, oldfidx;
	string chr, oldchr;

	// make name for index file
	string index_fname= pilefname+".cidx";
	pile_fname= pilefname;		// TO DO:
					// pilefname= local; pile_fname= global => rename later

	// check for pileup file existence
	fh_pile.open(pile_fname.c_str());
	if (fh_pile.bad())
	{
		cerr << "ERROR! couldn't open "<< pile_fname << endl;
		exit(1);
	}

	// check for index file existence
	cout << "loading "<< index_fname << endl;

	fileSlurp fh(index_fname);
	if (! fh.is_open() ){
		cout << "Error opening "<< index_fname << endl;
		return;
	}
	// get first record (e.g. chr10, 0, 0)
	line= fh.readLine();
	tokens.clear();
	Tokenize(line, tokens, "\t");
	oldchr= tokens[0];
	oldpos= atoi(tokens[1].c_str());
	oldfidx= atol(tokens[2].c_str());

	while (! fh.eof()){
		line= fh.readLine();
		tokens.clear();
		Tokenize(line, tokens, "\t");

		chr= tokens[0];
		pos= atoi(tokens[1].c_str());
		fidx= atol(tokens[2].c_str());

		tmpblock= (fidx- oldfidx);

		allpos2idx[oldchr][oldpos]= oldfidx;
		allpos2blocksize[oldchr][oldpos]= (tmpblock);

		if (tmpblock> maxblocksize)
			maxblocksize= tmpblock;

		oldchr= chr;
		oldpos= pos;
		oldfidx= fidx;
		// cout << tokens[0] << "\t" << pos << " -> " << allpos[chr][pos] << endl;

		// cout << fh.readLine() << endl;
	}
	fh.close();

	cout << "init. file chunk buffer.. size= "<< maxblocksize << endl;
	// free(buffer);
	buffer= new char [maxblocksize+1];		// remember to free
	buffer_str.reserve(maxblocksize);

}

/*
 * TO DO:
 * index the pileup file to speed up I/O and gene retrieveal
 */
void pileup::index(string fname){
	string line;
	string chr;
	int    pos;
	int    idxid;
	int    oldidxid= -1;
	ifstream myfile(fname.c_str());
	cout << "opening pileup file "<< fname << endl;

	// position in file for a chromosomal position (block size= 100000)
	// map <string, map<int, pos_type> > posinfile;
	streamoff fpointer;

	// open output index file
	string outfile= fname+".cidx";
	cout << "writing to index file file "<< outfile << endl;
	cout << "block size= "<< block_size << " bp" << endl;

	ofstream fout (outfile.c_str());
	if (! fout.is_open()){
		cerr << "Unable to open file" << endl;
		exit (1);
	}

	if (! myfile.is_open()){
		cerr << "Unable to open file "<< fname << endl;
		exit (1);
	}

	vector <string> tokens;

	myfile.seekg(0, ios::end);
	streamoff eof_pos= myfile.tellg();
	myfile.seekg(0, ios::beg);
	// cout << "eof pos= " << eof_pos << endl;

	while (! myfile.eof() ){
		fpointer= myfile.tellg();
		getline (myfile, line);

		tokens.clear();
		Tokenize(line, tokens, "\t");

		chr= tokens[0];
		pos= atoi(tokens[1].c_str());
		idxid= int(pos/block_size)*block_size;

		// if (posinfile.count(chr)==0 || posinfile[chr].count(idxid)==0)
		if (idxid != oldidxid){
			// posinfile[chr][idxid]= fpointer;
			// printf("adding %s %d %lu -> stored= %lu\n", chr.c_str(), idxid, fpointer, posinfile[chr][idxid]);
			fout << chr << "\t" << idxid << "\t" << fpointer << "\n";
		}
		oldidxid= idxid;
	}

	// record last position in file
	//
	// posinfile[chr][idxid]= myfile.tellg();
	// printf("adding %s %d %lu -> stored= %lu\n", chr, idxid+block_size, myfile.tellg(), posinfile[chr][idxid]);
	fout << chr << "\t" << (idxid+block_size) << "\t" << eof_pos << endl;
	fout.close();

}

/*
 * TO DO:
 * 	continue here
 * 		goals:
 * 			1. index pileup file
 * 			2. use file slurping to speed up
 * 			3. get expression using indexing
 */
void pileup::load(string fname){
	string line;
	ifstream myfile(fname.c_str());
	if (myfile.is_open())
	{
		int cnt= 0;
		vector <string> tokens;
		while (! myfile.eof() )
		{
			getline (myfile,line);
			tokens.clear();
			Tokenize(line, tokens, "\t");

			// load read information here

			// printf("%d elements\n", tokens.size());
			// printf("%s\t%s\t%s\n", tokens.at(0).c_str(), tokens.at(1).c_str(), tokens.at(2).c_str());
			cnt++;
			// cout << line << endl;
		}
		myfile.close();
		cout << "Total: "<< cnt << " lines\n";
	}
	else cout << "Unable to open file";
}




/*
 * generate an interaction plot for paired-end reads
 *
 * option to output only the simultaneously aligned reads
 *
 * 	bUniqRead: remove PCR duplicates
 * 	bUniqMap:  remove ambiguously mapped reads
 *
 * TO DO:
 * 	1. option to remove low-quality reads
 *
 * format:
 * 	track name=pairedReads description="Clone Paired Reads" useScore=1
 * 	chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
 * 	chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
 *
 */
void gbrowser::appPairedEndPlot(string file1, string file2, string outfile, int max_dist, bool bUniqRead){
	// read file1
	// 	store all aligned reads (id, position, strand)
	// read file2
	// 	load only those in file1
	//
	string line;
	string machine_code= "HWI-EAS261_0009:";

	// load the alignments (read positions, strand)
	cerr << "max_dist= "<< max_dist<< endl;
	cerr << "bUniqRead= "<< bUniqRead << endl;

	bool bUniqMap= true;			// filter out ambiguously mapped reads

	// to record the start position and strand of a read
	// 	chr -> pos -> strand (binarized)
	map <string, map <int, int> > readRegistry;
	const int bit_plus= 1;			// 01
	const int bit_minus= (1<< 1);		// 10
	const int bit_both= 3;			// 11

	// set up default parameters for bed graph output
	string id;
	int dist= -1, score= 1000;
	string chrname, name, color, blockSizes, blockStarts;
	int chrStart, chrEnd, thickStart, thickEnd;
	char strand;
	int blockCount;

	int read_len;
	map <string, map <int, int> > counts;
	map <string, read_str> null_reads;			// null read set for compatibility
	read_str rs1, rs2;					// structure to hold read information

	cerr << "loading SAM file 1"<< endl;
	map <string, read_str> id2readstr1= sam::loadAlignedReads(file1, counts, null_reads, bUniqMap, read_len);
	cerr << "loaded "<< id2readstr1.size() << " uniquely mapped reads\n";

	cerr << "loading SAM file 2"<< endl;
	map <string, read_str> id2readstr2= sam::loadAlignedReads(file2, counts, id2readstr1, bUniqMap, read_len);
	cerr << "loaded "<< id2readstr2.size() << " uniquely mapped reads\n";

	ofstream fout(outfile.c_str());
	if (! fout.is_open()){
		cerr << "Error opening "<< outfile << endl;
		exit(1);
	}

	// header is a problem for bedToBigBed, so now gotten rid of
	// fout << "track name=pairedReads description=\"chia_pet data\""<< endl;

	// constant stuff in the ucsc browser track
	score= 960;
	color= "255,0,0";
	blockCount= 2;
	blockSizes= castToString <int> (read_len) + ","+ castToString <int> (read_len);    // e.g. "76,76";

	// use iterator to traverse the keys
	map <string, read_str> :: iterator it;
	for (it= id2readstr1.begin(); it!= id2readstr1.end(); it++){
		id= (*it).first;
		rs1=  id2readstr1[id];
		rs2=  id2readstr2[id];

		if (rs1.chr == rs2.chr){
			// on the same chromosome
			dist = abs(rs1.start- rs2.start);
		}else{
			// not aligned to the same chromosome -> skipped
			dist = -1;
			continue;
		}

		// only output those that are <100Kbp apart
		if (dist > max_dist){
			// cout << "too far between two ends in "<< id<< " dist= "<< dist << endl;
			continue;
		}

		// continue processing the read pair
		// 	note: start: 0-based, end= 1-based
		// 	strand based on the 5' read
		chrname= rs1.chr;
		name= id;

		if (rs1.start < rs2.start){
			// rs1 upstream of rs2
			chrStart= rs1.start;
			chrEnd= rs2.end;
			thickStart= rs1.start;
			thickEnd= rs2.end;
			strand= rs1.strand;
			blockStarts= "0,";
			blockStarts+= castToString <int> (rs2.start-rs1.start);
		}else{
			chrStart= rs2.start;
			chrEnd= rs1.end;
			thickStart= rs2.start;
			thickEnd= rs1.end;
			strand= rs2.strand;
			blockStarts= "0,";
			blockStarts+= castToString <int> (rs1.start-rs2.start);
		}

		// get rid of PCR duplicates (bUniqRead option)
		// currently we use the 5' coordinate and orientation to mark duplicates (c.f. MarkDuplicates@Picard)
		// 	TO DO: maybe check both start and end are the same?
		if (bUniqRead){
			// chr -> position -> strand (00 01 10 11)
			if (readRegistry[chrname].count(chrStart)== 0){
				// not at this position before
				// use rs1.strand as a reference when determining unique reads
				if (strand == '-')
					readRegistry[chrname][chrStart]= bit_minus;		// 10	-
				else
					readRegistry[chrname][chrStart]= bit_plus;		// 01	+
			}else{
				// seen this positon before, check if strand is set already
				if (readRegistry[chrname][chrStart]== bit_plus && strand == '+'){
					continue;
				}else if (readRegistry[chrname][chrStart]== bit_minus && strand == '-'){
					continue;
				}else if (readRegistry[chrname][chrStart]== bit_both){
					continue;
				}else{
					// we'll get here if 10 and +, 01 and -, or 00
					if (strand== '+')
						readRegistry[chrname][chrStart] |= bit_plus;
					else if (strand == '-')
						readRegistry[chrname][chrStart] |= bit_minus;
				}
			}
		}

		fout << chrname << " "<< chrStart << " " << chrEnd << " " << name << " " << score << " " << strand << " " <<
				thickStart << " " << thickEnd << " " << color << " " << blockCount << " "<< blockSizes <<
				" " << blockStarts<< endl;
	}
	fout.close();

	cerr << "finished generating "<< outfile << endl;
	exit(0);
}

/*
 * output a bedGraph file from a dense read count object (vector <int> allpos)
 *
 */
void gbrowser::appendToBedGraph(ofstream &fout, string chr, int chrlen, vector <int> allpos){
	int cnt= -1;
	int oldcnt= 0;
	int start_pos= 0;

	// stringstream ss;

	for (int i=0; i< chrlen; i++){
		cnt= allpos[i];
		// cout << i << " " << cnt << endl;
		if (cnt != oldcnt || i== chrlen-1){
			if (start_pos != i)
				fout << chr << "\t"<< start_pos << "\t" << i << "\t" << oldcnt <<  endl;
			start_pos= i;
		}
		oldcnt= cnt;
	}
	// return ss.str();
}

/*
 * output a bedGraph file from a relatively sparse read count object (map <int, int> allpos)
 *
 * TO DO:
 * 	here we scan across all positions although the read counts are sparse.
 * 	should use sorted keys(read cnt) to speed up
 *
 */
void gbrowser::appendToBedGraph(ofstream &fout, string chr, int chrlen, map <int, int> allpos){
	int cnt= -1;
	int oldcnt= 0;
	int start_pos= 0;

	// stringstream ss;

	for (int i=0; i< chrlen; i++){
		if (allpos.count(i) > 0){
			cnt= allpos[i];
		}else{
			cnt= 0;
		}
		// cout << i << " " << cnt << endl;
		if (cnt != oldcnt || i== chrlen-1){
			if (start_pos != i)
				fout << chr << "\t"<< start_pos << "\t" << i << "\t" << oldcnt <<  endl;
			start_pos= i;
		}
		oldcnt= cnt;
	}
	// return ss.str();
}

/*
 * output a bedGraph file from a sparse read count object (map <int, int> allpos)
 *
 */
void gbrowser::appendToBedGraph(ofstream &fout, string chr, map <int, int> allpos){
	/*nt cnt= -1;
	int oldcnt= 0;
	int start_pos= 0;
	*/

	vector <int> sorted_pos= getKeys <int, int> (allpos);
	sort(sorted_pos.begin(), sorted_pos.end());

	for (int i=0; i< (int)sorted_pos.size(); i++){
		fout << chr << "\t"<< sorted_pos[i]<< "\t"<< (sorted_pos[i]+1)<< "\t"<< allpos[sorted_pos[i]]<< endl;
	}
}


string gbrowser::countToBedGraph(string chr, int chrlen, vector <int> &allpos){
	int cnt= -1;
	int oldcnt= 0;
	int start_pos= 0;

	stringstream ss;

	for (int i=0; i< chrlen; i++){
		cnt= allpos[i];
		// cout << i << " " << cnt << endl;
		if (cnt != oldcnt || i== chrlen-1){
			if (start_pos != i)
				ss << chr << "\t"<< start_pos << "\t" << i << "\t" << oldcnt <<  endl;
			start_pos= i;
		}
		oldcnt= cnt;
	}
	return ss.str();
}



/*
 * output bed graph
 *
 * option:
 * 	call bedTo
 *
 * TO DO:
 * 	add header
 * 	track type=bedGraph name="test - read coverage"
 */
string gbrowser::countToBedGraph(string chr, int chrlen, map <int, int> allpos){
	int cnt= -1;
	int oldcnt= 0;
	int start_pos= 0;

	stringstream ss;

	for (int i=0; i< chrlen; i++){
		if (allpos.count(i) > 0){
			cnt= allpos[i];
		}else{
			cnt= 0;
		}
		// cout << i << " " << cnt << endl;
		if (cnt != oldcnt || i== chrlen-1){
			if (start_pos != i)
				ss << chr << "\t"<< start_pos << "\t" << i << "\t" << oldcnt <<  endl;
			start_pos= i;
		}
		oldcnt= cnt;
	}
	return ss.str();
}

/*
 * translate codon to amino acid
 */
char codon2aa(const string &codon)
{
  if (codon ==  "TCA") return 'S';
  if (codon ==  "TCG") return 'S';
  if (codon ==  "TCC") return 'S';
  if (codon ==  "TCT") return 'S';
  if (codon ==  "TTT") return 'F';
  if (codon ==  "TTC") return 'F';
  if (codon ==  "TTA") return 'L';
  if (codon ==  "TTG") return 'L';
  if (codon ==  "TAT") return 'Y';
  if (codon ==  "TAC") return 'Y';
  if (codon ==  "TAA") return '*';
  if (codon ==  "TAG") return '*';
  if (codon ==  "TGT") return 'C';
  if (codon ==  "TGC") return 'C';
  if (codon ==  "TGA") return '*';
  if (codon ==  "TGG") return 'W';
  if (codon ==  "CTA") return 'L';
  if (codon ==  "CTG") return 'L';
  if (codon ==  "CTC") return 'L';
  if (codon ==  "CTT") return 'L';
  if (codon ==  "CCA") return 'P';
  if (codon ==  "CCG") return 'P';
  if (codon ==  "CCC") return 'P';
  if (codon ==  "CCT") return 'P';
  if (codon ==  "CAT") return 'H';
  if (codon ==  "CAC") return 'H';
  if (codon ==  "CAA") return 'Q';
  if (codon ==  "CAG") return 'Q';
  if (codon ==  "CGA") return 'R';
  if (codon ==  "CGG") return 'R';
  if (codon ==  "CGC") return 'R';
  if (codon ==  "CGT") return 'R';
  if (codon ==  "ATT") return 'I';
  if (codon ==  "ATC") return 'I';
  if (codon ==  "ATA") return 'I';
  if (codon ==  "ATG") return 'M';
  if (codon ==  "ACA") return 'T';
  if (codon ==  "ACG") return 'T';
  if (codon ==  "ACC") return 'T';
  if (codon ==  "ACT") return 'T';
  if (codon ==  "AAT") return 'N';
  if (codon ==  "AAC") return 'N';
  if (codon ==  "AAA") return 'K';
  if (codon ==  "AAG") return 'K';
  if (codon ==  "AGT") return 'S';
  if (codon ==  "AGC") return 'S';
  if (codon ==  "AGA") return 'R';
  if (codon ==  "AGG") return 'R';
  if (codon ==  "GTA") return 'V';
  if (codon ==  "GTG") return 'V';
  if (codon ==  "GTC") return 'V';
  if (codon ==  "GTT") return 'V';
  if (codon ==  "GCA") return 'A';
  if (codon ==  "GCG") return 'A';
  if (codon ==  "GCC") return 'A';
  if (codon ==  "GCT") return 'A';
  if (codon ==  "GAT") return 'D';
  if (codon ==  "GAC") return 'D';
  if (codon ==  "GAA") return 'E';
  if (codon ==  "GAG") return 'E';
  if (codon ==  "GGA") return 'G';
  if (codon ==  "GGG") return 'G';
  if (codon ==  "GGC") return 'G';
  if (codon ==  "GGT") return 'G';

  return -1;

}


