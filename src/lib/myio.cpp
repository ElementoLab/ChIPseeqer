/*
 * =====================================================================================
 *
 *       Filename:  myio.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/08/2010 02:01:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  David Soong (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "myio.h"

using namespace std;

/*
 * check the existence of a file
 */
bool fileExists(const string &fname){
	return ifstream(fname.c_str());
}


/*
 * convert integer (num) to 8-bit binary string
 */
string itob(const int &num){
	string out="";
	for(int i = 7; i>=0; i--) {
		if((1<<i)&num)
			out+= "1";
		else 
			out+= "0";
	}
	return out;
}

/*
 * read a gzipped file
 *
 * NOTE:
 * 	for more information on file compression, look at lib/download/zlib-1.2.5/example.c 
 */
void readGZfile(const string &fname){
	// char * infilename= "tmp.gz";
	gzFile infile = gzopen(fname.c_str(), "rb");
	if (!infile ) {
		exit (1);
	}

	int buff_size= 10* 1024 * 1024; // 10 MB;
	char buffer[buff_size];
	int num_read = 0;
	while ((num_read = gzread(infile, buffer, sizeof(buffer)-1)) > 0) {
		// fwrite(buffer, 1, num_read, outfile);
		buffer[num_read]= '\0';
		// cout << "buffer= |"<< buffer << "| num_read= "<< num_read << endl;
		cout << buffer;
	}
	gzclose(infile);
}

/*
 * lazy way to write content to a gzipped file
 */
void writeGZfile(const string &fname, const string &content){
	gzFile fout = gzopen(fname.c_str(), "w");

	if (! fout ) {
		cerr << "Error writing to "<< fname << endl;
		exit (1);
	}

	int ans= gzwrite(fout, content.c_str(), content.length());
	if (ans <= 0 ){
		cerr << "Error writing to "<< fname << endl;
		exit(1);
	}
	gzclose(fout);
}

string fileSlurpWrap::readLine(){
	if (bStdin){
		getline(cin,str);
		return str;
	}else{
		// cerr << "IMPLEMENT\n";
		return (*fs).readLine();
	}
}

/*
 * whether EOF has reached
 */
bool fileSlurpWrap::eof(){
	if (bStdin){
		return (cin.peek() < 0);
	}else{
		// cerr << "IMPLEMENT\n";
		return (*fs).eof();
	}
	// iseof= ( ! getline(cin,str) );
	// return iseof;	
}


/*
 * read in a huge chunk of data
 * 	ifstream -> use read(buffer, buffer_len)
 * 	gzFile -> use read(fpt, buffer, buffer_len)
 * 	
 */
void fileSlurp::readChunk(){
	// read in a chunk of data 
	// 	note: this is not formatted data so have to append '\0' by my self
	if ( isGzipped){
		// gzipped file
		num_read = gzread(gzf, buffer, buffer_len);
		buffer[num_read]= '\0';
	}else{
		fs.read(buffer, buffer_len);
		buffer[fs.gcount()]= '\0';
		num_read= fs.gcount();
	}

	if (buffer[num_read -1] == '\n')
		intact_last_line= true;
	else
		intact_last_line= false;
		
	// cout << "eof pointer= "<< ptr_eof  << endl;
	// cout << "current file ptr= "<< fs.tellg() << endl; 		// get end of file position 

	//	string ttt= string(buffer);
	//
	//	cerr << "\n";
	//	cerr << "isGzipped= "<< isGzipped<< endl;
	//	cerr << "buffer_len= "<< buffer_len << "\nnum_read= "<< num_read<< endl;
	//	cerr << "buffer (first 30)= "<< ttt.substr(0,30) << "|"<< endl;
	//	cerr << "buffer (last 30)= "<< ttt.substr( ttt.length()-30, 30)<< "|"<< endl;
	//	cerr << "buffer string= "<< ttt << "|"<< endl;
	//	cerr << "buffer string length= "<< ttt.length() << endl;
	//	cerr << "leftover= "<< leftover<< "|" << endl;
	
	tmpstr= leftover+ string(buffer);
	// cerr << "tmpstr string length= "<< tmpstr.length() << endl;
	// cout << "## reading in a new chunk # " << blockid << "\ndata (first 20 char): "<< tmpstr.substr(0, 20) << endl;

	// all_lines.clear();
	Tokenize_strict(tmpstr, all_lines, "\n");
	// Tokenize(tmpstr, all_lines, "\n");
	
	// cerr << "all_lines size= "<< all_lines.size()<< endl;

	last_line_id= all_lines.size()- 1;

	// DEBUG
	//if (all_lines.size()>0){
	//	cerr << "first line= "<< all_lines[0]<< endl;
	//	cerr << "last line= "<< all_lines[last_line_id]<< endl;
	//}

	current_line_id= 0;
	// cout << "\t"<< all_lines.size() << " lines\n";
	hasChunk= true;
	blockid++;
}

/*
 * TO DO:
 * 	modify to enable reading from gzipped files!!
 *
 */
string fileSlurp::readLine(){
	// if has a data chunk && not last line, output a line
	if (hasChunk){
		if (current_line_id < last_line_id){
			// cout << "now at line "<< current_line_id << " last= "<< last_line_id << endl;
			output= all_lines[current_line_id];
		}else{
			// cout << "	last line in chunk -> read in another chunk\n";
			leftover= all_lines[current_line_id];
			// cout << "	leftover= "<< leftover <<endl;
			//
			// problem if leftover happens to be intact (with \n)
			if (intact_last_line){
				leftover+= "\n";
			}

			// last line of block
			// 	1. check if we are at the very end of the file
			//	2. added support for gzipped files
			if ( ! isGzipped && streamoff(fs.tellg()) != -1){
				// if normal file -> read in next chunk only if we are not at eof
				readChunk();
				output= all_lines[0];
			}else if ( isGzipped && gztell(gzf) != -1 ){
				// if gzipped file -> read in next chunk only if we are not at eof
				readChunk();
				output= all_lines[0];
			}else{
				// cout << "EOF last line\n";
				output= all_lines[last_line_id];
			}
		}
	}else{
		// read in another chunk
		// concatenate with last line in previous chunk
		// reset line ids
		// return the first line
		readChunk();
		output= all_lines[0];
	}
	current_line_id++;
	return output;
}

/*
 * read in the next line from buffer. when reaching the end of buffer, read another block
 */
void fileSlurp::readAllLines(){
	cerr << "WARNING: not sure if readAllLines() implementation finished..." << endl;
	cout << "buffer_len= "<< buffer_len << endl;
	
	string tmpstr;
	tmpstr.reserve(buffer_len+1);
	vector <string> tokens;

	while (1){
		// do we need a new chunk?
		cout << "\n### reading in chunk # "<< blockid<< endl;
		fs.read(buffer, buffer_len);
		
		cout << "eof pointer= "<< ptr_eof  << endl;
		cout << "current file ptr= "<< fs.tellg() << endl; 		// get end of file position 
		tmpstr= string(buffer);
		cout << "data: "<< tmpstr.substr(0, 20) << endl;

		tokens.clear();
		Tokenize_strict(tmpstr, tokens, "\n");

		cout << "last line: "<< tokens[tokens.size()-1] << endl;
		
		// if last line && not end of file -> read another chunk

		blockid++;
		// if not at the end of the file, read in another chunk, concatenate with this last line 
		if (streamoff(fs.tellg()) != -1){
			cout << "have to read in another chunk...\n";
			cout << "leftover= "<< leftover << endl;

		}else{
			cout << "last line= "<< leftover << endl;
			break;
		}
	}
}


/*
 * load a hash from file
 */
map <string, string> loadHash(string fname, int column1, int column2){
	string line;
	vector <string> tokens;
	map <string, string> output;
	fileSlurp fs(fname);
	while(! fs.eof()){
		line= fs.readLine();
		if (line.at(0)== '#')
			continue;
		Tokenize_strict(line, tokens, "\t");
		output[tokens[column1]]= tokens[column2];
	}
	fs.close();
	return output;
}

bool isnumeric(const string & num){
	// cout << "checking: "<< num<< endl;
	for (int i=0; i< (int)num.length(); i++){
		if ( ! isdigit(num.at(i)) )
			return false;
	}
	return true;
}


/*
 * string fullstr starts with substr
 */
bool startsWith(const string & fullstr, const string & substr){
	size_t found= fullstr.find(substr);
	if (found != string::npos && found == 0){
		return true;
	}else{
		return false;
	}
}


/*
 * checks if a string (fullString) ends with (ending)
 */
bool hasEnding (std::string const &fullString, std::string const &ending)
{
	if (fullString.length() > ending.length()) {
		return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
}

/*
 * change each element of the string to upper case and return a new string
 */ 
string stringToUpper(string strToConvert){
	// strToConvert copied from input (pass by value)
	for(unsigned int i=0;i<strToConvert.length();i++)
	{
		strToConvert[i] = toupper(strToConvert[i]);
	}
	return strToConvert;
}


/*
 * change each element of the string to upper case and return a new string
 */ 
string stringToLower(string strToConvert){
	// strToConvert copied from input (pass by value)
	for(unsigned int i=0;i<strToConvert.length();i++)
	{
		strToConvert[i] = tolower(strToConvert[i]);
	}
	return strToConvert;
}

/*
 * Tokenize the string using ANY of the specified delimiters. Every encounter with a delimiter produces a split. No further split is 
 * created if we see a delimiter at the end of the string. 
 *
 * NOTE: 
 *	1. A run of n consecutive delimiters will produce n splits (n-1 tokens between the first and last delimiters).
 *	   c.f. Tokenize() produces only one split when encountering a stretch of n delimiters. 
 *
 * Examples:
 * 	Tokenize_strict with "\n" on: 
 * 		1111\n2222\n will produce [1111,2222]
 * 		1111\n2222   will produce [1111,2222]
 * 		\n1111\n2222 will produce [,1111,2222]
 * 		\n\n1111\n2222\n will produce [,,1111,2222]
 *
 * TO DO:
 * 	rename Tokenize() to Tokenize_soft()
 *
 */
void Tokenize_strict(const string& str, vector<string>& tokens, const string& delimiters ){
	tokens.clear();

	string::size_type pos= 0;
	string::size_type lastPos= str.find_first_of(delimiters, 0);

	// cerr << "first 30 chr= "<< str.substr(0, 30)<< "|"<< endl;

	// cerr << "::npos= "<< string::npos<< endl;
	// while (string::npos != pos || string::npos != lastPos)
	while (1)
	{
		// string toadd= str.substr(pos, lastPos- pos);
		// cerr << "pos= "<< pos << " lastpos= "<< lastPos<< " adding: "<< toadd << "|"<<endl;
		tokens.push_back(str.substr(pos, lastPos- pos));

		pos= lastPos+1;
		lastPos = str.find_first_of(delimiters, pos);
		if (lastPos == string::npos) {
			// no more delimiters found; add the remaining piece
			// cerr << "no more delimiters found, pos= "<< pos << ", lastPos= "<< lastPos << " len= "<< str.length() << endl;
			// cerr << "last piece: "<< str.substr(pos) <<endl;
			if (pos < str.length()){
				// cerr << "adding last piece "<< str.substr(pos)<< endl;
				tokens.push_back(str.substr(pos));
			}
			return;
		}
	}
}


/************************************************************************************ 
* slit string into tokens using delimiters (can be more than one character)
*
* NOTE:
* 	1. a run of multiple delimiters will produce NO token between them,
* 		e.g. "1,2,,,3" will give [1,2,3] instead of [1, 2, NULL, NULL, 3]
*		=> use Tokenize_strict to give one token per delimiter
* 	2. "tokens" is passed by reference and is reusable 
*
* NOTE:	we now automatically clear the content of tokens before use !!
*
*************************************************************************************/ 
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = "\t"){
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	tokens.clear();

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

/*
 * Expand the ~ notation to absolute home directory
 */
string expandHome (string path){
	if (path.empty()){
		return "";
	}

	if (path.at(0) == '~'){
		char *homedir= getenv("HOME");
		path.replace(0, 1, homedir);
		return (path);
		// cout << "expanding to "+ path << endl;
	}else{
		return path;
	}
}

/*
 * convert vector to string 
 */
string vector2string(vector <string> &v, string sep){
	int ii;
	string output= "";
	if (v.size()== 0)
		return output;
	else{	
		for(ii=0; ii < (int) v.size(); ii++)
		{
			output+= v[ii];
			output+= sep;
		}
		output.erase(output.length()-1);
		return output;
	}
}

/*
 * get command options
 *
 * E.g.
 *     char * filename = getCmdOption(argv, argv + argc, "-file");
 *
 */
string getCmdOption(char ** argv, int argc, const std::string & option){
	/*
	 * 
	 char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		return *itr;
	}
	return 0;
	*/
	string tmp;
	string::size_type found;
	// cout << "option= "<< option<< endl;
	for (int i=0; i< argc; i++){
		tmp= argv[i];
		found= tmp.find(option);
		// cout << tmp << " " << found  << " " << (found== string::npos)<< endl;
		// cout << tmp << endl;
		if (found != string::npos && found == 0){
			// cout << "\tfound" << endl;
			if (i+1 < argc){
				return string(argv[i+1]);
			}else{
				return string();
			}
		}
	}
	return string();
}

/*
 * get binary command options, e.g. --unique (without getting the next arg)
 */
bool getCmdOptionBin(char ** argv, int argc, const std::string & option){
	string tmp;
	string::size_type found;
	// cout << "option= "<< option<< endl;
	for (int i=0; i< argc; i++){
		tmp= argv[i];
		found= tmp.find(option);
		// cout << tmp << " " << found  << " " << (found== string::npos)<< endl;
		// cout << tmp << endl;
		if (found != string::npos && found == 0){
			return true;
		}
	}
	return false;
}

HG18::HG18(){
	chr2len["chr1"]=247249719;
	chr2len["chr2"]=242951149;
	chr2len["chr3"]=199501827;
	chr2len["chr4"]=191273063;
	chr2len["chr5"]=180857866;
	chr2len["chr6"]=170899992;
	chr2len["chr7"]=158821424;
	chr2len["chr8"]=146274826;
	chr2len["chr9"]=140273252;
	chr2len["chr10"]=135374737;
	chr2len["chr11"]=134452384;
	chr2len["chr12"]=132349534;
	chr2len["chr13"]=114142980;
	chr2len["chr14"]=106368585;
	chr2len["chr15"]=100338915;
	chr2len["chr16"]=88827254;
	chr2len["chr17"]=78774742;
	chr2len["chr18"]=76117153;
	chr2len["chr19"]=63811651;
	chr2len["chr20"]=62435964;
	chr2len["chr22"]=49691432;
	chr2len["chr21"]=46944323;
	chr2len["chrX"]=154913754;
	chr2len["chrY"]=57772954;
	/* 
	chr2len["chrM"]=16571;
	chr2len["chr6_cox_hap1"]=4731698;
	chr2len["chr6_qbl_hap2"]=4565931;
	chr2len["chr17_random"]=2617613;
	chr2len["chr6_random"]=1875562;
	chr2len["chr5_h2_hap1"]=1794870;
	chr2len["chrX_random"]=1719168;
	chr2len["chr21_random"]=1679693;
	chr2len["chr1_random"]=1663265;
	chr2len["chr9_random"]=1146434;
	chr2len["chr8_random"]=943810;
	chr2len["chr4_random"]=842648;
	chr2len["chr15_random"]=784346;
	chr2len["chr3_random"]=749256;
	chr2len["chr7_random"]=549659;
	chr2len["chr19_random"]=301858;
	chr2len["chr22_random"]=257318;
	chr2len["chr11_random"]=215294;
	chr2len["chr13_random"]=186858;
	chr2len["chr2_random"]=185571;
	chr2len["chr5_random"]=143687;
	chr2len["chr10_random"]=113275;
	chr2len["chr16_random"]=105485;
	chr2len["chr22_h2_hap1"]=63661;
	chr2len["chr18_random"]=4262;
	*/
}

int HG18::getChrLen(const string &chr){
	if (chr2len.count(chr)> 0){
		return chr2len[chr];
	}else{
		return 0;
	}
}

/*
 * return a vector of chromosome names (1..22, X, Y) 
 *
 * NOTE: _random chromosomes are excluced 
 * 
 */
vector <string> HG18::getChrs(){
	vector <string> out= getKeys <string, int> (chr2len); 
	sort(out.begin(), out.end());
	return out;
}


/*
 * get an ordered number array (start... end-1 )
 */
vector <int> get_num_array(int start, int end){
	vector <int> read_idx;
	for (int i= start; i< end; i++) 
		read_idx.push_back(i);
	return read_idx;
}

/*
 * get a randomly shuffled array (range between start and end-1 )
 */
vector <int> get_rnd_seq_array(int start, int end){
	vector <int> read_idx= get_num_array(start, end);
	random_shuffle ( read_idx.begin(), read_idx.end() );
	return read_idx;
}
		
/*
 * put content in a text file fname
 */
bool lazyWriteFile(const string &fname, const string &content){
	ofstream outfile;
	outfile.open (fname.c_str());
	if (outfile.is_open())
	{
		outfile << content;
		outfile.close();
		return true;
	}else{
		cerr << "Error opening file " << fname << endl;
		return false;
	}
}



/*
 * process_mem_usage(double &, double &) - takes two doubles by reference,
 * attempts to read the system-dependent data for a process' virtual memory
 * size and resident set size, and return the results in KB.
 *
 * On failure, returns 0.0, 0.0 
 *
 * TO DO:
 * 	this function is for Linux only. find solution for OSX and windows
 *
 */
void process_mem_usage(double& vm_usage, double& resident_set){
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage     = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	//
	unsigned long vsize= 0;
	long rss= 0;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
			>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
			>> utime >> stime >> cutime >> cstime >> priority >> nice
			>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage     = vsize / 1024.0;
	resident_set = rss * page_size_kb;
	if (vm_usage == 0.0 || resident_set == 0.0)
		cerr << "couldn't read /proc stats\n";
}

