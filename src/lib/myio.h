/*
 * =====================================================================================
 *
 *       Filename:  myio.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/08/2010 01:54:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  David Soong (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef MYIO_H
#define MYIO_H

#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include <time.h>
#include <set>

extern "C"{
	#include "third_party/zlib/zlib.h"
}

const int MB= 1024*1024;

using namespace std;

// TO DO:
// 	categorize these functions and give them namespaces
//	-> I/O, string, etc.
//
string itob(const int &num);
bool isnumeric(const string & num);
bool startsWith(const string & fullstr, const string & substr);
bool hasEnding (std::string const &fullString, std::string const &ending);

void Tokenize_strict(const string& str, vector<string>& tokens, const string& delimiters = "\t");
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
string expandHome (string path);
string vector2string(vector <string> &v, string sep);
string stringToUpper(string strToConvert);
string stringToLower(string strToConvert);

string getCmdOption(char ** argv, int argc, const std::string & option);
bool getCmdOptionBin(char ** argv, int argc, const std::string & option);
map <string, string> loadHash(string fname, int column1, int column2);

bool lazyWriteFile(const string &fname, const string &content);
void readGZfile(const string &fname);
void writeGZfile(const string &fname, const string &content);
bool fileExists(const string &fname);

vector <int> get_num_array(int start, int end);
vector <int> get_rnd_seq_array(int start, int end);

void process_mem_usage(double& vm_usage, double& resident_set);



/*
 * read in the file in huge chunks to reduce I/O system calls
 *
 * TO DO:
 * 	1. actually we can use pure C code to further speed up
 *	2. add option to set buffer size
 */
class fileSlurp {
	private: 
		int buffer_len;	// = 10000000;
		char * buffer;
		char * buffer_st;
		char * s;
		char * line;

		ifstream fs;
		gzFile   gzf;				// reading gzipped files 
		bool     isGzipped;				// gzipped file?
		int      num_read;

		int blockid;
		streamoff ptr_eof;
		bool bHasLeftover;
		vector <string> all_lines;
		bool hasChunk;
		int current_line_id;
		int last_line_id;
		bool intact_last_line;
		string tmpstr; 
		string leftover;
		string gFname;
		string output;				// allocated globally to avoid repetitive obj cons/destruction

	public:
		fileSlurp(int tmp_buff_len= 20* 1024* 1024){
			// TO DO:
			// 	test buffer_len to see if enough memory can be allocated
			init(tmp_buff_len);
		}
		inline void init(int tmp_buff_len){
			// initialization goes here
			// TO DO:
			// 	set buffer size
			// 	-> check if buffer exists -> free ? 
			// 	-> reallocate buffer and reset tmpstr size
			buffer_len= tmp_buff_len;		//20* 1024* 1024; 			// 10 MB ;
			buffer= new char [buffer_len+1];
			tmpstr.reserve(buffer_len+1);
			all_lines.reserve(buffer_len+1);

			blockid= 0;
			hasChunk= false;
			current_line_id= 0;
			last_line_id= 0;
			bHasLeftover= false;
			leftover= ""; 
			intact_last_line= false;
			gFname= "";

			// new addition: read gzipped files
			isGzipped= false;
			num_read= -1;
		}
		fileSlurp(string fname, int tmp_buff_len= 20* 1024* 1024){
			// cout << "in fileSlurp(fname), opening "<< fname << endl;
			init(tmp_buff_len);

			// determine it it's a gzipped file
			if (hasEnding(fname, ".gz")){
				isGzipped= true;				
				gzf= gzopen(fname.c_str(), "rb");
				if (! gzf){
					cerr << "Error: couldn't open "<< fname << endl;
					exit(1);
				}
			}else{
				// TO DO:
				// 	add code to read from STDIN
				isGzipped= false;
				fs.open(fname.c_str(), ios::binary);
				if (fs.fail()){
					cerr << "Error: couldn't open "<< fname << endl;
					exit(1);
				}
			}
			gFname= fname;

			fs.seekg(0, ios::end);		// go to end of file
			ptr_eof= fs.tellg(); 		// get end of file position 
			fs.seekg(0, ios::beg);		// back to beginning of file
		}

		void reset(){
			// cerr << "resetting fileSlurp()..."<< endl;
			tmpstr.clear();
			all_lines.clear();
			delete [] buffer;
			output.clear();
		}
		~fileSlurp(){
			// cerr << "closing fileSlurp "<< gFname << endl;
			if (fs.is_open()){
				fs.close();
				reset();
			}
		}
		void close();
		bool is_open();
		void readChunk();
		string readLine();
		void readAllLines();
		bool eof();
};

inline void fileSlurp::close(){
	if ( isGzipped ){
		gzclose(gzf);
		reset();
	}else{
		if (fs.is_open()){
			fs.close();
			reset();
		}
	}
}

inline bool fileSlurp::is_open(){
	return fs.is_open();
}

inline bool fileSlurp::eof(){
	if ( isGzipped ){
		// TO DO: 
		// 	double check
		//cout << "num_read= "<< num_read << " eof= " << gzeof(gzf)<< " gztell= "<< gztell(gzf) << 
		//	" current_line_id= "<< current_line_id << " last_line_id= "<< last_line_id << endl;
		return (gzeof(gzf)== 1 &&  (current_line_id > last_line_id) );
	}else{
		return (streamoff(fs.tellg())== -1 && (current_line_id > last_line_id));
	}
}

/*
 * wrapper for original fileSlurp that includes STDIN capability
 *
 * TO DO:
 * 	finish implementation
 */
class fileSlurpWrap {
	private:
		fileSlurp *fs;
		string str;
		bool bStdin;

	public:
		fileSlurpWrap(){
			cerr << "init fileSlurpWrap()"<< endl;
			cerr << "finish implementation to include gzip and slurp support\n";
			str= "";
			bStdin= false;
		};
		
		/*
		fileSlurpWrap(int tmp_buff_len= 20* 1024* 1024){
			// from FILE
			cerr << "does not seem to be used\n";
		} 
		*/

		// void init(int tmp_buff_len);
		
		fileSlurpWrap(string fname, int tmp_buff_len= 20* 1024* 1024){
			if (fname.empty() || fname == "-"){
				// use STDIN
				cerr << "init. use stdin"<< endl;
				str= "";
				bStdin= true;
			}else{
				cerr << "using file stream or gzopen to open "<< fname << endl;
				fs= new fileSlurp(fname, tmp_buff_len);
				bStdin= false;
			}

		}

		// void reset();
		//
		~fileSlurpWrap(){
			// cerr << "~fileSlurpWrap()"<< endl;	
			if ( ! bStdin){
				delete fs;
			}
		}
		void close(){
			// cerr << "closing STDIN (for compatability; do nothing)"<< endl;
			if ( ! bStdin ){
				(*fs).close();
			}
		}

		// bool is_open();
		// void readChunk();
		string readLine();
		// void readAllLines();
		bool eof();
};

/* 
 * data structure to hold aligned reads
 * 	note: 
 * 		start= 0-based
 * 		end= 1-based
 *
 * */
class read_str{
	//	static int chrid= 1;
	public:
		string name;
		string chr;
		char strand;
		int start;
		int end;

		read_str(){
			start= end= -1;
			name= chr= "";
			strand= 0;
		}

};


class timer {
	private:
		clock_t st;
		clock_t en;
	public:
		timer(){
			st= time(NULL);
			en= -1;
		}
		void start(){
			st= time(NULL);
		}
		void stop(){
			en= time(NULL);
		}
		int getElapsedSec(){
			// cout << st << "-" << en<< endl;
			if (en == -1){
				return -1;
			}else{
				return (en- st);
			}
		}
};

class HG18{
	private:
		map <string, int> chr2len;
	public:
		HG18();
		int getChrLen(const string &chr);
		vector <string> getChrs();
};


/**************************************************************************
 * template declarations go here
 **************************************************************************/
template <class T> string castToString(const T t);
template <class T> string vector2string2(const vector <T> &v, string sep);
template <class T1, class T2> vector <T1> getKeys( map<T1, T2> m );
template <class T1, class T2> bool hasKey(map <T1, T2> table, string key);
template <class T1, class T2> void deleteMapComponents( map < T1, T2 > &m );


/*
 * template implementations go here
 */

/*
 * convert vector to string 
 */
template <class T>
string vector2string2(const vector <T> &v, string sep){
	int ii;
	string output= "";
	if (v.size()== 0)
		return output;
	else{	
		for(ii=0; ii < (int) v.size(); ii++)
		{
			// cout << "position " << ii << " -> " << v[ii] << " = " << castToString <T> (v[ii]) << endl;
			output+= castToString <T> (v[ii]);
			output+= sep;
		}
		output.erase(output.length()- sep.length());
		return output;
	}
}

template <class T>
string castToString(const T t)
{
	stringstream s;
	s << t;
	return s.str();
}

/* 
 * get a vector of keys 
 */
template <class T1, class T2>
vector <T1> getKeys( map<T1, T2> m ){
	vector<T1> v;
	for (typename map<T1, T2>::iterator it = m.begin(); it != m.end(); ++it) {
		v.push_back(it->first);
		// cout << it->first << "\n";
	}
	return v;
}

/*
 * note: can use map.count(key) instead
 */
template <class T1, class T2>
bool hasKey(map <T1, T2> table, string key){
	typename map <T1, T2 > ::const_iterator it;
	it= table.find(key);
	return (it != table.end());
}

template <class T1, class T2>
void deleteMapComponents( map < T1, T2 > &m ){
	// iterate and delete 
	typename map<T1,T2>::iterator it;
	for ( it= m.begin() ; it != m.end(); it++ ){
		// cout << "deleting " << (*it).first << endl; // " => " << (*it).second << endl;
		(*it).second.clear();		// STL object -> use clear()
	}
}


#endif
