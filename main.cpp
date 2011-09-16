#include <fstream>
#include <iostream>
#include <string>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#define CHROM 0 
#define POS 1 
#define ID  2
#define REF 3
#define ALT 4
#define QUAL    5
#define FILTER  6
#define INFO    7
#define FORMAT  8

#define CHR 0
#define PS  1
#define RF  2
#define DP  3
#define EQ_DEPTH    4
#define A   5
#define C   6
#define G   7
#define T   8
#define N   9


using namespace std;
using namespace boost;

struct base_data {
        int dp;
        float bq,mq;
};

struct read_count {
    int chrom, pos, depth;
    base_data   a,c,g,t,n;
};

void split_line(vector<string> * v, const string& s);
void split_record(vector<string> * v,const string& s);
void process_read_count_line( vector<string> * v, read_count * r);

int main(int argc, char * argv[]) {

    string vcf_cols = argv[1];
    string readcount = argv[2];
    string sample_col = argv[3];

    ifstream vcfcols(vcf_cols.c_str(), ios_base::in | ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> vcf_cols_in;
    vcf_cols_in.push(boost::iostreams::gzip_decompressor());
    vcf_cols_in.push(vcfcols);

    ifstream rcount_is(readcount.c_str(), ios_base::in | ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> rcount_in;
    rcount_in.push(boost::iostreams::gzip_decompressor());
    rcount_in.push(rcount_is);

    ifstream sample_is(sample_col.c_str(), ios_base::in | ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> sample_col_in;
    sample_col_in.push(boost::iostreams::gzip_decompressor());
    sample_col_in.push(sample_is);

    string vcf_buff;
    string readcount_buff;
    string sample_buff;

    std::istream vcf_stream(&vcf_cols_in); 
    std::istream rcount_stream(&rcount_in); 
    std::istream sample_stream(&sample_col_in); 

    getline(vcf_stream, vcf_buff);
    getline(rcount_stream, readcount_buff);
    getline(sample_stream, sample_buff);


    vector<string> vcf;
    vector<string> rcount;
    vector<string> sample;

    bool eof=false;

    while( ! eof ){ 

        split_line(&vcf,vcf_buff);
        split_line(&rcount,readcount_buff);
        split_line(&sample,sample_buff); 
/*
        if(sample_buff == "."){
            cout << "Found it!"<< endl;
        }
        
        for(vector<string>::iterator it=sample.begin(); it!=sample.end();it++){
            cout << *it << endl;
        }
        cout << "size of sample vector is " << sample.size() << endl;
        return 0;
*/
        if( vcf[CHROM] != rcount[CHR] ){
            


        } else {


        }
            //if( to_int(vcf[POS]) != to_int(rcount[PS]) ) {

                                   
        


        for(vector<string>::iterator it=vcf.begin(); it!=vcf.end();it++){

        }

        getline(vcf_stream, vcf_buff);
        getline(rcount_stream, readcount_buff);
        getline(sample_stream, sample_buff);

        eof = (vcf_stream.eof() || rcount_stream.eof() || sample_stream.eof());

    }
}

void make_sample_field(const char& ref, const string& alt, const string& format, const int& depth, const read_count& rc){



}

void split_line(vector<string> * v,const string& s) {

    tokenizer<> tok(s);
    for(tokenizer<>::iterator it=tok.begin();it!=tok.end();it++){
        v->push_back(*it);
    }
}

void split_record(vector<string> * v,const string& s) {
    char_separator<char> sep(":");
    tokenizer<char_separator<char> > tok(s, sep);
    for(tokenizer<char_separator<char> >::iterator it=tok.begin();it!=tok.end();it++){
        v->push_back(*it);
    }
}

int to_int( const string& s){
    return lexical_cast<int>(s);
}
float to_float(const string& s){
    return lexical_cast<float>(s);
}

void process_read_count_line( vector<string> * v, read_count * r){

    r->chrom = lexical_cast<int> ( (*v)[CHR].c_str());
    r->pos = lexical_cast<int>  ((*v)[POS].c_str()); 
    r->depth = lexical_cast<int> ( (*v)[DP].c_str());
    vector<string> bases;

    split_record(&bases,(*v)[A]);
    r->a.dp = lexical_cast<int> (bases[1]);
    r->a.bq = lexical_cast<float>(bases[2]);
    r->a.mq = lexical_cast<float>(bases[3]);
    bases.clear();

    split_record(&bases,(*v)[C]);
    r->c.dp = lexical_cast<int> (bases[1]);
    r->c.bq = lexical_cast<float>(bases[2]);
    r->c.mq = lexical_cast<float>(bases[3]);
    bases.clear();

    split_record(&bases,(*v)[G]);
    r->g.dp = lexical_cast<int> (bases[1]);
    r->g.bq = lexical_cast<float>(bases[2]);
    r->g.mq = lexical_cast<float>(bases[3]);
    bases.clear();

    split_record(&bases,(*v)[T]);
    r->t.dp = lexical_cast<int> (bases[1]);
    r->t.bq = lexical_cast<float>(bases[2]);
    r->t.mq = lexical_cast<float>(bases[3]);
    bases.clear();
}
