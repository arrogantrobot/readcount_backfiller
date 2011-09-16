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
    base_data():dp(0),bq(0.0),mq(0.0){}
};

struct read_count {
    string chrom;
    int pos, depth;
    base_data   a,c,g,t,n;
};

void split_line(vector<string> * v, const string& s);
void split_record(vector<string> * v,const string& s);
void process_read_count_line( vector<string> * v, read_count * r);
void make_sample_field(const string& ref, const string& alt, const string& format, const int& depth, const read_count& rc, const string& sample_field, string * answer);
void split_record_comma(vector<string> * v,const string& s);
string to_string(const float& f);
string to_string(const int& f);

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

    bool rc_eof=false;
    bool done=false;
    bool get_vcf=1,get_rcount=1,get_sample=1;
    while( ! done ){ 

        if(sample_buff != ".") {
            cout << sample_buff << endl;
        } else {
            string answer;

            split_line(&vcf,vcf_buff);
            split_line(&rcount,readcount_buff);

            if(rc_eof){

                read_count rc;
               // rc.depth = "0";
                make_sample_field(vcf[REF],vcf[ALT],vcf[FORMAT],0,rc,sample_buff, &answer);
                cout << "rc\t" <<  answer << endl;

            }else{
                if( vcf[CHROM] != rcount[CHR] ){
                    
                    read_count rc;
                    rc.depth = 0;
                    make_sample_field(vcf[REF],vcf[ALT],vcf[FORMAT],0,rc,sample_buff, &answer);
                    get_rcount = 0;
                    cout << "1\t" <<  answer << endl;

                } else if ( vcf[POS] != rcount[PS] ) {
                    read_count rc;
                    rc.depth = 0;
                    make_sample_field(vcf[REF],vcf[ALT],vcf[FORMAT],0,rc,sample_buff, &answer);
                    get_rcount = 0;
                    cout << "2\t" << answer << endl;

                } else {
                    read_count rc;
                    process_read_count_line(&rcount,&rc);
                    make_sample_field(vcf[REF],vcf[ALT],vcf[FORMAT],0,rc,sample_buff, &answer);
                    get_rcount = 1;
                    cout << "3\t" <<  answer << endl;
                }
            }
        }                       
        


        for(vector<string>::iterator it=rcount.begin(); it!=rcount.end();it++){
            //cout << *it << endl;
        }

        if(get_vcf) getline(vcf_stream, vcf_buff);
        if(get_rcount) getline(rcount_stream, readcount_buff);
        if(get_sample) getline(sample_stream, sample_buff);
        if(rcount_stream.eof()) rc_eof=1;
        //eof = (vcf_stream.eof() || rcount_stream.eof() || sample_stream.eof());
        done = (vcf_stream.eof());
    }
    return 0;
}

void make_sample_field(const string& ref, const string& alt, const string& format, const int& depth, const read_count& rc,const string& sample_field, string * answer){
    string gt,mq,ad,bq;
    if( ref == "A"){
        mq = to_string(rc.a.mq);
    } else if ( ref == "C") {
        mq = to_string(rc.c.mq);
    } else if ( ref == "G") {
        mq = to_string(rc.g.mq);
    } else if ( ref == "T") {
        mq = to_string(rc.t.mq);
    } else {
        throw "wtf, ref is neither A nor C nor G nor T.";
    }
    vector<string> alts;
    alts.push_back(ref);
    split_record_comma(&alts, alt);
    for(vector<string>::iterator it=alts.begin();it!=alts.end();it++){
        if(bq != ""){
            bq+=",";
        }
        if(ad != ""){
            ad+=",";
        }
        if( *it == "A"){
            bq += to_string(rc.a.bq);
            ad += to_string(rc.a.dp);
        } else if ( *it == "C") {
            bq += to_string(rc.c.bq);
            ad += to_string(rc.c.dp);
        } else if ( *it == "G") {
            bq += to_string(rc.g.bq);
            ad += to_string(rc.g.dp);
        } else if ( *it == "T") {
            bq += to_string(rc.t.bq);
            ad += to_string(rc.t.dp);
        } else {
            throw "ref is neither A nor C nor G nor T.";
        }
    }

    vector<string> formats;
    split_record(&formats, format);
 
    if(rc.depth > 0){
        gt = "0/0";
    } else {
        gt = mq = ad = bq = ".";

    }

        
    for(vector<string>::iterator it=formats.begin();it!=formats.end();it++){
        if(*answer != ""){
            *answer+=":";
        }
        if( *it == "AD"){
            *answer += ad;
        } else if ( *it == "BQ") {
            *answer += bq;
        } else if ( *it == "GT") {
            *answer += gt;
        } else if ( *it == "MQ") {
            *answer += mq;
        } else if ( *it == "DP") {
            *answer += rc.depth;
        } else {
            *answer += ".";
        }
    }
}

/*
void split_line(vector<string> * v,const string& s) {
    if(! v->empty()){
        v->clear();
    }
    tokenizer<> tok(s);
    for(tokenizer<>::iterator it=tok.begin();it!=tok.end();it++){
        v->push_back(*it);
    }
}
*/
void split_record(vector<string> * v,const string& s) {
    char_separator<char> sep(":");
    tokenizer<char_separator<char> > tok(s, sep);
    for(tokenizer<char_separator<char> >::iterator it=tok.begin();it!=tok.end();it++){
        v->push_back(*it);
    }
}

void split_record_comma(vector<string> * v,const string& s) {
    char_separator<char> sep(",");
    tokenizer<char_separator<char> > tok(s, sep);
    for(tokenizer<char_separator<char> >::iterator it=tok.begin();it!=tok.end();it++){
        v->push_back(*it);
    }
}


void split_line(vector<string> * v,const string& s) {
    if(! v->empty()){
        v->clear();
    }
    char_separator<char> sep("\t");
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
string to_string(const float& f){
    return lexical_cast<string>(f);
}
string to_string(const int& f){
    return lexical_cast<string> (f);
}


void process_read_count_line( vector<string> * v, read_count * r){

    r->chrom = lexical_cast<int> ( (*v)[CHR].c_str());
    r->pos = lexical_cast<int>  ((*v)[PS].c_str()); 
    r->depth = lexical_cast<int> ( (*v)[DP].c_str());
    vector<string> bases;
    //cout << "here!" << endl;
    split_record(&bases,(*v)[A]);
    
    //cout << "A = " << (*v)[A] <<  endl;
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
