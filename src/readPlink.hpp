//
//  readPlink.cpp
//  ReadPlinkGit
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//
// [[Rcpp::depends(BH)]]
#ifndef readPlink_hpp
#define readPlink_hpp

#include<iostream>
#include <sstream>
#include <stdio.h>
#include <boost/algorithm/string.hpp>

#include "plinkfun.hpp"
#include "IGESS_aux.hpp"
//#include <libiomp/omp.h>
#include <map>

using namespace std;
using namespace arma;
using namespace boost;

//GenoInfo ReadDataFromFile(string stringname);
class Chroms;

class SNP{
public:
    SNP(){

    }
    SNP(string name, int idx, int type){
        this -> name = name;
        this -> idx = idx;
        this -> type = type;
    }
    SNP(const SNP& snp){
        this -> name = snp.name;
        this -> idx = snp.idx;
        this -> type = snp.type;
    }
    string name;
    int idx;
    int type;
    bool operator<(const SNP&) const;
    bool operator>(const SNP&) const;
    bool operator != (const SNP&) const;
    bool operator == (const SNP&) const;

    //   friend ostream& operator<<(ostream& os, const SNP& obj);
};


const int const_chrom_no = 22;
Chroms read_snpnames(string filename, int P);
class Chroms{
public:
    //      chromsomes[const_chrom_no];
    Col<int> chromsome;
    Col<int> A1Array;
    Col<int> A2Array;
    vector<SNP> snps;
  //  vector<string> snpidvector;
    int chrom_no;
    int P;
    arma::Col<int> index;
    //to allocate the chromsomes for the chroms, chr_idx = -1 for one chromsome, chr_idx != -1 for const_chrom_no 22
    Chroms(int P){
        chromsome.resize(P);
        A1Array.resize(P);
        A2Array.resize(P);
    //    snps = new SNP[P];
        this -> P = P;
    }

    Chroms(string bimfile, int P){
        *this = read_snpnames(bimfile, P);
    }

    Chroms(Chroms& chrom){
      this -> chromsome = chrom.chromsome;
      this -> A1Array = chrom.A1Array;
      this -> A2Array = chrom.A2Array;
      this -> snps = chrom.snps;
   //   this -> snpidvector = chrom.snpidvector;
      this -> index = chrom.index;
      this -> chrom_no = chrom.chrom_no;
      this -> P = chrom.P;

    //   Col<int> A1Array;
    //   Col<int> A2Array;
    //   vector<SNP> snps;
    //   vector<string> snps_string;
    //   int chrom_no;
    //   int P;
    //   arma::Col<int> index;
    // //  *this = read_snpnames(bimfile, P);
    }

    Chroms(vector<string>& snpnames){
      this -> P = snpnames.size();
      for(int i = 0; i < this -> P; i++){
        SNP snp(snpnames[i], i, 1);
        snps.push_back(snp);
      }
    }

    void clear(){
      chromsome.reset();
      A1Array.reset();
      A2Array.reset();
      index.reset();
    }

    Chroms(){}

};

#define pvalue_v 0
#define zvalue_v 1

#define beta_ss 0
#define pvalue_ss 1

#define from_ss 0
#define from_x 1

/*get the positions of the identifiers in the column names*/
Col<int> getPositions(vector <string> fields, vector<string>& identifiers);


vec read_phenotypes(string filename, int N, int phenopos = 5);

float normalCFD(float value);

class GenoInfo;

class Summary{
public:
    int P;
    int chrom_no;
    int type;
    Chroms* chroms;
    mat* lpsummary;
    mat* convert_lpsummary;
    vector<string> gwasidentifiers;
    map< string, vector <string> > config;
    void convert(){
        clock_t t1 = clock();
        if(this -> type == zvalue_v){
            float zvalue = 0;
            for(uword i = 0; i < lpsummary ->n_rows; i++){
                for(uword j = 0; j < lpsummary -> n_cols; j++){
                    zvalue = (*lpsummary).at(i,j);
                    zvalue = abs(zvalue);
                    (*lpsummary).at(i,j) = 2*(1 - normalCFD(zvalue));
                }
            }
        }
        cout << "Convert z-value to p-value, time elapsed " << (clock() - t1) / CLOCKS_PER_SEC << endl;
    }
    Summary();
    Summary(string summaryfile, string configfile){
        std::ifstream ifs(configfile.c_str());
        if(!ifs.is_open()){
            cout <<"Config file is not properly provided!" << endl;
        }
        std::string line;
        while(!ifs.eof()){
            std::getline(ifs, line);
            if(line.length() == 0) continue;
            if(line.find("=") == -1) continue;
            vector <string> fields;
            boost::split( fields, line, boost::is_any_of("="));
            string key = fields[0];
            string value = fields[1];
            cout <<"Key="<<key << endl;
            boost::split( fields, value, boost::is_any_of(", "));
            config[key] = fields;
        }
        ifs.close();
        vector<string> identifiers;
        if(config.find("zvalue") != config.end()){
            type = zvalue_v;
            identifiers = config["zvalue"];
        }else if(config.find("pvalue") != config.end()){
            type = pvalue_v;
            identifiers = config["pvalue"];
        }else{
            cout <<"Error of Config File! No p-value or z-value config!" << endl;
        }

        string snpidentifier = "SNP";
         if(config.find("snp") != config.end()){
             snpidentifier = config["snp"][0];
         }

        string chromidentifier = "CHR";
        if(config.find("chr") != config.end()){
            chromidentifier = config["chr"][0];
        }

        new (this)Summary(summaryfile, identifiers, snpidentifier, chromidentifier);


    }
    Summary(string summaryfile, vector<string> identifiers, string snpidentifier = "SNP",
            string chromidentifier = "CHR");
    Summary(vector<string>& snps, mat* lp_summary);
    void cal_overlap(GenoInfo& genoinfo);
    void clear(){
        if(lpsummary != NULL)
            delete lpsummary;
        if(convert_lpsummary != NULL)
            delete convert_lpsummary;
        this -> chroms -> clear();
    }

};


class GenoInfo{
public:
    GenoInfo(){

    }
    GenoInfo(string stringname);

    ~GenoInfo(){
    }
    arma::Mat<unsigned> X;
    arma::vec y;
    Col<uword> xindex;
    Chroms chroms;
    int N;  //sample size
    int P;  //snp number
};


int chroms_overlap(Chroms& chrom_x, Chroms& chrom_s,Col<uword>& xindex, Col<uword>& sindex);
int snps_overlap(vector<SNP>& chrom_x_i, vector<SNP>& chrom_s_i, Col<uword>& xindex, Col<uword>& sindex);


template <typename T>
void keep_indices(arma::Mat<T>& mat_obj, int type, uvec index);


#endif
