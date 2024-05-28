/*
 * Metodología de la Programación: Kmer1
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 27 de octubre de 2023, 11:03
 */

#include "KmerFreq.h"
#include "Kmer.h"
#include <string>
#include <iostream>
using namespace std;

KmerFreq::KmerFreq(){
    _kmer = Kmer();
    _frequency = 0;
}

Kmer KmerFreq::getKmer() const{
    return _kmer;
}

int KmerFreq::getFrequency() const{
    return _frequency;
}

void KmerFreq::setKmer(Kmer kmer){
    _kmer = kmer;
}

void KmerFreq::setFrequency(int frequency){
    if (frequency < 0)
    {
        throw out_of_range("Valor incorrecto de frecuencia.");
    }
    else
    {
        _frequency = frequency;
    }
}

string KmerFreq::toString(){
    string freq = to_string(_frequency);
    return _kmer.toString() + freq;                 
}

void KmerFreq::write(ostream &outputStream){
    outputStream << _kmer.toString() << " " << _frequency;
}

void KmerFreq::read(istream &inputStream){
    inputStream >> _kmer >> _frequency;
}

ostream operator<<(ostream &os, KmerFreq kmerFreq){
    kmerFreq.write(os);
    
    return os;
}

istream operator>>(istream &is, KmerFreq kmerFreq){   
    kmerFreq.read(is);
    
    return is;
}

bool operator>(KmerFreq kmerFreq1, KmerFreq kmerFreq2){
    if(kmerFreq1.getFrequency() > kmerFreq2.getFrequency()){
        return true;
    }
    else if(kmerFreq1.getFrequency() == kmerFreq2.getFrequency() && kmerFreq1.getKmer().toString() > kmerFreq2.getKmer().toString()){
        return true;
    }
    else
        return false;
}

bool operator<(KmerFreq kmerFreq1, KmerFreq kmerFreq2){
    if(!(kmerFreq1 > kmerFreq2)){
        return true;
    }
    else
        return false;
}

bool operator==(KmerFreq kmerFreq1, KmerFreq kmerFreq2){
    if(kmerFreq1.getFrequency() == kmerFreq2.getFrequency()){
        return true;
    }
    else
        return false;
}

bool operator!=(KmerFreq kmerFreq1, KmerFreq kmerFreq2){
    if(!(kmerFreq1 == kmerFreq2)){
        return true;
    }
    else
        return false;
}

bool operator<=(KmerFreq kmerFreq1, KmerFreq kmerFreq2){
    if(kmerFreq1 < kmerFreq2 || kmerFreq1 == kmerFreq2){
        return true;
    }
    else
        return false;
}

bool operator>=(KmerFreq kmerFreq1, KmerFreq kmerFreq2){
    if(kmerFreq1 > kmerFreq2 || kmerFreq1 == kmerFreq2){
        return true;
    }
    else
        return false;
}
 