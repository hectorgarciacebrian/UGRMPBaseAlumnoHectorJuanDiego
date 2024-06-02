/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Hector García Cebrian
 * @author Juan Diego Martín Payán
 * 
 * 
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

bool operator>(const KmerFreq &kmerFreq1, const KmerFreq &kmerFreq2){
    return (kmerFreq1.getFrequency() > kmerFreq2.getFrequency());
}

bool operator<(const KmerFreq &kmerFreq1, const KmerFreq &kmerFreq2){
    return (kmerFreq1.getFrequency() < kmerFreq2.getFrequency());
}

bool operator==(const KmerFreq &kmerFreq1, const KmerFreq &kmerFreq2){
    if(kmerFreq1.getFrequency() == kmerFreq2.getFrequency()){
        if(kmerFreq1.getKmer().toString() == kmerFreq2.getKmer().toString()){
            return true;
        }
    }
    
    return false;
}

bool operator!=(const KmerFreq &kmerFreq1, const KmerFreq &kmerFreq2){
    if(kmerFreq1.getFrequency() != kmerFreq2.getFrequency()){
        if(kmerFreq1.getKmer().toString() != kmerFreq2.getKmer().toString()){
            return true;
        }
    }
    
    return false;
}

bool operator<=(const KmerFreq &kmerFreq1, const KmerFreq &kmerFreq2){
    return (kmerFreq1.getFrequency() <= kmerFreq2.getFrequency());
}

bool operator>=(const KmerFreq &kmerFreq1, const KmerFreq &kmerFreq2){
    return (kmerFreq1.getFrequency() >= kmerFreq2.getFrequency());
}

ostream& operator<<(ostream& os, KmerFreq& kmerFreq) {
    os << kmerFreq.toString();
    return os;
}

istream& operator>>(istream& is, KmerFreq& kmerFreq) {
    string kmerText;
    int frequency;

    is >> kmerText >> frequency;

    if (!is) {
        // Por si hay algún error al leer
        return is;
    }

    Kmer kmer(kmerText);

    kmerFreq.setKmer(kmer);
    kmerFreq.setFrequency(frequency);

    return is;
}

