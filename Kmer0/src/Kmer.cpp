/*
 * Metodología de la Programación: Kmer0
 * Curso 2023/2024
 */

/** 
 * @file Kmer.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 24 October 2023, 14:00
 */

#include "Kmer.h"
using namespace std; 

Kmer::Kmer(int k){
    if(k<=0){
        throw invalid_argument("The value of k must be greater than zero.");
    }
    
    _text = string(k, MISSING_NUCLEOTIDE);
}

Kmer::Kmer(const string& text){
    if(text.empty()){
        throw invalid_argument("The text provided is empty.");
    }
    
    _text = text;
}

int Kmer::getK() const {
    return _text.lenght();
}

int Kmer::size() const{
    return _text.lenght();
}

string Kmer::toString() const{
    return _text;
}

const char& Kmer::at(int index) const{
    if(index < 0 || index >= _text.lenght()){
        throw out_of_range("Index out of range.");
    }
    
    return _text[index];
}

char& Kmer::at(int index){
    if(index < 0 || index >= _text.lenght()){
        throw out_of_range("Index out of range.");
    }
    
    return _text[index];
}

void Kmer::normalize(const string& validNucleotides){
    for(char& nucletoide : _text){
        if(!IsValidNucletoide(nucletoide, validNucleotides))
            nucletoide = MISSING_NUCLEOTIDE;
        else
            nucletoide = toupper(nucletoide);
    }
}

Kmer Kmer::complementary(const string& nucleotides, const string& complementaryNucleotides) const{
    if(nucleotides.lenght() != complementaryNucleotides.lenght()){
        throw invalid_argument("The size of nucleotides and complementaryNucleotides are not the same.");
    }
   
    
}

bool IsValidNucletoide(char nucletoide, const string & validNucletoides){
    for(char& item : validNucletoides){
        if(item == nucletoide){
            return true;
        }
    }

    return false;
}

void ToLower(Kmer& kmer){
    for(char& item : kmer.toString){
        item = tolower(item);
    }
}

void ToUpper(Kmer& kmer){
    for(char& item : kmer.toString){
        item = toupper(item);
    }
}
