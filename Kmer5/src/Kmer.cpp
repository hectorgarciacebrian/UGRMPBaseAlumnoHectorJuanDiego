

/** 
 * @file Kmer.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 24 October 2023, 14:00
 */

#include "../include/Kmer.h"
#include <string>
using namespace std; 

Kmer::Kmer(int k){
    if(k<1){
        throw invalid_argument(string("Kmer(int k): ") + 
                "invalid length " + to_string(k));
    }
    
    this->_text = string(k, MISSING_NUCLEOTIDE);
}

Kmer::Kmer(const string& text){
    if(text.empty()){
        throw invalid_argument("The text provided is empty.");
    }
    
    _text = text;
}

int Kmer::size() const{
    return _text.length();
}

int Kmer::getK() const {
    return size();
}

string Kmer::toString() const{
    return _text;
}

const char& Kmer::at(int index) const{
    if(index < 0 || index >= _text.length()){
        throw out_of_range("Index out of range.");
    }
    
    return _text[index];
}

char& Kmer::at(int index){
    if(index < 0 || index >= _text.length()){
        throw out_of_range("Index out of range.");
    }
    
    return _text[index];
}

void Kmer::normalize(const string& validNucleotides){
    for(char& nucletoide : _text){
        for(int i=0; i<validNucleotides.length(); i++){
            if(nucletoide != validNucleotides[i]){
                nucletoide = MISSING_NUCLEOTIDE;
            }
            else{
                nucletoide = toupper(nucletoide);
            }
        }
    }
}

Kmer Kmer::complementary(const string& nucleotides, const string& complementaryNucleotides) const{
    if(nucleotides.length() != complementaryNucleotides.length()){
        throw invalid_argument("The size of nucleotides and complementaryNucleotides are not the same.");
    }
    
    Kmer result(*this);
    
    for(char& nucletoide : result._text){
        size_t idx = nucleotides.find(toupper(nucletoide));
        if(idx != string::npos){
            nucletoide = complementaryNucleotides[idx];
        }
    }
    
    return result;
    
}

bool IsValidNucletoide(char nucletoide, const string & validNucletoides){
    for(int i=0; i<validNucletoides.length(); i++){
        if(nucletoide == validNucletoides[i]){
            return true;
        }
    }
    
    return false;
}

void ToLower(Kmer& kmer){
    for(char& item : kmer.toString()){
        item = tolower(item);
    }
}

void ToUpper(Kmer& kmer){
    for(char& item : kmer.toString()){
        item = toupper(item);
    }
}

ostream& operator<<(ostream& os, const Kmer& kmer) {
    os << kmer.toString();
    return os;
}

istream& operator>>(istream& is, Kmer& kmer) {
    string input;
    is >> input;
    kmer = Kmer(input); // Reasignar el objeto kmer
    return is;
}

