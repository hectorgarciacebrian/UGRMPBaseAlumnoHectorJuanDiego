/*
 * Metodología de la Programación: Kmer2
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 29 January 2023, 11:00
 */


#include <fstream>
#include <algorithm>

#include "Profile.h"

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";

Profile::Profile(){
    _profileId = "Unknown";
    _size = 0;
}

Profile::Profile(int size) : _size(size){
    _profileId = "Unknown";
    if (size < 0 || size >= DIM_VECTOR_KMER_FREQ)
    {
        throw out_of_range("El tamaño está fuera del rango permitido");
    }

    for (int i = 0; i < size; i++)
    {
        KmerFreq kmerFreq;
        _vectorKmerFreq[i] = kmerFreq;
    }
}

string Profile::getProfileId() const{
    return _profileId;
}

void Profile::setProfileId(std::string id){
    _profileId = id;
}

KmerFreq Profile::at(int index) const{
    if (index < 0 || index >= _size)
    {
        throw out_of_range("El tamaño está fuera del rango permitido");
    }

    return _vectorKmerFreq[index];
}

KmerFreq Profile::at(int index){
    if (index < 0 || index >= _size)
    {
        throw out_of_range("El tamaño está fuera del rango permitido");
    }

    return _vectorKmerFreq[index];
}

int Profile::getSize() const{
    return _size;
}

int Profile::getCapacity() const{
    return DIM_VECTOR_KMER_FREQ;
}


int Profile::findKmer(Kmer kmer, int initialPos, int finalPos){
    if (initialPos < 0 || initialPos >= _size || finalPos < initialPos || finalPos < 0 || finalPos > _size)
    {
        throw out_of_range("El tamaño está fuera del rango permitido");
    }

    for (int i = initialPos; i <= finalPos; i++)
    {
        if (_vectorKmerFreq[i].getKmer().toString() == kmer.toString())
        {
            return i;
        }
    }
    return -1;
}

int Profile::findKmer(Kmer kmer){

    for (int i = 0; i < this->getSize(); i++)
    {
        if (_vectorKmerFreq[i].getKmer().toString() == kmer.toString())
        {
            return i;
        }
    }
    return -1;
}

string Profile::toString(){
    string resultado = " ";
    resultado = resultado + _profileId + "\n";
    resultado = resultado + to_string(this->getSize()) + "\n";
    for (int i = 0; i < this->getSize(); i++)
    {
        resultado += _vectorKmerFreq[i].toString() + " ";
    }

    return resultado;
}

void Profile::sort(){
    for (int i = 0; i < this->getSize(); i++) {
        for (int j = 0; j < this->getSize() - i - 1; j++) {
            if (_vectorKmerFreq[j].getFrequency() < _vectorKmerFreq[j + 1].getFrequency()) 
            {
                swap(_vectorKmerFreq[j], _vectorKmerFreq[j + 1]);
            }
            else if (_vectorKmerFreq[j].getFrequency() == _vectorKmerFreq[j + 1].getFrequency())
            {
                if (_vectorKmerFreq[j].toString()[0] > _vectorKmerFreq[j + 1].toString()[0])
                {
                    swap(_vectorKmerFreq[j], _vectorKmerFreq[j + 1]);
                }  
            }  
        }
    }
}

void Profile::save(char fileName[]){
    ofstream file(fileName);
    if(!file.is_open()){
        throw ios_base::failure("Error opening file");
    }
    
    file << MAGIC_STRING_T << endl;
    file << _profileId << endl;
    file << _size << endl;
    
    for(int i=0; i<_size; i++){
        file << _vectorKmerFreq[i].getKmer().toString() << " " << endl;
    }
    
    file.close();
}

void Profile::load(char fileName[]){
    ifstream file(fileName);
    if(!file.is_open()){
        throw ios_base::failure("Error opening file");
    }
    
    string magic;
    getline(file, magic);
    if(magic != MAGIC_STRING_T){
        throw invalid_argument("Invalid magic string");
    }
    
    getline(file, _profileId);
    int newSize;
    file >> newSize;
    file.ignore();
    
    if(newSize > DIM_VECTOR_KMER_FREQ){
        throw out_of_range("Number of kmers in file exceeds maximum capacity");
        
    }
        
    if (newSize < 0) {
       throw out_of_range("Number of kmers in file is negative");
    }
    
    _size = newSize;

    for(int i=0; i<_size; i++){
        string kmerStr;
        int freq;
        file >> kmerStr >> freq;
        Kmer kmer(kmerStr);
        _vectorKmerFreq[i].setKmer(kmer);
        _vectorKmerFreq[i].setFrequency(freq);
    }
}

void Profile::append(KmerFreq kmerFreq){
    int index = -1;
    if(_size >= DIM_VECTOR_KMER_FREQ){
       throw std::out_of_range("Array capacity exceeded");
    }
    
    index = findKmer(kmerFreq.getKmer());
    
    if(index != -1){
        _vectorKmerFreq[index].setFrequency(_vectorKmerFreq[index].getFrequency() + kmerFreq.getFrequency());
    }else{
        _vectorKmerFreq[_size++] = kmerFreq;
    }
}

void Profile::normalize(string validNucleotides){
    for(int i=0; i<_size; i++){
        int index = -1;
        Kmer kmer = _vectorKmerFreq[i].getKmer();
        int freq = _vectorKmerFreq[i].getFrequency();
        
        kmer.normalize(validNucleotides);
        
        index = findKmer(kmer);
        if(index != -1){
            _vectorKmerFreq[index].setFrequency(_vectorKmerFreq[index].getFrequency() + freq);
            deletePos(i);
            i--;
        }else{
            _vectorKmerFreq[i].setKmer(kmer);
        }
    }
}

void Profile::deletePos(int pos){
    if(pos < 0 || pos >= _size){
        throw out_of_range("Invalid position");
    }
    
    for(int i=pos; i<_size-1; i++){
        _vectorKmerFreq[i] = _vectorKmerFreq[i+1];
    }
    
    _size--;
}

void Profile::zip(bool deleteMissing, int lowerBound){
    int newSize = 0;
    for(int i=0; i<_size; i++){
        if(!deleteMissing /*&& !array[i].getKmer().toString().find(UNKNOWN_NUCLEOTIDE)*/ || _vectorKmerFreq[i].getFrequency() > lowerBound){
            _vectorKmerFreq[newSize] = _vectorKmerFreq[i];
            newSize++;
        }
    }
    
    _size = newSize;
}

void Profile::join(Profile profile){
    for(int i=0; i<profile.getSize(); i++){
        append(profile.at(i));
    }
}

