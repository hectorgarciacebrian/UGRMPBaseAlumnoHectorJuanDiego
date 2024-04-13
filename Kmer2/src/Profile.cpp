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


#include "Profile.h"

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";

Profile::Profile(){
    _profileId = "Unknown";
    _size = 0;
    _vectorKmerFreq[DIM_VECTOR_KMER_FREQ];
}

Profile::Profile(int size) : _size(size){
    _profileId = "Unknown";
    if (size > 0 || size < DIM_VECTOR_KMER_FREQ)
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
    if (index > 0 || index < _size)
    {
        throw out_of_range("El tamaño está fuera del rango permitido");
    }

    return _vectorKmerFreq[index];
}

KmerFreq Profile::at(int index){
    if (index > 0 || index < _size)
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
    if (initialPos < 0 || initialPos > _size || finalPos < initialPos || finalPos < 0 || finalPos > _size)
    {
        throw out_of_range("El tamaño está fuera del rango permitido");
    }

    for (int i = initialPos; i < finalPos; i++)
    {
        if (_vectorKmerFreq[i].getKmer().toString() == kmer.toString())
        {
            return i;
        }
    }
    return -1;
}

int Profile::findKmer(Kmer kmer){

    for (int i = 0; i < this->getSize() - 1; i++)
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
    for (int i = 0; i < this->getSize() - 1; i++)
    {
        resultado = resultado + this->at(i).toString() + " ";
    }

    return resultado;
}

void Profile::sort(){
    for (int i = 0; i < this->getSize() - 1; ++i) {
        for (int j = 0; j < this->getSize() - i - 1; ++j) {
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

