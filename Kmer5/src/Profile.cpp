/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include <fstream>

#include "Profile.h"

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";
const string Profile::MAGIC_STRING_B="MP-KMER-B-1.0";

void Profile::allocate(int capacity){
    _vectorKmerFreq = new KmerFreq[capacity];
}
void Profile::deallocate(){
    if(_vectorKmerFreq != nullptr){
        delete[] _vectorKmerFreq;
        _vectorKmerFreq = nullptr;
    }
}

void Profile::reallocate(int newCapacity){
    deallocate();
    allocate(newCapacity);
}

void Profile::copy(const Profile &otro){
    _profileId = otro._profileId;
    _size = otro._size;
    _capacity = otro._capacity;
    
    deallocate();
    allocate(_capacity);
    
    for(int i=0; i<_size; i++){
        _vectorKmerFreq[i] = otro._vectorKmerFreq[i];
    }
}

Profile::Profile() : _profileId("Unknown"), _size(0), _capacity(INITIAL_CAPACITY){
    allocate(_capacity);
}

Profile::Profile(int size) : _size(size){
    _profileId = "Unknown";
    _capacity = size;
    if (size < 0)
    {
        throw out_of_range("El tamaño está fuera del rango permitido");
    }

    _vectorKmerFreq = new KmerFreq[size];
    for (int i = 0; i < _size; i++)
    {
        KmerFreq kmerFreq;
        _vectorKmerFreq[i] = kmerFreq;
    }
}

Profile::~Profile(){
    deallocate();
}

Profile& Profile::operator=(const Profile &orig){
    if(this!=&orig){
        copy(orig);
    }
    
    return *this;
}

Profile::Profile(const Profile &orig){
    copy(orig);
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
    return _capacity;
}

double Profile::getDistance(Profile otherProfile){
    double dist = 0.0;
    double distancia_final = 0.0;
    if( this->getSize() == 0 || otherProfile.getSize() == 0){
        throw invalid_argument("El objeto actual o el otro obheto están vacias");
    }
    
    for(int i=0; i<_size; i++){
        int pos1 = i;
        Kmer kmer1 = _vectorKmerFreq[i].getKmer();
        int pos2 = otherProfile.findKmer(kmer1);
        if(pos2 == -1){
            pos2 = 15;
        }
        
        double rank = abs(pos1 - pos2);
        dist += rank;
    }
    
    distancia_final = dist/(_size*otherProfile.getSize());
    
    return distancia_final;
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
            if (_vectorKmerFreq[j] < _vectorKmerFreq[j + 1]) 
            {
                swap(_vectorKmerFreq[j], _vectorKmerFreq[j + 1]);
            }
            else if (_vectorKmerFreq[j] == _vectorKmerFreq[j + 1])
            {
                if (_vectorKmerFreq[j] > _vectorKmerFreq[j + 1])
                {
                    swap(_vectorKmerFreq[j], _vectorKmerFreq[j + 1]);
                }  
            }  
        }
    }
}

void Profile::save(char fileName[], char mode){
    ofstream file(fileName);
    if(mode == 'b'){
        file.open(fileName, ios::binary);
        file << MAGIC_STRING_B << endl;
        file << _profileId << endl;
        file << _size << endl;
        
        for(int i=0; i<_size; i++){
            file.write(_vectorKmerFreq[i].toString(), _vectorKmerFreq[i].getKmer().size());
        }
    }
    else if(mode == 't'){
        file.open(fileName, ios::out);
        file << MAGIC_STRING_T << endl;
        file << _profileId << endl;
        file << _size << endl;
        
        for(int i=0; i<_size; i++){
            file << _vectorKmerFreq[i].toString();
        }
    }
    else{
        throw invalid_argument("El modo seleccionado no esta dentro de los permitidos");
    }
    
    if(!file.is_open()){
        throw ios_base::failure("Error opening file");
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
    if(magic != MAGIC_STRING_T && magic != MAGIC_STRING_B){
        throw invalid_argument("Invalid magic string");
    }
    
    getline(file, _profileId);
    int newSize;
    file >> newSize;
    file.ignore();
    
    if(newSize >= _capacity){
        reallocate(_capacity + BLOCK_SIZE);
        
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
    if(_size >= _capacity){
        reallocate(_capacity + BLOCK_SIZE);
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

void Profile::deprecated(Profile profile){
    for(int i=0; i<profile.getSize(); i++){
        append(profile.at(i));
    }
}

KmerFreq Profile::operator[](int index) const{
    return _vectorKmerFreq[index];
}

KmerFreq Profile::operator[](int index){
    return _vectorKmerFreq[index];
}

Profile Profile::operator+=(KmerFreq kmerFreq){
    bool encontrado = false;
    
    for(int i=0; i<_size; i++){
        if(_vectorKmerFreq[i] == kmerFreq){
            int frecuencia = _vectorKmerFreq[i].getFrequency() + kmerFreq.getFrequency();
            _vectorKmerFreq[i].setFrequency(frecuencia);
            encontrado = true;
        }
    }
    
    if(!encontrado){
        if(_size >= _capacity){
            reallocate(_capacity + BLOCK_SIZE);
        }
        
        _vectorKmerFreq[size++] = kmerFreq;
    }
    
    return *this;
}

Profile Profile::operator+=(Profile profile){
    for(int j=0; j<profile.getSize(); j++){
        bool encontrado = false;
        for(int i=0; i<_size; i++){
            if(_vectorKmerFreq[i].getKmer() == profile._vectorKmerFreq[j].getKmer()){
                int frecuencia = _vectorKmerFreq[i].getFrequency() + profile._vectorKmerFreq[j].getFrequency();
                _vectorKmerFreq[i].setFrequency(frecuencia);
                encontrado = true;
            }
        }
        if(!encontrado){
            if(_size>=_capacity){
                reallocate(_capacity + BLOCK_SIZE);
            }

            _vectorKmerFreq[size++] = profile._vectorKmerFreq[j];
        }
    }
    
    return *this;
}

ostream operator<<(ostream os, Profile &profile){
    os << profile.toString();
    
    return os;
}


istream operator>>(istream is, Profile &profile){
    profile.d
}