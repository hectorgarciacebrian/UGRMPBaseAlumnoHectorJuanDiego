/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerCounter.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "KmerCounter.h"
#include <string>
#include <fstream>
#include <stdexcept>

using namespace std;

/**
 * DEFAULT_VALID_NUCLEOTIDES is a c-string that contains the set of characters
 * that will be considered as valid nucleotides. 

 * The constructor of the class KmerCounter uses this c-string as a 
 * default parameter. It is possible to use a different c-string if that
 * constructor is used with a different c-string
 */
const char* const KmerCounter::DEFAULT_VALID_NUCLEOTIDES="ACGT";
//const string KmerCounter::_allNucleotides = "ACGT"; // Initializing _allNucleotides

KmerCounter::KmerCounter(int k, string validNucleotides) {
    _k = 5; //Estos valores se asignan según el guión de prácticas.
    _validNucleotides = DEFAULT_VALID_NUCLEOTIDES; //Estos valores se asignan según el guión de prácticas.
    _frequency = new int*[_k];
    for (int i = 0; i < _k; ++i) {
        _frequency[i] = new int[_k]();
    }
}

KmerCounter::KmerCounter(const KmerCounter &orig) {
    _k = orig._k;
    _validNucleotides = orig._validNucleotides;
    _frequency = new int*[_k];
    for (int i = 0; i < _k; ++i) {
        _frequency[i] = new int[_k];
    }
    for (int i = 0; i < _k; i++) {
        for (int j = 0; j < _k; j++) {
            _frequency[i][j] = orig._frequency[i][j];
        }
    }
}

KmerCounter::~KmerCounter() {
    for (int i = 0; i < _k; ++i) {
        delete[] _frequency[i];
    }
    delete[] _frequency;
    _frequency = nullptr;
    _validNucleotides = " ";
}

int KmerCounter::getNumNucleotides() {
    return _allNucleotides.size();
}

int KmerCounter::getK() {
    return _k;
}

int KmerCounter::getNumKmers() {
    int n = 1;
    for (int i = 0; i < _k; ++i)
        n *= _k;
    return n;
}

int KmerCounter::getNumRows() const {
    return _k;
}

int KmerCounter::getNumCols() const {
    return _k;
}

int KmerCounter::getNumberActiveKmers() {
    int n = 0;
    for (int i = 0; i < this->getNumRows(); ++i) {
        for (int j = 0; j < this->getNumCols(); ++j) {
            if (_frequency[i][j] > 0)
                n++;
        }
    }
    return n;
}

string KmerCounter::toString() {
    string outputString = _allNucleotides + " " + to_string(_k) + "\n";
    for (int row = 0; row < getNumRows(); row++) {
        for (int col = 0; col < getNumCols(); col++) {
            outputString += to_string((*this)(row, col)) + " ";
        }
        outputString += "\n";
    }
    return outputString;
}

void KmerCounter::getRowColumn(const Kmer &kmer, int &row, int &column) const {
    string kmerStr = kmer.toString();
    if (kmer.getK() != _k) {
        row = -1;
        column = -1;
        return;
    }

    int rowLength = _k / 2 + _k % 2;
    //int columnLength = _k / 2;
    
    string rowPart, columnPart;

    for (int i = 0; i < rowLength; ++i) {
        rowPart += kmerStr[i];
    }

    for (int i = rowLength; i < _k; ++i) {
        columnPart += kmerStr[i];
    }

    row = getIndex(rowPart);
    column = getIndex(columnPart);

    if (row == -1 || column == -1) {
        row = -1;
        column = -1;
    }
}

void KmerCounter::increaseFrequency(Kmer kmer, int frequency) {
    frequency = 1; //Valor asignado según el guión de prácticas.
    string kmerStr = kmer.toString();
    for (char nucleotido : kmerStr) {
        bool valido = false;
        for (char validNucleotide : _allNucleotides) {
            if (nucleotido == validNucleotide) {
                valido = true;
                break;
            }
        }
        if (!valido) {
            throw invalid_argument("Este Kmer contiene nucleótidos inválidos");
        }
    }
    
    int row, column;
    getRowColumn(kmer, row, column);

    if (row == -1 || column == -1) {
        throw invalid_argument("Este Kmer contiene nucleótidos no válidos en posiciones inválidas de la matriz.");
    }

    _frequency[row][column] += frequency;
}

KmerCounter& KmerCounter::operator=(const KmerCounter &orig) {
    if (this != &orig) {
        for (int i = 0; i < _k; ++i) {
            delete[] _frequency[i];
        }
        delete[] _frequency;
        
        _k = orig._k;
        _validNucleotides = orig._validNucleotides;
        _frequency = new int*[_k];
        for (int i = 0; i < _k; ++i) {
            _frequency[i] = new int[_k];
        }
        for (int i = 0; i < _k; i++) {
            for (int j = 0; j < _k; j++) {
                _frequency[i][j] = orig._frequency[i][j];
            }
        }
    }
    return *this;
}

KmerCounter& KmerCounter::operator+=(const KmerCounter &kc) {
    if (_validNucleotides != kc._validNucleotides || _allNucleotides != kc._allNucleotides) {
        throw invalid_argument("Los objetos tienen diferentes sets de nucleotidos.");
    }

    if (_k != kc._k) {
        throw invalid_argument("Los dos objetos tienen valores diferentes de K.");
    }
    
    int numRows = getNumRows();
    int numCols = getNumCols();

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            _frequency[i][j] += kc._frequency[i][j];
        }
    }

    return *this;
}

void KmerCounter::initFrequencies() {
    int numRows = getNumRows();
    int numCols = getNumCols();
    _frequency = new int*[numRows];
    for (int i = 0; i < numRows; ++i) {
        _frequency[i] = new int[numCols];
    }
    
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            _frequency[i][j] = 0;
        }
    }
}

void KmerCounter::calculateFrequencies(char* fileName) {
    initFrequencies();

    ifstream file(fileName);
    if (!file.is_open()) {
        throw ios_base::failure("Error al abrir el archivo");
    }

    string secuencia;
    while (file >> secuencia) {
        for (char &c : secuencia) {
            bool valido = false;
            for (char validNucleotide : _validNucleotides) {
                if (c == validNucleotide) {
                    valido = true;
                    break;
                }
            }
            if (!valido) {
                c = Kmer::MISSING_NUCLEOTIDE;
            }
        }

        for (size_t i = 0; i <= secuencia.length() - _k; ++i) {
            string kmerStr;
            for (int j = 0; j < _k; ++j) {
                kmerStr.push_back(secuencia[i + j]);
            }
            Kmer kmer(kmerStr);
            int row, col;
            getRowColumn(kmer, row, col);
            if (row != -1 && col != -1) {
                _frequency[row][col]++;
            }
        }
    }

    file.close();
}

Kmer KmerCounter::getKmer(int row, int column) {
    if (row < 0 || row >= getNumRows() || column < 0 || column >= getNumCols()) {
        throw invalid_argument("Coordenadas inválidas");
    }
    
    string rowKmer = getInvertedIndex(row, _k / 2 + _k % 2);
    string colKmer = getInvertedIndex(column, _k / 2);
    
    string kmerString = rowKmer + colKmer;
    
    Kmer kmer(kmerString);

    return kmer;
}

Profile KmerCounter::toProfile() {
    Profile profile;

    for (int row = 0; row < getNumRows(); ++row) {
        for (int column = 0; column < getNumCols(); ++column) {
            int frequency = _frequency[row][column];
            if (frequency > 0) {
                Kmer kmer = getKmer(row, column);
                KmerFreq kmerFreq;
                kmerFreq.setKmer(kmer);
                kmerFreq.setFrequency(frequency);
                profile.append(kmerFreq);
            }
        }
    }

    return profile;
}

int KmerCounter::operator()(int row, int column) const {
    return _frequency[row][column];
}

int KmerCounter::operator()(int row, int column) {
    return _frequency[row][column];
}

int KmerCounter::getIndex(const std::string &kmer) const {
    int index = 0;
    int base = 1;

    for (size_t i = 0; i < kmer.size(); i++) {
        size_t pos = _allNucleotides.find(kmer[kmer.size() - i - 1]);
        if (pos == string::npos)
            return -1;
        index += pos * base;
        base *= _allNucleotides.size();
    }
    return index;
}

string KmerCounter::getInvertedIndex(int index, int nCharacters) const {
    string result(nCharacters, Kmer::MISSING_NUCLEOTIDE);

    for (int i = result.size(); i > 0; i--) {
        result[i - 1] = _allNucleotides[index % _allNucleotides.size()];
        index = index / _allNucleotides.size();
    }
    return result;
}
