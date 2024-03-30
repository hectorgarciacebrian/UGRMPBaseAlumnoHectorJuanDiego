/*
 * Metodología de la Programación: Kmer1
 * Curso 2023/2024
 */

/** 
 * @file ArrayKmerFreqFunctions.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 27 October 2023, 12:00
 */


#include "ArrayKmerFreqFunctions.h"
#include "string.h"
#include <string>
#include <iostream>
using namespace std;



void ReadArrayKmerFreq(const KmerFreq array[], const int dim, int nElements){
      for (int i = 0; i < dim; i++)
      {
            if (!array[i].getKmer().toString().empty()) //Cambiar esto, condicion: si hay un kmer en la posicion dada.
            {
                 nElements++;
            }
            
      }
      
}

void PrintArrayKmerFreq(KmerFreq array[], int nElements){
      cout << "El número de elementos es: " << nElements << "\n";

      for (int i = 0; i < nElements; i++)
      {
            cout << array[i].toString() << "\n";
      }
      
}

void SwapElementsArrayKmerFreq(KmerFreq array[], int nElements, int first, int second){
      if ((first < 0 || first >= nElements) || (second < 0 || second >= nElements))
      {
            throw out_of_range("Valor incorrecto de frecuencia.");
      }
      else
      {
            KmerFreq aux = array[first];
            array[first] = array[second];
            array[second] = aux;

      }
}

int FindKmerInArrayKmerFreq(KmerFreq array[], Kmer kmer, int initialPos, int finalPos){
    for(int i=initialPos; i<=finalPos; i++){
        if(array[i].getKmer().toString() == kmer.toString()){
            return i;
        }
    }
    
    return -1;
}

void SortArrayKmerFreq(KmerFreq array[], int nElements){
    for(int i = 0; i < nElements; i++){
        for(int j = i + 1; j < nElements; j++){
            if (array[j].getFrequency() < array[i].getFrequency() || (array[i].getFrequency() == array[j].getFrequency() && array[i].getKmer().toString() > array[j].getKmer().toString())) {
                SwapElementsArrayKmerFreq(array, nElements, i, j);
            }
        }
    }
}

void NormalizeArrayKmerFreq(KmerFreq array[], int nElements, string validNucleotides){ 
    
    // Loop to traverse and normalize each one of the kmers in array
    for(int i=0; i<nElements; i++){
        Kmer currentKmer(array[i].getKmer());
        currentKmer.normalize(validNucleotides);
        array[i].setKmer(currentKmer.toString());
    }
    
    for(int i=1; i<nElements; i++){
        int index = -1;
        bool encontrado = false;
        
        for(int j=0; j<i && !encontrado; j++){
            if(array[i].getKmer().toString() == array[j].getKmer().toString()){
                index = j;
                encontrado = true;
            }
        }
        
        if(index!=-1){
            array[index].setFrequency(array[index].getFrequency()+array[i].getFrequency());
            for(int k=i; k<nElements-1; k++){
                array[k] = array[k+1];
            }
            
            nElements--;
            i--;
        }
    } 
}

void DeletePosArrayKmerFreq(KmerFreq array[], int nElements, int pos){
    if(pos < 0 || pos >= nElements){
        throw out_of_range("El elemento a borrar no se encuntra en el array");
    }
    
    for(int i=pos; i<nElements-1; i++){
        array[i] = array[i+1];
    }
    
    nElements--;
}void ZipArrayKmerFreq(KmerFreq array[], int nElements, bool deleteMissing, int lowerBound){
    int cont=0;
    
    for(int i=0; i<nElements; i++){
        if(!deleteMissing /*&& !array[i].getKmer().toString().find(UNKNOWN_NUCLEOTIDE)*/ || array[i].getFrequency() > lowerBound){
            array[cont] = array[i];
            cont++;
        }
    }
    
    nElements = cont;
}