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
#include <iostream>
using namespace std;



void ReadArrayKmerFreq(const KmerFreq array[], const int dim, int nElements){
      for (int i = 0; i < dim; i++)
      {
            if (i = 0) //Cambiar esto, condicion: si hay un kmer en la posicion dada.
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



void NormalizeArrayKmerFreq(KmerFreq array[], int nElements, string validNucleotides){ 
    
    // Loop to traverse and normalize each one of the kmers in array
          // Normalize kmer i
    
    
    // Loop to traverse the kmers in array from position 1 to position nElements-1
          // index = Position of array[i].getKmer() in the subarray that begins
          //         at position 0 and ends at position i-1
          // If array[i].getKmer() was found in the the subarray from 0 to i-1 
               // Accumulate the frequencies of the kmers at positions 
               //    index and i in the kmer at position index
               // Delete from the array, the kmer at position i 
}
