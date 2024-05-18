/*
 * Metodología de la Programación: Kmer4
 * Curso 2023/2024
 */

/* 
 * File:   main.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 *
 * Created on 17 November 2023, 12:45
 */

#include <iostream>

#include "Profile.h"

using namespace std;


/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(ostream& outputStream) {
    outputStream << "ERROR in Kmer4 parameters" << endl;
    outputStream << "Run with the following parameters:" << endl;
    outputStream << "kmer4 [-t min|max] <file1.prf> <file2.prf> [ ... <filen.prf>]" << endl;
    outputStream << endl;
    outputStream << "Parameters:" << endl;
    outputStream << "-t min | -t max: search for minimun distances or maximum distances (-t min by default)" << endl;
    outputStream << "<file1.prf>: source profile file for computing distances" << endl;
    outputStream << "<file2.prf> [ ... <filen.prf>]: target profile files for computing distances" << endl;  
    outputStream << endl;
    outputStream << "This program computes the distance from profile <file1.prf> to the rest" << endl;
    outputStream << endl;
}

/**
 * This program reads an undefined number of Profile objects from the set of 
 * files passed as parameters to main(). All the Profiles object, except the 
 * first one, must be stored in a dynamic array of Profile objects. Then, 
 * for each Profile in the dynamic array, this program prints to the 
 * standard output the name of the file of that Profile and the distance from 
 * the first Profile to the current Profile. 
 * Finally, the program should print in the standard output, the name of 
 * the file with the Profile with the minimum|maximum  distance to the Profile 
 * of the first file and its profile identifier.
 * 
 * At least, two Profile files are required to run this program.
 * 
 * This program assumes that the profile files are already normalized and 
 * sorted by frequency. This is not checked in this program. Unexpected results
 * will be obtained if those conditions are not met.
 * 
 * Running sintax:
 * > kmer4 [-t min|max] <file1.prf> <file2.prf> [  ... <filen.prf>] 
 * 
 * Running example:
 * > kmer4 ../Genomes/human1.prf ../Genomes/worm1.prf ../Genomes/mouse1.prf 
Distance to ../Genomes/worm1.prf: 0.330618
Distance to ../Genomes/mouse1.prf: 0.224901
Nearest profile file: ../Genomes/mouse1.prf
Identifier of the nearest profile: mus musculus
 * 
 * Running example:
 * > kmer4 -t max ../Genomes/human1.prf ../Genomes/worm1.prf ../Genomes/mouse1.prf 
Distance to ../Genomes/worm1.prf: 0.330618
Distance to ../Genomes/mouse1.prf: 0.224901
Farthest profile file: ../Genomes/worm1.prf
Identifier of the farthest profile: worm
 */
int main(int argc, char* argv[]) {
    
    // Verificar que se proporcionen al menos 2 archivos Profile
    if (argc < 3)
    {
        cerr << "Se necesitan al menos 2 archivos Profile\n";
        return -1;
    }
    
    // Cargar el primer perfil desde el primer archivo
    Profile firstProfile;
    try {
        firstProfile.load(argv[1]);
    } catch (const ios_base::failure& e) {
        cerr << "Error al cargar el primer perfil: " << e.what() << endl;
        return -1;
    }
    
    // Cargar el resto de los perfiles en un array dinámico
    Profile* profileArray = new Profile[argc - 2];
    for (int i = 2; i < argc; ++i)
    {
        try {
            profileArray[i - 2].load(argv[i]);
        } catch (const ios_base::failure& e) {
            cerr << "Error al cargar el perfil " << argv[i] << ": " << e.what() << endl;
            // Liberar la memoria antes de salir
            delete[] profileArray;
            return -1;
        }
    }
    
    // Calcular y mostrar la distancia desde el primer perfil al resto
    for (int i = 0; i < argc - 2; ++i)
    {
        double distance = firstProfile.getDistance(profileArray[i]);
        cout << "Distancia a " << argv[i + 2] << ": " << distance << endl;
    }
    
    // Encontrar el perfil con la distancia mínima o máxima
    double minDistance = firstProfile.getDistance(profileArray[0]);
    double maxDistance = firstProfile.getDistance(profileArray[0]);
    string minProfileName = argv[2];
    string maxProfileName = argv[2];
    for (int i = 1; i < argc - 2; ++i)
    {
        double distance = firstProfile.getDistance(profileArray[i]);
        if (distance < minDistance)
        {
            minDistance = distance;
            minProfileName = argv[i + 2];
        }
        if (distance > maxDistance)
        {
            maxDistance = distance;
            maxProfileName = argv[i + 2];
        }
    }

    // Mostrar el nombre del archivo y el identificador del perfil con la distancia mínima o máxima
    cout << "Perfil más cercano: " << minProfileName << " (distancia: " << minDistance << ")" << endl;
    cout << "Perfil más lejano: " << maxProfileName << " (distancia: " << maxDistance << ")" << endl;

    // Liberar el array dinámico de Profile
    delete[] profileArray;
    
    return 0;
}

