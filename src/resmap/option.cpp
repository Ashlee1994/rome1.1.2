/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include "./mpi.h"
#include "option.h"

static inline std::string trimKey(std::string key) {
    if (key.find('-') == std::string::npos) return key;
    else return key.substr(key.find_last_of('-'));
}

// add all options
void Option::addOption(std::string key,std::string comment,std::string default_value){
    Elem oneOption;
    key = trimKey(key);
    keyLength = keyLength < key.length()?key.length():keyLength;
    int current_commentLength = 0;// consider comment has '\n'
    for (const auto & p : comment) {
        current_commentLength++;
        if (p=='\n') {
            commentLength = commentLength < current_commentLength?current_commentLength:commentLength;
            current_commentLength = 0;
        }
    }
    commentLength = commentLength < current_commentLength?current_commentLength:commentLength;
    oneOption.key = key;oneOption.value = default_value;oneOption.comment = comment;
    AllOptions.push_back(oneOption);
    IOParser[key] = default_value;
}

// read command line
void Option::readCommandLine(int argc, char * argv[]){
    // add options one by one
    for (int i = 1; i < argc;) {
        auto key = trimKey(argv[i+0]);
        if (IOParser.find(key) == IOParser.end()) { // if it does not have this option,escape it
//            MPI_LOG<<"wrong option(command line argument) ,don't need "<<key<<" option."<<std::endl;
            i++;
        }
        else{ // if it has this option
            if (i+1 < argc) {
                std::string value = argv[i+1];
                if (IOParser.find(trimKey(value)) == IOParser.end()) { // if this option has value
                    IOParser[key] = value;
                    i += 2;
                }
                else{
                    IOParser[key] = "1";
                    i += 1;
                }
            }
            else{
                IOParser[key] = "1";
                i += 1;
            }
        }
    }
}

std::string Option::getOption(std::string key){
    // trim the key to format : "-v"
    key = trimKey(key);
    if (IOParser.find(key) == IOParser.end()) {
        std::cerr<<"make sure you have add "<<key<<" option."<<std::endl;
        ERROR_REPORT("missing some key in option.")
    }
    std::string value = IOParser[key];
    if (value == Option::unspecified()) {
        std::cerr<<"you need input "<<key<<" option , which means : "<<std::endl;
        std::for_each(AllOptions.begin(), AllOptions.end(),
                      [&](Elem& elem){if(elem.key == key) MPI_LOG<<elem.comment<<std::endl;});
        ERROR_REPORT("missing some input parameter.")
    }
    return value;
}

// read option
int Option::getIntOption(std::string key){
    std::string value = getOption(key);
    return atoi(value.c_str());
}

double Option::getFloatOption(std::string key){
    std::string value = getOption(key);
    float retval;
    sscanf(value.c_str(), "%f", &retval);
    return retval;
}

std::string Option::getStrOption(std::string key){
    std::string value = getOption(key);
    return value;
}

bool Option::getBoolOption(std::string key){
    std::string value = getOption(key);
    return bool(atoi(value.c_str()));
}

void Option::printValue(){
    MPI_LOG<<"--------------------------------  user's options    ------------------------------------------------------------"<<std::endl;
    std::for_each(AllOptions.begin(), AllOptions.end(),
                  [&](Elem& elem){
                      MPI_LOG<<std::setw(keyLength)<<std::left<<elem.key<<" , ";
                      MPI_LOG<<std::setw(commentLength)<<elem.comment;
                      MPI_LOG<<" , "<<IOParser[elem.key]<<std::endl;});
    MPI_LOG<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
}

void Option::printHelp(){
    MPI_LOG<<"----------------------  general option(need to set)      -------------------------------------------------------"<<std::endl;
    std::for_each(AllOptions.begin(), AllOptions.end(),
                  [&](Elem& elem){if(elem.value == Option::unspecified())
                      std::cerr<<std::setw(keyLength)<<std::left<<elem.key<<" , "<<elem.comment<<std::endl;});
    MPI_LOG<<"--------------------  advanced option(with default value)  -----------------------------------------------------"<<std::endl;
    std::for_each(AllOptions.begin(), AllOptions.end(),
                  [&](Elem& elem){if(elem.value != Option::unspecified())
                      std::cerr<<std::setw(keyLength)<<elem.key<<" , "<<std::setw(commentLength)<<elem.comment<<" , default : "<<elem.value<<std::endl;});
    MPI_LOG<<"----------------------------------------------------------------------------------------------------------------"<<std::endl;
}

