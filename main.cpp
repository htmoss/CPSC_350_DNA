#include<string>
#include<iostream>
#include<fstream>
#include "dna.h"
#include<cstdlib>

using namespace std;

int main(int argc, char** argv)
{
  char answer = 'Y';
  Dna dna;
  //loop so that the user can process multiple lists
  while (answer == 'Y' or answer == 'y')
  {
    dna.clearVariables();
    dna.readFile();
    dna.printSumStats();
    cout<<endl<<"Type 'Y' if you want to process again\n";
    cin>>answer;
  }
  return 0;
}
