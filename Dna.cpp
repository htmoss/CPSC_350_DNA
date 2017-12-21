#include<string>
#include<iostream>
#include "dna.h"
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<cmath>
using namespace std;

//constructor
Dna::Dna()
{

}
//deconstructor
Dna::~Dna()
{
  cout<<"Object Deleted"<<endl;
}
//clears the private variables from last calculations
void Dna::clearVariables()
{
  sum = 0;
  dnaCount = 0;
  varianceSum = 0;
  finalVariance = 0;
  standDev = 0;
  aCount = 0;
  cCount = 0;
  gCount = 0;
  tCount = 0;
}
//initially reads the file, counts sum of words and number of nucleotides
void Dna::readFile()
{
  cout<<"\nPlease enter the name of the file you would like to process:\n";
  cin>>inputFile;
  ifstream dnafile;

  dnafile.open(inputFile);
  if(!dnafile.is_open())
  {
    exit(EXIT_FAILURE);
  }

  string word;
  dnafile >> word;
  while(dnafile.good())
  {
    makeLowerCase(word);
    addSum(word); //counts the number of dna nucleotides in file
    addDnaCount();
    dnafile >> word;
  }
  dnafile.close();
}
//takes each DNA string and counts how many nucleotides
void Dna::addSum(string word)
{
  sum+= word.size();
}
//returns the sum
int Dna::getSum()
{
  return sum;
}
//keeps track of how many dna sequences
void Dna::addDnaCount()
{
  dnaCount++;
}
//returns the number of DNA strings
int Dna::getDnaCount()
{
  return dnaCount;
}
//returns the mean
double Dna::getMean()
{
  return (double)getSum()/getDnaCount();
}
//keeps track of sum of squares
void Dna::addVariance()
{
  ifstream dnafile;
  dnafile.open(inputFile);

  if(!dnafile.is_open())
  {
    exit(EXIT_FAILURE);
  }

  string wordv;
  dnafile >> wordv;
  while(dnafile.good())
  {
    makeLowerCase(wordv);
    double i = wordv.size()-getMean();
    i = pow(i, 2);
    varianceSum += i;
    dnafile >> wordv;
  }
  dnafile.close();
}
//calculates and returns the variance based on the sum of squares
double Dna::getVariance()
{
  addVariance();
  finalVariance =  ((double)varianceSum)/(getDnaCount()-1);
  return finalVariance;
}
//returns the standard deviation of the file
double Dna::getStandDev()
{
  return sqrt(finalVariance);
}
//prints data to hankmoss.out.txt
void Dna::printSumStats()
{
  //opening file
  string nameOfOutputFile = "hankmoss.out.txt";
  cout<<"\n\nYour data will be printed to: "<<nameOfOutputFile<<endl<<endl;
  ofstream outputFile(nameOfOutputFile);
  outputFile<<"Hank Moss\n2252167\n\n";
  outputFile <<"Sum of length of DNA Strings: "<<getSum()<<endl;
  outputFile <<"Mean of length of DNA Strings: "<<getMean()<<endl;
  outputFile <<"Variance of length of DNA Strings: "<<getVariance()<<endl;
  outputFile <<"Standard Deviation of length of DNA Strings: "<<getStandDev()<<endl<<endl;
  compProbSingle();
  compProbPairs();
  outputFile.close();
}
//makes all of the letters in a string lower case
void Dna::makeLowerCase(string& word)
{
  for(int i = 0; i<word.length(); i++)
  {
    word[i] = tolower(word[i]);
  }
}
//makes a Char lower case
void Dna::makeCharLowerCase(char& letter)
{
  letter = tolower(letter);
}
//computes the probability of each single nucleotide in the inputFile and prints it to the outputFile
void Dna::compProbSingle()
{
  ifstream dnafile;
  dnafile.open(inputFile);

  if(!dnafile.is_open())
  {
    exit(EXIT_FAILURE);
  }

  char nucleo;
  dnafile >> nucleo;
  while(dnafile.good())
  {
    makeCharLowerCase(nucleo);
    if (nucleo == 'a')
    {
      aCount++;
    }
    else if (nucleo == 'c')
    {
      cCount++;
    }
    else if (nucleo == 'g')
    {
      gCount++;
    }
    else if (nucleo == 't')
    {
      tCount++;
    }
    else
    {
      cout<<"Error, mutation"<<endl;
    }
    dnafile >> nucleo;
  }
  dnafile.close();

  ofstream outputFile;
  outputFile.open("hankmoss.out.txt", ios::app);

  //prints out the probability of each
  outputFile <<"Probability of A's: "<<(double)aCount/getSum()<<endl;
  outputFile <<"Probability of C's: "<<(double)cCount/getSum()<<endl;
  outputFile <<"Probability of G's: "<<(double)gCount/getSum()<<endl;
  outputFile <<"Probability of T's: "<<(double)tCount/getSum()<<endl;

  outputFile.close();
}
//computes the probability of each nucleotide pair in the intputFile and prints it to the outputFile
//also prints 1000 nucleotide strings with same probabilities
void Dna::compProbPairs()
{
  //reseting nucleotide pair counts to 0
  int AAcount = 0;
  int ACcount = 0;
  int AGcount = 0;
  int ATcount = 0;
  int CAcount = 0;
  int CCcount = 0;
  int CGcount = 0;
  int CTcount = 0;
  int GAcount = 0;
  int GCcount = 0;
  int GGcount = 0;
  int GTcount = 0;
  int TAcount = 0;
  int TCcount = 0;
  int TGcount = 0;
  int TTcount = 0;

  ifstream dnafile;
  dnafile.open(inputFile);

  if(!dnafile.is_open())
  {
    exit(EXIT_FAILURE);
  }

  //retrieves the first two nucleotides
  char nucleo1;
  char nucleo2;
  dnafile >> nucleo1;
  dnafile >> nucleo2;
  while(dnafile.good())
  {
    makeCharLowerCase(nucleo1);
    makeCharLowerCase(nucleo2);
    //sorts and counts how many of each nucleotide pairs there are
    if (nucleo1 == 'a')
    {
      if(nucleo2 == 'a')
        AAcount++;
      else if(nucleo2 == 'c')
        ACcount++;
      else if(nucleo2 == 'g')
        AGcount++;
      else if(nucleo2 == 't')
        ATcount++;
    }
    else if (nucleo1 == 'c')
    {
      if(nucleo2 == 'a')
        CAcount++;
      else if(nucleo2 == 'c')
        CCcount++;
      else if(nucleo2 == 'g')
        CGcount++;
      else if(nucleo2 == 't')
        CTcount++;
    }
    else if (nucleo1 == 'g')
    {
      if(nucleo2 == 'a')
        GAcount++;
      else if(nucleo2 == 'c')
        GCcount++;
      else if(nucleo2 == 'g')
        GGcount++;
      else if(nucleo2 == 't')
        GTcount++;
    }
    else if (nucleo1 == 't')
    {
      if(nucleo2 == 'a')
        TAcount++;
      else if(nucleo2 == 'c')
        TCcount++;
      else if(nucleo2 == 'g')
        TGcount++;
      else if(nucleo2 == 't')
        TTcount++;
    }
    else
    {
      cout<<"Error, mutation"<<endl;
    }
    nucleo1 = nucleo2;
    dnafile >> nucleo2;
  }
  dnafile.close();

  //calculates the probability of each pair based off of the total number of them
  double AAprob = (double)AAcount/(getSum()-1);
  double ACprob = (double)ACcount/(getSum()-1);
  double AGprob = (double)AGcount/(getSum()-1);
  double ATprob = (double)ATcount/(getSum()-1);
  double CAprob = (double)CAcount/(getSum()-1);
  double CCprob = (double)CCcount/(getSum()-1);
  double CGprob = (double)CGcount/(getSum()-1);
  double CTprob = (double)CTcount/(getSum()-1);
  double GAprob = (double)GAcount/(getSum()-1);
  double GCprob = (double)GCcount/(getSum()-1);
  double GGprob = (double)GGcount/(getSum()-1);
  double GTprob = (double)GTcount/(getSum()-1);
  double TAprob = (double)TAcount/(getSum()-1);
  double TCprob = (double)TCcount/(getSum()-1);
  double TGprob = (double)TGcount/(getSum()-1);
  double TTprob = (double)TTcount/(getSum()-1);

  //prints the probabilities to the outputFile
  ofstream outputFile;
  outputFile.open("hankmoss.out.txt", ios::app);
  outputFile << endl;
  outputFile <<"Probability of AA's: "<<AAprob<<endl;
  outputFile <<"Probability of AC's: "<<ACprob<<endl;
  outputFile <<"Probability of AG's: "<<AGprob<<endl;
  outputFile <<"Probability of AT's: "<<ATprob<<endl;
  outputFile <<"Probability of CA's: "<<CAprob<<endl;
  outputFile <<"Probability of CC's: "<<CCprob<<endl;
  outputFile <<"Probability of CG's: "<<CGprob<<endl;
  outputFile <<"Probability of CT's: "<<CTprob<<endl;
  outputFile <<"Probability of GA's: "<<GAprob<<endl;
  outputFile <<"Probability of GC's: "<<GCprob<<endl;
  outputFile <<"Probability of GG's: "<<GGprob<<endl;
  outputFile <<"Probability of GT's: "<<GTprob<<endl;
  outputFile <<"Probability of TA's: "<<TAprob<<endl;
  outputFile <<"Probability of TC's: "<<TCprob<<endl;
  outputFile <<"Probability of TG's: "<<TGprob<<endl;
  outputFile <<"Probability of TT's: "<<TTprob<<endl<<endl;

  //prints 1000 strings to outputFile
  int k = 0; //random constant
  for(int i = 0; i < 1000; ++i)
  {
    int length = findD(i);
    double randDecimal = 0;
    //chooses the pairs in the strings based on the probabilities already calculated
    for(int j = 0; j < length; ++j)
    {
      k++;
      srand(k);
      randDecimal = (rand() % 100000)/100000.0;

      double total = 0;
      if(randDecimal>=total && randDecimal<(total += AAprob))
      {
        outputFile<<"AA";
      }
      else if(randDecimal>=total && randDecimal<(total += ACprob))
      {
        outputFile<<"AC";
      }
      else if(randDecimal>=total && randDecimal<(total += AGprob))
      {
        outputFile<<"AG";
      }
      else if(randDecimal>=total && randDecimal<(total += ATprob))
      {
        outputFile<<"AT";
      }
      else if(randDecimal>=total && randDecimal<(total += CAprob))
      {
        outputFile<<"CA";
      }
      else if(randDecimal>=total && randDecimal<(total += CCprob))
      {
        outputFile<<"CC";
      }
      else if(randDecimal>=total && randDecimal<(total += CGprob))
      {
        outputFile<<"CG";
      }
      else if(randDecimal>=total && randDecimal<(total += CTprob))
      {
        outputFile<<"CT";
      }
      else if(randDecimal>=total && randDecimal<(total += GAprob))
      {
        outputFile<<"GA";
      }
      else if(randDecimal>=total && randDecimal<(total += GCprob))
      {
        outputFile<<"GC";
      }
      else if(randDecimal>=total && randDecimal<(total += GGprob))
      {
        outputFile<<"GG";
      }
      else if(randDecimal>=total && randDecimal<(total += GTprob))
      {
        outputFile<<"GT";
      }
      else if(randDecimal>=total && randDecimal<(total += TAprob))
      {
        outputFile<<"TA";
      }
      else if(randDecimal>=total && randDecimal<(total += TCprob))
      {
        outputFile<<"TC";
      }
      else if(randDecimal>=total && randDecimal<(total += TGprob))
      {
        outputFile<<"TG";
      }
      else if(randDecimal>=total && randDecimal<(total += TTprob))
      {
        outputFile<<"TT";
      }
      else
      {
        outputFile<<"Uh oh";
      }
    }
    outputFile<<endl;
  }
  outputFile.close();
}
//calculates and returns a postive D, the length of the string based on Gaussian distribution
double Dna::findD(int i)
{
  srand(i);
  double a = rand() / (RAND_MAX + 1.0);
  double b = rand() / (RAND_MAX + 1.0);
  double c = sqrt(-2*log(a))*cos(2*M_PI*b);
  double d = c*getStandDev()+getMean();

  return abs(d);
}
