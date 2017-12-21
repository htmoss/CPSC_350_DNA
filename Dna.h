#include<string>
using namespace std;

class Dna
{
private:
  int sum;
  int dnaCount;
  double varianceSum;
  double finalVariance;
  double standDev;
  string inputFile;
  int aCount;
  int cCount;
  int gCount;
  int tCount;
  double AAprob;
  double ACprob;
  double AGprob;
  double ATprob;
  double CAprob;
  double CCprob;
  double CGprob;
  double CTprob;
  double GAprob;
  double GCprob;
  double GGprob;
  double GTprob;
  double TAprob;
  double TCprob;
  double TGprob;
  double TTprob;

public:
  Dna();
  ~Dna();
  void compMean();
  void compVariance();
  void compStanDev();
  void readFile();
  void addSum(string word);
  int getSum();
  void addDnaCount();
  int getDnaCount();
  double getMean();
  void addVariance();
  double getVariance();
  void printSumStats();
  void makeLowerCase(string& word);
  void makeCharLowerCase(char& word);
  double getStandDev();
  void compProbSingle();
  void compProbPairs();
  double findD(int i);
  void clearVariables();

};
