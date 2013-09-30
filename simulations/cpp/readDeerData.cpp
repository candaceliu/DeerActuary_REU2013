#include <fstream>
#include <iostream>
#include <string.h>


#define FILE_NAME "threaded_trial.dat"

int main(int argc,char **argv)
{

  std::ifstream theFile(FILE_NAME,std::ios::in|std::ios::binary);

  struct output 
  {
    double theTime;
    double P;
    double alpha;
    double m1;
    double m2;
    double sumX;
    double sumX2;
    double sumM;
    double sumM2;
    int iterations;
  };
  output basicInfo;
  char   buffer[1024];


  while(!theFile.eof())
    {
      theFile.read(buffer,sizeof(basicInfo));
      memcpy(&basicInfo,buffer,sizeof(basicInfo));
      std::cout << "Time: "        << basicInfo.theTime << std::endl
                << " P: "          << basicInfo.P << std::endl
                << " alpha: "      << basicInfo.alpha << std::endl
                << " m1: "         << basicInfo.m1 << std::endl
                << " m2: "         << basicInfo.m2 << std::endl
                << " sum x: "      << basicInfo.sumX << std::endl
                << " sum x2: "     << basicInfo.sumX2 << std::endl
                << " sum m: "      << basicInfo.sumM << std::endl
                << " sum m2: "     << basicInfo.sumM2 << std::endl
                << " iterations: " << basicInfo.iterations << std::endl;
    }

  theFile.close();
  return(0);
}
