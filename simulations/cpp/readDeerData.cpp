#include <fstream>
#include <iostream>
#include <string.h>


#define FILE_NAME "threaded_trial.dat"

int main(int argc,char **argv)
{

  std::ifstream theFile(FILE_NAME,std::ios::in|std::ios::binary);
  int lupe = 0;

  struct output 
  {
    double theTime;
    double P;
    double alpha;
    int lupeAlpha;
    int lupeP;
    double sumX;
    double sumX2;
    double sumM;
    double sumM2;
    int iterations;
  };
  output basicInfo;
  char   buffer[1024];


  //std::cout << sizeof(basicInfo) << std::endl;
  theFile.seekg(0);
  while(!theFile.eof())
    {
      //theFile.seekg((lupe++)*80);
      theFile.read(buffer,sizeof(basicInfo));
      memcpy(&basicInfo,buffer,sizeof(basicInfo));
      std::cout << "Time: "        << basicInfo.theTime << std::endl
                << " P: "          << basicInfo.P << std::endl
                << " alpha: "      << basicInfo.alpha << std::endl
                << " Lalpha: "     << basicInfo.lupeAlpha << std::endl
                << " LP: "         << basicInfo.lupeP << std::endl
                << " sum x: "      << basicInfo.sumX << std::endl
                << " sum x2: "     << basicInfo.sumX2 << std::endl
                << " sum m: "      << basicInfo.sumM << std::endl
                << " sum m2: "     << basicInfo.sumM2 << std::endl
                << " iterations: " << basicInfo.iterations << std::endl
                << " loop: "       << lupe << std::endl << std::endl;
    }

  theFile.close();
  return(0);
}
