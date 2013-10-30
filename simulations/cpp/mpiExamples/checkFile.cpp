#include <fstream>
#include <iostream>
#include <string.h>


#define FILE_NAME "fileExample-01.dat"

int main(int argc,char **argv)
{

  int i;
  double x;
  std::ifstream theFile(FILE_NAME,std::ios::in|std::ios::binary);

  struct output 
  {
    double x;
    int    i;
  };
  output basicInfo;
  char   buffer[1024];


  /*
  while(!theFile.eof())
    {
      char c;
      theFile.read(&c,sizeof(c));
      std::cout << c;
    }
  std::cout << std::endl;
  */

  theFile.read(buffer,sizeof(basicInfo));
  memcpy(&basicInfo,buffer,sizeof(basicInfo));
  std::cout << "Int: " << basicInfo.i << "/" << sizeof(i) 
            << " Double: " << basicInfo.x << "/" << sizeof(x) << std::endl;

  theFile.read(buffer,sizeof(basicInfo));
  memcpy(&basicInfo,buffer,sizeof(basicInfo));
  std::cout << "Int: " << basicInfo.i << "/" << sizeof(i) 
            << " Double: " << basicInfo.x << "/" << sizeof(x) << std::endl;

  theFile.close();
  return(0);
}
