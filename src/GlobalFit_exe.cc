#include "GlobalFit.hpp"
using namespace std;

int main(){
  PrintLine("MAIN", blue);
  std::unique_ptr<GlobalFit> GF(new GlobalFit("Run2Test", "ratio", 0.0, 1.3));
  GF->Run();
  return 0;
}
