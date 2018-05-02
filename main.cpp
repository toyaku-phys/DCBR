#include "PDB.hpp"


int main(int argc, char* argv[])
{
   std::vector<std::string> argments;
   for(int i=0;i<argc;++i){ argments.push_back(argv[i]); }

   std::string input_file_name;

   for(size_t i=0,i_size=argments.size();i<i_size;++i)
   {
      std::vector<std::string> vs;
      boost::algorithm::split(vs,argments.at(i));
      if("in"==vs.at(0))
      {
         input_file_name=vs.at(1);
      }
   }

   return EXIT_SUCCESS;
}
