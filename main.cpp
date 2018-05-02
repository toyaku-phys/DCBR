#include "PDB.hpp"

Vector3D displacement_of_residue(const Protein& protein_a,const Protein& protein_b,const int& index_target_residue);
double correlation(const Vector3D& a, const Vector3D& b);

int main(int argc, char* argv[])
{
   std::vector<std::string> argments;
   for(int i=0;i<argc;++i){ argments.push_back(argv[i]); }

   std::string input_file_name;

   for(size_t i=0,i_size=argments.size();i<i_size;++i)
   {
      std::vector<std::string> vs;
      boost::algorithm::split(vs,argments.at(i),boost::is_any_of(" "));
      if("in"==vs.at(0))
      {
         input_file_name=vs.at(1);
      }
   }

   std::vector<Protein> proteins = load(input_file_name);

   for(size_t i=0,i_size=proteins.size();i<i_size;++i)
   {
      
   }

   return EXIT_SUCCESS;
}

Vector3D displacement_of_residue
(
   const Protein& protein_a,
   const Protein& protein_b,
   const int& index_target_residue
)
{
   return
   protein_b.get_center_of_residue(index_target_residue)
   -protein_a.get_center_of_residue(index_target_residue);
}

double correlation(const Vector3D& a, const Vector3D& b)
{
   return a*b;   
}
