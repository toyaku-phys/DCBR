#include "PDB.hpp"
#include "Kahan_summation_algorithm.hpp"

Vector3D displacement_of_residue(const Protein& protein_a,const Protein& protein_b,const int& index_target_residue);
double correlation(const Vector3D& a, const Vector3D& b);

int main(int argc, char* argv[])
{
   std::vector<std::string> argments;
   for(int i=0;i<argc;++i){ argments.push_back(argv[i]); }

   std::string input_file_name;
   std::string output_file_name;
   int target_residue_index;
   std::tuple<double,double> time_range(-DBL_MAX,+DBL_MAX);

   try
   {
      for(size_t i=0,i_size=argments.size();i<i_size;++i)
      {
         std::vector<std::string> vs;
         boost::algorithm::split(vs,argments.at(i),boost::is_any_of(" "));
         if("in"==vs.at(0))
         {
            input_file_name=vs.at(1);
         }
         if("out"==vs.at(0))
         {
            output_file_name=vs.at(1);
         }
         if("target"==vs.at(0))
         {
            target_residue_index=boost::lexical_cast<double>(vs.at(1));
         }
         if("time"==vs.at(0))
         {
            std::vector<std::string> vs_;
            boost::algotrithm::split(vs_,vs.at(1),boost::is_any_of("-"));
            time_range=std::tuple<double,double>
            (
               boost::lexical_cast<double>(vs_.at(0)),
               boost::lexical_cast<double>(vs_.at(1))
            ); 
         }
      }
   }catch(...)
   {
      std::cout<<"Failure @input argments"<<std::endl;   
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
