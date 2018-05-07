/*
   Copyright: Itoga and Takeuchi
   License: MIT (for all files)
*/


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
         boost::algorithm::split(vs,argments.at(i),boost::is_any_of("="));
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
            boost::algorithm::split(vs_,vs.at(1),boost::is_any_of("-"));
            time_range=std::tuple<double,double>
            (
               boost::lexical_cast<double>(vs_.at(0)),
               boost::lexical_cast<double>(vs_.at(1))
            ); 
         }
      }
   }catch(...)
   {
      std::cout<<"Failure @input argments"<<std::endl;exit(0);
   }

   std::vector<Protein> proteins = load(input_file_name,time_range);
   std::cout<<proteins.size()<<std::endl;
   pbc_setup(proteins.front().cryst);

   if(proteins.size()<2){std::cout<<"The number of samples is insufficient."<<std::endl;exit(0);}

   const int first_res = proteins.at(0).get_index_of_first_resudue();
   const int last_res  = proteins.at(0).get_index_of_last_resudue();
   const int pos_target = target_residue_index-first_res;
   std::vector<Kahan> correlations(last_res-first_res+1);
   for(size_t i=0,i_size=proteins.size();i<i_size;++i)
   {
      std::vector<Vector3D> us;//include disp. of each residue
      for(int r=first_res;r<=last_res;++r)
      {
         us.push_back(displacement_of_residue(proteins.at(i),proteins.at(i+1),r));
      }
      //const int pos_target = target_residue_index-first_res;
      const Vector3D& u_i = us.at(pos_target);
      for(size_t r=0,r_size=us.size();r<r_size;++r)
      {
         const Vector3D& u_j = us.at(r);
         correlations.at(r) += correlation(u_i,u_j);
      }
   }
   std::ofstream ofs(output_file_name,std::ios::trunc);
   ofs<<"# Res.index <u_i*u_j>_t"<<std::endl;
   for(size_t i=0,i_size=correlations.size();i<i_size;++i)
   {
      if(i!=target_residue_index)
      {
         ofs<<i+first_res<<" "<<correlations.at(i).get_av();
      }else
      {
         ofs<<std::endl<<std::endl;
      }
   }
   ofs.close();
   return EXIT_SUCCESS;
}

Vector3D displacement_of_residue
(
   const Protein& protein_a,
   const Protein& protein_b,
   const int& index_target_residue
)
{
   const Vector3D b = protein_b.get_center_of_residue(index_target_residue);
   const Vector3D a = pbc(protein_a.get_center_of_residue(index_target_residue),a);
   return b-a;
}

double correlation(const Vector3D& a, const Vector3D& b)
{
   return a*b;   
}

