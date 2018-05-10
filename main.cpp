/*
   Copyright: Itoga and Takeuchi
   License: MIT (for all files)
*/


#include "PDB.hpp"
#include "Kahan_summation_algorithm.hpp"

Vector3D displacement_of_residue(const Protein& protein_a,const Protein& protein_b,const int& index_target_residue);
double correlation(const Vector3D& a, const Vector3D& b);
void plot_ca_for_yuba(const std::string file_name, const int step, const Protein& pr);

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

   Getline gl(input_file_name);
   Protein master     = get_next(gl,time_range,true);
   Protein protein_m1 = master;
   Protein protein_m0;
   pbc_setup(protein_m1.cryst);


   const int first_res  = protein_m1.get_index_of_first_resudue();
   const int last_res   = protein_m1.get_index_of_last_resudue();
   const int pos_target = target_residue_index-first_res;
   std::vector<Kahan> correlations(last_res-first_res+1);
   std::vector<Kahan> thetas(last_res-first_res+1);
   int step=0;
   for(;;)
   {
      try{
      protein_m0 = get_next(gl,time_range);
      }catch(...){std::cout<<"end"<<std::endl;break;}

      std::vector<Vector3D> us;//include disp. of each residue
      for(int r=first_res;r<=last_res;++r)
      {
         const Vector3D b = protein_m0.get_center_of_residue(r);
         const Vector3D a = pbc(protein_m1.get_center_of_residue(r),b);
         const auto ab = b-a;
         us.push_back(ab);
      }
      const Vector3D& u_i = us.at(pos_target);
      for(size_t r=0,r_size=us.size();r<r_size;++r)
      {
         const Vector3D& u_j    = us.at(r);
         const auto      res_c  = correlation(u_i,u_j);
         correlations.at(r)    += res_c;
         thetas.at(r)          += std::acos(res_c/(u_i.norm()*u_j.norm()));
      }
      protein_m1 = protein_m0;
      //protein_m0.fit_to_(master,target_residue_index);
      //plot_ca_for_yuba("test.cas",step,protein_m0);
      std::cout<<"†"<<std::flush;
      ++step;
   }
   {
      std::ofstream ofs(output_file_name,std::ios::trunc);
      ofs<<"# Res.index <u_i*u_j>_t <theta>_t"<<std::endl;
      for(size_t i=0,i_size=correlations.size();i<i_size;++i)
      {
         if((i+first_res)!=target_residue_index)
         {
            ofs<<(i+first_res)<<" "<<correlations.at(i).get_av()<<" "<<thetas.at(i).get_av()<<std::endl;
         }else
         {
            ofs<<std::endl<<std::endl;
         }
      }
      ofs.close();
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
   const Vector3D b = protein_b.get_center_of_residue(index_target_residue);
   const Vector3D a = pbc(protein_a.get_center_of_residue(index_target_residue),b);
   return b-a;
}

double correlation(const Vector3D& u_i, const Vector3D& u_j)
{
   return u_i*u_j;   
}

void plot_ca_for_yuba(const std::string file_name, const int step, const Protein& pr)
{
   const std::vector<Vector3D> cas = pr.get_CAs();
   if(0==step)
   {
      std::ofstream ofs(file_name,std::ios::trunc);
      ofs<<"# Coordinate "<<step<<std::endl;
      for(size_t i=0,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i<<" "<<cas.at(i)<<std::endl;
      }
      ofs<<std::endl;
      ofs<<"# Bond "<<step<<std::endl;
      for(size_t i=1,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i-1<<" "<<i<<std::endl;
      }
      ofs<<std::endl<<std::endl;
      ofs.close();
   }
   else
   {
      std::ofstream ofs(file_name,std::ios::app);
      ofs<<"# Coordinate "<<step<<std::endl;
      for(size_t i=0,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i<<" "<<cas.at(i)<<std::endl;
      }
      ofs<<std::endl;
      ofs<<"# Bond "<<step<<std::endl;
      for(size_t i=1,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i-1<<" "<<i<<" "<<step<<std::endl;
      }
      ofs<<std::endl<<std::endl;
      ofs.close();
   }
}
