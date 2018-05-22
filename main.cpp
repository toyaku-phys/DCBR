/*
   Copyright: Itoga and Takeuchi
   License: MIT (for all files)
*/


#include "PDB.hpp"
#include "Kahan_summation_algorithm.hpp"

Vector3D displacement_of_residue(const Protein& protein_a,const Protein& protein_b,const int& index_target_residue);
Vector3D velocity_of_residue
(const Protein& protein_a,const Protein& protein_b,const int& index_target_residue);
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

   constexpr int MAX_INTERVAL=10;
   int first_res;
   std::vector<std::vector<Kahan> > correlationss;
   for(int interval=0;interval<=MAX_INTERVAL;++interval)
   {
      Getline gl(input_file_name);
      Protein master     = get_next(gl,time_range,true);
      Protein protein_m1 = master;
      Protein protein_m0;
      pbc_setup(protein_m1.cryst);

                first_res  = protein_m1.get_index_of_first_resudue();
      const int last_res   = protein_m1.get_index_of_last_resudue();
      const int pos_target = target_residue_index-first_res;
      std::vector<Kahan> correlations(last_res-first_res+1);
      int step=0;
      for(;;)
      {
         try{
            for(int s=0;s<interval;++s)
            {protein_m0 = get_next(gl,time_range);}
         }catch(...){std::cout<<"end"<<std::endl;break;}

         std::vector<Vector3D> vs;//include disp. of each residue
         for(int r=first_res;r<=last_res;++r)
         {
            vs.push_back(velocity_of_residue(protein_m0,protein_m1,r));
         }
         const Vector3D& v_i = vs.at(pos_target);
         for(size_t r=0,r_size=vs.size();r<r_size;++r)
         {
            const Vector3D& v_j    = vs.at(r);
            const auto      res_c  = correlation(v_i,v_j)/(v_i.norm()*v_j.norm());
            correlations.at(r)    += res_c;
         }
         protein_m1 = protein_m0;
         //protein_m0.fit_to_(master,target_residue_index);
         //plot_ca_for_yuba("test.cas",step,protein_m0);
         std::cout<<"†"<<std::flush;
         ++step;
      }//end of for(;;)
      correlationss.push_back(correlations);
   }//end of for(interval)
   {
      std::ofstream ofs(output_file_name,std::ios::trunc);
      ofs<<"# Res.index <(v_i(0)*v_j(t))/sqrt{v_i^2(0) vj^2(t)}>_t"<<std::endl;
      size_t size = correlationss.front().size();
      for(size_t r=0;r<size;++r)
      {
         ofs<<(r+first_res)<<" ";
         for(size_t i=0,i_size=correlationss.size();i<i_size;++i)
         {
            ofs<<correlationss.at(i).at(r).get_av()<<" ";
         }
         ofs<<std::endl;
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

Vector3D velocity_of_residue
(
   const Protein& protein_a,
   const Protein& protein_b,
   const int& index_target_residue
)
{
   const Vector3D b = protein_b.get_center_of_residue(index_target_residue);
   const Vector3D a = pbc(protein_a.get_center_of_residue(index_target_residue),b);
   return (b-a)/(protein_b.time-protein_a.time);
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
